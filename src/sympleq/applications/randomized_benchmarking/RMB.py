from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
import datetime

from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.backends import RMBBackend, SympleqBackend, QuantinuumBackend
from sympleq.applications.randomized_benchmarking.backends.base import backend_from_dict
from sympleq.applications.randomized_benchmarking.update_strategy import UpdateStrategy, default_update_strategy
from sympleq.core.bayesian_estimation import BayesianEstimator

DATA_DIR = Path(__file__).parent / "rmb_data"


class RMB:
    def __init__(self,
                 update_strategy: UpdateStrategy,
                 backend: RMBBackend,
                 rng: RNGGenerator,
                 threshold: float = 10**(-3),
                 min_runs: int = 100,
                 max_runs: int | None = None) -> None:
        self._update_strategy = update_strategy
        self.backend = backend
        self._threshold = threshold
        self._min_runs = min_runs
        self._max_runs = max_runs
        self.rng = rng
        self._data: RMBData = {}
        self._save_filename = ""

    def _default_bayesian_estimator(self) -> BayesianEstimator:
        return BayesianEstimator(threshold=self._threshold, min_runs=self._min_runs, max_runs=self._max_runs)

    def with_backend(self, backend: RMBBackend) -> RMB:
        """Set ``self.backend`` and return ``self`` for chaining."""
        self.backend = backend
        return self

    def with_update_strategy(self, update_strategy: UpdateStrategy) -> RMB:
        """Set ``self._update_strategy`` and return ``self`` for chaining."""
        self._update_strategy = update_strategy
        return self

    def with_bayesian_estimator(self,
                                threshold: float = 10**(-3),
                                min_runs: int = 100,
                                max_runs: int | None = None) -> RMB:
        """Set the estimator parameters and return ``self`` for chaining.

        Parameters
        ----------
        threshold : float
            Variance threshold under which an estimator is considered converged.
        min_runs : int
            Minimum number of samples before convergence may be declared.
        max_runs : int | None
            Hard cap on samples; ``None`` means unlimited.

        Returns
        -------
        RMB
            ``self``.
        """
        self._threshold = threshold
        self._min_runs = min_runs
        self._max_runs = max_runs
        return self

    @classmethod
    def default(cls, rng: RNGGenerator | None = None) -> RMB:
        """Return an RMB with the default update strategy and a noiseless ``SympleqBackend``.

        Parameters
        ----------
        rng : numpy.random.Generator | None
            RNG to attach. ``None`` uses ``default_rng()``.

        Returns
        -------
        RMB
            A new RMB instance.
        """
        if rng is None:
            rng = default_rng()
        return cls(update_strategy=default_update_strategy,
                   backend=SympleqBackend(),
                   rng=rng)

    def run(self,
            config: RMBConfig | None = None,
            verbose: bool = False,
            max_iterations: int = 100):
        """Run the RMB sweep for up to ``max_iterations`` configurations.

        Parameters
        ----------
        config : RMBConfig | None
            Starting configuration. ``None`` falls back to ``RMBConfig.default()``.
        verbose : bool
            If ``True``, the underlying estimator prints live progress.
        max_iterations : int
            Maximum number of estimator-convergence iterations to run in this
            call.
        """
        self._save_filename = f"{self.backend.type}-{datetime.datetime.now():%Y-%m-%dT%H-%M-%S}.json"

        if config is None:
            config = RMBConfig.default()

        for _ in range(max_iterations):
            seen = {config}
            while config in self._data and self._data[config].is_converged():
                config = self._update_strategy(self._data, config)
                if config in seen:
                    print("All reachable configs are already converged; stopping.")
                    RMB.print_data(self._data)
                    return
                seen.add(config)

            estimator = self._data.setdefault(config, self._default_bayesian_estimator())
            print("Running with config", config)
            estimator.run(
                lambda c=config: self.backend.fidelity_estimation(c, self.rng),
                verbose,
            )
            self.save()
            config = self._update_strategy(self._data, config)

        RMB.print_data(self._data)

    @classmethod
    def merge_close_configs(cls,
                            data: RMBData,
                            depth_bin: int = 10,
                            ratio_digits: int = 1) -> RMBData:
        """Return a copy of ``data`` with close configs collapsed and their estimators merged.

        Parameters
        ----------
        data : RMBData
            Mapping from configurations to their Bayesian estimators.
        depth_bin : int
            Configs' ``depth`` is snapped to the nearest multiple of this value
            (with a floor of ``depth_bin`` so the dataclass ``depth >= 1``
            validation is preserved).
        ratio_digits : int
            Number of decimal places to round the two-qudit ratio bounds to.

        Returns
        -------
        RMBData
            New mapping keyed by coarsened configs.
        """
        merged: RMBData = {}
        for config, estimator in data.items():
            coarse_depth = max(depth_bin, round(config.depth / depth_bin) * depth_bin)
            coarse = replace(
                config,
                depth=coarse_depth,
                min_two_qubit_gate_ratio=round(config.min_two_qubit_gate_ratio, ratio_digits),
                max_two_qubit_gate_ratio=round(config.max_two_qubit_gate_ratio, ratio_digits),
            )
            target = merged.setdefault(coarse, BayesianEstimator(
                threshold=estimator.threshold,
                min_runs=estimator.min_runs,
                max_runs=estimator.max_runs,
            ))
            target.merge(estimator)

        return merged

    @classmethod
    def print_data(cls, data: RMBData):
        """Print ``data`` as a markdown-style table of configs and fidelities."""
        headers = ["depth", "two_qubit_gate_ratio", "n_qubits", "random_elimination", "fidelity"]
        rows = []
        for config, estimator in data.items():
            ratio = (
                f"{config.min_two_qubit_gate_ratio}"
                if config.min_two_qubit_gate_ratio == config.max_two_qubit_gate_ratio
                else f"[{config.min_two_qubit_gate_ratio}, {config.max_two_qubit_gate_ratio}]"
            )
            rows.append([
                f"{config.depth}",
                ratio,
                f"{config.n_qubits}",
                f"{config.random_elimination}",
                f"{estimator.probability(True):.4f} ± {estimator.variance(True):.4f}",
            ])

        widths = [max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
                  for i, h in enumerate(headers)]

        def fmt(cells: list[str]) -> str:
            return "| " + " | ".join(c.center(w) for c, w in zip(cells, widths)) + " |"

        sep = "|-" + "-|-".join("-" * w for w in widths) + "-|"
        print(sep)
        print(fmt(headers))
        print(sep)
        for row in rows:
            print(fmt(row))
        print(sep)

    def save(self, path: str | Path | None = None):
        """Save this RMB run as JSON.

        Parameters
        ----------
        path : str | Path | None
            Output file. Bare file names (no separators) are resolved
            inside the package's ``rmb_data/`` directory; otherwise the given
            path is used as-is. ``None`` falls back to ``self._save_filename``.
        """
        import json
        if path is None:
            path = self._save_filename

        path = Path(path)
        if not path.is_absolute() and path.parent == Path("."):
            path = DATA_DIR / path
        path.parent.mkdir(parents=True, exist_ok=True)

        records = []
        for config, estimator in self._data.items():
            records.append({
                "depth": config.depth,
                "scrambling_probability": config.scrambling_probability,
                "min_two_qubit_gate_ratio": config.min_two_qubit_gate_ratio,
                "max_two_qubit_gate_ratio": config.max_two_qubit_gate_ratio,
                "n_qubits": config.n_qubits,
                "random_elimination": config.random_elimination,
                "gates_set": [g.name for g in config.gates_set],
                "results": [[outcome, count] for outcome, count in estimator.counts().items()],
            })

        payload = {
            "backend": self.backend.to_dict(),
            "estimator": {
                "threshold": self._threshold,
                "min_runs": self._min_runs,
                "max_runs": self._max_runs,
            },
            "data": records,
        }

        with open(path, "w") as f:
            json.dump(payload, f, indent=2)

    @classmethod
    def load(cls, path: str | Path = "rmb.json", rng: RNGGenerator | None = None) -> RMB:
        """Load an RMB run previously written by :meth:`save`.

        ``gates_set`` is not reconstructed and falls back to the
        :class:`RMBConfig` default.

        Parameters
        ----------
        path : str | Path
            Input file. Bare file names (no separators) are resolved
            inside the package's ``data/`` directory; otherwise the given
            path is used as-is.
        rng : numpy.random.Generator | None
            RNG to attach to the rebuilt RMB. ``None`` uses ``default_rng()``.

        Returns
        -------
        RMB
            The reconstructed RMB instance.
        """
        import json
        path = Path(path)
        if not path.is_absolute() and path.parent == Path("."):
            path = DATA_DIR / path
        with open(path) as f:
            payload = json.load(f)

        if rng is None:
            rng = default_rng()

        backend = backend_from_dict(payload["backend"])
        estimator_params = payload["estimator"]
        rmb = cls(update_strategy=default_update_strategy,
                  backend=backend,
                  rng=rng,
                  threshold=estimator_params["threshold"],
                  min_runs=estimator_params["min_runs"],
                  max_runs=estimator_params["max_runs"])

        for rec in payload["data"]:
            config = RMBConfig(
                depth=rec["depth"],
                scrambling_probability=rec["scrambling_probability"],
                min_two_qubit_gate_ratio=rec["min_two_qubit_gate_ratio"],
                max_two_qubit_gate_ratio=rec["max_two_qubit_gate_ratio"],
                n_qubits=rec["n_qubits"],
                random_elimination=rec["random_elimination"],
            )
            estimator = rmb._default_bayesian_estimator()
            for outcome, count in rec["results"]:
                estimator.record(outcome, count)
            rmb._data[config] = estimator
        return rmb

    @classmethod
    def plot_data(cls, data: RMBData, axes=None, show: bool = True, skip_incomplete: bool = True):
        """Scatter fidelity per ``n_qubits``: x = depth, y = two-qudit ratio, color = fidelity.

        Parameters
        ----------
        data : RMBData
            Mapping from configurations to their Bayesian estimators.
        axes : Sequence[matplotlib.axes.Axes] | None
            Axes to plot on, one per distinct ``n_qubits`` (sorted
            ascending). If ``None``, a new figure with one subplot per
            ``n_qubits`` is created.
        show : bool
            If ``True``, call ``plt.show()`` after building the plot.
        skip_incomplete : bool
            If ``True``, drop any estimator that has not converged.
            ``n_qubits`` groups left empty after filtering are not given a
            subplot.

        Returns
        -------
        list[matplotlib.axes.Axes]
            The axes the plots were drawn on (sorted by ``n_qubits``).
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap

        data = RMB.merge_close_configs(data, depth_bin=25, ratio_digits=1)

        groups: dict[int, RMBData] = {}
        for config, estimator in data.items():
            if config.n_qubits != 5:
                continue
            # estimator.min_runs = 6
            if skip_incomplete and not estimator.is_converged():
                continue
            groups.setdefault(config.n_qubits, {})[config] = estimator

        if not groups:
            return []

        sorted_groups = sorted(groups.items())
        n_groups = len(sorted_groups)

        if axes is None:
            _, axes_arr = plt.subplots(
                1, n_groups, figsize=(5 * n_groups, 4), squeeze=False)
            axes_list = list(axes_arr[0])
        else:
            axes_list = list(axes)
            if len(axes_list) < n_groups:
                raise ValueError(
                    f"Need at least {n_groups} axes for {n_groups} n_qubits "
                    f"groups, got {len(axes_list)}.")

        cmap = LinearSegmentedColormap.from_list(
            "darkred_to_lime",
            ["darkred", "red", "orange", "lime", "green"])

        for ax, (n_qubits, group) in zip(axes_list, sorted_groups):
            depths = np.array([c.depth for c in group])
            ratios = np.array(
                [0.5 * (c.min_two_qubit_gate_ratio + c.max_two_qubit_gate_ratio) for c in group])
            fidelities = np.array([e.probability(True) for e in group.values()])
            stds = np.array([np.sqrt(e.variance(True)) for e in group.values()])

            outer_color = np.clip(fidelities - stds, 0.0, 1.0)
            middle_color = fidelities
            inner_color = np.clip(fidelities + stds, 0.0, 1.0)

            ax.scatter(depths, ratios, c=outer_color, cmap=cmap,
                       vmin=0.0, vmax=1.0, s=200, edgecolors="none", zorder=1)
            sc = ax.scatter(depths, ratios, c=middle_color, cmap=cmap,
                            vmin=0.0, vmax=1.0, s=100, edgecolors="none", zorder=2)
            ax.scatter(depths, ratios, c=inner_color, cmap=cmap,
                       vmin=0.0, vmax=1.0, s=30, edgecolors="none", zorder=3)

            ax.set_xlabel("depth")
            ax.set_ylabel("two-qudit gate ratio")
            ax.set_title(f"n_qubits = {n_qubits}")
            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label("fidelity")

        if show:
            plt.show()

        return axes_list


if __name__ == "__main__":
    rng = default_rng()
    verbose = True

    # noise_model = GenericNoise.from_paulis([0.001, 0.001, 0.002], rng=rng)
    # two_qubit_noise_model = GenericNoise.from_paulis([0.0075, 0.0075, 0.0075], rng=rng)
    # rmb = RMB.default(rng)\
    #     .with_backend(SympleqBackend(noise_model, two_qubit_noise_model))\
    #     .with_bayesian_estimator(threshold=10**(-3), min_runs=100)

    backend = QuantinuumBackend(device_name="H2-Emulator")
    # rmb = RMB.default(rng)\
    #     .with_backend(backend)\
    #     .with_bayesian_estimator(threshold=10**(-1), min_runs=backend.batch_size)
    # data = backend.populate_from_recent_jobs(1000)
    # rmb._data = data
    # rmb.save("quantinuum_fetched.json")
    rmb = RMB.load("quantinuum_fetched.json")
    # data = RMB.merge_close_configs(rmb._data)
    # print(len(data))
    # for v in data.values():
    #     print(v._probabilities)
    #     print(v._variances)
    #     print(v.counts())
    #     print()
    # initial_config = rmb.backend.default_config().with_depth(90).with_n_qubits(5).with_two_qubit_gate_ratio(0.3, 0.4)
    # rmb.run(initial_config, verbose=True, max_iterations=40)

    RMB.plot_data(rmb._data, skip_incomplete=True)
