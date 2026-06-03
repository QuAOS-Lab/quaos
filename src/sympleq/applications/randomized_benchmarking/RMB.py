from __future__ import annotations
import copy
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
from sympleq.core.noise.noise_model import GenericNoise
from sympleq.integrations.quantinuum.utils import NATIVE_GATES_SET

DATA_DIR = Path(__file__).parent / "rmb_data"


class RMB:
    def __init__(self,
                 update_strategy: UpdateStrategy,
                 backend: RMBBackend,
                 rng: RNGGenerator,
                 default_estimator: BayesianEstimator | None = None) -> None:
        self._update_strategy = update_strategy
        self.backend = backend
        self._default_estimator = default_estimator

        self.rng = rng
        self._data: RMBData = {}
        self._save_filename = ""

    def new_default_estimator(self) -> BayesianEstimator:
        if self._default_estimator is None:
            return self.backend.default_estimator()
        return copy.deepcopy(self._default_estimator)

    def with_backend(self, backend: RMBBackend) -> RMB:
        """Set ``self.backend`` and return ``self`` for chaining."""
        self.backend = backend
        return self

    def with_update_strategy(self, update_strategy: UpdateStrategy) -> RMB:
        """Set ``self._update_strategy`` and return ``self`` for chaining."""
        self._update_strategy = update_strategy
        return self

    def with_bayesian_estimator(self, default_estimator: BayesianEstimator) -> RMB:
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
        self._default_estimator = default_estimator
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

            estimator = self._data.setdefault(config, self.new_default_estimator())
            estimator.run(
                lambda c=config: self.backend.fidelity_estimation(c, self.rng),
                verbose,
            )
            self.save()

            # Pick next config, i.e. pick next point in parameter space, based on the update_strategy
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

    def update(self, data: RMBData):
        """Merge ``data`` into ``self._data``.

        For configs already present, the incoming estimator is merged into
        the existing one; new configs are added with a fresh estimator
        seeded from the incoming one's parameters.

        Parameters
        ----------
        data : RMBData
            Mapping from configurations to their Bayesian estimators.
        """
        for config, estimator in data.items():
            target = self._data.setdefault(config, BayesianEstimator(
                threshold=estimator.threshold,
                min_runs=estimator.min_runs,
                max_runs=estimator.max_runs,
            ))
            target.merge(estimator)

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

        estimator = self.new_default_estimator()
        payload = {
            "backend": self.backend.to_dict(),
            "estimator": {
                "threshold": estimator.threshold,
                "min_runs": estimator.min_runs,
                "max_runs": estimator.max_runs,
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
        estimator = BayesianEstimator(estimator_params["threshold"],
                                      estimator_params["min_runs"], estimator_params["max_runs"])
        rmb = cls(update_strategy=default_update_strategy,
                  backend=backend,
                  rng=rng,
                  default_estimator=estimator)

        for rec in payload["data"]:
            config = RMBConfig(
                depth=rec["depth"],
                scrambling_probability=rec["scrambling_probability"],
                min_two_qubit_gate_ratio=rec["min_two_qubit_gate_ratio"],
                max_two_qubit_gate_ratio=rec["max_two_qubit_gate_ratio"],
                n_qubits=rec["n_qubits"],
                random_elimination=rec["random_elimination"],
            )
            estimator = rmb.backend.default_estimator()
            for outcome, count in rec["results"]:
                estimator.record(outcome, count)
            rmb._data[config] = estimator
        return rmb

    def plot(self, axes=None, show: bool = True, skip_incomplete: bool = True):
        RMB.plot_data(self._data, axes, show, skip_incomplete)

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

        data = RMB.merge_close_configs(data, depth_bin=20, ratio_digits=2)

        groups: dict[int, RMBData] = {}
        for config, estimator in data.items():
            # if config.n_qubits != 4:
            # continue
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

            ax.set_xlabel("# Gates")
            ax.set_ylabel("Two-qudit gate ratio")
            ax.set_title(f"# Qubits = {n_qubits}")
            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label("Fidelity")

        if show:
            plt.show()

        return axes_list


def sympleq_pipeline() -> RMB:
    rng = default_rng()
    noise_model = GenericNoise.from_paulis([0.000075, 0.000075, 0.000075], rng=rng)
    two_qubit_noise_model = GenericNoise.from_paulis([0.00039, 0.00039, 0.00039], rng=rng)
    backend = SympleqBackend(noise_model, two_qubit_noise_model)
    rmb = RMB.default(rng).with_backend(backend)

    for n_qubits in range(2, 9):
        if n_qubits != 5:
            continue
        folder_name = f"synthetic-{n_qubits}q"
        # for depth in range(120 - 6 * n_qubits, 170 - 6 * n_qubits, 2):
        #     initial_config = RMBConfig.default()\
        #         .with_n_qubits(n_qubits)\
        #         .with_depth(depth)\
        #         .with_two_qubit_gate_ratio(0.2, 0.4)\
        #         .with_scrambling_probability(0.5)\
        #         .with_gates_set(tuple(NATIVE_GATES_SET))
        #     generate_random_pytket_circuits(initial_config, 20, folder_name)
        data = backend.simulate_pytket_circuits(folder_name)
        rmb.update(data)
        rmb.save(f"{folder_name}_data.json")
        # _rmb = RMB.load(f"{folder_name}_data.json")

    # data = backend.populate_from_recent_jobs()
    rmb.save("sympleq_from_fetched.json")
    return rmb


def quantinuum_pipeline() -> RMB:
    backend = QuantinuumBackend(device_name="H2-Emulator")
    rmb = RMB.default().with_backend(backend)
    data = backend.populate_from_recent_jobs(1000)
    rmb.update(data)
    rmb.save("quantinuum_fetched.json")
    return rmb


def filter_data(data: RMBData) -> RMBData:
    merged_data = RMB.merge_close_configs(data, depth_bin=20, ratio_digits=1)
    filtered_data: RMBData = {}
    for config, estimator in sorted(
            merged_data.items(),
            key=lambda item: (item[0].depth, item[0].min_two_qubit_gate_ratio)):
        # if config.n_qubits != 5:
        #     continue
        # if config.depth >= 1000:
        #     continue
        if not (estimator.variance(True) <= 1e-1 and estimator._counts_tot > 10):
            continue
        filtered_data[config] = estimator
    return filtered_data


def compare():

    # qrmb = quantinuum_pipeline()
    # srmb = RMB.load("sympleq_from_fetched.json")
    qrmb = RMB.load("quantinuum_fetched.json")
    q_data = filter_data(qrmb._data)
    print("Quantinuum jobs loaded.")
    srmb = sympleq_pipeline()
    print("Sympleq jobs loaded.")
    s_data = filter_data(srmb._data)
    q_data = filter_data(qrmb._data)

    shared = s_data.keys() & q_data.keys()
    diffs = [
        (c, s_data[c].probability(True) - q_data[c].probability(True))
        for c in sorted(shared, key=lambda c: (c.depth, c.min_two_qubit_gate_ratio))
    ]
    fom = sum(d * d for _, d in diffs)
    print(f"shared configs: {len(shared)}; sum of squared fidelity differences: {fom:.6f}")
    for c, d in diffs:
        print(f"  depth={c.depth:>4} 2qr={c.min_two_qubit_gate_ratio:.2f} -> Δ={d:+.4f}")


if __name__ == "__main__":
    # srmb = sympleq_pipeline()

    rmb = RMB.default().with_backend(QuantinuumBackend(device_name="H2-Emulator", batch_size=1))
    initial_config = RMBConfig.default()\
        .with_n_qubits(2)\
        .with_depth(30)\
        .with_two_qubit_gate_ratio(0.4, 0.6)\
        .with_scrambling_probability(0.5)\
        .with_gates_set(tuple(NATIVE_GATES_SET))

    rmb.run(initial_config)
