from __future__ import annotations
from pathlib import Path
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData, UpdateStrategy
from sympleq.applications.randomized_benchmarking.backends import RMBBackend, SympleqBackend, QuantinuumBackend
from sympleq.applications.randomized_benchmarking.backends.base import backend_from_dict
from sympleq.core.bayesian_estimation import BayesianEstimator

DATA_DIR = Path(__file__).parent / "data"


class RMB:
    def __init__(self,
                 update_strategy: UpdateStrategy,
                 backend: RMBBackend,
                 rng: RNGGenerator,
                 *,
                 threshold: float = 10**(-3),
                 min_runs: int = 100,
                 max_runs: int | None = None) -> None:
        self._update_strategy = update_strategy
        self._backend = backend
        self._threshold = threshold
        self._min_runs = min_runs
        self._max_runs = max_runs
        self.rng = rng
        self._data: RMBData = {}

    def _default_bayesian_estimator(self) -> BayesianEstimator:
        return BayesianEstimator(threshold=self._threshold, min_runs=self._min_runs, max_runs=self._max_runs)

    def with_backend(self, backend: RMBBackend) -> RMB:
        """Set ``self._backend`` and return ``self`` for chaining."""
        self._backend = backend
        return self

    def with_update_strategy(self, update_strategy: UpdateStrategy) -> RMB:
        """Set ``self._update_strategy`` and return ``self`` for chaining."""
        self._update_strategy = update_strategy
        return self

    def with_bayesian_estimator(self,
                                threshold: float = 10**(-3),
                                min_runs: int = 100,
                                max_runs: int | None = None) -> RMB:
        """
        Set the estimator parameters and return ``self`` for chaining.

        These parameters are forwarded to :meth:`default_bayesian_estimator`
        when a fresh estimator is built for each configuration.
        """
        self._threshold = threshold
        self._min_runs = min_runs
        self._max_runs = max_runs
        return self

    @classmethod
    def default(cls, rng: RNGGenerator | None = None) -> RMB:
        """
        Build an RMB instance with default settings.

        Uses :meth:`RMBConfig.default_update_strategy` as the update
        strategy and a noiseless :class:`SympleqBackend` as the backend.
        """
        if rng is None:
            rng = default_rng()
        return cls(update_strategy=RMBConfig.default_update_strategy,
                   backend=SympleqBackend(),
                   rng=rng)

    def run(self, config: RMBConfig, verbose: bool = False):
        """
        Run the RMB sweep starting from ``config``.

        For each configuration, samples random circuits and feeds the
        success/failure outcomes (provided by the backend) into a
        Bayesian estimator until it converges.

        Parameters
        ----------
        config : RMBConfig
            Starting configuration. If ``None``, falls back to
            :meth:`RMBConfig.default`.
        verbose : bool
            If ``True``, the underlying estimator prints live progress.
        """

        if config is None:
            config = RMBConfig.default()

        while True:
            if config not in self._data:
                self._data[config] = self._default_bayesian_estimator()

            estimator = self._data[config]
            estimator.run(lambda: self._backend.fidelity_estimation(config, self.rng), verbose)

            # FIXME: add breaking conditions, e.g. times or number of samples
            if len(self._data) >= 8:
                break

            config = self._update_strategy(self._data, config)

        RMB.print_data(self._data)

    @classmethod
    def print_data(cls, data: RMBData):
        """Print ``data`` as a markdown-style table of configs and fidelities."""
        headers = ["depth", "two_qudit_gate_ratio", "n_qubits", "random_elimination", "fidelity"]
        rows = []
        for config, estimator in data.items():
            rows.append([
                f"{config.depth}",
                f"{config.two_qudit_gate_ratio}",
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

    def save(self, path: str | Path = "rmb.json"):
        """
        Save this RMB run as JSON.

        The output file contains the backend parameters, the Bayesian
        estimator parameters, and one record per ``(config, estimator)``
        pair in :attr:`_data`. Each record includes the config fields
        and the estimator's per-outcome counts so :meth:`load` can
        reconstruct the full state.

        Parameters
        ----------
        path : str | Path
            Output file. Bare file names (no separators) are resolved
            inside the package's ``data/`` directory; otherwise the
            given path is used as-is.
        """
        import json
        path = Path(path)
        if not path.is_absolute() and path.parent == Path("."):
            path = DATA_DIR / path
        path.parent.mkdir(parents=True, exist_ok=True)

        records = []
        for config, estimator in self._data.items():
            records.append({
                "depth": config.depth,
                "scrambling_probability": config.scrambling_probability,
                "two_qudit_gate_ratio": config.two_qudit_gate_ratio,
                "n_qubits": config.n_qubits,
                "random_elimination": config.random_elimination,
                "gates_set": [g.name for g in config.gates_set],
                "results": [[outcome, count] for outcome, count in estimator._results.items()],
            })

        payload = {
            "backend": self._backend.to_dict(),
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
        """
        Load an RMB run previously written by :meth:`save`.

        Reconstructs the backend and estimator parameters, then rebuilds
        each entry in :attr:`_data` by replaying the saved per-outcome
        counts into a fresh Bayesian estimator. ``gates_set`` is not
        reconstructed and falls back to the :class:`RMBConfig` default.

        Parameters
        ----------
        path : str | Path
            Input file. Bare file names (no separators) are resolved
            inside the package's ``data/`` directory; otherwise the
            given path is used as-is.
        rng : numpy.random.Generator | None
            RNG to attach to the rebuilt RMB. ``None`` uses ``default_rng()``.
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
        rmb = cls(update_strategy=RMBConfig.default_update_strategy,
                  backend=backend,
                  rng=rng,
                  threshold=estimator_params["threshold"],
                  min_runs=estimator_params["min_runs"],
                  max_runs=estimator_params["max_runs"])

        for rec in payload["data"]:
            config = RMBConfig(
                depth=rec["depth"],
                scrambling_probability=rec["scrambling_probability"],
                two_qudit_gate_ratio=rec["two_qudit_gate_ratio"],
                n_qubits=rec["n_qubits"],
                random_elimination=rec["random_elimination"],
            )
            estimator = rmb._default_bayesian_estimator()
            for outcome, count in rec["results"]:
                for _ in range(count):
                    estimator._record_result(outcome)
            rmb._data[config] = estimator
        return rmb

    @classmethod
    def plot_data(cls, data: RMBData, ax=None, show: bool = True):
        """
        Scatter ``data`` with depth on x, two-qudit gate ratio on y,
        and color encoding the fidelity (bright green = 1.0, dark
        red = 0.0).

        Each data point is drawn as three concentric disks whose colors
        encode ``fidelity - std``, ``fidelity``, and ``fidelity + std``
        from outer to inner (so the inner disk is greener and the outer
        disk is redder when the estimate is uncertain). ``std`` is the
        square root of the estimator variance.

        Parameters
        ----------
        data : RMBData
            Mapping from configurations to their Bayesian estimators.
        ax : matplotlib.axes.Axes | None
            Axes to plot on. If ``None``, a new figure and axes are
            created.
        show : bool
            If ``True``, call ``plt.show()`` after building the plot.

        Returns
        -------
        matplotlib.axes.Axes
            The axes the plot was drawn on.
        """
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap

        if ax is None:
            _, ax = plt.subplots()

        depths = np.array([c.depth for c in data])
        ratios = np.array([c.two_qudit_gate_ratio for c in data])
        fidelities = np.array([e.probability(True) for e in data.values()])
        stds = np.array([np.sqrt(e.variance(True)) for e in data.values()])

        cmap = LinearSegmentedColormap.from_list(
            "darkred_to_lime",
            ["darkred", "red", "orange", "lime", "green"])

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
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("fidelity")

        if show:
            plt.show()

        return ax


if __name__ == "__main__":
    rng = default_rng()
    verbose = True

    initial_config = RMBConfig.default()\
        .with_depth(20)\
        .with_random_elimination(0.25)\
        .with_n_qubits(4)\
        .with_two_qudit_gate_ratio(0.75)\
        .with_scrambling_probability(0.5)

    # noise_model = GenericNoise.from_paulis([0.001, 0.001, 0.002], rng=rng)
    # two_qubit_noise_model = GenericNoise.from_paulis([0.0075, 0.0075, 0.0075], rng=rng)
    # rmb = RMB.default(rng)\
    #     .with_backend(SympleqBackend(noise_model, two_qubit_noise_model))\
    #     .with_bayesian_estimator(threshold=10**(-4), min_runs=100)

    rmb = RMB.default(rng)\
        .with_backend(QuantinuumBackend(device_name="H2-2E"))\
        .with_bayesian_estimator(threshold=10**(-1), min_runs=10)

    rmb.run(initial_config, verbose)
    rmb.save("quantinuum.json")
    RMB.plot_data(rmb._data)
