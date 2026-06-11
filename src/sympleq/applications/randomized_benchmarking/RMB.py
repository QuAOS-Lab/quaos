from __future__ import annotations
import copy
import json
from pathlib import Path
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.backends import RMBBackend, SympleqBackend
from sympleq.applications.randomized_benchmarking.backends.base import backend_from_dict
from sympleq.core.bayesian_estimation import BayesianEstimator

DATA_DIR = Path(__file__).parent / "rmb_data"


def resolve_data_path(path: str | Path) -> Path:
    """Resolve bare file names (no separators) inside the package's ``rmb_data`` directory."""
    path = Path(path)
    if not path.is_absolute() and path.parent == Path("."):
        path = DATA_DIR / path
    return path


class RMB:
    """
    Container for one randomized-benchmarking run.

    Pairs the backend that produces fidelity outcomes with the recorded
    data (config -> estimator) and its (de)serialization.
    """

    def __init__(self,
                 backend: RMBBackend,
                 rng: RNGGenerator,
                 default_estimator: BayesianEstimator | None = None) -> None:
        self.backend = backend
        self._default_estimator = default_estimator

        self.rng = rng
        self._data: RMBData = {}

    def new_default_estimator(self) -> BayesianEstimator:
        if self._default_estimator is None:
            return self.backend.default_estimator()
        return copy.deepcopy(self._default_estimator)

    def with_backend(self, backend: RMBBackend) -> RMB:
        """Set ``self.backend`` and return ``self`` for chaining."""
        self.backend = backend
        return self

    @classmethod
    def default(cls, rng: RNGGenerator | None = None) -> RMB:
        """Return an RMB with a noiseless ``SympleqBackend``.

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
        return cls(backend=SympleqBackend(), rng=rng)

    def save(self, path: str | Path):
        """Save this RMB run as JSON.

        Parameters
        ----------
        path : str | Path
            Output file. Bare file names (no separators) are resolved
            inside the package's ``rmb_data/`` directory; otherwise the given
            path is used as-is.
        """
        path = resolve_data_path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        records = []
        for config, estimator in self._data.items():
            records.append({
                "n_1qb_gates": config.n_1qb_gates,
                "n_2qb_gates": config.n_2qb_gates,
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
        path = resolve_data_path(path)
        with open(path) as f:
            payload = json.load(f)

        if rng is None:
            rng = default_rng()

        backend = backend_from_dict(payload["backend"])
        estimator_params = payload["estimator"]
        estimator = BayesianEstimator(estimator_params["threshold"],
                                      estimator_params["min_runs"], estimator_params["max_runs"])
        rmb = cls(backend=backend, rng=rng, default_estimator=estimator)

        for rec in payload["data"]:
            config = RMBConfig(
                n_1qb_gates=rec["n_1qb_gates"],
                n_2qb_gates=rec["n_2qb_gates"],
                n_qubits=rec["n_qubits"],
                random_elimination=rec["random_elimination"],
            )
            estimator = rmb.backend.default_estimator()
            for outcome, count in rec["results"]:
                estimator.record(outcome, count)
            rmb._data[config] = estimator
        return rmb
