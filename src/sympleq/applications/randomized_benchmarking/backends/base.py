from __future__ import annotations
from abc import ABC, abstractmethod
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.bayesian_estimation import BayesianEstimator


class RMBBackend(ABC):
    """
    Abstract backend for randomized benchmarking.

    A backend pairs a fidelity-estimation procedure with the parameters
    it depends on (e.g. noise models for simulation, device name for
    hardware). It is also responsible for serializing those parameters
    via :meth:`to_dict` so that an :class:`RMB` run can be saved and
    fully reconstructed.
    """

    type: str

    @abstractmethod
    def fidelity_estimation(self, config: RMBConfig, rng: RNGGenerator) -> list[bool]:
        """Run a single fidelity-estimation trial for ``config``."""
        ...

    @classmethod
    @abstractmethod
    def default_config(cls) -> RMBConfig:
        """Returns the default RMB config for the backend."""
        ...

    @abstractmethod
    def default_estimator(self) -> BayesianEstimator:
        """Returns the default BayesianEstimator  for the backend."""
        ...

    @abstractmethod
    def to_dict(self) -> dict:
        """Serialize backend parameters to a JSON-friendly dict."""
        ...

    @classmethod
    @abstractmethod
    def from_dict(cls, payload: dict) -> RMBBackend:
        """Rebuild a backend from a dict produced by :meth:`to_dict`."""
        ...


def backend_from_dict(payload: dict) -> RMBBackend:
    """
    Construct an :class:`RMBBackend` from a serialized payload.

    Dispatches on ``payload["type"]`` to the matching backend subclass.
    """
    from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
    from sympleq.applications.randomized_benchmarking.backends.quantinuum import QuantinuumBackend

    type_name = payload["type"]
    if type_name == SympleqBackend.type:
        return SympleqBackend.from_dict(payload)
    if type_name == QuantinuumBackend.type:
        return QuantinuumBackend.from_dict(payload)
    raise ValueError(f"Unknown backend type: {type_name}")
