from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Callable
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.bayesian_estimation import BayesianEstimator


@dataclass(frozen=True)
class MeasurementRequest:
    """``shots`` independently drawn random circuits of one config."""
    config: RMBConfig
    shots: int


@dataclass(frozen=True)
class MeasurementOutcomes:
    """
    Boolean fidelity outcomes of one :meth:`RMBBackend.fidelity_estimation` call.

    ``outcomes`` lists each config's results in execution order. ``cost`` is
    what the backend actually charged for the call, in its native cost units
    (HQC for Quantinuum; local simulation is free). ``n_submissions`` counts
    the device submissions the backend packed the circuits into.
    """
    outcomes: dict[RMBConfig, list[bool]] = field(default_factory=dict)
    cost: float = 0.0
    n_submissions: int = 0


type ShotRNG = Callable[[RMBConfig, int], RNGGenerator]
"""Maps (config, per-call shot index) to the rng driving that shot."""


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
    def fidelity_estimation(self, requests: list[MeasurementRequest], rng: RNGGenerator,
                            shot_rng: ShotRNG | None = None) -> MeasurementOutcomes:
        """
        Run the requested circuits and return their Boolean outcomes.

        Each requested shot is one independently drawn random circuit whose
        outcome is ``True`` when the measured state matches the initial
        state. With ``shot_rng``, shot ``i`` of a config (counting across
        the call's requests) draws from ``shot_rng(config, i)``, so seeded
        callers record the same outcomes no matter how requests are grouped
        into calls; ``None`` draws everything from ``rng``. How the circuits
        are executed (e.g. stitched into device submissions) is an
        implementation detail of the backend.
        """
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
