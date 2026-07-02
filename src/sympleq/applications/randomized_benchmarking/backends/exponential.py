from __future__ import annotations

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementOutcomes,
    MeasurementRequest,
    RMBBackend,
    ShotRNG,
)
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.bayesian_estimation import BayesianEstimator


def asymptote_value(model: str | float, q: int | float) -> float:
    """Survival asymptote B(Q) for an RB exponential-link model."""
    if isinstance(model, str):
        if model == "depolarizing":
            return 2.0 ** (-int(q))
        if model == "zero":
            return 0.0
        raise ValueError(f"Unknown asymptote model: {model!r}")
    return float(model)


class ExponentialBackend(RMBBackend):
    """
    Direct Bernoulli simulator for the RB exponential link.

    This backend is a positive-control model for acquisition and fitting
    methods.  It does not build SympleQ circuits.  For each requested config it
    samples independent Bernoulli outcomes from

        p = V * (1 - B(Q)) * 2^{-n D(r,Q)} + B(Q),

    with

        D(r,Q) = ((1-r) L1(Q) + r L2(Q)) / log(2),
        Li(Q) = Li0 + mi (Q - Qref).

    The rates L1, L2 and slopes mi are natural-log decay rates per gate.  The
    fidelity-0.5 renormalized boundary is n*(r,Q) = 1 / D(r,Q).
    """

    type = "exponential"

    def __init__(
        self,
        one_qubit_rate: float,
        two_qubit_rate: float,
        *,
        one_qubit_q_slope: float = 0.0,
        two_qubit_q_slope: float = 0.0,
        q_reference: float = 0.0,
        visibility: float = 1.0,
        asymptote_model: str | float = "depolarizing",
    ) -> None:
        self.one_qubit_rate = float(one_qubit_rate)
        self.two_qubit_rate = float(two_qubit_rate)
        self.one_qubit_q_slope = float(one_qubit_q_slope)
        self.two_qubit_q_slope = float(two_qubit_q_slope)
        self.q_reference = float(q_reference)
        self.visibility = float(np.clip(visibility, 0.0, 1.0))
        self.asymptote_model = asymptote_model

    def _asymptote(self, q: int) -> float:
        return asymptote_value(self.asymptote_model, q)

    def _rate(self, config: RMBConfig) -> float:
        q = float(config.n_qubits)
        dq = q - self.q_reference
        l1 = self.one_qubit_rate + self.one_qubit_q_slope * dq
        l2 = self.two_qubit_rate + self.two_qubit_q_slope * dq
        ratio = float(config.ratio_2_qb_gates)
        natural_rate = (1.0 - ratio) * l1 + ratio * l2
        return float(max(natural_rate / np.log(2.0), 0.0))

    def probability(self, config: RMBConfig) -> float:
        b = float(np.clip(self._asymptote(config.n_qubits), 0.0, 1.0))
        d = self._rate(config)
        decay = 2.0 ** (-float(config.n_gates) * d)
        p = self.visibility * (1.0 - b) * decay + b
        return float(np.clip(p, 0.0, 1.0))

    def fidelity_estimation(
        self,
        requests: list[MeasurementRequest],
        rng: RNGGenerator,
        shot_rng: ShotRNG | None = None,
    ) -> MeasurementOutcomes:
        outcomes: dict[RMBConfig, list[bool]] = {}
        shot_counts: dict[RMBConfig, int] = {}
        for request in requests:
            p = self.probability(request.config)
            results = outcomes.setdefault(request.config, [])
            for _ in range(max(0, request.shots)):
                index = shot_counts.get(request.config, 0)
                shot_counts[request.config] = index + 1
                sample_rng = rng if shot_rng is None else shot_rng(request.config, index)
                results.append(bool(sample_rng.random() < p))
        return MeasurementOutcomes(outcomes=outcomes)

    @classmethod
    def default_config(cls) -> RMBConfig:
        return RMBConfig.default()

    def default_estimator(self) -> BayesianEstimator:
        return BayesianEstimator.default()

    def to_dict(self) -> dict:
        return {
            "type": self.type,
            "one_qubit_rate": self.one_qubit_rate,
            "two_qubit_rate": self.two_qubit_rate,
            "one_qubit_q_slope": self.one_qubit_q_slope,
            "two_qubit_q_slope": self.two_qubit_q_slope,
            "q_reference": self.q_reference,
            "visibility": self.visibility,
            "asymptote_model": self.asymptote_model,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> ExponentialBackend:
        return cls(
            one_qubit_rate=payload["one_qubit_rate"],
            two_qubit_rate=payload["two_qubit_rate"],
            one_qubit_q_slope=payload.get("one_qubit_q_slope", 0.0),
            two_qubit_q_slope=payload.get("two_qubit_q_slope", 0.0),
            q_reference=payload.get("q_reference", 0.0),
            visibility=payload.get("visibility", 1.0),
            asymptote_model=payload.get("asymptote_model", "depolarizing"),
        )
