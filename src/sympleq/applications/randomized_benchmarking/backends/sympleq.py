from __future__ import annotations
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.noise.noise_model import NoiseModel


class SympleqBackend(RMBBackend):
    """
    Native SympleQ simulation backend.

    Builds a random circuit from the config, applies the configured
    single- and two-qudit noise models, and returns whether the final
    Pauli-sum state matches the initial state.

    Parameters
    ----------
    noise_model : NoiseModel | None
        Noise model applied to all gates after circuit construction.
    two_qubit_noise_model : NoiseModel | None
        Noise model applied only to two-qudit gates after circuit
        construction.
    """

    type = "sympleq"

    def __init__(self,
                 noise_model: NoiseModel | None = None,
                 two_qubit_noise_model: NoiseModel | None = None) -> None:
        self.noise_model = noise_model
        self.two_qubit_noise_model = two_qubit_noise_model

    def fidelity_estimation(self, config: RMBConfig, rng: RNGGenerator) -> list[bool]:
        initial_state = config.initial_state()
        n_runs = 1
        results = []
        for _ in range(n_runs):
            circuit = config.random_circuit(rng=rng)
            if self.noise_model is not None:
                circuit = circuit.with_noise(self.noise_model)
            if self.two_qubit_noise_model is not None:
                circuit = circuit.with_two_qudit_noise(self.two_qubit_noise_model)
            final_state = circuit.act(initial_state)
            results.append(final_state == initial_state)

        return results

    def to_dict(self) -> dict:
        return {
            "type": self.type,
            "noise_model": self.noise_model.to_dict() if self.noise_model is not None else None,
            "two_qubit_noise_model": (self.two_qubit_noise_model.to_dict()
                                      if self.two_qubit_noise_model is not None else None),
        }

    @classmethod
    def from_dict(cls, payload: dict) -> SympleqBackend:
        noise_model = (NoiseModel.from_dict(payload["noise_model"])
                       if payload.get("noise_model") is not None else None)
        two_qubit_noise_model = (NoiseModel.from_dict(payload["two_qubit_noise_model"])
                                 if payload.get("two_qubit_noise_model") is not None else None)
        return cls(noise_model=noise_model, two_qubit_noise_model=two_qubit_noise_model)
