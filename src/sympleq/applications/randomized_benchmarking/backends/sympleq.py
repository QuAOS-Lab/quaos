from __future__ import annotations
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.backends.utils import data_from_pytket_circuit_results, \
    load_pytket_circuits
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.core.noise.noise_model import NoiseModel
from sympleq.integrations.quantinuum.utils import from_pytket_circuit


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

    @classmethod
    def default_config(cls) -> RMBConfig:
        return RMBConfig.default()\
            .with_depth(20)\
            .with_random_elimination(0.25)\
            .with_n_qubits(6)\
            .with_two_qubit_gate_ratio(0.5)\
            .with_scrambling_probability(0.5)

    def default_estimator(self) -> BayesianEstimator:
        """Return a fresh :class:`BayesianEstimator` with ``threshold=0.3`` and ``min_runs=100``."""
        return BayesianEstimator(threshold=10**(-3), min_runs=100)

    def simulate_pytket_circuits(self, folder_name: str) -> RMBData:
        """Build :type:`RMBData` from the last ``n`` execute jobs in ``self.project_name``.

        Each fetched circuit's ``n_qubits``, non-measurement gate count, and
        two-qudit gate ratio are used to synthesise an :class:`RMBConfig`;
        remaining fields fall back to defaults. The outcome recorded is
        ``True`` iff the most-common measured bitstring is all-zeros.

        Parameters
        ----------
        n : int
            Maximum number of execute jobs to fetch.

        Returns
        -------
        RMBData
            Mapping of inferred configs to their estimators.
        """
        circuits = load_pytket_circuits(folder_name)
        print(f"Simulating {len(circuits)} circuits from {folder_name}")

        pairs = []
        for idx, p_circuit in enumerate(circuits):
            print(f"{idx / len(circuits) * 100:.2f}%", end="\r")
            circuit = from_pytket_circuit(p_circuit)
            initial_state = RMBConfig.initial_state_for_n_qubits(circuit.n_qudits())
            if self.noise_model is not None:
                circuit = circuit.with_noise(self.noise_model)
            if self.two_qubit_noise_model is not None:
                circuit = circuit.with_two_qudit_noise(self.two_qubit_noise_model)

            for _ in range(11):
                final_state = circuit.act(initial_state)
                pairs.append((p_circuit, final_state == initial_state))

        return data_from_pytket_circuit_results(pairs, self.default_estimator)

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
