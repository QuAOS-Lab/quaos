"""
Noise models for quantum channels, based on the Kraus operator formalism.

Reference: "Enhancing Quantum Computation via Superposition of Quantum Gates",
Phys. Rev. A 108, 062604 (2023). DOI: 10.1103/PhysRevA.108.062604

Overview
--------
A quantum noise channel transforms a Pauli string PS as:

    PS_out = Σ_{m,n} λ_{mn} σ_m PS σ†_n

where σ_m are Pauli basis operators and λ is the process matrix.

This can equivalently be written using Kraus operators K_i:

    PS_out = Σ_i K_i PS K†_i

where each Kraus operator is a linear combination of Paulis:

    K_i = Σ_j α_{ij} σ_j

The process matrix λ and Kraus coefficients α are related by:

    λ_{mn} = Σ_i α_{im} α*_{in}

The probability of each Kraus operator (used for quantum trajectory sampling) is:

    p_i = Σ_j |α_{ij}|² = ||K_i||²

We take two simplifying assumptions:
**Uncorrelated noise**
For uncorrelated noise on multiple qudits, the multi-qudit quantities are
tensor products of single-qudit quantities.

**Clifford noise**
For now, we assume that Kraus operators are Clifford. This simplifies greatly how they act on PauliSums,
as they are basically Gates.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
import itertools
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.core.circuits.gates import GATES, Gate
from sympleq.core.circuits.utils import embed_unitary
from sympleq.core.paulis._typing import DimensionsType, HilbertOperator
from sympleq.core.paulis.pauli_object import PauliObject
from sympleq.core.paulis.pauli_sum import PauliSum


class NoiseModel(ABC):
    """
    Abstract base class for quantum noise models.

    Subclasses must implement:
    - n_kraus_operators(): number of single-qudit Kraus operators
    - kraus_probabilities(n_qudits): probabilities p_i = ||K_i||²
    - kraus_operators(dimensions, qudit_indices): list of Kraus operators as PauliSums
    - process_matrix(n_qudits): the λ_{mn} matrix
    """
    def __init__(self, rng: RNGGenerator) -> None:
        self.rng = rng

    @abstractmethod
    def n_qudits(self) -> int:
        pass

    @abstractmethod
    def kraus_gates(self) -> list[Gate]:
        """"
        The gate corresponding to the Kraus operator, acting as ps1 = G_i.act(ps, qudit_indices)
        """
        pass

    @abstractmethod
    def kraus_probabilities(self) -> list[float]:
        """
        Returns the Kraus probabilities for the noise channel.

        p_i = sum_j |α_{ij}|², i.e. the squared norm of each Kraus operator's coefficients.
        The probabilities sum to 1.
        """
        pass

    def apply_quantum_trajectory(self, pauli_sum: PauliObject, qudits: tuple[int, ...]) -> PauliObject:
        # NOTE: this default implementation is valid for uncorrelated noise acting on single qudits only
        n_qudits = len(qudits)
        kraus_n_qudits = self.n_qudits()
        if n_qudits < kraus_n_qudits:
            return pauli_sum
        # Get probabilities to select one possible quantum trajectory
        # Given probabilities [p0, p1, p2, p3], cumsum gives [p0, p0+p1, p0+p1+p2, 1.0].
        # This allows O(log n) sampling via searchsorted with a uniform random number
        # in _apply_gate_to_pauli_with_error.
        probs = np.cumsum(self.kraus_probabilities())
        idx = int(np.searchsorted(probs, self.rng.random()))
        kraus_gate = self.kraus_gates()[idx]
        if kraus_gate.n_qudits == n_qudits:
            pauli_sum = kraus_gate.act(pauli_sum, qudits)
        elif kraus_gate.n_qudits == 1:
            # Act independently on each qudit
            for qudit in qudits:
                pauli_sum = kraus_gate.act(pauli_sum, qudit)
        # FIXME: handle generic case, else:...

        return pauli_sum

    def act_in_hilbert_space(self, rho: HilbertOperator,
                             qudits: tuple[int, ...], dimensions: DimensionsType) -> HilbertOperator:

        def _apply_unitary(output_rho: HilbertOperator | None, unitary: HilbertOperator,
                           probability: float) -> HilbertOperator:
            if output_rho is None:
                output_rho = probability * (unitary @ rho @ unitary.conjugate().transpose())
            else:
                output_rho += probability * (unitary @ rho @ unitary.conjugate().transpose())
            assert output_rho is not None
            return output_rho

        dimension = dimensions[qudits[0]]
        n_qudits = len(qudits)

        kraus_n_qudits = self.n_qudits()
        if n_qudits < kraus_n_qudits:
            return rho

        output_rho: HilbertOperator | None = None

        # E.g.: single qudit noise on single qudit gate
        cum_prob = 0
        if n_qudits == kraus_n_qudits:
            for kraus_gate, probability in zip(self.kraus_gates(), self.kraus_probabilities()):
                unitary = embed_unitary(kraus_gate.local_unitary(dimension), qudits, dimensions)
                output_rho = _apply_unitary(output_rho, unitary, probability)
                cum_prob += probability

            assert np.abs(cum_prob - 1.0) < 10**(-5), f"cum prob {cum_prob}"
            assert output_rho is not None
            return output_rho

        if kraus_n_qudits == 1:
            # Act independently on each qudit.

            # Get all gates combinations
            all_gates_combinations = list(itertools.product(self.kraus_gates(), repeat=n_qudits))
            # Pray the order is correct, like, do it
            all_probabilities_combinations = list(itertools.product(self.kraus_probabilities(), repeat=n_qudits))
            for kraus_gates, probabilities in zip(all_gates_combinations, all_probabilities_combinations):
                unitary = embed_unitary(kraus_gates[0].local_unitary(dimension), (qudits[0],), dimensions)
                for kraus_gate, qudit in zip(kraus_gates[1:], qudits[1:]):
                    unitary = embed_unitary(kraus_gate.local_unitary(dimension), (qudit,), dimensions) @ unitary

                probability = float(np.prod(probabilities))
                output_rho = _apply_unitary(output_rho, unitary, probability)
                cum_prob += probability

            assert np.abs(cum_prob - 1.0) < 10**(-5), f"cum prob {cum_prob}"
            assert output_rho is not None
            return output_rho

        # FIXME: extend this for self.n_qudits > 1

        return rho


class Noiseless(NoiseModel):
    def __init__(self) -> None:
        super().__init__(rng=default_rng())
        self._probabilities = [1.0]

    def n_qudits(self) -> int:
        return 1

    def kraus_gates(self) -> list[Gate]:
        return [GATES.Id]

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities

    def apply_quantum_trajectory(self, pauli_sum: PauliSum, qudits: tuple[int, ...]) -> PauliSum:
        return pauli_sum

    def act_in_hilbert_space(self, rho: HilbertOperator,
                             qudits: tuple[int, ...], dimensions: DimensionsType) -> HilbertOperator:
        return rho


class DephasingNoise(NoiseModel):
    def __init__(self, error_rate: float, rng: RNGGenerator | None = None) -> None:
        # The dephasing channel has tow Klaus operators: the identity and Z.
        # A single parameter p0 models the noise probability,
        # representing the probability of having no error.
        # Kraus operators:
        # K0 = sqrt(p0) 1
        # K1 = sqrt(1 − p0) Z
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        self._probabilities = [self.p0, 1.0 - self.p0]

        if rng is None:
            rng = default_rng()
        super().__init__(rng)

    def n_qudits(self) -> int:
        return 1

    def kraus_gates(self) -> list[Gate]:
        return [GATES.Id, GATES.Z]

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities


class DepolarizingNoise(NoiseModel):
    def __init__(self, error_rate: float, rng: RNGGenerator | None = None) -> None:
        # The depolarizing channel has four Klaus operators: the 4 paulis.
        # A single parameter p0 models the noise probability,
        # representing the probability of having no error.
        # Kraus operators:
        # K0 = sqrt(p0) 1
        # Ki =  sqrt(1 − p0/3) σi
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        self._probabilities = [self.p0, (1.0 - self.p0) / 3, (1.0 - self.p0) / 3, (1.0 - self.p0) / 3]

        if rng is None:
            rng = default_rng()
        super().__init__(rng)

    def n_qudits(self) -> int:
        return 1

    def kraus_gates(self) -> list[Gate]:
        return [GATES.Id, GATES.X, GATES.Y, GATES.Z]

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities


class CompositeNoise(NoiseModel):
    """
    A noise model that combines multiple noise models by merging their Kraus operators.

    Operators representing the same Pauli are combined by summing their weights
    and normalizing by the total number of models.
    """

    def __init__(self, noise_models: list[NoiseModel], rng: RNGGenerator) -> None:
        self.noise_models = noise_models
        self._n_models = len(noise_models)
        combined_gates: dict[Gate, float] = self._combine_gates()
        self._probabilities = list(combined_gates.values())
        self._gates = list(combined_gates.keys())

        super().__init__(rng)

    @classmethod
    def from_noise_models(cls, noise_models: list[NoiseModel],
                          rng: RNGGenerator | None = None) -> CompositeNoise:
        if not noise_models:
            raise ValueError("noise_models list cannot be empty")
        if rng is None:
            rng = default_rng()
        return cls(noise_models, rng)

    def _combine_gates(self, compound_identity: bool = True) -> dict[Gate, float]:
        """
        Collect and combine Kraus operators from all underlying models.

        Returns a dict mapping gates to combined probabilities:
        Operators with the same gates have their probability summed and normalized.
        """
        gates_probabilities: dict[Gate, float] = {}

        for model in self.noise_models:
            for gate, probability in zip(model.kraus_gates(), model.kraus_probabilities()):
                if not compound_identity and gate == GATES.Id:
                    continue

                # If operator is unity, don't add to the weight.
                # The noiseless channel will have a weight equal to 1 - Sum toher_channels

                if gate in gates_probabilities:
                    gates_probabilities[gate] += probability
                else:
                    gates_probabilities[gate] = probability

        # Normalize weights so probabilities (|w|^2) sum to 1
        norm = sum(gates_probabilities.values())
        for gate in gates_probabilities:
            gates_probabilities[gate] /= norm

        return gates_probabilities

    def n_qudits(self) -> int:
        return self.noise_models[0].n_qudits()

    def kraus_gates(self) -> list[Gate]:
        return self._gates

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities
