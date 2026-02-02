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

For uncorrelated noise on multiple qudits, the multi-qudit quantities are
tensor products of single-qudit quantities.

For now, we assume that Kraus operators are Clifford. This simplifies greatly how they act on PauliSums,
as they are basically Gates.
"""

from __future__ import annotations
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from abc import ABC, abstractmethod
from itertools import product
import math

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
        # Initialize noise models probabilities.
        # In principle, we don't know if there are gates with n_qudits larger than 2,
        # good enough for now.
        self.cached_probabilities = {}

        for gate_n_qudits in (1, 2, 3):
            probs = self.kraus_probabilities(gate_n_qudits)
            # Given probabilities [p0, p1, p2, p3], cumsum gives [p0, p0+p1, p0+p1+p2, 1.0].
            # This allows O(log n) sampling via searchsorted with a uniform random number
            # in _apply_gate_to_pauli_with_error.
            self.cached_probabilities[gate_n_qudits] = np.cumsum(probs)

        self.cached_operators = {}

    @abstractmethod
    def n_kraus_operators(self) -> int:
        pass

    @classmethod
    @abstractmethod
    def kraus_operators(cls,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray,
                        weighted: bool = True) -> list[PauliSum]:
        """"
        K_i = sum_{j} alpha_{ij} sigma_j
        """
        pass

    @abstractmethod
    def kraus_probabilities(self, n_qudits: int) -> np.ndarray:
        """
        Returns the Kraus probabilities for the noise channel.

        p_i = sum_j |α_{ij}|², i.e. the squared norm of each Kraus operator's coefficients.
        Output is an array of length n^M where n = n_kraus_operators and M = n_qudits.
        The probabilities sum to 1.
        """
        return np.sum(self.process_matrix(n_qudits), axis=1)

    @abstractmethod
    def process_matrix(self, n_qudits: int) -> np.ndarray:
        """"
        Returns the full process matrix λ_{ij} of shape (n^M, n^M),
        where n = n_kraus_operators and M = n_qudits.

        lambda_{ij} = sum_k alpha_{ki} alpha^*_{kj} in Eq.(3), where i, j are multi-indices.
        """
        pass

    def act(self, pauli_sum: PauliSum, qudit_indices: np.ndarray) -> PauliSum:
        n_qudits = len(qudit_indices)
        # Get probabilities to select one possible quantum trajectory
        if n_qudits not in self.cached_probabilities:
            probs = self.kraus_probabilities(n_qudits)
            # Given probabilities [p0, p1, p2, p3], cumsum gives [p0, p0+p1, p0+p1+p2, 1.0].
            # This allows O(log n) sampling via searchsorted with a uniform random number
            # in _apply_gate_to_pauli_with_error.
            self.cached_probabilities[n_qudits] = np.cumsum(probs)

        probs = self.cached_probabilities[n_qudits]
        idx = int(np.searchsorted(probs, self.rng.random()))

        key = (pauli_sum.dimensions.tobytes(), qudit_indices.tobytes())
        if key not in self.cached_operators:
            self.cached_operators[key] = self.kraus_operators(
                pauli_sum.dimensions, qudit_indices, weighted=False)

        operators = self.cached_operators[key]
        k = operators[idx]
        return (k * pauli_sum * k.H())


class Noiseless(NoiseModel):
    def __init__(self) -> None:
        super().__init__(rng=default_rng())

    def n_kraus_operators(self) -> int:
        return 1

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray,
                        weighted: bool = True) -> list[PauliSum]:

        # The noiseless channel has a single Kraus operator (the identity), with probability 1.
        tableau = np.zeros(2 * len(dimensions), dtype=int)
        sigma_0 = PauliSum.from_tableau(tableau, dimensions)
        return [sigma_0]

    def process_matrix(self, n_qudits: int) -> np.ndarray:
        """
        Returns the process matrix for the noiseless channel.
        Always a 1x1 matrix with value 1.0, since n_kraus_operators=1
        and 1^M = 1 for any number of qudits M.
        """
        return np.array([[1.0]], dtype=float)


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
        self.alpha_00 = math.sqrt(self.p0)
        self.alpha_11 = math.sqrt(1.0 - self.p0)

        if rng is None:
            rng = default_rng()
        super().__init__(rng)

    def n_kraus_operators(self) -> int:
        return 2

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray,
                        weighted: bool = True) -> list[PauliSum]:

        n_qudits = len(dimensions)
        n_gate_qudits = len(qudit_indices)

        # Get coefficients from probabilities (sqrt of probabilities)
        if weighted:
            probabilities = self.kraus_probabilities(n_gate_qudits)
            coefficients = np.sqrt(probabilities)
        else:
            coefficients = [1 for _ in range(self.n_kraus_operators()**n_qudits)]

        tableau = np.zeros(2 * n_qudits, dtype=int)
        output = []

        # Build all combinations of the Kraus operators for each qudit in qudit_indices
        combinations = list(product(range(self.n_kraus_operators()), repeat=n_gate_qudits))
        for i, comb in enumerate(combinations):
            p_tableau = tableau.copy()
            for idx, exp in zip(qudit_indices, comb):
                # identity
                if exp == 0:
                    pass
                # sigma-z
                elif exp == 1:
                    p_tableau[n_qudits + idx] = 1

            sigma = PauliSum.from_tableau(p_tableau, dimensions)
            output.append(coefficients[i] * sigma)

        return output

    def process_matrix(self, n_qudits: int) -> np.ndarray:
        """
        Returns the process matrix for the dephasing channel.
        Shape is (2^M, 2^M) where M = n_qudits, since n_kraus_operators=2.
        For uncorrelated noise, the multi-qudit matrix is the tensor product
        of single-qudit process matrices.
        """
        single_qudit = np.array([
            [np.abs(self.alpha_00)**2, 0.0],
            [0.0, np.abs(self.alpha_11)**2]
        ], dtype=float)

        result = single_qudit
        for _ in range(n_qudits - 1):
            result = np.kron(result, single_qudit)

        return result


class DepolarizingNoise(NoiseModel):
    def __init__(self, error_rate: float, rng: RNGGenerator | None = None) -> None:
        # The depolarizing channel has tow Klaus operators: the identity and Z.
        # A single parameter p0 models the noise probability,
        # representing the probability of having no error.
        # Kraus operators:
        # K0 = sqrt(p0) 1
        # Ki =  sqrt(1 − p0/3) σi
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        self.alpha_00 = math.sqrt(self.p0)
        self.alpha_ii = math.sqrt((1.0 - self.p0) / 3)

        if rng is None:
            rng = default_rng()
        super().__init__(rng)

    def n_kraus_operators(self) -> int:
        return 4

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray,
                        weighted: bool = True) -> list[PauliSum]:

        n_qudits = len(dimensions)
        n_gate_qudits = len(qudit_indices)

        # Get coefficients from probabilities (sqrt of probabilities)
        if weighted:
            probabilities = self.kraus_probabilities(n_gate_qudits)
            coefficients = np.sqrt(probabilities)
        else:
            coefficients = [1 for _ in range(self.n_kraus_operators()**n_qudits)]

        tableau = np.zeros(2 * n_qudits, dtype=int)
        output = []

        combinations = list(product(range(self.n_kraus_operators()), repeat=n_gate_qudits))
        for i, comb in enumerate(combinations):
            p_tableau = tableau.copy()
            phases = np.zeros(1, dtype=int)
            for idx, exp in zip(qudit_indices, comb):
                # identity
                if exp == 0:
                    pass
                # sigma-x
                elif exp == 1:
                    p_tableau[idx] = 1
                # sigma-y
                elif exp == 2:
                    p_tableau[idx] = 1
                    p_tableau[n_qudits + idx] = 1
                    phases[0] = 1  # FIXME: generalize to qudit?
                # sigma-z
                elif exp == 3:
                    p_tableau[n_qudits + idx] = 1

            sigma = PauliSum.from_tableau(p_tableau, dimensions)
            output.append(coefficients[i] * sigma)

        return output

    def process_matrix(self, n_qudits: int) -> np.ndarray:
        """
        Returns the process matrix for the depolarizing channel.
        Shape is (4^M, 4^M) where M = n_qudits, since n_kraus_operators=4.
        For uncorrelated noise, the multi-qudit matrix is the tensor product
        of single-qudit process matrices.
        """
        single_qudit = np.diag([
            np.abs(self.alpha_00)**2,
            np.abs(self.alpha_ii)**2,
            np.abs(self.alpha_ii)**2,
            np.abs(self.alpha_ii)**2
        ]).astype(float)

        result = single_qudit
        for _ in range(n_qudits - 1):
            result = np.kron(result, single_qudit)

        return result


class CompositeNoise(NoiseModel):
    """
    A noise model that combines multiple noise models by merging their Kraus operators.

    Operators representing the same Pauli are combined by summing their weights
    and normalizing by the total number of models.
    """

    def __init__(self, noise_models: list[NoiseModel], rng: RNGGenerator) -> None:
        self.noise_models = noise_models
        self._n_models = len(noise_models)
        super().__init__(rng)

    @classmethod
    def from_noise_models(cls, noise_models: list[NoiseModel],
                          rng: RNGGenerator | None = None) -> CompositeNoise:
        if not noise_models:
            raise ValueError("noise_models list cannot be empty")
        if rng is None:
            rng = default_rng()
        return cls(noise_models, rng)

    @classmethod
    def from_process_matrix(cls, noise_models: list[NoiseModel],
                            rng: RNGGenerator | None = None) -> CompositeNoise:
        if not noise_models:
            raise ValueError("noise_models list cannot be empty")
        if rng is None:
            rng = default_rng()
        return cls(noise_models, rng)

    def n_kraus_operators(self) -> int:
        # Return the number of unique Pauli operators across all models
        # Use dummy dimensions to count unique operators
        dimensions = [2]
        qudit_indices = [0]
        combined = self._combine_operators(dimensions, qudit_indices)
        return len(combined)

    def _combine_operators(self,
                           dimensions: list[int] | np.ndarray,
                           qudit_indices: list[int] | np.ndarray) -> dict[PauliSum, complex]:
        """
        Collect and combine Kraus operators from all underlying models.

        Returns a dict mapping tableau bytes to combined PauliSum operators.
        Operators with the same Pauli (tableau) have their weights summed and normalized.
        """
        # Track weights and base operators separately
        operator_weights: dict[PauliSum, complex] = {}

        n_qudits = len(dimensions)

        for model in self.noise_models:
            operators = model.kraus_operators(dimensions, qudit_indices, weighted=False)
            probs = model.kraus_probabilities(n_qudits)
            weights = np.sqrt(probs)
            for idx, op in enumerate(operators):
                if op.is_identity:
                    continue
                weight = weights[idx]

                # If operator is unity, don't add to the weight.
                # The noiseless channel will have a weight equal to 1 - Sum toher_channels

                if op in operator_weights:
                    operator_weights[op] += weight
                else:
                    operator_weights[op] = weight

        # Normalize by number of models
        for op in operator_weights:
            operator_weights[op] /= self._n_models

        # Normalize weights so probabilities (|w|^2) sum to 1
        norm = np.sqrt(sum(np.abs(w)**2 for w in operator_weights.values()))
        for op in operator_weights:
            operator_weights[op] /= norm

        return operator_weights

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray,
                        weighted: bool = True) -> list[PauliSum]:

        combined = self._combine_operators(dimensions, qudit_indices)
        operators = []
        for op, weight in combined.items():
            if weighted:
                operators.append(weight * op)
            else:
                operators.append(op)

        return operators

    def process_matrix(self, n_qudits: int) -> np.ndarray:
        """
        Returns the combined process matrix.

        Averages the process matrices from all underlying models.
        """
        result = None
        for model in self.noise_models:
            matrix = model.process_matrix(n_qudits)
            if result is None:
                result = matrix
            else:
                # Matrices may have different sizes, need to handle that
                if result.shape == matrix.shape:
                    result = result + matrix
                else:
                    # Pad smaller matrix to match larger
                    max_size = max(result.shape[0], matrix.shape[0])
                    if result.shape[0] < max_size:
                        new_result = np.zeros((max_size, max_size), dtype=float)
                        new_result[:result.shape[0], :result.shape[1]] = result
                        result = new_result
                    if matrix.shape[0] < max_size:
                        new_matrix = np.zeros((max_size, max_size), dtype=float)
                        new_matrix[:matrix.shape[0], :matrix.shape[1]] = matrix
                        matrix = new_matrix
                    result = result + matrix

        assert result is not None

        return result / self._n_models
