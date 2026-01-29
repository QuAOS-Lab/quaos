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
"""

import numpy as np
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
    def __init__(self) -> None:
        pass

    @classmethod
    @abstractmethod
    def n_kraus_operators(cls) -> int:
        pass

    @classmethod
    @abstractmethod
    def kraus_operators(cls,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:
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
        pass

    @abstractmethod
    def process_matrix(self, n_qudits: int) -> np.ndarray:
        """"
        Returns the full process matrix λ_{ij} of shape (n^M, n^M),
        where n = n_kraus_operators and M = n_qudits.

        lambda_{ij} = sum_k alpha_{ki} alpha^*_{kj} in Eq.(3), where i, j are multi-indices.
        """
        pass


class Noiseless(NoiseModel):
    @classmethod
    def n_kraus_operators(cls) -> int:
        return 1

    def kraus_probabilities(self, n_qudits: int) -> np.ndarray:
        """
        Returns the Kraus probabilities for the noiseless channel.
        Always [1.0] since there's only one Kraus operator (identity).
        """
        return np.array([1.0], dtype=float)

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:

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
    def __init__(self, error_rate: float) -> None:
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
        super().__init__()

    @classmethod
    def n_kraus_operators(cls) -> int:
        return 2

    def kraus_probabilities(self, n_qudits: int) -> np.ndarray:
        """
        Returns the Kraus probabilities for the dephasing channel.
        Each probability is the product of single-qudit probabilities.
        """
        single_qudit_probs = np.array([
            np.abs(self.alpha_00)**2,
            np.abs(self.alpha_11)**2
        ], dtype=float)

        result = single_qudit_probs
        for _ in range(n_qudits - 1):
            result = np.kron(result, single_qudit_probs)

        return result

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:

        n_qudits = len(dimensions)
        n_gate_qudits = len(qudit_indices)

        # Get coefficients from probabilities (sqrt of probabilities)
        probabilities = self.kraus_probabilities(n_gate_qudits)
        coefficients = np.sqrt(probabilities)

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
    def __init__(self, error_rate: float) -> None:
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
        super().__init__()

    @classmethod
    def n_kraus_operators(cls) -> int:
        return 4

    def kraus_probabilities(self, n_qudits: int) -> np.ndarray:
        """
        Returns the Kraus probabilities for the depolarizing channel.
        Each probability is the product of single-qudit probabilities.
        """
        single_qudit_probs = np.array([
            np.abs(self.alpha_00)**2,
            np.abs(self.alpha_ii)**2,
            np.abs(self.alpha_ii)**2,
            np.abs(self.alpha_ii)**2
        ], dtype=float)

        result = single_qudit_probs
        for _ in range(n_qudits - 1):
            result = np.kron(result, single_qudit_probs)

        return result

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:

        n_qudits = len(dimensions)
        n_gate_qudits = len(qudit_indices)

        # Get coefficients from probabilities (sqrt of probabilities)
        probabilities = self.kraus_probabilities(n_gate_qudits)
        coefficients = np.sqrt(probabilities)

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
