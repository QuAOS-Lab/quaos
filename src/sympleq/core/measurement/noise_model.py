import numpy as np
from abc import ABC, abstractmethod
from itertools import product
import math

from sympleq.core.paulis.pauli_sum import PauliSum


class NoiseModel(ABC):
    """
    See 10.1103/PhysRevA.108.062604.
    """
    def __init__(self) -> None:
        pass

    @classmethod
    @abstractmethod
    def n_kraus_operators(cls) -> int:
        pass

    @abstractmethod
    def apply_kraus_operator(self, pauli: PauliSum, qudit_indices: list[int] | np.ndarray,
                             index: int) -> PauliSum:
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
    def process_matrix(self, pauli_i: list[int], pauli_j: list[int]) -> complex:
        """"
        lambda_{ij} in Eq.(3), where i, j are multi-indices.
        The inout pauli_indices is a list of pauli operators indices, one per qudit,
        which are mapped to multi-index to calculate the process matrix.
        """
        pass


class Noiseless(NoiseModel):
    @classmethod
    def n_kraus_operators(cls) -> int:
        return 1

    def apply_kraus_operator(self, pauli: PauliSum, qudit_indices: list[int] | np.ndarray,
                             index: int) -> PauliSum:

        return pauli

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:

        # The noiseless channel has a single Klaus operator (the identity), wuth probability 1.
        tableau = np.zeros(2 * len(dimensions), dtype=int)
        sigma_0 = PauliSum.from_tableau(tableau, dimensions)
        return [sigma_0]

    def process_matrix(self, pauli_i: list[int], pauli_j: list[int]) -> complex:
        # The only entry in the process matrix is lambda_{00} (all qudit indices are zero).
        if sum(pauli_i) == 0 and sum(pauli_j) == 0:
            return 1.0

        return 0.0


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
        super().__init__()

    @classmethod
    def n_kraus_operators(cls) -> int:
        return 2

    def apply_kraus_operator(self, pauli: PauliSum, qudit_indices: list[int] | np.ndarray,
                             index: int) -> PauliSum:

        n_qudits = len(pauli.dimensions)
        k = len(qudit_indices)
        kraus_tableau = np.zeros(2 * n_qudits, dtype=int)
        kraus_phases = np.zeros(1, dtype=int)

        # Interpret index as binary digits for each qudit, 0 = I, 1 = Z
        for qudit_pos, qudit_idx in enumerate(qudit_indices):
            bit = (index >> (k - qudit_pos - 1)) & 1
            kraus_tableau[n_qudits + qudit_idx] = bit

        p1 = pauli.phases
        p2 = kraus_phases

        kraus_tableau = kraus_tableau.reshape(1, -1)

        # Extract z- and x-parts from tableau
        n = pauli.n_qudits()
        a_x = pauli.tableau[:, :n]
        a_z = pauli.tableau[:, n:]
        b_x = kraus_tableau[:, :n]
        b_z = kraus_tableau[:, n:]

        # Compute acquired phases via symplectic form
        factors = (pauli.lcm // pauli.dimensions)
        left_acquired_phases = 2 * factors * a_z  @ b_x.T
        right_acquired_phases = 2 * factors * b_z  @ a_x.T

        # Combine with existing phases and flatten
        new_phases = (p1 + p2 + p2 + left_acquired_phases + right_acquired_phases) % (2 * pauli.lcm)
        new_phases = new_phases.reshape(-1)

        # Multiplication between PauliString corresponds to summing the tableaus
        new_tableau = np.asarray([ps1 + 2 * ps2 for ps1 in pauli.tableau for ps2 in kraus_tableau], dtype=int)

        return PauliSum(new_tableau, pauli.dimensions, pauli.weights, new_phases)

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:

        n_qudits = len(dimensions)

        tableau = np.zeros(2 * n_qudits, dtype=int)
        output = []

        # Build all combinations of the two kraus operators for each qudit in qudit_indices
        combinations = list(product([0, 1], repeat=len(qudit_indices)))
        for comb in combinations:
            p_tableau = tableau.copy()
            coefficient = 1
            for idx, exp in zip(qudit_indices, comb):
                p_tableau[n_qudits + idx] = exp
                coefficient *= self.p0 if exp == 0 else (1.0 - self.p0)
                sigma = PauliSum.from_tableau(p_tableau, dimensions)
            output.append(math.sqrt(coefficient) * sigma)

        return output

    def process_matrix(self, pauli_i: list[int], pauli_j: list[int]) -> complex:
        lambda_0 = math.sqrt(self.p0)
        lambda_1 = math.sqrt(1.0 - self.p0)

        def map_pauli_index_to_coefficient(idx: int) -> float:
            if idx == 0:
                return lambda_0

            if idx == 3:
                return lambda_1

            return 0.0

        return math.prod(list(map(map_pauli_index_to_coefficient, pauli_i))) * \
            math.prod(list(map(map_pauli_index_to_coefficient, pauli_j)))


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
        super().__init__()

    @classmethod
    def n_kraus_operators(cls) -> int:
        return 4

    def apply_kraus_operator(self, pauli: PauliSum, qudit_indices: list[int] | np.ndarray,
                             index: int) -> PauliSum:

        n_qudits = len(pauli.dimensions)
        k = len(qudit_indices)
        kraus_tableau = np.zeros(2 * n_qudits, dtype=int)
        kraus_phases = np.zeros(1, dtype=int)

        # Convert index to base-4 digits
        for qudit_pos, qudit_idx in enumerate(qudit_indices):
            # Extract base-4 digit
            digit = (index // (4 ** (k - qudit_pos - 1))) % 4
            # identity
            if digit == 0:
                pass
            # sigma-x
            elif digit == 1:
                kraus_tableau[qudit_idx] = 1
            # sigma-y
            elif digit == 2:
                kraus_tableau[qudit_idx] = 1
                kraus_tableau[n_qudits + qudit_idx] = 1
                kraus_phases[0] = 1
            # sigma-z
            elif digit == 3:
                kraus_tableau[n_qudits + qudit_idx] = 1

        p1 = pauli.phases[:, None]
        p2 = kraus_phases[None, :]

        # Extract z- and x-parts from tableau
        n = pauli.n_qudits()
        a_z = pauli.tableau[:, n:]
        b_x = kraus_tableau[:, :n]

        # Compute acquired phases via symplectic form
        factors = (pauli.lcm // pauli.dimensions)
        acquired_phases = 2 * factors * a_z  @ b_x.T

        # Combine with existing phases and flatten
        new_phases = (p1 + p2 + acquired_phases) % (2 * pauli.lcm)
        new_phases = new_phases.reshape(-1)

        # Multiplication between PauliString corresponds to summing the tableaus
        new_tableau = np.asarray([ps1 + ps2 for ps1 in pauli.tableau for ps2 in kraus_tableau], dtype=int)

        return PauliSum(new_tableau, pauli.dimensions, pauli.weights, new_phases)

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:

        n_qudits = len(dimensions)

        tableau = np.zeros(2 * n_qudits, dtype=int)
        output = []

        combinations = list(product([0, 1, 2, 3], repeat=len(qudit_indices)))
        for comb in combinations:
            p_tableau = tableau.copy()
            phases = np.zeros(1, dtype=int)
            coefficient = 1
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

                coefficient *= self.p0 if exp == 0 else (1.0 - self.p0) / 3.0
                sigma = PauliSum.from_tableau(p_tableau, dimensions, phases=phases)
            output.append(math.sqrt(coefficient) * sigma)

        return output

    def process_matrix(self, pauli_i: list[int], pauli_j: list[int]) -> complex:
        lambda_0 = math.sqrt(self.p0)
        lambda_j = math.sqrt((1.0 - self.p0) / 3)

        def map_pauli_index_to_coefficient(idx: int) -> float:
            if idx == 0:
                return lambda_0

            if idx in [1, 2, 3]:
                return lambda_j

            raise ValueError(f"Invalid pauli index {idx}")

        return math.prod(list(map(map_pauli_index_to_coefficient, pauli_i))) * \
            math.prod(list(map(map_pauli_index_to_coefficient, pauli_j)))
