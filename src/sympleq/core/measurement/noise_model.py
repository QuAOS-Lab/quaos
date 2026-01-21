import numpy as np
from abc import ABC, abstractmethod
from itertools import product

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
    def kraus_operator_probabilities(self, n_qudits: int) -> list[float]:
        pass

    @classmethod
    @abstractmethod
    def kraus_operator(cls,
                       dimensions: list[int] | np.ndarray,
                       qudit_indices: list[int] | np.ndarray,
                       index: int) -> PauliSum:
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
        pass


class Noiseless(NoiseModel):
    @classmethod
    def n_kraus_operators(cls) -> int:
        return 1

    def kraus_operator_probabilities(self, n_qudits: int) -> list[float]:
        return [1.0]

    def apply_kraus_operator(self, pauli: PauliSum, qudit_indices: list[int] | np.ndarray,
                             index: int) -> PauliSum:

        return pauli

    @classmethod
    def kraus_operator(cls,
                       dimensions: list[int] | np.ndarray,
                       qudit_indices: list[int] | np.ndarray,
                       index: int) -> PauliSum:

        tableau = np.zeros(2 * len(dimensions), dtype=int)
        return PauliSum.from_tableau(tableau, dimensions)

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[tuple[complex, PauliSum]]:

        # The noiseless channel has a single Klaus operator (the identity), wuth probability 1.
        tableau = np.zeros(2 * len(dimensions), dtype=int)
        return [(1.0, PauliSum.from_tableau(tableau, dimensions))]


class DephasingNoise(NoiseModel):
    def __init__(self, error_rate: float) -> None:
        # The dephasing channel has tow Klaus operators: the identity and Z.
        # A single parameter p0 models the noise probability,
        # representing the probability of having no error.
        # Kraus operators:
        # K0 = sqrt(p0) 1 ⊗ 1
        # K1 = sqrt(1 − p0) Z ⊗ Z
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        super().__init__()

    @classmethod
    def n_kraus_operators(cls) -> int:
        return 2

    def kraus_operator_probabilities(self, n_qudits: int) -> list[float]:
        sqrt_p0 = self.p0
        sqrt_p1 = (1.0 - self.p0)

        p0_pows = [sqrt_p0 ** i for i in range(n_qudits + 1)]
        p1_pows = [sqrt_p1 ** i for i in range(n_qudits + 1)]

        output = []

        for comb in product((0, 1), repeat=n_qudits):
            m = sum(exp != 0 for exp in comb)
            output.append(p0_pows[n_qudits - m] * p1_pows[m])

        assert sum(output) == 1

        return output

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

    @classmethod
    def kraus_operator(cls,
                       dimensions: list[int] | np.ndarray,
                       qudit_indices: list[int] | np.ndarray,
                       index: int) -> PauliSum:

        n_qudits = len(dimensions)
        k = len(qudit_indices)
        tableau = np.zeros(2 * n_qudits, dtype=int)

        # Interpret index as binary digits for each qudit, 0 = I, 1 = Z
        for qudit_pos, qudit_idx in enumerate(qudit_indices):
            bit = (index >> (k - qudit_pos - 1)) & 1
            tableau[n_qudits + qudit_idx] = bit

        return PauliSum.from_tableau(tableau, dimensions)

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[PauliSum]:

        n_qudits = len(dimensions)

        tableau = np.zeros(2 * n_qudits, dtype=int)
        output = []

        combinations = list(product([0, 1], repeat=len(qudit_indices)))
        for comb in combinations:
            p_tableau = tableau.copy()
            coefficient = 1
            for idx, exp in zip(qudit_indices, comb):
                p_tableau[n_qudits + idx] = exp
                coefficient *= np.sqrt(self.p0) if exp == 0 else np.sqrt(1.0 - self.p0)
            output.append((coefficient, PauliSum.from_tableau(p_tableau, dimensions)))

        return output


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

    def kraus_operator_probabilities(self, n_qudits: int) -> list[float]:
        sqrt_p0 = self.p0
        sqrt_p1 = (1.0 - self.p0) / 3.0

        p0_pows = [sqrt_p0 ** i for i in range(n_qudits + 1)]
        p1_pows = [sqrt_p1 ** i for i in range(n_qudits + 1)]

        output = []

        for comb in product((0, 1, 2, 3), repeat=n_qudits):
            m = sum(exp != 0 for exp in comb)
            output.append(p0_pows[n_qudits - m] * p1_pows[m])

        assert sum(output) == 1

        return output

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

    @classmethod
    def kraus_operator(cls,
                       dimensions: list[int] | np.ndarray,
                       qudit_indices: list[int] | np.ndarray,
                       index: int) -> PauliSum:

        n_qudits = len(dimensions)
        k = len(qudit_indices)
        tableau = np.zeros(2 * n_qudits, dtype=int)
        phases = np.zeros(1, dtype=int)

        # Convert index to base-4 digits
        for qudit_pos, qudit_idx in enumerate(qudit_indices):
            # Extract base-4 digit
            digit = (index // (4 ** (k - qudit_pos - 1))) % 4

            # identity
            if digit == 0:
                pass
            # sigma-x
            elif digit == 1:
                tableau[qudit_idx] = 1
            # sigma-y
            elif digit == 2:
                tableau[qudit_idx] = 1
                tableau[n_qudits + qudit_idx] = 1
                phases[0] = 1
            # sigma-z
            elif digit == 3:
                tableau[n_qudits + qudit_idx] = 1

        return PauliSum.from_tableau(tableau, dimensions, phases=phases)

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

                coefficient *= np.sqrt(self.p0) if exp == 0 else np.sqrt((1.0 - self.p0) / 3.0)
            output.append((coefficient, PauliSum.from_tableau(p_tableau, dimensions, phases=phases)))

        return output
