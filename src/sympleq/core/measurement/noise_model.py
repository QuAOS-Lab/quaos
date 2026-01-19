import numpy as np
from abc import ABC, abstractmethod
from itertools import product
# from typing import NewType

from sympleq.core.paulis.pauli_sum import PauliSum

# KrausOperator = NewType("KrausOperator", tuple[float, PauliSum])


class NoiseModel(ABC):
    """
    See 10.1103/PhysRevA.108.062604.
    """
    def __init__(self) -> None:
        pass

    @classmethod
    @abstractmethod
    def kraus_operators(cls,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[tuple[complex, PauliSum]]:
        raise NotImplementedError


class Noiseless(NoiseModel):
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
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        super().__init__()

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[tuple[complex, PauliSum]]:

        n_qudits = len(dimensions)

        # E.g.: With m = 2 one can construct the map with the Kraus operators:
        # K0 = p0 1 ⊗ 1
        # K1 = sqrt(p0 * (1 − p0)) 1 ⊗ Z
        # K2 = sqrt(p0 * (1 − p0)) Z ⊗ 1
        # K3 = (1 − p0) Z ⊗ Z

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
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        super().__init__()

    def kraus_operators(self,
                        dimensions: list[int] | np.ndarray,
                        qudit_indices: list[int] | np.ndarray) -> list[tuple[complex, PauliSum]]:

        n_qudits = len(dimensions)

        # E.g.: With m = 2 one can construct the map with the Kraus operators:
        # K0 = sqrt(p0) 1
        # Ki =  sqrt(1 − p0/3) σi

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
