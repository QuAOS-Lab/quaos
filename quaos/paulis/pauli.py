from __future__ import annotations
import numpy as np
from typing import Any


def get_sanitized_x_exp(x_exp: int | np.integer | None) -> int:
    if x_exp is None:
        return 0

    if not isinstance(x_exp, (int, np.integer)):
        raise TypeError("x_exp must be an integer type.")
    return int(x_exp)


def get_sanitized_z_exp(z_exp: int | np.integer | None) -> int:
    if z_exp is None:
        return 0

    if not isinstance(z_exp, (int, np.integer)):
        raise TypeError("z_exp must be an integer type.")
    return int(z_exp)


def get_sanitized_dimension(dimension: int | np.integer,
                            x_exp: int,
                            z_exp: int) -> int:
    if not isinstance(dimension, (int, np.integer)):
        raise TypeError("dimension must be an integer type.")

    if dimension - 1 < x_exp or dimension - 1 < z_exp:
        raise ValueError(f"Dimension {dimension} is too small for exponents {x_exp} and {z_exp}")

    return int(dimension)


class Pauli:
    def __init__(self,
                 x_exp: int | None = None,
                 z_exp: int | None = None,
                 dimension: int = 2):

        """
        Constructor for Pauli class.

        Parameters
        ----------
        x_exp : int or str
            Exponent of X part of Pauli in symplectic form. If str, this describes x and z parts in form
            'xnzm', where n and m are integers representing the exponents of x and z respectively.
        z_exp : int
            Exponent of Z part of Pauli in symplectic form. If None, this is set to 0.
        dimension : int
            The dimension of the qudit. Default is 2.
        """
        self.x_exp = get_sanitized_x_exp(x_exp)
        self.z_exp = get_sanitized_z_exp(z_exp)
        self.dimension = get_sanitized_dimension(dimension, self.x_exp, self.z_exp)

    @classmethod
    def from_string(cls, pauli_str: str, dimension: int = 2) -> Pauli:
        return cls(x_exp=int(pauli_str[1]), z_exp=int(pauli_str[3]), dimension=dimension)

    def __mul__(self, A: str | Pauli) -> Pauli:
        if isinstance(A, str):
            return self * Pauli.from_string(A)
        elif isinstance(A, Pauli):
            if A.dimension != self.dimension:
                raise Exception("To multiply two Paulis, their dimensions"
                                f" {A.dimension} and {self.dimension} must be equal")

            return Pauli(x_exp=(self.x_exp + A.x_exp) % self.dimension,
                         z_exp=(self.z_exp + A.z_exp) % self.dimension,
                         dimension=self.dimension)
        else:
            raise Exception(f"Cannot multiply Pauli with type {type(A)}")

    def __str__(self) -> str:
        return f'x{self.x_exp}z{self.z_exp}'

    def __eq__(self, other_pauli: Any) -> bool:
        if not isinstance(other_pauli, Pauli):
            return False
        return self.x_exp == other_pauli.x_exp and self.z_exp == other_pauli.z_exp and self.dimension == other_pauli.dimension

    def __ne__(self, other_pauli: Any) -> bool:
        return not self.__eq__(other_pauli)

    def __dict__(self) -> dict:
        return {'x_exp': self.x_exp, 'z_exp': self.z_exp, 'dimension': self.dimension}

    def __gt__(self, other_pauli: Pauli) -> bool:
        d = self.dimension
        x_measure = min(self.x_exp % d, (d - self.x_exp) % d)
        x_measure_new = min(other_pauli.x_exp % d, (d - other_pauli.x_exp) % d)
        z_measure = min(self.z_exp % d, (d - self.z_exp) % d)
        z_measure_new = min(other_pauli.z_exp % d, (d - other_pauli.z_exp) % d)

        if x_measure > x_measure_new:
            return True
        elif x_measure == x_measure_new:
            if z_measure > z_measure_new:
                return True
            elif z_measure == z_measure_new:
                return False

        return False

    def copy(self) -> Pauli:
        return Pauli(x_exp=self.x_exp, z_exp=self.z_exp, dimension=self.dimension)


class Xnd(Pauli):
    def __init__(self, x_exp: int, dimension: int):
        super().__init__(x_exp, 0, dimension)


class Ynd(Pauli):
    def __init__(self, y_exp: int, dimension: int):
        super().__init__(y_exp, y_exp, dimension)


class Znd(Pauli):
    def __init__(self, z_exp: int, dimension: int):
        super().__init__(0, z_exp, dimension)


class Id(Pauli):
    def __init__(self, dimension: int):
        super().__init__(0, 0, dimension)
