from __future__ import annotations
import numpy as np
import scipy.sparse as sp
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .pauli_sum import PauliSum
    from .pauli_string import PauliString

from .pauli_object import PauliObject
from .constants import DEFAULT_QUDIT_DIMENSION


class Pauli(PauliObject):
    @classmethod
    def from_tableau(cls, tableau: np.ndarray, dimension: int = DEFAULT_QUDIT_DIMENSION) -> Pauli:
        """
        Create a Pauli from its tableau.

        Parameters
        ----------
        tableau : inp.ndarray
            The tableau of the Pauli, a 1D array of length 2.
        dimension : int, optional
            Qudit dimension used to reduce exponents modulo `dimension`.

        Returns
        -------
        Pauli
            Pauli instance with tableau reduced modulo `dimension`.
        """
        P = cls(tableau, dimensions=dimension)
        P._sanity_check()

        return P

    @classmethod
    def from_exponents(cls, x_exp: int | None = None, z_exp: int | None = None,
                       dimension: int = DEFAULT_QUDIT_DIMENSION) -> Pauli:
        """
        Create a Pauli from integer X and Z exponents.

        Parameters
        ----------
        x_exp : int | None
            Exponent for the X part. If None, treated as 0.
        z_exp : int | None
            Exponent for the Z part. If None, treated as 0.
        dimension : int, optional
            Qudit dimension used to reduce exponents modulo `dimension`.

        Returns
        -------
        Pauli
            Pauli instance with exponents reduced modulo `dimension`.
        """

        if x_exp is None:
            x_exp = 0
        if z_exp is None:
            z_exp = 0

        tableau = np.asarray([x_exp, z_exp], dtype=int) % dimension
        P = cls(tableau, dimensions=dimension)
        P._sanity_check()

        return P

    @classmethod
    def from_string(cls, pauli_str: str, dimension: int = DEFAULT_QUDIT_DIMENSION) -> Pauli:
        """
        Create a Pauli from a string representation.

        Parameters
        ----------
        pauli_str : str
            String of the form 'x{int}z{int}', e.g. 'x1z0'. Only single-digit
            exponents are currently supported by this parser.
        dimension : int, optional
            Qudit dimension to interpret the exponents modulo. Defaults to
            DEFAULT_QUDIT_DIMENSION.

        Returns
        -------
        Pauli
            Instance constructed from the parsed exponents.

        Raises
        ------
        ValueError
            If the string cannot be parsed into two integer exponents.
        """
        tableau = np.empty(2, dtype=int)
        xz_exponents = re.split('x|z', pauli_str)[1:]
        x_exp = int(xz_exponents[0])
        z_exp = int(xz_exponents[1])
        tableau[0] = x_exp
        tableau[1] = z_exp

        P = cls(tableau, dimensions=dimension)
        P._sanity_check()

        return P

    @classmethod
    def Xnd(cls, x_exp: int, dimension: int) -> Pauli:
        """
        Create a Pauli X operator with given X exponent.

        Parameters
        ----------
        x_exp : int
            Exponent for the X part.
        dimension : int
            Qudit dimension.

        Returns
        -------
        Pauli
            Pauli instance corresponding to X^x_exp.
        """
        return cls.from_exponents(x_exp, 0, dimension)

    @classmethod
    def Ynd(cls, y_exp: int, dimension: int) -> Pauli:
        """
        Create a Pauli Y-like operator with equal X and Z exponents.

        Parameters
        ----------
        y_exp : int
            Exponent for both X and Z parts.
        dimension : int
            Qudit dimension.

        Returns
        -------
        Pauli
            Pauli instance with X and Z exponents equal to `y_exp`.
        """
        return cls.from_exponents(y_exp, y_exp, dimension)

    @classmethod
    def Znd(cls, z_exp: int, dimension: int) -> Pauli:
        """
        Create a Pauli Z operator with given Z exponent.

        Parameters
        ----------
        z_exp : int
            Exponent for the Z part.
        dimension : int
            Qudit dimension.

        Returns
        -------
        Pauli
            Pauli instance corresponding to Z^z_exp.
        """
        return cls.from_exponents(0, z_exp, dimension)

    @classmethod
    def Idnd(cls, dimension: int) -> Pauli:
        """
        Create the identity Pauli (zero exponents) for given dimension.

        Parameters
        ----------
        dimension : int
            Qudit dimension.

        Returns
        -------
        Pauli
            Identity Pauli (x_exp=0, z_exp=0).
        """
        return cls.from_exponents(0, 0, dimension)

    @property
    def phases(self) -> np.ndarray:
        """
        Returns the phases associated with the Pauli.
        For a Pauli operator, this is just the trivial phase.

        Returns
        -------
        np.ndarray
            The phases as a 1d-vector.
        """
        return np.asarray([0], dtype=int)

    @property
    def weights(self) -> np.ndarray:
        """
        Returns the weights associated with the Pauli.
        For a Pauli operator, this is just 1.

        Returns
        -------
        np.ndarray
            The weights as a 1d-vector.
        """
        return np.asarray([1], dtype=complex)

    @property
    def dimension(self) -> int:
        """
        Returns the dimension of the Pauli as an int.

        Returns
        -------
        int
            The Pauli dimension
        """
        return self._dimensions[0]

    def as_pauli_sum(self) -> PauliSum:
        """
        Converts the Pauli to a PauliSum.

        Returns
        -------
        PauliSum
            A PauliSum instance representing the given Pauli operator.
        """
        return PauliSum(self.tableau, self.dimensions, self.weights, self.phases)

    def as_pauli_string(self) -> PauliString:
        """
        Converts the Pauli to a PauliString.

        Returns
        -------
        PauliString
            A PauliString instance representing the given Pauli operator.
        """
        return PauliString(self.tableau, self.dimensions, self.weights, self.phases)

    def to_hilbert_space(self) -> sp.csr_matrix:
        """
        Get the matrix form of the Pauli as a sparse matrix.

        Returns
        -------
        scipy.sparse.csr_matrix
            Matrix representation of input Pauli.
        """
        return self.as_pauli_sum().to_hilbert_space()

    def __mul__(self, A: str | Pauli) -> Pauli:
        """
        Multiply this Pauli by another Pauli or parseable string.

        Parameters
        ----------
        A : str | Pauli
            Other Pauli or string representation (e.g. 'x1z0').

        Returns
        -------
        Pauli
            Product Pauli with exponents added modulo the dimension.

        Raises
        ------
        Exception
            If operand type is unsupported or dimensions mismatch.
        """
        if isinstance(A, str):
            return self * Pauli.from_string(A)

        if not isinstance(A, Pauli):
            raise Exception(f"Cannot multiply Pauli with type {type(A)}")

        if not np.array_equal(self.dimensions, A.dimensions):
            raise Exception("To multiply two Paulis, their dimensions"
                            f" {A.dimensions} and {self.dimensions} must be equal")

        new_tableau = (self.tableau + A.tableau) % self.lcm

        return Pauli(new_tableau, dimensions=self.dimensions)

    def __repr__(self) -> str:
        """
        Return the string representation of the Pauli.
        (in a format that is helpful for debugging).

        Returns
        -------
        str
            A string in the format "Pauli(x_exp=..., z_exp=..., dimensions=...)".
        """

        return f"PauliString(tableau={self.tableau}, dimensions={self.dimensions})"

    def __str__(self) -> str:
        """
        String representation in the form 'x{X}z{Z}'.

        Returns
        -------
        str
            Human-readable short string for the Pauli.
        """
        return f'x{self.x_exp}z{self.z_exp}'

    @property
    def x_exp(self) -> int:
        """
        X exponent of the Pauli.

        Returns
        -------
        int
            X exponent.
        """
        return self._tableau[0][0]

    @property
    def z_exp(self) -> int:
        """
        Z exponent of the Pauli.

        Returns
        -------
        int
            Z exponent.
        """
        return self._tableau[0][1]
