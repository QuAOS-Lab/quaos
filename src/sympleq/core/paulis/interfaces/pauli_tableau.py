from __future__ import annotations
from abc import ABC
import functools
import numpy as np
from typing import TypeVar, Self, Union, TYPE_CHECKING
if TYPE_CHECKING:
    from ..pauli_sum import PauliSum

from ..constants import DEFAULT_QUDIT_DIMENSION

P = TypeVar("P", bound="PauliTableau")

ScalarType = Union[float, complex, int]
PauliOrScalarType = Union['PauliTableau', ScalarType]


@functools.total_ordering
class PauliTableau(ABC):
    def __init__(self, tableau: np.ndarray, dimensions: int | list[int] | np.ndarray | None = None):
        """
        Initialize a PauliTableau represented in symplectic tableau form.

        Represents a sum of Pauli operators acting on multiple qudits.
        See references:
            - Phys. Rev. A 71, 042315 (2005)
            - Phys. Rev. A 70, 052328 (2004)

        Parameters
        ----------
        tableau : np.ndarray
            Symplectic tableau representation of the object with shape
            (n_paulis, 2 * n_qudits). The first half corresponds to X
            exponents, and the second half to Z exponents.
        dimensions : int | list[int] | np.ndarray | None, optional
            Qudit dimension(s). If int, all qudits share the same dimension.
            If list or array, its length must equal the number of qudits.
            Defaults to `DEFAULT_QUDIT_DIMENSION` if None.
        weights : int | float | complex | list | np.ndarray | None, optional
            Coefficients associated with each Pauli term. Defaults to 1.
        phases : int | list[int] | np.ndarray | None, optional
            Integer phase factors for each Pauli term, typically modulo
            `2 * lcm(dimensions)`. Defaults to 0.

        Attributes
        ----------
        tableau : np.ndarray
            Symplectic tableau representing the Pauli operators.
        dimensions : np.ndarray
            Dimensions of each qudit.
        weights : np.ndarray
            Coefficients (amplitudes) for each Pauli term.
        phases : np.ndarray
            Integer phases for each Pauli term.
        lcm : int
            Least common multiple of all qudit dimensions.
        """

        if tableau.ndim == 1:
            tableau = tableau.reshape(1, -1)

        if tableau.ndim != 2:
            raise ValueError(f"Invalid tableau shape ({tableau.shape}). Tableaus should be two dimensional.")

        n_qudits = tableau.shape[1] // 2

        if dimensions is None:
            dimensions = np.ones(n_qudits, dtype=int) * DEFAULT_QUDIT_DIMENSION
        else:  # Catches int but also list and arrays of length 1
            dimensions = np.asarray(dimensions, dtype=int)
            if dimensions.ndim == 0:
                dimensions = np.full(n_qudits, dimensions.item(), dtype=int)

        self._dimensions = dimensions
        # Dimensions is read-only, so any type of assignment will fail.
        self._dimensions.setflags(write=False)
        self._lcm = int(np.lcm.reduce(self.dimensions))

        self._tableau = tableau % np.tile(self.dimensions, 2)

    @property
    def tableau(self) -> np.ndarray:
        """
        Return the symplectic tableau representation of the Pauli object.

        The tableau is a 2D array of shape (n_paulis, 2 * n_qudits),
        where the first half represents X exponents and the second half Z exponents.

        Returns
        -------
        np.ndarray
            Tableau representation of the Pauli object.
        """
        return self._tableau

    @property
    def dimensions(self) -> np.ndarray:
        """
        Return the dimensions of each qudit.

        Returns
        -------
        np.ndarray
            A 1D array of qudit dimensions of length `n_qudits`.
        """
        return self._dimensions

    @dimensions.setter
    def dimensions(self, value: np.ndarray):
        # Dimensions is read-only, and this setter is not strictly required.
        # We keep it to raise with a meaningful error message.
        raise Exception("The dimensions of a PauliTableau cannot be set.\
                        If you want to change the PauliTableau dimensions, generate a new one.")

    @property
    def lcm(self) -> int:
        """
        Return the least common multiple (LCM) of all qudit dimensions.

        Returns
        -------
        int
            The least common multiple of the qudit dimensions.
        """
        return self._lcm

    @lcm.setter
    def lcm(self, value: int):
        raise Exception("The lcm of a PauliTableau cannot be set, as it is derived from its dimensions.\
                        If you want to change the PauliTableau dimensions, generate a new one.")

    def n_qudits(self) -> int:
        """
        Return the number of qudits represented by the Pauli object.

        Returns
        -------
        int
            Number of qudits.
        """
        return len(self.dimensions)

    def n_paulis(self) -> int:
        """
        Return the number of Pauli terms in the object.

        Returns
        -------
        int
            Number of Pauli operators (rows in the tableau).
        """
        return len(self.tableau)

    def shape(self) -> tuple[int, int]:
        """
        Return the shape of the Pauli object.

        Returns
        -------
        tuple[int, int]
            Tuple of (n_paulis, n_qudits).
        """
        return self.n_paulis(), self.n_qudits()

    def has_equal_tableau(self, other_pauli: PauliTableau) -> bool:
        """
        Check whether two Pauli objects have the same tableau and dimensions.

        Parameters
        ----------
        other_pauli : PauliTableau
            Pauli object to compare against.

        literal : bool, optional
            If True, compares objects literally in their current form. If False,
            the objects are first brought to standard form. Default is True.

        Returns
        -------
        bool
            True if all tableau entries and dimensions match; False otherwise.
        """

        if not np.array_equal(self.dimensions, other_pauli.dimensions):
            return False

        if not np.array_equal(self.tableau, other_pauli.tableau):
            return False

        return True

    def _sanity_check(self):
        """
        Validate internal consistency of the Pauli object.

        Raises
        ------
        ValueError
            If tableau, dimensions, or exponents are inconsistent or invalid.
        """

        if len(self.tableau[0]) != 2 * self.n_qudits():
            raise ValueError(f"Tableau ({len(self.tableau)}) should be twice as long as"
                             f"dimensions ({len(self.dimensions)}).")

        if np.any(self.dimensions < DEFAULT_QUDIT_DIMENSION):
            bad_dims = self.dimensions[self.dimensions < DEFAULT_QUDIT_DIMENSION]
            raise ValueError(f"Dimensions {bad_dims} are less than {DEFAULT_QUDIT_DIMENSION}")

        d = np.tile(self.dimensions, 2)
        if np.any((self.tableau >= d)):
            bad_indices = np.where((self.tableau >= d))[0]
            raise ValueError(
                f"Exponents at indices {bad_indices} are too large:"
                f"tableau={self.tableau[bad_indices]}"
            )

        if np.any((self.tableau < 0)):
            bad_indices = np.where((self.tableau < 0))[0]
            raise ValueError(
                f"Exponents at indices {bad_indices} are negative:"
                f"tableau={self.tableau[bad_indices]}"
            )

    def __repr__(self) -> str:
        """
        Returns an unambiguous string representation of the PauliTableau.

        Returns
        -------
        str
            A string representation of the PauliTableau.
        """
        return f'{self.__class__.__name__}({self.tableau}, {self.dimensions})'

    def __eq__(self, other_pauli: Self) -> bool:
        """
        Determine if two Pauli objects are equal.

        Parameters
        ----------
        other_pauli : PauliTableau
            Object to compare with.

        Returns
        -------
        bool
            True if tableau, weights, phases, and dimensions match exactly;
            False otherwise.
        """
        if not isinstance(other_pauli, self.__class__):
            return False

        if not np.array_equal(self.tableau, other_pauli.tableau):
            return False

        if not np.array_equal(self.dimensions, other_pauli.dimensions):
            return False

        return True

    def __ne__(self, other_pauli: PauliTableau) -> bool:
        """
        Determine if two Pauli objects are different.

        Parameters
        ----------
        other_pauli : PauliTableau
            Object to compare with.

        Returns
        -------
        bool
            True if objects are not equal; False otherwise.
        """
        return not self == other_pauli

    def __gt__(self, other_pauli: PauliTableau) -> bool:
        """
        Strict greater-than comparison for ordering single Pauli terms.

        Behavior
        --------
        This operator is intended to impose an ordering on single-term Pauli
        objects by interpreting their tableau as an integer (or comparable)
        representation. It is undefined for multi-term objects.

        Parameters
        ----------
        other_pauli : PauliTableau
            Other Pauli object to compare against. Both objects must represent
            a single Pauli term and share identical `dimensions`.

        Returns
        -------
        bool
            True if `self` is greater than `other_pauli` according to the
            implemented integer-like ordering; False otherwise.

        Raises
        ------
        ValueError
            If either object contains multiple Pauli terms or if dimensions differ.

        Examples
        --------
        >>> ps1 = PauliString.from_string("x0z1 x1z0", [2, 2])
        >>> ps2 = PauliString.from_string("x1z0 x0z1", [2, 2])
        >>> ps1 > ps2
        True
        """

        if self.n_paulis() > 1:
            raise Exception("A Pauli object with more than one Pauli objects cannot be ordered.")

        if not np.array_equal(self.dimensions, other_pauli.dimensions):
            raise Exception("Cannot compare Pauli objects with different dimensions.")

        # Flatten tableaus to 1D-vectors
        self_tableau = self.tableau.ravel()
        other_tableau = other_pauli.tableau.ravel()

        for i in range(len(self_tableau)):
            if self_tableau[i] == other_tableau[i]:
                continue
            if self_tableau[i] < other_tableau[i]:
                return True
            return False

        # They are equal
        return False

    def __lt__(self, other_pauli: Self) -> bool:
        """
        Strict less-than comparison for ordering single Pauli terms.

        Parameters
        ----------
        other_pauli : PauliTableau
            Other Pauli object to compare against. Both must represent single terms.

        Returns
        -------
        bool
            True if `self` is less than `other_pauli` according to the implemented ordering.

        Raises
        ------
        ValueError
            If preconditions (single-term objects, matching dimensions) are not met.

        Examples
        --------
        >>> ps1 = PauliString.from_string("x1z0 x0z1", [2, 2])
        >>> ps2 = PauliString.from_string("x0z1 x1z0", [2, 2])
        >>> ps1 > ps2
        False
        """
        return not self.__gt__(other_pauli) and not self.__eq__(other_pauli)

    def __add__(self, A: PauliTableau) -> PauliSum:
        """
        Implements the addition of Pauli objects.

        Parameters
        ----------
        A : PauliTableau
            The Pauli operator to add.

        Returns
        -------
        PauliSum
            A new PauliSum instance representing the sum of `self` and `A`.

        Examples
        --------
        >>> p1 = PauliSum.from_pauli_strings("x1z0 x0z1", [3, 2])
        >>> p2 = PauliSum.from_pauli_strings("x2z1 x1z1", [3, 2])
        >>> p1 + p2
        PauliSum(...)

        Raises
        ------
        ValueError
            If the dimensions of `self` and `A` do not match.

        Notes
        -----
        - Dimensions must agree!
        """

        if not np.array_equal(self.dimensions, A.dimensions):
            raise ValueError(f"The dimensions of the Pauli objects do not match ({self.dimensions}, {A.dimensions}).")

        new_tableau = np.vstack([self.tableau, A.tableau])

        from ..pauli_sum import PauliSum
        return PauliSum(new_tableau, self.dimensions.copy())

    def __radd__(self, A: PauliTableau) -> PauliSum:
        """
        Implements the addition of Pauli objects.

        Parameters
        ----------
        A : PauliTableau
            The Pauli operator to add.

        Returns
        -------
        PauliSum
            A new PauliSum instance representing the sum of `self` and `A`.

        Examples
        --------
        >>> p1 = PauliSum.from_pauli_strings("x1z0 x0z1", [3, 2])
        >>> p2 = PauliSum.from_pauli_strings("x2z1 x1z1", [3, 2])
        >>> p1 + p2
        PauliSum(...)

        Raises
        ------
        ValueError
            If the dimensions of `self` and `A` do not match.

        Notes
        -----
        - Dimensions must agree!
        """

        return self + A

    def __sub__(self, A: PauliTableau) -> PauliSum:
        """
        Implements the subtraction of Pauli objects.

        Parameters
        ----------
        A : PauliTableau
            The Pauli operator to subtract.

        Returns
        -------
        PauliSum
            A new PauliSum instance representing the difference of `self` and `A`.

        Examples
        --------
        >>> p1 = PauliSum.from_pauli_strings("x1z0 x0z1", [3, 2])
        >>> p2 = PauliSum.from_pauli_strings("x2z1 x1z1", [3, 2])
        >>> p1 - p2
        PauliSum(...)

        Raises
        ------
        ValueError
            If the dimensions of `self` and `A` do not match.

        Notes
        -----
        - Dimensions must agree!
        """

        if not np.array_equal(self.dimensions, A.dimensions):
            raise ValueError(f"The dimensions of the Pauli objects do not match ({self.dimensions}, {A.dimensions}).")

        new_tableau = np.vstack([self.tableau, A.tableau])

        from ..pauli_sum import PauliSum
        return PauliSum(new_tableau, self.dimensions.copy())

    def __rsub__(self, A: PauliTableau) -> PauliSum:
        """
        Implements the subtraction of Pauli objects.

        Parameters
        ----------
        A : PauliTableau
            The Pauli operator to subtract.

        Returns
        -------
        PauliSum
            A new PauliSum instance representing the difference of `self` and `A`.

        Examples
        --------
        >>> p1 = PauliSum.from_pauli_strings("x1z0 x0z1", [3, 2])
        >>> p2 = PauliSum.from_pauli_strings("x2z1 x1z1", [3, 2])
        >>> p1 - p2
        PauliSum(...)

        Raises
        ------
        ValueError
            If the dimensions of `self` and `A` do not match.

        Notes
        -----
        - Dimensions must agree!
        """

        return self - A

    def __pow__(self, A: int) -> Self:
        """
        Integer power of a Pauli object.

        Parameters
        ----------
        A : int
            Exponent to raise the Pauli object to. Typically only defined for
            single-term Pauli objects; behavior for sums depends on implementation.

        Returns
        -------
        Self
            Resulting PauliTableau after exponentiation.

        Raises
        ------
        ValueError
            If exponentiation is undefined for the current object (e.g., multi-term).

        Examples
        --------
        >>> ps = PauliString.from_exponents(x_exp, z_exp, dimensions)
        >>> ps_squared = ps ** 2
        """

        if self.n_paulis() > 1:
            raise Exception("A Pauli object with more than a PauliString cannot be exponentiated.")

        tableau = np.mod(self.tableau * A, np.tile(self.dimensions, 2))
        return self.__class__(tableau, self.dimensions.copy())

    def __hash__(self) -> int:
        """
        Return the hash value of the Pauli object. That is a unique identifier.

        Returns
        -------
        int
            The hash value of the Pauli object instance.
        """
        return hash(
            (tuple(self.tableau),
             tuple(self.dimensions))
        )

    def __dict__(self) -> dict:
        """
        Returns a dictionary representation of the object's attributes.

        Returns
        -------
        dict
            A dictionary containing the values of `tableau`, `weights`, `phases`, and `dimensions`.
        """
        return {'tableau': self.tableau,
                'dimensions': self.dimensions,
                }

    def copy(self) -> Self:
        """
        Creates a copy of the Pauli object.

        Returns
        -------
        Pauli object
            A copy of the Pauli object.
        """
        return self.__class__(self.tableau.copy(), self.dimensions.copy())
