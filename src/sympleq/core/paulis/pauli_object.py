from __future__ import annotations
from abc import ABC, abstractmethod
import functools
import numpy as np
from typing import TypeVar, Self, Union, TYPE_CHECKING
if TYPE_CHECKING:
    from .pauli_sum import PauliSum

from .constants import DEFAULT_QUDIT_DIMENSION

P = TypeVar("P", bound="PauliObject")

ScalarType = Union[float, complex, int]
PauliOrScalarType = Union['PauliObject', ScalarType]


@functools.total_ordering
class PauliObject(ABC):
    def __init__(self, tableau: np.ndarray, dimensions: int | list[int] | np.ndarray | None = None,
                 weights: int | float | complex | list[int | float | complex] | np.ndarray | None = None,
                 phases: int | list[int] | np.ndarray | None = None):
        """
        Initialize a PauliObject represented in symplectic tableau form.

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

        n_pauli_strings = tableau.shape[0]
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

        if weights is None:
            weights = np.ones(n_pauli_strings, dtype=complex)
        else:  # Catches scalars but also list and arrays of length 1
            weights = np.asarray(weights, dtype=complex)
            if weights.ndim == 0:
                weights = np.full(n_pauli_strings, weights.item(), dtype=complex)
        self._weights = weights

        if phases is None:
            phases = np.zeros(n_pauli_strings, dtype=int)
        else:  # Catches scalars but also list and arrays of length 1
            phases = np.asarray(phases, dtype=int)
            if phases.ndim == 0:
                phases = np.full(n_pauli_strings, phases.item(), dtype=int)
        self._phases = phases % (2 * self.lcm)

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
        raise Exception("The dimensions of a PauliObject cannot be set.\
                        If you want to change the PauliObject dimensions, generate a new one.")

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
        raise Exception("The lcm of a PauliObject cannot be set, as it is derived from its dimensions.\
                        If you want to change the PauliObject dimensions, generate a new one.")

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

    @property
    def phases(self) -> np.ndarray:
        """
        Return the integer phases associated with the Pauli object.

        Phases represent numerators, where denominators are `2 * self.lcm`.

        Returns
        -------
        np.ndarray
            1D array of integer phase values modulo `2 * lcm`.
        """
        return self._phases

    @phases.setter
    def phases(self, new_phases: list[int] | np.ndarray):
        new_phases = np.asarray(new_phases, dtype=int)

        if len(new_phases) != self.n_paulis():
            raise ValueError(
                f"New phases ({len(new_phases)}) length must equal the number of Pauli strings ({self.n_paulis()}).")

        self._phases = new_phases % (2 * self.lcm)

    def set_phases(self, new_phases: list[int] | np.ndarray):
        """
        Set new integer phases for the Pauli object.

        Parameters
        ----------
        new_phases : list[int] | np.ndarray
            1D array or list of new integer phase values.
        """
        self.phases = new_phases

    def reset_phases(self):
        """
        Reset all phase values to zero.

        This replaces the internal phase vector with an array of zeros
        of the same length, effectively removing all accumulated phase
        contributions.

        Returns
        -------
        None
            This method modifies the object in place.
        """
        self._phases = np.zeros(len(self.phases))

    @property
    def weights(self) -> np.ndarray:
        """
        Return the weights (coefficients) of the PauliString terms.

        Returns
        -------
        np.ndarray
            1D array of scalar complex coefficients.
        """
        return self._weights

    @weights.setter
    def weights(self, new_weights: list[complex] | np.ndarray):
        new_weights = np.asarray(new_weights, dtype=complex)

        if len(new_weights) != self.n_paulis():
            raise ValueError(
                f"New weights ({len(new_weights)}) length must equal the number of Pauli strings ({self.n_paulis()}).")

        self._weights = new_weights

    def set_weights(self, new_weights: list[complex] | np.ndarray):
        """
        Set new weights (coefficients) for the Pauli object.

        Parameters
        ----------
        new_weights : list[int] | np.ndarray
            1D array or list of new scalar coefficients.
        """
        self.weights = new_weights

    def reset_weights(self):
        """
        Reset all weight coefficients to one.

        This replaces the internal weight vector with an array of ones
        of the same length, restoring the default uniform weighting of
        all PauliString terms.

        Returns
        -------
        None
            This method modifies the object in place.
        """
        self._weights = np.ones(len(self.weights))

    def has_equal_tableau(self, other_pauli: PauliObject, literal: bool = True) -> bool:
        """
        Check whether two Pauli objects have the same tableau and dimensions.

        Parameters
        ----------
        other_pauli : PauliObject
            Pauli object to compare against.

        literal : bool, optional
            If True, compares objects literally in their current form. If False,
            the objects are first brought to standard form. Default is True.

        Returns
        -------
        bool
            True if all tableau entries and dimensions match; False otherwise.
        """

        ps1 = self
        ps2 = other_pauli
        if not literal:
            ps1 = ps1.to_standard_form()
            ps2 = ps2.to_standard_form()

        if not np.array_equal(ps1.dimensions, ps2.dimensions):
            return False

        if not np.array_equal(ps1.tableau, ps2.tableau):
            return False

        return True

    def is_close(self, other_pauli: Self, threshold: int = 10, literal: bool = True) -> bool:
        """
        Check whether two Pauli objects are approximately equal.

        Parameters
        ----------
        other_pauli : PauliObject
            Pauli object to compare against.
        threshold : int, optional
            Number of matching decimal digits required for equality. Default is 10.
        literal : bool, optional
            If True, compares objects literally in their current form. If False,
            the objects are first brought to standard form. Default is True.

        Returns
        -------
        bool
            True if all tableau entries, weights, phases, and dimensions
            match within tolerance; False otherwise.
        """
        if not isinstance(other_pauli, self.__class__):
            return False

        ps1 = self
        ps2 = other_pauli
        if not literal:
            ps1 = ps1.to_standard_form()
            ps2 = ps2.to_standard_form()

        if not np.all(np.isclose(ps1.weights, ps2.weights, 10**(-threshold))):
            return False

        if not np.array_equal(ps1.phases, ps2.phases):
            return False

        if not np.array_equal(ps1.dimensions, ps2.dimensions):
            return False

        if not np.array_equal(ps1.tableau, ps2.tableau):
            return False

        return True

    def hermitian_conjugate(self) -> Self:
        """
        Return the Hermitian conjugate of the Pauli object.

        The conjugate operation negates tableau exponents, conjugates weights,
        and adjusts phases to preserve physical equivalence.

        Returns
        -------
        PauliObject
            Hermitian conjugate of the Pauli object.
        """
        conjugate_weights = np.conj(self.weights)
        conjugate_tableau = (-self.tableau) % np.tile(self.dimensions, 2)

        acquired_phases = []
        for i in range(self.n_paulis()):
            hermitian_conjugate_phase = 0
            for j in range(self.n_qudits()):
                r = self.tableau[i, j]
                s = self.tableau[i, j + self.n_qudits()]
                hermitian_conjugate_phase += ((r * s) % self.dimensions[j]) * self.lcm / self.dimensions[j]
            acquired_phases.append(2 * hermitian_conjugate_phase)
        acquired_phases = np.asarray(acquired_phases, dtype=int)

        conjugate_initial_phases = (-self.phases) % (2 * self.lcm)
        conjugate_phases = (conjugate_initial_phases + acquired_phases) % (2 * self.lcm)

        return self.__class__(tableau=conjugate_tableau, dimensions=self.dimensions,
                              weights=conjugate_weights, phases=conjugate_phases)

    H = hermitian_conjugate

    def is_hermitian(self) -> bool:
        """
        Check if the Pauli object is Hermitian.

        Returns
        -------
        bool
            True if the object equals its Hermitian conjugate; False otherwise.
        """
        # NOTE: rounding errors could make this fail, hence we call the is_close function.
        return self.is_close(self.H(), literal=False)

    def _sanity_check(self):
        """
        Validate internal consistency of the Pauli object.

        Raises
        ------
        ValueError
            If tableau, dimensions, or exponents are inconsistent or invalid.
        """
        if len(self.weights) != self.n_paulis():
            # FIXME: Improve error message
            raise ValueError("The weights and tableau have inconsistent shapes.")

        if len(self.phases) != self.n_paulis():
            # FIXME: Improve error message
            raise ValueError("The phases and tableau have inconsistent shapes.")

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
        Returns an unambiguous string representation of the PauliObject.

        Returns
        -------
        str
            A string representation of the PauliObject with tableau, dimensions, weights, and phases.
        """
        return f'{self.__class__.__name__}({self.tableau}, {self.dimensions}, {self.weights}, {self.phases})'

    def __eq__(self, other_pauli: Self) -> bool:
        """
        Determine if two Pauli objects are equal.

        Parameters
        ----------
        other_pauli : PauliObject
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

        if not np.array_equal(self.weights, other_pauli.weights):
            return False

        if not np.array_equal(self.phases, other_pauli.phases):
            return False

        if not np.array_equal(self.dimensions, other_pauli.dimensions):
            return False

        return True

    def __ne__(self, other_pauli: PauliObject) -> bool:
        """
        Determine if two Pauli objects are different.

        Parameters
        ----------
        other_pauli : PauliObject
            Object to compare with.

        Returns
        -------
        bool
            True if objects are not equal; False otherwise.
        """
        return not self == other_pauli

    def __gt__(self, other_pauli: PauliObject) -> bool:
        """
        Strict greater-than comparison for ordering single Pauli terms.

        Behavior
        --------
        This operator is intended to impose an ordering on single-term Pauli
        objects by interpreting their tableau as an integer (or comparable)
        representation. It is undefined for multi-term objects.

        Parameters
        ----------
        other_pauli : PauliObject
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
        other_pauli : PauliObject
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

    def __add__(self, A: PauliObject) -> PauliSum:
        """
        Implements the addition of Pauli objects.

        Parameters
        ----------
        A : PauliObject
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
        new_weights = np.concatenate([self.weights, A.weights])
        new_phases = np.concatenate([self.phases, A.phases])

        from .pauli_sum import PauliSum
        return PauliSum(new_tableau, self.dimensions.copy(), new_weights, new_phases)

    def __radd__(self, A: PauliObject) -> PauliSum:
        """
        Implements the addition of Pauli objects.

        Parameters
        ----------
        A : PauliObject
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

    def __sub__(self, A: PauliObject) -> PauliSum:
        """
        Implements the subtraction of Pauli objects.

        Parameters
        ----------
        A : PauliObject
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
        new_weights = np.concatenate([self.weights, -np.array(A.weights)])
        new_phases = np.concatenate([self.phases, A.phases])

        from .pauli_sum import PauliSum
        return PauliSum(new_tableau, self.dimensions.copy(), new_weights, new_phases)

    def __rsub__(self, A: PauliObject) -> PauliSum:
        """
        Implements the subtraction of Pauli objects.

        Parameters
        ----------
        A : PauliObject
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
            Resulting PauliObject after exponentiation.

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
        return self.__class__(tableau, self.dimensions.copy(), self.weights.copy(), self.phases.copy())

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
             tuple(self.weights),
             tuple(self.phases),
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
                'weights': self.weights,
                'phases': self.phases}

    def copy(self) -> Self:
        """
        Creates a copy of the Pauli object.

        Returns
        -------
        Pauli object
            A copy of the Pauli object.
        """
        return self.__class__(self.tableau.copy(), self.dimensions.copy(), self.weights.copy(), self.phases.copy())

    def phase_to_weight(self):
        """
        Include the phases into the weights of the Pauli object.
        This method modifies the weights of the Pauli object by multiplying them with the phases,
        and reset the phases to all zeros.
        """
        new_weights = np.zeros(self.n_paulis(), dtype=np.complex128)
        for i in range(self.n_paulis()):
            phase = self.phases[i]
            omega = np.exp(2 * np.pi * 1j * phase / (2 * self.lcm))
            new_weights[i] = self.weights[i] * omega
        self._phases = np.zeros(self.n_paulis(), dtype=int)
        self._weights = new_weights

    def to_standard_form(self) -> Self:
        """
        Produce a standardized form of the Pauli object.

        The standard form consolidates equivalent Pauli terms, normalizes phases,
        and ensures a canonical ordering of terms.

        Returns
        -------
        Pauli object
            The Pauli object in standard form.
        """
        ps_out = self.copy()
        ps_out.standardise()
        return ps_out

    def standardise(self):
        """
        In-place standardisation of the Pauli object.

        Combines equivalent terms, absorbs phases into weights where appropriate,
        and normalizes the internal representation for deterministic comparisons.
        """
        self.phase_to_weight()
        T = self.tableau
        W = self.weights

        order = np.lexsort(T.T)

        self._tableau = T[order]
        self._weights = W[order]

    standardize = standardise

    @abstractmethod
    def to_hilbert_space(self, pauli_string_index: int | None = None) -> np.ndarray:
        """
        Get the matrix form of the PauliObject as a sparse matrix.

        Parameters
        ----------
        pauli_string_index : int | None, optional
            Index of a specific Pauli term to convert. If None, the full operator
            (sum of all terms) is returned.

        Returns
        -------
        scipy.sparse.csr_matrix
            Matrix representation of input PauliObject.
        """
        pass
