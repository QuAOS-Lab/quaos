from __future__ import annotations
from typing import overload, TYPE_CHECKING
import numpy as np
import re
import scipy.sparse as sp

from .pauli_object import PauliObject
from .pauli import Pauli
from .constants import DEFAULT_QUDIT_DIMENSION
from .bases_manipulation import bases_to_int

if TYPE_CHECKING:
    from .pauli_sum import PauliSum


class PauliString(PauliObject):
    @classmethod
    def from_exponents(cls, x_exp: list[int] | np.ndarray | str | int, z_exp: list[int] | np.ndarray | int,
                       dimensions: list[int] | np.ndarray | int | None = None) -> PauliString:
        """
        Create a PauliString instance from X and Z exponents.

        Parameters
        ----------
        x_exp : list[int] | np.ndarray | str | int
            The X exponents used to create the PauliString tableau.
            If no Z exponents are given, the X exponents can be also given as a parseable string.
            If an integer is provided, it is assumed the PauliString contains a single Pauli.
        z_exp: list[int] | np.ndarray | int
            The Z exponents used to create the PauliString tableau.
            If an integer is provided, it is assumed the PauliString contains a single Pauli.
        dimensions: list[int] | np.ndarray | int, optional
            The dimensions of each qudit. If an integer is provided,
            all qudits are assumed to have the same dimension.
            If no value is provided, the default is `DEFAULT_QUDIT_DIMENSION`.

        Returns
        -------
        PauliString
            A PauliString instance representing the given exponents.

        Raises
        ------
        ValueError
            If x_exp, z_exp, and dimensions don't have the same length.
        """
        if isinstance(x_exp, str):
            if z_exp is not None:
                raise Warning('If input string is provided, z_exp is unnecessary')
            xz_exponents = re.split('x|z', x_exp)[1:]
            z_exp = np.array(xz_exponents[1::2], dtype=int)
            x_exp = np.array(xz_exponents[0::2], dtype=int)
        elif isinstance(x_exp, list):
            x_exp = np.array(x_exp)
        elif isinstance(x_exp, int):
            x_exp = np.array([x_exp], dtype=int)

        if isinstance(z_exp, list):
            z_exp = np.array(z_exp)
        elif isinstance(z_exp, int):
            z_exp = np.array([z_exp], dtype=int)

        if dimensions is None:
            dimensions = DEFAULT_QUDIT_DIMENSION * np.ones(len(x_exp), dtype=int)
        elif type(dimensions) is int:
            dimensions = dimensions * np.ones(len(x_exp), dtype=int)
        else:
            dimensions = np.asarray(dimensions)

        if len(x_exp) != len(z_exp):
            raise ValueError(f"Number of x and z exponents ({len(x_exp)}"
                             f" and {len(z_exp)}) must be equal.")

        if len(x_exp) != len(dimensions):
            raise ValueError(f"Number of x exponents ({len(x_exp)})"
                             f" and dimensions ({len(dimensions)}) must be equal.")

        tableau = np.empty(2 * len(dimensions), dtype=int)
        tableau[:len(dimensions)] = x_exp % dimensions
        tableau[len(dimensions):] = z_exp % dimensions

        P = cls(tableau, dimensions)
        P._sanity_check()

        return P

    @classmethod
    def from_pauli(cls, pauli: Pauli) -> PauliString:
        """
        Create a PauliString instance from a single Pauli.

        Parameters
        ----------
        pauli : Pauli
            The Pauli to convert into a PauliString.

        Returns
        -------
        PauliString
            A PauliString instance representing the given Pauli operator.
        """

        P = cls(pauli.tableau, pauli.dimensions)
        P._sanity_check()

        return P

    @classmethod
    def from_tableau(
        cls,
        tableau: list[int] | np.ndarray,
        dimensions: int | list[int] | np.ndarray | None = None
    ) -> PauliString:
        """
        Create a PauliString instance from a tableau.

        Parameters
        ----------
        tableau : list[int] | np.ndarray
            The tableau to convert into a PauliString.
        dimensions : list[int] | np.ndarray | int, optional
            The dimensions of each qudit. If an integer is provided,
            all qudits are assumed to have the same dimension.
            If no value is provided, the default is `DEFAULT_QUDIT_DIMENSION`.

        Returns
        -------
        PauliString
            A PauliString instance initialized with the exponents and dimensions from the input tableau.
        """

        if isinstance(tableau, list):
            tableau = np.asarray(tableau, dtype=int)
        P = cls(tableau, dimensions)
        P._sanity_check()

        return P

    @classmethod
    def from_string(cls, pauli_str: str, dimensions: int | list[int] | np.ndarray) -> PauliString:
        """
        Create a PauliString instance from a string representation.

        Parameters
        ----------
        pauli_str : str
            The string representation of the Pauli string, where exponents are separated by 'x' and 'z'.
        dimensions : list[int] | np.ndarray
            The dimensions parameter to be passed to the PauliString constructor.

        Returns
        -------
        PauliString
            An instance of PauliString initialized with the exponents and dimensions parsed from the input string.

        Examples
        --------
        >>> PauliString.from_string("x2z3 x4z1 x0z0", [4, 5, 2]])
        <PauliString ...>
        """

        xz_exponents = re.split('x|z', pauli_str)[1:]
        z_exp = np.array(xz_exponents[1::2], dtype=int)
        x_exp = np.array(xz_exponents[0::2], dtype=int)

        tableau = np.concatenate([x_exp, z_exp])
        P = cls(tableau, dimensions)
        P._sanity_check()

        return P

    @classmethod
    def from_random(cls, dimensions: int | list[int] | np.ndarray, seed: int | None = None) -> PauliString:
        """
        Generate a random PauliString instance for a given number of qudits and their dimensions.

        Parameters
        ----------
        n_qudits : int
            The number of qudits in the Pauli string.
        dimensions : int | list[int] | np.ndarray
            The dimensions of each qudit. Should be a list or array of integers specifying the dimension for each qudit.
            The size of dimensions determines the number of qudits.
        seed : int or None, optional
            Seed for the random number generator to ensure reproducibility. Default is None.

        Returns
        -------
        PauliString
            A randomly generated PauliString instance with random x and z exponents for each qudit.
        """
        if seed is not None:
            np.random.seed(seed)

        if isinstance(dimensions, int):
            dimensions = [dimensions]

        tableau = np.concatenate([np.random.randint(dimensions, dtype=int), np.random.randint(dimensions, dtype=int)])
        return cls(tableau, dimensions)

    def __repr__(self) -> str:
        """
        Return the string representation of the PauliString.
        (in a format that is helpful for debugging).

        Returns
        -------
        str
            A string in the format "PauliString(tableau=..., dimensions=...)".
        """

        return f"PauliString(tableau={self.tableau}, dimensions={self.dimensions})"

    def __str__(self) -> str:
        """
        Return the string representation of the PauliString.
        (in the format also employed for constructing an instance of the class).

        Returns
        -------
        str
            The string representation of the PauliString, with each qudit's
            exponents shown as 'x{exp}z{exp} x{exp}z{exp} x{exp}z{exp} ...'.
        """

        p_string = ''
        for i in range(self.n_qudits()):
            p_string += 'x' + f'{self.x_exp[i]}' + 'z' + f'{self.z_exp[i]} '
        return p_string

    def __matmul__(self, A: PauliString) -> PauliString:
        """
        Implements the tensor product between PauliString objects.
        It corresponds to operator tensor product (`@`).
        The resulting PauliString has the exponents of both strings concatenated.

        Parameters
        ----------
        A : PauliString
            The PauliString instance to be tensored with `self`.

        Returns
        -------
        PauliString
            A new PauliString instance representing the tensor product of `self` and `A`.

        Examples
        --------
        >>> p1 = PauliString.from_string("x1z0 x0z1", [2, 2])
        >>> p2 = PauliString.from_string("x2z1", [3])
        >>> p1 @ p2
        Pauli(x_exp=[1 0 2], z_exp=[0 1 1], dimensions=[2 2 3])

        Notes
        -----
        - This is NOT a product between two PauliString objects!
        """

        new_n_qudits = self.n_qudits() + A.n_qudits()
        tableau = np.empty(2 * new_n_qudits, dtype=int)
        tableau[:self.n_qudits()] = self.x_exp
        tableau[self.n_qudits():new_n_qudits] = A.x_exp
        tableau[new_n_qudits:new_n_qudits + self.n_qudits()] = self.z_exp
        tableau[new_n_qudits + self.n_qudits():] = A.z_exp

        dimensions = np.empty(new_n_qudits, dtype=int)
        dimensions[:self.n_qudits()] = self.dimensions
        dimensions[self.n_qudits():] = A.dimensions
        return PauliString(tableau, dimensions)

    def __rmatmul__(self, A: Pauli) -> PauliString:
        """
        Implements the tensor product between a PauliString (`self`) and a Pauli (`A`) objects.
        It corresponds to operator tensor product (`@`).
        The resulting PauliString has the exponents of both strings concatenated.

        Parameters
        ----------
        A : Pauli
            The Pauli to be right-multiplied with this PauliString.

        Returns
        -------
        PauliString
            A new PauliString instance with updated `x_exp`, `z_exp`, and `dimensions` arrays,
            resulting from concatenating the corresponding attributes of the operands.

        Examples
        --------
        >>> ps = PauliString(...)
        >>> p = Pauli(...)
        >>> result = p @ ps
        """

        new_n_qudits = self.n_qudits() + A.n_qudits()
        tableau = np.empty(2 * new_n_qudits, dtype=int)
        tableau[:self.n_qudits()] = self.x_exp
        tableau[self.n_qudits():new_n_qudits] = A.x_exp
        tableau[new_n_qudits:new_n_qudits + self.n_qudits()] = self.z_exp
        tableau[new_n_qudits + self.n_qudits():] = A.z_exp

        dimensions = np.empty(new_n_qudits, dtype=int)
        dimensions[:self.n_qudits()] = self.dimensions
        dimensions[self.n_qudits():] = A.dimensions
        return PauliString(tableau, dimensions)

    def __mul__(self, A: PauliString) -> PauliString:
        """
        Multiply two PauliString objects element-wise. It corresponds to operator multiplication (`*`).
        It adds the tableaus of the two PauliStrings modulo their dimensions.

        Parameters
        ----------
        A : PauliString
            The PauliString to multiply with `self`.

        Returns
        -------
        PauliString
            A new PauliString representing the product of `self` and `A`.

        Raises
        ------
        Exception
            If the dimensions of the two PauliStrings are not equal.
        ValueError
            If `A` is not an instance of PauliString.

        Examples
        --------
        >>> ps1 = PauliString.from_exponents(x_exp1, z_exp1, dims)
        >>> ps2 = PauliString.from_exponents(x_exp2, z_exp2, dims)
        >>> ps3 = ps1 * ps2
        """
        if not isinstance(A, PauliString):
            raise ValueError(f"Cannot multiply PauliString with type {type(A)}")

        if not np.array_equal(self.dimensions, A.dimensions):
            raise Exception("To multiply two PauliStrings, their dimensions"
                            f" {self.dimensions} and {A.dimensions} must be equal.")

        tableau = np.mod(self.tableau + A.tableau, np.tile(self.dimensions, 2))

        w1 = self._weights[:, None]
        w2 = A._weights[None, :]
        new_weights = (w1 * w2).reshape(-1)

        p1 = self._phases[:, None]
        p2 = A._phases[None, :]

        # Extract z- and x-parts from tableau
        n1, n2 = self.n_qudits(), A.n_qudits()
        a_z = self.tableau[:, n1:]
        b_x = A.tableau[:, :n2]
        # Compute acquired phases via symplectic form
        factors = (self.lcm // self.dimensions)
        acquired_phases = 2 * factors * a_z  @ b_x.T

        # Combine with existing phases and flatten
        new_phases = (p1 + p2 + acquired_phases) % (2 * self.lcm)
        new_phases = new_phases.reshape(-1)
        return PauliString(tableau, self.dimensions.copy(), weights=new_weights, phases=new_phases)

    @property
    def x_exp(self) -> np.ndarray:
        """
        x_exp : np.ndarray
        Array of X exponents for each qudit.
        """
        return self._tableau[0][:self.n_qudits()]

    @property
    def z_exp(self) -> np.ndarray:
        """
        z_exp : np.ndarray
        Array of Z exponents for each qudit.
        """
        return self._tableau[0][self.n_qudits():]

    @property
    def phases(self) -> np.ndarray:
        """
        Returns the phases associated with the PauliString.

        Returns
        -------
        np.ndarray
            The phases as a 1d-vector.
        """
        return np.asarray([0], dtype=int)

    @property
    def weights(self) -> np.ndarray:
        """
        Returns the weights associated with the Pauli String.

        Returns
        -------
        np.ndarray
            The weights as a 1d-vector.
        """
        return np.asarray([1], dtype=complex)

    def as_pauli_sum(self) -> PauliSum:
        """
        Converts the PauliString to the embedding PauliSum.

        Returns
        -------
        PauliSum
            The PauliSum containing the PauliString as single entry.
        """
        # FIXME: import at the top. Currently we can't because of circular imports.
        from .pauli_sum import PauliSum
        return PauliSum(self._tableau, self._dimensions, self._weights, self._phases)

    def n_identities(self) -> int:
        """
        Returns the number of identities within this PauliString.

        Returns
        -------
        int
            The number of identities.
        """
        return np.sum(np.logical_and(self.x_exp == 0, self.z_exp == 0))

    def get_pauli(self, index: int) -> Pauli:
        """
        Returns a Pauli at the input index.

        Parameters
        ----------
        index : int
            The index of the Pauli.

        Returns
        -------
        Pauli
            The Pauli at the given index.
        """

        tableau = np.asarray([self.x_exp[index], self.z_exp[index]], dtype=int)
        return Pauli(tableau, int(self.dimensions[index]))

    def get_paulis(self) -> list[Pauli]:
        """
        Returns a list of Paulis corresponding to the PauliString tableau.

        Returns
        -------
        list[Pauli]
            A list of Pauli.
        """
        return [self.get_pauli(i) for i in range(self.n_qudits())]

    def symplectic_residues(self, A: PauliString) -> np.ndarray:
        """
        Per-qudit symplectic residues r_j = x_j z'_j - z_j x'_j  (mod d_j).
        Returns a length-n int array with entry-wise mod d_j.
        """
        if not np.array_equal(self.dimensions, A.dimensions):
            raise ValueError(
                f"Incompatible PauliStrings: must identical dimensions "
                f"(currently {self.dimensions} and {A.dimensions})."
            )

        dims = self.dimensions

        # Per-site residue, reduced mod the local dimension
        residues = (self.x_exp * A.z_exp - self.z_exp * A.x_exp) % dims
        return residues.astype(int)

    def symplectic_product(self, A: PauliString, *, as_scalar: bool = True) -> np.ndarray | int:
        """
        Compute the symplectic product between this PauliString and another.
        The symplectic product is defined as the sum over all qudits of the difference between
        the product of the X component of this PauliString and the Z component of the other,
        and the product of the Z component of this PauliString and the X component of the other,
        modulo the least common multiple (lcm) of the qudit dimensions. See
        `Phys. Rev. A 70, 052328 (2004) <https://doi.org/10.1103/PhysRevA.70.052328>`_
        for more details.

        We generalise to mixed dimensions as follows:

        If as_scalar=True (default): return a single integer r (mod LCM),
        defined by the phase-preserving lift:
        r = sum_j (L / d_j) * r_j   (mod L),   L = lcm(d_j),  r_j as in symplectic_residues().

        If as_scalar=False: return the per-qudit residue vector r_j (mod d_j).

        Returns
        -------
        int or np.ndarray
            - int in [0, L-1] when as_scalar=True
            - length-n array of residues when as_scalar=False
        """
        residues = self.symplectic_residues(A)           # length-n residues mod d_j
        if not as_scalar:
            return residues

        L = self.lcm
        dims = self.dimensions
        weights = (L // dims).astype(int)
        r = int(np.sum((weights * residues) % L) % L)    # this is the number for acquired phase
        return r

    def quditwise_product(self, A: PauliString) -> int:
        """
        Compute the quditwise product between this PauliString and another.
        The quditwise product is defined as the elementwise product of the X and Z exponents.

        Parameters
        ----------
        A : PauliString
            The other PauliString to compute the quditwise product with.

        Returns
        -------
        np.ndarray
            The quditwise product as a NumPy array.
        """
        if self.n_qudits() != A.n_qudits() or not np.array_equal(self.dimensions, A.dimensions):
            raise ValueError(
                (
                    f"Incompatible PauliStrings: must have the same number of qudits "
                    f"(currently {self.n_qudits()} and {A.n_qudits()}) and dimensions "
                    f"(currently {self.dimensions} and {A.dimensions})."
                )
            )
        if any(
            np.sum(self.x_exp[i] * A.z_exp[i] - self.z_exp[i] * A.x_exp[i]) % self.dimensions[i]
            for i in range(self.n_qudits())
        ):
            return 1

        return 0

    def amend(self, qudit_index: int, new_x: int, new_z: int) -> PauliString:
        """
        Change the X and Z exponents for a specific qudit in the Pauli string.

        Parameters
        ----------
        qudit_index : int
            The index of the qudit to change.
        new_x : int
            The new exponent for the X operator at the specified qudit.
        new_z : int
            The new exponent for the Z operator at the specified qudit.

        Returns
        -------
        PauliString
            The modified PauliString instance.

        Raises
        ------
        ValueError
            If either `new_x` or `new_z` is greater than the dimension of the specified qudit.
        """
        if new_x > self.dimensions[qudit_index] or new_z > self.dimensions[qudit_index]:
            raise ValueError(f"Exponents ({new_x, new_z}) cannot be larger than qudit dimension"
                             f" ({self.dimensions[qudit_index]})")
        self._tableau[0][qudit_index] = new_x
        self._tableau[0][self.n_qudits() + qudit_index] = new_z
        return self

    def acquired_phase(self, other_pauli: PauliString) -> int:
        """
        Compute the phase acquired when multiplying two Pauli strings self * other_pauli.
        This method calculates the phase factor resulting from the multiplication of two Pauli strings,
        using their symplectic representations. The phase is computed modulo twice the least common multiple (lcm)
        of the underlying dimensions. For details, see:
        `IEEE International Symposium on Information Theory (ISIT), pp. 791-795. IEEE (2018)
        <https://doi.org/10.1109/ISIT.2018.8437652>`_

        Parameters
        ----------
        other_pauli : PauliString
            The Pauli string to multiply with this Pauli string.

        Returns
        -------
        int
            The acquired phase, as an integer modulo ``2 * self.lcm``.
        """

        if not np.array_equal(self.dimensions, other_pauli.dimensions):
            raise Exception("To acquire phase from another PauliString, their dimensions"
                            f" {self.dimensions} and {other_pauli.dimensions} must be equal")

        n = self.n_qudits()
        a = self.tableau[0]
        b = other_pauli.tableau[0]

        # U is zeros with identity in lower-left n x n block
        # This is equivalent to sum over j of 2 * x'_j * z_j
        # U @ b selects b[:n] (x part) and puts it in lower half
        phase = 2 * np.dot(a[n:], b[:n])  # THIS ASSUMES [1 | 1] is XZ, NOT Y

        return int(phase % (2 * self.lcm))

    # TODO: not sure if here it is best to return a new object or not
    def _delete_qudits(self, mask: np.ndarray, return_new: bool = True) -> PauliString:
        """
        Delete specified qudits from the PauliString.
        Removes the qudits at the given indices from the internal representations
        (`x_exp`, `z_exp`, and `dimensions`). Optionally returns a new PauliString
        instance with the specified qudits removed, or modifies the current instance
        in place.

        Parameters
        ----------
        mask : np.ndarray
            A boolean mask where the value at the indices of the qudits to be deleted is False.
        return_new : bool, optional
            If True (default), returns a new PauliString instance with the specified
            qudits removed. If False, modifies the current instance in place and
            returns self.

        Returns
        -------
        PauliString
            A new PauliString instance with the specified qudits removed if
            `return_new` is True, otherwise returns self after modification.
        """
        dimensions = self.dimensions[mask]
        new_tableau = np.empty(2 * self.n_qudits(), dtype=int)
        new_tableau[:self.n_qudits()] = self.x_exp[mask]
        new_tableau[self.n_qudits():] = self.z_exp[mask]

        if return_new:
            return PauliString(new_tableau, dimensions)

        self._tableau = new_tableau
        self._dimensions = dimensions
        self._lcm = np.lcm.reduce(dimensions)
        return self

    @overload
    def __getitem__(self, key: int) -> Pauli:
        ...

    @overload
    def __getitem__(self, key: slice | np.ndarray | list[int]) -> PauliString:
        ...

    def __getitem__(self, key: int | slice | np.ndarray | list[int]) -> PauliString | Pauli:
        """
        Retrieve a Pauli or a (smaller) PauliString from the PauliString.

        Parameters
        ----------
        key : int | slice | np.ndarray | list[int]
            The index or indices specifying which Pauli(s) to retrieve. If an int, returns a single Pauli.
            If a slice, numpy array, or list, returns a new PauliString containing the selected Paulis.

        Returns
        -------
        PauliString or Pauli
            The selected Pauli operator(s). Returns a single Pauli if `key` is an int, otherwise returns a PauliString.

        Raises
        ------
        ValueError
            If `key` is not an int, slice, numpy.ndarray, or list.

        Examples
        --------
        >>> ps = PauliString(...)
        >>> ps[0]  # Returns a single Pauli
        >>> ps[1:3]  # Returns a PauliString with selected Paulis
        >>> ps[[0, 2]]  # Returns a PauliString with Paulis at indices 0 and 2
        """

        # Return a single Pauli
        if isinstance(key, int):
            return self.get_pauli(key)

        # Return a (smaller) PauliString
        if isinstance(key, slice) or isinstance(key, list):
            key = np.asarray(key, dtype=int)

        if isinstance(key, np.ndarray):
            tableau_mask = np.concatenate([key, key + self.n_qudits()])
            return PauliString(
                self.tableau[tableau_mask], self.dimensions[key], self._weights, self._phases)

        raise ValueError(f"Cannot get item with key {key}. Key must be aof type int, slice, np.ndarray, or list[int].")

    def __setitem__(self, key: int | slice | np.ndarray | list[int], value: Pauli | PauliString):
        """
        Set the value(s) of the PauliString at the specified index or slice.

        Parameters
        ----------
        key : int | slice | np.ndarray | list[int]
            The index or slice at which to set the value. If an integer, a single Pauli is set.
            If a slice, a PauliString is set over the specified range.
        value : Pauli | PauliString
            The value to set at the specified key. Must be a Pauli if key is an int,
            or a PauliString if key is a slice.

        Raises
        ------
        ValueError
            If the key and value types do not match the expected combinations.
        """
        # FIXME: merge this with PauliSum __setitem__ into PauliObject
        # TODO: is it necessary to distinguish the two cases in the if... elif... loop?

        if isinstance(key, int) and isinstance(value, Pauli):
            self._tableau[key] = value.x_exp
            self._tableau[key + self.n_qudits()] = value.z_exp
            self._dimensions[key] = value.dimensions[0]
            self._lcm = np.lcm.reduce(self.dimensions)

        elif isinstance(value, PauliString):
            if isinstance(key, slice):
                # Trick to convert slice to NumPy array.
                # This is necessary to be able to get the number of items in the slice.
                key = np.asarray(range(key.stop)[key], dtype=int)
            elif isinstance(key, list):
                key = np.asarray(key, dtype=int)
            elif not isinstance(key, np.ndarray):
                raise ValueError(f"Cannot set item with key {key} and value {value}.\
                                 Invalid key type.")

            if len(key) != value.n_qudits():
                raise ValueError(f"Cannot set item with key {key} and value {value}:\
                                 mismatching dimensions.")
            self._tableau[key] = value.x_exp
            self._tableau[key + self.n_qudits()] = value.z_exp
            self._dimensions[key] = value.dimensions
            self._lcm = np.lcm.reduce(self.dimensions)

        else:
            raise ValueError(f"Cannot set item with key {key} and value {value}.")

    def get_subspace(self, qudit_indices: int | list[int] | np.ndarray) -> PauliString:
        """
        Extracts a subspace of the PauliString corresponding to the specified qudit indices.

        Parameters
        ----------
        qudit_indices : list[int] or int
            The indices of the qudits to extract. Can be a single integer or a list of integers.

        Returns
        -------
        PauliString
            A new PauliString object representing the subspace defined by the selected qudit indices.

        Examples
        --------
        >>> ps = PauliString(...)
        >>> sub_ps = ps.get_subspace([0, 2])
        """

        dimensions = self.dimensions[qudit_indices]
        tableau = self._tableau[0][qudit_indices]

        return PauliString(tableau, dimensions)

    def commute(self, other_pauli: PauliString) -> bool:
        """
        Determine whether this Pauli string commutes with another Pauli string.
        It is efficiently checked via the symplectic product.

        Parameters
        ----------
        other_pauli : PauliString
            The other Pauli string to check commutation with.

        Returns
        -------
        bool
            True if the two Pauli strings commute, False otherwise.

        Notes
        -----
        Two Pauli strings commute if their symplectic product is zero.
        """
        return self.symplectic_product(other_pauli) == 0

    def is_identity(self) -> bool:
        """
        Check if the PauliString represents the identity operator.

        Returns
        -------
        bool
            True if the PauliString is the identity operator, False otherwise.
        """
        return bool(np.all(self.tableau == 0))

    def to_hilbert_space(self) -> sp.csr_matrix:
        """
        Get the matrix form of the PauliString as a sparse matrix.

        Returns
        -------
        scipy.sparse.csr_matrix
            Matrix representation of input Pauli.
        """
        return self.as_pauli_sum().to_hilbert_space()

    def _to_int(self, reverse=False):
        """
        Convert the Pauli string exponents to an integer representation.

        Parameters
        ----------
        reverse : bool, optional
            If True, the order of X and Z exponents is swapped and the resulting
            base and dimensions arrays are reversed before conversion. Default is False.

        Returns
        -------
        int
            The integer representation of the Pauli string.

        Notes
        -----
        - The conversion uses the `bases_to_int` function to map the base array. Look there for an example.
        - IMPORTANT: if too many qudits are provided, the output could be exceedingly large.
        """

        n = self.n_qudits()
        if n > 15:
            raise ValueError(f"Cannot convert PauliString with more than {self.n_qudits()} qudits to int, "
                             "as it may exceed the maximum integer size. Current max set to 15 qudits.")
        dims_double = np.tile(self.dimensions, 2)
        base = np.empty(2 * self.n_qudits(), dtype=int)

        if not reverse:
            base[:n] = self.x_exp
            base[n:] = self.z_exp
            return bases_to_int(base, dims_double)

        base[:n] = self.z_exp
        base[n:] = self.x_exp
        return bases_to_int(base[::-1], dims_double[::-1])
