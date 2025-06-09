from __future__ import annotations
import functools
import numpy as np
import re
from .pauli import Pauli
from quaos.core.prime_Functions_Andrew import bases_to_int   # should move somewhere else


@functools.total_ordering
class PauliString:
 
    def __init__(self,
                 x_exp: list[int] | np.ndarray | str,
                 z_exp: list[int] | np.ndarray | None = None,
                 dimensions: list[int] | np.ndarray | None = None):

        if isinstance(x_exp, str):
            if z_exp is not None:
                raise Warning('If input string is provided, z_exp is unnecessary')
            xz_exponents = re.split('x|z', x_exp)[1:]
            z_exp = np.array(xz_exponents[1::2], dtype=int)
            x_exp = np.array(xz_exponents[0::2], dtype=int)
        elif isinstance(x_exp, list):
            z_exp = np.array(z_exp, dtype=int)
            x_exp = np.array(x_exp, dtype=int)
        
        if dimensions is None:
            self.dimensions = (max(max(self.x_exp), max(self.z_exp)) + 1) * np.ones(len(self.x_exp))
        else:
            self.dimensions = np.asarray(dimensions)
        
        self.x_exp = x_exp % self.dimensions
        self.z_exp = z_exp % self.dimensions
        self.lcm = np.lcm.reduce(self.dimensions)
        self._sanity_check()

    def _sanity_check(self):
        if len(self.x_exp) != len(self.dimensions):
            raise ValueError(f"Number of x exponents ({len(self.x_exp)})"
                             f" and dimensions ({len(self.dimensions)}) must be equal.")

        if len(self.x_exp) != len(self.z_exp):
            raise ValueError(f"Number of x and z exponents ({len(self.x_exp)}"
                             f" and {len(self.z_exp)}) must be equal.")
        
        if len(self.dimensions) != len(self.z_exp):
            raise ValueError(f"Number of dimensions ({len(self.dimensions)})"
                             f" and z exponents ({len(self.z_exp)}) must be equal.")

        for i in range(len(self.x_exp)):
            if self.dimensions[i] - 1 < self.x_exp[i] or self.dimensions[i] - 1 < self.z_exp[i]:
                raise ValueError(f"Dimension {self.dimensions[i]} is too small for"
                                 f" exponents {self.x_exp[i]} and {self.z_exp[i]}")
    
    def __repr__(self) -> str:
        return f"Pauli(x_exp={self.x_exp}, z_exp={self.z_exp}, dimensions={self.dimensions})"

    def __str__(self) -> str:
        p_string = ''
        for i in range(self.n_qudits()):
            p_string += 'x' + f'{self.x_exp[i]}' + 'z' + f'{self.z_exp[i]} '
        return p_string
    
    def __matmul__(self, A: PauliString) -> PauliString:
        new_x_exp = np.concatenate((self.x_exp, A.x_exp))
        new_z_exp = np.concatenate((self.z_exp, A.z_exp))
        new_dims = np.concatenate((self.dimensions, A.dimensions))
        return PauliString(new_x_exp, new_z_exp, new_dims)

    def __mul__(self, A: PauliString) -> PauliString:
        if isinstance(A, PauliString):
            if np.any(self.dimensions != A.dimensions):
                raise Exception("To multiply two PauliStrings, their dimensions"
                                f" {self.dimensions} and {A.dimensions} must be equal")
            x_new = np.mod(self.x_exp + A.x_exp, (self.dimensions))
            z_new = np.mod(self.z_exp + A.z_exp, (self.dimensions))
            return PauliString(x_new, z_new, self.dimensions)
        else:
            raise ValueError(f"Cannot multiply PauliString with type {type(A)}")

    def __pow__(self, A: int) -> PauliString:
        return PauliString(self.x_exp * A, self.z_exp * A, self.dimensions)

    # def _to_pauli_sum(self) -> PauliSum:
    #     from . import PauliSum
    #     return PauliSum([self], weights=[1], phases=[0], dimensions=self.dimensions, standardise=False)

    def __eq__(self, other_pauli: 'PauliString') -> bool:
        if not isinstance(other_pauli, PauliString):
            return False
    
        return bool(np.all(self.x_exp == other_pauli.x_exp) and np.all(self.z_exp == other_pauli.z_exp) and np.all(self.dimensions == other_pauli.dimensions))

    def __ne__(self, other_pauli: PauliString) -> bool:
        return not self.__eq__(other_pauli)
    
    def __gt__(self, other_pauli: PauliString) -> bool:
        this_pauli = self._to_int(reverse=True)
        other_pauli = other_pauli._to_int(reverse=True)
        return this_pauli > other_pauli
    
    def _to_int(self, reverse=False):
        dims = self.dimensions
        dims_double = [d for d in dims for _ in range(2)]
        base = np.zeros(len(dims_double), dtype=int)
        base[:len(dims)] = self.x_exp
        base[len(dims):] = self.z_exp
        if not reverse:
            return bases_to_int(base, dims_double)
        else:
            base[:len(dims)] = self.z_exp
            base[len(dims):] = self.x_exp
            return bases_to_int(base[::-1], dims_double[::-1])

    def __hash__(self) -> int:
        return hash((tuple(self.x_exp), tuple(self.z_exp), tuple(self.dimensions)))
    
    def __dict__(self) -> dict:
        return {'x_exp': self.x_exp, 'z_exp': self.z_exp, 'dimensions': self.dimensions}

    def n_qudits(self) -> int:
        return len(self.x_exp)
        
    def n_identities(self) -> int:
        """
        Get the number of identities in the PauliString
        :return: The number of identities
        """
        return np.sum(np.logical_and(self.x_exp == 0, self.z_exp == 0))
    
    def get_paulis(self) -> list[Pauli]:
        """
        Get a list of Pauli objects from the PauliString
        :return: A list of Pauli objects
        """
        return [Pauli(x_exp=self.x_exp[i], z_exp=self.z_exp[i], dimension=self.dimensions[i]) for i in range(len(self.x_exp))]

    def symplectic(self) -> np.ndarray:
        symp = np.zeros(2 * self.n_qudits())
        symp[0:self.n_qudits()] = self.x_exp
        symp[self.n_qudits():2 * self.n_qudits()] = self.z_exp
        return symp

    def symplectic_product(self, A: PauliString) -> np.ndarray:
        n = self.n_qudits()
        symp = self.symplectic()
        symp_A = A.symplectic()
        prod = sum([symp[i] * symp_A[i + n] - symp[i + n] * symp_A[i] for i in range(n)]) % self.lcm
        return prod
    
    def amend(self, qudit_index: int, new_x: int, new_z: int) -> PauliString:
        if new_x > self.dimensions[qudit_index] or new_z > self.dimensions[qudit_index]:
            raise ValueError(f"Exponents ({new_x, new_z}) cannot be larger than qudit dimension"
                             f" ({self.dimensions[qudit_index]})")
        self.x_exp[qudit_index] = new_x
        self.z_exp[qudit_index] = new_z
        return self
    
    def acquired_phase(self, other_pauli: 'PauliString') -> int:
        # phases acquired when multiplying two Pauli strings
        phi = 1.  # / self.dimensions
        phase = phi * (np.dot(self.x_exp, other_pauli.z_exp) - np.dot(self.z_exp, other_pauli.x_exp))
        return phase % self.lcm
    
    def _replace_symplectic(self, symplectic: np.ndarray, qudit_indices: list[int]) -> PauliString:
        x_exp_replace = symplectic[0:len(qudit_indices)]
        z_exp_replace = symplectic[len(qudit_indices):2 * len(qudit_indices)]

        x_exp = self.x_exp.copy()
        z_exp = self.z_exp.copy()
        for i, index in enumerate(qudit_indices):
            x_exp[index] = x_exp_replace[i]
            z_exp[index] = z_exp_replace[i]

        return PauliString(x_exp=x_exp, z_exp=z_exp, dimensions=self.dimensions)
    
    def _delete_qudits(self, qudit_indices: list[int], return_new=True) -> PauliString:  # not sure if here it is best to return a new object or not
        x_exp = np.delete(self.x_exp, qudit_indices)
        z_exp = np.delete(self.z_exp, qudit_indices)
        dimensions = np.delete(self.dimensions, qudit_indices)
        if return_new:
            return PauliString(x_exp=x_exp, z_exp=z_exp, dimensions=dimensions)
        else:
            self.x_exp = x_exp
            self.z_exp = z_exp
            self.dimensions = dimensions
            self._sanity_check()
            return self
        
    def __getitem__(self, key: int | slice | list[int]) -> PauliString | Pauli:
        if isinstance(key, int):
            return self.get_paulis()[key]
        else:
            return PauliString(x_exp=self.x_exp[key], z_exp=self.z_exp[key], dimensions=self.dimensions[key])
        
    def __setitem__(self, key: int | slice, value: Pauli | PauliString):
        if isinstance(key, int) and isinstance(value, Pauli):
            self.x_exp[key] = value.x_exp
            self.z_exp[key] = value.z_exp
            self.dimensions[key] = value.dimension
        elif isinstance(value, PauliString):
            self.x_exp[key] = value.x_exp
            self.z_exp[key] = value.z_exp
            self.dimensions[key] = value.dimensions
        else:
            raise Exception(f"Cannot set PauliString with type {type(value)}")

    def get_subspace(self, qudit_indices: list[int] | int) -> PauliString:
        return PauliString(x_exp=self.x_exp[qudit_indices], z_exp=self.z_exp[qudit_indices],
                           dimensions=self.dimensions[qudit_indices])
    
    def copy(self) -> PauliString:
        return PauliString(x_exp=self.x_exp.copy(), z_exp=self.z_exp.copy(), dimensions=self.dimensions.copy())
    
    def commute(self, other_pauli: PauliString) -> bool:
        """
        Check if two Pauli strings commute
        :param other_pauli: The other Pauli string
        :return: True if they commute, False otherwise
        """
        return self.symplectic_product(other_pauli) == 0
    
    def hermitian(self) -> PauliString:
        """
        Return the Hermitian conjugate of the PauliString.

        The Hermitian conjugate is obtained by negating the exponents
        of the Pauli operations and taking the modulus with respect to
        their dimensions.

        :return: A new PauliString representing the Hermitian conjugate.
        """
        return PauliString(x_exp=(-self.x_exp)%self.dimensions, z_exp=(-self.z_exp)%self.dimensions, dimensions=self.dimensions)
    

        

        
