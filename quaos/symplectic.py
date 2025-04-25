from typing import Any
import numpy as np
import re
import functools
import scipy


class Pauli:
    def __init__(self,
                 x_exp: int | str,
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
        if isinstance(x_exp, str):
            if z_exp is not None:
                raise Warning('If input string is provided, z_exp is unnecessary')
            z_exp = int(x_exp[3])
            x_exp = int(x_exp[1])
        else:
            if (type(x_exp) is not int and type(x_exp) is not np.int64) or (type(z_exp) is not int and type(z_exp) is not np.int64):
                raise TypeError("x_exp and z_exp must be integers or x_exp must be a string of format 'xrzs'")

        self.x_exp = x_exp
        self.z_exp = z_exp
        self.dimension = dimension

        if self.dimension - 1 < x_exp or self.dimension - 1 < z_exp:
            raise ValueError(f"Dimension {self.dimension} is too small for exponents {self.x_exp} and {self.z_exp}")

    def __mul__(self, A: 'str | Pauli | float') -> 'Pauli | PauliSum':
        if isinstance(A, str):
            return self * Pauli(A)
        elif isinstance(A, Pauli):
            if A.dimension != self.dimension:
                raise Exception("To multiply two Paulis, their dimensions"
                                f" {A.dimension} and {self.dimension} must be equal")
            
            return Pauli(x_exp=(self.x_exp + A.x_exp) % self.dimension,
                         z_exp=(self.z_exp + A.z_exp) % self.dimension,
                         dimension=self.dimension)
        elif isinstance(A, float):
            return PauliSum(self, weights=A)
        else:
            raise Exception(f"Cannot multiply Pauli with type {type(A)}")
    
    def __str__(self) -> str:
        return f'x{self.x_exp}z{self.z_exp}'
    
    def __matmul__(self, A: 'Pauli') -> 'PauliString':
        return PauliString(x_exp=[self.x_exp] + [A.x_exp], z_exp=[self.z_exp] + [A.z_exp],
                           dimensions=[self.dimension] + [A.dimension])

    def __add__(self, A: 'Pauli | PauliString | PauliSum') -> 'PauliSum':
        ps1 = PauliString(x_exp=[self.x_exp], z_exp=[self.z_exp], dimensions=[self.dimension])
        if isinstance(A, Pauli):
            ps2 = PauliString(x_exp=[A.x_exp], z_exp=[A.z_exp], dimensions=[A.dimension])
        elif isinstance(A, PauliString) or isinstance(A, PauliSum):
            ps2 = A
        else:
            raise Exception(f"Cannot add Pauli with type {type(A)}")
        return ps1 + ps2
    
    def __sub__(self, A: 'Pauli | PauliString | PauliSum') -> 'PauliSum':
        ps1 = PauliString(x_exp=[self.x_exp], z_exp=[self.z_exp], dimensions=[self.dimension])
        ps1 = PauliSum([ps1])
        if isinstance(A, Pauli):
            ps2 = PauliString(x_exp=[A.x_exp], z_exp=[A.z_exp], dimensions=[A.dimension])
            ps2 = ps2._to_pauli_sum()
        elif isinstance(A, PauliString):
            ps2 = A._to_pauli_sum()
        elif isinstance(A, PauliSum):
            ps2 = A
        else:
            raise Exception(f"Cannot add Pauli with type {type(A)}")
        return ps1 - ps2

    def __eq__(self, other_pauli: Any) -> bool:
        if not isinstance(other_pauli, Pauli):
            return False
        return self.x_exp == other_pauli.x_exp and self.z_exp == other_pauli.z_exp and self.dimension == other_pauli.dimension
    
    def __ne__(self, other_pauli: Any) -> bool:
        return not self.__eq__(other_pauli)
    
    def __dict__(self) -> dict:
        return {'x_exp': self.x_exp, 'z_exp': self.z_exp, 'dimension': self.dimension}
    
    def _to_pauli_string(self) -> 'PauliString':
        return PauliString(x_exp=[self.x_exp], z_exp=[self.z_exp], dimensions=[self.dimension])
    
    def _to_pauli_sum(self) -> 'PauliSum':
        return PauliSum([self._to_pauli_string()], standardise=False)
    
    def __gt__(self, other_pauli: 'Pauli') -> bool:
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
    
    def copy(self) -> 'Pauli':
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
            z_exp = np.array(z_exp)
            x_exp = np.array(x_exp)
        
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
    
    def __matmul__(self, A: 'PauliString') -> 'PauliString':
        new_x_exp = np.concatenate((self.x_exp, A.x_exp))
        new_z_exp = np.concatenate((self.z_exp, A.z_exp))
        new_dims = np.concatenate((self.dimensions, A.dimensions))
        return PauliString(new_x_exp, new_z_exp, new_dims)

    def __mul__(self, A: 'PauliString | PauliSum | float | int | complex') -> 'PauliString | PauliSum':
        if isinstance(A, PauliString):
            if np.any(self.dimensions != A.dimensions):
                raise Exception("To multiply two PauliStrings, their dimensions"
                                f" {self.dimensions} and {A.dimensions} must be equal")
            x_new = np.mod(self.x_exp + A.x_exp, (self.dimensions))
            z_new = np.mod(self.z_exp + A.z_exp, (self.dimensions))
            return PauliString(x_new, z_new, self.dimensions)
        elif isinstance(A, PauliSum):
            return self._to_pauli_sum() * A
        elif isinstance(A, float) or isinstance(A, int) or isinstance(A, complex):
            return PauliSum([self], weights=[A])
        else:
            raise ValueError(f"Cannot multiply PauliString with type {type(A)}")
        
    def __rmul__(self, A: 'PauliString | PauliSum | float | int | complex') -> 'PauliSum':
        if not (isinstance(A, int) or isinstance(A, float) or isinstance(A, complex)):
            raise ValueError(f"Right multiply only changes weight - passed type = {type(A)}")
        return PauliSum([self], weights=[A])

    def _to_pauli_sum(self) -> 'PauliSum':
        return PauliSum([self], weights=[1], phases=[0], dimensions=self.dimensions, standardise=False)

    def __add__(self, A: 'PauliString | PauliSum') -> 'PauliSum':
        if np.all(self.dimensions != A.dimensions):
            raise Exception("To add two PauliStrings, their dimensions"
                            f" {self.dimensions} and {A.dimensions} must be equal")
        if isinstance(A, PauliString):
            return self._to_pauli_sum() + A._to_pauli_sum()
        elif isinstance(A, PauliSum):
            return self._to_pauli_sum() + A
        else:
            raise ValueError("A PauliString can only be added to another PauliString or a PauliSum,"
                             f" not a {type(A)}")
    
    def __eq__(self, other_pauli: 'PauliString') -> bool:
        if not isinstance(other_pauli, PauliString):
            return False
        
        return bool(np.all(self.x_exp == other_pauli.x_exp) and np.all(self.z_exp == other_pauli.z_exp) and np.all(self.dimensions == other_pauli.dimensions))

    def __ne__(self, other_pauli: 'PauliString') -> bool:
        return not self.__eq__(other_pauli)
    
    def __gt__(self, other_pauli: 'PauliString') -> bool:
        """
        Arbitrary but useful in making a standard form

        Note: Some other standard form may be preferable and can be determined here by defining a sorting

        higher powers of x are first, then higher powers of z,
        then sorted by weight, then phase
        """
        for i in range(self.n_qudits()):
            pauli_i = Pauli(x_exp=self.x_exp[i], z_exp=self.z_exp[i], dimension=self.dimensions[i])
            other_pauli_i = Pauli(x_exp=other_pauli.x_exp[i], z_exp=other_pauli.z_exp[i], dimension=other_pauli.dimensions[i])
            if pauli_i > other_pauli_i:
                return True
            elif pauli_i == other_pauli_i:
                continue
        return False
        
    def __hash__(self) -> int:
        return hash((tuple(self.x_exp), tuple(self.z_exp), tuple(self.dimensions)))
    
    def __dict__(self) -> dict:
        return {'x_exp': self.x_exp, 'z_exp': self.z_exp, 'dimensions': self.dimensions}

    def n_qudits(self) -> int:
        return len(self.x_exp)
        
    def n_identities(self) -> 'PauliString':
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

    def symplectic_product(self, A: 'PauliString') -> np.ndarray:
        n = self.n_qudits()
        symp = self.symplectic()
        symp_A = A.symplectic()
        prod = sum([symp[i] * symp_A[i + n] - symp[i + n] * symp_A[i] for i in range(n)]) % self.lcm
        return prod
    
    def amend(self, qudit_index: int, new_x: int, new_z: int) -> 'PauliString':
        if new_x > self.dimensions[qudit_index] or new_z > self.dimensions[qudit_index]:
            raise ValueError(f"Exponents ({new_x, new_z}) cannot be larger than qudit dimension"
                             f" ({self.dimensions[qudit_index]})")
        self.x_exp[qudit_index] = new_x
        self.z_exp[qudit_index] = new_z
        return self
    
    def acquired_phase(self, other_pauli: 'PauliString') -> int:
        # phases acquired when multiplying two Pauli strings
        phi = 1.  # / self.dimensions
        phase = 0
        for i in range(self.n_qudits()):
            phase += phi * (self.x_exp[i] * other_pauli.z_exp[i] + self.z_exp[i] * other_pauli.x_exp[i])
        return phase % self.lcm
    
    def _replace_symplectic(self, symplectic: np.ndarray, qudit_indices: list[int]) -> 'PauliString':
        x_exp_replace = symplectic[0:len(qudit_indices)]
        z_exp_replace = symplectic[len(qudit_indices):2 * len(qudit_indices)]

        x_exp = self.x_exp.copy()
        z_exp = self.z_exp.copy()
        for i, index in enumerate(qudit_indices):
            x_exp[index] = x_exp_replace[i]
            z_exp[index] = z_exp_replace[i]

        return PauliString(x_exp=x_exp, z_exp=z_exp, dimensions=self.dimensions)
    
    def _delete_qudits(self, qudit_indices: list[int], return_new=True) -> 'PauliString':  # not sure if here it is best to return a new object or not
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
        
    def __getitem__(self, key: int | slice) -> 'PauliString | Pauli':
        if isinstance(key, int):
            return self.get_paulis()[key]
        else:
            return PauliString(x_exp=self.x_exp[key], z_exp=self.z_exp[key], dimensions=self.dimensions[key])

    def get_subspace(self, qudit_indices: list[int] | int) -> 'PauliString':
        return PauliString(x_exp=self.x_exp[qudit_indices], z_exp=self.z_exp[qudit_indices],
                           dimensions=self.dimensions[qudit_indices])
    
    def copy(self) -> 'PauliString':
        return PauliString(x_exp=self.x_exp.copy(), z_exp=self.z_exp.copy(), dimensions=self.dimensions.copy())


class PauliSum:
    """
    Lower level class for performing calculations in the symplectic representation

    Represents a weighted sum of Pauli strings with arbitrary phases
    """
    def __init__(self,
                 pauli_list: list[PauliString] | list[Pauli] | list[str] | PauliString | Pauli,
                 weights: list[float | complex] | np.ndarray | float | complex | None = None,
                 phases: list[float] | np.ndarray | None = None,
                 dimensions: list[int] | np.ndarray | None = None,
                 standardise: bool = True):
        """
        Constructor for SymplecticPauli class.

        Parameters
        ----------
        pauli_list : list of Pauli
            The Pauli operators to be represented.
        weights : list of float, optional
            The weights of the Pauli operators.
        phases : list of float, optional
            The phases of the Pauli operators.
        dimensions : int or list of int, optional
            The dimensions for each qudit.

        Raises
        ------
        ValueError
            If the length of pauli_list and weights do not match.
        """
        sanitized_pauli_list, sanitized_dimensions, sanitized_phases, sanitized_weights = self._sanity_checks(
            pauli_list, weights, phases, dimensions
        )

        self.pauli_strings = sanitized_pauli_list
        self.weights = np.asarray(sanitized_weights, dtype=np.complex128)
        self.phases = np.asarray(sanitized_phases, dtype=int)
        self.dimensions = sanitized_dimensions
        
        x_exp = np.zeros((len(self.pauli_strings), len(self.dimensions)))  # ensures we can always index [pauli #, qudit #]
        z_exp = np.zeros((len(self.pauli_strings), len(self.dimensions)))  # ensures we can always index [pauli #, qudit #]

        for i, p in enumerate(self.pauli_strings):
            x_exp[i, :] = p.x_exp
            z_exp[i, :] = p.z_exp
        
        self.x_exp = x_exp
        self.z_exp = z_exp
        self.lcm = np.lcm.reduce(self.dimensions)

        if standardise:
            print('SORTING')  # Keeping this here for now for debug purposes as it can cause issues.
            self.standardise()

    @staticmethod
    def _sanitize_pauli_list(pauli_list: list[PauliString] | list[Pauli] | list[str] | PauliString | Pauli,
                             dimensions: list[int] | np.ndarray | None) -> list[PauliString]:
        
        if isinstance(pauli_list, Pauli):
            pauli_list = [pauli_list]
        if isinstance(pauli_list, PauliString):
            pauli_list = [pauli_list]
        if isinstance(pauli_list, str):
            pauli_list = [pauli_list]

        sanitised_pauli_list = []
        for p in pauli_list:
            if isinstance(p, PauliString):
                sanitised_pauli_list.append(p)
            elif isinstance(p, Pauli):
                sanitised_pauli_list.append(p._to_pauli_string())
            elif isinstance(p, str):
                if dimensions is None:
                    raise SyntaxError("Input of strings into PauliSum requires explicit dimensions input")
                sanitised_pauli_list.append(PauliString(p, dimensions=dimensions))
            else:
                raise TypeError("Pauli list must be a list of PauliString or Pauli objects or strings")

        return sanitised_pauli_list

    @staticmethod
    def _sanitize_dimensions(pauli_list: list[PauliString],
                             dimensions: list[int] | np.ndarray | None = None) -> np.ndarray:
        if dimensions is None and len(pauli_list) == 0:
            return np.empty(0, dtype=int)

        if dimensions is None:
            for i in range(1, len(pauli_list)):
                if not np.array_equal(pauli_list[i].dimensions, pauli_list[0].dimensions):
                    raise ValueError("The dimensions of all Pauli strings must be equal.")
            dimensions = pauli_list[0].dimensions
            
        return np.array(dimensions)

    @staticmethod
    def _sanitize_phases(pauli_list: list[PauliString],
                         phases: list[float] | np.ndarray | None) -> np.ndarray:
        if phases is None:
            return np.zeros(len(pauli_list), dtype=int)

        return np.asarray(phases, dtype=int)

    @staticmethod
    def _sanitize_weights(pauli_list: list[PauliString],
                          weights: list[float | complex] | np.ndarray | float | complex | None) -> np.ndarray:
        if weights is None:
            return np.ones(len(pauli_list))

        if isinstance(weights, (float, complex)):
            return np.ones(len(pauli_list)) * weights

        if isinstance(weights, (np.ndarray, list)) and len(pauli_list) != len(weights):
            raise ValueError(f"Length of Pauli list ({len(pauli_list)}) and weights ({len(weights)}) must be equal.")
        
        return np.asarray(weights, dtype=float)

    def _sanity_checks(self,
                       pauli_list: list[PauliString] | list[Pauli] | list[str] | PauliString | Pauli,
                       weights: list[float | complex] | np.ndarray | float | complex | None,
                       phases: list[float] | np.ndarray | None,
                       dimensions: list[int] | np.ndarray | None) -> tuple[list[PauliString], np.ndarray, np.ndarray, np.ndarray]:
        sanitized_pauli_list = self._sanitize_pauli_list(pauli_list, dimensions)
        sanitized_dimensions = self._sanitize_dimensions(sanitized_pauli_list, dimensions)
        sanitized_phases = self._sanitize_phases(sanitized_pauli_list, phases)
        sanitized_weights = self._sanitize_weights(sanitized_pauli_list, weights)

        return sanitized_pauli_list, sanitized_dimensions, sanitized_phases, sanitized_weights
    
    def n_paulis(self) -> int:
        return len(self.pauli_strings)
    
    def n_qudits(self) -> int:
        return len(self.dimensions)
    
    def n_identities(self):
        """
        Get the number of identities in the PauliSum
        :return: The number of identities
        """
        n_is = []
        for i in range(self.n_paulis()):
            n_is.append(self.pauli_strings[i].n_identities())

    def standardise(self):
        """
        Standardises the PauliSum object by combining equivalent Paulis and
        adding phase factors to the weights then resetting the phases.
        """
        # combine equivalent
        # self.combine_equivalent_paulis()
        # sort
        self.weights = [x for _, x in sorted(zip(self.pauli_strings, self.weights))]
        self.phases = [x for _, x in sorted(zip(self.pauli_strings, self.phases))]
        self.pauli_strings = sorted(self.pauli_strings)
        # add phase factors to weights then reset phases
        new_weights = np.zeros(self.n_paulis(), dtype=np.complex128)
        for i in range(self.n_paulis()):
            phase = self.phases[i]
            omega = np.exp(2 * np.pi * 1j * phase / self.lcm)
            new_weights[i] = self.weights[i] * omega
        self.phases = np.zeros(self.n_paulis(), dtype=int)
        self.weights = new_weights

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.pauli_strings[key]
        elif isinstance(key, slice):
            return PauliSum(self.pauli_strings[key], self.weights[key], self.phases[key], self.dimensions, False)
        elif isinstance(key, tuple):
            if len(key) != 2:
                raise ValueError("Tuple key must be of length 2")
            if isinstance(key[0], int):
                return self.pauli_strings[key[0]][key[1]]
            if isinstance(key[0], slice):
                pauli_strings_all_qubits = self.pauli_strings[key[0]]
                pauli_strings = [p[key[1]] for p in pauli_strings_all_qubits]
                return PauliSum(pauli_strings, self.weights[key[0]], self.phases[key[0]], self.dimensions[key[1]], False)
        else:
            raise TypeError(f"Key must be int or slice, not {type(key)}")

    def __add__(self, A: 'Pauli | PauliString | PauliSum') -> 'PauliSum':
        if isinstance(A, PauliString) or isinstance(A, Pauli):
            A = A._to_pauli_sum()
        elif isinstance(A, Pauli):
            A = A._to_pauli_sum()
            
        new_pauli_list = self.pauli_strings + A.pauli_strings
        new_weights = np.concatenate([self.weights, A.weights])
        new_phases = np.concatenate([self.phases, A.phases])
        return PauliSum(list(new_pauli_list), new_weights, new_phases, self.dimensions, False)

    def __sub__(self, A: 'PauliSum') -> 'PauliSum':
        new_pauli_list = self.pauli_strings + A.pauli_strings
        new_weights = np.concatenate([self.weights, -np.array(A.weights)])
        new_phases = np.concatenate([self.phases, A.phases])
        return PauliSum(list(new_pauli_list), new_weights, new_phases, self.dimensions, False)
    
    def __matmul__(self, A: 'Pauli | PauliString | PauliSum') -> 'PauliSum':
        """
        @ is the operator for tensor product
        """
        if isinstance(A, PauliString):
            A = A._to_pauli_sum()
        elif isinstance(A, Pauli):
            A = A._to_pauli_sum()
        
        new_dimensions = np.hstack((self.dimensions, A.dimensions))
        new_lcm = np.lcm.reduce(new_dimensions)
        new_pauli_list = []
        new_weights = []
        new_phases = []
        for i in range(self.n_paulis()):
            for j in range(A.n_paulis()):
                new_pauli_list.append(self.pauli_strings[i] @ A.pauli_strings[j])
                new_weights.append(self.weights[i] * A.weights[j])
                new_phases.append(((self.phases[i] + A.phases[j]) % new_lcm))
        output_pauli = PauliSum(new_pauli_list, new_weights, new_phases, self.dimensions, False)
        return output_pauli

    def __mul__(self, A: 'PauliSum | PauliString | Pauli | float | int | complex') -> 'PauliSum':
        """
        Operator multiplication on two SymplecticPauli objects or multiplication of weights by constant
        """

        if isinstance(A, (int, float)):
            return PauliSum(list(self.pauli_strings), np.array(self.weights) * A, self.phases)
        elif isinstance(A, PauliString):
            return self * A._to_pauli_sum()
        elif not isinstance(A, PauliSum):
            raise ValueError("Multiplication only supported with SymplecticPauli objects or scalar")

        new_p_sum = []
        new_weights = []
        new_phases = []
        for i in range(self.n_paulis()):
            for j in range(A.n_paulis()):
                new_p_sum.append(self.pauli_strings[i] * A.pauli_strings[j])
                new_weights.append(self.weights[i] * A.weights[j])
                acquired_phase = self.pauli_strings[i].acquired_phase(A.pauli_strings[j])
                new_phases.append((self.phases[i] + A.phases[j] + acquired_phase) % self.lcm)
        output_pauli = PauliSum(new_p_sum, new_weights, new_phases, self.dimensions, False)

        return output_pauli
    
    def __truediv__(self, A: 'PauliSum | PauliString | Pauli') -> 'PauliSum':
        if not isinstance(A, (int, float)):
            raise ValueError("Division only supported with scalar")
        return self * (1 / A)
    
    def __eq__(self, value: 'PauliSum') -> bool:
        if not isinstance(value, PauliSum):
            return False
        t1 = np.all(self.pauli_strings == value.pauli_strings)
        t2 = np.all(self.weights == value.weights)
        t3 = np.all(self.phases == value.phases)
        return bool(t1 and t2 and t3)
    
    def __ne__(self, value: 'PauliSum') -> bool:
        return not self == value
    
    def __hash__(self) -> int:
        return hash((tuple(self.pauli_strings), tuple(self.weights), tuple(self.phases), tuple(self.dimensions)))
    
    def __dict__(self) -> dict:
        return {'pauli_strings': self.pauli_strings, 'weights': self.weights, 'phases': self.phases}
    
    def _to_standard_form(self):
        """
        Alters self to a standard form

        Sorting could be improved but ok for now

        """
        new_weights = np.array([x for _, x in sorted(zip(self.pauli_strings, self.weights))])
        new_phases = np.array([x for _, x in sorted(zip(self.pauli_strings, self.phases))])
        new_strings = sorted(self.pauli_strings)
        self.pauli_strings = new_strings
        self.weights = new_weights
        self.phases = new_phases

    def combine_equivalent_paulis(self):
        self.standardise()  # makes sure all phases are 0
        # combine equivalent Paulis
        to_delete = []
        for i in reversed(range(self.n_paulis())):
            for j in range(i + 1, self.n_paulis()):
                if self.pauli_strings[i] == self.pauli_strings[j]:
                    self.weights[i] = self.weights[i] + self.weights[j]
                    to_delete.append(j)
        self._delete_paulis(to_delete)

        # remove zero weight Paulis
        to_delete = []
        for i in range(self.n_paulis()):
            if self.weights[i] == 0:
                to_delete.append(i)
        self._delete_paulis(to_delete)

    def remove_trivial_paulis(self):
        # If entire Pauli string is I, remove it
        to_delete = []
        for i in range(self.n_paulis()):
            if np.all(self.x_exp[i, :] == 0) and np.all(self.z_exp[i, :] == 0):
                to_delete.append(i)
        self._delete_paulis(to_delete)

    def remove_trivial_qudits(self):
        # If entire Pauli string is I, remove it
        to_delete = []
        for i in range(self.n_qudits()):
            if np.all(self.x_exp[:, i] == 0) and np.all(self.z_exp[:, i] == 0):
                to_delete.append(i)
        self._delete_qudits(to_delete)

    def symplectic_matrix(self) -> np.ndarray:
        symplectic = np.zeros([self.n_paulis(), 2 * self.n_qudits()])
        for i, p in enumerate(self.pauli_strings):
            symplectic[i, :] = p.symplectic()
        return symplectic
        
    def is_x(self) -> bool:
        # check whether self has only X component
        # Outputs: (bool) - True if self has only X component, False otherwise
        return not np.any(self.z_exp)

    def is_z(self) -> bool:
        # check whether self has only Z component
        # Outputs: (bool) - True if self has only Z component, False otherwise
        return not np.any(self.x_exp)

    def is_commuting(self, pauli_string_indexes: list[int] | None = None) -> bool:
        # check whether the set of Paulis are pairwise commuting
        # Outputs:  (bool) - True if self is pairwise commuting set of Paulis
        spm = self.symplectic_product_matrix()
        if pauli_string_indexes is None:
            return not np.any(spm)
        else:
            i, j = pauli_string_indexes[0], pauli_string_indexes[1]
            return not spm[i, j]

    def select_pauli_string(self, pauli_index: int) -> PauliString:
        # Inputs:
        #     pauli_index - (int) - index of Pauli to be returned
        # Outputs:
        #     (PauliString) - the indexed Pauli in self
        return self.pauli_strings[pauli_index]

    def _delete_paulis(self, pauli_indices: list[int] | int):
        # Inputs:
        #     pauli_indices - (list of int or int)
        if isinstance(pauli_indices, int):
            pauli_indices = [pauli_indices]

        new_weights = np.delete(self.weights, pauli_indices)
        new_phases = np.delete(self.phases, pauli_indices)
        new_x_exp = np.delete(self.x_exp, pauli_indices, axis=0)
        new_z_exp = np.delete(self.z_exp, pauli_indices, axis=0)

        for i in sorted(pauli_indices, reverse=True):  # sort in reverse order to avoid index shifting
            self.pauli_strings = list(self.pauli_strings)  # Convert to list
            del self.pauli_strings[i]
            self.pauli_strings = np.array(self.pauli_strings, dtype=object)  # Convert back to NumPy array

        self.weights = new_weights
        self.phases = new_phases
        self.x_exp = new_x_exp
        self.z_exp = new_z_exp

    def _delete_qudits(self, qudit_indices: list[int] | int):
        # Inputs:
        #     qudit_indices - (list of int)
        if isinstance(qudit_indices, int):
            qudit_indices = [qudit_indices]
        
        new_pauli_strings = []
        for p in self.pauli_strings:
            new_pauli_strings.append(p._delete_qudits(qudit_indices))

        self.pauli_strings = np.array(new_pauli_strings, dtype=object)
        self.x_exp = np.delete(self.x_exp, qudit_indices, axis=1)
        self.z_exp = np.delete(self.z_exp, qudit_indices, axis=1)
        self.dimensions = np.delete(self.dimensions, qudit_indices)

        self.lcm = np.lcm.reduce(self.dimensions)

    def copy(self) -> 'PauliSum':
        return PauliSum([ps for ps in self.pauli_strings], self.weights.copy(), self.phases.copy(), self.dimensions.copy(), False)

    def symplectic_product_matrix(self) -> np.ndarray:
        """
        An n x n matrix, n is the number of Paulis.
        The entry S[i, j] is the symplectic product of the ith Pauli and the jth Pauli.
        """
        n = self.n_paulis()
        # list_of_symplectics = self.symplectic_matrix()

        spm = np.zeros([n, n], dtype=int)
        for i in range(n):
            for j in range(n):
                if i > j:
                    spm[i, j] = symplectic_product(self.pauli_strings[i], self.pauli_strings[j])
        spm = spm + spm.T
        return spm
    
    def __str__(self) -> str:
        p_string = ''
        max_str_len = max([len(f'{self.weights[i]}') for i in range(self.n_paulis())])
        for i in range(self.n_paulis()):
            pauli_string = self.pauli_strings[i]
            qudit_string = ''.join(['x' + f'{pauli_string.x_exp[j]}' +
                                    'z' + f'{pauli_string.z_exp[j]} ' for j in range(self.n_qudits())])
            n_spaces = max_str_len - len(f'{self.weights[i]}')
            p_string += f'{self.weights[i]}' + ' ' * n_spaces + '|' + qudit_string + f'| {self.phases[i]} \n'
        return p_string
    
    def get_subspace(self, qudit_indices: list[int], pauli_indices: list | None = None):
        """
        Get the subspace of the PauliSum corresponding to the qudit indices for the given Paulis
        Not strictly a subspace if we restrict the Pauli indices, so we could rename but this is still quite clear

        :param qudit_indices: The indices of the qudits to get the subspace for
        :param pauli_indices: The indices of the Paulis to get the subspace for
        :return: The subspace of the PauliSum
        """
        if pauli_indices is None:
            pauli_indices = np.arange(self.n_paulis()).tolist()

        dimensions = self.dimensions[qudit_indices]
        pauli_list = []
        for i in pauli_indices:
            p = self.pauli_strings[i]
            p = p.get_subspace(qudit_indices)
            pauli_list.append(p)
        return PauliSum(pauli_list, self.weights[pauli_indices], self.phases[pauli_indices], dimensions, False)

    def matrix_form(self, pauli_string_index: int | None = None) -> scipy.sparse.csr_matrix:
        """
        Returns
        -------
        scipy.sparse.csr_matrix
            Matrix representation of input Pauli.
        """
        if pauli_string_index is not None:
            ps = self.select_pauli_string(pauli_string_index)
            ps = ps._to_pauli_sum()
            return ps.matrix_form()

        else:
            list_of_pauli_matrices = []
            for i in range(self.n_paulis()):
                X, Z, dim, phase = int(self.x_exp[i, 0]), int(self.z_exp[i, 0]), self.dimensions[0], self.phases[i]
                h = self.xz_mat(dim, X, Z)
                
                for n in range(1, self.n_qudits()):
                    X, Z, dim, phase = int(self.x_exp[i, n]), int(self.z_exp[i, n]), self.dimensions[n], self.phases[i]
                    h_next = self.xz_mat(dim, X, Z)

                    h = scipy.sparse.kron(h, h_next, format="csr")
                list_of_pauli_matrices.append(np.exp(phase * 2 * np.pi * 1j / self.lcm) * self.weights[i] * h)
            m = sum(list_of_pauli_matrices)
        return m
    
    def acquire_phase(self, phases: list[int], pauli_index: int | list[int] | None = None):
        if pauli_index is not None:
            if isinstance(pauli_index, int):
                pauli_index = [pauli_index]
            elif len(pauli_index) != len(phases):
                raise ValueError(f"Number of phases ({len(phases)}) must be equal to number of Paulis ({len(pauli_index)})")
            else:
                raise ValueError(f"pauli_index must be int, list, or np.ndarray, not {type(pauli_index)}")
            for i in pauli_index:
                self.phases[i] = (self.phases[i] + phases) % self.lcm
        else:
            if len(phases) != self.n_paulis():
                raise ValueError(f"Number of phases ({len(phases)}) must be equal to number of Paulis ({self.n_paulis()})")
            new_phase = (np.array(self.phases) + np.array(phases)) % self.lcm
        self.phases = new_phase

    @staticmethod
    def xz_mat(d: int, aX: int, aZ: int) -> scipy.sparse.csr_matrix:
        """
        Temporary function for pauli reduction.

        Function for creating generalized Pauli matrix.

        Parameters
        ----------
        d : int
            Dimension of the qudit
        aX : int
            X-part of the Pauli matrix
        aZ : int
            Z-part of the Pauli matrix

        Returns
        -------
        scipy.sparse.csr_matrix
            Generalized Pauli matrix
        """
        omega = np.exp(2 * np.pi * 1j / d)
        aa0 = np.array([1 for i in range(d)])
        aa1 = np.array([i for i in range(d)])
        aa2 = np.array([(i - aX) % d for i in range(d)])
        X = scipy.sparse.csr_matrix((aa0, (aa1, aa2)))
        aa0 = np.array([omega**(i * aZ) for i in range(d)])
        aa1 = np.array([i for i in range(d)])
        aa2 = np.array([i for i in range(d)])
        Z = scipy.sparse.csr_matrix((aa0, (aa1, aa2)))
        if (d == 2) and (aX % 2 == 1) and (aZ % 2 == 1):
            return 1j * (X @ Z)
        return X @ Z


# the symplectic inner product of two pauli objects (each with a single Pauli)
def symplectic_product(pauli_string: PauliString, pauli_string2: PauliString) -> bool:
    # Inputs:
    #     pauli_string - (PauliString)
    #     pauli_string2 - (PauliString)
    # Outputs:
    #     (bool) - quditwise inner product of Paulis

    if any(pauli_string.dimensions - pauli_string.dimensions):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    sp = 0
    for i in range(pauli_string.n_qudits()):
        sp += (pauli_string.x_exp[i] * pauli_string2.z_exp[i] - pauli_string.z_exp[i] * pauli_string2.x_exp[i])
    return sp % pauli_string.lcm


def string_to_symplectic(string: str) -> tuple[np.ndarray, int]:
    # split into single qubit paulis by spaces
    substrings = string.split()
    local_symplectics = []
    phases = []
    for s in substrings:
        match = re.match(r'x(\d+)z(\d+)(?:p(\d+))?', s)
        if not match:
            raise ValueError(f"Invalid Pauli string: {s}")
        else:
            x = int(match.group(1))
            z = int(match.group(2))
            p = int(match.group(3)) if match.group(3) is not None else 0
            local_symplectics.append((x, z))
            phases.append(p)
    
    symplectic = np.array(local_symplectics).T
    return symplectic.flatten(), sum(phases)


def symplectic_to_string(symplectic: np.ndarray, dimension: int) -> str:
    if dimension == 2:
        if symplectic[0] == 0 and symplectic[1] == 0:
            return 'x0z0'
        elif symplectic[0] == 1 and symplectic[1] == 0:
            return 'x1z0'
        elif symplectic[0] == 0 and symplectic[1] == 1:
            return 'x0z1'
        elif symplectic[0] == 1 and symplectic[1] == 1:
            return 'x1z1'
        else:
            raise Exception("Symplectic vector must be of the form (0, 0), (1, 0), (0, 1), or (1, 1)")
    else:
        return f'x{symplectic[0]}z{symplectic[1]}'


if __name__ == '__main__':
    
    import sys
    sys.path.append("./")
    dims = [2, 2]
    x1x1 = PauliSum(['x0z1 x0z0', 'x1z0 x0z0', 'x1z0 x0z0'], dimensions=dims)
    # x1x1.weights = [1., 10,]

    print(x1x1)
    print(x1x1.symplectic_matrix())
    print(x1x1.symplectic_product_matrix())
    print(x1x1.is_commuting())
    print(bool(x1x1.is_commuting(pauli_string_indexes=[0, 1])))
    print(bool(x1x1.is_commuting(pauli_string_indexes=[1, 2])))
