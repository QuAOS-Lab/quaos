# from __future__ import annotations
from typing import Union, overload
import numpy as np
import scipy
from numpy.typing import NDArray
from .pauli import Pauli
from .pauli_string import PauliString
# from .utils import to_pauli_sum, to_pauli_string


class PauliSum:
    """
    Lower level class for performing calculations in the symplectic representation

    Represents a weighted sum of Pauli strings with arbitrary phases
    """
    def __init__(self,
                 pauli_list: Union[list[PauliString], list[Pauli], list[str], PauliString, Pauli],
                 weights: list[float | complex] | np.ndarray | float | complex | None = None,
                 phases: list[int] | np.ndarray | None = None,
                 dimensions: list[int] | np.ndarray | None = None,
                 standardise: bool = True):
        """
        TODO: Change everything possible to numpy arrays.
        TODO: Remove self.xz_mat - should be in a utils module
        TODO: Add stack method to concatenate additional strings or sums (could use utils concatenate)

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
        self.dimensions = sanitized_dimensions
        self.lcm = np.lcm.reduce(self.dimensions)
        self.phases = np.asarray(sanitized_phases, dtype=int) % self.lcm
        
        self._set_exponents()

        if standardise:
            self.standardise()

    def _set_exponents(self):
        x_exp = np.zeros((len(self.pauli_strings), len(self.dimensions)))  # we can always index [pauli #, qudit #]
        z_exp = np.zeros((len(self.pauli_strings), len(self.dimensions)))  # we can always index [pauli #, qudit #]

        for i, p in enumerate(self.pauli_strings):
            x_exp[i, :] = p.x_exp
            z_exp[i, :] = p.z_exp
        
        self.x_exp = x_exp
        self.z_exp = z_exp

    @staticmethod
    def _sanitize_pauli_list(pauli_list: Union[list[PauliString], list[Pauli], list[str], PauliString, Pauli],
                             dimensions: list[int] | np.ndarray | None) -> list[PauliString]:
        
        if isinstance(pauli_list, Pauli) or isinstance(pauli_list, PauliString) or isinstance(pauli_list, str):
            pauli_list = [pauli_list]

        sanitised_pauli_list = []
        for p in pauli_list:
            if isinstance(p, PauliString):
                sanitised_pauli_list.append(p)
            elif isinstance(p, Pauli):
                p_ = PauliString([p.x_exp], [p.z_exp], dimensions=[p.dimension])
                sanitised_pauli_list.append(p_)
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
                         phases: list[int] | np.ndarray | None) -> np.ndarray:
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
        
        return np.asarray(weights, dtype=complex)

    def _sanity_checks(self,
                       pauli_list: Union[list[PauliString], list[Pauli], list[str], PauliString, Pauli],
                       weights: list[float | complex] | np.ndarray | float | complex | None,
                       phases: list[int] | np.ndarray | None,
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
    
    def shape(self) -> tuple[int, int]:
        return self.n_paulis(), self.n_qudits()
    
    def n_identities(self):
        """
        Get the number of identities in the PauliSum
        :return: The number of identities
        """
        n_is = []
        for i in range(self.n_paulis()):
            n_is.append(self.pauli_strings[i].n_identities())

    def phase_to_weight(self):
        new_weights = np.zeros(self.n_paulis(), dtype=np.complex128)
        for i in range(self.n_paulis()):
            phase = self.phases[i]
            omega = np.exp(2 * np.pi * 1j * phase / self.lcm)
            new_weights[i] = self.weights[i] * omega
        self.phases = np.zeros(self.n_paulis(), dtype=int)
        self.weights = new_weights

    @overload
    def __getitem__(self, key: int | tuple[int, slice] | tuple[int, list[int]] | tuple[int, NDArray[np.int_]]) -> 'PauliString':
        ...

    @overload
    def __getitem__(self, key: slice | tuple[slice, int] | tuple[slice, slice] | list[int] | NDArray[np.int_]) -> 'PauliSum':
        ...

    @overload
    def __getitem__(self, key: tuple[int, int]) -> 'Pauli':
        ...

    def __getitem__(self, key):
        # TODO: enable list for either input
        if isinstance(key, int):
            return self.pauli_strings[key]
        elif isinstance(key, slice):
            return PauliSum(self.pauli_strings[key], self.weights[key], self.phases[key], self.dimensions, False)
        elif isinstance(key, np.ndarray) or isinstance(key, list):
            ps = [self.pauli_strings[i] for i in key]
            w = np.asarray([self.weights[i] for i in key])
            p = np.asarray([self.phases[i] for i in key])
            return PauliSum(ps, w, p, self.dimensions, False)
        elif isinstance(key, tuple):
            if len(key) != 2:
                raise ValueError("Tuple key must be of length 2")
            if isinstance(key[0], int):
                return self.pauli_strings[key[0]][key[1]]
            elif isinstance(key[0], slice) or isinstance(key, np.ndarray) or isinstance(key, list):
                pauli_strings_all_qubits = self.pauli_strings[key[0]]
                pauli_strings = [p[key[1]] for p in pauli_strings_all_qubits]
                if isinstance(key[1], int):
                    return PauliSum(pauli_strings, self.weights[key[0]], self.phases[key[0]], np.asarray([self.dimensions[key[1]]]), False)
                elif isinstance(key[1], slice):
                    return PauliSum(pauli_strings, self.weights[key[0]], self.phases[key[0]], self.dimensions[key[1]], False)
        else:
            raise TypeError(f"Key must be int or slice, not {type(key)}")

    def __setitem__(self, key, value):
        if isinstance(key, int):  # key indexes the pauli_string to be replaced by value
            self.pauli_strings[key] = value
        elif isinstance(key, slice):
            self.pauli_strings[key] = value
        elif isinstance(key, tuple):
            if len(key) != 2:
                raise ValueError("Tuple key must be of length 2")
            if isinstance(key[0], int):
                if isinstance(key[1], int):  # key[0] indexes the pauli string, key[1] indexes the qudit
                    self.pauli_strings[key[0]][key[1]] = value
                elif isinstance(key[1], slice):  # key[0] indexes the pauli string, key[1] indexes the qudits
                    self.pauli_strings[key[0]][key[1]] = value
            if isinstance(key[0], slice):
                if isinstance(key[1], int):  # key[0] indexes the pauli strings, key[1] indexes the qudit
                    for i in np.arange(self.n_paulis())[key[0]]:
                        self.pauli_strings[i][key[1]] = value
                elif isinstance(key[1], slice):  # key[0] indexes the pauli strings, key[1] indexes the qudits
                    for i_val, i in enumerate(np.arange(self.n_paulis())[key[0]]):
                        print(i, value[int(i_val)])
                        self.pauli_strings[i][key[1]] = value[int(i_val)]
        self._set_exponents()  # update exponents x_exp and z_exp

    def __add__(self, A: 'Pauli | PauliString | PauliSum') -> 'PauliSum':
        if isinstance(A, PauliString) or isinstance(A, Pauli):
            A = PauliSum([A])
            
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
        if isinstance(A, PauliString) or isinstance(A, Pauli):
            A = PauliSum([A])
        
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
        output_pauli = PauliSum(new_pauli_list, new_weights, new_phases, new_dimensions, False)
        return output_pauli

    def __mul__(self, A: 'PauliSum | PauliString | Pauli | float | int | complex') -> 'PauliSum':
        """
        Operator multiplication on two SymplecticPauli objects or multiplication of weights by constant
        """
        if isinstance(A, (int, float)):
            return PauliSum(list(self.pauli_strings), np.array(self.weights) * A, self.phases)
        elif isinstance(A, PauliString):
            return self * PauliSum([A])
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
        t2 = np.all(np.asarray(self.weights) - np.asarray(value.weights) <= 1e-15)
        t3 = np.all(self.phases == value.phases)
        return bool(t1 and t2 and t3)
    
    def __ne__(self, value: 'PauliSum') -> bool:
        return not self == value
    
    def __hash__(self) -> int:
        return hash((tuple(self.pauli_strings), tuple(self.weights), tuple(self.phases), tuple(self.dimensions)))
    
    def __dict__(self) -> dict:
        return {'pauli_strings': self.pauli_strings, 'weights': self.weights, 'phases': self.phases}

    def standardise(self):
        """
        Standardises the PauliSum object by combining equivalent Paulis and
        adding phase factors to the weights then resetting the phases.
        """
        # combine equivalent
        # self.combine_equivalent_paulis()
        # sort
        self.phase_to_weight()
        self.weights = [x for _, x in sorted(zip(self.pauli_strings, self.weights))]
        # self.phases = [x for _, x in sorted(zip(self.pauli_strings, self.phases))]
        self.pauli_strings = sorted(self.pauli_strings)
        
    def combine_equivalent_paulis(self):
        self.standardise()  # makes sure all phases are 0
        # combine equivalent Paulis
        to_delete = []
        for i in reversed(range(self.n_paulis())):
            for j in range(i + 1, self.n_paulis()):
                if self.pauli_strings[i] == self.pauli_strings[j]:
                    self.weights[i] = self.weights[i] + self.weights[j]
                    if j not in to_delete:
                        to_delete.append(j)
        self.delete_paulis(to_delete)

        # remove zero weight Paulis
        to_delete = []
        for i in range(self.n_paulis()):
            if self.weights[i] == 0:
                to_delete.append(i)
        self.delete_paulis(to_delete)

    def remove_trivial_paulis(self):
        # If entire Pauli string is I, remove it
        to_delete = []
        for i in range(self.n_paulis()):
            if np.all(self.x_exp[i, :] == 0) and np.all(self.z_exp[i, :] == 0):
                to_delete.append(i)
        self.delete_paulis(to_delete)

    def remove_trivial_qudits(self):
        # If entire qudit is I, remove it
        to_delete = []
        for i in range(self.n_qudits()):
            if np.all(self.x_exp[:, i] == 0) and np.all(self.z_exp[:, i] == 0):
                to_delete.append(i)
        if len(to_delete) > 0:
            self.delete_qudits(to_delete)

    def symplectic(self) -> np.ndarray:
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

    def delete_paulis(self, pauli_indices: list[int] | int):
        # Inputs:
        #     pauli_indices - (list of int or int)
            
        if isinstance(pauli_indices, int):
            pauli_indices = [pauli_indices]

        new_weights = np.delete(self.weights, pauli_indices)
        new_phases = np.delete(self.phases, pauli_indices)
        new_x_exp = np.delete(self.x_exp, pauli_indices, axis=0)
        new_z_exp = np.delete(self.z_exp, pauli_indices, axis=0)

        for i in sorted(pauli_indices, reverse=True):  # sort in reverse order to avoid index shifting # Convert to list
            del self.pauli_strings[i]

        self.weights = new_weights
        self.phases = new_phases
        self.x_exp = new_x_exp
        self.z_exp = new_z_exp

    def delete_qudits(self, qudit_indices: list[int] | int):
        # Inputs:
        #     qudit_indices - (list of int)
        if isinstance(qudit_indices, int):
            qudit_indices = [qudit_indices]
        
        new_pauli_strings = []
        for p in self.pauli_strings:
            new_pauli_strings.append(p._delete_qudits(qudit_indices))

        self.pauli_strings = new_pauli_strings
        self.x_exp = np.delete(self.x_exp, qudit_indices, axis=1)
        self.z_exp = np.delete(self.z_exp, qudit_indices, axis=1)
        self.dimensions = np.delete(self.dimensions, qudit_indices)

        self.lcm = np.lcm.reduce(self.dimensions)

    def copy(self) -> 'PauliSum':
        return PauliSum([ps.copy() for ps in self.pauli_strings], self.weights.copy(), self.phases.copy(), self.dimensions.copy(), False)

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
                    spm[i, j] = self.pauli_strings[i].symplectic_product(self.pauli_strings[j])
        spm = spm + spm.T
        return spm
    
    def __str__(self) -> str:
        p_string = ''
        max_str_len = max([len(f'{self.weights[i]}') for i in range(self.n_paulis())])
        for i in range(self.n_paulis()):
            pauli_string = self.pauli_strings[i]
            qudit_string = ''.join(['x' + f'{pauli_string.x_exp[j]}' + 'z' + f'{pauli_string.z_exp[j]} ' for j in range(self.n_qudits())])
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
            ps = PauliSum([ps])
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

    def reorder(self, order: list[int]):
        """
        Reorder the Paulis in the PauliSum. If a set of indices are not in the list, they are
        appended to the end in the original order.

        e.g. reorder[10, 42] will put 10th Pauli first and 42nd Pauli second, followed by the remaining paulis
        in their original order
        """
        if len(order) != self.n_paulis():
            for i in range(self.n_paulis()):
                if i not in order:
                    order.append(i)
        self.pauli_strings = [self.pauli_strings[i] for i in order]
        self.weights = np.array([self.weights[i] for i in order])
        self.phases = np.array([self.phases[i] for i in order])
        self.x_exp = np.array([self.x_exp[i] for i in order])
        self.z_exp = np.array([self.z_exp[i] for i in order])

    @staticmethod
    def xz_mat(d: int, aX: int, aZ: int) -> scipy.sparse.csr_matrix:
        """
        TODO: Move this to a better location and amend where it is used in the Pauli reduction code

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
