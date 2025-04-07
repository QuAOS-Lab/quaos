import sys
import numpy as np
import itertools
import re

sys.path.append("./")


class Pauli:
    def __init__(self, x_exp, z_exp=None, dimension=2):

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
            x_exp = int(x_exp[1])
            z_exp = int(x_exp[3])
        self.x_exp = x_exp
        self.z_exp = z_exp
        self.dimension = dimension

        if self.dimension - 1 < x_exp or self.dimension - 1 < z_exp:
            raise ValueError(f"Dimension {self.dimension} is too small for exponents {self.x_exp} and {self.z_exp}")

    def __mul__(self, A):
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
    
    def __str__(self):
        return f'x{self.x_exp}z{self.z_exp}'
    
    def __matmul__(self, A):
        return PauliString(x_exp=[self.x_exp] + [A.x_exp], z_exp=[self.z_exp] + [A.z_exp],
                           dimensions=[self.dimension] + [A.dimension])

    def __add__(self, A):
        ps1 = PauliString(x_exp=[self.x_exp], z_exp=[self.z_exp], dimensions=[self.dimension])
        if isinstance(A, Pauli):
            ps2 = PauliString(x_exp=[A.x_exp], z_exp=[A.z_exp], dimensions=[A.dimension])
        elif isinstance(A, PauliString) or isinstance(A, PauliSum):
            ps2 = A
        else:
            raise Exception(f"Cannot add Pauli with type {type(A)}")
        return ps1 + ps2
    
    def to_pauli_string(self):
        return PauliString(x_exp=[self.x_exp], z_exp=[self.z_exp], dimensions=[self.dimension])
    
    def to_pauli_sum(self):
        return PauliSum([self.to_pauli_string()])


class Xnd(Pauli):
    def __init__(self, x_exp, dimension):
        super().__init__(x_exp, 0, dimension)


class Ynd(Pauli):
    def __init__(self, y_exp, dimension):
        super().__init__(y_exp, y_exp, dimension)


class Znd(Pauli):
    def __init__(self, z_exp, dimension):
        super().__init__(0, z_exp, dimension)


class Id(Pauli):
    def __init__(self, dimension):
        super().__init__(0, 0, dimension)


class PauliString:
 
    def __init__(self, x_exp, z_exp=None, dimensions=None):

        # NOTE: There is an important generalisation to do for when dimensions are not equal
        #       this I assume is where the lcm is used

        # NOTE: Either we also define phases here or we need an extra rule in the PauliSum class to do multiplication.
        #       It would be neater here as then we can do multiplication of two PauliStrings and get the phase right

        if isinstance(x_exp, str):
            if z_exp is not None:
                raise Warning('If input string is provided, z_exp is unnecessary')
            xz_exponents = re.split('x|z', x_exp)[1:]
            x_exp = np.array(xz_exponents[0::2], dtype=int)
            z_exp = np.array(xz_exponents[1::2], dtype=int)
        elif isinstance(x_exp, list):
            x_exp = np.array(x_exp)
            z_exp = np.array(z_exp)
        self.x_exp = x_exp
        self.z_exp = z_exp
        if dimensions is None:
            dimensions = (max(max(self.x_exp), max(self.z_exp)) + 1) * np.ones(len(self.x_exp))
        self.dimensions = np.asarray(dimensions)
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
    
    def __repr__(self):
        return f"Pauli(x_exp={self.x_exp}, z_exp={self.z_exp}, dimension={self.dimension})"
    
    def n_qudits(self):
        return len(self.x_exp)
    
    def __str__(self):
        p_string = ''
        for i in range(self.n_qudits()):
            p_string += 'x' + f'{self.x_exp[i]}' + 'z' + f'{self.z_exp[i]} '
        return p_string
    
    def __matmul__(self, A):
        new_x_exp = np.concatenate((self.x_exp, A.x_exp))
        new_z_exp = np.concatenate((self.z_exp, A.z_exp))
        new_dims = np.concatenate((self.dimensions, A.dimensions))
        return PauliString(new_x_exp, new_z_exp, new_dims)

    def __mul__(self, A):
        # returns SymplecticPauli * SymplecticPauli
        if isinstance(A, PauliString):
            if np.any(self.dimensions != A.dimensions):
                raise Exception("To multiply two PauliStrings, their dimensions"
                                f" {self.dimensions} and {A.dimensions} must be equal")
            x_new = np.mod(self.x_exp + A.x_exp, (self.dimensions))
            z_new = np.mod(self.z_exp + A.z_exp, (self.dimensions))
            return PauliString(x_new, z_new, self.dimensions)
        elif isinstance(A, PauliSum):
            return self.to_symplectic() * A
        elif isinstance(A, str):
            return self.to_symplectic() * PauliSum(A)
        elif isinstance(A, float) or isinstance(A, int):
            return PauliSum(self, weights=A)
        else:
            raise ValueError(f"Cannot multiply PauliString with type {type(A)}")
        
    def to_pauli_sum(self, weight=None, phase=None):
        return PauliSum(self, weight, phase)
    
    def _to_pauli_sum(self):
        return PauliSum([self], weights=[1], phases=[0])

    def __add__(self, A):
        if np.all(self.dimensions != A.dimensions):
            raise Exception("To add two PauliStrings, their dimensions"
                            f" {self.dimensions} and {A.dimensions} must be equal")
        return self._to_pauli_sum() + A._to_pauli_sum()
    
    def get_paulis(self):
        """
        Get a list of Pauli objects from the PauliString
        :return: A list of Pauli objects
        """
        return [Pauli(x_exp=self.x_exp[i], z_exp=self.z_exp[i], dimension=self.dimensions[i]) for i in range(len(self.x_exp))]

    def symplectic(self):
        symp = np.zeros(2 * self.n_qudits())
        symp[0:self.n_qudits()] = self.x_exp
        symp[self.n_qudits():2 * self.n_qudits()] = self.z_exp
        return symp

    def symplectic_product(self, A):
        n = self.n_qudits()
        symp = self.symplectic()
        symp_A = A.symplectic()
        prod = sum([symp[i] * symp_A[i + n] - symp[i + n] * symp_A[i] for i in range(n)]) % self.dimensions
        return self * A - A * self


class PauliSum:
    """
    Lower level class for performing calculations in the symplectic representation
    """
    def __init__(self, pauli_list, weights=None, phases=None, dimensions=None):
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
        # sanity checks
        if weights is None:
            weights = np.ones(len(pauli_list))
        if phases is None:
            phases = np.zeros(len(pauli_list))
        if len(pauli_list) != len(weights):
            raise ValueError(f"Length of Pauli list ({len(pauli_list)}) and weights ({len(weights)}) must be equal.")
        if not isinstance(pauli_list, list):
            pauli_list = [pauli_list]
        # check all elements of pauli_list are PauliString objects
        if not all(isinstance(p, PauliString) for p in pauli_list):
            pauli_list = [PauliString(pauli_list[i].x_exp, dimensions=dimensions) for i in range(len(pauli_list))]
        
        # define attributes
        self.pauli_strings = pauli_list
        self.weights = weights
        self.phases = phases
        self.n_paulis = len(pauli_list)
        self.n_qudits = pauli_list[0].n_qudits()
        # My guess is that speed will be better if we store the symplectic matrix here - could be a method instead
        self.symplectic = self.symplectic_matrix()

        if dimensions is None:
            if isinstance(pauli_list[0], PauliString):
                dimensions = np.array([pauli_list[i].dimensions for i in range(len(pauli_list))])
            else:
                dimensions = 2 * np.ones(self.n_qudits)
            dimensions = 2 * np.ones(self.n_qudits)

    def x(self):
        return self.symplectic[:, 1:self.n_qudits + 1]
    
    def z(self):
        return self.symplectic[:, self.n_qudits + 1:-1]
        
    def phases(self):
        return self.symplectic[:, -1]
    
    def weights(self):
        return self.symplectic[:, 0]
    
    def change_weight(self, index, weight):
        self.symplectic[index, 0] = weight

    def change_phase(self, index, phase):
        self.symplectic[index, -1] = phase

    def change_pauli(self, index, pauli):
        self.pauli_strings[index] = pauli
        self.symplectic[index, :] = pauli_to_symplectic(pauli, self.weights[index])
    
    def symplectic_structure_matrix(self):
        return self.symplectic[:, 1:-1]

    def __add__(self, symplectic_pauli):
        if isinstance(symplectic_pauli, PauliString) or isinstance(symplectic_pauli, Pauli):
            symplectic_pauli = symplectic_pauli.to_pauli_sum()
        new_pauli_list = self.pauli_strings + symplectic_pauli.pauli_strings
        new_weights = np.concatenate([self.weights, symplectic_pauli.weights])
        new_phases = np.concatenate([self.phases, symplectic_pauli.phases])
        return PauliSum(new_pauli_list, new_weights, new_phases)

    def __sub__(self, symplectic_pauli):
        new_pauli_list = self.pauli_strings + symplectic_pauli.pauli_string
        new_weights = np.concatenate([self.weights, -symplectic_pauli.weights])
        new_phases = np.concatenate([self.phases, symplectic_pauli.phases])
        return PauliSum(new_pauli_list, new_weights, new_phases)
    
    def __matmul__(self, symplectic_pauli):
        """
        @ is the operator for tensor product of two SymplecticPauli objects
        """
        new_pauli_list = []
        new_weights = []
        new_phases = []
        for i in range(self.n_paulis):
            for j in range(symplectic_pauli.n_paulis):
                next_pauli = ''
                next_weight = 1
                next_phase = 0
                for n in range(self.n_qudits):
                    next_pauli += self.pauli_strings[i][n] + symplectic_pauli.pauli_string[j][n]
                    next_phase += (self.phases()[i] + symplectic_pauli.phases()[j]) % 2
                    next_weight *= self.weights()[i] * symplectic_pauli.weights()[j]
                new_pauli_list.append(next_pauli)
                new_weights.append(next_weight)
                new_phases.append(((self.phases()[i] + symplectic_pauli.phases()[j] + next_phase) % 2))
        output_pauli = PauliSum(new_pauli_list, new_weights, new_phases)
        return output_pauli

    def __mul__(self, A):
        """
        Operator multiplication on two SymplecticPauli objects or multiplication of weights by constant
        """
        # Uses the Pauli relations
        # XX = ZZ = = YY = I
        # ZX = -XZ = Y with phase +1
        # YZ = -ZY = X with phase +1
        # XY = -YX = Z with phase +1

        # for qudits making the phase increase by 1/d and

        # there is probably a better way of doing this with the symplectic representation

        if isinstance(A, (int, float)):
            return PauliSum(self.pauli_strings, self.weights() * A, self.phases())
        elif not isinstance(A, PauliSum):
            raise ValueError("Multiplication only supported with SymplecticPauli objects or scalar")
        
        pauli_multiplication_rules = {
            ('I', 'I'): ('I', 0, 1),
            ('X', 'I'): ('X', 0, 1),
            ('Y', 'I'): ('Y', 0, 1),
            ('Z', 'I'): ('Z', 0, 1),
            ('I', 'X'): ('X', 0, 1),
            ('I', 'Y'): ('Y', 0, 1),
            ('I', 'Z'): ('Z', 0, 1),
            ('X', 'X'): ('I', 0, 1),
            ('Y', 'Y'): ('I', 0, 1),
            ('Z', 'Z'): ('I', 0, 1),
            ('X', 'Z'): ('Y', 1, -1),
            ('Z', 'X'): ('Y', 1, 1),
            ('Y', 'Z'): ('X', 1, 1),
            ('Z', 'Y'): ('X', 1, -1),
            ('X', 'Y'): ('Z', 1, 1),
            ('Y', 'X'): ('Z', 1, -1),
        }

        def multiply_pauli_string(ps1, ps2):
            if (ps1, ps2) in pauli_multiplication_rules:
                return pauli_multiplication_rules[(ps1, ps2)]
            else:
                raise Exception(f"Invalid Pauli strings {ps1} and {ps2}")
            
        # Runs through all the Pauli strings and weights, multiplies with the above rules
        # and returns the result as a new SymplecticPauli object
        weights = self.weights
        phases = self.phases
        new_pauli_list = []
        new_weights = []
        new_phases = []
        for i in range(self.n_paulis):
            for j in range(A.n_paulis):
                next_pauli = ''
                next_phase = 0
                next_weight = 1

                for n in range(self.n_qudits):
                    
                    # pp, phase, weight = multiply_pauli_string(self.pauli_strings[i][n], A.pauli_strings[j][n]) 
                    next_pauli += [self.pauli_strings[i] * A.pauli_strings[j]]
                    # next_phase += phase
                    next_weight *= weights[i] * A.weights[j]
                new_pauli_list.append(next_pauli)
                new_weights.append(next_weight * weights[i] * A.weights[j])
                new_phases.append(((phases[i] + A.phases[j] + next_phase) % 2))
        # clean up
        output_pauli = PauliSum(new_pauli_list, new_weights, new_phases)
        # output_pauli._remove_trivial_paulis()
        # output_pauli._combine_equivalent_paulis()
        return output_pauli
    
    def __truediv__(self, A):
        if not isinstance(A, (int, float)):
            raise ValueError("Division only supported with scalar")
        return self * (1 / A)
    
    def gate_operation(self, gate_symplectic):
        """
        Pre and post multiply by gate_symplectic

        """
        output = gate_symplectic * self * gate_symplectic
        output.combine_equivalent_paulis()
        output.remove_trivial_paulis()
        return output

    def combine_equivalent_paulis(self):
        # combine equivalent Paulis
        for i in range(self.n_paulis):
            for j in range(i + 1, self.n_paulis):
                if self.pauli_strings[i] == self.pauli_strings[j]:
                    self. change_weight(i, self.weights()[i] + self.weights()[j])
                    self.delete_paulis_(j)
                    break
        # remove zero weight Paulis
        to_delete = []
        for i in range(self.n_paulis):
            if self.weights()[i] == 0:
                to_delete.append(i)
        self.delete_paulis_(to_delete)

    def remove_trivial_paulis(self):
        # If entire Pauli string is I, remove it
        to_delete = []
        for i in range(self.n_paulis):
            if self.pauli_strings[i] == 'I' * self.n_qudits:
                to_delete.append(i)
        self.delete_paulis_(to_delete)

    def symplectic_matrix(self):
        symplectic = np.zeros([self.n_paulis, 2 * self.n_qudits])
        for i, p in enumerate(self.pauli_strings):
            symplectic[i, :] = p.symplectic()
        return symplectic

        # check whether self has only X component
    def is_IX(self):
        # Outputs:
        #     (bool) - True if self has only X component, False otherwise
        return not np.any(self.z())

    # check whether self has only Z component 
    def is_IZ(self):
        # Outputs:
        #     (bool) - True if self has only Z component, False otherwise
        return not np.any(self.x())

    # check whether the set of Paulis are pairwise commuting
    def is_commuting(self):
        # Outputs:
        #     (bool) - True if self is pairwise commuting set of Paulis
        spm = self.symplectic_product_matrix()
        return not np.any(spm)

    # pull out the ath Pauli from self
    def select_pauli(self, a):
        # Inputs:
        #     a - (int) - index of Pauli to be returned
        # Outputs:
        #     (pauli) - the ath Pauli in self
        return self.pauli_strings[a]

    # count the number of Paulis in self
    def paulis(self):
        # Output: (int)
        return self.X.shape[0]

    # count the number of qudits in self
    def qudits(self):
        # Outputs: (int)
        return self.X.shape[1]

    # delete Paulis indexed by aa
    def delete_paulis_(self, aa):
        # Inputs: 
        #     aa - (list of int)
        if type(aa) is int:
            aa = [aa]
        self.symplectic = np.delete(self.symplectic, aa, axis=0)
        self.pauli_strings = np.delete(self.pauli_strings, aa)
        self.n_paulis -= 1

    # return self after deletion of qudits indexed by aa
    def delete_qudits_(self, aa):
        # Inputs: 
        #     aa - (list of int)
        if type(aa) is int:
            self.X = np.delete(self.X, aa, axis=1)
            self.Z = np.delete(self.Z, aa, axis=1)
            self.dims = np.delete(self.dims, aa)
        else:
            for a in sorted(aa, reverse=True):
                self.X = np.delete(self.X, a, axis=1)
                self.Z = np.delete(self.Z, a, axis=1)
                self.dims = np.delete(self.dims, a)
        self.lcm = np.lcm.reduce(self.dims)   # was just dims... doesn't exist

    # return deep copy of self
    def copy(self):
        # Outputs: (SymplecticPauli)
        return PauliSum(self.X.copy(), self.Z.copy(), self.dims.copy(), self.phases.copy())

    def symplectic_product_matrix(self):
        """
        The symplectic product matrix, S, is an n x n matrix, n is the number
        of Paulis. The entry S[i, j] is the symplectic product of the
        ith Pauli and the jth Pauli.
        """
        n = self.n_paulis

        list_of_symplectics = self.symplectic_structure_matrix()

        # then we can create the full symplectic product matrix
        spm = np.zeros([n, n], dtype=int)
        for i in range(n):
            for j in range(n):
                if i > j:
                    spm[i, j] = symplectic_product(list_of_symplectics[i], list_of_symplectics[j])
        spm = spm + spm.T
        return spm
    
    def __str__(self):
        p_string = ''
        for i in range(self.n_paulis):
            pauli_string = self.pauli_strings[i]
            qudit_string = ''.join(['x' + f'{pauli_string.x_exp[j]}' + 'z' + f'{pauli_string.z_exp[j]} ' for j in range(self.n_qudits)])
            p_string += f'{self.weights[i]} |' + qudit_string + f'| {self.phases[i]} \n'
        return p_string


class Gate(PauliSum):
    def __init__(self, name, index, generalised_pauli_list, weights=None, phases=None, dims=None):
        self.name = name
        self.index = index
        super().__init__(generalised_pauli_list, weights, phases, dims)


class Hadamard(Gate):
    def __init__(self, n_qudits, index, dims):
        str1 = ''
        str2 = ''
        for i in range(n_qudits):
            if i == index:
                str1 += 'X'
                str2 += 'Z'
            else:
                str1 += 'I'
                str2 += 'I'
        pauli_string = [str1, str2]
        weights = 1. / np.sqrt(2) * np.ones(2)
        super().__init__(name='H', generalised_pauli_list=pauli_string, weights=weights, phases=None, dims=dims)


class CNOT(Gate):
    """
    |0><0|_control I_all + |1><1|_control X_target I_rest

    Uses:
        |0><0| = (I - Z) / 2
        |1><1| = (I + Z) / 2
    """

    def __init__(self, control, target, n_qubits, dims=2):
        # strings depend on state of the control - each state has two contributions I +- Z
        control_0_str1 = ''
        control_0_str2 = ''
        control_1_str1 = ''
        control_1_str2 = ''
        for i in range(n_qubits):
            if i == control:
                control_0_str1 += 'I'
                control_0_str2 += 'Z'
                control_1_str1 += 'I'
                control_1_str2 += 'Z'

            elif i == target:
                control_0_str1 += 'I'
                control_0_str2 += 'I'
                control_1_str1 += 'X'
                control_1_str2 += 'X'

            else:
                control_0_str1 += 'I'
                control_0_str2 += 'I'
                control_1_str1 += 'I'
                control_1_str2 += 'I'

        w01 = 1 / 2
        w02 = -1 / 2
        w11 = 1 / 2
        w12 = 1 / 2
        weights = [w01, w02, w11, w12]
        pauli_string = [control_0_str1, control_0_str2, control_1_str1, control_1_str2]
        super().__init__(name='CNOT', index=(control, target), generalised_pauli_list=pauli_string,
                         weights=weights, phases=None, dims=dims)

        
class Circuit:
    def __init__(self, gates, indexes):
        """
        Initialize the Circuit with gates, indexes, and targets.

        If a multi-qubit gate has a target, the targets should be at the ent of the tuple of indexes
        e.g. a CNOT with control 1, target 3 is

        gate = 'CNOT'
        indexes = (1, 3)
        

        Parameters:
            gates (list): A list of Gate objects representing the gates in the circuit.
            indexes (list): A list of integers or tuples for multi-qudit gates,
              indicating the indexes of qudits the gates act upon.
        """
        if len(gates) != len(indexes):
            raise ValueError(f"Length of gates ({len(gates)}) and indexes ({len(indexes)}) must be equal.")

        self.gates = gates
        self.indexes = indexes

    def add_gate(self, gate, indexes):
        """
        Appends a gate to qudit index with specified target (if relevant)
        """
        self.gates.append(gate)
        self.indexes.append(indexes)

    def remove_gate(self, index):
        """
        Removes a gate from the circuit at the specified index
        """
        self.gates.pop(index)
        self.indexes.pop(index)


def pauli_to_symplectic(pauli_string, weight):
    """
    The symplectic form of a pauli operator is
    I = [0, 0]
    X = [1, 0]
    Y = [1, 1]
    Z = [0, 1]
    :param pauli_string: list of I, X, Y, Z
    :return: symplectic vector
    """
    symplectic = np.zeros(2 * len(pauli_string) + 2)
    symplectic[0] = weight
    phase = 0
    for i, p in enumerate(pauli_string):
        if p == 'X':
            symplectic[i + 1] = 1
        if p == 'Z':
            symplectic[i + len(pauli_string) + 1] = 1
        if p == 'Y':
            symplectic[i + 1] = 1
            symplectic[i + len(pauli_string) + 1] = 1
            phase += 1
    symplectic[-1] = phase % 2
    return symplectic



def symplectic_to_pauli(symplectic):
    """
    Convert a symplectic vector to a Pauli string.

    :param symplectic: symplectic vector
    :return: Pauli string
    """
    pauli_string = ""
    n = (len(symplectic) - 2) // 2
    weight = symplectic[0]
    symplectic = symplectic[1:-1]
    phase = symplectic[-1]  # remove the weight
    for i in range(n):
        if symplectic[i] == 0 and symplectic[i + n] == 0:
            pauli_string += "I"
        elif symplectic[i] == 1 and symplectic[i + n] == 0:
            pauli_string += "X"
        elif symplectic[i] == 0 and symplectic[i + n] == 1:
            pauli_string += "Z"
        else:
            pauli_string += "Y"
    return pauli_string, weight, phase


def symplectic_product(p1, p2):
    """Takes the symplectic product of the binary vectors p1 and p2

    If 1 the corresponding Paulis anti-commute
    if 0 the corresponding Paulis commute"""
    s_prod = 0
    n = int(len(p1) / 2)
    for i in range(n):
        s_prod += (p1[i] * p2[i + n] - p1[i + n] * p2[i])
    return s_prod % 2


# the symplectic inner product of two pauli objects (each with a single Pauli)
def quditwise_inner_product(symplectic_pauli1, symplectic_pauli2):
    # Inputs:
    #     symplectic_pauli1 - (SymplecticPauli) - must have n_paulis = 1
    #     symplectic_pauli2 - (SymplecticPauli) - must have n_paulis = 1
    # Outputs:
    #     (bool) - quditwise inner product of Paulis
    if (symplectic_pauli1.n_paulis != 1) or (symplectic_pauli2.n_paulis != 1):
        raise Exception("Qubitwise inner product only works with pair of single Paulis")
    if any(symplectic_pauli1.dims - symplectic_pauli1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    return any(np.sum(symplectic_pauli1.X[0, i] * symplectic_pauli2.Z[0, i]
                      - symplectic_pauli2.Z[0, i] * symplectic_pauli1.X[0, i])
                      % symplectic_pauli1.dims[i] for i in range(symplectic_pauli1.qudits()))


if __name__ == "__main__":
    dims = 3
    x1 = Xnd(1, dims)
    y1 = Ynd(1, dims)
    z1 = Znd(1, dims)

    xx = x1 + y1 + z1

    print(xx)

    # xz = x1 @ z1
    # print(xz.x_exp)
    # print(xz.z_exp)
    # print(xz.dimensions)

    # print(xz)

    ps1_in = 'x1z0 x2z0 x1z0'
    ps2_in = 'x1z1 x2z2 x0z2'
    ps1 = PauliString(ps1_in, dimensions=[dims]*3)
    ps2 = PauliString(ps2_in, dimensions=[dims]*3)
    ps3 = PauliString([1, 2, 0], [1, 2, 2], dimensions=[dims]*3)

    # print(ps1.x_exp)
    # print(ps1.z_exp)
    # print(ps1)
    # print(ps2)
    # print(ps3)

    # Multiply Pauli strings

    print('          ' + str(ps1 * ps2))
    print('Should be x2z1 x1z2 x1z2')
    print('          ' + str(ps1 * ps3))
    print('Should be x2z1 x1z2 x1z2')
    print('          ' + str(ps2 * ps3))
    print('Should be x2z2 x1z1 x0z1')

    # Tensor product Pauli strings

    # print('          ' + str(ps1 @ ps2))
    # print('Should be x1z0 x2z0 x1z0 x1z1 x2z2 x0z2' + '\n')
    # print('          ' + str(ps1 @ ps3))
    # print('Should be x1z0 x2z0 x1z0 x1z1 x2z2 x0z2' + '\n')
    # print('          ' + str(ps2 @ ps3))
    # print('Should be x1z1 x2z2 x0z2 x1z1 x2z2 x0z2' + '\n')

    # Add Pauli strings

    ps12 = ps1 + ps2
    ps13 = ps1 + ps3
    
    # # Multiply Pauli sums

    ps1213 = ps12 * ps13

    # # Tensor product Pauli sums

    # ps12_13 = ps12 @ ps13

    # # Add Pauli sums

    # ps12_13_1213 = ps12_13 + ps1213

    # I_s = ['I']
    # X_s = ['X']
    # Y_s = ['Y']
    # Z_s = ['Z']

    # X = PauliSum(X_s)
    # Y = PauliSum(Y_s)
    # Z = PauliSum(Z_s)
    # I = PauliSum(I_s)

    # CNOT1 = (I - Z) / 2. @ I + (I + Z) / 2 @ X
    # CNOT2 = CNOT(0, 1, 2)

    # print(CNOT1)
