import sys
import numpy as np
import itertools

sys.path.append("./")
import quaos.gaussian_elimination as ge


class SymplecticPauli:
    def __init__(self, pauli_list, weights=None, phases=None):
        if not isinstance(pauli_list, list):
            pauli_list = [pauli_list]
        if weights is None:
            weights = np.ones(len(pauli_list))
        if phases is None:
            phases = np.zeros(len(pauli_list))
        if len(pauli_list) != len(weights):
            raise ValueError(f"Length of Pauli list ({len(pauli_list)}) and weights ({len(weights)}) must be equal.")

        self.n_paulis = len(pauli_list)
        self.n_qubits = len(pauli_list[0])
        self.pauli_string = pauli_list

        self.symplectic = self.symplectic_matrix(weights)

    def x(self):
        return self.symplectic[:, 1:self.n_qubits + 1]
    
    def z(self):
        return self.symplectic[:, self.n_qubits + 1:-1]
    
    def phases(self):
        return self.symplectic[:, -1]
    
    def weights(self):
        return self.symplectic[:, 0]
    
    def change_weight(self, index, weight):
        self.symplectic[index, 0] = weight

    def change_phase(self, index, phase):
        self.symplectic[index, -1] = phase
    
    def change_pauli(self, index, pauli):
        self.pauli_string[index] = pauli
        self.symplectic[index, :] = pauli_to_symplectic(pauli, self.weights[index])
    
    def symplectic_structure_matrix(self):
        return self.symplectic[:, 1:-1]

    def __add__(self, symplectic_pauli):
        new_pauli_list = self.pauli_string + symplectic_pauli.pauli_string
        new_weights = np.concatenate([self.weights(), symplectic_pauli.weights()])
        new_phases = np.concatenate([self.phases(), symplectic_pauli.phases()])
        return SymplecticPauli(new_pauli_list, new_weights, new_phases)

    def __sub__(self, symplectic_pauli):
        new_pauli_list = self.pauli_string + symplectic_pauli.pauli_string
        new_weights = np.concatenate([self.weights(), -symplectic_pauli.weights()])
        new_phases = np.concatenate([self.phases(), symplectic_pauli.phases()])
        return SymplecticPauli(new_pauli_list, new_weights, new_phases)
    
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
                for n in range(self.n_qubits):
                    next_pauli += self.pauli_string[i][n] + symplectic_pauli.pauli_string[j][n]
                    next_phase += (self.phases()[i] + symplectic_pauli.phases()[j]) % 2
                    next_weight *= self.weights()[i] * symplectic_pauli.weights()[j]
                new_pauli_list.append(next_pauli)
                new_weights.append(next_weight)
                new_phases.append(((self.phases()[i] + symplectic_pauli.phases()[j] + next_phase) % 2))
        output_pauli = SymplecticPauli(new_pauli_list, new_weights, new_phases)
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

        # there is probably a better way of doing this with the symplectic representation

        if isinstance(A, (int, float)):
            return SymplecticPauli(self.pauli_string, self.weights() * A, self.phases())
        elif not isinstance(A, SymplecticPauli):
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
        weights = self.weights()
        phases = self.phases()
        new_pauli_list = []
        new_weights = []
        new_phases = []
        for i in range(self.n_paulis):
            for j in range(A.n_paulis):
                next_pauli = ''
                next_phase = 0
                next_weight = 1
                for n in range(self.n_qubits):
                    pp, phase, weight = multiply_pauli_string(self.pauli_string[i][n], A.pauli_string[j][n]) 
                    next_pauli += pp
                    next_phase += phase
                    next_weight *= weight
                new_pauli_list.append(next_pauli)
                new_weights.append(next_weight * weights[i] * A.weights()[j])
                new_phases.append(((phases[i] + A.phases()[j] + next_phase) % 2))
        # clean up
        output_pauli = SymplecticPauli(new_pauli_list, new_weights, new_phases)
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
                if self.pauli_string[i] == self.pauli_string[j]:
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
            if self.pauli_string[i] == 'I' * self.n_qubits:
                to_delete.append(i)
        self.delete_paulis_(to_delete)

    def symplectic_matrix(self, weights=None):
        symplectic = np.zeros([self.n_paulis, 2 * self.n_qubits + 2])
        for i, p in enumerate(self.pauli_string):
            symplectic[i, :] = pauli_to_symplectic(p, weights[i])
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
        symplectic = self.symplectic[a, :]
        return symplectic_to_pauli(symplectic)

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
        self.pauli_string = np.delete(self.pauli_string, aa)
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
        return SymplecticPauli(self.X.copy(), self.Z.copy(), self.dims.copy(), self.phases.copy())

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
    
    def diagonalize(self):
        pass

    def reduced_form(self, return_copy=False):
        # should be altered to track weights as well
        """
        Reduce a list of Paulis to a minimal set of generators.

        Takes the full symplectic product matrix of the Paulis and performs
        Gaussian elimination to determine the rank of the matrix. The rank
        determines the number of generators required to generate the original
        set of Paulis. The generators are then determined by applying the
        inverse of the lower triangular matrix to the original symplectic
        vectors. The weights of the generators are determined by the weights of
        the original Paulis.

        Args:
            return_copy (bool): If True, return new instance of SymplecticPauli
                with the reduced form, otherwise modifies the current instance.

        Returns:
            If return_copy is True, a new instance of SymplecticPauli with the
            reduced form. Otherwise, None.

        ################# Future  #################
            
        SHOULD ALSO OUTPUT THE CIRCUIT WHICH REDUCES THE PAULIS

        ###########################################
        """

        spm = self.symplectic_product_matrix()
        L, rank = ge.symmetric_gaussian_elimination(spm)
        dim = len(spm)
        registers = dim - rank // 2
        fund_p = np.zeros((dim, 2 * registers), dtype=int)
        for i in range(dim - rank):
            fund_p[i + 2 * rank // 2][i + registers + rank // 2] = 1
        for i in range(rank // 2):
            fund_p[2 * i][i] = 1
            fund_p[2 * i + 1][i + registers] = 1
        L_inv = np.linalg.inv(L) % 2
        minimal_rep = np.zeros((dim, 2 * registers))
        for i in range(dim):
            for j in range(dim):
                minimal_rep[i] += (L_inv[i, j] * fund_p[j]) % 2

        pauli_forms = [self.symplectic_to_pauli(minimal_rep[i])
                       for i in range(dim)]
        if return_copy:
            return SymplecticPauli(pauli_forms, self.weights)
        else:
            # update current instance to reduced form
            self.pauli_string = pauli_forms
            self.n_paulis = len(pauli_forms)
            self.symplectic = self.symplectic_matrix()
            self.symplectic_structure = self.symplectic_structure_matrix()
            self.phases = self.symplectic[-1, :]


class SymplecticQuditPauli(SymplecticPauli):
    def __init__(self, generalised_pauli_list, weights=None, phases=None, dims=None):

        if dims is None:
            dims = 2 * np.ones(len(generalised_pauli_list), dtype=int)
        if len(generalised_pauli_list) != len(weights):
            raise ValueError("Length of Pauli list and weights must be equal.")

        self.dims = dims
        self.lcm = np.lcm.reduce(dims)
        super().__init__(generalised_pauli_list, weights, phases)

    # check whether the set of Paulis are pairwise commuting on every qudit
    def is_quditwise_commuting(self):
        # Outputs:
        #     (bool) - True if self is pairwise quditwise commuting set of Paulis
        p = self.paulis()
        PP = [self.a_pauli(i) for i in range(p)]
        return not any(quditwise_inner_product(PP[i0], PP[i1]) for i0, i1 in itertools.combinations(range(p), 2))
    

class Gate(SymplecticQuditPauli):
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


def symplectic_product(p1, p2):
    """Takes the symplectic product of the binary vectors p1 and p2

    If 1 the corresponding Paulis anti-commute
    if 0 the corresponding Paulis commute"""
    s_prod = 0
    n = int(len(p1) / 2)
    for i in range(n):
        s_prod += (p1[i] * p2[i + n] - p1[i + n] * p2[i])
    return s_prod % 2


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


# the symplectic inner product of two pauli objects (each with a single Pauli)
def quditwise_inner_product(P0, P1):
    # Inputs:
    #     P0 - (pauli) - must have shape (1,q)
    #     P1 - (pauli) - must have shape (1,q)
    # Outputs:
    #     (bool) - quditwise inner product of Paulis
    if (P0.paulis() != 1) or (P1.paulis() != 1):
        raise Exception("Qubitwise inner product only works with pair of single Paulis")
    if any(P0.dims - P1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    return any(np.sum(P0.X[0, i] * P1.Z[0, i] - P0.Z[0, i] * P1.X[0, i]) % P0.dims[i] for i in range(P0.qudits()))


if __name__ == "__main__":

    I_s = ['I']
    X_s = ['X']
    Y_s = ['Y']
    Z_s = ['Z']

    X = SymplecticPauli(X_s)
    Y = SymplecticPauli(Y_s)
    Z = SymplecticPauli(Z_s)
    I = SymplecticPauli(I_s)

    CNOT1 = (I - Z) / 2. @ I + (I + Z) / 2 @ X
    CNOT2 = CNOT(0, 1, 2)

    print(CNOT1.pauli_string)
    print(CNOT2.pauli_string)
    print(CNOT1.weights())
    print(CNOT2.weights())
    print(CNOT1.phases())
    print(CNOT2.phases())
