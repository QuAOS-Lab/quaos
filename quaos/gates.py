import sys
import numpy as np
from qiskit import QuantumCircuit
sys.path.append("./")

import time
from quaos.symplectic import PauliSum, PauliString, Pauli, Xnd, Ynd, Znd, Id, string_to_symplectic, symplectic_to_string


# class Gate(PauliSum):
#     def __init__(self, name, index, generalised_pauli_list, weights=None, phases=None, dims=None):
#         self.name = name
#         self.index = index
#         super().__init__(generalised_pauli_list, weights, phases, dims)

#     def act(self, P):
#         return self * P * self


# class Hadamard(Gate):
#     def __init__(self, n_qudits, index, dims):
#         str1 = ''
#         str2 = ''
#         for i in range(n_qudits):
#             if i == index:
#                 str1 += 'X'
#                 str2 += 'Z'
#             else:
#                 str1 += 'I'
#                 str2 += 'I'
#         pauli_string = [str1, str2]
#         weights = 1. / np.sqrt(2) * np.ones(2)
#         super().__init__(name='H', generalised_pauli_list=pauli_string, weights=weights, phases=None, dims=dims)


# class CNOT(Gate):
#     """
#     |0><0|_control I_all + |1><1|_control X_target I_rest

#     Uses:
#         |0><0| = (I - Z) / 2
#         |1><1| = (I + Z) / 2
#     """

#     def __init__(self, control, target, n_qubits):
#         dims = [2] * n_qubits
#         # strings depend on state of the control - each state has two contributions I +- Z
#         control_0_str1 = ''
#         control_0_str2 = ''
#         control_1_str1 = ''
#         control_1_str2 = ''
#         for i in range(n_qubits):
#             if i == control:
#                 control_0_str1 += 'x0z0'
#                 control_0_str2 += 'x0z1'
#                 control_1_str1 += 'x0z0'
#                 control_1_str2 += 'x0z1'

#             elif i == target:
#                 control_0_str1 += 'x0z0'
#                 control_0_str2 += 'x0z0'
#                 control_1_str1 += 'x1z0'
#                 control_1_str2 += 'x1z0'

#             else:
#                 control_0_str1 += 'x0z0'
#                 control_0_str2 += 'x0z0'
#                 control_1_str1 += 'x0z0'
#                 control_1_str2 += 'x0z0'

#         w01 = 1 / 2
#         w02 = -1 / 2
#         w11 = 1 / 2
#         w12 = 1 / 2
#         weights = [w01, w02, w11, w12]
#         pauli_string = [control_0_str1, control_0_str2, control_1_str1, control_1_str2]
#         super().__init__(name='CNOT', index=(control, target), generalised_pauli_list=pauli_string,
#                          weights=weights, phases=None, dims=dims)


# class SUM(Gate):
#     def __init__(self, index, generalised_pauli_list, weights=None, phases=None, dims=None):
#         super().__init__('SUM', index, generalised_pauli_list, weights, phases, dims)

        
class GateOperation:
    """
    Mapping can be written as set of rules,
    
    e.g. for CNOT(control=3, target=1, dimensions=2)
                x1z0*x0z0 -> x1z0*x1z0
                x0z0*x1z0 -> x0z0*x1z0  # doesn't need specifying
                x0z1*x0z0 -> x0z1*x0z0  # doesn't need specifying
                x0z0*x0z1 -> -x0z1*x0z1 # note phase important here

    inputs are: 

    qudit_indices = (3, 1)
    mapping = ['x1z0*x0z0 -> x1z0*x1z0', 'x0z0*x0z1 -> -x0z1*x0z1']  # (control*target -> control*target)

    or for SUM(target=1, control=3)

    qudit_indices = (3, 1)
    mapping = ['X*I -> X*X', 'I*Z -> Z-1*Z']   # Interpreting Xn = X^n, Z-n = Z^{-n}
    
    on all mappings (output) % lcm is assumed

    Should act on a PauliSum to return another PauliSum

    """
    def __init__(self, name, qudit_indices, mapping, dimension):
        
        self.dimension = dimension
        self.name = name
        self.qudit_indices = qudit_indices  # number of total qudits including those not acted upon - could remove the need for this entirely...
        self.map_from, self.map_to, self.acquired_phase = self._interpret_mapping(mapping)
    
    def _interpret_mapping(self, map_string):

        map_from, map_to = zip(*[map_string[i].split('->') for i in range(len(map_string))])

        n_maps = len(map_from)

        symplectic_looked_for = []
        for i in range(n_maps):
            s, _ = string_to_symplectic(map_from[i])
            symplectic_looked_for.append(s)

        symplectic_mapped_to = []
        acquired_phase = []
        for i in range(len(map_from)):
            s, p = string_to_symplectic(map_to[i])
            symplectic_mapped_to.append(s)
            acquired_phase.append(p)
        
        # For an n qubit operation there should be a set of 2n mappings.
        # If the rest are unspecified, they should be identity mappings (those with equal looked_for and mapped_to)
        # If the user specifies less than 2n mappings, then the remaining ones will be identity mappings
        # These remaining identity mappings must all be linearly independent of the specified ones
        # note that this will fail if the input mappings are not of the form X*I*...*I, Z*I*...*I, (having only a single
        # X or Z in the string - checked for below) - generalisable if needed
        # This will only be needed to obtain the symplectic matrix for the gate operation, not for current method...
        # if len(symplectic_looked_for) < 2 * self.n_qudits:
        #     for s in symplectic_looked_for:
        #         if np.sum(s) != 1:
        #             raise Exception('Unable to automate completion of mappings. Specify 2*n_qudits linearly '
        #                             'independent mappings to fully define the gate operation.')
        #     for i in range(2 * self.n_qudits):
        #         mapping = np.zeros(2 * self.n_qudits, dtype=int)
        #         mapping[i] = 1
        #         if not any(np.array_equal(mapping, arr) for arr in symplectic_looked_for):
        #             symplectic_looked_for.append(mapping)
        #             symplectic_mapped_to.append(mapping)

        ### remove once debugged
        # assert len(symplectic_looked_for) == self.dimension ** (self.n_qudits), (len(symplectic_looked_for),
        #                                                                              self.dimension ** (self.n_qudits))
        # assert len(symplectic_mapped_to) == self.dimension ** (self.n_qudits)
        ###

        symplectic_looked_for = np.array(symplectic_looked_for)
        symplectic_mapped_to = np.array(symplectic_mapped_to)

        # reorder such that the symplectic looked for is always the identity
        # (this would need to be altered to generalise)
        # perm = np.argmax(symplectic_looked_for, axis=1)
        # inverse_perm = np.argsort(perm)
        # symplectic_looked_for = symplectic_looked_for[inverse_perm]
        # symplectic_mapped_to = symplectic_mapped_to[inverse_perm]

        # symplectic = self.build_symplectic_from_mappings(symplectic_looked_for, symplectic_mapped_to)

        return symplectic_looked_for, symplectic_mapped_to, acquired_phase
    
    # def mod_inv(self, a):
    #     """Modular inverse of a mod d, where d = self.dimension."""
    #     a = a % self.dimension
    #     if a == 0:
    #         raise ValueError("0 has no inverse modulo d")
    #     return pow(int(a), -1, self.dimension)

    # def mod_mat_inv(self, A):
    #     """Modular inverse of matrix A over Z_d using Gauss-Jordan elimination."""
    #     A = np.array(A, dtype=int) % self.dimension
    #     n = A.shape[0]
    #     I = np.eye(n, dtype=int)
    #     AI = np.hstack([A, I])  # Augmented matrix [A | I]

    #     for i in range(n):
    #         # Find pivot
    #         pivot = AI[i, i]
    #         if pivot == 0:
    #             # Try to swap with a lower row
    #             for j in range(i + 1, n):
    #                 if AI[j, i] != 0:
    #                     AI[[i, j]] = AI[[j, i]]
    #                     pivot = AI[i, i]
    #                     break
    #             else:
    #                 raise ValueError("Matrix not invertible")

    #         # Normalize pivot row
    #         inv_pivot = self.mod_inv(pivot)
    #         AI[i] = (AI[i] * inv_pivot) % self.dimension

    #         # Eliminate other rows
    #         for j in range(n):
    #             if j != i:
    #                 factor = AI[j, i]
    #                 AI[j] = (AI[j] - factor * AI[i]) % self.dimension

    #     A_inv = AI[:, n:]  # Extract right half
    #     return A_inv

    # def build_symplectic_from_mappings(self, looked_for, mapped_to):
    #     """Build symplectic matrix from Pauli generator mappings."""
    #     looked_for = np.array(looked_for, dtype=int) % self.dimension
    #     mapped_to = np.array(mapped_to, dtype=int) % self.dimension

    #     if looked_for.shape != mapped_to.shape:
    #         raise ValueError("Shape mismatch between looked_for and mapped_to.")
    #     # if looked_for.shape[0] != looked_for.shape[1]:
    #     #     raise ValueError("Expect square matrix of 2n symplectic vectors.", looked_for)

    #     inv_looked_for = self.mod_mat_inv(looked_for)
    #     S = (mapped_to @ inv_looked_for) % self.dimension
    #     return S

    def _act_on_pauli_string(self, P, return_phase=False):
        # Extract symplectic of PauliString
        # Extract the symplectic of the relevant qudit numbers in self.qudit_indices
        # Check if these correspond to any of the self.symplectics_looked_for
        # If so, replace with symplectic_mapped_to
        # Replace the mapped symplectics in the positions in self.qudit_indices
        # Return new PauliString
        # Phase ignored by default as it is a global phase for only a single PauliString
        symplectic = P.symplectic()
        local_symplectic = np.zeros(2 * len(self.qudit_indices))
        for i, index in enumerate(self.qudit_indices):
            local_symplectic[i] = symplectic[index]
            local_symplectic[i + len(self.qudit_indices)] = symplectic[index + P.n_qudits()]

        acquired_phase = 0
        for i, symplectic_looked_for in enumerate(self.map_from):
            symplectic_mapped_to = self.map_to[i]
            if (local_symplectic == symplectic_looked_for).all():
                local_symplectic = symplectic_mapped_to
                # print(self.acquired_phase[i])
                acquired_phase += self.acquired_phase[i]
                break

        P = P._replace_symplectic(local_symplectic, self.qudit_indices)
        if return_phase:
            return P, acquired_phase
        else:
            return P

    def _act_on_pauli_sum(self, P):
        if np.all(self.acquired_phase == 0):

            # In this case each output of _act_on_pauli_string is a PauliString
            P = PauliSum([self._act_on_pauli_string(p) for p in P.pauli_strings], P.weights, P.phases, P.dims, False)
        else:
            # in this case at least one is a PauliSum
            pauli_list = []
            acquired_phase = []
            for p in P.pauli_strings:
                new_pauli_string, additional_phase = self._act_on_pauli_string(p, return_phase=True)
                acquired_phase.append(additional_phase)
                pauli_list.append(new_pauli_string)
            weights = P.weights
            P = PauliSum(pauli_list, weights, P.phases, P.dimensions, False)
            P.acquire_phase(acquired_phase)
        # P.combine_equivalent_paulis()
        return P

    def act(self, P):
        if isinstance(P, Pauli):
            P = P.to_pauli_string()
        if isinstance(P, PauliString):
            return self._act_on_pauli_string(P)
        elif isinstance(P, PauliSum):
            return self._act_on_pauli_sum(P)
        else:
            raise ValueError(f"TwoQuditOperation cannot act on type {type(P)}")
        
    def __mul__(self, gate):
        circuit = Circuit([self + gate])
        return circuit
    

class CNOT(GateOperation):
    def __init__(self, control, target, n_qudits):
        CNOT_operations = ['x1z0 x0z0 -> x1z0 x1z0', 'x0z0 x0z1 -> x0z1 x0z1']
        super().__init__("CNOT", [control, target], CNOT_operations, n_qudits=n_qudits, dimension=[2, 2])
    

class Hadamard(GateOperation):
    def __init__(self, index, dimension):
        Hadamard_operations = self.hadamard_gate_operations(dimension)
        super().__init__("H", [index], Hadamard_operations, dimension=[dimension])

    @staticmethod
    def hadamard_gate_operations(dimension):
        operations = []
        for r in range(dimension):
            for s in range(dimension):
                operations.append(f"x{r}z{s} -> x{-s % dimension}z{r}p{r * s}")
        return operations
    

class PHASE(GateOperation):
    def __init__(self, index, dimension):
        SGate_operations = self.s_gate_operations(dimension)  # ['x1z0 -> x1z1', 'x0z1 -> x1z0']
        super().__init__("S", [index], SGate_operations, dimension=dimension)

    @staticmethod
    def s_gate_operations(dimension):
        operations = []
        for r in range(dimension):
            for s in range(dimension):
                operations.append(f"x{r}z{s} -> x{r}z{s + r}p{r * (r - 1) // 2}")
        return operations


class SUM(GateOperation):
    def __init__(self, control, target, dimension):
        SGate_operations = self.sum_gate_operations(dimension)
        super().__init__("SUM", [control, target], SGate_operations, dimension=dimension)
   
    @staticmethod
    def sum_gate_operations(dimension):
        operations = []
        for r1 in range(dimension):
            for s1 in range(dimension):
                for r2 in range(dimension):
                    for s2 in range(dimension):
                        new_r1 = r1
                        new_s1 = (s1 - s2) % dimension
                        new_r2 = (r2 + r1) % dimension
                        new_s2 = s2
                        phase = (r1 * s2) % dimension
                        operations.append(f"x{r1}z{s1} x{r2}z{s2} -> x{new_r1}z{new_s1} x{new_r2}z{new_s2}p{phase}")
        return operations


class Circuit:
    def __init__(self, dimensions, gates=None):
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
        if gates is None:
            gates = []
        self.dimensions = dimensions
        self.gates = gates
        self.indexes = [gate.qudit_indices for gate in gates]  # indexes accessible at the Circuit level

    def add_gate(self, gate):
        """
        Appends a gate to qudit index with specified target (if relevant)

        If gate is a list indexes should be a list of integers or tuples
        """
        if isinstance(gate, list) or isinstance(gate, np.ndarray):
            for i, g in enumerate(gate):
                self.gates.append(g)
                self.indexes.append(g.qudit_indices)
        else:
            self.gates.append(gate)
            self.indexes.append(gate.qudit_indices)

    def remove_gate(self, index):
        """
        Removes a gate from the circuit at the specified index
        """
        self.gates.pop(index)
        self.indexes.pop(index)

    def __str__(self):
        str_out = ''
        for gate in self.gates:
            str_out += gate.name + ' ' + str(gate.qudit_indices) + '\n'
        return str_out

    def act(self, pauli):
        for gate in self.gates:
            pauli = gate.act(pauli)
        return pauli
    
    def show(self):
        circuit = QuantumCircuit(len(self.dimensions))
        dict = {'X': circuit.x, 'H': circuit.h, 'S': circuit.s, 'SUM': circuit.cx, 'CNOT': circuit.cx}

        for gate in self.gates:
            name = gate.name
            if len(gate.qudit_indices) == 2:
                dict[name](gate.qudit_indices[0], gate.qudit_indices[1])
            else:
                dict[name](gate.qudit_indices[0])

        print(circuit)
        return circuit


if __name__ == "__main__":
    from quaos.hamiltonian import cancel_X, random_pauli_hamiltonian
    import random
    random.seed(27)

    # # testing gates

    # X = Xnd(1, 2)
    # Y = Ynd(1, 2)
    # Z = Znd(1, 2)
    # I = Id(2)

    # # print((I + Z) @ X)
    # # print(I @ X + Z @ X)
    
    # CNOT1 = ((I - Z) @ I + (I + Z) @ X) / 2.
    # CNOT2 = CNOT(0, 1, 2)
    # print(CNOT1.symplectic_matrix())
    # print(CNOT2.symplectic_matrix())

    # CNOT on two qubits

    # CNOT_operations = ['x1z0 x0z0 -> x1z0 x1z0', 'x0z0 x0z1 -> x0z1 x0z1']
    # CNOT3 = GateOperation('CNOT', [0, 1], CNOT_operations, 2)

    # ps1 = PauliString('x1z0 x0z0', dimensions=[2, 2])
    # ps2 = PauliString('x0z0 x0z1', dimensions=[2, 2])
    # ps3 = PauliString('x1z0 x1z0', dimensions=[2, 2])
    # ps4 = PauliString('x0z1 x0z1', dimensions=[2, 2])
    # ps5 = PauliString('x1z1 x0z0', dimensions=[2, 2])

    # print(ps1, '->', CNOT3.act(ps1))
    # print(ps2, '->', CNOT3.act(ps2))
    # print(ps3, '->', CNOT3.act(ps3))
    # print(ps4, '->', CNOT3.act(ps4))

    # psum1 = ps1 + 0.5 * ps2 + 1j * ps3 - 0.5j * ps4
    # print(psum1, '\n -> \n', CNOT3.act(psum1))

    # Hg = Hadamard(0, 2)
    # print(ps1, '->', H.act(ps1))
    # print(ps2, '->', H.act(ps2))
    # print(ps3, '->', H.act(ps3))
    # print(ps4, '->', Hg.act(ps4))
    # print(ps5, '->', Hg.act(ps5))

    ps = random_pauli_hamiltonian(5, [2] * 5, mode='uniform')
    Sum01 = SUM(0, 1, 2)

    print(Sum01.act(ps))
    # print(Sum03.act(Sum02.act(Sum01.act(ps))))RP
    # c = Circuit([2 * 5])
    # ps2, c = cancel_X(ps, 0, 5, c, 5)

    # print(ps2)
