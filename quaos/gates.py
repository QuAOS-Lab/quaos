import sys
import numpy as np
import itertools
import re
import functools

sys.path.append("./")

from quaos.symplectic import PauliSum, PauliString, Pauli, Xnd, Ynd, Znd, Id, string_to_symplectic, symplectic_to_string


class Gate(PauliSum):
    def __init__(self, name, index, generalised_pauli_list, weights=None, phases=None, dims=None):
        self.name = name
        self.index = index
        super().__init__(generalised_pauli_list, weights, phases, dims)

    def act(self, P):
        return self * P * self


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

    def __init__(self, control, target, n_qubits):
        dims = [2] * n_qubits
        # strings depend on state of the control - each state has two contributions I +- Z
        control_0_str1 = ''
        control_0_str2 = ''
        control_1_str1 = ''
        control_1_str2 = ''
        for i in range(n_qubits):
            if i == control:
                control_0_str1 += 'x0z0'
                control_0_str2 += 'x0z1'
                control_1_str1 += 'x0z0'
                control_1_str2 += 'x0z1'

            elif i == target:
                control_0_str1 += 'x0z0'
                control_0_str2 += 'x0z0'
                control_1_str1 += 'x1z0'
                control_1_str2 += 'x1z0'

            else:
                control_0_str1 += 'x0z0'
                control_0_str2 += 'x0z0'
                control_1_str1 += 'x0z0'
                control_1_str2 += 'x0z0'

        w01 = 1 / 2
        w02 = -1 / 2
        w11 = 1 / 2
        w12 = 1 / 2
        weights = [w01, w02, w11, w12]
        pauli_string = [control_0_str1, control_0_str2, control_1_str1, control_1_str2]
        super().__init__(name='CNOT', index=(control, target), generalised_pauli_list=pauli_string,
                         weights=weights, phases=None, dims=dims)


class SUM(Gate):
    def __init__(self, index, generalised_pauli_list, weights=None, phases=None, dims=None):
        super().__init__('SUM', index, generalised_pauli_list, weights, phases, dims)

        
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
    def __init__(self, name, qudit_indices, mapping, n_qudits):
        if len(qudit_indices) > n_qudits:
            raise ValueError("Number of qudits acted upon by GateOperation larger than number of qudits")
        self.name = name
        self.qudit_indices = qudit_indices
        self.map_from, self.map_to, self.specific_map = self._interpret_mapping(mapping)

    def _str_to_symplectic(self, string):
        if string == 'I':
            return (0, 0)
        elif string == 'X':
            return (1, 0)
        elif string == 'Z':
            return (0, 1)
    
    def _interpret_string(self, string):
        string = string.split('*')
        symplectic_list = []
        for i, qudit_state in enumerate(string):
            qudit_index = self.qudit_indices[i]
            # needs something here to extract phase if present and include it
            symplectic = string_to_symplectic(qudit_state, self.n_qudits)
            symplectic_list.append((symplectic, qudit_index))
        return symplectic_list
    
    def _interpret_mapping(self, map_string):
        map_string = "".join(map_string.split())  # remove spaces
        map_from, map_to = map_string.split('->')
        map_from = map_from.split(',')
        map_to = map_to.split(',')
        if len(map_from) != len(map_to):
            raise ValueError("Length of map_from and map_to must be equal")
        n_maps = len(map_from)

        symplectic_looked_for = []
        for i in range(n_maps):
            s = np.zeros(2 * self.n_qudits)
            


        symplectic_mapped_to = []
        for i in range(len(map_from)):
            map_from[i] = map_from[i].split('*')
        return map_from, map_to

    def act(self, P):
        if isinstance(P, PauliString):
            return self._act_on_pauli_string(P)
        elif isinstance(P, PauliSum):
            return self._act_on_pauli_sum(P)
        else:
            raise ValueError(f"TwoQuditOperation cannot act on type {type(P)}")
            


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

    def show(self):
        # uses scipy to print circuit
        pass


if __name__ == "__main__":
    # testing gates

    X = Xnd(1, 2)
    Y = Ynd(1, 2)
    Z = Znd(1, 2)
    I = Id(2)

    # print((I + Z) @ X)
    # print(I @ X + Z @ X)
    
    CNOT1 = ((I - Z) @ I + (I + Z) @ X) / 2.
    CNOT2 = CNOT(0, 1, 2)
    print(CNOT1)
    print(CNOT2)
