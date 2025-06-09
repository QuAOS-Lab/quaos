from quaos.gates import GateOperation
import numpy as np
from qiskit import QuantumCircuit
from quaos.paulis import PauliSum, PauliString, Pauli

# TODO: Replace GateOperation with Gate in quaos.circuits


class Circuit:
    def __init__(self, dimensions: list[int] | np.ndarray, gates: list[GateOperation] | None = None):
        """
        Initialize the Circuit with gates, indexes, and targets.

        If a multi-qubit gate has a target, the targets should be at the ent of the tuple of indexes
        e.g. a CNOT with control 1, target 3 is

        gate = 'CNOT'
        indexes = (1, 3)
        

        Parameters:
            dimensions (list[int] | np.ndarray): A list or array of integers representing the dimensions of the qudits.
            gates (list): A list of Gate objects representing the gates in the circuit.
            
        """
        if gates is None:
            gates = []
        self.dimensions = dimensions
        self.gates = gates
        self.indexes = [gate.qudit_indices for gate in gates]  # indexes accessible at the Circuit level

    def add_gate(self, gate: GateOperation | list[GateOperation]):
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

    def remove_gate(self, index: int):
        """
        Removes a gate from the circuit at the specified index
        """
        self.gates.pop(index)
        self.indexes.pop(index)

    def n_qudits(self) -> int:
        """
        Returns the number of qudits in the circuit.
        """
        return len(self.dimensions)

    def __add__(self, other: 'Circuit') -> 'Circuit':
        """
        Adds two circuits together by concatenating their gates and indexes.
        """
        if not isinstance(other, Circuit):
            raise TypeError("Can only add another Circuit object.")
        new_gates = self.gates + other.gates
        return Circuit(self.dimensions, new_gates)
    
    def __mul__(self, other: 'Circuit') -> 'Circuit':
        """
        THIS IS THE SAME FUNCTION AS ADDITION  -  PROBABLY WANT TO CHOOSE WHICH ONE DOES THIS

        Adds two circuits together by concatenating their gates and indexes.
        """
        if not isinstance(other, Circuit):
            raise TypeError("Can only add another Circuit object.")
        new_gates = self.gates + other.gates
        return Circuit(self.dimensions, new_gates)
    
    def __eq__(self, other: 'Circuit') -> bool:
        if not isinstance(other, Circuit):
            return False
        if len(self.gates) != len(other.gates):
            return False
        for i in range(len(self.gates)):
            if self.gates[i] != other.gates[i]:
                return False
        return True
    
    def __getitem__(self, index: int) -> GateOperation:
        return self.gates[index]
    
    def __setitem__(self, index: int, value: GateOperation):
        self.gates[index] = value
        self.indexes[index] = value.qudit_indices

    def __len__(self) -> int:
        return len(self.gates)

    def __str__(self) -> str:
        str_out = ''
        for gate in self.gates:
            str_out += gate.name + ' ' + str(gate.qudit_indices) + '\n'
        return str_out

    def act(self, pauli: Pauli | PauliString | PauliSum) -> PauliString | PauliSum:
        if isinstance(pauli, Pauli):
            if self.dimensions[0] != pauli.dimension or len(self.dimensions) != 1:
                raise ValueError("Pauli dimension does not match circuit dimensions")
        elif np.any(self.dimensions != pauli.dimensions):
            raise ValueError("Pauli dimensions do not match circuit dimensions")
        for gate in self.gates:
            pauli = gate.act(pauli)
        return pauli   # Why this warning? - it is a PauliSum or PauliString, never a Pauli

    def show(self) -> QuantumCircuit:
        if np.all(np.array(self.dimensions) == 2):
            print("Circuit dimensions are all 2, using Qiskit QuantumCircuit, some gates may not be supported")
        circuit = QuantumCircuit(len(self.dimensions))
        dict = {'X': circuit.x, 'H': circuit.h, 'S': circuit.s, 'SUM': circuit.cx, 'CNOT': circuit.cx,
                'Hdag': circuit.h}

        for gate in self.gates:
            name = gate.name
            if len(gate.qudit_indices) == 2:
                dict[name](gate.qudit_indices[0], gate.qudit_indices[1])
            else:
                dict[name](gate.qudit_indices[0])

        print(circuit)
        return circuit
    
    def copy(self) -> 'Circuit':
        return Circuit(self.dimensions, self.gates.copy())

    def embed_circuit(self, circuit: 'Circuit', qudit_indices: list[int] | np.ndarray | None = None):
        """
        Embed a circuit into current circuit at the specified qudit indices.
        """

        if qudit_indices is not None:
            if len(qudit_indices) != circuit.n_qudits():
                raise ValueError("Number of qudit indices does not match number of qudits in circuit to embed")
            
        for gate in circuit.gates:
            new_gate = gate.copy()
            if qudit_indices is not None:
                new_indexes = [qudit_indices[j] for j in gate.qudit_indices]
                new_gate.qudit_indices = new_indexes
            self.add_gate(new_gate)
