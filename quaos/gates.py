
import sys
sys.path.append("./")
import numpy as np
from typing import overload
from qiskit import QuantumCircuit
from quaos.paulis import (
    PauliSum, PauliString, Pauli,
    #  Xnd, Ynd, Znd, Id, symplectic_to_string,
    string_to_symplectic,
)


class GateOperation:
    """
    Mapping can be written as set of rules,
    
    e.g. for CNOT
                x1z0*x0z0 -> x1z0*x1z0
                x0z0*x1z0 -> x0z0*x1z0  # doesn't need specifying
                x0z1*x0z0 -> x0z1*x0z0  # doesn't need specifying
                x0z0*x0z1 -> x0z-1*x0z1 # 

    inputs are:

    mapping = ['x1z0*x0z0 -> x1z0*x1z0', 'x0z0*x0z1 -> -x0z1*x0z1']  # (control*target -> control*target)

    """
    def __init__(self, name: str, qudit_indices: list[int], mapping: list[str], dimension: list[int]):
        self.dimension = dimension
        self.name = name
        self.qudit_indices = qudit_indices
        self.mapping = mapping
        self.map_from, self.map_to, self.acquired_phase = self._interpret_mapping(mapping)
    
    def _interpret_mapping(self, map_string: list[str]) -> tuple[np.ndarray, np.ndarray, list[int]]:
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
        
        symplectic_looked_for = np.array(symplectic_looked_for)
        symplectic_mapped_to = np.array(symplectic_mapped_to)

        return symplectic_looked_for, symplectic_mapped_to, acquired_phase

    def _act_on_pauli_string(self, P: PauliString) -> tuple[PauliString, float]:
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
        return P, acquired_phase

    def _act_on_pauli_sum(self, P: PauliSum) -> PauliSum:
        if np.all(self.acquired_phase == 0):
            # In this case each output of _act_on_pauli_string is a PauliString
            P = PauliSum([self._act_on_pauli_string(p)[0] for p in P.pauli_strings], P.weights, P.phases, P.dimensions, False)
        else:
            # in this case at least one is a PauliSum
            pauli_list = []
            acquired_phase = []
            for p in P.pauli_strings:
                new_pauli_string, additional_phase = self._act_on_pauli_string(p)
                acquired_phase.append(additional_phase)
                pauli_list.append(new_pauli_string)
            weights = P.weights
            P = PauliSum(pauli_list, weights, P.phases, P.dimensions, False)
            P.acquire_phase(acquired_phase)
        # P.combine_equivalent_paulis()
        return P
    
    def copy(self) -> 'GateOperation':
        new_gate = GateOperation(self.name, self.qudit_indices, self.mapping, self.dimension)
        new_gate.acquired_phase = self.acquired_phase
        return new_gate

    @overload
    def act(self, P: Pauli) -> PauliString:
        ...

    @overload
    def act(self, P: PauliString) -> PauliString:
        ...

    @overload
    def act(self, P: PauliSum) -> PauliSum:
        ...

    def act(self, P: Pauli | PauliString | PauliSum) -> PauliString | PauliSum:
        if isinstance(P, Pauli):
            P = P._to_pauli_string()

        if isinstance(P, PauliString):
            return self._act_on_pauli_string(P)[0]
        elif isinstance(P, PauliSum):
            return self._act_on_pauli_sum(P)
        else:
            raise ValueError(f"GateOperation cannot act on type {type(P)}")
    
    def __mul__(self, gate: 'GateOperation') -> 'Circuit':
        # TODO: check if the two gates are compatible, set dimensions accordingly
        # circuit = Circuit([self + gate])  # TODO: add support for gate summation
        if self.dimension != gate.dimension:
            # this is a choice for the moment, that we select the dimensions of the entire circuit from the beginning
            # when defining the gates. We could instead define gates only locally and create the circuit from these 
            # plus the indexes on which they act.
            raise ValueError("Cannot compile Circuit from gates with different dimensions")
        circuit = Circuit(self.dimension, [self, gate])
        return circuit
    
    def __eq__(self, other_gate: 'GateOperation') -> bool:
        if self.name != other_gate.name:
            return False
        if self.qudit_indices != other_gate.qudit_indices:
            return False
        if self.mapping != other_gate.mapping:
            return False
        if self.dimension != other_gate.dimension:
            return False
        if self.acquired_phase != other_gate.acquired_phase:
            return False
        return True


class CNOT(GateOperation):
    def __init__(self, control: int, target: int):
        CNOT_operations = ['x1z0 x0z0 -> x1z0 x1z0', 'x0z0 x0z1 -> x0z1 x0z1']
        super().__init__("CNOT", [control, target], CNOT_operations, dimension=[2, 2])


class Hadamard(GateOperation):
    def __init__(self, index: int, dimension: int, inverse: bool = False):
        if inverse:
            Hadamard_operations = self.inverse_fourier_gate_operations(dimension)
        else:
            Hadamard_operations = self.hadamard_gate_operations(dimension)
        name = "H" if not inverse else "Hdag"
        super().__init__(name, [index], Hadamard_operations, dimension=[dimension])

    @staticmethod
    def hadamard_gate_operations(dimension: int) -> list[str]:
        operations = []
        for r in range(dimension):
            for s in range(dimension):
                phase = (r * s) % dimension
                operations.append(f"x{r}z{s} -> x{-s % dimension}z{r}p{phase}")
        return operations

    @staticmethod
    def inverse_fourier_gate_operations(dimension: int) -> list[str]:
        operations = []
        for r in range(dimension):
            for s in range(dimension):
                phase = (-r * s) % dimension
                operations.append(f"x{r}z{s} -> x{s}z{-r % dimension}p{phase}")
        return operations


class PHASE(GateOperation):
    def __init__(self, index: int, dimension: int):
        SGate_operations = self.s_gate_operations(dimension)  # ['x1z0 -> x1z1', 'x0z1 -> x1z0']
        super().__init__("S", [index], SGate_operations, dimension=[dimension])

    @staticmethod
    def s_gate_operations(dimension: int) -> list[str]:
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
    def sum_gate_operations(dimension: int) -> list[str]:
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


class SWAP(GateOperation):
    def __init__(self, index1, index2, dimension):
        SGate_operations = self.swap_gate_operations(dimension)
        super().__init__("SWAP", [index1, index2], SGate_operations, dimension=dimension)
   
    @staticmethod
    def swap_gate_operations(dimension):
        operations = []
        for r1 in range(dimension):
            for s1 in range(dimension):
                for r2 in range(dimension):
                    for s2 in range(dimension):
                        new_r1 = r2
                        new_s1 = s2
                        new_r2 = r1
                        new_s2 = s1
                        operations.append(f"x{r1}z{s1} x{r2}z{s2} -> x{new_r1}z{new_s1} x{new_r2}z{new_s2}")
        return operations


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
    
    @overload
    def act(self, pauli: Pauli) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliString) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliSum) -> PauliSum:
        ...

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
                new_indices = [qudit_indices[j] for j in gate.qudit_indices]
                new_gate.qudit_indices = new_indices
            self.add_gate(new_gate)


if __name__ == "__main__":
    import sys
    from pathlib import Path
    root_path = Path(__file__).parent.parent
    print(root_path)
    sys.path.append(str(root_path))
    from hamiltonian import random_pauli_hamiltonian
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
    CNOT3 = SUM(0, 1, 2)

    ps1 = PauliString('x1z0 x0z0', dimensions=[2, 2])
    ps2 = PauliString('x0z0 x0z1', dimensions=[2, 2])
    ps3 = PauliString('x1z0 x1z0', dimensions=[2, 2])
    ps4 = PauliString('x0z1 x0z1', dimensions=[2, 2])
    ps5 = PauliString('x1z1 x0z0', dimensions=[2, 2])

    print(ps1, '->', CNOT3.act(ps1))
    print(ps2, '->', CNOT3.act(ps2))
    print(ps3, '->', CNOT3.act(ps3))
    print(ps4, '->', CNOT3.act(ps4))

    # psum1 = ps1 + 0.5 * ps2 + 1j * ps3 - 0.5j * ps4
    # print(psum1, '\n -> \n', CNOT3.act(psum1))

    Hg = Hadamard(0, 2)
    # print(ps1, '->', H.act(ps1))
    # print(ps2, '->', H.act(ps2))
    # print(ps3, '->', H.act(ps3))
    # print(ps4, '->', Hg.act(ps4))
    # print(ps5, '->', Hg.act(ps5))

    # ps = random_pauli_hamiltonian(5, [2] * 5, mode='uniform')
    # Sum01 = SUM(0, 1, 2)

    # print(Sum01.act(ps))
    # print(Sum03.act(Sum02.act(Sum01.act(ps))))RP
    # c = Circuit([2 * 5])
    # ps2, c = cancel_X(ps, 0, 5, c, 5)

    # print(ps2)

