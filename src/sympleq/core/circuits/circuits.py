from __future__ import annotations
from typing import Generator, overload, TypeVar
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from .gates import Gate, Hadamard as H, PHASE as S, SUM as CX, SWAP, CNOT, PauliGate
from .utils import embed_symplectic
import scipy.sparse as sp
from collections import defaultdict

from sympleq.core.paulis import PauliSum, PauliString, Pauli, PauliObject


# We define a type using TypeVar to let the type checker know that
# the input and output of the `act` function share the same type.
P = TypeVar("P", bound="PauliObject")


class Circuit:
    def __init__(self, dimensions: list[int] | np.ndarray,
                 gates: list[Gate] = []):
        """
        Initialize the Circuit with gates, indexes, and targets.

        If a multi-qubit gate has a target, the targets should be at the end of the tuple of indexes
        e.g. a CNOT with control 1, target 3 is

        gate = 'CNOT'
        indexes = (1, 3)


        Parameters:
            dimensions (list[int] | np.ndarray): A list or array of integers representing the dimensions of the qudits.
            gates (list): A list of Gate objects representing the gates in the circuit.


        TODO: Remove dimensions as input - this can be obtained from the gates only - make this a method not attribute

        TODO: Perhaps store the composite gate as an attribute - it will allow gate.act to be significantly faster
        """
        if gates is None:
            gates = []

        dimensions = np.asarray(dimensions, dtype=int)
        self.dimensions = dimensions
        self.gates = gates
        self.indexes = [gate.qudit_indices for gate in gates]  # indexes accessible at the Circuit level

    @classmethod
    def from_random(cls, n_gates: int,
                    dimensions: list[int] | np.ndarray,
                    two_qudit_gate_ratio: float = 0.3,
                    rng: RNGGenerator | None = None) -> Circuit:
        """
        Creates a random circuit with the given number of qudits and depth.

        NOTE: It may be nice to have depth rather than n_gates, and a filling factor to control the number of gates
              per layer? Not too important.

        Parameters:
            n_qudits (int): The number of qudits in the circuit.
            depth (int): The depth of the circuit.
            dimensions (list[int] | np.ndarray): A list or array of integers representing the dimensions of the qudits.
            gate_list (list): A list of Gate objects representing the gates in the circuit.
            two_qudit_gate_ratio (float): The ratio of two-qudit gates to single-qudit gates.

        Returns:
            Circuit: A new Circuit object.
        """

        def index_lists(lst):
            groups = defaultdict(list)
            for i, val in enumerate(lst):
                groups[val].append(i)
            return list(groups.values())

        if rng is None:
            rng = default_rng()

        index_sets = index_lists(dimensions)  # list of lists of indexes for each dimension
        n_dims = len(index_sets)  # number of different dimensions

        single_qudit_gates = [H, S]
        two_qudit_gates = [CX, SWAP]
        gg = []
        for _ in range(n_gates):
            set_idx = rng.integers(0, n_dims)
            dim = dimensions[index_sets[set_idx][0]]
            if rng.random() < two_qudit_gate_ratio and len(index_sets[set_idx]) > 1:
                indices = rng.choice(index_sets[set_idx], 2, replace=False)
                gate_cls = rng.choice(two_qudit_gates)
                assert indices[0] != indices[1]
                gg.append(gate_cls(indices[0], indices[1], dim))
            else:
                index = rng.choice(index_sets[set_idx])
                gate = rng.choice(single_qudit_gates)
                gg.append(gate(index, dim))

        return cls(dimensions, gg)

    def add_gate(self, gate: Gate | list[Gate]):
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

    def n_gates(self) -> int:
        """
        Returns the number of gates in the circuit.
        """
        return len(self.gates)

    def n_qudits(self) -> int:
        """
        Returns the number of qudits in the circuit.
        """
        return len(self.dimensions)

    def __add__(self, other: Circuit | Gate) -> Circuit:
        """
        Adds two circuits together by concatenating their gates and indexes.
        """
        if not isinstance(other, Circuit) and not isinstance(other, Gate):
            raise TypeError("Can only add another Circuit or Gate object.")

        if not np.array_equal(self.dimensions, other.dimensions):
            raise ValueError("Cannot add circuits with different dimensions.")

        if isinstance(other, Gate):
            new_gates = self.gates + [other]
        else:
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

    def __getitem__(self, index: int) -> Gate:
        return self.gates[index]

    def __setitem__(self, index: int, value: Gate):
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
    def act(self, pauli: Pauli) -> Pauli:
        ...

    @overload
    def act(self, pauli: PauliString) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliSum) -> PauliSum:
        ...

    def act(self, pauli: Pauli | PauliString | PauliSum) -> Pauli | PauliString | PauliSum:
        for gate in self.gates:
            pauli = gate.act(pauli)

        return pauli

    @overload
    def act_iter(self, pauli: Pauli) -> Generator[Pauli, None, None]:
        ...

    @overload
    def act_iter(self, pauli: PauliString) -> Generator[PauliString, None, None]:
        ...

    @overload
    def act_iter(self, pauli: PauliSum) -> Generator[PauliSum, None, None]:
        ...

    def act_iter(self, pauli: Pauli | PauliString | PauliSum) -> Generator[Pauli | PauliString | PauliSum, None, None]:
        for gate in self.gates:
            pauli_sum = gate.act(pauli)
            yield pauli_sum

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
                new_gate.qudit_indices = np.ndarray(new_indexes)
            self.add_gate(new_gate)

    def _composite_phase_vector(self, F_1: np.ndarray, F_2: np.ndarray, h_2: np.ndarray, lcm: int) -> np.ndarray:
        """
        Returns the vector to add to h_1 to obtain h'' in PHYSICAL REVIEW A 71, 042315 (2005) - Eq. (8)

        New phase vector is h_1 + h_c

        """
        U = np.zeros((2 * self.n_qudits(), 2 * self.n_qudits()), dtype=int)
        U[self.n_qudits():, :self.n_qudits()] = np.eye(self.n_qudits(), dtype=int)

        U_conjugated = F_2.T @ U @ F_2

        p1 = np.dot(F_1, h_2)
        # negative sign in below as definition in paper is strictly upper diagonal, not including diagonal part
        p2 = np.diag(np.dot(F_1, np.dot((2 * np.triu(U_conjugated) - np.diag(np.diag(U_conjugated))), F_1.T)))
        p3 = np.dot(F_1, np.diag(U_conjugated))

        h_c = (p1 + p2 - p3) % (2 * lcm)

        return h_c

    def composite_gate(self) -> Gate:
        """Composes the list of symplectics acting on all qudits to a single symplectic"""

        n_qudits = self.n_qudits()
        total_symplectic = np.eye(2 * n_qudits, dtype=np.uint8)
        lcm = np.lcm.reduce(self.dimensions)
        total_phase_vector = np.zeros(2 * n_qudits, dtype=int)

        for i, gate in enumerate(self.gates):
            symplectic = gate.symplectic
            indexes = gate.qudit_indices
            phase_vector = gate.phase_vector

            F, h = embed_symplectic(symplectic, phase_vector, indexes, self.n_qudits())  #
            if i == 0:
                total_phase_vector = h
            else:
                total_phase_vector = np.mod(total_phase_vector + self._composite_phase_vector(total_symplectic, F, h,
                                                                                              lcm),
                                            2 * lcm)

            total_symplectic = np.mod(total_symplectic @ F.T, lcm)

        total_indexes = list(range(n_qudits))
        total_symplectic = total_symplectic.T
        return Gate('CompositeGate', total_indexes, total_symplectic, self.dimensions, total_phase_vector)

    def unitary(self):
        known_unitaries = (H, S, CX, SWAP, CNOT, PauliGate)
        if not np.all([isinstance(gate, known_unitaries) for gate in self.gates]):
            print([(gate.name, isinstance(gate, known_unitaries)) for gate in self.gates])
            raise NotImplementedError("Unitary not implemented for all gates in the circuit.")

        q = self.dimensions
        m = sp.csr_matrix(([1] * (np.prod(q)), (range(np.prod(q)), range(np.prod(q)))))
        for g in self.gates:
            m = g.unitary(dims=self.dimensions) @ m

        return m

    def inv(self):
        C_inv = Circuit(self.dimensions, [g.inv() for g in self.gates[::-1]])
        return C_inv

    def full_symplectic(self):
        return self.composite_gate().full_symplectic(self.n_qudits())

    def cleanup(self):
        # TODO If two gates are the inverse of each other and next to each other, remove them both. This happens
        # in a few algorithms
        raise NotImplementedError

    def gates_layout(self,
                     with_qudit_indices: bool = False,
                     with_input: PauliSum | None = None,
                     with_output: PauliSum | None = None,
                     wires: str | list[str] | None = None,
                     wrap: bool = True) -> str:
        """
        Returns a visual circuit diagram of the RMB.

        Renders the circuit as ASCII art with gates displayed as boxes
        connected by wires.

        Parameters
        ----------
        with_qudit_indices : bool, default False
            If True, display qudit indices on the left of each wire.
        with_input : PauliSum | None, default None
            If provided, display input phases on the left.
        with_output : PauliSum | None, default None
            If provided, display output phases on the right.
        wires : str | list[str] | None, default None
            If provided, overrides default wires strin. If a list is provided, its length must match circuit n_qudits.
        wrap : bool, default True
            If True, wrap the output to fit the terminal width by splitting
            at gate boundaries.

        Returns
        -------
        str
            A string representation of the circuit diagram.
        """

        def gate_name(gate: Gate) -> str:
            return gate.name.replace("-inv", "*")[:gate_name_len].center(gate_name_len)

        n_qudits = self.n_qudits()
        lines: list[str] = ["" for _ in range(3 * n_qudits)]
        gate_num: list[int] = [0 for _ in range(n_qudits)]

        if wires is None:
            wires = ["="] * n_qudits
        elif isinstance(wires, str):
            wires = wires * n_qudits
        elif len(wires) != n_qudits:
            raise ValueError("Wires list length must match the circuit number of qudits.")

        gate_name_len = 5
        gate_len = gate_name_len + 4

        if with_qudit_indices:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] = " " * 4
                lines[3 * l_idx + 1] = f"{l_idx:>2}: "
                lines[3 * l_idx + 2] = " " * 4
        else:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] = ""
                lines[3 * l_idx + 1] = ""
                lines[3 * l_idx + 2] = ""

        if with_input is None:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += " " * 2 + wires[l_idx] * 2
                lines[3 * l_idx + 2] += " " * 4
        else:
            # Put initial state phase on the left
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += f"{with_input.phases[l_idx]} " + wires[l_idx] * 2
                lines[3 * l_idx + 2] += " " * 4

        for gate in self.gates:
            if gate.n_qudits == 1:
                l_idx = gate.qudit_indices[0]
                gate_num[l_idx] += 1

                lines[3 * l_idx + 0] += " ┌" + "─" * gate_name_len + "┐ "
                lines[3 * l_idx + 1] += wires[l_idx] + "│" + f"{gate_name(gate)}" + "│" + wires[l_idx]
                lines[3 * l_idx + 2] += " └" + "─" * gate_name_len + "┘ "
            # 2-qudit gate
            else:
                # Get max line length of affected qudits
                max_num_gate_affected_qudits = max(
                    [gate_num[idx] for idx in gate.qudit_indices])

                for l_idx in gate.qudit_indices:
                    while gate_num[l_idx] < max_num_gate_affected_qudits:
                        gate_num[l_idx] += 1
                        lines[3 * l_idx + 0] += " " * gate_len
                        lines[3 * l_idx + 1] += wires[l_idx] * gate_len
                        lines[3 * l_idx + 2] += " " * gate_len

                    gate_num[l_idx] += 1

                    is_top_qudit = l_idx == min(gate.qudit_indices)
                    is_btm_qudit = l_idx == max(gate.qudit_indices)

                    if is_top_qudit:
                        lines[3 * l_idx + 0] += " ┌" + "─" * gate_name_len + "┐ "
                    else:
                        lines[3 * l_idx + 0] += " ┌" + "─" * (gate_name_len // 2) + \
                            "┴" + "─" * (gate_name_len // 2) + "┐ "

                    lines[3 * l_idx + 1] += wires[l_idx] + "│" + f"{gate_name(gate)}" + "│" + wires[l_idx]
                    if is_btm_qudit:
                        lines[3 * l_idx + 2] += " └" + "─" * gate_name_len + "┘ "
                    else:
                        lines[3 * l_idx + 2] += " └" + "─" * (gate_name_len // 2) + \
                            "┬" + "─" * (gate_name_len // 2) + "┘ "

        max_num_gate = max(gate_num)
        for l_idx in range(n_qudits):
            while gate_num[l_idx] < max_num_gate:
                gate_num[l_idx] += 1
                lines[3 * l_idx + 0] += " " * gate_len
                lines[3 * l_idx + 1] += wires[l_idx] * gate_len
                lines[3 * l_idx + 2] += " " * gate_len

        if with_output is None:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += wires[l_idx] * 2 + " " * 2
                lines[3 * l_idx + 2] += " " * 4
        else:
            # Put final state phase on the right
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += wires[l_idx] * 2 + " " + f"{with_output.phases[l_idx]}"
                lines[3 * l_idx + 2] += " " * 4

        if not wrap:
            return "\n".join(lines)

        # Wrap output to fit terminal width
        import shutil
        import re
        term_width = shutil.get_terminal_size().columns
        line_len = len(lines[0])  # Top border has no ANSI codes

        if line_len <= term_width:
            return "\n".join(lines)

        # Strip ANSI codes from wire lines for correct slicing
        ansi_pattern = re.compile(r'\033\[[0-9;]*m')
        plain_wire_lines = {
            l_idx: ansi_pattern.sub('', lines[3 * l_idx + 1])
            for l_idx in range(n_qudits)
        }

        sections = []
        pos = 0

        while pos < line_len:
            if pos == 0:
                # First section includes prefix from input state
                prefix_len = 8 if with_qudit_indices else 4
                gates_per_section = max(1, (term_width - prefix_len - 2) // gate_len)
                end = min(prefix_len + gates_per_section * gate_len, line_len)
                section_lines = []
                for l_idx in range(n_qudits):
                    section_lines.append(lines[3 * l_idx + 0][pos:end])
                    wire_slice = plain_wire_lines[l_idx][pos:end]
                    wire_slice = wire_slice.replace("=", wires[l_idx])
                    section_lines.append(wire_slice)
                    section_lines.append(lines[3 * l_idx + 2][pos:end])
            else:
                # Continuation sections: add wire prefix for visual continuity
                extra_prefix_side = " " * 4 if with_qudit_indices else ""
                prefix_len = len(extra_prefix_side)
                gates_per_section = max(1, (term_width - prefix_len - 2) // gate_len)
                end = min(pos + gates_per_section * gate_len, line_len)
                section_lines = []

                for l_idx in range(n_qudits):
                    extra_prefix_central = f"{l_idx:>2}: " if with_qudit_indices else ""
                    section_lines.append(extra_prefix_side + "  " + lines[3 * l_idx + 0][pos:end])
                    wire_slice = plain_wire_lines[l_idx][pos:end]
                    wire_slice = wire_slice.replace("=", wires[l_idx])
                    section_lines.append(extra_prefix_central + wires[l_idx] * 2 + wire_slice)
                    section_lines.append(extra_prefix_side + "  " + lines[3 * l_idx + 2][pos:end])

            sections.append("\n".join(section_lines))
            pos = end

        return "\n\n\n".join(sections)
