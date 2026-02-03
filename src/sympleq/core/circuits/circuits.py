from typing import Generator, overload, TypeVar, Sequence
import numpy as np
from .gates import Gate, Hadamard as H, PHASE as S, SUM as CX, SWAP, CNOT, PauliGate
from sympleq.core.paulis import PauliSum, PauliString, Pauli, PauliObject
from sympleq.core.states.state import State
from .utils import embed_symplectic
import scipy.sparse as sp
from collections import defaultdict
import random


# We define a type using TypeVar to let the type checker know that
# the input and output of the `act` function share the same type.
P = TypeVar("P", bound="PauliObject")


class Circuit:
    def __init__(self, dimensions: list[int] | np.ndarray,
                 gates: Sequence[Gate] | None = None):
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
        self.dimensions = dimensions
        self.gates = list(gates) if gates is not None else []
        if gates is not None:
            self.indexes = [gate.qudit_indices for gate in gates]
            # indexes accessible at the Circuit level - note really necessary,
            #  we can see how much this is used in practice
        else:
            self.indexes = []

    @classmethod
    def from_random(
        cls,
        depth: int,
        dimensions: list[int] | np.ndarray,
        single_qudit_fill: float = 0.5,
        two_qudit_fill: float = 0.3,
    ) -> 'Circuit':
        """
        Creates a random circuit with the given depth.

        Each layer independently chooses (probabilistically) to place a 1-qudit gate,
        a 2-qudit gate, or nothing:
          - probability of single-qudit gate = single_qudit_fill
          - probability of two-qudit gate    = two_qudit_fill
          - probability of no gate           = 1 - (single_qudit_fill + two_qudit_fill)
        If a two-qudit gate is chosen but the selected dimension set only has one qudit,
        the layer falls back to a single-qudit gate when possible, otherwise it is skipped.

        Args:
            depth: Number of layers.
            dimensions: List/array of qudit dimensions.
            single_qudit_fill: Expected fraction of layers with a 1-qudit gate.
            two_qudit_fill: Expected fraction of layers with a 2-qudit gate.
        """
        def index_lists(lst):
            groups = defaultdict(list)
            for i, val in enumerate(lst):
                groups[val].append(i)
            return list(groups.values())
        index_sets = index_lists(dimensions)  # list of lists of indexes for each dimension
        n_dims = len(index_sets)  # number of different dimensions

        single_qudit_gates = [H, S]
        two_qudit_gates = [CX, SWAP]
        gg = []
        for _ in range(depth):
            set_idx = np.random.randint(n_dims)
            dim = dimensions[index_sets[set_idx][0]]
            r = np.random.rand()
            gate_added = False
            if r < two_qudit_fill:
                if len(index_sets[set_idx]) > 1:
                    i0, i1 = random.sample(index_sets[set_idx], 2)
                    gate_cls = random.choice(two_qudit_gates)
                    gg.append(gate_cls(i0, i1, dim))
                    gate_added = True
            elif r < two_qudit_fill + single_qudit_fill:
                index = random.choice(index_sets[set_idx])
                gate_cls = random.choice(single_qudit_gates)
                gg.append(gate_cls(index, dim))
                gate_added = True

            # If a two-qudit gate was requested but impossible, fall back to single if allowed.
            if not gate_added and len(index_sets[set_idx]) == 1 and single_qudit_fill > 0:
                index = random.choice(index_sets[set_idx])
                gate_cls = random.choice(single_qudit_gates)
                gg.append(gate_cls(index, dim))

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

    def n_qudits(self) -> int:
        """
        Returns the number of qudits in the circuit.
        """
        return len(self.dimensions)

    def __add__(self, other: "Circuit | Gate") -> "Circuit":
        """
        Adds two circuits together by concatenating their gates and indexes.
        """
        if not isinstance(other, Circuit) and not isinstance(other, Gate):
            raise TypeError("Can only add another Circuit or Gate object.")
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

    # def _composite_phase_vector(self, F_1: np.ndarray, F_2: np.ndarray, h_2: np.ndarray, lcm: int) -> np.ndarray:
    #     """
    #     Returns the vector to add to h_1 to obtain h'' in PHYSICAL REVIEW A 71, 042315 (2005) - Eq. (8)

    #     New phase vector is h_1 + h_c

    #     """
    #     U = np.zeros((2 * self.n_qudits(), 2 * self.n_qudits()), dtype=int)
    #     U[self.n_qudits():, :self.n_qudits()] = np.eye(self.n_qudits(), dtype=int)

    #     U_conjugated = F_2.T @ U @ F_2

    #     p1 = np.dot(F_1, h_2)
    #     # negative sign in below as definition in paper is strictly upper diagonal, not including diagonal part
    #     p2 = np.diag(np.dot(F_1, np.dot((2 * np.triu(U_conjugated) - np.diag(np.diag(U_conjugated))), F_1.T)))
    #     p3 = np.dot(F_1, np.diag(U_conjugated))

    #     h_c = (p1 + p2 - p3) % (2 * lcm)

        # return h_c

    def _composite_phase_vector(self, F_1: np.ndarray, F_2: np.ndarray, h_2: np.ndarray, lcm: int) -> np.ndarray:
        """
        Faster equivalent of:
            U = [[0,0],[I,0]]
            Uc = F2.T @ U @ F2
            p2 = diag(F1 @ ( (2*triu(Uc)-diag(diag(Uc))) @ F1.T ))
        Uses:
            p2 = sum( (F1 @ Q) * F1, axis=1 )   (one matmul)
        and block construction of Uc.
        """
        F1 = np.ascontiguousarray(F_1, dtype=np.int64)
        F2 = np.ascontiguousarray(F_2, dtype=np.int64)
        h2 = np.asarray(h_2, dtype=np.int64)

        n2 = F1.shape[0]
        assert n2 % 2 == 0 and F1.shape == (n2, n2) and F2.shape == (n2, n2)
        n = n2 // 2
        mod = 2 * int(lcm)

        # ---- Build U_conjugated = F2^T U F2 using blocks (A,B,C,D) ----
        A = F2[:n, :n]
        B = F2[:n, n:]
        C = F2[n:, :n]
        D = F2[n:, n:]

        # Uc = [[C^T A, C^T B],
        #       [D^T A, D^T B]]
        Uc = np.empty((n2, n2), dtype=np.int64)
        Uc[:n, :n] = C.T @ A
        Uc[:n, n:] = C.T @ B
        Uc[n:, :n] = D.T @ A
        Uc[n:, n:] = D.T @ B

        # diag(Uc) without forming np.diag(Uc)
        diag_uc = np.empty(n2, dtype=np.int64)
        diag_uc[:n] = np.sum(C * A, axis=0)
        diag_uc[n:] = np.sum(D * B, axis=0)

        # ---- p1 and p3 ----
        p1 = F1 @ h2
        p3 = F1 @ diag_uc

        # ---- Build Q = 2*triu(Uc) - diag(diag(Uc)) ----
        # After this, Uc == Q
        tril_i, tril_j = np.tril_indices(n2, -1)
        Uc[tril_i, tril_j] = 0          # keep only upper triangle + diag
        Uc *= 2                         # doubles diag too (we will restore)
        np.fill_diagonal(Uc, diag_uc)   # restore diagonal to original (not doubled)

        # ---- p2 with ONE matmul ----
        F1Q = F1 @ Uc                   # only one matmul
        p2 = np.sum(F1Q * F1, axis=1)   # diag(F1 Q F1^T)

        return (p1 + p2 - p3) % mod

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
        C_inv = Circuit(self.dimensions, [g.inv() for g in self.gates])
        return C_inv

    def full_symplectic(self):
        return self.composite_gate().full_symplectic(self.n_qudits())

    def cleanup(self):
        # TODO If two gates are the inverse of each other and next to each other, remove them both. This happens
        # in a few algorithms
        raise NotImplementedError

    def local_circuit(self, indices: list[int] | np.ndarray) -> 'Circuit':
        """
        Returns a Circuit object containing only the gates that act only on the specified qudit indices.
        """
        indices = list(indices)
        index_map = {old: new for new, old in enumerate(indices)}

        def _remap_gate(gate: Gate) -> Gate | None:
            # Skip gates that touch qudits outside the requested subset
            if not all(idx in index_map for idx in gate.qudit_indices):
                return None

            # Handle PauliGate separately to keep its pauli_string consistent with the new local qudit ordering
            if isinstance(gate, PauliGate):
                ps = gate.pauli_string
                x_local = []
                z_local = []
                dims_local = []
                for old_idx in indices:
                    x_local.append(int(ps.x_exp[old_idx]))
                    z_local.append(int(ps.z_exp[old_idx]))
                    dims_local.append(int(ps.dimensions[old_idx]))
                ps_local = PauliString.from_exponents(x_local, z_local, dims_local)
                try:
                    return PauliGate(ps_local, name=gate.name)
                except ValueError:
                    # PauliGate requires at least one non-trivial component; skip if trivial on this subset.
                    return None

            # Generic gate: copy and remap indices/dimensions to the local numbering
            g_copy = gate.copy()
            g_copy.qudit_indices = np.asarray([index_map[idx] for idx in gate.qudit_indices], dtype=int)
            g_copy.dimensions = np.asarray([self.dimensions[idx] for idx in gate.qudit_indices], dtype=int)
            return g_copy

        local_gates = []
        for gate in self.gates:
            remapped = _remap_gate(gate)
            if remapped is not None:
                local_gates.append(remapped)

        local_dimensions = [self.dimensions[i] for i in indices]
        return Circuit(local_dimensions, local_gates)

    def act_on_state(self, state: State) -> State:
        """
        Apply the entire circuit to a State object, gate by gate.
        """
        if not np.array_equal(state.dimensions, np.asarray(self.dimensions, dtype=int)):
            raise ValueError(
                "State dimensions do not match Circuit.dimensions: "
                f"{state.dimensions} vs {self.dimensions}"
            )

        out = state
        for gate in self.gates:
            out = gate.act_on_state(out)
        return out

    def apply_to_statevector(self, psi: np.ndarray) -> np.ndarray:
        """
        Convenience wrapper: apply circuit directly to a bare statevector psi.

        psi must be a 1D array of length prod(self.dimensions).
        Returns the updated statevector.
        """
        state = State(psi, self.dimensions)
        state_out = self.act_on_state(state)
        return state_out.as_array()
