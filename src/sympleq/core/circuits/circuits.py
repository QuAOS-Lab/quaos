from __future__ import annotations
from typing import Generator, overload, TypeVar, TypeAlias
import numpy as np
from qiskit import QuantumCircuit

from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from .utils import embed_symplectic
import scipy.sparse as sp
import random

from .gates import GATES, Hadamard as H, SUM, PHASE, Gate, SWAP
from sympleq.core.paulis import PauliSum, PauliString, Pauli, PauliObject

GateTuple: TypeAlias = tuple[Gate, *tuple[int, ...]]


# We define a type using TypeVar to let the type checker know that
# the input and output of the `act` function share the same type.
P = TypeVar("P", bound="PauliObject")


class Circuit:
    def __init__(self, n_qudits: int | None = None,
                 gates: np.ndarray | list[Gate] = [], qudits: np.ndarray | list[list[int]] = []):
        """
        Initialize the Circuit with gates, qudits.

        If a multi-qubit gate has a target, the targets should be at the ent of the tuple of indexes
        e.g. a CNOT with control 1, target 3 is

        gate = 'CNOT'
        qudits = (1, 3)


        Parameters:
            dimensions (list[int] | np.ndarray): A list or array of integers representing the dimensions of the qudits.
            gates (list): A list of Gate objects representing the gates in the circuit.


        TODO: Remove dimensions as input - this can be obtained from the gates only - make this a method not attribute

        TODO: Perhaps store the composite gate as an attribute - it will allow gate.act to be significantly faster
        """
        if n_qudits is None:
            n_qudits = int(max([idx for q_indices in qudits for idx in q_indices]) + 1)
        self._n_qudits = n_qudits
        self._gates = np.asarray(gates, dtype=object)
        self._qudits = np.asarray(qudits, dtype=object)

    @classmethod
    def empty(cls) -> Circuit:
        return Circuit(n_qudits=0, gates=[], qudits=[])

    @classmethod
    def from_random(cls, n_qudits: int, depth: int) -> Circuit:
        """
        Creates a random circuit with the given number of qudits and depth.

        Parameters:
            n_qudits (int): The number of qudits in the circuit.
            depth (int): The depth of the circuit.

        Returns:
            Circuit: A new Circuit object.
        """
        # FIXME: add weight of 2 qubits gates
        if n_qudits > 1:
            available_gates = np.asarray(
                [GATES.H, GATES.S, GATES.swap, GATES.cnot, GATES.sum])
        else:
            available_gates = np.asarray(
                [GATES.H, GATES.S])

        gates = np.random.choice(available_gates, size=depth, replace=True)
        qudits = [[random.randint(0, n_qudits - 1) for _ in range(g.n_qudits())] for g in gates]
        C = cls(n_qudits, gates, qudits)
        C._sanity_check()
        return C

    @classmethod
    def from_data(cls, data: GateTuple | list[GateTuple]) -> Circuit:
        """
        Creates a circuit from the given gates and qudits data.

        Parameters:
            data: list[tuple[Gate, *tuple[int, ...]]]
                The circuit gates and qudits, given as a list of tuples.
                Each tuple contains the Gate at first position and one or more qudit index as integer.

        Returns:
            Circuit: A new Circuit object.
        """

        if isinstance(data, tuple):
            data = [data]

        gates = [d[0] for d in data]
        qudits = [list(d[1:]) for d in data]
        C = cls(gates=gates, qudits=qudits)
        C._sanity_check()
        return C

    def n_qudits(self) -> int:
        return self._n_qudits

    size = n_qudits

    def gates(self) -> np.ndarray:
        return self._gates

    def qudits(self) -> np.ndarray:
        return self._qudits

    def add_gate(self, gate: Gate, *qudits: int):
        """
        Appends a gate to qudit index with specified target (if relevant).
        If the number of qudits of the Circuit is smaller than any qudit index,
        it is increased as to match it.
        """

        if len(qudits) != gate.n_qudits():
            raise ValueError(f"Gate {gate} acts on {gate.n_qudits()} qudits, but {len(qudits)} qudits were passed.")

        self._gates = np.concatenate(self._gates, np.asarray(gate), dtype=object)
        self._qudits = np.concatenate(self._qudits, np.asarray(list[qudits]), dtype=object)
        self._n_qudits = int(max([idx for q_indices in self._qudits for idx in q_indices]) + 1)

    def add_gates(self, gates: list[Gate], *qudits: list[int]):
        """
        Appends a gate to qudit index with specified target (if relevant).
        If the number of qudits of the Circuit is smaller than any qudit index,
        it is increased as to match it.
        """

        for (gate_qudits, gate) in zip(qudits, gates):
            if len(gate_qudits) != gate.n_qudits():
                raise ValueError(f"Gate {gate} acts on {gate.n_qudits()} qudits, but {len(qudits)} qudits were passed.")

        self._gates = np.concatenate(self._gates, np.asarray(gates), dtype=object)
        self._qudits = np.concatenate(self._qudits, np.asarray(qudits), dtype=object)
        self._n_qudits = int(max([idx for q_indices in self._qudits for idx in q_indices]) + 1)

    def remove_gate(self, index: int):
        """
        Removes a gate from the circuit at the specified index
        """
        np.delete(self._gates, index)
        np.delete(self._qudits, index)

    def __add__(self, other: Circuit) -> Circuit:
        """
        Adds two circuits together by concatenating their gates and indexes.
        """
        if not isinstance(other, Circuit):
            raise TypeError("Can only add another Circuit object.")

        new_gates = np.concatenate(self.gates(), other.gates())
        new_qudits = np.concatenate(self.qudits(), other.qudits())

        return Circuit(gates=new_gates, qudits=new_qudits)

    def __eq__(self, other: Circuit) -> bool:
        if not isinstance(other, Circuit):
            return False
        if len(self.gates()) != len(other.gates()):
            return False
        for i in range(len(self.gates())):
            if self.gates()[i] != other.gates()[i]:
                return False
        if len(self.qudits()) != len(other.qudits()):
            return False
        for i in range(len(self.qudits())):
            if self.qudits()[i] != other.qudits()[i]:
                return False
        return True

    def __getitem__(self, index: int) -> tuple[Gate, tuple[int, ...]]:
        return self.gates()[index], self.qudits()[index]

    def __setitem__(self, index: int, value: Gate):
        self.gates()[index] = value

    def __len__(self) -> int:
        return len(self.gates())

    def __str__(self) -> str:
        return "\n".join([f"{gate.name} {' '.join(qudits)}" for qudits, gate in zip(self.qudits(), self.gates())])

    @overload
    def act(self, pauli: Pauli) -> Pauli:
        ...

    @overload
    def act(self, pauli: PauliString) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliSum) -> PauliSum:
        ...

    @overload
    def act(self, pauli: Pauli | PauliString | PauliSum) -> Pauli | PauliString | PauliSum:
        ...

    def act(self, pauli: P) -> P:
        for (qudits, gate) in zip(self.qudits(), self.gates()):
            pauli = gate.act(pauli, qudits)

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

    def act_iter(self, pauli: P) -> Generator[P, None, None]:
        for (qudits, gate) in zip(self.qudits(), self.gates()):
            pauli_sum = gate.act(pauli, qudits)
            yield pauli_sum

    def show(self, n_qudits: int):
        circuit = QuantumCircuit(n_qudits)
        dict = {'X': circuit.x, 'H': circuit.h, 'S': circuit.s, 'SUM': circuit.cx, 'CNOT': circuit.cx,
                'Hdag': circuit.h}

        for (qudits, gate) in zip(self.qudits(), self.gates()):
            name = gate.name
            if gate.n_qudits() == 2:
                dict[name](qudits[0], qudits[1])
            else:
                dict[name](qudits[0])

        print(circuit)
        # return circuit

    def copy(self) -> Circuit:
        return Circuit(self.n_qudits(), self.gates().copy(), self.qudits().copy())

    def embed_circuit(self, circuit: Circuit, qudits: list[int] | np.ndarray | None = None):
        """
        Embed a circuit into current circuit at the specified qudit indices.
        """
        # FIXME: ask what this is supposed to do

        # for (qudits, gate) in zip(self.qudits(), self.gates()):
        #     new_qudits = [qudits[j] for j in gate.qudits]
        #     new_gate.qudits = np.ndarray(new_indexes)
        #     self.add_gate(new_gate)

    def _composite_phase_vector(self, F_1: np.ndarray, F_2: np.ndarray, h_2: np.ndarray) -> np.ndarray:
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

        # NOTE: we do not take modulo 2*lcm here.
        h_c = (p1 + p2 - p3)

        return h_c

    def composite_gate(self) -> Gate:
        """Composes the list of symplectics acting on all qudits to a single symplectic"""

        total_indexes = []
        total_symplectic = np.eye(2 * self.n_qudits(), dtype=np.uint8)
        for i, (qudits, gate) in enumerate(zip(self.qudits(), self.gates())):
            symplectic = gate.symplectic
            phase_vector = gate.phase_vector

            F, h = embed_symplectic(symplectic, phase_vector, qudits, self.n_qudits())
            if i == 0:
                total_phase_vector = h
            else:
                total_phase_vector = total_phase_vector + self._composite_phase_vector(total_symplectic, F, h)

            # NOTE: we do not take modulo 2*lcm here.
            total_symplectic = total_symplectic @ F.T
            total_indexes.extend(qudits)

        total_indexes = list(set(np.sort(total_indexes)))
        total_symplectic = total_symplectic.T
        return Gate('CompositeGate', total_symplectic, total_phase_vector)

    def unitary(self, dimensions: int | list[int] | np.ndarray | None = None) -> sp.csr_matrix:
        if dimensions is None:
            dimensions = np.ones(self.n_qudits()) * DEFAULT_QUDIT_DIMENSION
        else:  # Catches int but also list and arrays of length 1
            dimensions = np.asarray(dimensions, dtype=int)
            if dimensions.ndim == 0:
                dimensions = np.full(self.n_qudits(), dimensions.item(), dtype=int)

        known_unitaries = (H, PHASE, SUM, SWAP)
        if not np.all([isinstance(gate, known_unitaries) for gate in self.gates()]):
            raise NotImplementedError("Unitary not implemented for all gates in the circuit.")

        m = sp.csr_matrix(([1] * (np.prod(dimensions)), (range(np.prod(dimensions)), range(np.prod(dimensions)))))
        for g in self.gates():
            m = g.unitary(dimensions) @ m

        return m

    def _sanity_check(self):
        """
        Validates the consistency of the Circuit internal representation.

        Raises
        ------
        ValueError
            If any qudit index is large equal than the number of qudits.
        """

        if len(self.gates()) != len(self.qudits()):
            raise ValueError(
                f"There should be the same number of gates ({len(self.gates())}) and qudits ({len(self.qudits())}).")

        q = self.qudits()

        # Check that each element in qudit is a non-empty list
        lengths = np.fromiter((len(lst) for lst in q), dtype=int)
        empty_mask = lengths == 0
        if np.any(empty_mask):
            empty_indices = np.where(empty_mask)[0]
            raise ValueError(f"Gates and indices {empty_indices} must be applied to at least one qudit.")

        all_qudits = np.fromiter((x for lst in q for x in lst), dtype=int, count=sum(lengths))
        if not np.all((all_qudits >= 0) & (all_qudits < self.n_qudits())):
            raise ValueError(f"Qudits should be between 0 and {self.n_qudits()} (number of qudits in the circuit).")
