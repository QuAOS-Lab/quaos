import numpy as np
import pytest

pytest.importorskip("pytket")
from pytket.circuit import Circuit as PytketCircuit, OpType

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES
from sympleq.integrations.quantinuum.utils import to_pytket_circuit, from_pytket_circuit


class TestToPytketCircuit:

    def test_single_qubit_gates(self):
        """Convert circuit with H and S gates."""
        circuit = Circuit.from_tuples([2, 2], [
            (GATES.H, 0),
            (GATES.S, 1),
            (GATES.S_inv, 0),
        ])
        tk = to_pytket_circuit(circuit)

        assert tk.n_qubits == 2
        commands = tk.get_commands()
        assert len(commands) == 3
        assert commands[0].op.type == OpType.H
        assert commands[1].op.type == OpType.S
        assert commands[2].op.type == OpType.Sdg

    def test_two_qubit_gates(self):
        """Convert circuit with CX, SWAP, and CZ gates."""
        circuit = Circuit.from_tuples([2, 2, 2], [
            (GATES.CX, 0, 1),
            (GATES.SWAP, 1, 2),
            (GATES.CZ, 0, 2),
        ])
        tk = to_pytket_circuit(circuit)

        assert tk.n_qubits == 3
        commands = tk.get_commands()
        assert len(commands) == 3
        assert commands[0].op.type == OpType.CX
        assert commands[1].op.type == OpType.SWAP
        assert commands[2].op.type == OpType.CZ

    def test_qudit_indices_preserved(self):
        """Qudit indices are correctly mapped to pytket qubits."""
        circuit = Circuit.from_tuples([2, 2, 2], [
            (GATES.CX, 2, 0),
            (GATES.H, 1),
        ])
        tk = to_pytket_circuit(circuit)

        commands = tk.get_commands()
        assert [q.index[0] for q in commands[0].qubits] == [2, 0]
        assert [q.index[0] for q in commands[1].qubits] == [1]

    def test_empty_circuit(self):
        """Convert an empty circuit."""
        circuit = Circuit.empty([2, 2])
        tk = to_pytket_circuit(circuit)

        assert tk.n_qubits == 2
        assert len(tk.get_commands()) == 0

    def test_non_qubit_raises(self):
        """Raise ValueError for non-qubit dimensions."""
        circuit = Circuit.from_tuples([2, 3], [(GATES.H, 0)])
        with pytest.raises(ValueError, match="dimension=2"):
            to_pytket_circuit(circuit)

    def test_pauli_gates(self):
        """Convert circuit with X, Y, Z gates."""
        circuit = Circuit.from_tuples([2, 2], [
            (GATES.X, 0),
            (GATES.Y, 1),
            (GATES.Z, 0),
        ])
        tk = to_pytket_circuit(circuit)

        commands = tk.get_commands()
        assert len(commands) == 3
        assert commands[0].op.type == OpType.X
        assert commands[1].op.type == OpType.Y
        assert commands[2].op.type == OpType.Z

    def test_identity_gate(self):
        """Convert circuit with Id gate."""
        circuit = Circuit.from_tuples([2], [(GATES.Id, 0)])
        tk = to_pytket_circuit(circuit)

        commands = tk.get_commands()
        assert len(commands) == 1
        assert commands[0].op.type == OpType.noop

    def test_inverse_gates(self):
        """H_inv and CX_inv map to their pytket equivalents."""
        circuit = Circuit.from_tuples([2, 2], [
            (GATES.H_inv, 0),
            (GATES.CX_inv, 0, 1),
        ])
        tk = to_pytket_circuit(circuit)

        commands = tk.get_commands()
        assert commands[0].op.type == OpType.H
        assert commands[1].op.type == OpType.CX

    def test_zzphase_gates(self):
        """GATES.ZZPhase maps to OpType.ZZPhase with angle +0.5; inverse with -0.5."""
        circuit = Circuit.from_tuples([2, 2], [
            (GATES.ZZPhase, 0, 1),
            (GATES.ZZPhase_inv, 0, 1),
        ])
        tk = to_pytket_circuit(circuit)

        commands = tk.get_commands()
        assert len(commands) == 2
        assert commands[0].op.type == OpType.ZZPhase
        assert np.isclose(float(commands[0].op.params[0]), 0.5)
        assert commands[1].op.type == OpType.ZZPhase
        # pytket normalizes -0.5 modulo the gate's 4-half-turn period, so we accept either form.
        inv_angle = float(commands[1].op.params[0]) % 4.0
        assert np.isclose(inv_angle, 3.5)


class TestFromPytketCircuit:

    def test_single_qubit_gates(self):
        """Convert pytket circuit with H, S, Sdg."""
        tk = PytketCircuit(2)
        tk.H(0)
        tk.S(1)
        tk.add_gate(OpType.Sdg, [0])

        circuit = from_pytket_circuit(tk)

        assert circuit.n_qudits() == 2
        assert np.all(circuit.dimensions == 2)
        assert circuit.n_gates() == 3
        assert circuit.gates[0] is GATES.H
        assert circuit.gates[1] is GATES.S
        assert circuit.gates[2] is GATES.S_inv

    def test_pauli_gates(self):
        """Convert pytket circuit with X, Y, Z."""
        tk = PytketCircuit(2)
        tk.X(0)
        tk.Y(1)
        tk.Z(0)

        circuit = from_pytket_circuit(tk)

        assert circuit.n_gates() == 3
        assert circuit.gates[0] is GATES.X
        assert circuit.gates[1] is GATES.Y
        assert circuit.gates[2] is GATES.Z

    def test_two_qubit_gates(self):
        """Convert pytket circuit with CX, SWAP, CZ."""
        tk = PytketCircuit(3)
        tk.CX(0, 1)
        tk.add_gate(OpType.SWAP, [1, 2])
        tk.CZ(0, 2)

        circuit = from_pytket_circuit(tk)

        assert circuit.n_gates() == 3
        assert circuit.gates[0] is GATES.CX
        assert circuit.gates[1] is GATES.SWAP
        assert circuit.gates[2] is GATES.CZ
        assert circuit.qudit_indices[0] == (0, 1)
        assert circuit.qudit_indices[1] == (1, 2)
        assert circuit.qudit_indices[2] == (0, 2)

    def test_empty_circuit(self):
        """Convert an empty pytket circuit."""
        tk = PytketCircuit(3)
        circuit = from_pytket_circuit(tk)

        assert circuit.n_qudits() == 3
        assert circuit.n_gates() == 0

    def test_unsupported_gate_raises(self):
        """Raise ValueError for unsupported pytket gate."""
        tk = PytketCircuit(1)
        tk.T(0)

        with pytest.raises(ValueError, match="No SympleQ mapping"):
            from_pytket_circuit(tk)

    def test_zzphase_gates(self):
        """ZZPhase(±0.5) should map to GATES.ZZPhase and GATES.ZZPhase_inv."""
        tk = PytketCircuit(2)
        tk.add_gate(OpType.ZZPhase, [0.5], [0, 1])
        tk.add_gate(OpType.ZZPhase, [-0.5], [0, 1])

        circuit = from_pytket_circuit(tk)

        assert circuit.n_gates() == 2
        assert circuit.gates[0] is GATES.ZZPhase
        assert circuit.gates[1] is GATES.ZZPhase_inv

    def test_zzphase_non_clifford_angle_raises(self):
        """Non-Clifford ZZPhase angles should raise ValueError."""
        tk = PytketCircuit(2)
        tk.add_gate(OpType.ZZPhase, [0.25], [0, 1])

        with pytest.raises(ValueError, match="Clifford ZZPhase"):
            from_pytket_circuit(tk)


class TestRoundtrip:

    def test_sympleq_to_pytket_and_back(self):
        """Roundtrip: SympleQ -> pytket -> SympleQ preserves gates."""
        original = Circuit.from_tuples([2, 2, 2], [
            (GATES.H, 0),
            (GATES.S, 1),
            (GATES.CX, 0, 1),
            (GATES.SWAP, 1, 2),
            (GATES.CZ, 0, 2),
            (GATES.S_inv, 2),
        ])

        restored = from_pytket_circuit(to_pytket_circuit(original))

        assert original.n_gates() == restored.n_gates()
        for i in range(original.n_gates()):
            assert original.gates[i] is restored.gates[i]
            assert original.qudit_indices[i] == restored.qudit_indices[i]

    def test_pytket_to_sympleq_and_back(self):
        """Roundtrip: pytket -> SympleQ -> pytket preserves structure."""
        tk_original = PytketCircuit(3)
        tk_original.H(0)
        tk_original.S(1)
        tk_original.CX(0, 2)
        tk_original.add_gate(OpType.Sdg, [1])

        tk_restored = to_pytket_circuit(from_pytket_circuit(tk_original))

        orig_cmds = tk_original.get_commands()
        rest_cmds = tk_restored.get_commands()
        assert len(orig_cmds) == len(rest_cmds)
        for oc, rc in zip(orig_cmds, rest_cmds):
            assert oc.op.type == rc.op.type
            assert [q.index[0] for q in oc.qubits] == [q.index[0] for q in rc.qubits]

    def test_zzphase_roundtrip(self):
        """SympleQ -> pytket -> SympleQ roundtrip preserves ZZPhase and its inverse."""
        original = Circuit.from_tuples([2, 2], [
            (GATES.ZZPhase, 0, 1),
            (GATES.ZZPhase_inv, 0, 1),
        ])
        restored = from_pytket_circuit(to_pytket_circuit(original))

        assert restored.n_gates() == 2
        assert restored.gates[0] is GATES.ZZPhase
        assert restored.gates[1] is GATES.ZZPhase_inv

    def test_random_circuit_roundtrip(self):
        """Roundtrip random qubit circuits preserves gate count and gate set."""
        for _ in range(20):
            n_qudits = np.random.randint(2, 6)
            n_gates = np.random.randint(0, 15)
            original = Circuit.from_random(n_gates, [2] * n_qudits)

            restored = from_pytket_circuit(to_pytket_circuit(original))

            # pytket may reorder independent gates, so just check counts match
            assert original.n_gates() == restored.n_gates()
            assert original.n_qudits() == restored.n_qudits()
