import numpy as np
import pytest

from pytket.circuit import Circuit as PytketCircuit, OpType

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES
from sympleq.core.paulis import PauliString
from sympleq.integrations.quantinuum.utils import to_pytket_circuit, from_pytket_circuit, NATIVE_GATES_SET


def _action_signature(circuit: Circuit) -> list[PauliString]:
    """Conjugation action of a qubit Clifford circuit on all X/Z generators."""
    n = circuit.n_qudits()
    signature = []
    for i in range(n):
        for label in ("x1z0", "x0z1"):
            s = " ".join(label if j == i else "x0z0" for j in range(n))
            signature.append(circuit.act(PauliString.from_string(s, dimensions=[2] * n)))
    return signature


_H_SERIES_NATIVE_OPS = {OpType.Rz, OpType.PhasedX, OpType.ZZPhase, OpType.Measure, OpType.Barrier}


def _assert_native(tk: PytketCircuit):
    """Assert that a compiled circuit contains only H-series native ops."""
    for command in tk.get_commands():
        assert command.op.type in _H_SERIES_NATIVE_OPS, f"non-native op {command.op.type}"


class TestToPytketCircuit:
    """``to_pytket_circuit`` measures all qubits and compiles with the Nexus
    default pass, so outputs contain only H-series native ops."""

    def test_single_qubit_gates(self):
        """H/S/Sdg compile to native Rz/PhasedX ops."""
        circuit = Circuit.from_tuples([2, 2], [
            (GATES.H, 0),
            (GATES.S, 1),
            (GATES.S_inv, 0),
        ])
        tk = to_pytket_circuit(circuit)

        assert tk.n_qubits == 2
        _assert_native(tk)

    def test_two_qubit_gates(self):
        """CX/SWAP/CZ compile to native ZZPhase-based ops."""
        circuit = Circuit.from_tuples([2, 2, 2], [
            (GATES.CX, 0, 1),
            (GATES.SWAP, 1, 2),
            (GATES.CZ, 0, 2),
        ])
        tk = to_pytket_circuit(circuit)

        assert tk.n_qubits == 3
        _assert_native(tk)

    def test_all_qubits_measured(self):
        """measure_all() runs before compilation: every qubit gets a Measure."""
        circuit = Circuit.from_tuples([2, 2, 2], [(GATES.H, 0)])
        tk = to_pytket_circuit(circuit)

        measured = sorted(
            cmd.qubits[0].index[0]
            for cmd in tk.get_commands()
            if cmd.op.type == OpType.Measure
        )
        assert measured == [0, 1, 2]

    def test_empty_circuit(self):
        """An empty circuit compiles to measurements only."""
        circuit = Circuit.empty([2, 2])
        tk = to_pytket_circuit(circuit)

        assert tk.n_qubits == 2
        assert all(cmd.op.type == OpType.Measure for cmd in tk.get_commands())

    def test_non_qubit_raises(self):
        """Raise ValueError for non-qubit dimensions."""
        circuit = Circuit.from_tuples([2, 3], [(GATES.H, 0)])
        with pytest.raises(ValueError, match="dimension=2"):
            to_pytket_circuit(circuit)

    def test_ZZMax_gates(self):
        """GATES.ZZMax compiles to ZZPhase(0.5); GATES.ZZMax_inv to ZZPhase(3.5)."""
        circuit = Circuit.from_tuples([2, 2], [
            (GATES.ZZMax, 0, 1),
            (GATES.ZZMax_inv, 0, 1),
        ])
        tk = to_pytket_circuit(circuit)

        zz = [cmd for cmd in tk.get_commands() if cmd.op.type == OpType.ZZPhase]
        assert len(zz) == 2
        assert np.isclose(float(zz[0].op.params[0]) % 4.0, 0.5)
        assert np.isclose(float(zz[1].op.params[0]) % 4.0, 3.5)

    def test_V_gates(self):
        """GATES.V compiles to PhasedX(0.5, 0); GATES.V_inv to PhasedX(3.5, 0)."""
        circuit = Circuit.from_tuples([2], [
            (GATES.V, 0),
            (GATES.V_inv, 0),
        ])
        tk = to_pytket_circuit(circuit)

        px = [cmd for cmd in tk.get_commands() if cmd.op.type == OpType.PhasedX]
        assert len(px) == 2
        assert np.isclose(float(px[0].op.params[0]) % 4.0, 0.5)
        assert np.isclose(float(px[0].op.params[1]) % 4.0, 0.0)
        assert np.isclose(float(px[1].op.params[0]) % 4.0, 3.5)
        assert np.isclose(float(px[1].op.params[1]) % 4.0, 0.0)

    def test_NATIVE_GATES_SET_all_convert(self):
        """Every gate in NATIVE_GATES_SET compiles to native ops without error."""
        for gate in NATIVE_GATES_SET:
            qudits = tuple(range(gate.n_qudits))
            dims = [2] * max(2, gate.n_qudits)
            circuit = Circuit.from_tuples(dims, [(gate, *qudits)])
            tk = to_pytket_circuit(circuit)
            _assert_native(tk)


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

    def test_ZZMax_gates(self):
        """OpType.ZZMax (no params) maps to GATES.ZZMax; ZZPhase(-0.5) maps to GATES.ZZMax_inv."""
        tk = PytketCircuit(2)
        tk.add_gate(OpType.ZZMax, [0, 1])
        tk.add_gate(OpType.ZZPhase, [-0.5], [0, 1])

        circuit = from_pytket_circuit(tk)

        assert circuit.n_gates() == 2
        assert circuit.gates[0] is GATES.ZZMax
        assert circuit.gates[1] is GATES.ZZMax_inv

    def test_ZZPhase_non_clifford_angle_raises(self):
        """Non-Clifford ZZPhase angles should raise ValueError."""
        tk = PytketCircuit(2)
        tk.add_gate(OpType.ZZPhase, [0.25], [0, 1])

        with pytest.raises(ValueError, match="Clifford ZZPhase"):
            from_pytket_circuit(tk)


class TestRoundtrip:

    def test_sympleq_to_pytket_and_back(self):
        """Roundtrip: SympleQ -> pytket -> SympleQ preserves the Clifford action."""
        original = Circuit.from_tuples([2, 2, 2], [
            (GATES.H, 0),
            (GATES.S, 1),
            (GATES.CX, 0, 1),
            (GATES.SWAP, 1, 2),
            (GATES.CZ, 0, 2),
            (GATES.S_inv, 2),
        ])

        restored = from_pytket_circuit(to_pytket_circuit(original))

        assert restored.n_qudits() == original.n_qudits()
        assert _action_signature(restored) == _action_signature(original)

    def test_pytket_to_sympleq_and_back(self):
        """Roundtrip: pytket -> SympleQ -> pytket preserves the Clifford action."""
        tk_original = PytketCircuit(3)
        tk_original.H(0)
        tk_original.S(1)
        tk_original.CX(0, 2)
        tk_original.add_gate(OpType.Sdg, [1])

        s_original = from_pytket_circuit(tk_original)
        s_restored = from_pytket_circuit(to_pytket_circuit(s_original))

        assert _action_signature(s_restored) == _action_signature(s_original)

    def test_every_mapped_gate_roundtrips(self):
        """Each supported SympleQ gate survives compile + back-conversion."""
        one_qubit = [GATES.Id, GATES.H, GATES.H_inv, GATES.S, GATES.S_inv,
                     GATES.X, GATES.X_inv, GATES.Y, GATES.Y_inv,
                     GATES.Z, GATES.Z_inv, GATES.V, GATES.V_inv]
        two_qubit = [GATES.CX, GATES.CX_inv, GATES.SWAP, GATES.CZ,
                     GATES.ZZMax, GATES.ZZMax_inv]
        for gate in one_qubit:
            original = Circuit.from_tuples([2, 2], [(gate, 0)])
            restored = from_pytket_circuit(to_pytket_circuit(original))
            assert _action_signature(restored) == _action_signature(original), gate.name
        for gate in two_qubit:
            original = Circuit.from_tuples([2, 2], [(gate, 0, 1)])
            restored = from_pytket_circuit(to_pytket_circuit(original))
            assert _action_signature(restored) == _action_signature(original), gate.name

    def test_native_gate_counts_preserved(self):
        """NATIVE_GATES_SET gates (plus Id) rebase 1:1, so 1- and 2-qubit gate
        counts survive compilation and the roundtrip."""
        rng = np.random.default_rng(11)
        gate_pool = NATIVE_GATES_SET + [GATES.Id]
        for _ in range(5):
            n_qubits = int(rng.integers(2, 5))
            tuples = []
            for _ in range(15):
                gate = gate_pool[rng.integers(len(gate_pool))]
                if gate.n_qudits == 2:
                    a, b = rng.choice(n_qubits, size=2, replace=False)
                    tuples.append((gate, int(a), int(b)))
                else:
                    tuples.append((gate, int(rng.integers(n_qubits))))
            original = Circuit.from_tuples([2] * n_qubits, tuples)

            tk = to_pytket_circuit(original)
            assert tk.n_1qb_gates() == original.n_1qd_gates()
            assert tk.n_2qb_gates() == original.n_2qd_gates()

            restored = from_pytket_circuit(tk)
            assert restored.n_1qd_gates() == original.n_1qd_gates()
            assert restored.n_2qd_gates() == original.n_2qd_gates()

    def test_ZZMax_roundtrip(self):
        """SympleQ -> pytket -> SympleQ roundtrip preserves ZZMax and its inverse."""
        original = Circuit.from_tuples([2, 2], [
            (GATES.ZZMax, 0, 1),
            (GATES.ZZMax_inv, 0, 1),
        ])
        restored = from_pytket_circuit(to_pytket_circuit(original))

        assert restored.n_gates() == 2
        assert restored.gates[0] is GATES.ZZMax
        assert restored.gates[1] is GATES.ZZMax_inv

    def test_random_circuit_roundtrip(self):
        """Roundtrip random qubit circuits preserves the Clifford action."""
        for _ in range(20):
            n_qudits = np.random.randint(2, 6)
            n_gates = np.random.randint(0, 15)
            original = Circuit.from_depth(n_gates, [2] * n_qudits)

            restored = from_pytket_circuit(to_pytket_circuit(original))

            assert restored.n_qudits() == original.n_qudits()
            assert _action_signature(restored) == _action_signature(original)


class TestImplicitPermutation:
    """Compilation with ``allow_swaps=True`` absorbs SWAPs into an implicit
    wire permutation; ``from_pytket_circuit`` must restore it explicitly."""

    def test_swap_only_circuit(self):
        """A bare SWAP is fully absorbed into the implicit permutation."""
        original = Circuit.from_tuples([2, 2], [(GATES.SWAP, 0, 1)])
        restored = from_pytket_circuit(to_pytket_circuit(original))

        assert _action_signature(restored) == _action_signature(original)

    def test_swap_cycle(self):
        """Two chained SWAPs produce a 3-cycle wire permutation."""
        original = Circuit.from_tuples([2, 2, 2], [
            (GATES.SWAP, 0, 1),
            (GATES.SWAP, 1, 2),
        ])
        restored = from_pytket_circuit(to_pytket_circuit(original))

        assert _action_signature(restored) == _action_signature(original)

    def test_random_circuits_with_swaps(self):
        """Compiled circuits with interleaved SWAPs stay Clifford-equivalent."""
        rng = np.random.default_rng(7)
        one_qubit = [GATES.H, GATES.S, GATES.S_inv, GATES.X, GATES.Z, GATES.V, GATES.V_inv]
        two_qubit = [GATES.CX, GATES.CZ, GATES.SWAP, GATES.ZZMax]
        for _ in range(10):
            n_qubits = int(rng.integers(2, 5))
            tuples = []
            for _ in range(12):
                if rng.random() < 0.5:
                    gate = two_qubit[rng.integers(len(two_qubit))]
                    a, b = rng.choice(n_qubits, size=2, replace=False)
                    tuples.append((gate, int(a), int(b)))
                else:
                    gate = one_qubit[rng.integers(len(one_qubit))]
                    tuples.append((gate, int(rng.integers(n_qubits))))
            original = Circuit.from_tuples([2] * n_qubits, tuples)

            restored = from_pytket_circuit(to_pytket_circuit(original))

            assert _action_signature(restored) == _action_signature(original)
