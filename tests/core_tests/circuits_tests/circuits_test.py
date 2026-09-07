import json
from pathlib import Path
import tempfile
import numpy as np
import pytest
from scipy.sparse import issparse
from sympleq.core.circuits import Circuit, GATES
from sympleq.core.noise.noise_model import CompositeNoise, DepolarizingNoise, DephasingNoise, Noiseless, NoiseModel
from sympleq.core.paulis import PauliSum


class TestCircuits():

    def test_circuit_composition(self):
        # TODO: Full test for mixed dimensions
        for _ in range(10):
            n_qudits = 3
            dimensions = [2] * n_qudits
            n_paulis = 5
            # make a random circuit
            circuit = Circuit.from_random_symplectic(dimensions)

            # make a random pauli sum
            pauli_sum = PauliSum.from_random(n_paulis, dimensions)
            # compose the circuit and pauli sum
            composed_gate = circuit.composite_gate()
            # For composite gate, act on all qudits
            output_composite = composed_gate.act(pauli_sum, tuple(range(n_qudits)))
            output_sequential = circuit.act(pauli_sum)

            # show that the composed gate returns the same thing as the circuit when acting on the pauli sum
            assert output_composite == output_sequential, (
                f'Input: \n {pauli_sum}\n'
                f'Composed gate:\n{output_composite} \n'
                f'Sequential gate:\n{output_sequential}'
            )

    def test_hadamard_composition(self):
        # simple case of two Hadamards on different qubits. Known symplectic in this case.

        n_qudits = 1
        dimensions = [2] * n_qudits
        n_paulis = 2
        # make a random circuit
        circuit = Circuit.empty(dimensions)
        circuit.add_gate(GATES.H, 0)
        circuit.add_gate(GATES.S, 0)

        # make a random pauli sum
        pauli_sum = PauliSum.from_random(n_paulis, dimensions)

        # compose the circuit and pauli sum
        composed_gate = circuit.composite_gate()
        output_composite = composed_gate.act(pauli_sum, tuple(range(n_qudits)))
        output_sequential = circuit.act(pauli_sum)

        print(output_composite)
        print(output_sequential)

        # show that the composed gate returns the same thing as the circuit when acting on the pauli sum
        assert output_composite == output_sequential

    def test_random_circuit(self):
        # test that a random circuit can be generated with the correct dimensions on mixed qudits
        for _ in range(100):
            n_qudits = np.random.randint(2, 10)
            dimensions = np.random.choice([2, 3, 5], size=n_qudits)
            C = Circuit.from_random_symplectic(dimensions=dimensions)
            ps = PauliSum.from_random(10, dimensions)
            out = C.act(ps)
            assert np.all(out.dimensions == dimensions)

    def test_from_random_preserves_mixed_dimensions(self):
        """Random mixed-dimension circuit should preserve the input qudit register dimensions."""
        rng = np.random.default_rng(123)
        dimensions = np.array([2, 3, 2, 5, 3], dtype=int)

        circuit = Circuit.from_random_symplectic(dimensions=dimensions, rng=rng)

        assert np.array_equal(circuit.dimensions, dimensions)

    @pytest.mark.parametrize("d", [2, 3, 5, 11])
    def test_single_hadamard_unitary(self, d: int):
        # For a single-qudit circuit with one Hadamard, the circuit unitary
        # should equal the gate's local unitary.
        circuit = Circuit.from_gates_and_qudits([d], [GATES.H], [(0,)])
        U_circ = circuit.unitary()
        assert issparse(U_circ)
        U_gate = GATES.H.local_unitary(d)
        assert U_circ.shape == U_gate.shape
        assert np.allclose(U_circ.toarray(), U_gate.toarray())

    def test_mixed_qudits_phase_with_unitary(self):
        N = 100
        dimensions = [2, 3, 5]
        n_paulis = 1
        for _ in range(N):
            P = PauliSum.from_random(n_paulis, dimensions, rand_weights=False)
            C = Circuit.from_random_symplectic(dimensions=dimensions)
            U = C.unitary()

            ps_m = P.to_hilbert_space()
            ps_res = C.act(P)
            ps_res_m = ps_res.to_hilbert_space()
            phase_symplectic = ps_res.phases[0]

            ps_res.reset_phases()
            ps_res_m = ps_res.to_hilbert_space().toarray()
            ps_m_res = (U @ ps_m @ U.conj().T).toarray()
            mask = (ps_res_m != 0)
            factors = np.unique(np.around(ps_m_res[mask] / ps_res_m[mask], 10))
            assert len(factors) == 1
            factor = factors[0]
            d = P.lcm
            phase_unitary = int(np.around((d * np.angle(factor) / (np.pi)) % (2 * d), 1))
            assert np.array_equal(phase_symplectic, phase_unitary)

    def test_phase_mixed_species(self):
        def debug_steps(C: Circuit, P: PauliSum):
            print(f"Initial phases: {P.phases} -- exponents: {P.tableau}")
            for i, partial_p in enumerate(C.act_iter(P)):
                gate = C.gates[i]
                print(f"Phases after {gate.name}: {partial_p.phases} -- exponents: {partial_p.tableau}")

        # Test 1: Simple qutrit + qubit
        P = PauliSum.from_string(['x2z0 x0z0'],
                                 dimensions=[3, 2],
                                 weights=[1], phases=[0])
        idx = 0
        C = Circuit.from_tuples(dimensions=P.dimensions, data=[(GATES.S, idx)])
        debug_steps(C, P)
        P = C.act(P)
        assert P.phases[0] == 4

        # Test 2: Simple ququint + qubit
        P = PauliSum.from_string(['x2z0 x0z0'],
                                 dimensions=[5, 2],
                                 weights=[1], phases=[0])

        idx = 0
        C = Circuit.from_tuples(dimensions=P.dimensions, data=[(GATES.S, idx)])
        debug_steps(C, P)
        P = C.act(P)
        assert P.phases[0] == 4

        # Test 3: More complex ququint + qubit
        P = PauliSum.from_string(['x3z0 x0z0'],
                                 dimensions=[5, 2],
                                 weights=[1], phases=[0])
        idx = 0
        C = Circuit.from_tuples(dimensions=P.dimensions, data=[(GATES.S, idx)])
        debug_steps(C, P)
        P = C.act(P)
        assert P.phases[0] == 12

        # Test 4: Simple ququint + qutrit
        P = PauliSum.from_string(['x2z0 x0z0'],
                                 dimensions=[5, 3],
                                 weights=[1], phases=[0])
        idx = 0
        C = Circuit.from_tuples(dimensions=P.dimensions, data=[(GATES.S, idx)])
        debug_steps(C, P)
        P = C.act(P)
        assert P.phases[0] == 6

        # Test 5: Simple qutrit + qubit but action on qubit
        P = PauliSum.from_string(['x0z0 x1z0'],
                                 dimensions=[3, 2],
                                 weights=[1], phases=[0])
        idx = 1
        C = Circuit.from_tuples(dimensions=P.dimensions, data=[(GATES.S, idx)])
        debug_steps(C, P)
        P = C.act(P)
        assert P.phases[0] == 3

        # Test 6: Simple ququint + qutrit + qubit
        P = PauliSum.from_string(['x2z0 x0z0 x0z0'],
                                 dimensions=[5, 3, 2],
                                 weights=[1], phases=[0])
        idx = 0
        C = Circuit.from_tuples(dimensions=P.dimensions, data=[(GATES.S, idx)])
        debug_steps(C, P)
        P = C.act(P)
        assert P.phases[0] == 12

        # Test 7: composite circuit
        P = PauliSum.from_string(['x2z2 x0z0'],
                                 dimensions=[3, 2],
                                 weights=[1], phases=[0])
        idx = 0
        C = Circuit.from_tuples(dimensions=P.dimensions, data=[(GATES.S, idx), (GATES.S, idx), (GATES.H, idx)])
        debug_steps(C, P)
        P = C.act(P)
        assert P.phases[0] == 8

    @staticmethod
    def _linear_index(dims, idxs):
        # Row-major: idx = sum_k idxs[k] * prod_{l>k} dims[l]
        strides = [1] * len(dims)
        for k in range(len(dims) - 2, -1, -1):
            strides[k] = strides[k + 1] * dims[k + 1]
        return sum(idxs[k] * strides[k] for k in range(len(dims)))

    def test_swap_embedding_on_equal_dims(self):
        # Verify SWAP on qudits (0,1) within a 2-qudit system with equal dimensions.
        dims = [3, 3]
        c = Circuit.from_gates_and_qudits(dims, [GATES.SWAP], [(0, 1)])
        U = c.unitary()

        # Start in |i,j> with i=1, j=2
        i, j = 1, 2
        D = np.prod(dims)
        psi = np.zeros(D, dtype=complex)
        psi[self._linear_index(dims, [i, j])] = 1.0

        # Expected after SWAP: |j,i>
        phi = U @ psi
        expected = np.zeros(D, dtype=complex)
        expected[self._linear_index(dims, [j, i])] = 1.0

        assert np.allclose(phi, expected), f"Expected:\n{expected}\nGot:\n{phi}"

    def test_cx_embedding_on_three_qudits(self):
        # Verify CX on qudits (1,2) inside a 3-qudit system.
        d = 5
        dims = [d, d, d]
        c = Circuit.from_gates_and_qudits(dims, [GATES.CX], [(1, 2)])
        U = c.unitary()

        # Start in |i,j,k> = |3,1,4>
        i, j, k = 3, 1, 4
        D = np.prod(dims)
        psi = np.zeros(D, dtype=complex)
        psi[self._linear_index(dims, [i, j, k])] = 1.0

        # After CX(1->2): |i, j, k+j mod d>
        phi = U @ psi
        expected = np.zeros(D, dtype=complex)
        expected[self._linear_index(dims, [i, j, (k + j) % d])] = 1.0

        assert np.allclose(phi, expected)

    def test_phase_embedding_on_middle_qudit(self):
        # Verify PHASE acting on middle qudit multiplies amplitude appropriately.
        d0, d1, d2 = 3, 5, 2
        dims = [d0, d1, d2]
        c = Circuit.from_tuples(dims, (GATES.S, 1))
        U = c.unitary().toarray()

        # Basis |i,j,k> = |2,3,1>
        i, j, k = 2, 3, 1
        D = np.prod(dims)
        psi = np.zeros(D, dtype=complex)
        psi[self._linear_index(dims, [i, j, k])] = 1.0

        phi = U @ psi

        # PHASE unitary is diag(omega^{j(j-1)/2}) on that qudit, with omega = exp(2πi/d).
        omega = np.exp(1j * 2 * np.pi / d1)
        factor = omega ** (j * (j - 1) / 2)
        expected = np.zeros(D, dtype=complex)
        expected[self._linear_index(dims, [i, j, k])] = factor

        assert np.allclose(phi, expected)

    @pytest.mark.parametrize("dimension", [2, 3, 5])
    @pytest.mark.parametrize("n_qudits", [2, 3, 4])
    @pytest.mark.parametrize("c_depth", [1, 3, 5])
    def test_circuit_unitary(self, dimension: int, n_qudits: int, c_depth: int):
        n_tests = 100
        n_paulis = 10
        dimensions = [dimension] * n_qudits
        for _ in range(n_tests):
            ps = PauliSum.from_random(n_paulis, dimensions, False)
            c = Circuit.from_depth(c_depth, dimensions)
            U_c = c.unitary()
            P_from_conjugation = U_c @ ps.to_hilbert_space() @ U_c.conj().T
            P_from_act = c.act(ps).to_hilbert_space()
            diff_m = np.around(P_from_conjugation - P_from_act.toarray(), 10)
            assert not np.any(diff_m), f'failed for dimension {_, dimension, c.__str__()}'

    def test_to_string_basic(self):
        """Test basic to_string functionality."""
        dimensions = [2] * 3
        circuit = Circuit.from_tuples(dimensions, [(GATES.H, 0), (GATES.CX, 0, 1)])
        s = circuit.to_string()
        data = json.loads(s)

        assert data["dimensions"] == dimensions
        assert data["data"] == [["H", [0], "None"], ["CX", [0, 1], "None"]]

    def test_from_string_basic(self):
        """Test basic from_string functionality."""
        s = '{"dimensions": [2, 2], "data": [["H", [0]], ["CX", [0, 1]]]}'
        circuit = Circuit.from_string(s)

        assert len(circuit.gates) == 2
        assert circuit.gates[0] is GATES.H
        assert circuit.gates[1] is GATES.CX
        assert circuit.qudit_indices[0] == (0,)
        assert circuit.qudit_indices[1] == (0, 1)

    def test_roundtrip_simple(self):
        """Test that to_string -> from_string produces equivalent circuit."""
        dimensions = [2] * 3
        original = Circuit.from_tuples(dimensions, [(GATES.H, 0), (GATES.S, 1), (GATES.CX, 0, 1), (GATES.SWAP, 1, 2)])
        s = original.to_string()
        restored = Circuit.from_string(s)

        assert original == restored

    def test_roundtrip_all_gates(self):
        """Test roundtrip with all supported gate types."""
        dimensions = [2] * 3
        original = Circuit.from_tuples(
            dimensions,
            [
                (GATES.H, 0),
                (GATES.H_inv, 1),
                (GATES.S, 0),
                (GATES.S_inv, 1),
                (GATES.CX, 0, 1),
                (GATES.CX_inv, 1, 0),
                (GATES.SWAP, 0, 1),
                (GATES.CZ, 0, 1),
            ]
        )
        s = original.to_string()
        restored = Circuit.from_string(s)

        assert original == restored

    def test_roundtrip_mixed_dimensions(self):
        """Test roundtrip with mixed qudit dimensions."""
        dimensions = [2] * 3
        original = Circuit.from_tuples(
            dimensions,
            [(GATES.H, 0), (GATES.S, 1), (GATES.H, 2)]
        )
        s = original.to_string()
        restored = Circuit.from_string(s)

        assert original == restored

    def test_roundtrip_empty_circuit(self):
        """Test roundtrip with empty circuit."""
        dimensions = [2] * 3
        original = Circuit.empty(dimensions)
        s = original.to_string()
        restored = Circuit.from_string(s)

        assert original == restored
        assert len(restored.gates) == 0

    def test_file_roundtrip(self):
        """Test save_to_file and from_file."""
        dimensions = [2] * 3
        original = Circuit.from_tuples(dimensions, [(GATES.H, 0), (GATES.S, 1), (GATES.CX, 0, 1)])

        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            temp_path = Path(f.name)

        try:
            original.save_to_file(temp_path)
            restored = Circuit.from_file(temp_path)
            assert original == restored
        finally:
            temp_path.unlink()

    def test_from_string_unknown_gate_raises(self):
        """Test that unknown gate names raise ValueError."""
        s = '{"dimensions": [2], "data": [["UNKNOWN_GATE", [0]]]}'
        with pytest.raises(ValueError, match="Unknown gate name"):
            Circuit.from_string(s)

    def test_roundtrip_random_circuits(self):
        """Test roundtrip with random circuits."""
        for _ in range(50):
            depth = np.random.randint(0, 20)
            n_qudits = np.random.randint(2, 10)
            dimensions = np.random.choice([2, 3, 5], size=n_qudits)
            original = Circuit.from_depth(depth, dimensions)

            s = original.to_string()
            restored = Circuit.from_string(s)

            assert original == restored

    def test_roundtrip_preserves_behavior(self):
        """Test that restored circuit produces same results when acting on Paulis."""
        dimensions = [2, 3, 5]
        original = Circuit.from_depth(10, dimensions)
        restored = Circuit.from_string(original.to_string())

        pauli_sum = PauliSum.from_random(5, dimensions)
        original_result = original.act(pauli_sum)
        restored_result = restored.act(pauli_sum)

        assert original_result == restored_result

    def test_roundtrip_with_global_noise(self):
        """Test roundtrip with a single noise model applied to all gates."""
        dimensions = [2, 2]
        circuit = Circuit.from_tuples(
            dimensions, [(GATES.H, 0), (GATES.CX, 0, 1)]
        ).with_noise(DepolarizingNoise(0.05))

        restored = Circuit.from_string(circuit.to_string())
        assert circuit == restored

    def test_roundtrip_with_per_gate_noise(self):
        """Test roundtrip with different noise models per gate."""
        dimensions = [2, 2]
        circuit = Circuit.from_tuples(
            dimensions, [(GATES.H, 0), (GATES.S, 1), (GATES.CX, 0, 1)]
        ).with_noise([DepolarizingNoise(0.05), DephasingNoise(0.1), None])

        restored = Circuit.from_string(circuit.to_string())
        assert circuit == restored

    def test_roundtrip_with_noiseless(self):
        """Test roundtrip with Noiseless noise model."""
        dimensions = [2, 2]
        circuit = Circuit.from_tuples(
            dimensions, [(GATES.H, 0), (GATES.CX, 0, 1)]
        ).with_noise(Noiseless())

        restored = Circuit.from_string(circuit.to_string())
        assert circuit == restored

    def test_noise_model_from_string_none(self):
        """Test NoiseModel.from_string with None."""
        assert NoiseModel.from_string("None") is None

    def test_noise_model_from_string_unknown(self):
        """Test NoiseModel.from_string raises on unknown model."""
        with pytest.raises(ValueError, match="Unknown noise model"):
            NoiseModel.from_string("UnknownNoise(error_rate=0.1)")

    def test_noise_model_from_string_invalid(self):
        """Test NoiseModel.from_string raises on unparseable string."""
        with pytest.raises(ValueError, match="Cannot parse"):
            NoiseModel.from_string("garbage")

    def test_roundtrip_with_composite_noise(self):
        """Test roundtrip with CompositeNoise applied to all gates."""
        dimensions = [2, 2]
        noise = CompositeNoise.from_noise_models([
            DepolarizingNoise(0.05), DephasingNoise(0.1)
        ])
        circuit = Circuit.from_tuples(
            dimensions, [(GATES.H, 0), (GATES.CX, 0, 1)]
        ).with_noise(noise)

        restored = Circuit.from_string(circuit.to_string())
        assert circuit == restored

    def test_composite_noise_str(self):
        """Test CompositeNoise.__str__ output."""
        noise = CompositeNoise.from_noise_models([
            DepolarizingNoise(0.05), DephasingNoise(0.1)
        ])
        s = str(noise)
        assert s.startswith("CompositeNoise([")
        assert "DepolarizingNoise" in s
        assert "DephasingNoise" in s

    def test_gates_layout(self):
        dimensions = [3] * 4
        circuit = Circuit.from_depth(np.random.randint(4, 12), dimensions)

        # Just check that it doesn't crash
        _ = circuit.gates_layout()
        _ = circuit.gates_layout(with_qudit_indices=True)
        ps_in = PauliSum.from_random(len(dimensions), dimensions)
        _ = circuit.gates_layout(with_qudit_indices=True, with_input=ps_in)
        ps_out = PauliSum.from_random(len(dimensions), dimensions)
        _ = circuit.gates_layout(with_qudit_indices=True, with_input=ps_in, with_output=ps_out)
        _ = circuit.gates_layout(with_qudit_indices=True, with_input=ps_in, with_output=ps_out, wrap=False)

    def test_gates_layout_exception(self):
        dimensions = [3] * 4
        circuit = Circuit.from_depth(np.random.randint(4, 12), dimensions)

        # Just check that it doesn't crash
        _ = circuit.gates_layout(wires="-")
        _ = circuit.gates_layout(wires=["x"] * circuit.n_qudits())

        with pytest.raises(ValueError):
            _ = circuit.gates_layout(wires=["x"] * (circuit.n_qudits() + 1))

    def test_circuit_sanity_checks_exception(self):
        # Inconsistent dimensions
        dimensions = [2, 3, 3, 1]
        with pytest.raises(ValueError):
            _ = Circuit.from_depth(4, dimensions)

        # Inconsistent dimensions for CX
        dimensions = [2, 3, 3, 5]
        with pytest.raises(ValueError):
            _ = Circuit.from_gates_and_qudits(dimensions, [GATES.CX, GATES.H], [(0, 1), (2,)])

        with pytest.raises(ValueError):
            _ = Circuit.from_gates_and_qudits(dimensions, [GATES.CX, GATES.H], [(0,), (2,)])

        with pytest.raises(ValueError):
            _ = Circuit.from_gates_and_qudits(dimensions, [GATES.CX, GATES.H], [(0, 0), (2,)])

        # Inconsistent gates and qudit indices
        dimensions = [2, 3, 3, 5]
        with pytest.raises(ValueError):
            _ = Circuit.from_gates_and_qudits(dimensions, [GATES.CX, GATES.H], [(0, 1), (2,), (1,)])
        with pytest.raises(ValueError):
            _ = Circuit.from_gates_and_qudits(dimensions, [GATES.CX, GATES.H, GATES.S], [(0, 1), (2,)])

    def test_from_random_gates_set_default(self):
        """Default gate set is {H, S, CX, CZ, SWAP}; no ZZMax should appear."""
        rng = np.random.default_rng(0)
        c = Circuit.from_depth(50, [2, 2, 2], rng=rng)
        allowed = {GATES.H, GATES.S, GATES.CX, GATES.CZ, GATES.SWAP}
        for g in c.gates:
            assert g in allowed, f"Default gate set produced unexpected gate {g.name}"

    def test_from_random_gates_set_restricts_sampling(self):
        """Only the gates in gates_set should appear in the produced circuit."""
        gates_set = [GATES.S, GATES.ZZMax]
        rng = np.random.default_rng(42)
        c = Circuit.from_depth(100, [2] * 4, gates_set=gates_set,
                               two_qudit_gate_ratio=0.5, rng=rng)
        seen = set(c.gates)
        assert seen.issubset({GATES.S, GATES.ZZMax})
        # With 100 gates and 50/50 ratio both should actually appear.
        assert GATES.S in seen
        assert GATES.ZZMax in seen

    def test_from_random_gates_set_only_single_qudit(self):
        """With only single-qudit gates and ratio=0, only single-qudit gates should appear."""
        gates_set = [GATES.H, GATES.S]
        rng = np.random.default_rng(0)
        c = Circuit.from_depth(40, [2, 2, 2], gates_set=gates_set,
                               two_qudit_gate_ratio=0.0, rng=rng)
        for g in c.gates:
            assert g.n_qudits == 1
            assert g in {GATES.H, GATES.S}

    def test_from_random_gates_set_warns_on_multi_qudit(self):
        """Passing a gate that acts on >2 qudits should emit a warning (and be ignored)."""
        from sympleq.core.circuits import Gate
        # A 3-qudit random Clifford as a stand-in for an unsupported multi-qudit gate
        three_qudit_gate = Gate.from_random(n_qudits=3, dimension=2)
        gates_set = [GATES.H, GATES.CX, three_qudit_gate]
        with pytest.warns(UserWarning, match="single qudit and 2-qudits"):
            c = Circuit.from_depth(10, [2, 2, 2], gates_set=gates_set,
                                   rng=np.random.default_rng(0))
        for g in c.gates:
            assert g.n_qudits <= 2
            assert g is not three_qudit_gate

    def test_from_random_gates_set_invalid_type(self):
        """A non-container gates_set should raise ValueError."""
        with pytest.raises(ValueError):
            _ = Circuit.from_depth(5, [2, 2], gates_set=GATES.H)  # type: ignore[arg-type]

    def test_from_depth_default_gate_set(self):
        """Default from_depth gate set is {H, S, CX, CZ, SWAP}."""
        rng = np.random.default_rng()
        c = Circuit.from_depth(depth=10, dimensions=[2] * 4, rng=rng)
        allowed = {GATES.H, GATES.S, GATES.CX, GATES.CZ, GATES.SWAP}
        for g in c.gates:
            assert g in allowed, f"Default gate set produced unexpected gate {g.name}"

    def test_from_depth_single_qudit_only_covers_full_grid(self):
        """With ratio=0 every layer is filled with single-qudit gates -> n_gates == depth * n_qudits."""
        depth = 4
        dimensions = [2, 2, 2]
        c = Circuit.from_depth(depth=depth, dimensions=dimensions,
                               two_qudit_gate_ratio=0.0,
                               rng=np.random.default_rng(0))
        assert c.n_gates() == depth * len(dimensions)
        for g in c.gates:
            assert g.n_qudits == 1

    def test_from_depth_qudit_positions_cover_full_grid(self):
        """Every layer fully tiles the qudit register: total qudit-slots == depth * n_qudits."""
        depth = 6
        dimensions = [2, 2, 2, 2]
        c = Circuit.from_depth(depth=depth, dimensions=dimensions,
                               rng=np.random.default_rng(7))
        total_positions = sum(len(idxs) for idxs in c.qudit_indices)
        assert total_positions == depth * len(dimensions)

    def test_from_depth_gates_set_restricts_sampling(self):
        """Only gates in gates_set should appear."""
        gates_set = [GATES.H, GATES.CX]
        c = Circuit.from_depth(depth=10, dimensions=[2] * 3, gates_set=gates_set,
                               rng=np.random.default_rng(0))
        seen = set(c.gates)
        assert seen.issubset({GATES.H, GATES.CX})

    def test_from_depth_only_single_qudit_in_gates_set(self):
        """If gates_set contains only single-qudit gates, the produced circuit has only those."""
        gates_set = [GATES.H, GATES.S]
        c = Circuit.from_depth(depth=5, dimensions=[2, 2, 2], gates_set=gates_set,
                               two_qudit_gate_ratio=0.0,
                               rng=np.random.default_rng(0))
        for g in c.gates:
            assert g.n_qudits == 1
            assert g in {GATES.H, GATES.S}

    def test_from_depth_warns_on_multi_qudit(self):
        """Passing a gate acting on >2 qudits emits a warning and ignores it."""
        from sympleq.core.circuits import Gate
        three_qudit_gate = Gate.from_random(n_qudits=3, dimension=2)
        gates_set = [GATES.H, GATES.CX, three_qudit_gate]
        with pytest.warns(UserWarning, match="single qudit and 2-qudits"):
            c = Circuit.from_depth(depth=4, dimensions=[2, 2, 2], gates_set=gates_set,
                                   rng=np.random.default_rng(0))
        for g in c.gates:
            assert g.n_qudits <= 2
            assert g is not three_qudit_gate

    def test_from_depth_invalid_gates_set_type(self):
        """A non-container gates_set should raise ValueError."""
        with pytest.raises(ValueError):
            _ = Circuit.from_depth(depth=3, dimensions=[2, 2], gates_set=GATES.H)  # type: ignore[arg-type]

    def test_from_depth_invalid_ratio(self):
        """two_qudit_gate_ratio outside [0, 1] should raise ValueError."""
        with pytest.raises(ValueError):
            _ = Circuit.from_depth(depth=3, dimensions=[2, 2], two_qudit_gate_ratio=-0.1)
        with pytest.raises(ValueError):
            _ = Circuit.from_depth(depth=3, dimensions=[2, 2], two_qudit_gate_ratio=1.5)

    def test_from_depth_rng_reproducible(self):
        """Same seed yields the same circuit."""
        c1 = Circuit.from_depth(depth=5, dimensions=[2, 2, 2],
                                rng=np.random.default_rng(123))
        c2 = Circuit.from_depth(depth=5, dimensions=[2, 2, 2],
                                rng=np.random.default_rng(123))
        assert c1.n_gates() == c2.n_gates()
        for g1, g2 in zip(c1.gates, c2.gates):
            assert g1 is g2
        assert list(c1.qudit_indices) == list(c2.qudit_indices)


class TestCircuitLayers():

    @staticmethod
    def _as_tuples(circuit: Circuit) -> list[tuple[str, tuple[int, ...]]]:
        return [(g.name, idxs) for g, idxs in zip(circuit.gates, circuit.qudit_indices)]

    def test_layer_basic_partition(self):
        c = Circuit.from_tuples([2, 2, 2], [(GATES.H, 0), (GATES.H, 1), (GATES.CX, 0, 1), (GATES.H, 2)])

        assert c.n_layers() == 2
        assert self._as_tuples(c.layer(0)) == [("H", (0,)), ("H", (1,)), ("H", (2,))]
        assert self._as_tuples(c.layer(1)) == [("CX", (0, 1))]

    def test_layer_single_int_index_returns_circuit(self):
        c = Circuit.from_tuples([2, 2], [(GATES.H, 0), (GATES.S, 1)])

        L = c.layer(0)
        assert isinstance(L, Circuit)
        assert np.array_equal(L.dimensions, c.dimensions)
        assert self._as_tuples(L) == [("H", (0,)), ("S", (1,))]

    def test_layer_multiple_indices_preserve_original_order(self):
        c = Circuit.from_tuples([2, 2, 2], [(GATES.H, 0), (GATES.H, 1), (GATES.CX, 0, 1), (GATES.H, 2)])

        # Passing indices out of order should not change the resulting gate order (original array order).
        combo = c.layer([1, 0])
        assert self._as_tuples(combo) == [("H", (0,)), ("H", (1,)), ("CX", (0, 1)), ("H", (2,))]

    def test_layer_out_of_range_raises(self):
        c = Circuit.from_tuples([2, 2], [(GATES.H, 0), (GATES.CX, 0, 1)])

        with pytest.raises(IndexError):
            c.layer(-1)
        with pytest.raises(IndexError):
            c.layer(c.n_layers())

    def test_layer_preserves_noise_models(self):
        c = Circuit.from_tuples(
            [2, 2], [(GATES.H, 0), (GATES.S, 1)]
        ).with_noise([DepolarizingNoise(0.05), DephasingNoise(0.1)])

        L = c.layer(0)
        assert L.noise_model_per_gate == c.noise_model_per_gate[:L.n_gates()]

    def test_layer_empty_circuit(self):
        c = Circuit.empty([2, 2])

        assert c.n_layers() == 0
        with pytest.raises(IndexError):
            c.layer(0)

    def test_n_layers_matches_depth_for_single_qudit_only(self):
        depth = 4
        c = Circuit.from_depth(depth=depth, dimensions=[2, 2, 2], two_qudit_gate_ratio=0.0,
                               rng=np.random.default_rng(0))

        assert c.n_layers() == depth

    def test_insert_layer_shifts_existing_layers(self):
        # Both qudits are busy in layer 0, so a SWAP touching them cannot merge backwards
        # and genuinely forces a new layer, pushing CX to layer 2.
        c = Circuit.from_tuples([2, 2], [(GATES.H, 0), (GATES.H, 1), (GATES.CX, 0, 1)])
        assert c.n_layers() == 2

        c.insert_layer(1, [GATES.SWAP], [(0, 1)])

        assert c.n_layers() == 3
        assert self._as_tuples(c.layer(0)) == [("H", (0,)), ("H", (1,))]
        assert self._as_tuples(c.layer(1)) == [("SWAP", (0, 1))]
        assert self._as_tuples(c.layer(2)) == [("CX", (0, 1))]

    def test_insert_layer_at_end_equals_append(self):
        c = Circuit.from_tuples([2, 2], [(GATES.H, 0), (GATES.H, 1)])
        n_before = c.n_layers()

        c.insert_layer(n_before, [GATES.S, GATES.S], [(0,), (1,)])

        assert c.n_layers() == n_before + 1
        assert self._as_tuples(c.layer(n_before)) == [("S", (0,)), ("S", (1,))]

    def test_insert_layer_mismatched_lengths_raises(self):
        c = Circuit.from_tuples([2, 2], [(GATES.H, 0)])

        with pytest.raises(ValueError):
            c.insert_layer(0, [GATES.H, GATES.S], [(0,)])

    def test_insert_layer_invalid_position_raises(self):
        c = Circuit.from_tuples([2, 2], [(GATES.H, 0)])

        with pytest.raises(ValueError):
            c.insert_layer(-1, [GATES.S], [(1,)])
        with pytest.raises(ValueError):
            c.insert_layer(c.n_layers() + 1, [GATES.S], [(1,)])

    def test_insert_layer_with_noise_models(self):
        c = Circuit.from_tuples([2, 2], [(GATES.H, 0), (GATES.CX, 0, 1)])

        c.insert_layer(1, [GATES.S], [(1,)], noise_models=[DepolarizingNoise(0.05)])

        assert self._as_tuples(c.layer(0)) == [("H", (0,)), ("S", (1,))]
        s_gate_index = c.gates.index(GATES.S)
        assert isinstance(c.noise_model_per_gate[s_gate_index], DepolarizingNoise)

    def test_layer_partition_reconstructs_full_circuit(self):
        c = Circuit.from_depth(depth=6, dimensions=[2, 2, 2, 2], rng=np.random.default_rng(7))

        reconstructed: list[tuple[str, tuple[int, ...]]] = []
        for i in range(c.n_layers()):
            reconstructed.extend(self._as_tuples(c.layer(i)))

        assert reconstructed == self._as_tuples(c)
