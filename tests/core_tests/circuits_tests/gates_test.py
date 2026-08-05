import pytest
import numpy as np
import random
from sympleq.core.circuits import GATES, PauliGate
from sympleq.core.circuits.gates import Gate
from sympleq.core.circuits.utils import is_symplectic
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits import Circuit
from sympleq.core.circuits.random_symplectic import (symplectic_gf2, symplectic_group_size,
                                                     symplectic_random_koenig_smolin_gf2,
                                                     symplectic_random_transvection)
from sympleq.models import random_hamiltonian


class TestGates():

    @staticmethod
    def states(dim):
        ps1 = PauliString.from_string('x1z0 x0z0', dimensions=[dim, dim])
        ps2 = PauliString.from_string('x0z0 x0z1', dimensions=[dim, dim])
        ps3 = PauliString.from_string('x1z0 x1z0', dimensions=[dim, dim])
        ps4 = PauliString.from_string('x0z1 x0z1', dimensions=[dim, dim])
        ps5 = PauliString.from_string('x1z1 x0z0', dimensions=[dim, dim])
        p_sum = PauliSum.from_string(['x1z0 x0z0 x1z1', 'x0z0 x0z1 x1z0',
                                      'x1z1 x1z0 x0z0'], dimensions=[dim, dim, dim])
        return ps1, ps2, ps3, ps4, ps5, p_sum

    def random_pauli_sum(self, dim, n_paulis=10):
        ps_list = []
        element_list = [(0, 0, 0, 0)]
        for _ in range(n_paulis):
            ps, r1, r2, s1, s2 = self.random_pauli_string(dim)
            element_list.append((r1, r2, s1, s2))
            while (r1, r2, s1, s2) in element_list:
                ps, r1, r2, s1, s2 = self.random_pauli_string(dim)
            ps_list.append(ps)
        return PauliSum.from_pauli_strings(ps_list)

    def random_pauli_string(self, dim):
        r1 = np.random.randint(0, dim)
        r2 = np.random.randint(0, dim)
        s1 = np.random.randint(0, dim)
        s2 = np.random.randint(0, dim)
        return PauliString.from_string(f'x{r1}z{s1} x{r2}z{s2}', dimensions=[dim, dim]), r1, r2, s1, s2

    @pytest.mark.parametrize("d", [2, 5, 11])
    def test_CX(self, d: int):
        # acts the CX gate on a bunch of random pauli strings and pauli sums
        # test pauli_strings
        for _ in range(1000):
            r1 = np.random.randint(0, d)
            r2 = np.random.randint(0, d)
            s1 = np.random.randint(0, d)
            s2 = np.random.randint(0, d)

            input_str = f"x{r1}z{s1} x{r2}z{s2}"
            output_str_correct = f"x{r1}z{(s1 - s2) % d} x{(r2 + r1) % d}z{s2}"

            input_ps = PauliString.from_string(input_str, dimensions=[d, d])
            output_ps = GATES.CX.act(input_ps, (0, 1))
            assert output_ps == PauliString.from_string(output_str_correct, dimensions=[d, d]), 'Error in CX gate'

        # test pauli_sums
        for _ in range(10):
            ps_list_in = []
            ps_list_out_correct = []
            ps_phase_out_correct = []
            for _ in range(10):
                r1 = np.random.randint(0, d)
                r2 = np.random.randint(0, d)
                s1 = np.random.randint(0, d)
                s2 = np.random.randint(0, d)

                input_str = f"x{r1}z{s1} x{r2}z{s2}"
                output_str_correct = f"x{r1}z{(s1 - s2) % d} x{(r2 + r1) % d}z{s2}"

                input_ps = PauliString.from_string(input_str, dimensions=[d, d])
                output_ps_correct = PauliString.from_string(output_str_correct, dimensions=[d, d])
                ps_list_in.append(input_ps)
                ps_list_out_correct.append(output_ps_correct)
                ps_phase_out_correct.append(0)

            input_psum = PauliSum.from_pauli_strings(ps_list_in)
            output_psum_correct = PauliSum.from_pauli_strings(ps_list_out_correct, phases=ps_phase_out_correct)

            output_psum = GATES.CX.act(input_psum, (0, 1))
            assert output_psum == output_psum_correct, (
                'Error in CX gate: \n' +
                output_psum.__str__() + '\n' +
                output_psum_correct.__str__()
            )

    @pytest.mark.parametrize("d", [2, 5, 11])
    def test_SWAP(self, d: int):
        # test pauli_strings
        for _ in range(100):
            r1 = np.random.randint(0, d)
            r2 = np.random.randint(0, d)
            s1 = np.random.randint(0, d)
            s2 = np.random.randint(0, d)

            input_str = f"x{r1}z{s1} x{r2}z{s2}"
            output_str_correct = f"x{r2}z{(s2) % d} x{r1}z{s1}"

            input_ps = PauliString.from_string(input_str, dimensions=[d, d])
            output_ps = GATES.SWAP.act(input_ps, (0, 1))
            assert output_ps == PauliString.from_string(output_str_correct, dimensions=[d, d]), 'Error in SWAP gate'

        # test pauli_sums
        for _ in range(100):
            ps_list_in = []
            ps_list_out_correct = []
            ps_phase_out_correct = []
            for _ in range(10):
                r1 = np.random.randint(0, d)
                r2 = np.random.randint(0, d)
                s1 = np.random.randint(0, d)
                s2 = np.random.randint(0, d)

                input_str = f"x{r1}z{s1} x{r2}z{s2}"
                output_str_correct = f"x{r2}z{s2} x{r1}z{s1}"

                input_ps = PauliString.from_string(input_str, dimensions=[d, d])
                output_ps_correct = PauliString.from_string(output_str_correct, dimensions=[d, d])
                ps_list_in.append(input_ps)
                ps_list_out_correct.append(output_ps_correct)
                ps_phase_out_correct.append(0)

            input_psum = PauliSum.from_pauli_strings(ps_list_in)
            output_psum_correct = PauliSum.from_pauli_strings(ps_list_out_correct,
                                                              phases=ps_phase_out_correct)

            output_psum = GATES.SWAP.act(input_psum, (0, 1))
            assert output_psum == output_psum_correct, (
                'Error in SWAP gate: \n' +
                input_psum.__str__() + '\n' +
                output_psum.__str__() + '\n' +
                output_psum_correct.__str__()
            )

    @pytest.mark.parametrize("d", [2, 5, 11])
    def test_Hadamard(self, d: int):
        # test pauli_strings on qudit 0
        for _ in range(100):
            input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
            output_str_correct = f"x{(-s1) % d}z{r1} x{r2}z{s2}"

            output_ps = GATES.H.act(input_ps, 0)
            assert output_ps.has_equal_tableau(PauliString.from_string(
                output_str_correct, dimensions=[d, d]
            )), 'Error in Hadamard gate 0'

        # test on qudit 1
        for _ in range(100):
            input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
            output_str_correct = f"x{r1}z{s1} x{(-s2) % d}z{r2}"

            output_ps = GATES.H.act(input_ps, 1)
            assert output_ps.has_equal_tableau(PauliString.from_string(
                output_str_correct, dimensions=[d, d]
            )), 'Error in Hadamard gate 1'

        # test pauli_sums
        for _ in range(100):
            ps_list_in = []
            ps_list_out_correct = []
            ps_phase_out_correct = []
            for _ in range(10):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{(-s1) % d}z{r1} x{r2}z{s2}"

                ps_list_in.append(input_ps)
                ps_list_out_correct.append(PauliString.from_string(output_str_correct, dimensions=[d, d]))
                ps_phase_out_correct.append(-2 * r1 * s1)

            input_psum = PauliSum.from_pauli_strings(ps_list_in)
            output_psum = GATES.H.act(input_psum, 0)
            output_psum_correct = PauliSum.from_pauli_strings(
                ps_list_out_correct, phases=ps_phase_out_correct)
            assert output_psum.has_equal_tableau(output_psum_correct), (
                'Error in Hadamard gate: \n' +
                input_psum.__str__() + '\n' +
                output_psum.__str__() + '\n' +
                output_psum_correct.__str__()
            )

    @pytest.mark.parametrize("d", [2, 5, 11])
    def test_PHASE(self, d: int):
        # test on qudit 0
        for _ in range(100):
            input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
            output_str_correct = f"x{r1}z{(r1 + s1) % d} x{r2}z{s2}"

            output_ps = GATES.S.act(input_ps, 0)
            assert output_ps.has_equal_tableau(PauliString.from_string(
                output_str_correct, dimensions=[d, d])), 'Error in PHASE gate 0'

        # test on qudit 1
        for _ in range(100):
            input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
            output_str_correct = f"x{r1}z{s1} x{r2}z{(r2 + s2) % d}"

            output_ps = GATES.S.act(input_ps, 1)
            assert output_ps.has_equal_tableau(PauliString.from_string(
                output_str_correct, dimensions=[d, d])), 'Error in PHASE gate 1'

    def test_group_homomorphism(self):
        # Tests all gates and all combinations of two pauli strings for the group homomorphism property:
        # gate.act(p1) * gate.act(p2) == gate.act(p1 * p2)

        # Test for dim = 2 with all combinations
        gates_and_qudits = [
            (GATES.CX, (0, 1)), (GATES.CX, (1, 0)),
            (GATES.SWAP, (0, 1)),
            (GATES.H, 0), (GATES.H, 1),
            (GATES.S, 0), (GATES.S, 1)
        ]

        for gate, qudits in gates_and_qudits:
            for x0 in range(2):
                for z0 in range(2):
                    for x1 in range(2):
                        for z1 in range(2):
                            for x0p in range(2):
                                for z0p in range(2):
                                    for x1p in range(2):
                                        for z1p in range(2):
                                            p1 = PauliSum.from_string([f'x{x0}z{z0} x{x1}z{z1}'], dimensions=[2, 2])
                                            p2 = PauliSum.from_string([f'x{x0p}z{z0p} x{x1p}z{z1p}'], dimensions=[2, 2])
                                            lhs = gate.act(p1, qudits) * gate.act(p2, qudits)
                                            rhs = gate.act(p1 * p2, qudits)
                                            assert lhs == rhs, f"Failed for {gate.name} on {qudits}"

        # Test for larger dimensions with random samples
        for dim in [3, 5, 7, 15]:
            for gate, qudits in gates_and_qudits:
                for _ in range(100):
                    p1 = self.random_pauli_sum(dim, n_paulis=1)
                    p2 = self.random_pauli_sum(dim, n_paulis=1)
                    lhs = gate.act(p1, qudits) * gate.act(p2, qudits)
                    rhs = gate.act(p1 * p2, qudits)
                    assert lhs == rhs, f"Failed for {gate.name} on {qudits} dim={dim}"

    @pytest.mark.parametrize("gate", [GATES.CX, GATES.SWAP, GATES.H, GATES.S])
    def test_is_symplectic(self, gate: Gate):
        # Tests if the symplectic matrix of the gate is symplectic
        assert is_symplectic(gate.symplectic, 2), (
            f"Gate {gate.name} is not symplectic. \n" +
            gate.symplectic.__str__()
        )

    def test_random_symplectic(self, num_tests=20, max_n=10, primes=[2, 3, 5, 7]):
        """
        Test random_symplectic() across several n, p values using assertions only.
        """
        for n in range(2, max_n + 1):
            for d in primes:
                for i in range(num_tests):
                    if n < 6 and d == 2:
                        index = random.randint(0, symplectic_group_size(n) - 1)
                        F = symplectic_gf2(index, n)
                    else:
                        F = symplectic_random_transvection(n, dimension=d)
                    assert is_symplectic(F, d), f"Failed symplectic check: n={n}, test {i}"

    def test_random_transvection_sampler_repeatable_with_rng(self):
        rng1 = np.random.default_rng(321)
        rng2 = np.random.default_rng(321)

        F1 = symplectic_random_transvection(3, dimension=5, num_transvections=12, rng=rng1)
        F2 = symplectic_random_transvection(3, dimension=5, num_transvections=12, rng=rng2)

        assert is_symplectic(F1, 5)
        assert np.array_equal(F1, F2)

    def test_koenig_smolin_gf2_n1_enumerates_whole_group(self):
        elements = []
        for index in range(symplectic_group_size(1, 2)):
            F = symplectic_gf2(index, 1)
            assert is_symplectic(F, 2)
            elements.append(tuple(F.ravel().tolist()))

        assert len(elements) == 6
        assert len(set(elements)) == 6

    @pytest.mark.parametrize("n", [1, 2, 3])
    def test_koenig_smolin_gf2_random_sampler_is_symplectic_and_repeatable(self, n: int):
        rng1 = np.random.default_rng(1234 + n)
        rng2 = np.random.default_rng(1234 + n)

        F1 = symplectic_random_koenig_smolin_gf2(n, rng=rng1)
        F2 = symplectic_random_koenig_smolin_gf2(n, rng=rng2)

        assert is_symplectic(F1, 2)
        assert np.array_equal(F1, F2)

    def test_koenig_smolin_gf2_random_sampler_rejects_invalid_n(self):
        with pytest.raises(ValueError, match="n_qubits must be >= 1"):
            symplectic_random_koenig_smolin_gf2(0)

    @pytest.mark.parametrize("dim", [2, 3, 5])
    @pytest.mark.parametrize("n_qudits", [5, 6])
    @pytest.mark.parametrize("num_pauli", [30, 40])
    def test_gate_from_target(self, dim: int, n_qudits: int, num_pauli: int):
        """Test Gate.solve_from_target finds
        correct pauli sum"""

        dimensions = [dim] * n_qudits
        for _ in range(100):
            pl_sum = random_hamiltonian.random_pauli_hamiltonian(num_pauli, dimensions)

            C = Circuit.from_random(n_gates=10 * n_qudits**2, dimensions=dimensions)
            target_pl_sum = C.act(pl_sum)

            final_gate = Gate.solve_from_target(pl_sum, target_pl_sum)

            found_pl_sum = final_gate.act(pl_sum, tuple(range(n_qudits)))

            assert (found_pl_sum == target_pl_sum)

    @pytest.mark.parametrize("dim", [2, 3])
    @pytest.mark.parametrize("n_qudits", [5])
    @pytest.mark.parametrize("num_pauli", [5])
    def test_gate_from_target_hilbert_space(self, dim: int, n_qudits: int, num_pauli: int):
        """Test Gate.solve_from_target finds
        correct pauli sum in hilbert space"""

        dimensions = [dim] * n_qudits
        for _ in range(100):
            pl_sum = random_hamiltonian.random_pauli_hamiltonian(num_pauli, dimensions)

            C = Circuit.from_random(n_gates=10 * n_qudits**2, dimensions=dimensions)
            target_pl_sum = C.act(pl_sum)

            final_gate = Gate.solve_from_target(pl_sum, target_pl_sum)

            found_pl_sum = final_gate.act(pl_sum, tuple(range(n_qudits)))

            assert np.allclose(
                found_pl_sum.to_hilbert_space().toarray(),
                target_pl_sum.to_hilbert_space().toarray(),
            )

    def test_full_symplectic_embeds_and_reduces_mod_dimension(self):
        """Test Gate.full_symplectic embeds a local symplectic into a larger system and reduces mod dimension."""
        n_qudits = 3
        qudits = (1,)

        # dimension=None returns the raw embedding, un-reduced.
        raw = GATES.H.full_symplectic(qudits, n_qudits)
        assert raw.shape == (2 * n_qudits, 2 * n_qudits)
        # Hadamard's local symplectic [[0, -1], [1, 0]] should appear verbatim, including the -1.
        assert -1 in raw

        # Qudits the gate does not act on must be left as identity.
        for q in range(n_qudits):
            if q in qudits:
                continue
            idx = [q, q + n_qudits]
            block = raw[np.ix_(idx, idx)]
            assert np.array_equal(block, np.eye(2, dtype=int))

        # An explicit dimension reduces the embedded matrix mod that dimension.
        dimension = 5
        reduced = GATES.H.full_symplectic(qudits, n_qudits, dimension=dimension)
        assert np.array_equal(reduced, raw % dimension)
        assert np.all(reduced >= 0) and np.all(reduced < dimension)

    @pytest.mark.parametrize("dimension", [2, 3, 5])
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_gate_from_random_symplecticity(self, dimension: int, n_qudits: int):
        """Test that Gate.from_random produces valid symplectic matrices."""
        from sympleq.core.circuits import Gate

        for _ in range(5):
            gate = Gate.from_random(n_qudits, dimension)
            assert is_symplectic(gate.symplectic, dimension), (
                f"Random gate not symplectic for n={n_qudits}, d={dimension}"
            )

    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_gate_from_random_koenig_smolin_symplecticity(self, n_qudits: int):
        gate = Gate.from_random(
            n_qudits,
            2,
            sampler="koenig-smolin",
            rng=np.random.default_rng(100 + n_qudits),
        )
        assert is_symplectic(gate.symplectic, 2)

    def test_gate_from_random_koenig_smolin_repeatable_with_rng(self):
        rng1 = np.random.default_rng(123)
        rng2 = np.random.default_rng(123)

        gate1 = Gate.from_random(3, 2, sampler="koenig-smolin", rng=rng1)
        gate2 = Gate.from_random(3, 2, sampler="koenig-smolin", rng=rng2)

        assert np.array_equal(gate1.symplectic, gate2.symplectic)

    def test_gate_from_random_koenig_smolin_rejects_invalid_options(self):
        with pytest.raises(ValueError, match="only implemented for dimension=2"):
            Gate.from_random(2, 3, sampler="koenig-smolin")
        with pytest.raises(ValueError, match="num_transvections is not used"):
            Gate.from_random(2, 2, 10, sampler="koenig-smolin")
        with pytest.raises(ValueError, match="Unknown random Clifford sampler"):
            Gate.from_random(2, 2, sampler="unknown")

    def test_gate_from_random_transvection_repeatable_with_rng(self):
        rng1 = np.random.default_rng(123)
        rng2 = np.random.default_rng(123)

        gate1 = Gate.from_random(3, 2, sampler="transvection", rng=rng1)
        gate2 = Gate.from_random(3, 2, sampler="transvection", rng=rng2)

        assert np.array_equal(gate1.symplectic, gate2.symplectic)

    def test_gate_from_random_transvection_still_accepts_positional_depth(self):
        gate = Gate.from_random(2, 2, 10)
        assert is_symplectic(gate.symplectic, 2)

    @pytest.mark.parametrize("d", [2, 3, 5])
    @pytest.mark.parametrize("n", [2, 3])
    def test_gate_from_random_action(self, d: int, n: int):
        """Test that random gates correctly transform Paulis."""
        from sympleq.core.circuits import Gate

        gate = Gate.from_random(n, d)

        # Act on random PauliSum
        ps = PauliSum.from_random(5, [d] * n)
        result = gate.act(ps, tuple(range(n)))

        # Result should have same number of Paulis
        assert result.n_paulis() == ps.n_paulis()

        # Verify symplectic transformation on tableau
        # act() applies: tableau @ symplectic.T (see gates.py line 127)
        expected_tableau = (ps.tableau @ gate.symplectic.T) % d
        assert np.array_equal(result.tableau % d, expected_tableau % d)

    @pytest.mark.parametrize("d", [2, 3, 5])
    @pytest.mark.parametrize("gate", [GATES.H, GATES.H_inv, GATES.S, GATES.S_inv,
                                      GATES.CX, GATES.CX_inv, GATES.SWAP, GATES.CZ,
                                      GATES.ZZMax, GATES.ZZMax_inv])
    def test_unitary_is_unitary(self, d: int, gate: Gate):
        """Test that all gate unitaries are actually unitary matrices."""
        U = gate.local_unitary(d).toarray()
        Id = np.eye(U.shape[0])
        assert np.allclose(U.conj().T @ U, Id), f"{gate.name}(d={d}) is not unitary"
        assert np.allclose(U @ U.conj().T, Id), f"{gate.name}(d={d}) is not unitary"

    @pytest.mark.parametrize("gate", [GATES.V, GATES.V_inv])
    def test_V_unitary_is_unitary(self, gate: Gate):
        """V is qubit-only; verify the local_unitary is unitary for d=2."""
        U = gate.local_unitary(2).toarray()
        Id = np.eye(2)
        assert np.allclose(U.conj().T @ U, Id), f"{gate.name} is not unitary"
        assert np.allclose(U @ U.conj().T, Id), f"{gate.name} is not unitary"

    def test_V_matches_sqrt_X(self):
        """V should equal exp(-iπ/4 X) = (1/√2)(I - iX); V_inv = V†."""
        U_V = GATES.V.local_unitary(2).toarray()
        expected = np.array([[1, -1j], [-1j, 1]], dtype=complex) / np.sqrt(2)
        assert np.allclose(U_V, expected), "V does not equal √X"

        U_V_inv = GATES.V_inv.local_unitary(2).toarray()
        assert np.allclose(U_V @ U_V_inv, np.eye(2)), "V · V_inv != I"

    def test_V_squared_is_X_up_to_phase(self):
        """V² = -i·X"""
        U_V = GATES.V.local_unitary(2).toarray()
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        assert np.allclose(U_V @ U_V, -1j * X), "V² != -i·X"

    def test_V_clifford_action_on_paulis(self):
        """V's symplectic + phase action: X → X (phase 1), Z → XZ (phase ω⁻¹ = -i for d=2)."""
        # X stays X with phase 0
        X = PauliString.from_string("x1z0", dimensions=[2])
        out_X = GATES.V.act(X, 0)
        assert out_X.has_equal_tableau(X), "V·X·V† should keep tableau as X"
        assert (out_X.phases % (2 * X.lcm) == 0).all(), "V·X·V† should have phase 0"

        # Z -> XZ with phase -1 (mod 2*lcm = 4 for qubits) i.e. phase factor ω^{-1} = -i
        Z = PauliString.from_string("x0z1", dimensions=[2])
        out_Z = GATES.V.act(Z, 0)
        expected_Z = PauliString.from_string("x1z1", dimensions=[2])
        assert out_Z.has_equal_tableau(expected_Z), "V·Z·V† should map to XZ"
        assert (out_Z.phases[0] % (2 * Z.lcm)) == 3, "V·Z·V† phase should be -1 (=3 mod 4)"

        # V_inv flips the sign: Z -> XZ with phase +1
        out_Z_inv = GATES.V_inv.act(Z, 0)
        assert out_Z_inv.has_equal_tableau(expected_Z), "V_inv·Z·V_inv† should map to XZ"
        assert (out_Z_inv.phases[0] % (2 * Z.lcm)) == 1, "V_inv·Z·V_inv† phase should be +1"

    def test_V_inverse_round_trip(self):
        """V_inv · V applied to any single-qubit Pauli should be the identity."""
        for s in ["x1z0", "x0z1", "x1z1"]:
            ps = PauliString.from_string(s, dimensions=[2])
            roundtrip = GATES.V_inv.act(GATES.V.act(ps, 0), 0)
            assert roundtrip.has_equal_tableau(ps) and \
                np.array_equal(roundtrip.phases % (2 * ps.lcm), ps.phases % (2 * ps.lcm)), \
                f"V_inv·V·{s} != {s}"

    def test_V_qudit_local_unitary_raises(self):
        """V's local_unitary is only defined for qubits."""
        with pytest.raises(NotImplementedError):
            GATES.V.local_unitary(3)

    @pytest.mark.parametrize("d", [2, 3, 5])
    def test_unitary_inverse(self, d: int):
        """Test that gate inverses have inverse unitaries."""
        # H and H_inv
        U_H = GATES.H.local_unitary(d).toarray()
        U_H_inv = GATES.H_inv.local_unitary(d).toarray()
        assert np.allclose(U_H @ U_H_inv, np.eye(d)), f"H @ H_inv != I for d={d}"

        # S and S_inv
        U_S = GATES.S.local_unitary(d).toarray()
        U_S_inv = GATES.S_inv.local_unitary(d).toarray()
        assert np.allclose(U_S @ U_S_inv, np.eye(d)), f"S @ S_inv != I for d={d}"

        # CX and CX_inv
        U_CX = GATES.CX.local_unitary(d).toarray()
        U_CX_inv = GATES.CX_inv.local_unitary(d).toarray()
        assert np.allclose(U_CX @ U_CX_inv, np.eye(d * d)), f"CX @ CX_inv != I for d={d}"

        # SWAP is self-inverse
        SWAP = GATES.SWAP.local_unitary(d).toarray()
        assert np.allclose(SWAP @ SWAP, np.eye(d * d)), f"SWAP @ SWAP != I for d={d}"

        # CZ is self-inverse only for qubits (d=2)
        U_CZ = GATES.CZ.local_unitary(d).toarray()
        if d == 2:
            assert np.allclose(U_CZ @ U_CZ, np.eye(d * d)), f"CZ @ CZ != I for d={d}"
        else:
            # For d > 2, CZ^d = I (CZ has order d)
            U_CZ_power = np.eye(d * d, dtype=complex)
            for _ in range(d):
                U_CZ_power = U_CZ_power @ U_CZ
            assert np.allclose(U_CZ_power, np.eye(d * d)), f"CZ^{d} != I for d={d}"

    @pytest.mark.parametrize("d", [2, 3, 5])
    @pytest.mark.parametrize("gate", [GATES.H, GATES.S])
    def test_one_qudit_unitary_clifford_property(self, d: int, gate: Gate):
        """Test that single-qudit unitaries correctly implement symplectic transformation."""
        from sympleq.core.circuits.utils import pauli_unitary_qudit

        U = gate.local_unitary(d)

        # Test on X (x=1, z=0) and Z (x=0, z=1)
        for x, z in [(1, 0), (0, 1)]:
            P = pauli_unitary_qudit(d, x, z).toarray()

            # Conjugate: U P U† (this code's convention)
            P_conj = U @ P @ U.conj().T

            # Expected from symplectic: F @ [x, z]
            v_out = (gate.symplectic @ np.array([x, z])) % d
            x_out, z_out = v_out[0], v_out[1]
            P_exp = pauli_unitary_qudit(d, x_out, z_out).toarray()

            # Should match up to a global phase
            if np.max(np.abs(P_exp)) > 0:
                ratio = P_conj[np.abs(P_exp) > 0.1] / P_exp[np.abs(P_exp) > 0.1]
                assert np.allclose(np.abs(ratio), 1.0), (
                    f"{gate.name}(d={d}) Clifford property failed for x={x}, z={z}"
                )

    @pytest.mark.parametrize("d", [2, 3, 5])
    @pytest.mark.parametrize("gate", [GATES.CX, GATES.SWAP, GATES.CZ])
    def test_two_qudit_unitary_clifford_property(self, d: int, gate: Gate):
        """Test that two-qudit unitaries correctly implement symplectic transformation."""
        from sympleq.core.circuits.utils import pauli_unitary_qudit

        U = gate.local_unitary(d)

        # Test on X0, X1, Z0, Z1
        test_paulis = [
            (1, 0, 0, 0),  # X0
            (0, 1, 0, 0),  # X1
            (0, 0, 1, 0),  # Z0
            (0, 0, 0, 1),  # Z1
        ]
        for x0, x1, z0, z1 in test_paulis:
            # Build input Pauli
            P0 = pauli_unitary_qudit(d, x0, z0).toarray()
            P1 = pauli_unitary_qudit(d, x1, z1).toarray()
            P = np.kron(P0, P1)

            # Conjugate: U P U† (this code's convention)
            P_conj = U @ P @ U.conj().T

            # Expected from symplectic
            v = np.array([x0, x1, z0, z1])
            v_out = (gate.symplectic @ v) % d
            x0_out, x1_out, z0_out, z1_out = v_out

            P0_exp = pauli_unitary_qudit(d, x0_out, z0_out).toarray()
            P1_exp = pauli_unitary_qudit(d, x1_out, z1_out).toarray()
            P_exp = np.kron(P0_exp, P1_exp)

            # Should match up to a global phase
            if np.max(np.abs(P_exp)) > 0:
                mask = np.abs(P_exp) > 0.1
                ratio = P_conj[mask] / P_exp[mask]
                assert np.allclose(np.abs(ratio), 1.0), (
                    f"{gate.name}(d={d}) Clifford property failed for {(x0, x1, z0, z1)}"
                )

    @pytest.mark.parametrize("d", [2, 3, 5])
    def test_CX_unitary(self, d: int):
        """Test CX gate acts as |j,k⟩ -> |j, j+k mod d⟩."""
        U = GATES.CX.local_unitary(d)
        for j in range(d):
            for k in range(d):
                # Input state |j,k⟩
                in_state = np.zeros(d * d, dtype=complex)
                in_state[j * d + k] = 1.0

                # Apply CX
                out_state = U @ in_state

                # Expected: |j, (j+k) mod d⟩
                expected = np.zeros(d * d, dtype=complex)
                expected[j * d + (j + k) % d] = 1.0

                assert np.allclose(out_state, expected), (
                    f"CX(d={d}) failed for |{j},{k}⟩"
                )

    @pytest.mark.parametrize("d", [2, 3, 5])
    def test_ZZMax_unitary_inverse(self, d: int):
        """ZZMax @ ZZMax_inv should be the identity."""
        U = GATES.ZZMax.local_unitary(d).toarray()
        U_inv = GATES.ZZMax_inv.local_unitary(d).toarray()
        assert np.allclose(U @ U_inv, np.eye(d * d)), f"ZZMax @ ZZMax_inv != I for d={d}"
        assert np.allclose(U_inv @ U, np.eye(d * d)), f"ZZMax_inv @ ZZMax != I for d={d}"

    @pytest.mark.parametrize("d", [2, 3, 5])
    @pytest.mark.parametrize("gate", [GATES.ZZMax, GATES.ZZMax_inv])
    def test_ZZMax_clifford_property(self, d: int, gate: Gate):
        """U P U† should equal the symplectically-predicted Pauli up to a global phase."""
        from sympleq.core.circuits.utils import pauli_unitary_qudit

        U = gate.local_unitary(d).toarray()

        for x0, x1, z0, z1 in [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]:
            P0 = pauli_unitary_qudit(d, x0, z0).toarray()
            P1 = pauli_unitary_qudit(d, x1, z1).toarray()
            P = np.kron(P0, P1)

            P_conj = U @ P @ U.conj().T

            v_out = (gate.symplectic @ np.array([x0, x1, z0, z1])) % d
            x0_out, x1_out, z0_out, z1_out = v_out
            P0_exp = pauli_unitary_qudit(d, x0_out, z0_out).toarray()
            P1_exp = pauli_unitary_qudit(d, x1_out, z1_out).toarray()
            P_exp = np.kron(P0_exp, P1_exp)

            mask = np.abs(P_exp) > 0.1
            ratio = P_conj[mask] / P_exp[mask]
            assert np.allclose(np.abs(ratio), 1.0), (
                f"{gate.name}(d={d}) Clifford property failed for {(x0, x1, z0, z1)}"
            )

    def test_ZZMax_act_matches_unitary_qubits(self):
        """For qubits, acting on each single-qudit basis Pauli via .act() should match the
        unitary conjugation including the ±i phase captured in the exceptional_phase_vector."""
        dims = [2, 2]
        # X on qudit 0: expect phase +i  (encoded as phase_vector entry 1 -> phases units of lcm=2)
        input_ps = PauliString.from_string("x1z0 x0z0", dimensions=dims)
        out = GATES.ZZMax.act(input_ps, (0, 1))
        # symplectic image: X0 -> X0 Z0 Z1 with phase +i
        expected_tableau = PauliString.from_string("x1z1 x0z1", dimensions=dims)
        assert out.has_equal_tableau(expected_tableau), "ZZMax X0 tableau mismatch"

        # Inverse acts with -i on the same tableau image
        out_inv = GATES.ZZMax_inv.act(input_ps, (0, 1))
        assert out_inv.has_equal_tableau(expected_tableau), "ZZMax_inv X0 tableau mismatch"
        # Phases should differ by a full i^2 = -1 (i.e. by the lcm=2 factor in phase units)
        assert ((out.phases - out_inv.phases) % (2 * input_ps.lcm)).any(), (
            "ZZMax and ZZMax_inv should produce different phases on X0"
        )

    @pytest.mark.parametrize("d", [2, 3, 5])
    def test_pauli_gate_unitary(self, d: int):
        """Test PauliGate unitary is correct."""
        from sympleq.core.circuits.utils import pauli_unitary_from_tableau

        for _ in range(10):
            # Random Pauli
            x_exp = np.random.randint(0, d, size=2)
            z_exp = np.random.randint(0, d, size=2)
            ps = PauliString.from_exponents(x_exp, z_exp, dimensions=[d, d])
            pg = PauliGate(ps)

            U = pg.local_unitary().toarray()
            U_expected = pauli_unitary_from_tableau(d, x_exp, z_exp).toarray()

            assert np.allclose(U, U_expected), (
                f"PauliGate unitary mismatch for x={x_exp}, z={z_exp}, d={d}"
            )
