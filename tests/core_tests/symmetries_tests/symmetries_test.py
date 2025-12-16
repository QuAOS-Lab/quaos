from sympleq.models.random_hamiltonian import random_gate_symmetric_hamiltonian
from sympleq.core.circuits import SWAP
from sympleq.core.symmetries.clifford import find_clifford_symmetries, qudit_cost, min_qudit_clifford_symmetry
from sympleq.core.circuits import Circuit
from sympleq.core.symmetries.block_decomposition import (
    symplectic_basis_from_span,
    independent_columns,
    complete_symplectic_local,
)
from sympleq.core.symmetries.modular_helpers import omega_matrix, mod_p, inv_mod_scalar
import numpy as np


class TestSymmetryFinder:
    def test_symplectic_basis_from_span_handles_overlapping_pairs(self):
        """
        Regression: spans where chosen u vectors overlap previous pairs used to fail the
        symplectic completion, raising RuntimeError. Ensure we now build a canonical basis.
        """
        p = 2
        B = np.array(
            [
                [1, 0, 1, 0],
                [1, 0, 0, 1],
                [1, 1, 0, 1],
                [1, 0, 0, 0],
            ],
            dtype=np.int64,
        )
        S = independent_columns(B, p)
        Omega = omega_matrix(S.shape[0] // 2, p)

        T_new = symplectic_basis_from_span(S, p)
        G_new = mod_p(T_new.T @ Omega @ T_new, p)
        assert np.array_equal(G_new % p, omega_matrix(T_new.shape[1] // 2, p))

    def test_symplectic_basis_from_span_random_nondegenerate(self):
        """Random full-dimension spans produce a symplectic basis."""
        p = 2
        rng = np.random.default_rng(0)
        for n_modes in [1, 2, 3]:
            n2 = 2 * n_modes
            Omega = omega_matrix(n_modes, p)
            for _ in range(10):
                # Build a full-rank span; retry until the symplectic form is non-degenerate on it
                for _ in range(100):
                    B = independent_columns(rng.integers(0, p, size=(n2, n2), dtype=np.int64), p)
                    if B.shape[1] != n2:
                        continue
                    if np.linalg.matrix_rank(mod_p(B.T @ Omega @ B, p)) == n2:
                        break
                else:
                    raise AssertionError("Failed to sample a non-degenerate span")

                T = symplectic_basis_from_span(B, p)
                G = mod_p(T.T @ Omega @ T, p)
                assert np.array_equal(G % p, omega_matrix(T.shape[1] // 2, p))

    def test_complete_symplectic_local_extends_basis(self):
        """A non-degenerate subspace can be completed to a full symplectic basis."""
        p = 2
        rng = np.random.default_rng(1)
        n_modes = 3  # ambient dimension 6
        n2 = 2 * n_modes
        Omega = omega_matrix(n_modes, p)

        # Sample a 4D non-degenerate subspace inside the 6D space
        for _ in range(200):
            B = independent_columns(rng.integers(0, p, size=(n2, n2), dtype=np.int64), p)
            if B.shape[1] < 4:
                continue
            B = B[:, :4]  # enforce even dimension 4
            if np.linalg.matrix_rank(mod_p(B.T @ Omega @ B, p)) == 4:
                break
        else:
            raise AssertionError("Failed to sample a non-degenerate 4D subspace")

        T_blk = symplectic_basis_from_span(B, p)
        T_local = complete_symplectic_local(T_blk, p)

        # First k columns should preserve the original block's U part
        k = T_blk.shape[1] // 2
        assert np.array_equal(T_local[:, :k] % p, T_blk[:, :k] % p)

        G_local = mod_p(T_local.T @ Omega @ T_local, p)
        assert np.array_equal(G_local % p, omega_matrix(n_modes, p))

    def test_random_SWAP_symmetry(self):
        n_tests = 30
        p = 2
        n_qudits = 10
        n_paulis = 58
        for _ in range(n_tests):
            sym = SWAP(0, 1, p)
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(sym, n_qudits, n_paulis, scrambled=False)
            C = Circuit.from_random(100, H.dimensions).composite_gate()  # scrambling circuit
            H = C.act(H)
            H.weight_to_phase()
            scrambled_sym = Circuit(H.dimensions, [C.inv(), sym, C]).composite_gate()
            assert H.to_standard_form() == scrambled_sym.act(H).to_standard_form(
            ), f"\n{H.to_standard_form().__str__()}\n{sym.act(H).to_standard_form().__str__()}"

            known_F = scrambled_sym.symplectic
            symmetries = find_clifford_symmetries(H)

            assert len(symmetries) != 0

            for c in symmetries:
                print(np.all(c.symplectic == known_F) and np.all(
                    c.phase_vector == scrambled_sym.phase_vector))
                H_s = H.to_standard_form()
                H_out = c.act(H).to_standard_form()
                H_s.weight_to_phase()
                H_out.weight_to_phase()
                assert np.all(H_s.tableau == H_out.tableau)
                assert np.all(H_s.phases == H_out.phases)
                assert np.all(H_s.weights == H_out.weights)

                assert c.act(H).to_standard_form() == H.to_standard_form()

    def test_random_SWAP_symmetry_with_block_decomposition(self):
        n_tests = 30
        p = 2
        n_qudits = 10
        n_paulis = 58
        for _ in range(n_tests):
            sym = SWAP(0, 1, p)
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(sym, n_qudits, n_paulis, scrambled=False)
            C = Circuit.from_random(100, H.dimensions).composite_gate()  # scrambling circuit
            H = C.act(H)
            H.weight_to_phase()
            scrambled_sym = Circuit(H.dimensions, [C.inv(), sym, C]).composite_gate()
            assert H.to_standard_form() == scrambled_sym.act(H).to_standard_form(
            ), f"\n{H.to_standard_form().__str__()}\n{sym.act(H).to_standard_form().__str__()}"

            F, S, T = min_qudit_clifford_symmetry(H, check_symmetry=False)

            assert np.all(F.symplectic == scrambled_sym.symplectic)
            assert np.all(F.phase_vector == scrambled_sym.phase_vector)
            assert F == Circuit(F.dimensions, [T.inv(), S, T]).composite_gate()

            assert H.to_standard_form() == F.act(H).to_standard_form()
            assert T.act(S.act(T.inv().act(H))).to_standard_form() == H.to_standard_form()

            assert S.act(T.inv().act(H)).to_standard_form() == T.inv().act(H).to_standard_form()
            assert qudit_cost(S) == 2

    def test_random_multi_SWAP_symmetry(self):

        n_tests = 100
        p = 2
        n_qudits = 3
        n_paulis = 7
        for _ in range(n_tests):
            sym = Circuit([p] * n_qudits, [SWAP(0, 1, p), SWAP(1, 2, p)])  #
            sym = sym.composite_gate()
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(sym, n_qudits, n_paulis, scrambled=False)
            C = Circuit.from_random(100, H.dimensions).composite_gate()  # scrambling circuit
            H = C.act(H)
            H.weight_to_phase()
            scrambled_sym = Circuit(H.dimensions, [C.inv(), sym, C]).composite_gate()
            assert H.to_standard_form() == scrambled_sym.act(H).to_standard_form(
            ), f"\n{H.to_standard_form().__str__()}\n{sym.act(H).to_standard_form().__str__()}"

            known_F = scrambled_sym.symplectic
            circ = find_clifford_symmetries(H)

            assert len(circ) != 0

            for c in circ:
                print(np.all(c.symplectic == known_F) and np.all(
                    c.phase_vector == scrambled_sym.phase_vector))
                H_s = H.to_standard_form()
                H_out = c.act(H).to_standard_form()
                H_s.weight_to_phase()
                H_out.weight_to_phase()
                assert np.all(H_s.tableau == H_out.tableau)
                assert np.all(H_s.phases == H_out.phases)
                assert np.all(H_s.weights == H_out.weights)

                assert c.act(H).to_standard_form() == H.to_standard_form()

    def test_random_multi_SWAP_symmetry_with_block_decomposition(self):

        n_tests = 100
        p = 2
        n_qudits = 6
        n_paulis = 15
        for _ in range(n_tests):
            sym = Circuit([p] * n_qudits, [SWAP(0, 1, p), SWAP(1, 2, p)])  #
            sym = sym.composite_gate()

            H = random_gate_symmetric_hamiltonian(sym, n_qudits, n_paulis, scrambled=False)
            C = Circuit.from_random(100, H.dimensions).composite_gate()  # scrambling circuit
            H = C.act(H)
            H.weight_to_phase()
            scrambled_sym = Circuit(H.dimensions, [C.inv(), sym, C]).composite_gate()
            assert H.to_standard_form() == scrambled_sym.act(H).to_standard_form(
            ), f"\n{H.to_standard_form().__str__()}\n{sym.act(H).to_standard_form().__str__()}"

            F, S, T = min_qudit_clifford_symmetry(H, check_symmetry=False)

            # there may be multiple expressions of the symmetry so these are too harsh
            # assert np.all(F.symplectic == scrambled_sym.symplectic)
            # assert np.all(F.phase_vector == scrambled_sym.phase_vector)

            assert F == Circuit(F.dimensions, [T.inv(), S, T]).composite_gate()

            assert H.to_standard_form() == F.act(H).to_standard_form()
            assert T.act(S.act(T.inv().act(H))).to_standard_form() == H.to_standard_form()

            assert S.act(T.inv().act(H)).to_standard_form() == T.inv().act(H).to_standard_form()
            assert qudit_cost(S) <= 3

    # def random_Hadamard_symmetry(self):
        # pass

    # def test_random_multi_gate_symmetry(self):
    #     """ Selects a random set of gates and injects that as a symmetry."""
    #     pass
