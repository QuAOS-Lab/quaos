from sympleq.models.random_hamiltonian import random_gate_symmetric_hamiltonian
from sympleq.core.circuits import SWAP
from sympleq.core.symmetries.clifford import find_clifford_symmetries, qudit_cost, min_qudit_clifford_symmetry
from sympleq.core.circuits import Circuit
import numpy as np


class TestSymmetryFinder:

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

            F, S, T = min_qudit_clifford_symmetry(H)

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

            F, S, T = min_qudit_clifford_symmetry(H)

            # there may be multiple expressions of the symmetry so these are too harsh
            # assert np.all(F.symplectic == scrambled_sym.symplectic)
            # assert np.all(F.phase_vector == scrambled_sym.phase_vector)

            assert F == Circuit(F.dimensions, [T.inv(), S, T]).composite_gate()

            assert H.to_standard_form() == F.act(H).to_standard_form()
            assert T.act(S.act(T.inv().act(H))).is_close(H, literal=False)

            assert S.act(T.inv().act(H)).to_standard_form() == T.inv().act(H).to_standard_form()
            assert qudit_cost(S) <= 3

    def test_random_arbitrary_symmetry(self):

        n_tests = 2000
        p = 2
        n_qudits = 3
        n_paulis = 7

        for _ in range(n_tests):
            sym = Circuit.from_random(10, [p] * n_qudits)  #
            sym = sym.composite_gate()
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(sym, n_qudits, n_paulis, scrambled=False)
            C = Circuit.from_random(100, H.dimensions).composite_gate()  # scrambling circuit
            H = C.act(H)
            H.weight_to_phase()
            H.weights = np.round(H.weights, 2)
            scrambled_sym = Circuit(H.dimensions, [C.inv(), sym, C]).composite_gate()
            # , f"\n{H.to_standard_form().__str__()}\n{sym.act(H).to_standard_form().__str__()}"
            assert H.is_close(scrambled_sym.act(H), literal=False), "Scrambled Hamiltonian not symmetric."

            known_F = scrambled_sym.symplectic
            if np.array_equal(known_F, np.eye(known_F.shape[0], dtype=known_F.dtype)) or H.n_paulis() <= 2 * n_qudits:
                # Trivial symmetry, or incomplete basis, skipping test
                continue
            else:
                circ = find_clifford_symmetries(H)

                assert len(
                    circ) != 0, (f"No symmetries found for run {_} \n F:{known_F}"
                                 " \n H:\n{H.tableau}\nHF:\n{scrambled_sym.act(H).tableau}")

                for c in circ:

                    H_s = H.to_standard_form()
                    H_out = c.act(H).to_standard_form()
                    H_s.weight_to_phase()
                    H_out.weight_to_phase()
                    assert np.all(H_s.tableau == H_out.tableau)
                    assert np.all(H_s.phases == H_out.phases)
                    assert np.all(H_s.weights == H_out.weights)

                    assert c.act(H).to_standard_form() == H.to_standard_form()

    def test_random_arbitrary_symmetry_with_block_decomposition(self):

        n_tests = 50
        p = 2
        n_qudits = 10
        n_paulis = 25

        for _ in range(n_tests):
            sym = Circuit.from_random(10, [p] * n_qudits)  #
            sym = sym.composite_gate()
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(sym, n_qudits, n_paulis, scrambled=False)
            C = Circuit.from_random(100, H.dimensions).composite_gate()  # scrambling circuit
            qc = qudit_cost(sym)
            H = C.act(H)
            H.weight_to_phase()
            H.weights = np.round(H.weights, 2)
            scrambled_sym = Circuit(H.dimensions, [C.inv(), sym, C]).composite_gate()
            # , f"\n{H.to_standard_form().__str__()}\n{sym.act(H).to_standard_form().__str__()}"
            assert H.is_close(scrambled_sym.act(H), literal=False), "Scrambled Hamiltonian not symmetric."

            known_F = scrambled_sym.symplectic
            if np.array_equal(known_F, np.eye(known_F.shape[0], dtype=known_F.dtype)) or H.n_paulis() <= 2 * n_qudits:
                # Trivial symmetry, or incomplete basis, skipping test
                continue
            else:
                F, S, T = min_qudit_clifford_symmetry(H)

                assert F == Circuit(F.dimensions, [T.inv(), S, T]).composite_gate()

                assert H.to_standard_form() == F.act(H).to_standard_form()
                assert T.act(S.act(T.inv().act(H))).to_standard_form() == H.to_standard_form()

                assert S.act(T.inv().act(H)).to_standard_form() == T.inv().act(H).to_standard_form()
                assert qudit_cost(S) <= qc

    def test_generate_symmetric_hamiltonian(self):
        n_qudits = 3
        n_paulis = 8
        p = 2
        n_tests = 1000
        for _ in range(n_tests):
            sym = Circuit.from_random(10, [p] * n_qudits)  #
            sym = sym.composite_gate()
            # sym = SWAP(0, 1, p)
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(sym, n_qudits, n_paulis, scrambled=False)
            # print(H.tableau)
            # H = H.to_standard_form()
            # print(H.to_standard_form().tableau)
            assert H.is_close(sym.act(H), literal=False), "Hamiltonian not symmetric. \n H: \n" + \
                str(H.to_standard_form().tableau) + "\n sym: \n" + str(sym.act(H).to_standard_form().tableau)
            C = Circuit.from_random(100, H.dimensions).composite_gate()  # scrambling circuit
            H = C.act(H)
            H.weight_to_phase()
            H.weights = np.round(H.weights, 2)
            scrambled_sym = Circuit(H.dimensions, [C.inv(), sym, C]).composite_gate()
            assert H.is_close(scrambled_sym.act(H), literal=False), "Scrambled Hamiltonian not symmetric."
