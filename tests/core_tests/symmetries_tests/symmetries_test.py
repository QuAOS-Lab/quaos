from sympleq.core.circuits.gates import GATES, Gate
from sympleq.models.random_hamiltonian import random_gate_symmetric_hamiltonian
from sympleq.core.symmetries.clifford import find_clifford_symmetries  # , qudit_cost, min_qudit_clifford_symmetry
from sympleq.core.circuits import Circuit
import numpy as np


class TestSymmetryFinder:

    def test_random_SWAP_symmetry(self):
        n_tests = 30
        dimension = 2
        n_qudits = 15
        # Need enough terms to determine a non-trivial automorphism robustly.
        n_paulis = 6
        for _ in range(n_tests):
            qudit_indices = (0, 1)
            all_qudit_indices = tuple(range(n_qudits))
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(
                GATES.SWAP, dimension, qudit_indices, n_qudits, n_paulis, scrambled=False)

            C_gate = Gate.from_random(n_qudits, dimension, 100)  # scrambling gate
            H = C_gate.act(H, all_qudit_indices)
            H.weight_to_phase()

            scrambled_sym = Circuit.from_gates_and_qudits(
                H.dimensions,
                [C_gate.inverse(), GATES.SWAP, C_gate],
                [all_qudit_indices, qudit_indices, all_qudit_indices]).composite_gate()

            check = scrambled_sym.act(H, all_qudit_indices)
            assert H.is_close(check, literal=False), f"\n{H}\n{check}"

            symmetries = find_clifford_symmetries(H)

            assert len(symmetries) != 0

            for sym_gate in symmetries:
                H_s = H.to_standard_form()
                H_out = sym_gate.act(H, all_qudit_indices).to_standard_form()
                H_s.weight_to_phase()
                H_out.weight_to_phase()
                assert np.all(H_s.tableau == H_out.tableau)
                assert np.all(H_s.phases == H_out.phases)
                assert np.all(H_s.weights == H_out.weights)

                assert sym_gate.act(H, all_qudit_indices).is_close(H, literal=False)

    def test_random_multi_SWAP_symmetry(self):

        n_tests = 100
        dimension = 2
        n_qudits = 3
        n_paulis = 7
        all_qudit_indices = tuple(range(n_qudits))
        for _ in range(n_tests):
            sym = Circuit.from_gates_and_qudits(
                [dimension] * n_qudits,
                [GATES.SWAP, GATES.SWAP],
                [(0, 1), (1, 2)])
            sym = sym.composite_gate()
            # unscrambled H
            H = random_gate_symmetric_hamiltonian(sym, dimension, all_qudit_indices, n_qudits, n_paulis,
                                                  scrambled=False)
            C = Gate.from_random(n_qudits, dimension, 100)  # scrambling gate
            H = C.act(H, all_qudit_indices)
            H.weight_to_phase()
            scrambled_sym = Circuit.from_gates_and_qudits(H.dimensions,
                                                          [C.inverse(), sym, C],
                                                          [all_qudit_indices, all_qudit_indices, all_qudit_indices])
            scrambled_sym = scrambled_sym.composite_gate()

            assert H.to_standard_form() == scrambled_sym.act(H, all_qudit_indices).to_standard_form(
            ), f"\n{H.to_standard_form().__str__()}\n{sym.act(H, all_qudit_indices).to_standard_form().__str__()}"
            circ = find_clifford_symmetries(H)

            assert len(circ) != 0

            for c in circ:
                H_s = H.to_standard_form()
                H_out = c.act(H, all_qudit_indices).to_standard_form()
                H_s.weight_to_phase()
                H_out.weight_to_phase()
                assert np.all(H_s.tableau == H_out.tableau)
                assert np.all(H_s.phases == H_out.phases)
                assert np.all(H_s.weights == H_out.weights)

                assert c.act(H, all_qudit_indices).to_standard_form() == H.to_standard_form()

    def test_generate_symmetric_hamiltonian(self):
        n_qudits = 5
        n_paulis = 12
        dimension = 2
        all_qudit_indices = tuple(range(n_qudits))
        n_tests = 100
        for _ in range(n_tests):
            C1 = Circuit.from_random(10, [dimension] * n_qudits)
            C1_gate = C1.composite_gate()

            H = random_gate_symmetric_hamiltonian(C1_gate, dimension, all_qudit_indices,
                                                  n_qudits, n_paulis, scrambled=False)
            check = C1_gate.act(H, all_qudit_indices)
            assert H.is_close(check, literal=False), "Hamiltonian not symmetric. \n H: \n" + \
                H.to_standard_form().tableau + "\n sym: \n" + \
                check.to_standard_form().tableau

            C2_gate = Gate.from_random(n_qudits, dimension, 100)  # scrambling gate
            H = C2_gate.act(H, all_qudit_indices)
            H.weight_to_phase()
            H.weights = np.round(H.weights, 2)
            scrambled_C = Circuit.from_gates_and_qudits(H.dimensions,
                                                        [C2_gate.inverse(), C1_gate, C2_gate],
                                                        [all_qudit_indices, all_qudit_indices, all_qudit_indices])
            scrambled_C_gate = scrambled_C.composite_gate()
            assert H.is_close(scrambled_C_gate.act(H, all_qudit_indices),
                              literal=False), "Scrambled Hamiltonian not symmetric."

    def test_random_arbitrary_symmetry(self):
        n_tests = 10
        dimension = 2
        n_qudits = 10
        n_paulis = 30
        all_qudit_indices = tuple(range(n_qudits))

        for _ in range(n_tests):
            C1 = Circuit.from_random(10, [dimension] * n_qudits)
            C1_gate = C1.composite_gate()

            H = random_gate_symmetric_hamiltonian(C1_gate, dimension, all_qudit_indices,
                                                  n_qudits, n_paulis, scrambled=False)
            check = C1_gate.act(H, all_qudit_indices)
            assert H.is_close(check, literal=False), "Hamiltonian not symmetric. \n H: \n" + \
                H.to_standard_form().tableau + "\n sym: \n" + \
                check.to_standard_form().tableau

            C2_gate = Gate.from_random(n_qudits, dimension, 100)  # scrambling gate

            H = C2_gate.act(H, all_qudit_indices)
            H.weight_to_phase()
            H.weights = np.round(H.weights, 2)
            scrambled_C = Circuit.from_gates_and_qudits(H.dimensions,
                                                        [C2_gate.inverse(), C1_gate, C2_gate],
                                                        [all_qudit_indices, all_qudit_indices, all_qudit_indices])
            scrambled_C_gate = scrambled_C.composite_gate()
            assert H.is_close(scrambled_C_gate.act(H, all_qudit_indices),
                              literal=False), "Scrambled Hamiltonian not symmetric."
            circ = find_clifford_symmetries(H)

            assert len(circ) != 0

            for c in circ:
                H_s = H.to_standard_form()
                H_out = c.act(H, all_qudit_indices).to_standard_form()
                H_s.weight_to_phase()
                H_out.weight_to_phase()
                assert np.all(H_s.tableau == H_out.tableau)
                assert np.all(H_s.phases == H_out.phases)
                assert np.all(H_s.weights == H_out.weights)

                assert c.act(H, all_qudit_indices).to_standard_form() == H.to_standard_form()
