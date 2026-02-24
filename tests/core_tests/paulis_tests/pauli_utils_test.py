import numpy as np
import random
import pytest
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.paulis import PauliSum, check_mappable_via_clifford, mod_inv, hamiltonian_mean, covariance_matrix
# from sympleq.core.paulis.utils import row_reduce_mod_d
from tests import choose_random_dimensions


N_tests = 30


class TestUtils:

    def test_check_mappable_via_clifford(self):
        for _ in range(N_tests):
            dimensions = choose_random_dimensions(250)
            n_qudits = len(dimensions)
            n_paulis = random.randint(1, max(2, 2 * n_qudits ** 2))

            P1 = PauliSum.from_random(n_paulis=n_paulis, dimensions=dimensions)
            C = Circuit.from_random(n_gates=10 * n_qudits ** 2, dimensions=dimensions)
            P2 = C.act(P1)

            assert check_mappable_via_clifford(P1, P2), (
                f"Expected mapped PauliSum to be Clifford-mappable.\nP1:\n{P1}\nP2:\n{P2}"
            )

            # Build a modified target by adding a new PauliString that is not
            # already in P2; this should break mappability.
            existing_tableau = P2.tableau
            while True:
                candidate = PauliSum.from_random(1, dimensions)
                candidate_tableau = candidate.tableau
                if not any(np.array_equal(candidate_tableau[0], row) for row in existing_tableau):
                    extra_pauli = candidate
                    break

            P2_not_mappable = P2 + extra_pauli

            assert not check_mappable_via_clifford(P1, P2_not_mappable), (
                "Expected non-mappability after adding a new PauliString to P2. "
                f"Added: {extra_pauli}"
            )

    def test_mod_inv(self):
        for _ in range(N_tests):
            d = random.randint(2, 250)
            a = random.randint(1, d - 1)
            inv_1 = None
            for i in range(1, d):
                if (a * i) % d == 1:
                    inv_1 = i
                    break

            if inv_1 is None:
                with pytest.raises(ValueError):
                    mod_inv(a, d)
            else:
                inv_2 = mod_inv(a, d)
                assert inv_1 == inv_2, (f"Expected modular inverse to yield "
                                        f"{inv_1}, yet we got {inv_2} for a={a}, d={d}.")

    def test_hamiltonian_mean(self):
        for _ in range(N_tests):

            dimensions = choose_random_dimensions(25)
            m_size = int(np.prod(dimensions))

            matrix = np.random.rand(m_size, m_size) + 1j * np.random.rand(m_size, m_size) + \
                - (1 / 2) * (1 + 1j) * np.ones((m_size, m_size))
            matrix = matrix + matrix.conj().T

            P = PauliSum.from_hilbert_space(matrix, dimensions=dimensions)

            state = np.random.rand(m_size) + 1j * np.random.rand(m_size) + \
                - (1 / 2) * (1 + 1j) * np.ones((m_size))
            state /= np.linalg.norm(state)

            mean_1 = state.conj().T @ matrix @ state
            mean_2 = hamiltonian_mean(P, state)

            assert np.isclose(mean_1, mean_2), (
                f"Expected Hamiltonian mean to match direct computation. "
                f"Got {mean_1} vs {mean_2} for state {state} and PauliSum:\n{P}"
            )

    def test_covariance_matrix(self):
        for _ in range(N_tests):

            dimensions = choose_random_dimensions(25)
            m_size = int(np.prod(dimensions))

            matrix = np.random.rand(m_size, m_size) + 1j * np.random.rand(m_size, m_size) + \
                - (1 / 2) * (1 + 1j) * np.ones((m_size, m_size))
            matrix = matrix + matrix.conj().T

            P = PauliSum.from_hilbert_space(matrix, dimensions=dimensions)

            state = np.random.rand(m_size) + 1j * np.random.rand(m_size) + \
                - (1 / 2) * (1 + 1j) * np.ones((m_size))
            state /= np.linalg.norm(state)
            state_dag = state.conj().T

            cov_1 = state_dag @ matrix @ matrix @ state - (state_dag @ matrix @ state) ** 2
            cov_mat = covariance_matrix(P, state)
            cov_2 = np.sum(cov_mat)

            assert np.isclose(cov_1, cov_2), (
                f"Expected covariance matrix to match direct computation. "
                f"Got {cov_1} vs {cov_2} for state {state} and PauliSum:\n{P}"
            )
