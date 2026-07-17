import numpy as np
from sympleq.models.random_hamiltonian import (random_pauli_hamiltonian, random_gate_symmetric_hamiltonian,
                                               random_pauli_symmetry_hamiltonian)
from sympleq.core.circuits.gates import GATES
from sympleq.core.paulis import PauliSum


class TestModel:

    def test_random_pauli_hamiltonian(self):
        n_paulis = 5
        dimensions = [2, 2, 2]
        P = random_pauli_hamiltonian(n_paulis, dimensions)
        assert isinstance(P, PauliSum)

        dimensions = [2, 2, 3]
        P = random_pauli_hamiltonian(n_paulis, dimensions)
        assert isinstance(P, PauliSum)

        dimensions = [3, 3, 3]
        P = random_pauli_hamiltonian(n_paulis, dimensions)
        assert isinstance(P, PauliSum)

        dimensions = [3, 5, 7]
        P = random_pauli_hamiltonian(n_paulis, dimensions)
        assert isinstance(P, PauliSum)

    def test_random_gate_symmetric_hamiltonian(self):
        dimension = 2
        n_qudits = 15
        n_paulis = 6
        qudit_indices = (0, 1)
        P = random_gate_symmetric_hamiltonian(
            GATES.SWAP, dimension, qudit_indices, n_qudits, n_paulis, scrambled=False)
        assert isinstance(P, PauliSum)

    def test_random_pauli_symmetry_hamiltonian(self):
        n_tests = 10
        n_qudits = 5
        n_paulis = 30
        for _ in range(n_tests):
            n_redundant = np.random.randint(0, n_qudits - 3)
            n_conditional = np.random.randint(0, n_qudits - n_redundant - 1)
            ham = random_pauli_symmetry_hamiltonian(n_qudits, n_paulis, n_redundant=n_redundant,
                                                    n_conditional=n_conditional)
            assert isinstance(ham, PauliSum)