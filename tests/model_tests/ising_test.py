from sympleq.models.Ising import ising_chain_hamiltonian, ising_2d_hamiltonian
from sympleq.core.paulis import PauliSum


class TestModel:
    def test_ising_chain(self):
        J_zz = 1.0
        h_x = 1.0
        periodic = False
        for n_spins in range(1, 5):
            hamiltonian = ising_chain_hamiltonian(n_spins, J_zz, h_x, periodic)
            assert isinstance(hamiltonian, PauliSum)

            # Additional checks can be added here based on expected properties of the model

        periodic = True
        for n_spins in range(2, 5):
            hamiltonian = ising_chain_hamiltonian(n_spins, J_zz, h_x, periodic)
            assert isinstance(hamiltonian, PauliSum)

    def test_ising_2d(self):
        J_zz = 1.0
        h_x = 1.0
        periodic = False
        for n_x in range(1, 4):
            for n_y in range(1, 4):
                hamiltonian = ising_2d_hamiltonian(n_x, n_y, J_zz, h_x, periodic)
                assert isinstance(hamiltonian, PauliSum)

        periodic = True
        for n_x in range(2, 4):
            for n_y in range(2, 4):
                hamiltonian = ising_2d_hamiltonian(n_x, n_y, J_zz, h_x, periodic)
                assert isinstance(hamiltonian, PauliSum)
