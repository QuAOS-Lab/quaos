from sympleq.models.toric_code import ToricCode
from sympleq.core.paulis import PauliSum


class TestModel:
    def test_toric_code(self):
        for Nx in range(2, 6):
            for Ny in range(2, 6):
                periodic = True
                c_x = 1.0
                c_z = 1.0
                c_g = 1.0

                TC = ToricCode(Nx, Ny, c_x, c_z, c_g, periodic)
                hamiltonian = TC.hamiltonian()
                assert isinstance(hamiltonian, PauliSum)

                # Additional checks can be added here based on expected properties of the model

        for Nx in range(2, 6):
            for Ny in range(2, 6):
                periodic = False
                c_x = 1.0
                c_z = 1.0
                c_g = 1.0

                TC = ToricCode(Nx, Ny, c_x, c_z, c_g, periodic)
                hamiltonian = TC.hamiltonian()
                assert isinstance(hamiltonian, PauliSum)

                # Additional checks can be added here based on expected properties of the model
