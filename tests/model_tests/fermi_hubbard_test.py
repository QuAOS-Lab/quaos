from sympleq.models.fermi_hubbard import fermi_hubbard_model
from sympleq.core.paulis import PauliSum


class TestModel:
    def test_fermi_hubbard_model(self):
        t, U, mu = 1.0, 4.0, 0.0
        periodic = False
        spinless = False
        for Lx in range(1, 5):
            for Ly in range(1, 5):
                tableau = fermi_hubbard_model(x_dimension=Lx, y_dimension=Ly,
                                              tunneling=t, coulomb=U, chemical_potential=mu,
                                              periodic=periodic, spinless=spinless)
                assert isinstance(tableau, PauliSum)

                # Additional checks can be added here based on expected properties of the model
