import numpy as np

from sympleq.models.Ising import ising_2d_hamiltonian


def test_ising_2d_hamiltonian_accepts_site_resolved_transverse_field():
    h_x = np.array([0.1, 0.2, 0.3, 0.4], dtype=float)
    H = ising_2d_hamiltonian(2, 2, J_zz=1.0, h_x=h_x, periodic=False).to_standard_form()

    x_terms = H.tableau[:, :4]
    z_terms = H.tableau[:, 4:]

    local_x_mask = np.sum(x_terms, axis=1) == 1
    local_x_mask &= np.sum(z_terms, axis=1) == 0

    local_x_terms = H.weights[local_x_mask]
    np.testing.assert_allclose(np.sort(local_x_terms.real), np.sort(h_x))
