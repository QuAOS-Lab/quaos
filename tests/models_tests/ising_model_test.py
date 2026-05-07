import numpy as np

from sympleq.core.paulis import PauliSum, PauliString
from sympleq.models.Ising import ising_2d_hamiltonian, modified_ising_ladder_hamiltonian


def test_ising_2d_hamiltonian_accepts_site_resolved_transverse_field():
    h_x = np.array([0.1, 0.2, 0.3, 0.4], dtype=float)
    H = ising_2d_hamiltonian(2, 2, J_zz=1.0, h_x=h_x, periodic=False).to_standard_form()

    x_terms = H.tableau[:, :4]
    z_terms = H.tableau[:, 4:]

    local_x_mask = np.sum(x_terms, axis=1) == 1
    local_x_mask &= np.sum(z_terms, axis=1) == 0

    local_x_terms = H.weights[local_x_mask]
    np.testing.assert_allclose(np.sort(local_x_terms.real), np.sort(h_x))


def test_modified_ising_ladder_hamiltonian_adds_plaquette_diagonals():
    n_x = 3
    n_y = 2
    J_zz = 1.7
    h_x = 0.25

    H_base = ising_2d_hamiltonian(n_x, n_y, J_zz=J_zz, h_x=h_x, periodic=False)
    H_modified = modified_ising_ladder_hamiltonian(n_x, n_y, J_zz=J_zz, h_x=h_x)

    dims = [2] * (n_x * n_y)

    def site_index(x: int, y: int) -> int:
        return y * n_x + x

    diagonal_terms: list[PauliString] = []
    diagonal_weights: list[float] = []
    z0 = np.zeros(n_x * n_y, dtype=int)

    for x in range(n_x - 1):
        zz_a = z0.copy()
        zz_a[site_index(x, 0)] = 1
        zz_a[site_index(x + 1, 1)] = 1
        diagonal_terms.append(PauliString.from_exponents(z0, zz_a, dims))
        diagonal_weights.append(J_zz)

        zz_b = z0.copy()
        zz_b[site_index(x + 1, 0)] = 1
        zz_b[site_index(x, 1)] = 1
        diagonal_terms.append(PauliString.from_exponents(z0, zz_b, dims))
        diagonal_weights.append(J_zz)

    expected = H_base + PauliSum.from_pauli_strings(diagonal_terms, weights=diagonal_weights, phases=None)
    expected.combine_equivalent_paulis()
    expected.remove_zero_weight_paulis()

    assert H_modified.to_standard_form().is_close(expected.to_standard_form(), literal=False)
