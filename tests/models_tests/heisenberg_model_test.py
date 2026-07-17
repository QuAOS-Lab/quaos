import numpy as np

from sympleq.models.heisenberg import (
    heisenberg_2d_hamiltonian,
    modified_heisenberg_ladder_hamiltonian,
)


def _embed_single_qubit(pauli: np.ndarray, site: int, n_qubits: int) -> np.ndarray:
    ops = [np.eye(2, dtype=complex) for _ in range(n_qubits)]
    ops[site] = pauli
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def _embed_two_qubit(pauli_a: np.ndarray, site_a: int,
                     pauli_b: np.ndarray, site_b: int,
                     n_qubits: int) -> np.ndarray:
    ops = [np.eye(2, dtype=complex) for _ in range(n_qubits)]
    ops[site_a] = pauli_a
    ops[site_b] = pauli_b
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def test_heisenberg_2d_hamiltonian_matches_open_2x2_reference():
    n_x = 2
    n_y = 2
    J = 1.7
    h_z = np.array([0.3, -0.4, 0.2, 0.1], dtype=float)

    H = heisenberg_2d_hamiltonian(n_x, n_y, J, h_z=h_z, periodic=False)
    H_dense = H.to_hilbert_space().toarray()

    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    expected = np.zeros_like(H_dense)
    bonds = [(0, 1), (0, 2), (1, 3), (2, 3)]
    for i, j in bonds:
        expected += J * _embed_two_qubit(X, i, X, j, 4)
        expected += J * _embed_two_qubit(Y, i, Y, j, 4)
        expected += J * _embed_two_qubit(Z, i, Z, j, 4)

    for i, field in enumerate(h_z):
        expected += field * _embed_single_qubit(Z, i, 4)

    np.testing.assert_allclose(H_dense, expected, atol=1e-12)


def test_heisenberg_2d_hamiltonian_periodic_2x1_counts_two_wrap_bonds():
    H = heisenberg_2d_hamiltonian(2, 1, 1.0, periodic=True).to_standard_form()

    expected_weights = np.array([2.0, 2.0, -2.0], dtype=complex)
    expected_tableau = np.array([
        [1, 1, 0, 0],
        [0, 0, 1, 1],
        [1, 1, 1, 1],
    ], dtype=np.uint8)
    expected_phases = np.array([0, 0, 0], dtype=int)

    assert np.array_equal(H.tableau, expected_tableau)
    assert np.array_equal(H.phases, expected_phases)
    np.testing.assert_allclose(H.weights, expected_weights)


def test_modified_heisenberg_ladder_hamiltonian_adds_plaquette_diagonals():
    n_x = 2
    n_y = 2
    J = 1.3
    h_z = np.array([0.2, -0.1, 0.4, -0.3], dtype=float)

    H = modified_heisenberg_ladder_hamiltonian(n_x, n_y, J, h_z=h_z)
    H_dense = H.to_hilbert_space().toarray()

    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    expected = np.zeros_like(H_dense)
    bonds = [(0, 1), (0, 2), (1, 3), (2, 3), (0, 3), (1, 2)]
    for i, j in bonds:
        expected += J * _embed_two_qubit(X, i, X, j, 4)
        expected += J * _embed_two_qubit(Y, i, Y, j, 4)
        expected += J * _embed_two_qubit(Z, i, Z, j, 4)

    for i, field in enumerate(h_z):
        expected += field * _embed_single_qubit(Z, i, 4)

    np.testing.assert_allclose(H_dense, expected, atol=1e-12)
