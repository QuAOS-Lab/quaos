import numpy as np
from typing import Tuple, Optional
from sympleq.core.paulis import PauliSum
from .utils import _zeros_xz, _add_pauli_term, _finalize_terms


def _all_to_all_heisenberg_hamiltonian(
    n_qubits: int,
    J: float | np.ndarray,
    h_z: Optional[np.ndarray] = None,
    tol: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Build a Pauli-tableau representation of the all-to-all Heisenberg model

        H = Σ_{i<j} J_{ij} (X_i X_j + Y_i Y_j + Z_i Z_j) + Σ_i h_i Z_i

    If J is a scalar, uses uniform J_{ij}=J on the complete graph.
    If J is an array, it must be shape (n,n) and only the upper triangle i<j is used.

    Output convention matches your code:
        P = (1j)**phase * Π_k X_k^{x_k} Z_k^{z_k}
    """
    n = int(n_qubits)
    if n <= 0:
        raise ValueError("n_qubits must be positive.")

    # Parse couplings
    if np.isscalar(J):
        J_mat = np.full((n, n), J, dtype=float)
        np.fill_diagonal(J_mat, 0.0)
    else:
        J_mat = np.asarray(J, dtype=float)
        if J_mat.shape != (n, n):
            raise ValueError(f"J matrix must have shape ({n},{n}), got {J_mat.shape}.")
        # symmetrize defensively
        J_mat = 0.5 * (J_mat + J_mat.T)
        np.fill_diagonal(J_mat, 0.0)

    if h_z is None:
        h_vec = np.zeros(n, dtype=float)
    else:
        h_vec = np.asarray(h_z, dtype=float).reshape(-1)
        if h_vec.shape[0] != n:
            raise ValueError(f"h_z must have length {n}, got {h_vec.shape[0]}.")

    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], float] = {}

    # Pair terms
    for i in range(n):
        for j in range(i + 1, n):
            Jij = J_mat[i, j]
            if Jij == 0.0:
                continue

            # XX
            x, z = _zeros_xz(n)
            x[i] = 1
            x[j] = 1
            _add_pauli_term(terms, x, z, phase=0, coeff=Jij)

            # YY = (iXZ)(iXZ) = i^2 XZ XZ  -> phase=2
            x, z = _zeros_xz(n)
            x[i] = 1
            z[i] = 1
            x[j] = 1
            z[j] = 1
            _add_pauli_term(terms, x, z, phase=2, coeff=Jij)

            # ZZ
            x, z = _zeros_xz(n)
            z[i] = 1
            z[j] = 1
            _add_pauli_term(terms, x, z, phase=0, coeff=Jij)

    # Local z-fields
    for i, hi in enumerate(h_vec):
        if hi == 0.0:
            continue
        x, z = _zeros_xz(n)
        z[i] = 1
        _add_pauli_term(terms, x, z, phase=0, coeff=hi)

    return _finalize_terms(n, terms, tol=tol)


def all_to_all_heisenberg_hamiltonian(
    n: int,
    J: float | np.ndarray,
    h_z: Optional[np.ndarray] = None,
    tol: float = 0.0,
) -> 'PauliSum':
    tableau, coeffs, phases = _all_to_all_heisenberg_hamiltonian(n, J, h_z=h_z, tol=tol)
    return PauliSum.from_tableau(tableau, weights=coeffs, phases=phases)


def _heisenberg_2d_hamiltonian(
    n_x: int,
    n_y: int,
    J: float | np.ndarray,
    h_z: float | np.ndarray = 0.0,
    periodic: bool = False,
    tol: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Build a Pauli-tableau representation of a nearest-neighbour 2D Heisenberg model
    on an ``n_x`` by ``n_y`` square lattice:

        H = \sum_{<i,j>} J_{ij} (X_i X_j + Y_i Y_j + Z_i Z_j) + \sum_i h_i Z_i

    Sites are indexed in row-major order, matching ``ising_2d_hamiltonian``:
        site_index(x, y) = y * n_x + x

    If ``J`` is a scalar, a uniform coupling is used on every nearest-neighbour bond.
    If ``J`` is an array, it must be a full ``(n_x*n_y, n_x*n_y)`` coupling matrix;
    only nearest-neighbour entries are used.
    """
    n_x = int(n_x)
    n_y = int(n_y)
    if n_x <= 0 or n_y <= 0:
        raise ValueError("n_x and n_y must be positive.")

    n = n_x * n_y

    if np.isscalar(J):
        J_mat = np.zeros((n, n), dtype=float)
        use_uniform_J = True
        uniform_J = float(J)
    else:
        J_mat = np.asarray(J, dtype=float)
        if J_mat.shape != (n, n):
            raise ValueError(f"J matrix must have shape ({n},{n}), got {J_mat.shape}.")
        J_mat = 0.5 * (J_mat + J_mat.T)
        np.fill_diagonal(J_mat, 0.0)
        use_uniform_J = False
        uniform_J = 0.0

    if np.isscalar(h_z):
        h_vec = np.full(n, float(h_z), dtype=float)
    else:
        h_vec = np.asarray(h_z, dtype=float).reshape(-1)
        if h_vec.shape[0] != n:
            raise ValueError(f"h_z must have length {n}, got {h_vec.shape[0]}.")

    def site_index(x: int, y: int) -> int:
        return y * n_x + x

    def coupling(i: int, j: int) -> float:
        return uniform_J if use_uniform_J else float(J_mat[i, j])

    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], float] = {}

    for x in range(n_x):
        for y in range(n_y):
            i = site_index(x, y)

            neighbours: list[int] = []
            if x + 1 < n_x:
                neighbours.append(site_index(x + 1, y))
            elif periodic and n_x > 1:
                neighbours.append(site_index(0, y))

            if y + 1 < n_y:
                neighbours.append(site_index(x, y + 1))
            elif periodic and n_y > 1:
                neighbours.append(site_index(x, 0))

            for j in neighbours:
                Jij = coupling(i, j)
                if Jij == 0.0:
                    continue

                x_term, z_term = _zeros_xz(n)
                x_term[i] = 1
                x_term[j] = 1
                _add_pauli_term(terms, x_term, z_term, phase=0, coeff=Jij)

                x_term, z_term = _zeros_xz(n)
                x_term[i] = 1
                z_term[i] = 1
                x_term[j] = 1
                z_term[j] = 1
                _add_pauli_term(terms, x_term, z_term, phase=2, coeff=Jij)

                x_term, z_term = _zeros_xz(n)
                z_term[i] = 1
                z_term[j] = 1
                _add_pauli_term(terms, x_term, z_term, phase=0, coeff=Jij)

    for i, hi in enumerate(h_vec):
        if hi == 0.0:
            continue
        x_term, z_term = _zeros_xz(n)
        z_term[i] = 1
        _add_pauli_term(terms, x_term, z_term, phase=0, coeff=hi)

    return _finalize_terms(n, terms, tol=tol)


def heisenberg_2d_hamiltonian(
    n_x: int,
    n_y: int,
    J: float | np.ndarray,
    h_z: float | np.ndarray = 0.0,
    periodic: bool = False,
    tol: float = 0.0,
) -> 'PauliSum':
    tableau, coeffs, phases = _heisenberg_2d_hamiltonian(
        n_x,
        n_y,
        J,
        h_z=h_z,
        periodic=periodic,
        tol=tol,
    )
    return PauliSum.from_tableau(tableau, weights=coeffs, phases=phases)
