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
