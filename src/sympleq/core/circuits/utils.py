from __future__ import annotations
import numpy as np
import scipy.sparse as sp
from functools import reduce
from sympleq._typing import IntNDArray
from sympleq.core.paulis._typing import TableauType, PhasesType, HilbertOperator, DimensionsType


def is_symplectic(F: TableauType, p: int) -> bool:
    """
    Check if matrix F is symplectic over GF(p).

    Args:
        F: (2n x 2n) numpy array, entries in {0, 1, ..., p-1}.
        p: prime modulus.

    Returns:
        True if F is symplectic over GF(p), False otherwise.
    """
    n = F.shape[0] // 2
    Omega = np.zeros((2 * n, 2 * n), dtype=int)
    Omega[:n, n:] = np.eye(n, dtype=int)
    Omega[n:, :n] = -np.eye(n, dtype=int)

    Omega = Omega % p

    lhs = (F.T @ Omega @ F) % p
    return np.array_equal(lhs, Omega)


def symplectic_product_arrays(u: TableauType, v: TableauType, p: int = 2) -> int:
    """
    Compute the symplectic inner product of two binary vectors.

    Args:
        u, v: Binary vectors of length 2n

    Returns:
        Symplectic inner product modulo p
    """
    n = len(u) // 2
    return int(np.sum(u[:n] * v[n:] - u[n:] * v[:n])) % p


def symplectic_product_matrix(pauli_sum_tableau: TableauType, p: int = 2) -> TableauType:
    m = len(pauli_sum_tableau)
    spm = np.zeros((m, m), dtype=int)

    for i in range(m):
        for j in range(m):
            spm[i, j] = symplectic_product_arrays(pauli_sum_tableau[i], pauli_sum_tableau[j], p)

    return spm


def symplectic_form(n: int, p: int = 2) -> TableauType:
    """
    Construct the symplectic matrix Omega for a given dimension n over GF(p).

    Args:
        n: Half the length of the vectors (2n is the full length)
        p: Modulus (default 2)

    Returns:
        Omega matrix as a 2n x 2n numpy array
    """
    Id = np.eye(n, dtype=int)
    Z = np.zeros((n, n), dtype=int)

    if p == 2:
        return np.block([[Z, Id], [Id, Z]])
    else:
        return np.block([[Z, Id], [-Id, Z]])


def transvection_matrix(h: IntNDArray, p: int = 2, multiplier: int = 1) -> TableauType:
    """
    Compute the transvection matrix corresponding to the vector h.

    Args:
        h: Binary vector of length 2n
        p: Modulus (default 2)

    Returns:
        The transvection matrix as a 2n x 2n matrix over integers modulo p
    """
    n = len(h) // 2
    Omega = symplectic_form(n, p)

    F_h = (np.eye(2 * n, dtype=int) + multiplier * (Omega @ np.outer(h.T, h))) % p
    return F_h


def transvection(h, x, p=2):
    return (x + symplectic_product_arrays(x, h.T, p) * h) % p


def embed_symplectic(symplectic_local: TableauType, phase_vector_local: PhasesType,
                     qudit_indices: tuple[int, ...] | list[int] | int, n_qudits: int) -> tuple[TableauType, PhasesType]:
    """
    Embed a local Clifford (F_local, h_local) into a larger 2n-dimensional space,
    correctly handling arbitrary qudit index ordering.
    """
    _qudit_indices = np.asarray(qudit_indices, dtype=int).reshape(-1)

    m = len(_qudit_indices)
    if symplectic_local.shape != (2 * m, 2 * m):
        raise ValueError("symplectic_local must be 2m x 2m")
    if len(phase_vector_local) != 2 * m:
        raise ValueError("phase_vector_local must have length 2m")

    # Full 2n x 2n identity
    F_full = np.eye(2 * n_qudits, dtype=int)
    h_full = np.zeros(2 * n_qudits, dtype=int)

    # Build row/column index mapping for the full space
    # First X rows/columns
    row_indices = np.concatenate([_qudit_indices, n_qudits + _qudit_indices])
    col_indices = np.concatenate([_qudit_indices, n_qudits + _qudit_indices])

    # Place the full local symplectic block into the full system
    F_full[np.ix_(row_indices, col_indices)] = symplectic_local

    # Embed phase vector
    h_full[row_indices] = phase_vector_local

    return F_full, h_full


def _multi_index_to_linear(index: list[int] | np.ndarray, dims: list[int] | np.ndarray) -> int:
    """Convert a mixed-radix index to a linear index using row-major order.

    idx(i0, i1, ..., iN-1) = sum_k i_k * prod_{l>k} dims[l]
    """
    dims = list(map(int, dims))
    idx = 0
    # Compute strides from right to left
    strides = [1] * len(dims)
    for k in range(len(dims) - 2, -1, -1):
        strides[k] = strides[k + 1] * dims[k + 1]
    for k, ik in enumerate(index):
        idx += int(ik) * strides[k]
    return idx


def embed_unitary(U_local: HilbertOperator,
                  qudit_indices: tuple[int, ...] | list[int] | np.ndarray,
                  total_dimensions: DimensionsType) -> HilbertOperator:
    """
    Embed a local unitary acting on a subset of qudits into the full Hilbert space.

    - Basis ordering: |q0> ⊗ |q1> ⊗ ... ⊗ |qN-1>
    - Linear index mapping: idx(q) = sum_k q[k] * prod_{l>k} d[l]

    The local unitary is assumed to act on qudits in the order given by
    `qudit_indices`.

    Parameters
    ----------
    U_local : HilbertOperator
        Local unitary of shape ``(D_loc, D_loc)`` where
        ``D_loc = prod(d[qudit_indices])``.
    qudit_indices : tuple[int, ...] | list[int] | np.ndarray
        Indices of the qudits the local unitary acts on.
    total_dimensions : DimensionsType
        Dimensions of each qudit in the full system.

    Returns
    -------
    HilbertOperator
        Full unitary of shape ``(D_total, D_total)`` with
        ``D_total = prod(total_dimensions)``.
    """
    dims = np.asarray(total_dimensions)
    N = len(dims)
    sel = list(qudit_indices)
    rest = [k for k in range(N) if k not in sel]
    perm_order = sel + rest

    D_rest = int(np.prod(dims[rest])) if rest else 1
    D_total = int(np.prod(dims))

    ## Build permutation matrix P that reorders tensor factors to [sel..., rest...].
    # Construct in COO-style triplets and convert once to CSR to avoid expensive
    # repeated structural updates on CSR.
    dims_perm = [dims[k] for k in sel + rest]
    n_states = D_loc_expected * D_rest
    rows = np.empty(n_states, dtype=int)
    cols = np.empty(n_states, dtype=int)
    data = np.ones(n_states, dtype=complex)

    for idx, q in enumerate(np.ndindex(*dims)):
        q = list(q)
        old_idx = _multi_index_to_linear(q, dims)
        q_perm = [q[k] for k in (sel + rest)]
        new_idx = _multi_index_to_linear(q_perm, dims_perm)
        rows[idx] = new_idx
        cols[idx] = old_idx
    P = sp.csr_matrix((data, (rows, cols)), shape=(n_states, D_total))

    # Construct full operator: P^T (U_local ⊗ I_rest) P
    U_kron = sp.kron(U_local, sp.eye(D_rest))
    return P.conj().T @ U_kron @ P


def tensor(mm: list[HilbertOperator]) -> HilbertOperator:
    # Inputs:
    #     mm - (list{scipy.sparse.csr_matrix}) - matrices to tensor
    # Outputs:
    #     (scipy.sparse.csr_matrix) - tensor product of matrices
    if len(mm) == 0:
        return sp.csr_matrix([])

    if len(mm) == 1:
        return mm[0]

    return sp.csr_matrix(sp.kron(mm[0], tensor(mm[1:]), format="csr"))


def H_mat(d: int) -> sp.csr_matrix:
    omega = np.exp(2 * np.pi * 1j / d)
    return sp.csr_matrix(1 / np.sqrt(d) * np.array([[omega ** (i0 * i1) for i0 in range(d)] for i1 in range(d)]))


def S_mat(d: int) -> HilbertOperator:
    if d == 2:
        return sp.csr_matrix(np.diag([1, 1j]))

    omega = np.exp(2 * np.pi * 1j / d)
    return sp.csr_matrix(np.diag([omega ** (i * (i - 1) / 2) for i in range(d)]))


def _X_power(d: int, a: int) -> HilbertOperator:
    """
    Sparse X^a on a single qudit (dimension d).
    Places ones at rows ((j + a) % d) and columns j.
    """
    a %= d
    cols = np.arange(d, dtype=int)
    rows = (cols + a) % d
    data = np.ones(d, dtype=complex)
    return sp.csr_matrix((data, (rows, cols)), shape=(d, d))


def _Z_power(d: int, b: int) -> HilbertOperator:
    """
    Sparse Z^b on a single qudit (dimension d).
    Diagonal with entries ω^{b*j}, ω = exp(2πi/d).
    """
    b %= d
    j = np.arange(d)
    omega = np.exp(2j * np.pi / d)
    diag = omega ** (b * j)
    return sp.csr_matrix(sp.diags(diag, offsets=0, dtype=complex, format="csr"))


def pauli_unitary_qudit(d: int, x: int, z: int, convention: str = "bare") -> HilbertOperator:
    """
    Sparse unitary for single-qudit Pauli specified by tableau [x | z] over Z_d.

    Conventions:
      - "bare": U = Z^z X^x
      - "weyl": U = τ^{x z} Z^z X^x,  τ = exp(iπ(d+1)/d)
    """
    Xx = _X_power(d, x)
    Zz = _Z_power(d, z)
    U = Zz @ Xx  # Z then X

    if convention.lower() == "weyl":
        tau = np.exp(1j * np.pi * (d + 1) / d)
        U = (tau ** (x * z)) * U
    elif convention.lower() != "bare":
        raise ValueError("convention must be 'bare' or 'weyl'.")

    return sp.csr_matrix(U)


def pauli_unitary_from_tableau(
    d: int, x: TableauType, z: TableauType, convention: str = "bare"
) -> HilbertOperator:
    """
    Sparse multi-qudit unitary. x, z are length-n integer arrays (mod d),
    representing a tableau row [x0..x_{n-1} | z0..z_{n-1}] on n qudits.

    Returns a (d^n)×(d^n) CSR sparse matrix:  ⊗_k (Z^{z_k} X^{x_k})
    with optional Weyl phase τ^{x_k z_k} per local factor.
    """
    x = np.asarray(x, dtype=int)
    z = np.asarray(z, dtype=int)
    assert x.shape == z.shape and x.ndim == 1, "x and z must be 1D arrays of same length"

    locals_ = [
        pauli_unitary_qudit(d, int(xk), int(zk), convention=convention)
        for xk, zk in zip(x, z)
    ]
    # Tensor product (left-to-right order matches locals_ order)
    U = reduce(lambda A, B: sp.kron(A, B, format="csr"), locals_)
    return sp.csr_matrix(U)
