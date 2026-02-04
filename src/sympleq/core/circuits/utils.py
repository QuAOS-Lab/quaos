from __future__ import annotations
import numpy as np
import scipy.sparse as sp
from sympleq.utils import int_to_bases, bases_to_int
from functools import reduce


def is_symplectic(F, p: int) -> bool:
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


def symplectic_product_arrays(u: np.ndarray, v: np.ndarray, p: int = 2) -> int:
    """
    Compute the symplectic inner product of two binary vectors.

    Args:
        u, v: Binary vectors of length 2n

    Returns:
        Symplectic inner product modulo p
    """
    n = len(u) // 2
    return (np.sum(u[:n] * v[n:] - u[n:] * v[:n])) % p


def symplectic_product_matrix(pauli_sum_tableau: np.ndarray, p: int = 2) -> np.ndarray:
    m = len(pauli_sum_tableau)
    spm = np.zeros((m, m), dtype=int)

    for i in range(m):
        for j in range(m):
            spm[i, j] = symplectic_product_arrays(pauli_sum_tableau[i], pauli_sum_tableau[j], p)

    return spm


def symplectic_form(n: int, p: int = 2) -> np.ndarray:
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


def transvection_matrix(h: np.ndarray, p=2, multiplier=1):
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


def embed_symplectic(symplectic_local, phase_vector_local, qudit_indices, n_qudits):
    """
    Embed a local Clifford (F_local, h_local) into a larger 2n-dimensional space,
    correctly handling arbitrary qudit index ordering.
    """
    m = len(qudit_indices)
    if symplectic_local.shape != (2 * m, 2 * m):
        raise ValueError("symplectic_local must be 2m x 2m")
    if len(phase_vector_local) != 2 * m:
        raise ValueError("phase_vector_local must have length 2m")

    # Full 2n x 2n identity
    F_full = np.eye(2 * n_qudits, dtype=int)
    h_full = np.zeros(2 * n_qudits, dtype=int)

    qudit_indices = np.array(qudit_indices, dtype=int)

    # Build row/column index mapping for the full space
    # First X rows/columns
    row_indices = np.concatenate([qudit_indices, n_qudits + qudit_indices])
    col_indices = np.concatenate([qudit_indices, n_qudits + qudit_indices])

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


def embed_unitary(U_local: sp.csr_matrix,
                  qudit_indices: list[int] | np.ndarray,
                  total_dimensions: list[int] | np.ndarray) -> sp.csr_matrix:
    """
    Embed a local unitary acting on a subset of qudits into the full Hilbert space.

    - Basis ordering: |q0> ⊗ |q1> ⊗ ... ⊗ |qN-1>
    - Linear index mapping: idx(q) = sum_k q[k] * prod_{l>k} d[l]

    The local unitary is assumed to act on qudits in the order given by `qudit_indices`.

    Args:
        U_local: Local unitary of shape (D_loc, D_loc) where D_loc = prod(d[qudit_indices]).
        qudit_indices: Indices of the qudits the local unitary acts on.
        total_dimensions: Dimensions of each qudit in the full system.

    Returns:
        Full unitary of shape (D_total, D_total) with D_total = prod(total_dimensions).
    """
    qudit_indices = list(map(int, qudit_indices))
    dims = list(map(int, total_dimensions))
    N = len(dims)
    sel = qudit_indices
    rest = [k for k in range(N) if k not in sel]

    # Validate local size matches product of selected dims
    D_loc_expected = int(np.prod([dims[k] for k in sel]))
    if U_local.shape != (D_loc_expected, D_loc_expected):
        raise ValueError(
            f"U_local shape {U_local.shape} doesn't match selected dims product {D_loc_expected}"
        )

    D_rest = int(np.prod([dims[k] for k in rest]) if rest else 1)
    D_total = int(np.prod(dims))

    # Build permutation matrix P that reorders tensor factors to [sel..., rest...]
    dims_perm = [dims[k] for k in sel + rest]
    P = sp.csr_matrix(np.zeros((D_loc_expected * D_rest, D_total), dtype=complex))

    # Iterate over all basis states
    for q in np.ndindex(*dims):
        q = list(q)
        old_idx = _multi_index_to_linear(q, dims)
        q_perm = [q[k] for k in (sel + rest)]
        new_idx = _multi_index_to_linear(q_perm, dims_perm)
        P[new_idx, old_idx] = 1.0

    # Construct full operator: P^T (U_local ⊗ I_rest) P
    U_kron = sp.kron(U_local, sp.eye(D_rest))
    return P.conj().T @ U_kron @ P


def tensor(mm: list[sp.csr_matrix]) -> sp.csr_matrix:
    # Inputs:
    #     mm - (list{scipy.sparse.csr_matrix}) - matrices to tensor
    # Outputs:
    #     (scipy.sparse.csr_matrix) - tensor product of matrices
    if len(mm) == 0:
        return sp.csr_matrix([])

    if len(mm) == 1:
        return mm[0]

    return sp.csr_matrix(sp.kron(mm[0], tensor(mm[1:]), format="csr"))


def I_mat(d: int) -> sp.csr_matrix:
    return sp.csr_matrix(np.diag([1] * d))


def H_mat(d: int) -> sp.csr_matrix:
    omega = np.exp(2 * np.pi * 1j / d)
    return sp.csr_matrix(1 / np.sqrt(d) * np.array([[omega ** (i0 * i1) for i0 in range(d)] for i1 in range(d)]))


def S_mat(d: int) -> sp.csr_matrix:
    if d == 2:
        return sp.csr_matrix(np.diag([1, 1j]))

    omega = np.exp(2 * np.pi * 1j / d)
    return sp.csr_matrix(np.diag([omega ** (i * (i - 1) / 2) for i in range(d)]))


def _X_power(d: int, a: int) -> sp.csr_matrix:
    """
    Sparse X^a on a single qudit (dimension d).
    Places ones at rows ((j + a) % d) and columns j.
    """
    a %= d
    cols = np.arange(d, dtype=int)
    rows = (cols + a) % d
    data = np.ones(d, dtype=complex)
    return sp.csr_matrix((data, (rows, cols)), shape=(d, d))


def _Z_power(d: int, b: int) -> sp.csr_matrix:
    """
    Sparse Z^b on a single qudit (dimension d).
    Diagonal with entries ω^{b*j}, ω = exp(2πi/d).
    """
    b %= d
    j = np.arange(d)
    omega = np.exp(2j * np.pi / d)
    diag = omega ** (b * j)
    return sp.csr_matrix(sp.diags(diag, offsets=0, dtype=complex, format="csr"))


def pauli_unitary_qudit(d: int, x: int, z: int, convention: str = "bare") -> sp.csr_matrix:
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
    d: int, x: np.ndarray, z: np.ndarray, convention: str = "bare"
) -> sp.csr_matrix:
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


def CX_func(i, a0, a1, dims):
    aa = int_to_bases(i, dims)
    aa[a1] = (aa[a1] + aa[a0]) % dims[a1]
    return bases_to_int(aa, dims)


def _mixed_radix_strides(dims: np.ndarray) -> np.ndarray:
    """
    strides[k] = product of dims[k+1:], with strides[-1] = 1.
    Convention: linear index i = sum_k digits[k] * strides[k],
    where digits[k] in [0, dims[k]-1], qudit 0 is most significant.
    """
    dims = np.asarray(dims, dtype=int)
    q = len(dims)
    strides = np.empty(q, dtype=int)
    strides[-1] = 1
    for k in range(q - 2, -1, -1):
        strides[k] = strides[k + 1] * dims[k + 1]
    return strides


def _int_to_digits(i: int, dims: np.ndarray, strides: np.ndarray) -> np.ndarray:
    """Decode linear index i to mixed-radix digits using the given strides."""
    # digits[k] = (i // strides[k]) % dims[k]
    return (i // strides) % dims


def _digits_to_int(digits: np.ndarray, strides: np.ndarray) -> int:
    """Encode mixed-radix digits to linear index using the given strides."""
    return int(np.dot(digits, strides))


def SWAP_func(i: int, a0: int, a1: int, dims: np.ndarray) -> int:
    """
    Map a basis index i -> f(i) by swapping qudit positions a0 <-> a1
    in a mixed-radix register with local dimensions `dims`.
    """
    dims = np.asarray(dims, dtype=int)
    q = len(dims)
    if not (0 <= a0 < q and 0 <= a1 < q and a0 != a1):
        raise ValueError("Invalid swap positions a0, a1.")
    strides = _mixed_radix_strides(dims)
    digits = _int_to_digits(i, dims, strides).astype(int)
    digits[a0], digits[a1] = digits[a1], digits[a0]
    return _digits_to_int(digits, strides)
