"""Implements a random symplectic with the approach from """
from __future__ import annotations

import numpy as np

from sympleq.core.paulis._typing import TableauType


def symplectic_group_size(n: int, p: int = 2) -> int:
    """
    |Sp(2n, p)| = p^(n^2) * prod_{j=1..n} (p^(2j) - 1)
    Default p=2.
    """
    total = p ** (n * n)
    for j in range(1, n + 1):
        total *= (p ** (2 * j) - 1)
    return total


def direct_sum(m1, m2):
    n1, n2 = m1.shape[0], m2.shape[0]
    out = np.zeros((n1 + n2, n1 + n2), dtype=np.int8)
    out[:n1, :n1] = m1
    out[n1:, n1:] = m2
    return out


def int2bits(i: int, n: int) -> np.ndarray:
    """LSB-first length-n bit vector (dtype int8)."""
    out = np.zeros(n, dtype=np.int8)
    for j in range(n):
        out[j] = i & 1
        i >>= 1
    return out


# ----------------------------
# Symplectic arithmetic over GF(2)
# (interleaved ordering: [x0,z0,x1,z1,...])
# ----------------------------
def inner(v: np.ndarray, w: np.ndarray) -> int:
    """Symplectic inner product over GF(2) for interleaved ordering.
       For each pair (x,z): x_i * w_z_i + z_i * w_x_i (== subtraction in GF(2))."""
    nn = v.size
    assert nn == w.size and (nn % 2 == 0)
    t = 0
    for i in range(0, nn, 2):
        t ^= (int(v[i]) & int(w[i + 1]))
        t ^= (int(v[i + 1]) & int(w[i]))
    return t & 1


def transvection(k: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Z_k(v) = v + <k,v> k  (mod 2)."""
    coeff = inner(k, v)
    if coeff == 0:
        return v.copy()
    return (v + k) & 1


def find_transvection(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Find up to two transvection vectors h1,h2 (returned shape (2, len(x))) such that
    y = Z_{h1} Z_{h2} x. Follows the full-case logic (robust) from the paper.
    All arrays dtype=int8.
    """
    n = x.size
    out = np.zeros((2, n), dtype=np.int8)

    if np.array_equal(x, y):
        return out

    if inner(x, y) == 1:
        out[0] = (x + y) & 1
        return out

    # Try to find a single-block index i where both x and y have non-zero 2-block
    z = np.zeros(n, dtype=np.int8)
    for b in range(0, n // 2):
        ii = 2 * b
        xsum = (int(x[ii]) + int(x[ii + 1]))
        ysum = (int(y[ii]) + int(y[ii + 1]))
        if xsum != 0 and ysum != 0:
            # set z_block = x_block + y_block mod 2
            z[ii] = (x[ii] + y[ii]) & 1
            z[ii + 1] = (x[ii + 1] + y[ii + 1]) & 1
            # if z_block == 00, fix it as in the paper
            if (z[ii] + z[ii + 1]) == 0:
                z[ii + 1] = 1
                if x[ii] != x[ii + 1]:
                    z[ii] = 1
            out[0] = (x + z) & 1
            out[1] = (y + z) & 1
            return out

    # Otherwise find two blocks: one where x != 00 and y == 00, and one where x == 00 and y != 00
    # find block where x != 00 and y == 00
    for b in range(0, n // 2):
        ii = 2 * b
        if ((x[ii] + x[ii + 1]) != 0) and ((y[ii] + y[ii + 1]) == 0):
            if x[ii] == x[ii + 1]:
                z[ii + 1] = 1
            else:
                z[ii + 1] = x[ii]
                z[ii] = x[ii + 1]
            break

    # find block where x == 00 and y != 00
    for b in range(0, n // 2):
        ii = 2 * b
        if ((x[ii] + x[ii + 1]) == 0) and ((y[ii] + y[ii + 1]) != 0):
            if y[ii] == y[ii + 1]:
                z[ii + 1] = 1
            else:
                z[ii + 1] = y[ii]
                z[ii] = y[ii + 1]
            break

    out[0] = (x + z) & 1
    out[1] = (y + z) & 1
    return out


def symplectic_gf2_interleaved(index: int, n: int) -> TableauType:
    """
    Deterministic canonical enumeration of Sp(2n,2) per Koenig/Smolin appendix.
    Returns 2n x 2n numpy array dtype=int8 in INTERLEAVED ordering [x0,z0,x1,z1,...].
    index must be 0 <= index < symplectic_group_size(n,2).
    """
    if n <= 0:
        raise ValueError("n must be >= 1")

    total = symplectic_group_size(n, p=2)
    if index < 0 or index >= total:
        raise ValueError(f"index out of range: should be in [0, {total - 1}]")

    # Working copy of index that we peel off at each recursion level
    i = int(index)

    def _symplectic_recursive(i_local: int, n_local: int) -> np.ndarray:
        nn_local = 2 * n_local
        # Step 1: choose k in 1..(2^{nn}-1)
        s = (1 << nn_local) - 1
        k = (i_local % s) + 1
        i_new = i_local // s

        # Step 2: f1 is k as nn bits (LSB-first)
        f1 = int2bits(k, nn_local)

        # Step 3: find T mapping e1 -> f1
        e1 = np.zeros(nn_local, dtype=np.int8)
        e1[0] = 1
        T = find_transvection(e1, f1)  # shape (2, nn_local)

        # Step 4: read next nn_local-1 bits for e'
        bits = int2bits(i_new % (1 << (nn_local - 1)), nn_local - 1)
        i_new //= (1 << (nn_local - 1))

        # Step 5: construct e' (bits[0] used for step 6 later; bits[1:] fill higher coords)
        e_prime = e1.copy()
        # fill coordinates j = 2..nn_local-1 with bits[1..] (paper's indexing; using LSB-first)
        for j in range(2, nn_local):
            e_prime[j] = bits[j - 1]

        # Step 6: h0 = T(e')
        h0 = transvection(T[0], e_prime)
        h0 = transvection(T[1], h0)
        # if bits[0]==1 then h0 = h0 + f1  (this is the correct GF(2) action)
        if bits[0] == 1:
            h0 = (h0 + f1) & 1

        # Step 7: recursive call for remaining block
        id2 = np.zeros((2, 2), dtype=np.int8)
        id2[0, 0] = 1
        id2[1, 1] = 1

        if n_local > 1:
            g_small = _symplectic_recursive(i_new, n_local - 1)
            g = direct_sum(id2, g_small)
        else:
            g = id2.copy()

        # Apply transvections (left multiplication) by transforming columns
        for col_idx in range(nn_local):
            col = g[:, col_idx].copy()
            col = transvection(T[0], col)
            col = transvection(T[1], col)
            col = transvection(h0, col)
            col = transvection(f1, col)
            g[:, col_idx] = col

        return g

    return _symplectic_recursive(i, n)


def symplectic_gf2(index: int, n: int) -> TableauType:
    return interleaved_to_grouped(symplectic_gf2_interleaved(index, n))


def interleaved_to_grouped(F_inter: np.ndarray) -> np.ndarray:
    """
    Convert interleaved [x0,z0,x1,z1,...] => grouped [x0,x1,...,z0,z1,...]
    by similarity transform F_group = S^T F_inter S (implemented by selecting rows/cols
    with the inverse permutation).
    """
    n2 = F_inter.shape[0]
    assert n2 % 2 == 0
    n = n2 // 2

    # interleaved -> grouped mapping
    i2g = [i // 2 + (i % 2) * n for i in range(n2)]
    # invert
    perm_inv = [0] * n2
    for i, g in enumerate(i2g):
        perm_inv[g] = i

    # similarity: S^T @ F @ S  <=> select rows perm_inv, cols perm_inv
    return F_inter[np.ix_(perm_inv, perm_inv)].astype(np.int8)


def is_symplectic_interleaved(F: np.ndarray) -> bool:
    """Check symplectic in interleaved ordering (Ω = diag(J,J,...), J=[[0,1],[1,0]])."""
    nn = F.shape[0]
    assert nn % 2 == 0
    Omega = np.zeros((nn, nn), dtype=np.int8)
    for i in range(nn // 2):
        ii = 2 * i
        Omega[ii, ii + 1] = 1
        Omega[ii + 1, ii] = 1
    lhs = (F.T @ Omega @ F) & 1
    return np.array_equal(lhs, Omega)


def _isotropic_vector(n, d):
    """
    Sample an isotropic vector v = (a|b) in Z_d^{2n}.
    For d=2 (qubits), every vector is isotropic.
    For prime d, ensures <v,v> = 0 mod d.
    """
    if d == 2:
        # Any vector works
        return np.random.randint(0, 2, size=(2 * n,), dtype=int)

    # d prime
    while True:
        a = np.random.randint(0, d, size=(n,), dtype=int)
        if np.all(a == 0):
            continue  # avoid trivial a
        # Find random b orthogonal to a
        while True:
            b = np.random.randint(0, d, size=(n,), dtype=int)
            if (a @ b) % d == 0:
                return np.concatenate([a, b])


def _vector_to_transvection(v, J, d):
    """
    Return the symplectic transvection matrix T_v over Z_d:
        T_v(w) = w + <v, w> * v (mod d).
    """
    v = v.reshape(-1, 1)
    return (np.identity(len(v), dtype=int) + (J @ v) @ v.T) % d


def symplectic_random_transvection(n_qudits: int, dimension: int = 2,
                                   num_transvections: int | None = None) -> TableauType:
    """
    Return a random 2n x 2n symplectic matrix over Z_d by composing
    num_transvections random transvections.

    Parameters
    ----------
    n_qudits : int
        Number of qudits (i.e. pairs of rows/cols).
    dimension : int
        Dimension of each qudit (>=2).
    num_transvections : int or None
        Number of transvections to compose. If None, defaults to 2 * (2n).

    Returns
    -------
    M : (2n x 2n) integer matrix
        Random symplectic matrix over Z_d.
    """
    Id_n = np.identity(n_qudits, dtype=int)
    Zero_n = np.zeros((n_qudits, n_qudits), dtype=int)
    J = np.block([[Zero_n, Id_n], [-Id_n, Zero_n]]) % dimension

    dim = 2 * n_qudits
    M = np.identity(dim, dtype=int)

    if num_transvections is None:
        num_transvections = 2 * dim

    for _ in range(num_transvections):
        v = _isotropic_vector(n_qudits, dimension)
        Mv = _vector_to_transvection(v, J, dimension)
        M = (M @ Mv) % dimension

    return M
