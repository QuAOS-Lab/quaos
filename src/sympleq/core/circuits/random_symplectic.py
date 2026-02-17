"""Implements a random symplectic with the approach from """

import numpy as np


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
# (grouped ordering: [x0,x1,...,z0,z1,...])
# ----------------------------
def inner(v: np.ndarray, w: np.ndarray) -> int:
    """Symplectic inner product over GF(2) in grouped ordering."""
    n2 = v.size
    assert n2 == w.size and (n2 % 2 == 0)
    n = n2 // 2
    x_dot_zw = int(np.dot(v[:n].astype(np.int64), w[n:].astype(np.int64)))
    z_dot_xw = int(np.dot(v[n:].astype(np.int64), w[:n].astype(np.int64)))
    return (x_dot_zw + z_dot_xw) & 1


def transvection(k: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Z_k(v) = v + <k,v> k  (mod 2)."""
    coeff = inner(k, v)
    if coeff == 0:
        return v.copy()
    return (v + k) & 1


def find_transvection(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Find up to two transvection vectors h1,h2 (returned shape (2, len(x))) such that
    y = Z_{h1} Z_{h2} x.
    All arrays dtype=int8.
    """
    n2 = x.size
    n = n2 // 2
    out = np.zeros((2, n2), dtype=np.int8)

    if np.array_equal(x, y):
        return out

    if inner(x, y) == 1:
        out[0] = (x + y) & 1
        return out

    # Try to find a single-block index i where both x and y have non-zero 2-block
    z = np.zeros(n2, dtype=np.int8)
    for q in range(n):
        xi = q
        zi = q + n
        xsum = int(x[xi]) + int(x[zi])
        ysum = int(y[xi]) + int(y[zi])
        if xsum != 0 and ysum != 0:
            # set z_block = x_block + y_block mod 2
            z[xi] = (x[xi] + y[xi]) & 1
            z[zi] = (x[zi] + y[zi]) & 1
            # if z_block == 00, fix it as in the paper
            if (z[xi] + z[zi]) == 0:
                z[zi] = 1
                if x[xi] != x[zi]:
                    z[xi] = 1
            out[0] = (x + z) & 1
            out[1] = (y + z) & 1
            return out

    # Otherwise find two blocks: one where x != 00 and y == 00, and one where x == 00 and y != 00
    # find block where x != 00 and y == 00
    for q in range(n):
        xi = q
        zi = q + n
        if ((x[xi] + x[zi]) != 0) and ((y[xi] + y[zi]) == 0):
            if x[xi] == x[zi]:
                z[zi] = 1
            else:
                z[zi] = x[xi]
                z[xi] = x[zi]
            break

    # find block where x == 00 and y != 00
    for q in range(n):
        xi = q
        zi = q + n
        if ((x[xi] + x[zi]) == 0) and ((y[xi] + y[zi]) != 0):
            if y[xi] == y[zi]:
                z[zi] = 1
            else:
                z[zi] = y[xi]
                z[xi] = y[zi]
            break

    out[0] = (x + z) & 1
    out[1] = (y + z) & 1
    return out


def symplectic_gf2(index: int, n: int) -> np.ndarray:
    """
    Deterministic canonical enumeration of Sp(2n,2) per Koenig/Smolin appendix.
    Returns 2n x 2n numpy array dtype=int8 in grouped ordering [x0,x1,...,z0,z1,...].
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
        # In grouped ordering, x0 is at 0 and z0 is at n_local; exclude both.
        free_coords = [idx for idx in range(nn_local) if idx not in (0, n_local)]
        for bit_idx, coord in enumerate(free_coords, start=1):
            e_prime[coord] = bits[bit_idx]

        # Step 6: h0 = T(e')
        h0 = transvection(T[0], e_prime)
        h0 = transvection(T[1], h0)
        # if bits[0]==1 then h0 = h0 + f1  (this is the correct GF(2) action)
        if bits[0] == 1:
            h0 = (h0 + f1) & 1

        # Step 7: recursive call for remaining block
        if n_local > 1:
            g_small = _symplectic_recursive(i_new, n_local - 1)
            g = np.eye(nn_local, dtype=np.int8)
            rest = list(range(1, n_local)) + list(range(n_local + 1, nn_local))
            g[np.ix_(rest, rest)] = g_small
        else:
            g = np.eye(nn_local, dtype=np.int8)

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


def symplectic_gf2_interleaved(index: int, n: int) -> np.ndarray:
    """
    Backward-compatible alias.

    Despite the historical name, this returns grouped ordering:
    [x0,x1,...,z0,z1,...].
    """
    return symplectic_gf2(index, n)


def interleaved_to_grouped(F_inter: np.ndarray) -> np.ndarray:
    """
    Backward-compatible no-op for grouped-order matrices.

    This module now uses grouped ordering everywhere:
    [x0,x1,...,z0,z1,...].
    """
    if F_inter.ndim != 2 or F_inter.shape[0] != F_inter.shape[1]:
        raise ValueError("Input must be a square matrix.")
    if F_inter.shape[0] % 2 != 0:
        raise ValueError("Matrix size must be even (2n x 2n).")
    return np.asarray(F_inter, dtype=np.int8)


def is_symplectic_interleaved(F: np.ndarray) -> bool:
    """Check symplectic in grouped ordering (legacy function name)."""
    nn = F.shape[0]
    assert nn % 2 == 0
    n = nn // 2
    I = np.eye(n, dtype=np.int8)
    Z = np.zeros((n, n), dtype=np.int8)
    Omega = np.block([[Z, I], [I, Z]])
    lhs = (F.T @ Omega @ F) & 1
    return np.array_equal(lhs, Omega)


def random_isotropic_vector(n, d, rng=None):
    """
    Sample an isotropic vector v = (a|b) in Z_d^{2n}.
    For d=2 (qubits), every vector is isotropic.
    For prime d, ensures <v,v> = 0 mod d.
    """
    if rng is None:
        rng = np.random.default_rng()
    if d == 2:
        # Any vector works
        return rng.integers(0, 2, size=(2 * n,), dtype=int)

    # d prime
    while True:
        a = rng.integers(0, d, size=(n,), dtype=int)
        if np.all(a == 0):
            continue  # avoid trivial a
        # Find random b orthogonal to a
        while True:
            b = rng.integers(0, d, size=(n,), dtype=int)
            if (a @ b) % d == 0:
                return np.concatenate([a, b])


def _vector_to_transvection(v, J, d):
    """
    Return the symplectic transvection matrix T_v over Z_d:
        T_v(w) = w + <v, w> * v (mod d).
    """
    v = v.reshape(-1, 1)
    return (np.identity(len(v), dtype=int) + (J @ v) @ v.T) % d


def symplectic_random_transvection(n_qudits, dimension=2, num_transvections=None, rng=None):
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

    if rng is None:
        rng = np.random.default_rng()

    Id_n = np.identity(n_qudits, dtype=int)
    Zero_n = np.zeros((n_qudits, n_qudits), dtype=int)
    J = np.block([[Zero_n, Id_n], [-Id_n, Zero_n]]) % dimension

    dim = 2 * n_qudits
    M = np.identity(dim, dtype=int)

    if num_transvections is None:
        num_transvections = 2 * dim
    for _ in range(num_transvections):
        v = random_isotropic_vector(n_qudits, dimension, rng=rng)
        Mv = _vector_to_transvection(v, J, dimension)
        M = (M @ Mv) % dimension

    return M


def random_symplectic(n_qubits: int, dimension: int = 2, num_transvections: int | None = None) -> np.ndarray:
    return symplectic_random_transvection(n_qubits, dimension, num_transvections)
