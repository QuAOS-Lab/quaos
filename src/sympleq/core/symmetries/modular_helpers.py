import numpy as np


def independent_columns(B: np.ndarray, p: int) -> np.ndarray:
    """Return a column-subset of B with independent columns over GF(p) (single RREF pass)."""
    if B.size == 0:
        return B
    _, piv = rref_mod(mod_p(B, p), p)
    piv = [pc for pc in piv if pc < B.shape[1]]
    return B[:, piv] if piv else np.zeros((B.shape[0], 0), dtype=np.int64)


def omega_matrix(n: int, p: int) -> np.ndarray:
    Id = np.eye(n, dtype=np.int64)
    zeros = np.zeros((n, n), dtype=np.int64)
    top = np.concatenate([zeros, Id], axis=1)
    bot = np.concatenate([mod_p(-Id, p), zeros], axis=1)
    return np.concatenate([top, bot], axis=0)


def is_symplectic(F: np.ndarray, p: int) -> bool:
    n2 = F.shape[0]
    assert n2 % 2 == 0 and F.shape[1] == n2
    Omega = omega_matrix(n2 // 2, p)
    return np.array_equal(mod_p(F.T @ Omega @ F, p), Omega % p)


def basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    """
    Canonical greedy basis extension over GF(p) (single shared implementation).

    Deterministically pick ``want`` columns from ``candidates`` that extend
    ``span(base)``, returning them as an (n x want) array. Raises if fewer than
    ``want`` independent extending columns exist.

    Incremental implementation: a running forward-echelon of span(base ∪ picked)
    is maintained, and each candidate is accepted iff reducing it against the
    echelon leaves a nonzero pivot (i.e. it strictly increases the rank). This is
    the same accept/reject decision as the previous per-candidate ``rank_mod``,
    so the selected columns are identical, but it avoids a full RREF per
    candidate.
    """
    n = candidates.shape[0]
    echelon: list[tuple[int, np.ndarray]] = []  # (pivot_pos, normalized row), pivot entry == 1

    def _reduce_and_add(v: np.ndarray) -> bool:
        v = mod_p(v, p).reshape(-1).copy()
        for pos, pv in echelon:
            if v[pos]:
                v = mod_p(v - int(v[pos]) * pv, p)
        nz = np.flatnonzero(v)
        if nz.size == 0:
            return False
        pos = int(nz[0])
        inv = inv_mod_scalar(int(v[pos]), p)
        echelon.append((pos, mod_p(v * inv, p)))
        return True

    if base.size:
        bred = independent_columns(mod_p(base, p), p)
        for j in range(bred.shape[1]):
            _reduce_and_add(bred[:, j])

    picked_cols: list[np.ndarray] = []
    for j in range(candidates.shape[1]):
        col = mod_p(candidates[:, j:j + 1], p)
        if _reduce_and_add(col[:, 0]):
            picked_cols.append(col)
            if len(picked_cols) == want:
                return np.concatenate(picked_cols, axis=1)
    raise RuntimeError("basis_extend: could not extend by required amount.")


def mod_p(A: np.ndarray, p: int) -> np.ndarray:
    A = np.asarray(A)
    if p == 2:
        return (A.astype(np.int64, copy=False) & 1)
    return np.asarray(A % p, dtype=np.int64)


def rref_mod(aug: np.ndarray, p: int) -> tuple[np.ndarray, list[int]]:
    if p == 2:
        return rref_mod2(aug)
    return rref_modp(aug, p)


def rref_mod2(aug: np.ndarray) -> tuple[np.ndarray, list[int]]:
    A = mod_p(aug, 2).astype(np.int64, copy=False)
    m, n = A.shape
    if n == 0 or m == 0 or n > 62:
        return _rref_mod2_dense(A)
    # Pack each row into one int64 (bit (n-1-j) = A[i, j]); GF(2) row reduction
    # is then scalar integer XOR, far cheaper than per-column numpy slicing for
    # the small widths used here. Pack/unpack are vectorized. Output is the
    # identical reduced matrix and pivot columns (RREF is canonical).
    weights = (1 << np.arange(n - 1, -1, -1, dtype=np.int64))
    packed = [int(x) for x in (A @ weights)]
    piv_cols: list[int] = []
    r = 0
    for c in range(n):
        if r >= m:
            break
        bit = 1 << (n - 1 - c)
        piv = -1
        for i in range(r, m):
            if packed[i] & bit:
                piv = i
                break
        if piv < 0:
            continue
        if piv != r:
            packed[r], packed[piv] = packed[piv], packed[r]
        pr = packed[r]
        for i in range(m):
            if i != r and (packed[i] & bit):
                packed[i] ^= pr
        piv_cols.append(c)
        r += 1
    parr = np.array(packed, dtype=np.int64).reshape(-1, 1)
    R = ((parr >> np.arange(n - 1, -1, -1, dtype=np.int64)) & 1).astype(np.int64)
    return R, piv_cols


def _rref_mod2_dense(A: np.ndarray) -> tuple[np.ndarray, list[int]]:
    A = mod_p(A.copy(), 2).astype(np.int64, copy=False)
    m, n = A.shape
    r = 0
    piv_cols: list[int] = []
    for c in range(n):
        if r >= m:
            break
        nz = np.flatnonzero(A[r:, c] & 1)
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        mask = (A[:, c] & 1).astype(bool)
        mask[r] = False
        if np.any(mask):
            A[mask, :] ^= A[r, :]
        piv_cols.append(c)
        r += 1
    return A, piv_cols


def rref_modp(aug: np.ndarray, p: int) -> tuple[np.ndarray, list[int]]:
    A = mod_p(aug.copy(), p)
    m, n = A.shape
    r = 0
    c = 0
    piv_cols: list[int] = []
    while r < m and c < n:
        piv = None
        for i in range(r, m):
            if A[i, c] % p != 0:
                piv = i
                break
        if piv is None:
            c += 1
            continue
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = inv_mod_scalar(A[r, c], p)
        A[r, :] = mod_p(A[r, :] * inv, p)
        # Eliminate column c from every other row in one vectorized step:
        # subtract fac_i * (pivot row) from row i, with fac_r forced to 0.
        facs = (A[:, c] % p).copy()
        facs[r] = 0
        if np.any(facs):
            A = mod_p(A - np.outer(facs, A[r, :]), p)
        piv_cols.append(c)
        r += 1
        c += 1
    return A, piv_cols


def rank_mod(A: np.ndarray, p: int) -> int:
    _, piv = rref_mod(mod_p(A, p), p)
    # safer than “count nonzero rows” when you start passing augmented matrices around
    return len([pc for pc in piv if pc < A.shape[1]])


def nullspace_mod(A: np.ndarray, p: int) -> np.ndarray:
    """Right nullspace basis of A over GF(p); columns form a basis."""
    A = mod_p(A, p)
    m, n = A.shape
    aug = np.concatenate([A, np.zeros((m, 1), dtype=np.int64)], axis=1)
    R, piv_cols = rref_mod(aug, p)
    piv_set = set(piv_cols)
    free = [j for j in range(n) if j not in piv_set]
    if not free:
        return np.zeros((n, 0), dtype=np.int64)
    # For the basis vector of free column f, the free part is the indicator e_f,
    # so each pivot coordinate is simply x[pc] = -R[row, f]. This is the closed
    # form of the old triple loop (which summed R[row, j]*x[j] over free j with
    # x = e_f), produced here in one vectorized assignment per pivot row.
    free_arr = np.asarray(free, dtype=np.int64)
    basis = np.zeros((n, free_arr.size), dtype=np.int64)
    basis[free_arr, np.arange(free_arr.size)] = 1
    pivs = [pc for pc in piv_cols if pc < n]
    for row_idx, pc in enumerate(pivs):
        basis[pc, :] = (-R[row_idx, free_arr]) % p
    return mod_p(basis, p)


def _solve_linear(A: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    """Solve A x = b over GF(p); returns one particular solution (free vars = 0)."""
    A = mod_p(A, p)
    b = mod_p(b.reshape(-1, 1), p)
    aug = np.concatenate([A, b], axis=1)
    R, piv_cols = rref_mod(aug, p)
    m, n = A.shape
    # Consistency
    for i in range(m):
        if np.all(R[i, :n] % p == 0) and (R[i, n] % p != 0):
            raise RuntimeError("No solution to linear system over GF(p)")
    x = np.zeros((n, 1), dtype=np.int64)
    row_idx = 0
    for pc in piv_cols:
        if pc < n:
            x[pc, 0] = R[row_idx, n] % p
            row_idx += 1
    return x


# Public name for cross-module use (the leading-underscore alias is retained for
# backward compatibility with existing imports).
solve_linear = _solve_linear


def solve_linear_many(A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """
    Solve A X = B over GF(p) for multiple RHS (columns of B) in ONE RREF.
    Returns one particular solution (free vars = 0).
    """
    A = mod_p(A, p)
    B = mod_p(B, p)
    if B.ndim == 1:
        B = B.reshape(-1, 1)
    m, n = A.shape
    if B.shape[0] != m:
        raise ValueError("solve_linear_many: incompatible shapes")

    aug = np.concatenate([A, B], axis=1)
    R, piv_cols = rref_mod(aug, p)

    # Consistency: rows with 0...0 | nonzero RHS
    left = R[:, :n]
    right = R[:, n:]
    bad = np.where(np.all(left % p == 0, axis=1) & np.any(right % p != 0, axis=1))[0]
    if bad.size:
        raise RuntimeError("No solution to linear system over GF(p)")

    X = np.zeros((n, B.shape[1]), dtype=np.int64)
    row = 0
    for pc in piv_cols:
        if pc < n:
            X[pc, :] = right[row, :] % p
            row += 1
    return mod_p(X, p)


def matmul_mod(A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    return mod_p(A @ B, p)


def inv_mod_scalar(a: int | np.integer, p: int) -> int:
    a = int(a) % p
    if a == 0:
        raise ZeroDivisionError("inv_mod_scalar: cannot invert 0 mod p")
    if p == 2:
        return 1  # only 1 is invertible
    return pow(a, p - 2, p)


def inv_mod_mat(A: np.ndarray, p: int) -> np.ndarray:
    """Gauss-Jordan inverse over GF(p). Raises if singular."""
    n = A.shape[0]
    aug = np.concatenate([mod_p(A, p), np.eye(n, dtype=np.int64)], axis=1)
    R, _ = rref_mod(aug, p)
    left = R[:, :n]
    right = R[:, n:]
    if not np.array_equal(left % p, np.eye(n, dtype=np.int64)):
        raise ValueError("Matrix not invertible mod p")
    return mod_p(right, p)



def solve_in_span(S: np.ndarray, A: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    """
    Find x in span(S) such that A x = b over GF(p).
    Inputs:
      S: (n x k) basis columns spanning the allowed subspace
      A: (m x n) constraint matrix
      b: (m,) or (m x 1)
    Returns:
      x: (n x 1) one solution in span(S)
    Raises if no solution exists.
    """
    S = mod_p(S, p)
    A = mod_p(A, p)
    b = mod_p(np.asarray(b).reshape(-1, 1), p)

    AS = matmul_mod(A, S, p)          # m × k
    y = _solve_linear(AS, b, p)       # k × 1
    x = matmul_mod(S, y, p)           # n × 1
    return x


def mat_pow_mod(A: np.ndarray, e: int, p: int) -> np.ndarray:
    """Compute A^e mod p by fast exponentiation."""
    if e < 0:
        raise ValueError("mat_pow_mod: e must be >= 0")
    A = mod_p(A, p)
    n = A.shape[0]
    R = np.eye(n, dtype=np.int64)
    B = A
    ee = int(e)
    while ee:
        if ee & 1:
            R = matmul_mod(R, B, p)
        ee >>= 1
        if ee:
            B = matmul_mod(B, B, p)
    return R


def column_rank_mod(B: np.ndarray, p: int) -> int:
    """Rank of the column span of B over GF(p)."""
    return rank_mod(B, p)

 
def solve_mod(mat: np.ndarray, vec: np.ndarray, modulus: int) -> np.ndarray:
    """Gaussian elimination over Z_mod, requiring unit pivots (gcd=1).
    Not necessarily prime_modulus.

    Returns a particular solution x (length = number of columns) with free
    variables set to 0. The solution is reconstructed from the recorded pivot
    columns, so it is correct even when pivots skip columns or when m < ncols
    (the previous ``aug[:ncols, -1]`` slice silently assumed diagonal pivots).
    """
    mat = mat.copy().astype(int)
    vec = vec.copy().astype(int)
    m, ncols = mat.shape
    aug = np.concatenate([mat, vec.reshape(-1, 1)], axis=1) % modulus
    row = 0
    pivot_cols: list[int] = []
    for col in range(ncols):
        if row >= m:
            break
        pivot = None
        for r in range(row, m):
            if np.gcd(int(aug[r, col]), modulus) == 1:
                pivot = r
                break
        if pivot is None:
            continue
        if pivot != row:
            aug[[row, pivot]] = aug[[pivot, row]]
        inv = pow(int(aug[row, col]) % modulus, -1, modulus)
        aug[row] = (aug[row] * inv) % modulus
        for r in range(m):
            if r == row:
                continue
            factor = aug[r, col]
            aug[r] = (aug[r] - factor * aug[row]) % modulus
        pivot_cols.append(col)
        row += 1

    # Reconstruct the solution aligned to pivot columns (free vars = 0).
    x = np.zeros(ncols, dtype=int)
    for r, col in enumerate(pivot_cols):
        x[col] = int(aug[r, -1]) % modulus
    return x % modulus
