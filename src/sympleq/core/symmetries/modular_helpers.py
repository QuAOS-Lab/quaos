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
    Ω = omega_matrix(n2 // 2, p)
    return np.array_equal(mod_p(F.T @ Ω @ F, p), Ω % p)


def mod_p(A: np.ndarray, p: int) -> np.ndarray:
    return np.asarray(A % p, dtype=np.int64)


def rank_mod(A: np.ndarray, p: int) -> int:
    R, _ = rref_mod(mod_p(A, p), p)
    # Count nonzero rows
    return int(np.sum(np.any(R % p != 0, axis=1)))


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
    basis = []
    for f in free:
        x = np.zeros((n, 1), dtype=np.int64)
        x[f, 0] = 1
        row_idx = 0
        for pc in piv_cols:
            if pc < n:
                s = 0
                for j in free:
                    s = (s + (R[row_idx, j] % p) * (x[j, 0] % p)) % p
                x[pc, 0] = (-s) % p
                row_idx += 1
        basis.append(x.reshape(-1))
    return np.stack(basis, axis=1)


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


def rref_mod(aug: np.ndarray, p: int) -> tuple[np.ndarray, list[int]]:
    """RREF over GF(p). Returns (RREF_augmented, pivot_cols)."""
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
        for i in range(m):
            if i != r and A[i, c] % p != 0:
                fac = A[i, c] % p
                A[i, :] = mod_p(A[i, :] - fac * A[r, :], p)
        piv_cols.append(c)
        r += 1
        c += 1
    return A, piv_cols


def matmul_mod(A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    return mod_p(A @ B, p)


def inv_mod_scalar(a: int | np.integer, p: int) -> int:
    return pow(int(a) % p, p - 2, p)


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
