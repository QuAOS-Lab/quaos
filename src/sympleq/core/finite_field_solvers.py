import numpy as np
from math import gcd
import galois
from collections import defaultdict
# TODO: Move tests to test suite (and write more!)


def solve_linear_system_over_gf(A: np.ndarray, b: np.ndarray, GF: type | int) -> np.ndarray:
    """
    Solve the system A @ x = b over a finite field.
    Returns one particular solution with free variables set to zero.
    Raises ValueError if the system is inconsistent.
    """
    if isinstance(GF, int):
        GF = galois.GF(GF)

    A_gf = GF(A)
    b_gf = GF(b).reshape(-1, 1)

    if A_gf.ndim != 2:
        raise ValueError("Coefficient matrix must be 2-dimensional.")
    if A_gf.shape[0] != b_gf.shape[0]:
        raise ValueError(
            f"Incompatible shapes for A ({A_gf.shape}) and b ({b_gf.shape})."
        )

    m, n = A_gf.shape
    augmented = np.hstack((A_gf, b_gf))
    R = augmented.copy()

    pivot_rows: list[tuple[int, int]] = []
    row = 0
    for col in range(n):
        pivot = None
        for r in range(row, m):
            if R[r, col] != GF(0):
                pivot = r
                break
        if pivot is None:
            continue

        if pivot != row:
            R[[row, pivot], :] = R[[pivot, row], :]

        pivot_val = R[row, col]
        R[row, :] /= pivot_val

        for r in range(m):
            if r != row and R[r, col] != GF(0):
                R[r, :] -= R[r, col] * R[row, :]

        pivot_rows.append((row, col))
        row += 1
        if row == m:
            break

    for r in range(m):
        if all(R[r, c] == GF(0) for c in range(n)) and R[r, n] != GF(0):
            raise ValueError("Inconsistent linear system over GF(p).")

    x = GF.Zeros(n)
    for row_idx, pivot_col in reversed(pivot_rows):
        value = R[row_idx, n]
        for j in range(pivot_col + 1, n):
            if R[row_idx, j] != GF(0):
                value -= R[row_idx, j] * x[j]
        x[pivot_col] = value

    return x.copy()


def solve_gf2(A: np.ndarray, b: np.ndarray) -> np.ndarray | None:
    """
    Solve the linear system Ax = b over GF(2) using Gaussian elimination.

    Args:
        A: Coefficient matrix
        b: Right-hand side vector

    Returns:
        Solution vector or None if no solution exists
    """
    GF2 = galois.GF(2)
    try:
        solution = solve_linear_system_over_gf(
            np.asarray(A, dtype=int) % 2, np.asarray(b, dtype=int) % 2, GF2
        )
    except ValueError:
        return None

    return solution.view(np.ndarray).astype(int, copy=False)


def solve_modular_linear_additive(x: int, z: int, d: int) -> int:
    """
    Find smallest non-negative integer n such that (x + z*n) % d == 0.
    Returns None if no solution exists
    """
    x, z, d = int(x), int(z), int(d)
    g = gcd(z, d)
    if x % g != 0:
        raise ValueError(f"No solution of (x + z*n) % d == 0 exists for x ={x}, z ={z}, d ={d}")

    # Reduce the equation modulo d // g
    z_ = z // g
    d_ = d // g
    x_ = (-x // g) % d_

    # compute modular inverse
    try:
        z_inv = pow(z_, -1, d_)
    except ValueError:
        raise ValueError(f"No solution of (x + z*n) % d == 0 exists for x ={x}, z ={z}, d ={d}")

    n = (x_ * z_inv) % d_
    return n


def solve_modular_linear_system(B, v):
    """
    Solve x @ B = v over GF(p) using row-reduction (RREF)
    """
    GF = type(B)
    solution = solve_linear_system_over_gf(B.T, v, GF)
    if not np.array_equal(solution @ B, v):
        raise ValueError("Failed to solve linear system over GF(p).")
    return solution


def gf_solve(A: np.ndarray, b: np.ndarray, GF: type) -> np.ndarray:
    """
    Solve A x = b over GF(p). Returns one particular solution or raises ValueError if inconsistent.
    """
    solution = solve_linear_system_over_gf(A, b, GF)
    return solution.reshape(-1, 1)


def get_linear_dependencies(vectors: np.ndarray,
                            p: int | list[int] | np.ndarray
                            ) -> tuple[list[int], dict[int, list[tuple[int, int]]]]:
    """
    Find linear dependencies among rows of `vectors` over finite fields.

    Parameters
    ----------
    vectors : np.ndarray
        Matrix of shape (m, n). Each row is a vector.
    p : int | list[int] | np.ndarray
        - If int: all columns are over GF(p).
        - If list/array of length n: p[j] is the prime for column j.
        - If list/array of length n//2: per-qudit primes; each value is duplicated for X/Z columns.
        - If list/array of length m: p[i] is the prime field per row (legacy mode).

    Returns
    -------
    pivot_indices : list[int]
        Indices of linearly independent rows.
    dependencies : dict[int, list[tuple[int, int]]]
        Mapping: row index -> list of (pivot_index, coefficient) such that
        row = sum(coeff * pivot).
    """
    m, n = vectors.shape  # m = rows (vectors), n = columns

    # Normalize p into per-column primes when possible.
    p_cols: list[int]
    mode: str
    if isinstance(p, int):
        p_cols = [int(p)] * n
        mode = "per-column"
    elif isinstance(p, (list, np.ndarray)):
        lp = len(p)
        if lp == n:
            p_cols = [int(x) for x in p]
            mode = "per-column"
        elif lp == n // 2 and n % 2 == 0:
            # Per-qudit list provided; duplicate for X/Z halves
            p_cols = [int(x) for x in np.repeat(p, 2)]
            mode = "per-column"
        elif lp == m:
            # Legacy behavior: per-row primes
            ps = [int(x) for x in p]
            mode = "per-row"
        else:
            raise AssertionError(
                f"Length of p must be either rows={m}, cols={n}, or cols/2={n // 2} (for qudits). Got {lp}"
            )
    else:
        raise TypeError(f"p must be int or list/np.ndarray of ints, it is {type(p)}")

    pivot_indices: list[int] = []
    dependencies: dict[int, list[tuple[int, int]]] = {}

    if mode == "per-row":
        # Group rows by prime and process each group independently
        prime_groups = defaultdict(list)
        for i, prime in enumerate(ps):
            prime_groups[prime].append(i)

        for prime, indices in prime_groups.items():
            if len(indices) == 1:
                pivot_indices.append(indices[0])
                continue

            GF = galois.GF(prime)
            V = GF(vectors[indices, :])

            # Identify group pivots incrementally
            group_pivots = []
            seen = GF.Zeros((0, V.shape[1]))
            for idx_in_group, i in enumerate(indices):
                candidate = V[idx_in_group]
                A = GF(np.vstack([seen, candidate]))
                R = A.row_reduce()
                rank = sum(1 for row in R if np.any(row != 0))
                if rank > seen.shape[0]:
                    group_pivots.append(i)
                    seen = A
            pivot_indices.extend(group_pivots)

            # Dependencies within the group
            if len(group_pivots) > 0:
                B = V[[indices.index(j) for j in group_pivots], :]
                for idx_in_group, i in enumerate(indices):
                    if i in group_pivots:
                        continue
                    x = solve_modular_linear_system(B, V[idx_in_group])
                    deps = [(group_pivots[j], int(x[j])) for j in range(len(x)) if x[j] != 0]
                    dependencies[i] = deps

        return pivot_indices, dependencies

    # Per-column primes path
    # Build prime -> column indices map
    prime_cols: dict[int, list[int]] = defaultdict(list)
    for j, prime in enumerate(p_cols):
        prime_cols[int(prime)].append(j)

    unique_primes = list(prime_cols.keys())

    # Track seen rows per-prime on the respective column subsets
    seen_per_prime: dict[int, np.ndarray] = {}
    for q in unique_primes:
        GFq = galois.GF(q)
        seen_per_prime[q] = GFq.Zeros((0, len(prime_cols[q])))

    # Select pivots: a row is independent if it increases rank for any prime on its column subset
    for i in range(m):
        increases_any = False
        updated_seen = {}
        for q in unique_primes:
            cols = prime_cols[q]
            if len(cols) == 0:
                continue
            GFq = galois.GF(q)
            candidate = GFq(vectors[i, cols])
            A_prev = seen_per_prime[q]
            A = GFq(np.vstack([A_prev, candidate]))
            R = A.row_reduce()
            rank_new = sum(1 for row in R if np.any(row != 0))
            if rank_new > A_prev.shape[0]:
                increases_any = True
                updated_seen[q] = A
            else:
                updated_seen[q] = A_prev

        if increases_any:
            pivot_indices.append(i)
            for q in unique_primes:
                seen_per_prime[q] = updated_seen[q]
        else:
            # Mark as dependent; coefficients are computed later as best-effort
            continue

    # Compute dependencies best-effort:
    # If all columns share the same prime, compute exact coefficients as before.
    if len(unique_primes) == 1:
        q = unique_primes[0]
        GFq = galois.GF(q)
        cols = prime_cols[q]
        Vq = GFq(vectors[:, cols])
        Bq = Vq[pivot_indices, :]
        for i in range(m):
            if i in pivot_indices:
                continue
            x = solve_modular_linear_system(Bq, Vq[i])
            deps = [(pivot_indices[j], int(x[j])) for j in range(len(x)) if x[j] != 0]
            dependencies[i] = deps
    else:
        # Mixed primes: return independent set; dependency coefficients are non-unique across fields.
        dependencies = {}

    return pivot_indices, dependencies


def gf_inv(A, p: int = 2):
    """
    Compute the inverse of a square matrix over GF(p) for a prime p.
    Defaults to GF(2) for backwards compatibility.
    """
    A = np.asarray(A, dtype=int) % p
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("Matrix must be square to compute an inverse over GF(p).")

    n = A.shape[0]
    Id = np.eye(n, dtype=int) % p
    AI = np.concatenate((A.copy(), Id), axis=1)

    for i in range(n):
        pivot_row = None
        for r in range(i, n):
            if AI[r, i] % p != 0:
                pivot_row = r
                break
        if pivot_row is None:
            raise ValueError("Matrix is singular over GF(p); inverse does not exist.")
        if pivot_row != i:
            AI[[i, pivot_row]] = AI[[pivot_row, i]]

        pivot_val = int(AI[i, i] % p)
        pivot_inv = pow(pivot_val, -1, p)
        AI[i, :] = (AI[i, :] * pivot_inv) % p

        for r in range(n):
            if r == i:
                continue
            factor = AI[r, i] % p
            if factor != 0:
                AI[r, :] = (AI[r, :] - factor * AI[i, :]) % p

    return AI[:, n:] % p


def gf_rref(A, p: int = 2):
    """
    Compute the reduced row echelon form of a matrix over GF(p) for a prime p.
    Returns the transformed matrix, the left transformation matrix M (row ops),
    the right transformation matrix N (column ops), and the rank.
    Defaults to GF(2) for backwards compatibility.
    """
    A = np.asarray(A, dtype=int) % p
    if A.ndim != 2:
        raise ValueError("Input matrix must be 2-dimensional for RREF.")

    A = A.copy()
    m, n = A.shape
    i = j = 0
    M = np.eye(m, dtype=int) % p
    N = np.eye(n, dtype=int) % p

    while i < m and j < n:
        pivot_row = None
        for r in range(i, m):
            if A[r, j] % p != 0:
                pivot_row = r
                break

        if pivot_row is None:
            j += 1
            continue

        if pivot_row != i:
            A[[i, pivot_row]] = A[[pivot_row, i]]
            M[[i, pivot_row]] = M[[pivot_row, i]]

        pivot_val = int(A[i, j] % p)
        pivot_inv = pow(pivot_val, -1, p)
        A[i, :] = (A[i, :] * pivot_inv) % p
        M[i, :] = (M[i, :] * pivot_inv) % p

        for r in range(m):
            if r == i:
                continue
            factor = A[r, j] % p
            if factor != 0:
                A[r, :] = (A[r, :] - factor * A[i, :]) % p
                M[r, :] = (M[r, :] - factor * M[i, :]) % p

        for c in range(n):
            if c == j:
                continue
            factor = A[i, c] % p
            if factor != 0:
                A[:, c] = (A[:, c] - factor * A[:, j]) % p
                N[:, c] = (N[:, c] - factor * N[:, j]) % p

        i += 1
        j += 1

    rank = i
    return A % p, M % p, N % p, rank


def gf_lu(A, p: int = 2):
    """
    Perform LU decomposition of a matrix over GF(p) for prime p.
    Returns L, U, P such that P @ A = L @ U with unit diagonal L.
    Defaults to GF(2) for backwards compatibility.
    """
    A = np.asarray(A, dtype=int) % p
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("LU decomposition requires a square matrix over GF(p).")

    m = A.shape[0]
    U = A.copy()
    L = np.eye(m, dtype=int) % p
    P = np.eye(m, dtype=int)

    for k in range(m):
        pivot_candidates = np.where(U[k:, k] % p != 0)[0]
        if pivot_candidates.size == 0:
            continue
        pivot = pivot_candidates[0] + k
        if pivot != k:
            U[[k, pivot], k:] = U[[pivot, k], k:]
            L[[k, pivot], :k] = L[[pivot, k], :k]
            P[[k, pivot]] = P[[pivot, k]]

        pivot_val = int(U[k, k] % p)
        pivot_inv = pow(pivot_val, -1, p)

        for j in range(k + 1, m):
            factor = (U[j, k] * pivot_inv) % p
            L[j, k] = factor
            if factor != 0:
                U[j, k:] = (U[j, k:] - factor * U[k, k:]) % p

    return L % p, U % p, P


def _select_row_basis_indices(A_int: np.ndarray, p: int, max_rows: int) -> np.ndarray:
    """
    Return indices of a greedily chosen row basis of `A_int` over GF(p).

    The routine performs a light Gauss-Jordan sweep, moving left-to-right across
    columns and picking the first unused row with a non-zero entry as the pivot.
    Each pivot row is normalized and used to eliminate the pivot column from the
    other unused rows. Collection stops once either all columns are processed or
    `max_rows` independent rows have been gathered.

    Parameters
    ----------
    A_int : np.ndarray
        Integer matrix whose rows are candidate equations.
    p : int
        Prime for the finite field GF(p).
        max_rows : int
        Maximum number of independent rows to return (often the number of columns).

    Returns
    -------
    np.ndarray
        1-D array of row indices that are linearly independent over GF(p).
    """
    GF = galois.GF(p)
    A = GF(A_int % p).copy()
    N, M = A.shape
    used = np.zeros(N, dtype=bool)
    basis = []
    col = 0
    for _ in range(N):
        if col >= M:
            break
        pick = None
        for r in range(N):
            if used[r]:
                continue
            if A[r, col] != GF(0):
                pick = r
                break
        if pick is None:
            col += 1
            continue
        basis.append(pick)
        used[pick] = True
        inv = GF(1) / A[pick, col]
        A[pick, :] = A[pick, :] * inv
        for r in range(N):
            if r == pick or used[r]:
                continue
            if A[r, col] != GF(0):
                A[r, :] = A[r, :] - A[r, col] * A[pick, :]
        col += 1
        if len(basis) >= max_rows:
            break
    return np.array(basis, dtype=int)


def _random_invertible_matrix(p: int, size: int, rng: np.random.Generator) -> np.ndarray:
    """Generate a random invertible matrix over GF(p) with the given dimension."""
    GFp = galois.GF(p)
    while True:
        mat = GFp(rng.integers(0, p, size=(size, size)))
        if np.linalg.matrix_rank(mat) == size:
            return np.asarray(mat, dtype=int) % p


def _random_matrix(p: int, shape: tuple[int, int], rng: np.random.Generator) -> np.ndarray:
    """Generate a random matrix over GF(p) with the given shape."""
    return rng.integers(0, p, size=shape, dtype=int) % p


def _is_permutation_matrix(P: np.ndarray) -> bool:
    """Check whether a matrix is a permutation matrix."""
    if P.ndim != 2 or P.shape[0] != P.shape[1]:
        return False
    return bool(np.all((P == 0) | (P == 1)) and np.all(P.sum(axis=0) == 1) and np.all(P.sum(axis=1) == 1))


def _test_gf_inv() -> None:
    rng = np.random.default_rng(1234)
    for p in (2, 3, 5, 7):
        for n in (1, 2, 4):
            A = _random_invertible_matrix(p, n, rng)
            inv = gf_inv(A, p=p)
            prod = (A @ inv) % p
            assert np.array_equal(prod, np.eye(n, dtype=int) % p), f"Inverse failed for p={p}, n={n}"


def _test_gf_rref() -> None:
    rng = np.random.default_rng(5678)
    for p in (2, 3, 5):
        m, n = 4, 6
        for _ in range(5):
            A = _random_matrix(p, (m, n), rng)
            R, M, N, rank = gf_rref(A, p=p)
            left = (M @ A) % p
            recon = (left @ N) % p
            assert np.array_equal(recon, R), f"Reconstruction failed for p={p}"
            pivots = []
            for row_idx in range(m):
                row = R[row_idx]
                nz = np.nonzero(row)[0]
                if nz.size == 0:
                    assert np.all(row % p == 0)
                    continue
                pivot_col = nz[0]
                pivots.append(pivot_col)
                assert row[pivot_col] % p == 1
                assert np.all(row[:pivot_col] % p == 0)
                assert np.all(row[pivot_col + 1:] % p == 0)
            assert rank == len(pivots), f"Rank mismatch for p={p}"


def _test_gflu() -> None:
    rng = np.random.default_rng(91011)
    for p in (2, 5, 11):
        for n in (2, 3, 5):
            A = _random_matrix(p, (n, n), rng)
            L, U, P = gf_lu(A, p=p)
            assert _is_permutation_matrix(P), f"P is not a permutation matrix for p={p}"
            PA = (P @ A) % p
            LU = (L @ U) % p
            assert np.array_equal(PA, LU), f"LU factorization failed for p={p}"
            assert np.array_equal(np.diag(L) % p, np.ones(n, dtype=int)), f"Diagonal of L not unit for p={p}"
            assert np.all((np.triu(L, k=1) % p) == 0), f"L not lower-triangular for p={p}"


if __name__ == "__main__":
    _test_gf_inv()
    _test_gf_rref()
    _test_gflu()
    print("All finite field solver self-tests passed.")
