
from typing import Optional
import numpy as np
from math import gcd
import galois
from collections import defaultdict


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


def get_linear_dependencies(
    vectors: np.ndarray,
    p: int | list[int] | np.ndarray,
    compute_dependencies: bool = True,
) -> tuple[list[int], dict[int, list[tuple[int, int]]]]:
    """
    Fast replacement for get_linear_dependencies.

    - For single prime p (especially p=2): returns exact pivot rows and exact dependencies.
    - For per-column primes: returns correct pivot rows (same criterion as your code),
      but dependency coefficients only if all primes are identical (same as your code).
    - For per-row primes: processes each group with the fast single-prime path.

    Notes:
    - For p=2, dependencies are coefficients in GF(2) (0/1).
    - For odd prime p, dependencies are coefficients in GF(p).
    """
    V = np.asarray(vectors, dtype=np.int64)
    m, n = V.shape

    # Normalize p into mode and per-column primes when possible.
    if isinstance(p, int):
        p_int = int(p)
        return _get_deps_single_prime(V, p_int, compute_dependencies)

    if isinstance(p, (list, np.ndarray)):
        p = np.asarray(p, dtype=int)
        lp = len(p)

        if lp == m:
            # per-row legacy: group and process each group
            pivot_indices: list[int] = []
            dependencies: dict[int, list[tuple[int, int]]] = {}
            prime_groups = defaultdict(list)
            for i, prime in enumerate(p.tolist()):
                prime_groups[int(prime)].append(i)
            for prime, idxs in prime_groups.items():
                piv, deps = _get_deps_single_prime(V[idxs, :], int(prime), compute_dependencies)
                pivot_indices.extend([idxs[j] for j in piv])
                if compute_dependencies:
                    for local_row, expr in deps.items():
                        dependencies[idxs[local_row]] = [(idxs[piv_j], c) for piv_j, c in expr]
            return pivot_indices, dependencies

        if lp == n:
            p_cols = p.tolist()
        elif lp == n // 2 and n % 2 == 0:
            p_cols = np.repeat(p, 2).tolist()
        else:
            raise AssertionError(f"Length of p must be rows={m}, cols={n}, or cols/2={n // 2}. Got {lp}")

        # Per-column primes: pivot criterion is “independent if increases rank for any prime subset”.
        # We can do this by running elimination separately per distinct prime on its column subset,
        # but WITHOUT galois row_reduce. Just incremental elimination mod p.

        prime_cols = defaultdict(list)
        for j, prime in enumerate(p_cols):
            prime_cols[int(prime)].append(j)
        primes = list(prime_cols.keys())

        # Maintain separate eliminators per prime
        eliminators = {q: _IncrementalElim(mod=q, ncols=len(prime_cols[q])) for q in primes}

        pivot_indices: list[int] = []
        # We can only produce dependencies if there's a single prime overall (same as your code)
        dependencies: dict[int, list[tuple[int, int]]] = {} if compute_dependencies else {}

        # We also need a mapping from global pivot list to each prime’s pivot indices if we want deps;
        # but we only do deps when single prime.
        for i in range(m):
            independent = False
            # test rank increase on any prime subset
            for q in primes:
                cols = prime_cols[q]
                row = V[i, cols] % q
                if eliminators[q].would_increase_rank(row):
                    independent = True
                    break

            if independent:
                pivot_indices.append(i)
                # actually add to all eliminators (mirrors your “updated_seen” logic)
                for q in primes:
                    cols = prime_cols[q]
                    row = V[i, cols] % q
                    eliminators[q].add_row(row, i, track_combo=False)

        # dependencies only if exactly one prime
        if compute_dependencies and len(primes) == 1:
            q = primes[0]
            cols = prime_cols[q]
            piv, deps = _get_deps_single_prime(V[:, cols], q, compute_dependencies=True)
            # piv returned are *row indices* already, but might differ slightly from the “any prime” criterion
            # Only return deps for those not in pivot_indices according to our pivot selection:
            pivot_set = set(pivot_indices)
            # Build basis rows from pivot_indices in this field:
            B = (V[pivot_indices, :][:, cols] % q).astype(np.int64)
            elim = _IncrementalElim(mod=q, ncols=B.shape[1])
            for k, ridx in enumerate(pivot_indices):
                elim.add_row(B[k], ridx, track_combo=True)

            for i in range(m):
                if i in pivot_set:
                    continue
                row = (V[i, cols] % q).astype(np.int64)
                combo = elim.solve_in_span(row)
                if combo is not None:
                    dependencies[i] = [(ridx, int(coeff)) for ridx, coeff in combo.items() if coeff % q != 0]

        return pivot_indices, dependencies

    raise TypeError(f"p must be int or list/np.ndarray of ints, got {type(p)}")


# ----------------------------
# Core single-prime routines
# ----------------------------

def _get_deps_single_prime(
    V: np.ndarray, p: int, compute_dependencies: bool
) -> tuple[list[int], dict[int, list[tuple[int, int]]]]:
    """
    Exact pivots + dependencies for a single prime field GF(p).
    """
    m, n = V.shape
    p = int(p)

    elim = _IncrementalElim(mod=p, ncols=n)

    pivots: list[int] = []
    deps: dict[int, list[tuple[int, int]]] = {} if compute_dependencies else {}

    for i in range(m):
        row = (V[i] % p).astype(np.int64, copy=False)
        if p == 2:
            row = (row & 1).astype(np.uint8, copy=False)

        if elim.would_increase_rank(row):
            pivots.append(i)
            elim.add_row(row, i, track_combo=compute_dependencies)
        else:
            if compute_dependencies:
                combo = elim.solve_in_span(row)
                if combo is None:
                    # should not happen if would_increase_rank returned False
                    continue
                # return as list of (pivot_row_index, coeff)
                deps[i] = [(ridx, int(coeff)) for ridx, coeff in combo.items() if coeff % p != 0]

    return pivots, deps


class _IncrementalElim:
    """
    Incremental row-space basis over GF(p) for prime p.

    Stores:
      - piv_col_to_row[piv]   : the reduced basis row with pivot at column piv (pivot value = 1)
      - piv_col_to_combo[piv] : dict mapping ORIGINAL pivot-row ids -> coefficients,
                               such that (basis row) = sum coeff * V[row_id]  (mod p)

    Then for a dependent row `v`, solve_in_span(v) returns a dict combo with
        v = sum combo[row_id] * V[row_id]  (mod p)
    """

    def __init__(self, mod: int, ncols: int):
        self.p = int(mod)
        self.ncols = int(ncols)

        self.piv_col_to_row: dict[int, np.ndarray] = {}
        self.piv_col_to_combo: dict[int, dict[int, int]] = {}
        self.pivot_cols: list[int] = []  # keep sorted for deterministic reduction order

    def would_increase_rank(self, row: np.ndarray) -> bool:
        r = self._reduce_vec_only(row)
        return self._first_nonzero_col(r) is not None

    def add_row(self, row: np.ndarray, row_id: int, track_combo: bool):
        """
        Add `row` to the basis if independent.

        If track_combo=True, maintain combinations in terms of pivot-row IDs.
        """
        if self.p == 2:
            r = (row & 1).astype(np.uint8, copy=True)
        else:
            r = (row % self.p).astype(np.int64, copy=True)

        combo: dict[int, int] | None
        if track_combo:
            combo = {int(row_id): 1 if self.p != 2 else 1}
        else:
            combo = None

        # Reduce using existing basis; for basis-building we update combo with SUBTRACTION
        for piv in self.pivot_cols:
            coeff = (r[piv] & 1) if self.p == 2 else int(r[piv] % self.p)
            if not coeff:
                continue
            prow = self.piv_col_to_row[piv]
            if self.p == 2:
                r ^= prow
                if combo is not None:
                    pc = self.piv_col_to_combo[piv]
                    # combo := combo - coeff*pc ; in GF(2), "-" == "+"
                    for k, v in pc.items():
                        combo[k] = combo.get(k, 0) ^ (v & 1)
            else:
                r = (r - coeff * prow) % self.p
                if combo is not None:
                    pc = self.piv_col_to_combo[piv]
                    for k, v in pc.items():
                        combo[k] = (combo.get(k, 0) - coeff * v) % self.p

        piv_new = self._first_nonzero_col(r)
        if piv_new is None:
            # dependent row; nothing to add
            return

        # Normalize new pivot row so pivot entry is 1, and scale combo accordingly
        if self.p != 2:
            inv = pow(int(r[piv_new]), -1, self.p)
            r = (r * inv) % self.p
            if combo is not None:
                for k in list(combo.keys()):
                    combo[k] = (combo[k] * inv) % self.p
        # p==2 already normalized

        # Insert pivot col into sorted pivot list (deterministic reduction order)
        insert_at = np.searchsorted(self.pivot_cols, piv_new)
        self.pivot_cols.insert(int(insert_at), int(piv_new))

        # Store new basis row
        self.piv_col_to_row[int(piv_new)] = r.copy()
        if track_combo:
            assert combo is not None
            # ensure mod p
            if self.p != 2:
                combo = {k: (v % self.p) for k, v in combo.items() if (v % self.p) != 0}
            else:
                combo = {k: (v & 1) for k, v in combo.items() if (v & 1) != 0}
            self.piv_col_to_combo[int(piv_new)] = combo

        # Eliminate this pivot from the other basis rows to keep RREF-like form,
        # and apply same ops to combos.
        for pc in list(self.piv_col_to_row.keys()):
            if pc == piv_new:
                continue
            basis_row = self.piv_col_to_row[pc]
            factor = (basis_row[piv_new] & 1) if self.p == 2 else int(basis_row[piv_new] % self.p)
            if not factor:
                continue

            if self.p == 2:
                self.piv_col_to_row[pc] = basis_row ^ r
                if track_combo:
                    c = self.piv_col_to_combo[pc]
                    newc = dict(c)
                    # c := c - factor*combo ; in GF(2) subtraction is XOR
                    for k, v in self.piv_col_to_combo[piv_new].items():
                        newc[k] = newc.get(k, 0) ^ (v & 1)
                    # cleanup
                    newc = {k: (v & 1) for k, v in newc.items() if (v & 1) != 0}
                    self.piv_col_to_combo[pc] = newc
            else:
                self.piv_col_to_row[pc] = (basis_row - factor * r) % self.p
                if track_combo:
                    c = self.piv_col_to_combo[pc]
                    newc = dict(c)
                    for k, v in self.piv_col_to_combo[piv_new].items():
                        newc[k] = (newc.get(k, 0) - factor * v) % self.p
                    newc = {k: (v % self.p) for k, v in newc.items() if (v % self.p) != 0}
                    self.piv_col_to_combo[pc] = newc

    def solve_in_span(self, row: np.ndarray) -> Optional[dict[int, int]]:
        """
        If row is in span(basis), return coefficients on pivot-row IDs:
            row = sum combo[row_id] * V[row_id]  (mod p)
        """
        if self.p == 2:
            r = (row & 1).astype(np.uint8, copy=True)
        else:
            r = (row % self.p).astype(np.int64, copy=True)

        combo: dict[int, int] = {}

        # Reduce; for solving, we update combo with ADDITION
        for piv in self.pivot_cols:
            coeff = (r[piv] & 1) if self.p == 2 else int(r[piv] % self.p)
            if not coeff:
                continue
            prow = self.piv_col_to_row[piv]
            if self.p == 2:
                r ^= prow
                pc = self.piv_col_to_combo[piv]
                for k, v in pc.items():
                    combo[k] = combo.get(k, 0) ^ (v & 1)
            else:
                r = (r - coeff * prow) % self.p
                pc = self.piv_col_to_combo[piv]
                for k, v in pc.items():
                    combo[k] = (combo.get(k, 0) + coeff * v) % self.p

        if self._first_nonzero_col(r) is not None:
            return None

        # cleanup zeros
        if self.p == 2:
            combo = {k: (v & 1) for k, v in combo.items() if (v & 1) != 0}
        else:
            combo = {k: (v % self.p) for k, v in combo.items() if (v % self.p) != 0}
        return combo

    def _reduce_vec_only(self, row: np.ndarray) -> np.ndarray:
        """Reduce a row using basis rows, without tracking combos."""
        if self.p == 2:
            r = (row & 1).astype(np.uint8, copy=True)
            for piv in self.pivot_cols:
                if r[piv] & 1:
                    r ^= self.piv_col_to_row[piv]
            return r
        else:
            r = (row % self.p).astype(np.int64, copy=True)
            for piv in self.pivot_cols:
                coeff = int(r[piv] % self.p)
                if coeff:
                    r = (r - coeff * self.piv_col_to_row[piv]) % self.p
            return r

    def _first_nonzero_col(self, r: np.ndarray) -> Optional[int]:
        nz = np.flatnonzero(r) if self.p == 2 else np.flatnonzero(r % self.p)
        return int(nz[0]) if nz.size else None


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
    rng = np.random.default_rng()
    for p in (2, 3, 5, 7):
        for n in (1, 2, 4):
            A = _random_invertible_matrix(p, n, rng)
            inv = gf_inv(A, p=p)
            prod = (A @ inv) % p
            assert np.array_equal(prod, np.eye(n, dtype=int) % p), f"Inverse failed for p={p}, n={n}"


def _test_gf_rref() -> None:
    rng = np.random.default_rng()
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
    rng = np.random.default_rng()
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
