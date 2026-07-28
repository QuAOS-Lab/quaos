# TODO: understand where these functions are used and either build the algebraic module
#       we need (instead of using galois) or move them where relevant.
from typing import Optional
import numpy as np
from math import gcd
import galois
from collections import defaultdict

from sympleq._typing import IntNDArray


def solve_linear_system_over_gf(coefficient_matrix: IntNDArray, right_hand_side: IntNDArray,
                                field: type | int) -> IntNDArray:
    """
    Solve the system coefficient_matrix @ solution = right_hand_side over a finite field.

    Returns one particular solution with free variables set to zero.

    Parameters
    ----------
    coefficient_matrix : IntNDArray
        Coefficient matrix.
    right_hand_side : IntNDArray
        Right-hand side vector.
    field : type or int
        The Galois field to solve over: either a `galois.GF` field class, or an
        integer prime dimension (uses a GF(2)-specialized fast path when dimension == 2).

    Returns
    -------
    IntNDArray
        One particular solution with free variables set to zero.

    Raises
    ------
    ValueError
        If the system is inconsistent, or if `coefficient_matrix` and `right_hand_side`
        have incompatible shapes.
    """
    # Fast path for GF(2) using uint8 XOR Gaussian elimination
    if field == 2 or field == galois.GF(2):
        coefficient_matrix_gf2 = (np.asarray(coefficient_matrix, dtype=np.uint8) & 1)
        right_hand_side_gf2 = (np.asarray(right_hand_side, dtype=np.uint8).reshape(-1, 1) & 1)
        m, n = coefficient_matrix_gf2.shape
        if right_hand_side_gf2.shape[0] != m:
            raise ValueError(
                f"Incompatible shapes for coefficient_matrix ({coefficient_matrix_gf2.shape}) "
                f"and right_hand_side ({right_hand_side_gf2.shape})."
            )

        R = np.hstack((coefficient_matrix_gf2, right_hand_side_gf2))
        row = 0
        pivots: list[tuple[int, int]] = []
        for col in range(n):
            if row >= m:
                break
            # find pivot
            nz = np.flatnonzero(R[row:, col])
            if nz.size == 0:
                continue
            piv = row + nz[0]
            if piv != row:
                R[[row, piv]] = R[[piv, row]]
            # eliminate other rows
            mask = R[:, col].astype(bool)
            mask[row] = False
            R[mask] ^= R[row]
            pivots.append((row, col))
            row += 1

        # inconsistency check
        for r in range(m):
            if not R[r, :n].any() and R[r, n]:
                raise ValueError("Inconsistent linear system over GF(2).")

        solution = np.zeros(n, dtype=np.uint8)
        for r, c in pivots:
            solution[c] = R[r, n]
        return solution.astype(int)

    if isinstance(field, int):
        field = galois.GF(field)

    coefficient_matrix_gf = field(coefficient_matrix)
    right_hand_side_gf = field(right_hand_side).reshape(-1, 1)

    if coefficient_matrix_gf.ndim != 2:
        raise ValueError("Coefficient matrix must be 2-dimensional.")
    if coefficient_matrix_gf.shape[0] != right_hand_side_gf.shape[0]:
        raise ValueError(
            f"Incompatible shapes for coefficient_matrix ({coefficient_matrix_gf.shape}) "
            f"and right_hand_side ({right_hand_side_gf.shape})."
        )

    m, n = coefficient_matrix_gf.shape
    augmented = np.hstack((coefficient_matrix_gf, right_hand_side_gf))
    R = augmented.copy()

    pivot_rows: list[tuple[int, int]] = []
    row = 0
    for col in range(n):
        pivot = None
        for r in range(row, m):
            if R[r, col] != field(0):
                pivot = r
                break
        if pivot is None:
            continue

        if pivot != row:
            R[[row, pivot], :] = R[[pivot, row], :]

        pivot_val = R[row, col]
        R[row, :] /= pivot_val

        for r in range(m):
            if r != row and R[r, col] != field(0):
                R[r, :] -= R[r, col] * R[row, :]

        pivot_rows.append((row, col))
        row += 1
        if row == m:
            break

    for r in range(m):
        if all(R[r, c] == field(0) for c in range(n)) and R[r, n] != field(0):
            raise ValueError("Inconsistent linear system over GF(dimension).")

    solution = field.Zeros(n)
    for row_idx, pivot_col in reversed(pivot_rows):
        value = R[row_idx, n]
        for j in range(pivot_col + 1, n):
            if R[row_idx, j] != field(0):
                value -= R[row_idx, j] * solution[j]
        solution[pivot_col] = value

    return solution.copy()


def solve_gf2(coefficient_matrix: IntNDArray, right_hand_side: IntNDArray) -> IntNDArray | None:
    """
    Solve the linear system coefficient_matrix @ solution = right_hand_side over GF(2)
    using Gaussian elimination.

    Parameters
    ----------
    coefficient_matrix : IntNDArray
        Coefficient matrix.
    right_hand_side : IntNDArray
        Right-hand side vector.

    Returns
    -------
    IntNDArray or None
        Solution vector, or None if no solution exists.
    """
    field = galois.GF(2)
    try:
        solution = solve_linear_system_over_gf(
            np.asarray(coefficient_matrix, dtype=int) % 2, np.asarray(right_hand_side, dtype=int) % 2, field
        )
    except ValueError:
        return None

    return solution.view(np.ndarray).astype(int, copy=False)


def solve_modular_linear_additive(x: int, z: int, d: int) -> int:
    """
    Find smallest non-negative integer n such that (x + z*n) % d == 0.

    Parameters
    ----------
    x : int
        Additive offset.
    z : int
        Multiplicative step.
    d : int
        Modulus.

    Returns
    -------
    int
        The smallest non-negative integer `n` satisfying (x + z*n) % d == 0.

    Raises
    ------
    ValueError
        If no such `n` exists.
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


def solve_modular_linear_system(coefficient_matrix: IntNDArray, right_hand_side: IntNDArray) -> IntNDArray:
    """
    Solve solution @ coefficient_matrix = right_hand_side over GF(dimension) using row-reduction (RREF).

    Parameters
    ----------
    coefficient_matrix : IntNDArray
        Coefficient matrix (a `galois.FieldArray`); its type determines the field.
    right_hand_side : IntNDArray
        Right-hand side vector.

    Returns
    -------
    IntNDArray
        A solution such that solution @ coefficient_matrix == right_hand_side.

    Raises
    ------
    ValueError
        If the computed solution does not satisfy solution @ coefficient_matrix == right_hand_side.
    """
    field = type(coefficient_matrix)
    solution = solve_linear_system_over_gf(coefficient_matrix.T, right_hand_side, field)
    if not np.array_equal(solution @ coefficient_matrix, right_hand_side):
        raise ValueError("Failed to solve linear system over GF(dimension).")
    return solution


def gf_solve(coefficient_matrix: IntNDArray, right_hand_side: IntNDArray, field: type) -> IntNDArray:
    """
    Solve coefficient_matrix @ solution = right_hand_side over GF(dimension).

    Parameters
    ----------
    coefficient_matrix : IntNDArray
        Coefficient matrix.
    right_hand_side : IntNDArray
        Right-hand side vector.
    field : type
        The `galois.GF` field class to solve over.

    Returns
    -------
    IntNDArray
        One particular solution, as a column vector of shape (n, 1).

    Raises
    ------
    ValueError
        If the system is inconsistent.
    """
    solution = solve_linear_system_over_gf(coefficient_matrix, right_hand_side, field)
    return solution.reshape(-1, 1)


def get_linear_dependencies(
    vectors: IntNDArray,
    dimension: int | list[int] | IntNDArray,
    compute_dependencies: bool = True,
) -> tuple[list[int], dict[int, list[tuple[int, int]]]]:
    """
    Find the linearly dependent rows of `vectors` over one or more finite fields.

    - For a single dimension (especially dimension=2): returns exact pivot rows and exact dependencies.
    - For per-column dimensions: returns correct pivot rows (same criterion as the single-dimension
      case), but dependency coefficients only if all dimensions are identical.
    - For per-row dimensions: processes each group with the fast single-dimension path.

    Parameters
    ----------
    vectors : IntNDArray
        Matrix whose rows are the candidate vectors to test for linear dependence.
    dimension : int or list of int or IntNDArray
        The prime dimension(s) defining the finite field(s). A single int applies
        GF(dimension) to all rows. A list/array applies per-row or per-column
        dimensions, depending on its length (see raised `AssertionError` below for
        the accepted lengths).
    compute_dependencies : bool
        If True, also compute the dependency coefficients for each dependent row.

    Returns
    -------
    pivot_indices : list of int
        Row indices of a linearly independent basis.
    dependencies : dict of int to list of tuple of (int, int)
        For each dependent row index, a list of (pivot_row_index, coefficient)
        pairs expressing that row as a linear combination of the pivot rows.
        Empty if `compute_dependencies` is False.

    Raises
    ------
    TypeError
        If `dimension` is neither an int nor a list/array of ints.
    AssertionError
        If `dimension` is a list/array whose length matches neither the number of
        rows, the number of columns, nor half the number of columns.

    Notes
    -----
    For dimension=2, dependencies are coefficients in GF(2) (0/1).
    For odd prime dimension, dependencies are coefficients in GF(dimension).
    """
    vectors = np.asarray(vectors, dtype=np.int64)
    m, n = vectors.shape

    # Normalize dimension into mode and per-column dimensions when possible.
    if isinstance(dimension, int):
        dimension_int = int(dimension)
        return _get_deps_single_prime(vectors, dimension_int, compute_dependencies)

    if isinstance(dimension, (list, np.ndarray)):
        dimension = np.asarray(dimension, dtype=int)
        lp = len(dimension)

        if lp == m:
            # per-row legacy: group and process each group
            pivot_indices: list[int] = []
            dependencies: dict[int, list[tuple[int, int]]] = {}
            prime_groups = defaultdict(list)
            for i, prime in enumerate(dimension.tolist()):
                prime_groups[int(prime)].append(i)
            for prime, idxs in prime_groups.items():
                piv, deps = _get_deps_single_prime(vectors[idxs, :], int(prime), compute_dependencies)
                pivot_indices.extend([idxs[j] for j in piv])
                if compute_dependencies:
                    for local_row, expr in deps.items():
                        dependencies[idxs[local_row]] = [(idxs[piv_j], c) for piv_j, c in expr]
            return pivot_indices, dependencies

        if lp == n:
            dimension_cols = dimension.tolist()
        elif lp == n // 2 and n % 2 == 0:
            dimension_cols = np.repeat(dimension, 2).tolist()
        else:
            raise AssertionError(f"Length of dimension must be rows={m}, cols={n}, or cols/2={n // 2}. Got {lp}")

        # Per-column dimensions: pivot criterion is “independent if increases rank for any prime subset”.
        # We can do this by running elimination separately per distinct prime on its column subset,
        # but WITHOUT galois row_reduce. Just incremental elimination mod dimension.

        prime_cols = defaultdict(list)
        for j, prime in enumerate(dimension_cols):
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
                row = vectors[i, cols] % q
                if eliminators[q].would_increase_rank(row):
                    independent = True
                    break

            if independent:
                pivot_indices.append(i)
                # actually add to all eliminators (mirrors your “updated_seen” logic)
                for q in primes:
                    cols = prime_cols[q]
                    row = vectors[i, cols] % q
                    eliminators[q].add_row(row, i, track_combo=False)

        # dependencies only if exactly one prime
        if compute_dependencies and len(primes) == 1:
            q = primes[0]
            cols = prime_cols[q]
            piv, deps = _get_deps_single_prime(vectors[:, cols], q, compute_dependencies=True)
            # piv returned are *row indices* already, but might differ slightly from the “any prime” criterion
            # Only return deps for those not in pivot_indices according to our pivot selection:
            pivot_set = set(pivot_indices)
            # Build basis rows from pivot_indices in this field:
            basis_rows = (vectors[pivot_indices, :][:, cols] % q).astype(np.int64)
            elim = _IncrementalElim(mod=q, ncols=basis_rows.shape[1])
            for k, ridx in enumerate(pivot_indices):
                elim.add_row(basis_rows[k], ridx, track_combo=True)

            for i in range(m):
                if i in pivot_set:
                    continue
                row = (vectors[i, cols] % q).astype(np.int64)
                combo = elim.solve_in_span(row)
                if combo is not None:
                    dependencies[i] = [(ridx, int(coeff)) for ridx, coeff in combo.items() if coeff % q != 0]

        return pivot_indices, dependencies

    raise TypeError(f"dimension must be int or list/np.ndarray of ints, got {type(dimension)}")


# ----------------------------
# Core single-prime routines
# ----------------------------

def _get_deps_single_prime(
    vectors: IntNDArray, dimension: int, compute_dependencies: bool
) -> tuple[list[int], dict[int, list[tuple[int, int]]]]:
    """
    Exact pivots + dependencies for a single prime field GF(dimension).

    Parameters
    ----------
    vectors : IntNDArray
        Matrix whose rows are the candidate vectors to test for linear dependence.
    dimension : int
        Prime defining the finite field GF(dimension).
    compute_dependencies : bool
        If True, also compute the dependency coefficients for each dependent row.

    Returns
    -------
    pivots : list of int
        Row indices of a linearly independent basis.
    deps : dict of int to list of tuple of (int, int)
        For each dependent row index, a list of (pivot_row_index, coefficient)
        pairs expressing that row as a linear combination of the pivot rows.
        Empty if `compute_dependencies` is False.
    """
    m, n = vectors.shape
    dimension = int(dimension)

    elim = _IncrementalElim(mod=dimension, ncols=n)

    pivots: list[int] = []
    deps: dict[int, list[tuple[int, int]]] = {} if compute_dependencies else {}

    for i in range(m):
        row = (vectors[i] % dimension).astype(np.int64, copy=False)
        if dimension == 2:
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
                deps[i] = [(ridx, int(coeff)) for ridx, coeff in combo.items() if coeff % dimension != 0]

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

        self.piv_col_to_row: dict[int, IntNDArray] = {}
        self.piv_col_to_combo: dict[int, dict[int, int]] = {}
        self.pivot_cols: list[int] = []  # keep sorted for deterministic reduction order

    def would_increase_rank(self, row: IntNDArray) -> bool:
        r = self._reduce_vec_only(row)
        return self._first_nonzero_col(r) is not None

    def add_row(self, row: IntNDArray, row_id: int, track_combo: bool):
        """
        Add `row` to the basis if independent.

        Parameters
        ----------
        row : IntNDArray
            Candidate row to add to the basis.
        row_id : int
            Original row index of `row`, used to key its coefficients when
            `track_combo` is True.
        track_combo : bool
            If True, maintain combinations in terms of pivot-row IDs.
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

    def solve_in_span(self, row: IntNDArray) -> Optional[dict[int, int]]:
        """
        If row is in span(basis), return coefficients on pivot-row IDs:
            row = sum combo[row_id] * V[row_id]  (mod p)

        Parameters
        ----------
        row : IntNDArray
            Row to express in terms of the basis.

        Returns
        -------
        dict of int to int or None
            Mapping from pivot-row ID to coefficient, or None if `row` is not
            in the span of the current basis.
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

    def _reduce_vec_only(self, row: IntNDArray) -> IntNDArray:
        """
        Reduce a row using basis rows, without tracking combos.

        Parameters
        ----------
        row : IntNDArray
            Row to reduce.

        Returns
        -------
        IntNDArray
            The row reduced against the current basis.
        """
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

    def _first_nonzero_col(self, r: IntNDArray) -> Optional[int]:
        nz = np.flatnonzero(r) if self.p == 2 else np.flatnonzero(r % self.p)
        return int(nz[0]) if nz.size else None


def gf_inv(matrix: IntNDArray, dimension: int = 2) -> IntNDArray:
    """
    Compute the inverse of a square matrix over GF(dimension) for a prime dimension.

    Defaults to GF(2) for backwards compatibility.

    Parameters
    ----------
    matrix : IntNDArray
        Square matrix to invert.
    dimension : int
        Prime defining the finite field GF(dimension).

    Returns
    -------
    IntNDArray
        The inverse of `matrix` over GF(dimension).

    Raises
    ------
    ValueError
        If `matrix` is not square, or is singular over GF(dimension).
    """
    matrix = np.asarray(matrix, dtype=int) % dimension
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Matrix must be square to compute an inverse over GF(dimension).")

    n = matrix.shape[0]
    Id = np.eye(n, dtype=int) % dimension
    AI = np.concatenate((matrix.copy(), Id), axis=1)

    for i in range(n):
        pivot_row = None
        for r in range(i, n):
            if AI[r, i] % dimension != 0:
                pivot_row = r
                break
        if pivot_row is None:
            raise ValueError("Matrix is singular over GF(dimension); inverse does not exist.")
        if pivot_row != i:
            AI[[i, pivot_row]] = AI[[pivot_row, i]]

        pivot_val = int(AI[i, i] % dimension)
        pivot_inv = pow(pivot_val, -1, dimension)
        AI[i, :] = (AI[i, :] * pivot_inv) % dimension

        for r in range(n):
            if r == i:
                continue
            factor = AI[r, i] % dimension
            if factor != 0:
                AI[r, :] = (AI[r, :] - factor * AI[i, :]) % dimension

    return AI[:, n:] % dimension


def gf2_inv(matrix: IntNDArray) -> IntNDArray:
    """Invert a square GF(2) matrix using XOR elimination.

    Parameters
    ----------
    matrix : IntNDArray
        Square matrix with entries in {0,1}. Any integer dtype is accepted.

    Returns
    -------
    IntNDArray
        Inverse matrix over GF(2), dtype=uint8.

    Raises
    ------
    np.linalg.LinAlgError
        If the matrix is not square or is singular over GF(2).
    """
    matrix = np.asarray(matrix, dtype=np.uint8).copy() & 1
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise np.linalg.LinAlgError("GF2 inverse requires square matrix")

    n = matrix.shape[0]
    aug = np.hstack((matrix, np.eye(n, dtype=np.uint8)))

    row = 0
    for col in range(n):
        if row >= n:
            break
        nz = np.flatnonzero(aug[row:, col])
        if nz.size == 0:
            continue
        piv = row + int(nz[0])
        if piv != row:
            aug[[row, piv]] = aug[[piv, row]]
        mask = aug[:, col].astype(bool)
        mask[row] = False
        aug[mask] ^= aug[row]
        row += 1

    if not np.array_equal(aug[:, :n] & 1, np.eye(n, dtype=np.uint8)):
        raise np.linalg.LinAlgError("matrix is singular over GF(2)")

    return aug[:, n:] & 1


def gf_rref(matrix: IntNDArray, dimension: int = 2) -> tuple[IntNDArray, IntNDArray, IntNDArray, int]:
    """
    Compute the reduced row echelon form of a matrix over GF(dimension) for a prime dimension.

    Defaults to GF(2) for backwards compatibility.

    Parameters
    ----------
    matrix : IntNDArray
        Matrix to row-reduce.
    dimension : int
        Prime defining the finite field GF(dimension).

    Returns
    -------
    IntNDArray
        The reduced row echelon form of `matrix`.
    IntNDArray
        The row_transform_matrix (row operations), such that row_transform_matrix @ matrix ==
        the RREF (before column ops).
    IntNDArray
        The column_transform_matrix (column operations).
    int
        The rank of `matrix`.

    Raises
    ------
    ValueError
        If `matrix` is not 2-dimensional.
    """
    matrix = np.asarray(matrix, dtype=int) % dimension
    if matrix.ndim != 2:
        raise ValueError("Input matrix must be 2-dimensional for RREF.")

    matrix = matrix.copy()
    m, n = matrix.shape
    i = j = 0
    row_transform_matrix = np.eye(m, dtype=int) % dimension
    column_transform_matrix = np.eye(n, dtype=int) % dimension

    while i < m and j < n:
        pivot_row = None
        for r in range(i, m):
            if matrix[r, j] % dimension != 0:
                pivot_row = r
                break

        if pivot_row is None:
            j += 1
            continue

        if pivot_row != i:
            matrix[[i, pivot_row]] = matrix[[pivot_row, i]]
            row_transform_matrix[[i, pivot_row]] = row_transform_matrix[[pivot_row, i]]

        pivot_val = int(matrix[i, j] % dimension)
        pivot_inv = pow(pivot_val, -1, dimension)
        matrix[i, :] = (matrix[i, :] * pivot_inv) % dimension
        row_transform_matrix[i, :] = (row_transform_matrix[i, :] * pivot_inv) % dimension

        for r in range(m):
            if r == i:
                continue
            factor = matrix[r, j] % dimension
            if factor != 0:
                matrix[r, :] = (matrix[r, :] - factor * matrix[i, :]) % dimension
                row_transform_matrix[r, :] = (
                    row_transform_matrix[r, :] - factor * row_transform_matrix[i, :]
                ) % dimension

        for c in range(n):
            if c == j:
                continue
            factor = matrix[i, c] % dimension
            if factor != 0:
                matrix[:, c] = (matrix[:, c] - factor * matrix[:, j]) % dimension
                column_transform_matrix[:, c] = (
                    column_transform_matrix[:, c] - factor * column_transform_matrix[:, j]
                ) % dimension

        i += 1
        j += 1

    rank = i
    return matrix % dimension, row_transform_matrix % dimension, column_transform_matrix % dimension, rank


def gf_lu(matrix: IntNDArray, dimension: int = 2) -> tuple[IntNDArray, IntNDArray, IntNDArray]:
    """
    Perform LU decomposition of a matrix over GF(dimension) for prime dimension.

    Defaults to GF(2) for backwards compatibility.

    Parameters
    ----------
    matrix : IntNDArray
        Square matrix to decompose.
    dimension : int
        Prime defining the finite field GF(dimension).

    Returns
    -------
    IntNDArray
        lower_triangular_matrix, the lower-triangular factor with unit diagonal.
    IntNDArray
        upper_triangular_matrix, the upper-triangular factor.
    IntNDArray
        permutation_matrix, such that permutation_matrix @ matrix ==
        lower_triangular_matrix @ upper_triangular_matrix.

    Raises
    ------
    ValueError
        If `matrix` is not square.
    """
    matrix = np.asarray(matrix, dtype=int) % dimension
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("LU decomposition requires a square matrix over GF(dimension).")

    m = matrix.shape[0]
    upper_triangular_matrix = matrix.copy()
    lower_triangular_matrix = np.eye(m, dtype=int) % dimension
    permutation_matrix = np.eye(m, dtype=int)

    for k in range(m):
        pivot_candidates = np.where(upper_triangular_matrix[k:, k] % dimension != 0)[0]
        if pivot_candidates.size == 0:
            continue
        pivot = pivot_candidates[0] + k
        if pivot != k:
            upper_triangular_matrix[[k, pivot], k:] = upper_triangular_matrix[[pivot, k], k:]
            lower_triangular_matrix[[k, pivot], :k] = lower_triangular_matrix[[pivot, k], :k]
            permutation_matrix[[k, pivot]] = permutation_matrix[[pivot, k]]

        pivot_val = int(upper_triangular_matrix[k, k] % dimension)
        pivot_inv = pow(pivot_val, -1, dimension)

        for j in range(k + 1, m):
            factor = (upper_triangular_matrix[j, k] * pivot_inv) % dimension
            lower_triangular_matrix[j, k] = factor
            if factor != 0:
                upper_triangular_matrix[j, k:] = (
                    upper_triangular_matrix[j, k:] - factor * upper_triangular_matrix[k, k:]
                ) % dimension

    return lower_triangular_matrix % dimension, upper_triangular_matrix % dimension, permutation_matrix


def _select_row_basis_indices(matrix: IntNDArray, dimension: int, max_rows: int) -> IntNDArray:
    """
    Return indices of a greedily chosen row basis of `matrix` over GF(dimension).

    The routine performs a light Gauss-Jordan sweep, moving left-to-right across
    columns and picking the first unused row with a non-zero entry as the pivot.
    Each pivot row is normalized and used to eliminate the pivot column from the
    other unused rows. Collection stops once either all columns are processed or
    `max_rows` independent rows have been gathered.

    Parameters
    ----------
    matrix : IntNDArray
        Integer matrix whose rows are candidate equations.
    dimension : int
        Prime for the finite field GF(dimension).
    max_rows : int
        Maximum number of independent rows to return (often the number of columns).

    Returns
    -------
    IntNDArray
        1-D array of row indices that are linearly independent over GF(dimension).
    """
    field = galois.GF(dimension)
    matrix = field(matrix % dimension).copy()
    num_rows, num_cols = matrix.shape
    used = np.zeros(num_rows, dtype=bool)
    basis = []
    col = 0
    for _ in range(num_rows):
        if col >= num_cols:
            break
        pick = None
        for r in range(num_rows):
            if used[r]:
                continue
            if matrix[r, col] != field(0):
                pick = r
                break
        if pick is None:
            col += 1
            continue
        basis.append(pick)
        used[pick] = True
        inv = field(1) / matrix[pick, col]
        matrix[pick, :] = matrix[pick, :] * inv
        for r in range(num_rows):
            if r == pick or used[r]:
                continue
            if matrix[r, col] != field(0):
                matrix[r, :] = matrix[r, :] - matrix[r, col] * matrix[pick, :]
        col += 1
        if len(basis) >= max_rows:
            break
    return np.array(basis, dtype=int)
