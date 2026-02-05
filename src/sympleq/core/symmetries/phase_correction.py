import numpy as np
from typing import Optional
from sympleq.core.finite_field_solvers import solve_linear_system_over_gf, _select_row_basis_indices
from sympleq.core.symmetries.modular_helpers import solve_mod
from sympleq.core.circuits import Gate, Circuit


def _gf_solve_one_solution(A_int: np.ndarray, b_int: np.ndarray, p: int) -> Optional[np.ndarray]:
    """
    Wrapper around solve_linear_system_over_gf returning one solution or None if inconsistent.
    """
    try:
        sol = solve_linear_system_over_gf(A_int % p, b_int % p, p)
    except ValueError:
        return None
    return np.asarray(sol, dtype=int)


def solve_phase_vector_h_from_residual(
    tableau_in: np.ndarray,   # (N, 2n) ints; rows of input Paulis
    delta_2L: np.ndarray,     # (N,) ints mod 2L; desired phase corrections
    dimensions: np.ndarray | list[int],
    debug: bool = False,
    row_basis_cache: dict[str, np.ndarray] | None = None,
) -> Optional[np.ndarray]:
    """
    Solve A @ h ≡ delta (mod 2L) for the *linear* phase correction h.

    Qubits (p=2) are the tricky case because the modulus is 4, i.e. we're solving a linear
    system over Z_4 (a ring, not a field). We therefore use a dedicated Z_4 solver that
    reduces the problem to a pair of GF(2) systems via a standard lifting construction.
    """
    A = np.asarray(tableau_in, dtype=int)
    b = np.asarray(delta_2L, dtype=int)
    dims = np.asarray(dimensions, dtype=int)

    if A.ndim != 2:
        raise ValueError("tableau_in must be a 2D array")
    N, M = A.shape
    if M != 2 * dims.size:
        raise ValueError("tableau_in columns must equal 2 * number of qudits")

    L = int(np.lcm.reduce(dims))
    modulus = 2 * L

    # Uniform dimension fast paths
    if dims.size and np.all(dims == dims[0]):
        p_uni = int(dims[0])

        def _gf2_rref_pivots(A2: np.ndarray) -> tuple[np.ndarray, list[tuple[int, int]]]:
            """
            Return (RREF(A2), pivots) over GF(2) using XOR elimination.
            pivots is a list of (row, col) pairs in the produced RREF.
            """
            R = (np.asarray(A2, dtype=np.uint8) & 1).copy()
            m, n = R.shape
            row = 0
            pivots: list[tuple[int, int]] = []
            for col in range(n):
                if row >= m:
                    break
                nz = np.flatnonzero(R[row:, col])
                if nz.size == 0:
                    continue
                piv = row + int(nz[0])
                if piv != row:
                    R[[row, piv]] = R[[piv, row]]
                mask = R[:, col].astype(bool)
                mask[row] = False
                R[mask] ^= R[row]
                pivots.append((row, col))
                row += 1
            return R, pivots

        def _gf2_nullspace_basis(A2: np.ndarray) -> np.ndarray:
            """
            Return a GF(2) nullspace basis N with shape (n_vars, d),
            i.e. columns span {x : A2 x = 0 (mod 2)}.
            """
            R, pivots = _gf2_rref_pivots(A2)
            m, n = R.shape
            piv_cols = [c for _, c in pivots]
            free_cols = [j for j in range(n) if j not in set(piv_cols)]
            if not free_cols:
                return np.zeros((n, 0), dtype=np.uint8)

            N = np.zeros((n, len(free_cols)), dtype=np.uint8)
            # In RREF, pivot col c has a 1 in its pivot row r, so x[c] = sum_j R[r,j] x[j] for free cols.
            for k, f in enumerate(free_cols):
                N[f, k] = 1
                for r, c in pivots:
                    N[c, k] = R[r, f] & 1
            return N

        def _solve_z4_by_enumerating_mod2_solutions(A01: np.ndarray, b4: np.ndarray) -> Optional[np.ndarray]:
            """
            Solve A x = b over Z_4 for A entries in {0,1}.

            We solve the mod-2 system A x0 = b (mod 2), then attempt to lift to mod 4 by writing
            x = x0 + 2 y, where y is binary. This lift depends on the *specific* mod-2 solution x0
            (because A x0 is an integer count, not a GF(2) product), so we may need to search over
            the (typically small) GF(2) nullspace.
            """
            A2 = (np.asarray(A01, dtype=np.uint8) & 1)
            b4 = (np.asarray(b4, dtype=np.int64).reshape(-1) % 4)
            m, n = A2.shape
            if b4.shape[0] != m:
                raise ValueError("Incompatible shapes for A and b.")

            b0 = (b4 & 1).astype(np.uint8)
            x_p = _gf_solve_one_solution(A2, b0, 2)
            if x_p is None:
                return None
            x_p = (np.asarray(x_p, dtype=np.uint8) & 1)

            Nbasis = _gf2_nullspace_basis(A2)  # (n, d)
            d = int(Nbasis.shape[1])

            A_int = A2.astype(np.int64, copy=False)

            # Enumerate all nullspace assignments up to a cap (d is usually tiny in our callers).
            max_enum = 1 << 14  # 16384
            total = 1 << d if d < 31 else (max_enum + 1)
            if total > max_enum:
                # Deterministic pseudo-random subset of the solution space.
                rng = np.random.default_rng(0)
                masks = [0]
                masks.extend(rng.integers(0, 1 << min(d, 30), size=512).tolist())
            else:
                masks = range(total)

            for mask in masks:
                if d == 0:
                    x0 = x_p
                else:
                    t = np.array([(mask >> k) & 1 for k in range(d)], dtype=np.uint8)
                    x0 = (x_p ^ ((Nbasis @ t) & 1)) & 1

                Ax0_4 = (A_int @ x0.astype(np.int64)) % 4
                # b - A x0 is even because x0 solves mod2, so division by 2 is well-defined in Z2.
                r = (((b4 - Ax0_4) % 4) // 2).astype(np.uint8) & 1
                y = _gf_solve_one_solution(A2, r, 2)
                if y is None:
                    continue
                y = (np.asarray(y, dtype=np.uint8) & 1)
                x = (x0.astype(np.int64) + 2 * y.astype(np.int64)) % 4
                if np.all((A_int @ x) % 4 == b4):
                    return x

            return None

        # -------------------------
        # QUBITS: modulus=4.
        #
        # Important: the linear correction h lives in Z_4^(2n) in general (it is NOT restricted to even
        # entries). Some valid symmetries require odd components in h.
        # -------------------------
        if p_uni == 2:
            A01 = (A & 1).astype(np.uint8, copy=False)
            b4_full = (b % 4).astype(np.int64, copy=False)

            # The lifting solver is robust but slightly heavier than a plain GF(2) solve.
            # We still keep the row-basis cache for other callers (and for potential future
            # optimizations), but correctness comes first here.
            sol4 = _solve_z4_by_enumerating_mod2_solutions(A01, b4_full)
            if sol4 is None:
                return None
            if not np.all((A01.astype(np.int64) @ sol4) % 4 == b4_full):
                # Should not happen; keep a hard guard to avoid returning incorrect corrections.
                if debug:
                    print("[phase] qubit: internal solver produced non-solution; rejecting")
                return None
            return np.asarray(sol4, dtype=int) % modulus

        # ODD PRIME
        else:
            # raise NotImplementedError("Odd-prime uniform dimensions not yet implemented for phase correction.")
            A_p_full = (A % p_uni).astype(int, copy=False)
            b_p_full = (b % p_uni).astype(int, copy=False)

            rows = None
            if row_basis_cache is not None:
                rows = row_basis_cache.get("gfp")
            if rows is None or rows.size == 0:
                rows = _select_row_basis_indices(A_p_full, p_uni, A_p_full.shape[1])
                if row_basis_cache is not None:
                    row_basis_cache["gfp"] = rows

            if rows is not None and rows.size > 0:
                sol_p = _gf_solve_one_solution(A_p_full[rows], b_p_full[rows], p_uni)
                if sol_p is not None:
                    if np.all((A_p_full @ sol_p) % p_uni == b_p_full):
                        return sol_p.astype(int) % modulus
                    if debug:
                        print("[phase] GF(p): basis-row solution failed full verification; retry full")

            sol_p = _gf_solve_one_solution(A_p_full, b_p_full, p_uni)
            if sol_p is None:
                return None
            return sol_p.astype(int) % modulus

    # -------------------------
    # General composite/mixed dims: DO NOT call GF(modulus) here.
    # If you need this branch, you should implement a Z_(2L) solver (CRT / prime-power).
    # For now: return None (or raise) rather than silently doing the wrong thing.
    # -------------------------
    if debug:
        print("[phase] mixed/composite dimensions: no safe solver implemented for Z_(2L)")
    return None


def solve_linear_system_mod_prime(A: np.ndarray, b: np.ndarray, p: int) -> Optional[np.ndarray]:
    """
    Solve A x = b over Z_p where p is prime.
    Returns one particular solution (free vars = 0) or None if inconsistent.
    """
    A = np.asarray(A, dtype=np.int64) % p
    b = np.asarray(b, dtype=np.int64).reshape(-1) % p
    m, n = A.shape
    if b.shape[0] != m:
        raise ValueError("Incompatible shapes for A and b.")

    row = 0
    piv_cols = []

    if p == 2:
        # GF(2) fast path using XOR
        A2 = (A & 1).astype(np.uint8, copy=True)
        b2 = (b & 1).astype(np.uint8, copy=True)

        for col in range(n):
            if row >= m:
                break
            # find pivot
            piv = row + np.argmax(A2[row:, col])
            if piv >= m or A2[piv, col] == 0:
                continue
            if piv != row:
                A2[[row, piv]] = A2[[piv, row]]
                b2[[row, piv]] = b2[[piv, row]]

            # eliminate
            for r in range(m):
                if r != row and A2[r, col]:
                    A2[r, :] ^= A2[row, :]
                    b2[r] ^= b2[row]

            piv_cols.append(col)
            row += 1

        # inconsistency check
        for r in range(m):
            if not A2[r].any() and b2[r]:
                return None

        x = np.zeros(n, dtype=np.uint8)
        # in RREF, pivot variable equals RHS
        for r, c in enumerate(piv_cols):
            x[c] = b2[r]
        return x.astype(np.int64)

    # odd prime p
    for col in range(n):
        if row >= m:
            break
        piv = None
        for r in range(row, m):
            if A[r, col] % p != 0:
                piv = r
                break
        if piv is None:
            continue
        if piv != row:
            A[[row, piv]] = A[[piv, row]]
            b[[row, piv]] = b[[piv, row]]

        inv = pow(int(A[row, col]), -1, p)
        A[row, :] = (A[row, :] * inv) % p
        b[row] = (b[row] * inv) % p

        for r in range(m):
            if r == row:
                continue
            factor = A[r, col] % p
            if factor:
                A[r, :] = (A[r, :] - factor * A[row, :]) % p
                b[r] = (b[r] - factor * b[row]) % p

        piv_cols.append(col)
        row += 1

    # inconsistency
    for r in range(m):
        if np.all(A[r, :] == 0) and (b[r] % p) != 0:
            return None

    x = np.zeros(n, dtype=np.int64)
    for r, c in enumerate(piv_cols):
        x[c] = b[r]
    return x


def clifford_phase_decomposition(F: np.ndarray, h_F: np.ndarray,
                                 S: np.ndarray, T: np.ndarray, d: int,
                                 l_T: np.ndarray | None = None):
    """
    Inputs:
      F,h_F : target Clifford (symplectic F, phase vector h_F) with phases mod 2d
      S,T   : symplectics satisfying F = T S T^{-1}
      d     : qudit dimension
      l_T   : optional gauge vector added to h_T (default 0)

    Outputs:
      h_S, h_T : phase vectors of S and T (mod 2d)

    Conventions:
      - Pauli exponent rows update as a' = a @ F.T.
    """
    mod = 2 * d
    n2 = F.shape[0]
    n = n2 // 2
    dims = [d] * n

    def diag_U(C: np.ndarray) -> np.ndarray:
        """Return diag(C^T U C) in the same convention used by Gate.act."""
        U = np.zeros((n2, n2), dtype=int)
        U[n:, :n] = np.eye(n, dtype=int)
        return np.diag(C.T @ U @ C) % mod

    # Gauge choice for T: default ℓ_T = 0  ⇒  h_T = diag(T^T U T) + ℓ_T
    if l_T is None:
        l_T = np.zeros(n2, dtype=int)
    h_T = (diag_U(T) + l_T) % mod

    # Helper to build the composite phase for given h_S, with h_T fixed above
    def composite_phase(h_S: np.ndarray) -> np.ndarray:
        T_gate = Gate('T', list(range(n)), T, dims, h_T)
        S_gate = Gate('S', list(range(n)), S, dims, h_S)
        # Order matters: [T^{-1}, S, T] gives composite symplectic T S T^{-1}
        circuit = Circuit(dims, [T_gate.inv(), S_gate, T_gate])
        return circuit.composite_gate().phase_vector % mod

    base = composite_phase(np.zeros(n2, dtype=int))

    # Build the linear map from h_S to the resulting phase vector
    cols = []
    eye = np.eye(n2, dtype=int)
    for i in range(n2):
        cols.append((composite_phase(eye[i]) - base) % mod)
    A = np.stack(cols, axis=1) % mod
    b = (h_F % mod - base) % mod

    h_S = solve_mod(A, b, mod)
    return h_S.astype(int), h_T.astype(int)
