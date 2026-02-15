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
    tableau_in: np.ndarray,
    delta_2L: np.ndarray,
    dimensions: np.ndarray | list[int],
    debug: bool = False,
    row_basis_cache: dict[str, np.ndarray] | None = None,
) -> Optional[np.ndarray]:
    """
    Solve A @ h ≡ delta (mod 2L) for the *linear* phase correction h.

    For qubits (p=2): h must be even mod 4, so solve over GF(2) for t where h = 2 t (mod 4):
        (A mod 2) t = (delta/2 mod 2)
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

        # -------------------------
        # QUBITS: modulus=4, solve in GF(2) on delta/2
        # -------------------------
        if p_uni == 2:
            # delta must be even for a Pauli-frame correction to exist
            if np.any(b & 1):
                if debug:
                    print("[phase] qubit: residual has odd entries -> no linear correction")
                return None

            A2_full = (A & 1).astype(int, copy=False)
            b2_full = ((b // 2) & 1).astype(int, copy=False)

            # Try cached/basis rows first
            rows = None
            if row_basis_cache is not None:
                rows = row_basis_cache.get("gf2")
            if rows is None or rows.size == 0:
                rows = _select_row_basis_indices(A2_full, 2, A2_full.shape[1])
                if row_basis_cache is not None:
                    row_basis_cache["gf2"] = rows

            # If we got a small row set, solve it; else fall back to full
            if rows is not None and rows.size > 0:
                sol2 = solve_linear_system_mod_prime(A2_full[rows], b2_full[rows], 2)
                if sol2 is not None:
                    # verify
                    if np.all((A2_full @ sol2) % 2 == b2_full):
                        return (2 * sol2.astype(int)) % modulus
                    # Otherwise redo with full system
                    if debug:
                        print("[phase] qubit: basis-row solution failed full verification; retry full")

            sol2 = solve_linear_system_mod_prime(A2_full, b2_full, 2)
            if sol2 is None:
                return None
            return (2 * sol2.astype(int)) % modulus

        # ODD PRIME
        else:
            raise NotImplementedError("Odd-prime uniform dimensions not yet implemented for phase correction.")
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
    all_qudit_indices = tuple(range(n))

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
        T_gate = Gate('T', T, h_T)
        S_gate = Gate('S', S, h_S)
        # Order matters: [T^{-1}, S, T] gives composite symplectic T S T^{-1}
        circuit = Circuit.from_gates_and_qudits(dims, [T_gate.inverse(), S_gate, T_gate],
                                                      [all_qudit_indices, all_qudit_indices, all_qudit_indices])
        return circuit.composite_gate().phase_vector() % mod

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
