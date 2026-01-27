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
    Solve tableau_in @ h ≡ delta_2L (mod 2L), where h is the additional phase vector.

    Returns h or None if inconsistent.
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

    if dims.size and np.all(dims == dims[0]):
        p_uni = int(dims[0])
        if p_uni == 2:
            if np.any(b % 2 != 0):
                if debug:
                    print("[phase] qubit fallback: residual has odd entries, cannot fix")
                return None
            rows = None
            if row_basis_cache is not None:
                rows = row_basis_cache.get("gf2")
            if rows is None or rows.size == 0:
                rows = _select_row_basis_indices(A % 2, 2, A.shape[1])
                if row_basis_cache is not None:
                    row_basis_cache["gf2"] = rows
            if rows.size == 0:
                return None
            A2 = (A[rows] % 2).astype(int, copy=False)
            b2 = ((b[rows] // 2) % 2).astype(int, copy=False)
            sol2 = _gf_solve_one_solution(A2, b2, 2)
            if sol2 is not None:
                return (2 * sol2.astype(int)) % modulus
        else:
            rows = None
            if row_basis_cache is not None:
                rows = row_basis_cache.get("gfp")
            if rows is None or rows.size == 0:
                rows = _select_row_basis_indices(A % p_uni, p_uni, A.shape[1])
                if row_basis_cache is not None:
                    row_basis_cache["gfp"] = rows
            if rows.size == 0:
                return None
            A_p = (A[rows] % p_uni).astype(int, copy=False)
            b_p = (b[rows] % p_uni).astype(int, copy=False)
            sol_p = _gf_solve_one_solution(A_p, b_p, p_uni)
            if sol_p is not None:
                if debug:
                    print(f"[phase] GF({p_uni}) fallback succeeded")
                return sol_p.astype(int) % modulus

    sol = solve_linear_system_over_gf(A % modulus, b % modulus, modulus)
    if sol is not None:
        if debug:
            print("[phase] direct solve mod", modulus, "succeeded")
        return sol % modulus

    if debug:
        print("[phase] direct mod", modulus, "solve failed")

    return None


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
