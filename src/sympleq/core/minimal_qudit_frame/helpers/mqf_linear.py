from __future__ import annotations
import numpy as np
from sympleq.core.symmetries import modular_helpers as _modular_helpers
from sympleq.core.symmetries.modular_helpers import (
    independent_columns,
    inv_mod_mat,
    inv_mod_scalar,
    mod_p,
    nullspace_mod,
    omega_matrix,
    rank_mod,
)
# ``mat_pow_mod`` is re-exported from modular_helpers so existing
# ``from .mqf_linear import mat_pow_mod`` callers keep working (B7: single
# shared implementation rather than a per-module copy).
mat_pow_mod = _modular_helpers.mat_pow_mod


def symplectic_left_inverse(T: np.ndarray, p: int) -> np.ndarray:
    """
    For T with T^T Ω T invertible (typically Ω_k), return L such that L T = I on span(T).
    L = (T^T Ω T)^{-1} T^T Ω
    """
    n2 = T.shape[0]
    Omega = omega_matrix(n2 // 2, p)
    J = mod_p(T.T @ Omega @ T, p)
    if rank_mod(J, p) != J.shape[0]:
        raise RuntimeError("symplectic_left_inverse: T^T Ω T is singular; span(T) is degenerate.")
    Jinv = inv_mod_mat(J, p)
    return mod_p(Jinv @ T.T @ Omega, p)


def restrict_operator(F: np.ndarray, T: np.ndarray, p: int) -> np.ndarray:
    """
    Return F_restricted in the T-coordinates:  F_T = T^{-1} F T
    where T^{-1} means symplectic left inverse above.
    """
    L = symplectic_left_inverse(T, p)
    return mod_p(L @ F @ T, p)


def symplectic_orthogonal_complement_in_ambient(T: np.ndarray, p: int) -> np.ndarray:
    """
    Return a basis (columns) for W^⊥ in the ambient space, where W = span(T).
    W^⊥ = {x : T^T Omega x = 0}.
    """
    n2 = T.shape[0]
    Omega = omega_matrix(n2 // 2, p)
    A = mod_p(T.T @ Omega, p)
    N = nullspace_mod(A, p)
    return independent_columns(N, p)


def kernel_in_span(A: np.ndarray, span_basis: np.ndarray, p: int) -> np.ndarray:
    """
    Return a basis (ambient columns) for { x in span(span_basis) : A x = 0 }.

    If span_basis is n×d (columns spanning subspace U), then any x ∈ U is x = span_basis c.
    Constraint A x = 0 becomes (A span_basis) c = 0, so c ∈ ker(A span_basis).
    """
    span_basis = independent_columns(mod_p(span_basis, p), p)
    if span_basis.size == 0:
        return span_basis
    AB = mod_p(A @ span_basis, p)
    C = nullspace_mod(AB, p)            # d×k coefficients
    X = mod_p(span_basis @ C, p)        # n×k ambient vectors
    return independent_columns(X, p)


def span_columns(cols: np.ndarray, p: int) -> np.ndarray:
    """Return an independent column basis of span(cols)."""
    return independent_columns(mod_p(cols, p), p)


def symp_pairing_matrix(Omega: np.ndarray, A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """Return A^T Ω B mod p."""
    return mod_p(A.T @ Omega @ B, p)


def is_nondegenerate(Omega: np.ndarray, B: np.ndarray, p: int) -> bool:
    """True iff the restricted form on span(B) is nondegenerate."""
    B = independent_columns(mod_p(B, p), p)
    d = B.shape[1]
    if d == 0:
        return True
    G = mod_p(B.T @ Omega @ B, p)
    return rank_mod(G, p) == d


def symplectic_orthogonal_complement(Omega: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """
    Return a basis C for B^⊥ in the ambient space:
        B^T Ω C = 0.
    """
    B = independent_columns(mod_p(B, p), p)
    if B.size == 0:
        return np.eye(Omega.shape[0], dtype=np.int64)
    A = mod_p(B.T @ Omega, p)  # (d×n)
    C = nullspace_mod(A, p)    # (n×k)
    return independent_columns(C, p)


def symplectic_orthogonal_complement_in_span(
    Omega: np.ndarray, B: np.ndarray, span_basis: np.ndarray, p: int
) -> np.ndarray:
    """
    Return a basis (ambient) for (span(span_basis) ∩ B^⊥).
    """
    span_basis = independent_columns(mod_p(span_basis, p), p)
    if span_basis.size == 0:
        return span_basis
    A = mod_p(B.T @ Omega, p)  # constraints A x = 0
    return kernel_in_span(A, span_basis, p)


def darboux_basis_from_span(Omega: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """
    Given B (nx2h) whose columns span a nondegenerate symplectic subspace,
    produce T (nx2h) such that T^T Ω T = Ω_h.

    Deterministic symplectic Gram-Schmidt on columns.
    """
    B = independent_columns(mod_p(B, p), p)
    d = B.shape[1]
    if d % 2 != 0:
        raise ValueError("darboux_basis_from_span: span dimension must be even.")
    if d == 0:
        return B

    # Working set of remaining vectors
    S = B.copy()

    U_cols: list[np.ndarray] = []
    V_cols: list[np.ndarray] = []

    def pair(u: np.ndarray, v: np.ndarray) -> int:
        return int(mod_p(u.T @ Omega @ v, p).reshape(()))

    while S.shape[1] > 0:
        # pick a nonzero u (first column)
        u = S[:, 0:1]
        S = S[:, 1:]

        # Find v in span(S) with <u,v> != 0.
        # If none of the basis columns work, solve for a linear combination.
        v: np.ndarray | None = None
        for j in range(S.shape[1]):
            cand = S[:, j:j + 1]
            if pair(u, cand) % p != 0:
                v = cand
                # remove that column from S
                S = np.concatenate([S[:, :j], S[:, j + 1:]], axis=1) if j + 1 <= S.shape[1] else S[:, :j]
                break

        if v is None:
            # Solve (u^T Ω S) c = 1 for c
            if S.shape[1] == 0:
                raise RuntimeError("darboux_basis_from_span: nondegeneracy violated (no partner).")
            r = mod_p(u.T @ Omega @ S, p)     # 1×m
            # Try pick a column with r_j != 0 (fast)
            nz = np.where(r.reshape(-1) % p != 0)[0]
            if nz.size:
                j = int(nz[0])
                c = np.zeros((S.shape[1], 1), dtype=np.int64)
                c[j, 0] = inv_mod_scalar(int(r[0, j]), p)
            else:
                # If r is identically zero, then u pairs with nothing in span(S): contradiction
                raise RuntimeError("darboux_basis_from_span: nondegeneracy violated (u pairs with nothing).")
            v = mod_p(S @ c, p)

        # Orthogonalize u,v against existing pairs
        for u_prev, v_prev in zip(U_cols, V_cols):
            # u <- u - <u,v_prev> u_prev + <u,u_prev> v_prev
            a = pair(u, v_prev) % p
            g = pair(u, u_prev) % p
            if a or g:
                u = mod_p(u - a * u_prev + g * v_prev, p)

            # v <- v - <v,v_prev> u_prev + <v,u_prev> v_prev
            cu = pair(v, v_prev) % p
            cv = pair(v, u_prev) % p
            if cu or cv:
                v = mod_p(v - cu * u_prev + cv * v_prev, p)

        if np.all(u % p == 0) or np.all(v % p == 0):
            raise RuntimeError("darboux_basis_from_span: orthogonalization collapsed vectors.")

        # Normalize so that <u,v> = 1
        beta = pair(u, v) % p
        if beta == 0:
            raise RuntimeError("darboux_basis_from_span: failed to create a hyperbolic pair.")
        v = mod_p(v * inv_mod_scalar(beta, p), p)

        U_cols.append(u)
        V_cols.append(v)

        # Make remaining S orthogonal to the new pair:
        # w <- w - <w,v> u + <w,u> v
        if S.shape[1] > 0:
            alpha = mod_p(S.T @ Omega @ v, p)  # m×1 with alpha_j = <w_j,v>
            beta2 = mod_p(S.T @ Omega @ u, p)  # m×1 with beta_j  = <w_j,u>
            S = mod_p(S - u @ alpha.T + v @ beta2.T, p)
            S = independent_columns(S, p)

    T = np.concatenate(U_cols + V_cols, axis=1)
    return mod_p(T, p)


def split_uv(T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Given canonical [U | V], return (U, V)."""
    k2 = T.shape[1]
    if k2 % 2 != 0:
        raise ValueError("split_uv: expected even number of columns")
    k = k2 // 2
    return T[:, :k], T[:, k:]


def symplectic_completion_from_block(T_blk: np.ndarray, p: int) -> np.ndarray:
    """
    Given T_blk (n2x2k) whose columns form a canonical symplectic basis of a nondegenerate subspace W,
    complete to a full symplectic basis T_full (n2xn2) with column order:
        [U_W | U_perp | V_W | V_perp].
    """
    n2 = T_blk.shape[0]
    Omega = omega_matrix(n2 // 2, p)

    # Complement basis (ambient columns)
    N = symplectic_orthogonal_complement(Omega, T_blk, p)  # n2×(n2-2k)
    if N.shape[1] % 2 != 0:
        raise RuntimeError("symplectic_completion_from_block: complement has odd dimension")

    T_perp = darboux_basis_from_span(Omega, N, p) if N.shape[1] else np.zeros((n2, 0), dtype=np.int64)
    U_W, V_W = split_uv(T_blk)
    U_p, V_p = split_uv(T_perp) if T_perp.shape[1] else (T_perp, T_perp)

    T_full = np.concatenate([U_W, U_p, V_W, V_p], axis=1)
    G = mod_p(T_full.T @ Omega @ T_full, p)
    if not np.array_equal(G % p, Omega % p):
        raise RuntimeError("symplectic_completion_from_block: completion is not symplectic")
    return T_full
