import numpy as np
from ..modular_helpers import mod_p, solve_linear, independent_columns, rank_mod, solve_linear_many, mat_pow_mod
from .module_invariants import cyclic_submodule_basis
from .atomic_extension import top_quotient_E_basis


def _col(v: np.ndarray) -> np.ndarray:
    """Ensure v is a column vector (n×1)."""
    v = np.asarray(v, dtype=np.int64)
    return v.reshape((-1, 1))


def krylov_chain(op: np.ndarray, v: np.ndarray, L: int, p: int) -> np.ndarray:
    """
    Return [v, op v, op^2 v, ..., op^(L-1) v] as an n×L matrix (columns).
    Column-action convention.
    """
    op = mod_p(op, p)
    v = _col(v)
    n = op.shape[0]
    L = int(L)
    out = np.zeros((n, L), dtype=np.int64)
    w = v.copy()
    for i in range(L):
        out[:, i:i + 1] = mod_p(w, p)
        w = mod_p(op @ w, p)
    return out


def krylov_closure(op: np.ndarray, v: np.ndarray, p: int, cap: int | None = None) -> np.ndarray:
    """
    Stabilized Krylov span K = span{v, op v, op^2 v, ...} until rank stops increasing.
    Returns an n×k basis matrix (columns), independent_columns-selected.
    """
    op = mod_p(op, p)
    v = _col(v)
    n = op.shape[0]
    cap_local = n if cap is None else int(cap)

    K = np.zeros((n, 0), dtype=np.int64)
    w = v.copy()
    r_prev = 0

    for _ in range(cap_local):
        trial = np.concatenate([K, w], axis=1)
        r = rank_mod(trial, p)
        if r > r_prev:
            K = trial
            r_prev = r
            w = mod_p(op @ w, p)
        else:
            break

    return independent_columns(K, p)


def krylov_closure_in_span(
    op: np.ndarray, v: np.ndarray, span_basis: np.ndarray, p: int, cap: int | None = None
) -> np.ndarray:
    """
    Stabilized Krylov span inside span(span_basis).

    Contract: with v in span(span_basis) and span_basis op-invariant, every
    Krylov vector op^t v lies in span(span_basis). This is now *enforced* rather
    than assumed: each Krylov vector is expressed in the span_basis coordinates
    (solve_linear_many) and the reconstruction is checked. A containment failure
    means op does not leave span(span_basis) invariant (or v is outside it), and
    raises a clear error instead of silently returning ambient vectors.

    Returns ambient columns (n x k).
    """
    span_basis = independent_columns(mod_p(span_basis, p), p)
    if span_basis.shape[1] == 0:
        return span_basis
    K = krylov_closure(op, v, p, cap=cap)
    if K.shape[1] == 0:
        return K

    # Express K in span_basis coordinates: span_basis @ C ?= K.
    try:
        C = solve_linear_many(span_basis, mod_p(K, p), p)
    except RuntimeError as exc:
        raise RuntimeError(
            "krylov_closure_in_span: Krylov vectors are not contained in span(span_basis); "
            "op does not leave the span invariant, or v lies outside it."
        ) from exc

    K_in = mod_p(span_basis @ C, p)
    if not np.array_equal(K_in % p, mod_p(K, p) % p):
        raise RuntimeError(
            "krylov_closure_in_span: span containment check failed; "
            "span(span_basis) is not op-invariant for this seed."
        )
    return independent_columns(K_in, p)


def build_partner_in_span(
    K: np.ndarray,
    Omega: np.ndarray,
    span_basis: np.ndarray,
    p: int,
    *,
    return_none: bool = False,
) -> np.ndarray | None:
    """
    Find z in span(span_basis) such that:
        <K[:,b], z> = δ_{b,k-1}  for b=0..k-1,
    where K is nxk (columns) and <u,v> := u^T Omega v.

    Returns:
        z (nx1) column vector in the ambient space.

    If return_none=True, returns None when no solution exists in span(span_basis).
    """
    K = mod_p(K, p)
    Omega = mod_p(Omega, p)
    span_basis = independent_columns(mod_p(span_basis, p), p)

    n = K.shape[0]
    if Omega.shape != (n, n):
        raise ValueError("build_partner_in_span: Omega has incompatible shape.")
    if span_basis.shape[0] != n:
        raise ValueError("build_partner_in_span: span_basis has incompatible shape.")

    k = K.shape[1]
    if k == 0:
        raise ValueError("build_partner_in_span: K must have at least one column.")
    if span_basis.shape[1] == 0:
        if return_none:
            return None
        raise RuntimeError("build_partner_in_span: span_basis is empty.")

    # Solve (K^T Ω span_basis) c = e_{k-1}
    A = mod_p(K.T @ Omega @ span_basis, p)  # k×d
    b = np.zeros((k, 1), dtype=np.int64)
    b[-1, 0] = 1

    try:
        coeff = solve_linear(A, b, p)      # d×1
    except RuntimeError:
        if return_none:
            return None
        raise

    z = mod_p(span_basis @ coeff, p)        # n×1
    return z


def select_module_generators_from_top_quotient(
    Fp: np.ndarray,
    Np: np.ndarray,
    top_candidates: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
    *,
    denom: np.ndarray | None = None,
) -> np.ndarray:
    """
    Select one GF(p)[F]/(q)-module generator per q^L block from a top quotient.

    This is the quotient-orbit version of the chain-head selection.  The input
    ``top_candidates`` is a GF(p)-basis of representatives for

        K_L / (K_{L-1} + N K_{L+1}).

    Its base-field dimension is deg_q * multiplicity.  We select only one
    representative from each deg_q-dimensional F-orbit in this quotient.  This
    avoids the common mistake of treating all deg_q*multiplicity basis vectors
    as separate cyclic blocks.
    """
    try:
        tq = top_quotient_E_basis(
            Fp, Np, int(L), int(deg_q), int(p),
            top_reps=top_candidates,
            denom=denom,
        )
    except RuntimeError:
        # Preserve the old public behaviour: a sector builder may treat an empty
        # return as a failed candidate and continue to another length.  Strict
        # callers should use the sector diagnostics to distinguish failure modes.
        return np.zeros((Fp.shape[0], 0), dtype=np.int64)

    if not tq.e_basis_reps:
        return np.zeros((Fp.shape[0], 0), dtype=np.int64)
    return np.concatenate(tq.e_basis_reps, axis=1)
