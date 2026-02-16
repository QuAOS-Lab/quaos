import numpy as np
from .modular_helpers import mod_p, _solve_linear, independent_columns, rank_mod
from .module_invariants import cyclic_submodule_basis


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
    We iteratively add op^t v, but finally project the resulting basis back into the span
    using independent_columns (the caller should ensure v ∈ span_basis and op-invariance
    of span_basis when required).
    Returns ambient columns (nxk).
    """
    span_basis = independent_columns(mod_p(span_basis, p), p)
    if span_basis.shape[1] == 0:
        return span_basis
    K = krylov_closure(op, v, p, cap=cap)
    # ensure we stay within span(span_basis) numerically (and for rank bookkeeping)
    # The span restriction is: x ∈ span_basis iff x is in the column span.
    # Enforce by expressing K as span_basis * coeff (solve_linear_many later if you want),
    # but a cheap, robust way is just to intersect spans by rank filtering:
    # The intersection basis itself is handled elsewhere; here we just return K
    # and rely on invariance assumptions. Keep as-is to avoid extra solves.
    return independent_columns(K, p)


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
        coeff = _solve_linear(A, b, p)      # d×1
    except RuntimeError:
        if return_none:
            return None
        raise

    z = mod_p(span_basis @ coeff, p)        # n×1
    return z


def _select_module_generators_from_top_space(
    Fp: np.ndarray,
    Np: np.ndarray,
    top_candidates: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
) -> np.ndarray:
    """
    Convert a GF(p)-basis of the 'top quotient space' at length L into a set of
    *module generators* (one per indecomposable q^L-block).

    A vector v is accepted iff its cyclic submodule basis C(v) increases the current
    module-span by exactly deg_q*L dimensions.

    Robustification:
      - We may be handed an arbitrary basis of the top space; its columns need not
        individually be good cyclic generators. After a greedy pass over raw columns,
        we deterministically search linear combinations within the same top space.
    """
    top_candidates = independent_columns(mod_p(top_candidates, p), p)
    if top_candidates.shape[1] == 0:
        return top_candidates

    deg_q = int(deg_q)
    L = int(L)
    target = deg_q * L
    if target <= 0:
        return np.zeros((Fp.shape[0], 0), dtype=np.int64)

    # Expected number of indecomposable q^L-blocks in this top space:
    # top space dimension is typically deg_q * mult_L (over GF(p)).
    # If this isn't divisible, we fall back to "as many as we can certify".
    top_dim = int(top_candidates.shape[1])
    want = top_dim // deg_q if deg_q > 0 else 0

    # Local N^e mod p (avoid importing extra helpers)
    def _mat_pow_mod_local(A: np.ndarray, e: int) -> np.ndarray:
        e = int(e)
        n = A.shape[0]
        if e < 0:
            raise ValueError("negative exponent")
        if e == 0:
            return np.eye(n, dtype=np.int64)
        A = mod_p(A, p)
        res = np.eye(n, dtype=np.int64)
        base = A
        while e:
            if e & 1:
                res = mod_p(res @ base, p)
            base = mod_p(base @ base, p)
            e >>= 1
        return res

    Np = mod_p(Np, p)
    # Length-L sanity: v ∈ ker(N^L) but v ∉ ker(N^(L-1)) for L>1 (and for L=1 require N v = 0)
    NL = _mat_pow_mod_local(Np, L)
    NLm1 = _mat_pow_mod_local(Np, L - 1) if L - 1 >= 0 else np.eye(Np.shape[0], dtype=np.int64)

    span = np.zeros((Fp.shape[0], 0), dtype=np.int64)  # module-span accumulated so far
    gens: list[np.ndarray] = []

    def _try_accept(v: np.ndarray) -> bool:
        nonlocal span, gens
        v = mod_p(v, p)

        # Enforce the "length L" condition as a hard filter.
        # (This stops accidental shorter-length reps being accepted.)
        if np.any(mod_p(NL @ v, p) % p != 0):
            return False
        if L >= 2:
            if np.all(mod_p(NLm1 @ v, p) % p == 0):
                return False
        else:
            # L == 1: require v ∈ ker(N) (NL == N here).
            # Already enforced by NL @ v == 0 above.
            pass

        try:
            C = cyclic_submodule_basis(Fp, Np, v, int(deg_q), int(L), p)  # d × (deg_q*L)
        except Exception:
            return False

        C = independent_columns(mod_p(C, p), p)
        if C.shape[1] != target:
            return False

        new_span = independent_columns(np.concatenate([span, C], axis=1), p)
        if new_span.shape[1] == span.shape[1] + target:
            gens.append(v)
            span = new_span
            return True
        return False

    # --- Pass 1: greedy over raw basis columns
    for j in range(top_candidates.shape[1]):
        if want and len(gens) >= want:
            break
        _try_accept(top_candidates[:, j:j + 1])

    # --- Pass 2: if not enough, try deterministic linear combinations inside the top space
    # This is the crucial robustness step: the given basis may not contain clean generators.
    if want and len(gens) < want:
        # Use the same candidate set; deterministic nested loops.
        # Try v + a*u where u runs over earlier columns; a runs over 1..p-1.
        # (We do NOT use randomization here.)
        for j in range(top_candidates.shape[1]):
            if len(gens) >= want:
                break
            v = top_candidates[:, j:j + 1]
            if _try_accept(v):
                continue

            # Mix with other top-space vectors deterministically
            for i in range(top_candidates.shape[1]):
                if len(gens) >= want:
                    break
                if i == j:
                    continue
                u = top_candidates[:, i:i + 1]
                # coefficients in ascending order, skipping 0
                for a in range(1, int(p)):
                    cand = mod_p(v + a * u, p)
                    if _try_accept(cand):
                        break
                if len(gens) >= want:
                    break

    if not gens:
        return np.zeros((Fp.shape[0], 0), dtype=np.int64)
    return np.concatenate(gens, axis=1)


# def _select_module_generators_from_top_space(
#     Fp: np.ndarray,
#     Np: np.ndarray,
#     top_candidates: np.ndarray,
#     deg_q: int,
#     L: int,
#     p: int,
# ) -> np.ndarray:
#     """
#     Convert a GF(p)-basis of the 'top quotient space' at length L into a set of
#     *module generators* (one per indecomposable q^L-block).

#     A vector v is accepted iff its cyclic submodule basis C(v):
#       (i) has full expected dimension deg_q*L, and
#       (ii) increases the current module-span by exactly deg_q*L dimensions.
#     """
#     top_candidates = independent_columns(mod_p(top_candidates, p), p)
#     if top_candidates.shape[1] == 0:
#         return top_candidates

#     target = int(deg_q) * int(L)
#     gens: list[np.ndarray] = []
#     span = np.zeros((Fp.shape[0], 0), dtype=np.int64)

#     for j in range(top_candidates.shape[1]):
#         v = top_candidates[:, j:j + 1]

#         try:
#             C = cyclic_submodule_basis(Fp, Np, v, int(deg_q), int(L), p)  # d × (deg_q*L) nominally
#         except Exception:
#             # Not a valid length-L generator; skip deterministically.
#             continue

#         # IMPORTANT: enforce the *actual* generated module has full size target
#         C = independent_columns(mod_p(C, p), p)
#         if C.shape[1] != target:
#             continue

#         # Accept iff it adds exactly one whole module
#         new_span = independent_columns(np.concatenate([span, C], axis=1), p)
#         if new_span.shape[1] == span.shape[1] + target:
#             gens.append(v)
#             span = new_span

#     if not gens:
#         return np.zeros((Fp.shape[0], 0), dtype=np.int64)
#     return np.concatenate(gens, axis=1)
