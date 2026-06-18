# sympleq/core/symmetries/atomic_decomposition_helpers/cost_bounds.py
"""
Invariant-derived (search-independent) cost lower bounds.

This is the Phase-2 ``decoupled cost certificate``: the minimal achievable
half-dimension of a sector is computed from conjugacy invariants alone
(Theorems 3.3 and 4.13 of the review notes), *before* any block extraction, so
that ``certified_minimal`` becomes the genuine statement

    lower_bound (from invariants)  ==  attained (verified construction)

rather than being true by construction.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from ..polynomials_fp import poly_monic
from .atomic_filtration import build_nilpotent_filtration
from .module_invariants import restrict_operator_invariant, q_of_F_restricted


def is_x_pm_1(q: np.ndarray, p: int) -> bool:
    """True iff monic ``q`` is ``x-1`` or ``x+1`` over GF(p)."""
    q = poly_monic(np.asarray(q), p)
    if q.size != 2 or int(q[1]) % p != 1:
        return False
    a0 = int(q[0]) % p
    return (a0 == 1 % p) or (a0 == (p - 1) % p)


def sector_lower_bound(
    sector_type: str,
    q: np.ndarray,
    p: int,
    deg_q: int,
    lengths_present: List[int],
    p2_diag_flags: Optional[Dict[int, bool]] = None,
) -> int:
    """
    Exact minimal half-dimension achievable on this sector.

    ``lengths_present`` : chain lengths L of N = q(F) on the sector (as
        F_p[x]/(q)-module lengths, i.e. the set of L with T_L != 0).
    ``p2_diag_flags`` : for p=2, q=x+1 only -- {even L: diag(b_L) != 0}
        (the negation of ``B_alt_ok`` from ``_p2_length_form_invariants``).

    Justification: Thm. 3.3 (paired / Hermitian / odd-p) and Thm. 4.13 (p=2
    unipotent) of the accompanying notes.
    """
    d = int(deg_q)
    Lmax = max(lengths_present, default=0)
    if Lmax == 0:
        return 0

    if sector_type == "paired":
        return d * Lmax

    if not is_x_pm_1(q, p):  # self-reciprocal, q != x +/- 1 (deg even)
        return (d * Lmax) // 2

    if p != 2:  # q = x +/- 1, odd p
        odd = max((L for L in lengths_present if L % 2 == 1), default=0)
        even = max((L // 2 for L in lengths_present if L % 2 == 0), default=0)
        return max(1, odd, even)

    # p = 2, q = x + 1: Theorem 4.13
    flags = p2_diag_flags or {}
    c = 1
    for L in lengths_present:
        if L % 2 == 1:
            c = max(c, L)
        else:
            c = max(c, (L // 2) if flags.get(L, False) else L)
    return c


def _nilpotent_on_sector(F: np.ndarray, ctx, p: int) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Return (N, Omega_sec) for the sector's nilpotent N = q(F).

    For self sectors this is the cached ``ctx.N_sec`` / ``ctx.Omega_sec``.  For
    paired sectors we form N on the V_q half from the primary basis (Omega is
    not needed there -- the paired bound uses only chain lengths).
    """
    q = np.asarray(ctx.poly_key, dtype=np.int64)
    if ctx.sector_type == "self":
        return ctx.N_sec, ctx.Omega_sec

    # paired: build N on the V_q primary component
    meta = ctx.meta or {}
    prim = meta.get("primary_info") or {}
    Vb = prim.get("V_basis")
    if Vb is None or np.asarray(Vb).size == 0:
        return None, None
    Fr = restrict_operator_invariant(F, np.asarray(Vb, dtype=np.int64), p)
    N = q_of_F_restricted(Fr, q, p)
    return N, None


def compute_sector_profile(F: np.ndarray, ctx, p: int) -> Dict[str, object]:
    """
    Compute, from invariants alone, the data feeding ``sector_lower_bound``:
    the set of chain lengths present and (for p=2, q=x+1) the per-length
    diagonal flags.  Returns a dict with keys ``lengths_present``,
    ``p2_diag_flags``, ``lower_bound`` (int or None), and ``complete`` (bool).
    """
    q = np.asarray(ctx.poly_key, dtype=np.int64)
    deg_q = int(ctx.deg_q)
    try:
        N, Omega_sec = _nilpotent_on_sector(F, ctx, p)
        if N is None or np.asarray(N).size == 0:
            # Empty/trivial sector: no chains, lower bound 0 is exact.
            return {"lengths_present": [], "p2_diag_flags": {}, "lower_bound": 0, "complete": True}

        N = np.asarray(N, dtype=np.int64)
        dim = N.shape[0]
        max_exp = int(ctx.max_exp) if int(ctx.max_exp) > 0 else dim
        filt = build_nilpotent_filtration(N, np.eye(dim, dtype=np.int64), max_exp, p)

        lengths: List[int] = [
            L for L in range(1, max_exp + 1)
            if (filt.tops.get(L) is not None and filt.tops[L].shape[1] > 0)
        ]

        flags: Dict[int, bool] = {}
        if p == 2 and ctx.sector_type == "self" and is_x_pm_1(q, p) and Omega_sec is not None:
            # Lazy import avoids a heavy import at prepass load time and any cycle.
            from .atomic_unipotent_p2 import _p2_length_form_invariants
            for L in lengths:
                if L % 2 == 0:
                    info = _p2_length_form_invariants(filt.tops[L], Omega_sec, N, int(L))
                    flags[L] = (not bool(info.get("B_alt_ok", True)))

        lb = sector_lower_bound(ctx.sector_type, q, p, deg_q, lengths, flags)
        return {"lengths_present": lengths, "p2_diag_flags": flags, "lower_bound": int(lb), "complete": True}
    except Exception as exc:  # never let bound computation break the prepass
        return {
            "lengths_present": [],
            "p2_diag_flags": {},
            "lower_bound": None,
            "complete": False,
            "error": f"{type(exc).__name__}: {exc}",
        }
