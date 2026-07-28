from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple

from sympleq.core.symmetries.modular_helpers import mod_p, independent_columns, rank_mod, omega_matrix
from sympleq.core.symmetries.polynomials_fp import (
    poly_monic, poly_pow, poly_divmod, extended_euclidean, poly_is_zero, poly_reciprocal, poly_mul, poly_eval_matrix
)
from sympleq.core.symmetries.minpoly import minimal_polynomial, factor_poly_over_fp
from .mqf_types import PrepassContext, SectorContext
from .mqf_linear import darboux_basis_from_span, restrict_operator, symplectic_left_inverse
from .module_invariants import q_of_F_restricted
from .cost_bounds import is_x_pm_1 as _is_x_pm_1, compute_sector_profile


def _half_dim_floor_from_minpoly_factor(q: np.ndarray, e: int, p: int) -> int:
    """
    A *safe* (possibly loose) lower bound on the half-dimension of any
    nontrivial symplectic invariant block that can live inside this primary.

    IMPORTANT:
      - For p=2 and q=x+1 (unipotent), minimal polynomial data alone is NOT
        sufficient to certify the true minimum block size. We therefore return 1
        (meaning: no certified >1 lower bound from this prepass alone).
      - For paired sectors (q != q*), the usual e*deg bound is safe.
      - For self-reciprocal q != x±1, e*deg is a reasonable bound.
      - For x±1 in odd p, ceil(e/2)*deg appears; for p=2, treat as uncertified.
    """
    q = poly_monic(q, p)
    deg = len(q) - 1

    if p == 2 and _is_x_pm_1(q, p):
        # This is the hard unipotent/self-reciprocal corner; do not "certify" from minpoly.
        return 1

    if _is_x_pm_1(q, p):
        # classical floor for ±1 sectors when p is odd (or when you accept this as heuristic)
        return ((int(e) + 1) // 2) * deg

    return int(e) * deg


def _compute_Lmin(sectors: list[dict], n: int) -> int:
    """
    Global lower bound on the optimal qudit cost.
    We deliberately keep this certified given the information rcf_prepass computes.
    """
    if not sectors:
        return 1
    L = max(int(sec.get("half_dim_floor", 1)) for sec in sectors)
    # Always at least 1 (identity cost convention), at most n.
    return max(1, min(int(L), int(n)))


def _normalize_factorization_output(factors, p: int) -> List[Tuple[np.ndarray, int]]:
    """
    Normalize factorization output into [(q, e)] with monic q and total exponent e.
    Accepts:
      A) [q, q, r, ...]                   (repeated factors)
      B) [(q, e), (q, e2), (r, f), ...]   (pairs; possibly repeated)
    Returns:
      [(q1, e1), (q2, e2), ...] with distinct qi.
    """
    mult: Dict[Tuple[int, ...], int] = {}
    reps: Dict[Tuple[int, ...], np.ndarray] = {}

    for item in factors:
        if isinstance(item, (tuple, list)) and len(item) == 2 and isinstance(item[1], (int, np.integer)):
            q, e = item[0], int(item[1])
            q = poly_monic(q, p)
            key = tuple(q.tolist())
            reps[key] = q
            mult[key] = mult.get(key, 0) + e   # NOTE: sum, not max
        else:
            q = poly_monic(item, p)
            key = tuple(q.tolist())
            reps[key] = q
            mult[key] = mult.get(key, 0) + 1

    return [(reps[k], mult[k]) for k in reps.keys()]


def _safe_sector_coordinates(
    F: np.ndarray,
    W: np.ndarray,
    Omega: np.ndarray,
    p: int,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray, np.ndarray | None, str]:
    """
    Build a deterministic sector coordinate system.

    If W is a nondegenerate symplectic sector, T_sec is returned as a Darboux
    basis, so T_sec.T @ Omega @ T_sec is the standard symplectic form.  If W is
    empty or unexpectedly degenerate, the raw independent W basis is returned
    and the error is recorded in the final string instead of hiding it.

    Returns (T_sec, T_sec_leftinv, Omega_sec, F_sec, note).
    """
    F = mod_p(F, p)
    W = independent_columns(mod_p(W, p), p)
    if W.shape[1] == 0:
        return W, None, np.zeros((0, 0), dtype=np.int64), None, "empty sector basis"

    try:
        Gram = mod_p(W.T @ Omega @ W, p)
        if W.shape[1] % 2 != 0 or rank_mod(Gram, p) != W.shape[1]:
            note = "sector span is not nondegenerate; using raw sector basis"
            return W, None, Gram, None, note

        T_sec = darboux_basis_from_span(Omega, W, p)
        T_sec_leftinv = symplectic_left_inverse(T_sec, p)
        Omega_sec = mod_p(T_sec.T @ Omega @ T_sec, p)
        F_sec = restrict_operator(F, T_sec, p)
        return T_sec, T_sec_leftinv, Omega_sec, F_sec, ""
    except Exception as exc:
        Omega_sec = mod_p(W.T @ Omega @ W, p)
        note = f"failed to construct Darboux sector coordinates: {type(exc).__name__}: {exc}"
        return W, None, Omega_sec, None, note


def _build_prepass_context(
    *,
    F: np.ndarray,
    p: int,
    n: int,
    primaries: Dict[Tuple[int, ...], Dict],
    sectors: list[dict],
    mF: np.ndarray,
    factors: List[Tuple[np.ndarray, int]],
    Lmin_star: int,
) -> PrepassContext:
    """Build typed sector contexts from primary-sector records."""
    Omega = omega_matrix(n, p)
    sector_contexts: list[SectorContext] = []

    for index, sec in enumerate(sectors):
        sector_type = "paired" if sec.get("type") == "paired" else "self"
        key = tuple(sec["key"])
        key_star = tuple(sec["key_star"]) if sec.get("key_star") is not None else None
        poly_key = key

        W = independent_columns(mod_p(sec.get("W_basis", np.zeros((2 * n, 0), dtype=np.int64)), p), p)
        T_sec, T_sec_leftinv, Omega_sec, F_sec, coord_note = _safe_sector_coordinates(F, W, Omega, p)

        N_sec = None
        if sector_type == "self" and F_sec is not None:
            q = primaries[key]["poly"]
            N_sec = q_of_F_restricted(F_sec, q, p)

        ctx_meta: Dict[str, object] = {
            "index": int(index),
            "W_basis": W,
            "primaries": primaries,
            "primary_info": primaries.get(key),
            "dim2": int(W.shape[1]),
            "half_dim": int(W.shape[1] // 2),
            "half_dim_floor": int(sec.get("half_dim_floor", 1)),
            "floor_certified": bool(sec.get("floor_certified", True)),
            "sector_note": str(sec.get("note", "")),
            "mF": mF,
            "factors": factors,
            "Lmin_star": int(Lmin_star),
            "coordinate_note": coord_note,
        }
        if sector_type == "paired" and key_star is not None:
            ctx_meta["primary_info_star"] = primaries.get(key_star)

        ctx = SectorContext(
            sector_key=key,
            sector_type=sector_type,
            poly_key=poly_key,
            sector_key_star=key_star,
            p=int(p),
            deg_q=int(sec.get("deg", primaries.get(key, {}).get("deg", 0))),
            max_exp=int(sec.get("exponent", primaries.get(key, {}).get("exponent", 0))),
            T_sec=T_sec,
            T_sec_leftinv=T_sec_leftinv,
            F_sec=F_sec,
            Omega_sec=Omega_sec,
            N_sec=N_sec,
            meta=ctx_meta,
        )

        # Compute the invariant-derived (search-independent) cost lower bound
        # before any extraction, and stash it on the context meta.
        profile = compute_sector_profile(F, ctx, p)
        ctx_meta["lengths_present"] = profile.get("lengths_present", [])
        ctx_meta["p2_diag_flags"] = profile.get("p2_diag_flags", {})
        ctx_meta["cost_lower_bound"] = profile.get("lower_bound")
        ctx_meta["cost_lower_bound_complete"] = bool(profile.get("complete", False))

        sector_contexts.append(ctx)

    # Prefer the exact invariant bound, obtained by taking the maximum over
    # independent sector bounds, over the loose minpoly floor when every sector
    # produced a certified bound.
    exact_lbs: list[int] = []
    for ctx in sector_contexts:
        lower_bound = (ctx.meta or {}).get("cost_lower_bound")
        if not isinstance(lower_bound, int):
            exact_lbs = []
            break
        exact_lbs.append(lower_bound)

    if exact_lbs:
        Lmin_exact = max(1, min(max(exact_lbs), int(n)))
    else:
        Lmin_exact = int(Lmin_star)

    return PrepassContext(
        p=int(p),
        n=int(n),
        Omega=Omega,
        sectors=sector_contexts,
        meta={
            "mF": mF,
            "factors": factors,
            "Lmin_star": int(Lmin_exact),
            "Lmin_star_minpoly_floor": int(Lmin_star),
            "primaries": primaries,
        },
    )


def primary_components_crt(F: np.ndarray, p: int) -> Dict:
    """
    Build primary components using coprime-polynomial projectors.
    """
    F = mod_p(F, p)
    n2 = F.shape[0]

    mF = minimal_polynomial(F, p)
    raw = factor_poly_over_fp(mF, p)
    factors = _normalize_factorization_output(raw, p)

    # Build the list of pairwise coprime moduli f_i = q_i^{e_i}
    moduli: List[np.ndarray] = [poly_pow(q, e, p) for (q, e) in factors]

    # Sanity: product of moduli should equal mF up to monic scaling
    prod = np.array([1], dtype=np.int64)
    for fi in moduli:
        # multiply in poly space
        prod = poly_mul(prod, fi, p)
    prod = poly_monic(prod, p)
    if not np.array_equal(prod, poly_monic(mF, p)):
        raise RuntimeError(
            "Factorization mismatch: product of (q^e) != mF. "
            "Your factor_poly_over_fp output is inconsistent."
        )

    primaries: Dict[Tuple[int, ...], Dict] = {}

    # Project onto each primary subspace V_i = ker(fi(F)) using a Bezout
    # coefficient for the complementary factor Mi = mF / fi.
    for (q, e), fi in zip(factors, moduli):
        Mi, r = poly_divmod(mF, fi, p)
        if not poly_is_zero(r):
            raise RuntimeError("CRT setup failed: fi does not divide mF (should not happen).")

        s, t, g = extended_euclidean(Mi, fi, p)   # s*Mi + t*fi = g
        # g must be 1 (monic constant)
        if not (len(g) == 1 and int(g[0]) % p == 1):
            raise RuntimeError("CRT xgcd failed: Mi and fi not coprime.")

        Pi = mod_p(poly_eval_matrix(F, s, p) @ poly_eval_matrix(F, Mi, p), p)
        V_basis = independent_columns(Pi, p)

        key = tuple(poly_monic(q, p).tolist())
        primaries[key] = {
            "poly": poly_monic(q, p),
            "deg": len(q) - 1,
            "exponent": int(e),
            "V_basis": V_basis,
            "dim": int(V_basis.shape[1]),
        }

    # Group into symplectic sectors by reciprocity
    for key, data in primaries.items():
        q = data["poly"]
        q_star = poly_reciprocal(q, p)
        data["reciprocal_key"] = tuple(q_star.tolist())
        data["self_reciprocal"] = (data["reciprocal_key"] == key)

    # Build sector bases W by combining V_basis of reciprocal pairs,
    #  then independent_columns to clean up any linear dependencies.
    used = set()
    sectors = []
    for key, data in primaries.items():
        if key in used:
            continue
        if data["self_reciprocal"]:
            e = int(data["exponent"])
            half_floor = _half_dim_floor_from_minpoly_factor(data["poly"], e, p)
            W = independent_columns(data["V_basis"], p)
            sectors.append({
                "type": "self",
                "key": key,
                "W_basis": W,
                "deg": int(data["deg"]),
                "exponent": int(e),
                "dim2": int(W.shape[1]),
                "half_dim": int(W.shape[1] // 2),
                "half_dim_floor": int(half_floor),
                "floor_certified": bool(not (p == 2 and _is_x_pm_1(data["poly"], p))),
                "note": ("p=2, q=x±1: minpoly does not certify block size"
                         if (p == 2 and _is_x_pm_1(data["poly"], p)) else "")
            })
            used.add(key)
        else:
            k_star = data["reciprocal_key"]
            if k_star not in primaries:
                # Should not happen; be robust
                W = independent_columns(data["V_basis"], p)
                sectors.append({"type": "self", "key": key, "W_basis": W})
                used.add(key)
            else:
                W = np.concatenate([data["V_basis"], primaries[k_star]["V_basis"]], axis=1)
                W = independent_columns(W, p)
                e = int(data["exponent"])
                half_floor = _half_dim_floor_from_minpoly_factor(data["poly"], e, p)  # for paired: returns e*deg
                sectors.append({
                    "type": "paired",
                    "key": key,
                    "key_star": k_star,
                    "W_basis": W,
                    "deg": int(data["deg"]),
                    "exponent": int(e),
                    "dim2": int(W.shape[1]),
                    "half_dim": int(W.shape[1] // 2),
                    "half_dim_floor": int(half_floor),
                    "floor_certified": True,
                    "note": ""
                })
                used.add(key)
                used.add(k_star)

    # Sanity: sectors should span V (columns total rank = n2)
    all_cols = np.concatenate([sec["W_basis"] for sec in sectors], axis=1) if sectors else np.zeros((n2, 0),
                                                                                                    dtype=np.int64)
    if rank_mod(all_cols, p) != n2:
        raise RuntimeError("Primary sectorization failed: sectors do not span V.")

    n = n2 // 2
    Lmin_star = _compute_Lmin(sectors, n)

    prepass_context = _build_prepass_context(
        F=F,
        p=p,
        n=n,
        primaries=primaries,
        sectors=sectors,
        mF=mF,
        factors=factors,
        Lmin_star=Lmin_star,
    )

    # The prepass context recomputes an exact Lmin from invariant sector bounds;
    # surface it, falling back to the minpoly floor.
    Lmin_exact = int(prepass_context.meta.get("Lmin_star", Lmin_star))

    return {
        "p": int(p),
        "n": int(n),
        "mF": mF,
        "primaries": primaries,
        "sector_contexts": prepass_context.sectors,
        "prepass_context": prepass_context,
        "factors": factors,
        "Lmin_star": int(Lmin_exact),
        "Lmin_star_minpoly_floor": int(Lmin_star),
    }


def rcf_prepass(F: np.ndarray, p: int) -> Dict:
    return primary_components_crt(F, p)
