# sympleq/core/symmetries/rcf_prepass.py
from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple

from sympleq.core.symmetries.modular_helpers import mod_p, independent_columns, rank_mod
from sympleq.core.symmetries.polynomials_fp import (
    poly_monic, poly_pow, poly_divmod, poly_xgcd, poly_is_zero, poly_reciprocal, poly_is_zero, poly_mul
)
from sympleq.core.symmetries.minpoly import minimal_polynomial, factor_poly_over_fp


def _is_x_pm_1(q: np.ndarray, p: int) -> bool:
    q = poly_monic(q, p)
    if q.size != 2 or int(q[1]) % p != 1:
        return False
    a0 = int(q[0]) % p
    # x - 1  => [p-1, 1]
    # x + 1  => [1, 1]
    return (a0 == 1 % p) or (a0 == (p - 1) % p)


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



def poly_eval_matrix(F: np.ndarray, poly: np.ndarray, p: int) -> np.ndarray:
    """
    Evaluate poly(F) for column-action convention, coeffs low->high.
    """
    n2 = F.shape[0]
    poly = mod_p(poly, p)
    M = np.zeros((n2, n2), dtype=np.int64)
    P = np.eye(n2, dtype=np.int64)
    for a in poly:
        a = int(a) % p
        if a:
            M = mod_p(M + a * P, p)
        P = mod_p(P @ F, p)
    return M


def primary_components_crt(F: np.ndarray, p: int) -> Dict:
    """
    Provably-correct primary decomposition via CRT projectors.
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

    # CRT projector for each modulus fi:
    # Let Mi = mF/fi. Find s,t with s*Mi + t*fi = 1. Then Pi = s(F) Mi(F).
    for (q, e), fi in zip(factors, moduli):
        Mi, r = poly_divmod(mF, fi, p)
        if not poly_is_zero(r):
            raise RuntimeError("CRT setup failed: fi does not divide mF (should not happen).")

        s, t, g = poly_xgcd(Mi, fi, p)   # s*Mi + t*fi = g
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
        }

    # Group into symplectic sectors by reciprocity
    for key, data in primaries.items():
        q = data["poly"]
        q_star = poly_reciprocal(q, p)
        data["reciprocal_key"] = tuple(q_star.tolist())
        data["self_reciprocal"] = (data["reciprocal_key"] == key)

    used = set()
    sectors = []
    for key, data in primaries.items():
        if key in used:
            continue
        if data["self_reciprocal"]:
            deg = int(data["deg"])
            e = int(data["exponent"])
            if _is_x_pm_1(data["poly"], p):
                half_floor = ((e + 1) // 2) * deg
            else:
                half_floor = e * deg
            W = independent_columns(data["V_basis"], p)
            sectors.append({"type": "self", "key": key, "W_basis": W, "half_dim_floor": int(half_floor)})
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
                deg = int(data["deg"])
                e = int(data["exponent"])
                half_floor = e * deg

                sectors.append({"type": "paired", "key": key, "key_star": k_star, "W_basis": W,
                                "half_dim_floor": int(half_floor)})
                used.add(key)
                used.add(k_star)

    # Sanity: sectors should span V (columns total rank = n2)
    all_cols = np.concatenate([sec["W_basis"] for sec in sectors], axis=1) if sectors else np.zeros((n2, 0),
                                                                                                    dtype=np.int64)
    if rank_mod(all_cols, p) != n2:
        raise RuntimeError("Primary sectorization failed: sectors do not span V.")

    return {"mF": mF, "primaries": primaries, "sectors": sectors, "factors": factors}


def rcf_prepass(F: np.ndarray, p: int) -> Dict:
    return primary_components_crt(F, p)


if __name__ == "__main__":
    # Minimal sanity: identity => only one primary (x-1), self sector spans whole space.
    p = 2
    F = np.eye(6, dtype=np.int64)
    meta = rcf_prepass(F, p)
    assert len(meta["sectors"]) == 1
    W = meta["sectors"][0]["W_basis"]
    assert W.shape[0] == 6
    print("rcf_prepass.py tests passed")
