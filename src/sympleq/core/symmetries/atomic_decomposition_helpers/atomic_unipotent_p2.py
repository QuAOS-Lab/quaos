# sympleq/core/symmetries/atomic_unipotent_p2.py
from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Tuple

from ..modular_helpers import mod_p, independent_columns, rank_mod, omega_matrix
from .atomic_types import AtomicBlock, AtomicInvariant
from .atomic_linear import (
    restrict_operator,
    is_nondegenerate,
    darboux_basis_from_span,
    mat_pow_mod,
    kernel_in_span,
    symplectic_orthogonal_complement_in_span,
)
from .module_invariants import (
    jordan_chain_tops_nilpotent,
    cyclic_submodule_basis,
    q_of_F_restricted,
)
from .atomic_krylov import _select_module_generators_from_top_space


# ----------------------------
# Small local helpers (p=2)
# ----------------------------

def _basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    """Deterministically pick 'want' columns from 'candidates' extending span(base)."""
    base = independent_columns(mod_p(base, p), p) if base.size else base
    picked = np.zeros((candidates.shape[0], 0), dtype=np.int64)
    r_base = rank_mod(base, p) if base.size else 0

    for j in range(candidates.shape[1]):
        c = candidates[:, j:j + 1]
        r_try = rank_mod(np.concatenate([base, picked, c], axis=1), p)
        if r_try > r_base + picked.shape[1]:
            picked = np.concatenate([picked, c], axis=1)
            if picked.shape[1] == want:
                return picked

    raise RuntimeError("_basis_extend: could not extend by required amount.")


def jordan_chain_tops_nilpotent_in_span(
    N: np.ndarray, space_basis: np.ndarray, max_exp: int, p: int
) -> Dict[int, np.ndarray]:
    """
    Restricted analogue:
      tops[L] represents K_L / (K_{L-1} + N K_{L+1}) inside span(space_basis),
    with K_j := ker(N^j) ∩ span(space_basis).

    Deterministic: uses _basis_extend.
    """
    N = mod_p(N, p)
    space_basis = independent_columns(mod_p(space_basis, p), p)
    d = N.shape[0]
    if space_basis.shape[1] == 0:
        return {}

    # K[0]=0, K[1..max_exp], and K[max_exp+1]=K[max_exp]
    K: List[np.ndarray] = [np.zeros((d, 0), dtype=np.int64)]
    for j in range(1, int(max_exp) + 1):
        Kj = kernel_in_span(mat_pow_mod(N, j, p), space_basis, p)
        K.append(independent_columns(Kj, p))
    K.append(K[int(max_exp)])

    tops: Dict[int, np.ndarray] = {}
    for L in range(1, int(max_exp) + 1):
        KL = K[L]
        if KL.shape[1] == 0:
            continue

        S = K[L - 1]
        NKLp1 = mod_p(N @ K[L + 1], p) if K[L + 1].shape[1] else np.zeros((d, 0), dtype=np.int64)
        if NKLp1.shape[1]:
            S = np.concatenate([S, NKLp1], axis=1) if S.shape[1] else NKLp1
        S = independent_columns(S, p) if S.shape[1] else S

        rS = rank_mod(S, p) if S.shape[1] else 0
        need = KL.shape[1] - rS
        if need <= 0:
            continue

        chosen = _basis_extend(S, KL, need, p)
        tops[L] = chosen

    return tops


def _pair_value(v: np.ndarray, w: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int, p: int) -> int:
    """Compute <v, N^{L-1} w> mod p (scalar)."""
    v = mod_p(v.reshape(-1, 1), p)
    w = mod_p(w.reshape(-1, 1), p)
    Nr = np.eye(N.shape[0], dtype=np.int64) if (L - 1) == 0 else mat_pow_mod(N, L - 1, p)
    return int(mod_p(v.T @ Omega @ (Nr @ w), p).reshape(()))


def _find_partner_in_top_span_p2(
    v: np.ndarray, A: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int
) -> np.ndarray | None:
    """
    p=2: find w among columns of A such that <v, N^{L-1} w> = 1 and rank([v,w])=2.
    Returns None if impossible (v is in the radical of this pairing on span(A)).
    Deterministic: picks the first suitable column.
    """
    p = 2
    v = mod_p(v.reshape(-1, 1), p)
    A = independent_columns(mod_p(A, p), p)
    if A.shape[1] == 0:
        return None

    Nr = np.eye(N.shape[0], dtype=np.int64) if (L - 1) == 0 else mat_pow_mod(N, L - 1, p)
    r = mod_p(v.T @ Omega @ (Nr @ A), p).reshape(-1)  # length m

    nz = np.where(r % 2 != 0)[0]
    if nz.size == 0:
        return None

    for j in nz:
        w = A[:, int(j):int(j) + 1]
        # exclude w == v
        if np.array_equal(mod_p(w, p), mod_p(v, p)):
            continue
        # ensure independent in the top quotient (rank 2 for {v,w})
        if rank_mod(np.concatenate([v, w], axis=1), p) != 2:
            continue
        return w

    return None

def _beta_from_top_generators_p2(
    gens: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int
) -> int:
    """
    Chapter-5-style beta on the 'top quotient' for this extracted indecomposable,
    computed from the quadratic invariant evaluated on a symplectic top basis.

    (p=2): the label is detected at mid-chain, not at N^{L-1}.
      - V_beta(2k):     use q(v) = <v, N^{k-1} v> = <v, N^{L/2 - 1} v>
      - W_beta(2l+1):   use q(v) = <v, N^{l}   v> = <v, N^{(L-1)/2} v>
      - W(even): label should come out 0 in “honest” cases; we still compute it.
    """
    p = 2
    gens = mod_p(gens, p)
    if gens.shape[1] == 0:
        return 0

    # Mid-chain exponent for the quadratic refinement (Chapter 5, p=2 case)
    if L % 2 == 0:
        # L = 2k  -> exponent k-1
        exp = max(L // 2 - 1, 0)
    else:
        # L = 2l+1 -> exponent l
        exp = (L - 1) // 2

    Nr = np.eye(N.shape[0], dtype=np.int64) if exp == 0 else mat_pow_mod(N, exp, p)

    def q(vcol: np.ndarray) -> int:
        vcol = mod_p(vcol.reshape(-1, 1), p)
        return int(mod_p(vcol.T @ Omega @ (Nr @ vcol), p).reshape(()))

    if gens.shape[1] == 1:
        return q(gens[:, 0]) & 1

    if gens.shape[1] == 2:
        qv = q(gens[:, 0]) & 1
        qw = q(gens[:, 1]) & 1
        return (qv & qw)  # product in GF(2)

    return 0



def _p2_hyperbolic_pairs_from_alternating_form(B: np.ndarray, p: int = 2) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Deterministic Gram–Schmidt for an alternating form B over GF(2).

    Returns (e_list, f_list, rad_list), where each entry is a *coefficient vector*
    (length m) in the original coordinate basis such that:
      - e_i^T B f_i = 1, and e_i,f_i are mutually orthogonal to all other pairs
      - rad_list spans the radical {v: v^T B = 0}
    """
    if p != 2:
        raise ValueError("_p2_hyperbolic_pairs_from_alternating_form is p=2 only")
    B = mod_p(np.asarray(B, dtype=np.int64), 2)
    m = B.shape[0]
    R = [np.eye(m, dtype=np.int64)[:, i] for i in range(m)]  # remaining vectors (standard basis), as 1D arrays

    e_list: List[np.ndarray] = []
    f_list: List[np.ndarray] = []
    rad_list: List[np.ndarray] = []

    def pair(u: np.ndarray, v: np.ndarray) -> int:
        return int((u.reshape(1, -1) @ B @ v.reshape(-1, 1)) % 2)

    # Work on a mutable list of remaining vectors, orthogonalizing as we go
    while R:
        a = R.pop(0)
        # find b with <a,b>=1
        found = None
        for idx, cand in enumerate(R):
            if pair(a, cand) == 1:
                found = idx
                break

        if found is None:
            # a is radical (w.r.t. current remaining set, which is already orthogonal to previous pairs)
            rad_list.append(a.copy())
            continue

        b = R.pop(found)

        # Orthogonalize remaining vectors to the new pair (a,b)
        newR = []
        for x in R:
            ax = pair(x, b)  # <x,b>
            bx = pair(x, a)  # <x,a>
            if ax:
                x = (x ^ a)  # x <- x + <x,b> a
            if bx:
                x = (x ^ b)  # x <- x + <x,a> b
            newR.append(x)
        R = newR

        e_list.append(a.copy())
        f_list.append(b.copy())

    return e_list, f_list, rad_list


def _p2_length_form_invariants(A_top, Omega, N, L):
    p = 2
    A_top = independent_columns(mod_p(A_top, p), p)
    m = int(A_top.shape[1])
    if m == 0:
        return {
            "top_dim": 0, "B_rank": 0, "rad_dim": 0,
            "B_sym_ok": True, "B_alt_ok": True,
            "q_defined": False, "q_mid_polar_ok": None, "q_witness_polar_ok": True,
            "arf": None,
        }

    Nr = np.eye(N.shape[0], dtype=np.int64) if (L - 1) == 0 else mat_pow_mod(N, L - 1, p)
    B = mod_p(A_top.T @ Omega @ (Nr @ A_top), p)
    B_rank = int(rank_mod(B, p))
    rad_dim = int(m - B_rank)

    B_sym_ok = bool(np.array_equal(B, B.T))
    B_alt_ok = bool(np.all(np.diag(B) % 2 == 0))

    # "mid" q: the one you were *trying* to use (debug only)
    q_defined = bool(L > 1)
    q_mid_polar_ok = None
    if q_defined:
        if L % 2 == 0:
            exp = L // 2          # IMPORTANT: don't use k-1 if you intend L=2 to be meaningful
        else:
            exp = (L - 1) // 2
        Nm = np.eye(N.shape[0], dtype=np.int64) if exp == 0 else mat_pow_mod(N, exp, p)
        Q_mid = mod_p(A_top.T @ Omega @ (Nm @ A_top), p)
        Bq_mid = mod_p(Q_mid + Q_mid.T, p)
        q_mid_polar_ok = bool(np.array_equal(Bq_mid, B))

    # "witness" quadratic refinement: always exists when B is alternating.
    # Choose Q_wit as strict upper-triangular part of B so that Q_wit + Q_wit^T = B (diag=0).
    q_witness_polar_ok = bool(B_alt_ok)
    # (We don't compute Arf from this witness: it isn't canonical.)

    return {
        "top_dim": int(m),
        "B_rank": int(B_rank),
        "rad_dim": int(rad_dim),
        "B_sym_ok": bool(B_sym_ok),
        "B_alt_ok": bool(B_alt_ok),
        "q_defined": bool(q_defined),
        "q_mid_polar_ok": q_mid_polar_ok,
        "q_witness_polar_ok": bool(q_witness_polar_ok),
        "arf": None,  # leave None unless/until you implement the *canonical* quadratic
    }



# ----------------------------
# Main builder (p=2 unipotent self sector)
# ----------------------------

def atomic_blocks_in_unipotent_self_sector_p2(
    F: np.ndarray, key: Tuple[int, ...], primaries: dict
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    p=2, self sector with q(x)=x±1 (unipotent in bad characteristic).

    Construct atomic indecomposables in the Chapter-5 sense:
      - Try to peel off V_beta(2k) blocks when a single cyclic module is nondegenerate (L even).
      - Otherwise peel off W-like blocks by pairing two cyclic modules of the same length L.
      - Attach 'beta' from the quadratic top invariant q_L on the corresponding top quotient.

    Returns blocks + invariant record. Marks status="OK" iff blocks span the sector and
    each extracted block span is nondegenerate in the sector symplectic form.
    """
    p = 2
    F = mod_p(F, p)

    q = primaries[key]["poly"]
    deg_q = int(primaries[key]["deg"])       # should be 1 here
    max_exp0 = int(primaries[key]["exponent"])

    V = independent_columns(mod_p(primaries[key]["V_basis"], p), p)
    if V.shape[1] == 0:
        inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data={"status": "empty"})
        return [], inv

    # Canonicalize sector basis to Darboux so Ω becomes standard in sector coords
    n2 = F.shape[0]
    Ω_amb = omega_matrix(n2 // 2, p)
    if not is_nondegenerate(Ω_amb, V, p):
        raise RuntimeError("Unipotent self sector basis is degenerate; cannot proceed.")

    T_sec = darboux_basis_from_span(Ω_amb, V, p)  # (2n × 2m)
    m2 = T_sec.shape[1]
    if m2 % 2 != 0:
        raise RuntimeError("Unipotent self sector dimension must be even.")
    m = m2 // 2

    F_sec = restrict_operator(F, T_sec, p)  # (2m × 2m)
    Ω = omega_matrix(m, p)

    # Nilpotent N = q(F) on this primary (for x±1 over p=2 this is F+I)
    N = q_of_F_restricted(F_sec, q, p)
    max_exp = max_exp0

    # -------------------------
    # Invariant summary (full sector)
    # -------------------------
    tops_full = jordan_chain_tops_nilpotent(N, max_exp, p)

    # Conjugacy-invariant kernel profile for nilpotent N (dims of ker N^k)
    kernel_profile = [
        int((2 * m) - rank_mod(mat_pow_mod(N, k, p), p)) for k in range(1, int(max_exp0) + 1)
    ]

    length_summary: Dict[int, Dict[str, Any]] = {}
    length_invariants: Dict[int, Dict[str, Any]] = {}

    for L, Araw in tops_full.items():
        Araw = independent_columns(mod_p(Araw, p), p)
        invL = _p2_length_form_invariants(Araw, Ω, N, int(L))
        length_invariants[int(L)] = invL

        # Deterministic selection of genuine module generators (guardrail).
        A = _select_module_generators_from_top_space(F_sec, N, Araw, deg_q, int(L), p)

        q1_count = 0
        for j in range(A.shape[1]):
            q1_count += (_beta_from_top_generators_p2(A[:, j:j + 1], Ω, N, int(L)) & 1)

        length_summary[int(L)] = {
            "mult": int(invL["top_dim"]),
            "gen_dim": int(A.shape[1]),
            "gen_dim_ok": bool(int(A.shape[1]) == int(invL["top_dim"])),
            "q1_count": int(q1_count),
            **invL,
        }

    # -------------------------
    # Deterministic atomic extraction (sector coords) (sector coords)
    # -------------------------
    blocks: List[AtomicBlock] = []
    blocks_meta: List[Dict[str, Any]] = []
    built_cols_sec: List[np.ndarray] = []

    space_basis = np.eye(2 * m, dtype=np.int64)

    while space_basis.shape[1] > 0:
        tops = jordan_chain_tops_nilpotent_in_span(N, space_basis, max_exp, p)
        if not tops:
            break

        L = max(tops.keys())
        A_full = independent_columns(mod_p(tops[L], p), p)  # full top space basis (quotient reps)
        A = _select_module_generators_from_top_space(F_sec, N, A_full, deg_q, int(L), p)

        if A.shape[1] == 0:
            max_exp = int(L) - 1
            continue

        progressed = False

        # Try candidates in a deterministic order until we peel off one indecomposable.
        for j in range(A.shape[1]):
            v_top = A[:, j:j + 1]

            # Build cyclic module Cv (deg=1 => expected dim L)
            try:
                Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, int(L), p)
            except Exception:
                continue
            Cv = independent_columns(mod_p(Cv, p), p)
            if Cv.shape[1] == 0:
                continue

            # --- Attempt V_beta(2k): single cyclic module already nondegenerate (requires L even) ---
            took_V = False
            if (int(L) % 2 == 0) and (Cv.shape[1] % 2 == 0) and is_nondegenerate(Ω, Cv, p):
                T_blk = darboux_basis_from_span(Ω, Cv, p)  # sector coords
                built_cols_sec.append(T_blk)

                rem = symplectic_orthogonal_complement_in_span(Ω, T_blk, space_basis, p)

                beta = _beta_from_top_generators_p2(v_top, Ω, N, int(L))
                # Chapter-5 naming: V_beta(2k) is the “single Jordan block” indecomposable
                typ = "V_alpha" if beta == 1 else "V"

                T_blk_amb = mod_p(T_sec @ T_blk, p)
                blocks.append(AtomicBlock(T_blk=T_blk_amb, half_dim=int(T_blk_amb.shape[1] // 2), sector_key=key, inv=None))
                blocks_meta.append({"type": typ, "L": int(L), "half_dim": int(T_blk_amb.shape[1] // 2), "beta": int(beta)})

                space_basis = rem
                progressed = True
                took_V = True

            if took_V:
                break

            # --- Otherwise attempt W-like block: pair with a partner in the top space ---
            w_top = _find_partner_in_top_span_p2(v_top, A_full, Ω, N, int(L))

            if w_top is None:
                # v_top is “unpairable” at this L in the current top space; try next candidate.
                continue

            try:
                Cw = cyclic_submodule_basis(F_sec, N, w_top, deg_q, int(L), p)
            except Exception:
                continue
            Cw = independent_columns(mod_p(Cw, p), p)

            span = independent_columns(np.concatenate([Cv, Cw], axis=1), p)
            if span.shape[1] % 2 != 0 or not is_nondegenerate(Ω, span, p):
                continue

            T_blk = darboux_basis_from_span(Ω, span, p)  # sector coords
            built_cols_sec.append(T_blk)
            rem = symplectic_orthogonal_complement_in_span(Ω, T_blk, space_basis, p)

            # beta from the 2D top space of this indecomposable
            gens = np.concatenate([v_top, w_top], axis=1)
            beta = _beta_from_top_generators_p2(gens, Ω, N, int(L))

            # Chapter-5 naming: W_alpha exists only for odd L >= 3; otherwise treat as W
            if (int(L) % 2 == 1) and (int(L) >= 3) and beta == 1:
                typ = "W_alpha"
            else:
                typ = "W"

            T_blk_amb = mod_p(T_sec @ T_blk, p)
            blocks.append(AtomicBlock(T_blk=T_blk_amb, half_dim=int(T_blk_amb.shape[1] // 2), sector_key=key, inv=None))
            blocks_meta.append({"type": typ, "L": int(L), "half_dim": int(T_blk_amb.shape[1] // 2), "beta": int(beta)})

            space_basis = rem
            progressed = True
            break

        if not progressed:
            # If we get stuck at this L, deterministically drop L and continue.
            # This prevents infinite loops in pathological p=2 cases.
            max_exp = int(L) - 1
            if max_exp <= 0:
                break

    # Sanity: blocks should span the whole sector in sector coordinates
    dim_sector = 2 * m
    if built_cols_sec:
        all_cols_sec = np.concatenate(built_cols_sec, axis=1)
        dim_blocks = rank_mod(all_cols_sec, p)
    else:
        dim_blocks = 0

    if dim_blocks != dim_sector:
        raise RuntimeError(
            f"Unipotent p=2 self sector: blocks do not span sector "
            f"(dim_blocks={dim_blocks}, dim_sector={dim_sector}). "
            f"key={key}, deg={deg_q}, exp0={max_exp0}"
        )

    if built_cols_sec:
        all_cols_sec = np.concatenate(built_cols_sec, axis=1)
        if rank_mod(all_cols_sec, p) != all_cols_sec.shape[1]:
            raise RuntimeError("Unipotent p=2 self sector: extracted block bases overlap (not a direct sum).")

    # Kernel profile cross-check: sum of block-restricted profiles must match the whole-sector profile.
    N_blks = [restrict_operator(N, T_blk, p) for T_blk in built_cols_sec] if built_cols_sec else []
    kernel_profile_blocks: List[int] = []
    for k in range(1, int(max_exp0) + 1):
        tot = 0
        for Nb in N_blks:
            Nk = mat_pow_mod(Nb, k, p)
            tot += int(Nb.shape[0] - rank_mod(Nk, p))
        kernel_profile_blocks.append(int(tot))

    if kernel_profile_blocks != kernel_profile:
        raise RuntimeError(
            "Unipotent p=2 self sector: kernel profile mismatch (block sum != sector). "
            f"key={key}, profile={kernel_profile}, block_profile={kernel_profile_blocks}"
        )

    # Beta bookkeeping: counts by length and type.
    beta_counts_by_L: Dict[int, Dict[str, int]] = {}
    for bm in blocks_meta:
        L_raw = bm.get("L")
        if L_raw is None:
            raise RuntimeError("Missing 'L' in unipotent block metadata.")
        L = int(L_raw)

        beta_raw = bm.get("beta", 0)
        beta = int(0 if beta_raw is None else beta_raw) & 1

        typ_raw = bm.get("type", "")
        typ = "" if typ_raw is None else str(typ_raw)
        d = beta_counts_by_L.setdefault(L, {"V0": 0, "V1": 0, "W0": 0, "W1": 0})
        if typ.startswith("V"):
            d["V1" if beta else "V0"] += 1
        else:
            d["W1" if (typ.endswith("alpha") or beta) else "W0"] += 1

    inv_data: Dict[str, Any] = {
        "status": "OK",
        "deg": int(deg_q),
        "exponent": int(max_exp0),
        "length_summary": length_summary,
        "blocks": blocks_meta,
        "p2_unipotent": {
            "kernel_profile": kernel_profile,
            "kernel_profile_blocks": kernel_profile_blocks,
            "length_invariants": length_invariants,          # per-L top-space invariants (Arf, ranks, ...)
            "beta_counts_by_L": beta_counts_by_L,            # per-L V/W beta splits as extracted
        },
    }

    inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data=inv_data)

    blocks = [AtomicBlock(b.T_blk, b.half_dim, key, inv) for b in blocks]
    return blocks, inv
