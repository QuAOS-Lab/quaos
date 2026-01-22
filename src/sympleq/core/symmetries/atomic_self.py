# sympleq/core/symmetries/atomic_self.py
from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Tuple

from .modular_helpers import mod_p, independent_columns, rank_mod, _solve_linear, omega_matrix
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
    q_of_F_restricted,
    jordan_chain_tops_nilpotent,   # full-space tops (used for invariant summary)
    cyclic_submodule_basis,
)

from .atomic_krylov import _select_module_generators_from_top_space


def _is_alternating(B: np.ndarray, p: int) -> bool:
    """
    Alternating bilinear form matrix test over GF(p).
    For odd p: B^T = -B and diag=0.
    For p=2:  -B = B, so "alternating" reduces to symmetric with diag=0.
    """
    B = mod_p(B, p)
    if np.any(np.diag(B) % p != 0):
        return False
    return np.array_equal(B.T % p, (-B) % p)


def _basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    """
    Deterministically pick 'want' columns from 'candidates' that extend span(base).
    Returns picked columns (ambient coords).
    """
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


def jordan_chain_tops_nilpotent_in_span(N: np.ndarray, space_basis: np.ndarray, max_exp: int, p: int) -> Dict[int, np.ndarray]:
    """
    Restricted analogue of jordan_chain_tops_nilpotent:
      tops[L] columns represent K_L / (K_{L-1} + N K_{L+1}) inside span(space_basis),
    where K_j := ker(N^j) ∩ span(space_basis).

    Deterministic: uses the same extension rule as module_invariants.
    """
    N = mod_p(N, p)
    space_basis = independent_columns(mod_p(space_basis, p), p)
    d = N.shape[0]
    if space_basis.shape[1] == 0:
        return {}

    # K[0] = 0, K[j] for j=1..max_exp, and K[max_exp+1] := K[max_exp]
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


def _induced_form_matrix_on_tops(A: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int, p: int) -> np.ndarray:
    """
    Induced (quotient) form on L-tops:
        B_L(u,v) = <u, N^{L-1} v>  (mod p)
    for columns u,v of A.
    """
    if L <= 0:
        raise ValueError("L must be >= 1")
    if A.shape[1] == 0:
        return np.zeros((0, 0), dtype=np.int64)
    Nr = np.eye(N.shape[0], dtype=np.int64) if (L - 1) == 0 else mat_pow_mod(N, L - 1, p)
    return mod_p(A.T @ Omega @ (Nr @ A), p)


def _pair_value(v: np.ndarray, w: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int, p: int) -> int:
    """
    Compute <v, N^{L-1} w> mod p (scalar).
    """
    Nr = np.eye(N.shape[0], dtype=np.int64) if (L - 1) == 0 else mat_pow_mod(N, L - 1, p)
    return int(mod_p(v.T @ Omega @ (Nr @ w), p).reshape(()))


def _find_partner_in_top_span(
    v: np.ndarray, A: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int, p: int
) -> np.ndarray:
    """
    Given v (a top vector) and A whose columns span the top space at length L,
    find w in span(A) such that <v, N^{L-1} w> = 1.

    Deterministic:
      - try basis columns with scaling,
      - else solve a 1×m linear system.
    """
    v = mod_p(v.reshape(-1, 1), p)
    A = independent_columns(mod_p(A, p), p)
    if A.shape[1] == 0:
        raise RuntimeError("_find_partner_in_top_span: empty top space.")

    Nr = np.eye(N.shape[0], dtype=np.int64) if (L - 1) == 0 else mat_pow_mod(N, L - 1, p)
    r = mod_p(v.T @ Omega @ (Nr @ A), p)  # 1×m

    # fast: pick first nonzero entry and scale that basis vector
    nz = np.where(r.reshape(-1) % p != 0)[0]
    if nz.size:
        j = int(nz[0])
        a = int(r[0, j]) % p
        inv = pow(a, p - 2, p) if p != 2 else 1  # in p=2, a=1 anyway
        return mod_p(A[:, j:j + 1] * inv, p)

    # else: solve r * c = 1
    Aeq = mod_p(r, p)  # 1×m
    beq = np.array([[1]], dtype=np.int64)
    c = _solve_linear(Aeq, beq, p)  # m×1
    w = mod_p(A @ c, p)
    # sanity
    if _pair_value(v, w, Omega, N, L, p) % p != 1 % p:
        raise RuntimeError("_find_partner_in_top_span: failed to build partner with pairing=1.")
    return w

def atomic_blocks_in_self_sector_nonunipotent(
    F: np.ndarray, p: int, key: Tuple[int, ...], primaries: dict
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Self-reciprocal sector q=q* with q != x±1.

    Pipeline:
      1) take ambient sector span V = primaries[key]["V_basis"]
      2) canonicalize to a Darboux basis T (so restricted Ω is standard)
      3) restrict F to sector coordinates, compute N=q(F) nilpotent
      4) compute invariant summaries from the induced forms on top spaces
      5) deterministically extract atomic blocks by largest-L top pairing,
         removing each block by symplectic orthogonal complement.
      6) lift each block back to ambient via T

    Notes:
      - For p=2 this routine is still heuristic; if block extraction fails,
        we degrade gracefully to a single sector-sized block so callers can proceed.
    """
    F = mod_p(F, p)

    q = primaries[key]["poly"]
    deg_q = int(primaries[key]["deg"])
    max_exp_orig = int(primaries[key]["exponent"])   # keep original for reporting
    exp_work = int(max_exp_orig)                     # may be decreased locally during extraction

    V = independent_columns(mod_p(primaries[key]["V_basis"], p), p)
    if V.shape[1] == 0:
        inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data={"status": "empty"})
        return [], inv

    # Canonicalize sector basis to Darboux so Ω becomes standard in sector coords
    n2 = F.shape[0]
    Ω_amb = omega_matrix(n2 // 2, p)
    if not is_nondegenerate(Ω_amb, V, p):
        raise RuntimeError("Self sector basis is degenerate; cannot proceed.")

    T_sec = darboux_basis_from_span(Ω_amb, V, p)  # (2n × 2m)
    m2 = T_sec.shape[1]
    if m2 % 2 != 0:
        raise RuntimeError("Self sector dimension must be even.")
    m = m2 // 2

    F_sec = restrict_operator(F, T_sec, p)  # (2m × 2m), in canonical symplectic coords
    Ω = omega_matrix(m, p)                  # standard Ω on sector

    N = q_of_F_restricted(F_sec, q, p)      # nilpotent on this primary

    # -------------------------
    # Invariant summary (full sector)
    # -------------------------
    tops_full = jordan_chain_tops_nilpotent(N, max_exp_orig, p)  # Dict[L] -> (2m × mult_L)

    top_multiplicities: Dict[int, int] = {}
    top_form: Dict[int, Dict[str, Any]] = {}

    for L, A in tops_full.items():
        A = independent_columns(mod_p(A, p), p)
        top_multiplicities[int(L)] = int(A.shape[1])
        B = _induced_form_matrix_on_tops(A, Ω, N, int(L), p)
        rB = rank_mod(B, p)
        top_form[int(L)] = {
            "rank": int(rB),
            "alternating": bool(_is_alternating(B, p)),
        }

    # -------------------------
    # Deterministic block extraction in sector coords
    # -------------------------
    blocks_meta: List[Dict[str, Any]] = []
    blocks: List[AtomicBlock] = []

    # keep accepted block bases in sector coordinates for span sanity-check
    built_cols_sec: List[np.ndarray] = []

    space_basis = np.eye(2 * m, dtype=np.int64)

    extraction_error: Exception | None = None
    try:
        while space_basis.shape[1] > 0:
            # tops in current space (use exp_work, which may shrink deterministically)
            if exp_work <= 0:
                break
            tops = jordan_chain_tops_nilpotent_in_span(N, space_basis, exp_work, p)
            if not tops:
                break

            L = max(tops.keys())

            A_raw = independent_columns(mod_p(tops[L], p), p)
            A = _select_module_generators_from_top_space(F_sec, N, A_raw, deg_q, int(L), p)

            if A.shape[1] == 0:
                # nothing usable at this L; deterministically drop this length
                exp_work = int(L) - 1
                continue

            v_top = A[:, 0:1]

            # Candidate 1: self-dual cyclic module
            bvv = _pair_value(v_top, v_top, Ω, N, int(L), p)

            Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, int(L), p)  # (2m × deg*L)
            Cv = independent_columns(mod_p(Cv, p), p)

            if Cv.shape[1] > 0 and (Cv.shape[1] % 2 == 0) and is_nondegenerate(Ω, Cv, p):
                T_blk = darboux_basis_from_span(Ω, Cv, p)  # SECTOR COORDS
                built_cols_sec.append(T_blk)

                rem = symplectic_orthogonal_complement_in_span(Ω, T_blk, space_basis, p)

                T_blk_amb = mod_p(T_sec @ T_blk, p)
                blocks.append(
                    AtomicBlock(
                        T_blk=T_blk_amb,
                        half_dim=int(T_blk_amb.shape[1] // 2),
                        sector_key=key,
                        inv=None,
                    )
                )
                blocks_meta.append(
                    {"type": "self", "L": int(L), "half_dim": int(T_blk_amb.shape[1] // 2), "bvv": int(bvv)}
                )
                space_basis = rem
                continue

            # Candidate 2: hyperbolic pairing block from two cyclic modules
            w_top = _find_partner_in_top_span(v_top, A, Ω, N, int(L), p)
            Cw = cyclic_submodule_basis(F_sec, N, w_top, deg_q, int(L), p)
            Cw = independent_columns(mod_p(Cw, p), p)

            span = independent_columns(np.concatenate([Cv, Cw], axis=1), p)
            if span.shape[1] % 2 != 0 or not is_nondegenerate(Ω, span, p):
                found = False
                for j in range(1, A.shape[1]):
                    cand = A[:, j:j + 1]
                    if _pair_value(v_top, cand, Ω, N, int(L), p) % p == 0:
                        continue
                    a = _pair_value(v_top, cand, Ω, N, int(L), p) % p
                    inva = pow(int(a), p - 2, p) if p != 2 else 1
                    cand = mod_p(cand * inva, p)
                    Ccand = cyclic_submodule_basis(F_sec, N, cand, deg_q, int(L), p)
                    Ccand = independent_columns(mod_p(Ccand, p), p)
                    span2 = independent_columns(np.concatenate([Cv, Ccand], axis=1), p)
                    if span2.shape[1] % 2 == 0 and is_nondegenerate(Ω, span2, p):
                        span = span2
                        w_top = cand
                        found = True
                        break
                if not found:
                    raise RuntimeError("Self sector: failed to form a nondegenerate hyperbolic block from top pairing.")

            T_blk = darboux_basis_from_span(Ω, span, p)  # SECTOR COORDS
            built_cols_sec.append(T_blk)

            rem = symplectic_orthogonal_complement_in_span(Ω, T_blk, space_basis, p)

            T_blk_amb = mod_p(T_sec @ T_blk, p)
            blocks.append(
                AtomicBlock(
                    T_blk=T_blk_amb,
                    half_dim=int(T_blk_amb.shape[1] // 2),
                    sector_key=key,
                    inv=None,
                )
            )
            blocks_meta.append(
                {"type": "hyperbolic", "L": int(L), "half_dim": int(T_blk_amb.shape[1] // 2), "bvv": int(bvv)}
            )
            space_basis = rem

    except Exception as e:
        extraction_error = e

    # If extraction failed in p=2, degrade gracefully to a single sector-sized block.
    if extraction_error is not None and p == 2:
        blocks = [AtomicBlock(T_blk=mod_p(T_sec, p), half_dim=m, sector_key=key, inv=None)]
        blocks_meta = [{
            "type": "fallback_sector",
            "half_dim": int(m),
            "note": f"{type(extraction_error).__name__}: {extraction_error}",
        }]
        built_cols_sec = []  # skip span-check in degraded mode

    # Span sanity-check (only if we did not degrade)
    if extraction_error is None:
        dim_sector = 2 * m
        if built_cols_sec:
            all_cols_sec = np.concatenate(built_cols_sec, axis=1)
            dim_blocks = rank_mod(all_cols_sec, p)
        else:
            dim_blocks = 0

        if dim_blocks != dim_sector:
            raise RuntimeError(
                f"Self sector: blocks do not span sector "
                f"(dim_blocks={dim_blocks}, dim_sector={dim_sector}). "
                f"key={key}, deg={deg_q}, exp={max_exp_orig}"
            )

    # Status / invariant payload
    if p != 2:
        status = "OK"
    else:
        status = "heuristic_p2" if extraction_error is None else "DEGRADED_p2"

    inv_data: Dict[str, Any] = {
        "status": status,
        "deg": int(deg_q),
        "exponent": int(max_exp_orig),          # report the true primary exponent
        "exp_work_final": int(exp_work),        # optional: helps debugging
        "top_multiplicities": top_multiplicities,
        "top_form": top_form,
        "blocks": blocks_meta,
    }
    if extraction_error is not None:
        inv_data["note"] = f"{type(extraction_error).__name__}: {extraction_error}"

    inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data=inv_data)
    blocks = [AtomicBlock(b.T_blk, b.half_dim, key, inv) for b in blocks]
    return blocks, inv
