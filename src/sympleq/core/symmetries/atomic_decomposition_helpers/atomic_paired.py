# sympleq/core/symmetries/atomic_paired.py
from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Tuple

from ..modular_helpers import (
    mod_p,
    independent_columns,
    omega_matrix,
    rank_mod,
    nullspace_mod,
    inv_mod_scalar,
)
from .atomic_types import AtomicBlock, AtomicInvariant
from .atomic_linear import (
    is_nondegenerate,
    darboux_basis_from_span,
    symplectic_orthogonal_complement_in_span,
)
from .module_invariants import (
    restrict_operator_invariant,
    q_of_F_restricted,
    jordan_chain_tops_nilpotent,
    cyclic_submodule_basis,
)
from .atomic_krylov import _select_module_generators_from_top_space


def _mat_pow_mod(A: np.ndarray, e: int, p: int) -> np.ndarray:
    """Matrix power A^e mod p (e>=0)."""
    e = int(e)
    n = A.shape[0]
    if e < 0:
        raise ValueError("_mat_pow_mod: e must be >= 0")
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


def _pairing_matrix_between(V_left: np.ndarray, V_right: np.ndarray, p: int) -> np.ndarray:
    """P = V_left^T Ω V_right in GF(p)."""
    n2 = V_left.shape[0]
    Ω = omega_matrix(n2 // 2, p)
    return mod_p(V_left.T @ Ω @ V_right, p)


def _intersection_basis(A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """Return an independent column basis for span(A) ∩ span(B) over GF(p)."""
    A = independent_columns(mod_p(A, p), p)
    B = independent_columns(mod_p(B, p), p)
    if A.size == 0 or B.size == 0:
        return np.zeros((A.shape[0], 0), dtype=np.int64)

    # Solve A x = B y  <=>  [A | -B] [x;y] = 0
    M = np.concatenate([A, mod_p(-B, p)], axis=1)
    N = nullspace_mod(M, p)  # (a+b)×k
    if N.size == 0:
        return np.zeros((A.shape[0], 0), dtype=np.int64)
    X = N[: A.shape[1], :]
    Id = mod_p(A @ X, p)
    return independent_columns(Id, p)


def _cyclic_module_has_full_rank(
    Fp: np.ndarray,
    Np: np.ndarray,
    v: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
) -> bool:
    """
    True iff cyclic_submodule_basis produces rank exactly deg_q*L (after independent filtering).
    """
    v = mod_p(v.reshape(-1, 1), p)
    target = int(deg_q) * int(L)
    try:
        C = cyclic_submodule_basis(Fp, Np, v, int(deg_q), int(L), p)
    except Exception:
        return False
    C = independent_columns(mod_p(C, p), p)
    return int(C.shape[1]) == target


def _candidate_generators_from_top_space(
    Fp: np.ndarray,
    Np: np.ndarray,
    top_basis: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
) -> np.ndarray:
    """
    Deterministic pool of candidate generators drawn from the provided top-space basis.
    Keeps ONLY vectors v such that cyclic_submodule_basis(Fp,Np,v) has full rank deg_q*L.
    """
    top_basis = independent_columns(mod_p(top_basis, p), p)
    if top_basis.shape[1] == 0:
        return top_basis

    cols: List[np.ndarray] = []
    for j in range(top_basis.shape[1]):
        v = top_basis[:, j:j + 1]
        if _cyclic_module_has_full_rank(Fp, Np, v, deg_q, L, p):
            cols.append(v)

    if not cols:
        return np.zeros((Fp.shape[0], 0), dtype=np.int64)
    return np.concatenate(cols, axis=1)


def _select_right_generators_with_full_pairing(
    *,
    NA: np.ndarray,                # (dimVq × m)
    P: np.ndarray,                 # (dimVq × dimVqs)
    right_pool: np.ndarray,        # (dimVqs × t)
    p: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Deterministically pick m columns W from right_pool such that M = NA^T P W is invertible.
    Returns (W, M).
    """
    NA = mod_p(NA, p)
    P = mod_p(P, p)
    right_pool = independent_columns(mod_p(right_pool, p), p)

    m = NA.shape[1]
    if m == 0:
        return np.zeros((right_pool.shape[0], 0), dtype=np.int64), np.zeros((0, 0), dtype=np.int64)
    if right_pool.shape[1] == 0:
        raise RuntimeError("Paired sector: empty right generator pool; cannot form dual pairing.")

    W_list: List[np.ndarray] = []
    M_cols: List[np.ndarray] = []
    rM = 0

    for j in range(right_pool.shape[1]):
        w = right_pool[:, j:j + 1]
        col = mod_p(NA.T @ (P @ w), p)  # (m×1)

        # rank test on the growing pairing matrix
        if not M_cols:
            M_try = col
        else:
            M_try = np.concatenate(M_cols + [col], axis=1)

        r_try = rank_mod(M_try, p)
        if r_try > rM:
            W_list.append(w)
            M_cols.append(col)
            rM = r_try
            if len(W_list) == m:
                break

    if len(W_list) != m:
        raise RuntimeError(
            f"Paired sector: could not select {m} right generators with full-rank pairing; "
            f"got {len(W_list)} (rank={rM})."
        )

    W = mod_p(np.concatenate(W_list, axis=1), p)          # (dimVqs × m)
    M = mod_p(np.concatenate(M_cols, axis=1), p)          # (m × m)
    if rank_mod(M, p) != m:
        raise RuntimeError("Paired sector: internal error, selected pairing matrix is not full rank.")
    return W, M


def atomic_blocks_in_paired_sector(
    F: np.ndarray,
    p: int,
    key: tuple[int, ...],
    key_star: tuple[int, ...],
    primaries: dict,
    allow_fallback: bool = False,
) -> tuple[list[AtomicBlock], AtomicInvariant]:
    """
    Atomic block construction for paired sector W = V_q ⊕ V_{q*} (q != q*).

    Certified ("OK") iff:
      - for every length L, the induced chain-level pairing between the chosen generators
        is non-singular (enforced by construction),
      - each constructed block span is nondegenerate,
      - and the blocks span the whole paired sector.
    """
    F = mod_p(F, p)

    q = primaries[key]["poly"]
    q_star = primaries[key_star]["poly"]

    Vq = independent_columns(mod_p(primaries[key]["V_basis"], p), p)
    Vqs = independent_columns(mod_p(primaries[key_star]["V_basis"], p), p)

    inv_data: Dict[str, Any] = {
        "status": "PENDING",
        "deg": int(primaries[key]["deg"]),
        "exponent": int(primaries[key]["exponent"]),
        "length_multiplicities": {},
        "pairing_rank": {},
        "checks_passed": [],
        "note": "",
    }

    if Vq.shape[1] == 0 or Vqs.shape[1] == 0:
        msg = "one side of paired sector is empty"
        if not allow_fallback:
            raise RuntimeError(f"Paired sector: {msg}.")
        inv_data["status"] = "DEGRADED"
        inv_data["note"] = msg
        inv = AtomicInvariant(sector_key=key, sector_type="paired", poly_key=key, data=inv_data)
        return [], inv

    # --- Certified-by-construction algorithm (extract-and-remove) ---
    # The previous "bulk" builder can produce invariant nondegenerate summands that are not
    # symplectically orthogonal to each other, which later makes the global basis non-symplectic.
    # Here we instead extract one atomic block at a time and remove its symplectic orthogonal
    # complement inside the paired sector. This guarantees pairwise orthogonality of blocks.

    n2 = F.shape[0]
    Ω_amb = omega_matrix(n2 // 2, p)

    max_exp = int(primaries[key]["exponent"])
    deg_q = int(primaries[key]["deg"])
    deg_qs = int(primaries[key_star]["deg"])
    if deg_q != deg_qs:
        raise RuntimeError("Paired sector: deg(q) != deg(q*) (unexpected).")

    # The full paired sector span W = Vq ⊕ Vq*
    W_sector = independent_columns(np.concatenate([Vq, Vqs], axis=1), p)
    if not is_nondegenerate(Ω_amb, W_sector, p):
        raise RuntimeError("Paired sector: sector span is degenerate (unexpected).")

    inv_data.setdefault("progress", [])  # list of dicts per extracted block

    blocks: List[AtomicBlock] = []

    try:
        rem = W_sector
        rem_dim = rank_mod(rem, p)

        # Safety: prevent infinite loops if something goes wrong.
        max_iters = rem_dim // (2 * deg_q) + 5
        it = 0

        while rem_dim > 0:
            it += 1
            if it > max_iters:
                raise RuntimeError("Paired sector: extraction stuck (too many iterations).")

            # Recompute the current left/right parts inside the remaining invariant subspace.
            Vq_r = _intersection_basis(Vq, rem, p)
            Vqs_r = _intersection_basis(Vqs, rem, p)
            if Vq_r.shape[1] == 0 or Vqs_r.shape[1] == 0:
                raise RuntimeError("Paired sector: remaining subspace lost one side (unexpected).")

            # Restrict to each side in its own coordinates.
            Fq = restrict_operator_invariant(F, Vq_r, p)
            Fqs = restrict_operator_invariant(F, Vqs_r, p)
            Nq = q_of_F_restricted(Fq, q, p)
            Nqs = q_of_F_restricted(Fqs, q_star, p)

            tops_left = jordan_chain_tops_nilpotent(Nq, max_exp, p)
            tops_right = jordan_chain_tops_nilpotent(Nqs, max_exp, p)
            if not tops_left and max_exp >= 1:
                tops_left = {1: np.eye(Fq.shape[0], dtype=np.int64)}
            if not tops_right and max_exp >= 1:
                tops_right = {1: np.eye(Fqs.shape[0], dtype=np.int64)}

            lengths = sorted(set(tops_left.keys()) | set(tops_right.keys()), reverse=True)
            if not lengths and max_exp >= 1:
                lengths = [1]

            # Pairing between current left/right bases.
            P = _pairing_matrix_between(Vq_r, Vqs_r, p)

            extracted = False
            last_err: str | None = None
            for L in lengths:
                A_raw = tops_left.get(L, np.zeros((Fq.shape[0], 0), dtype=np.int64))
                B_raw = tops_right.get(L, np.zeros((Fqs.shape[0], 0), dtype=np.int64))

                A = _select_module_generators_from_top_space(Fq, Nq, A_raw, deg_q, int(L), p)
                pool = _candidate_generators_from_top_space(Fqs, Nqs, B_raw, deg_q, int(L), p)

                # Bookkeeping only (not used by the algorithm)
                inv_data["length_multiplicities"][int(L)] = (int(A.shape[1]), int(pool.shape[1]))

                if A.shape[1] == 0 or pool.shape[1] == 0:
                    continue

                # Pick a single generator pair (v,w) with nonzero chain-level top pairing.
                N_pow = _mat_pow_mod(Nq, int(L) - 1, p)
                v_top = None
                w_top = None
                s_val = 0

                for ai in range(A.shape[1]):
                    v_cand = A[:, ai:ai + 1]
                    Nv = mod_p(N_pow @ v_cand, p)
                    for j in range(pool.shape[1]):
                        w_cand = pool[:, j:j + 1]
                        s = int(mod_p(Nv.T @ (P @ w_cand), p).reshape(())) % p
                        if s != 0:
                            v_top = v_cand
                            w_top = w_cand
                            s_val = s
                            break
                    if w_top is not None:
                        break

                if w_top is None or v_top is None:
                    last_err = f"no generator pair with nonzero top pairing at L={L}"
                    continue

                # Normalize so that Nv^T P w = 1.
                w_top = mod_p(w_top * inv_mod_scalar(s_val, p), p)

                # Build the cyclic submodules and lift to ambient.
                C_left = cyclic_submodule_basis(Fq, Nq, v_top, deg_q, int(L), p)
                C_right = cyclic_submodule_basis(Fqs, Nqs, w_top, deg_q, int(L), p)
                C_left = independent_columns(mod_p(C_left, p), p)
                C_right = independent_columns(mod_p(C_right, p), p)

                W_left = mod_p(Vq_r @ C_left, p)
                W_right = mod_p(Vqs_r @ C_right, p)
                span = independent_columns(np.concatenate([W_left, W_right], axis=1), p)

                if not is_nondegenerate(Ω_amb, span, p):
                    last_err = f"constructed span degenerate at L={L}"
                    continue

                # Ensure span is inside rem.
                if rank_mod(np.concatenate([rem, span], axis=1), p) != rem_dim:
                    last_err = f"constructed span not contained in remaining subspace at L={L}"
                    continue

                T_blk = darboux_basis_from_span(Ω_amb, span, p)
                blocks.append(
                    AtomicBlock(
                        T_blk=mod_p(T_blk, p),
                        half_dim=int(T_blk.shape[1] // 2),
                        sector_key=key,
                        inv=None,
                    )
                )

                # Remove its symplectic orthogonal complement within the remaining subspace.
                rem2 = symplectic_orthogonal_complement_in_span(Ω_amb, span, rem, p)
                rem2 = independent_columns(mod_p(rem2, p), p)
                rem2_dim = rank_mod(rem2, p)
                drop = rem_dim - rem2_dim
                if drop != span.shape[1]:
                    raise RuntimeError(
                        f"Paired sector: rank-drop mismatch when removing block (expected {span.shape[1]}, got {drop})."
                    )
                inv_data["progress"].append(
                    {
                        "L": int(L),
                        "block_dim": int(span.shape[1]),
                        "rem_dim_before": int(rem_dim),
                        "rem_dim_after": int(rem2_dim),
                    }
                )

                rem = rem2
                rem_dim = rem2_dim
                extracted = True
                break

            if not extracted:
                raise RuntimeError(
                    ("Paired sector: could not extract a valid block from remaining subspace" +
                     f" (last_err={last_err})" if last_err else "")
                )

        # Final span check
        all_cols = np.concatenate([b.T_blk for b in blocks], axis=1) if blocks else np.zeros((n2, 0), dtype=np.int64)
        if rank_mod(np.concatenate([W_sector, all_cols], axis=1), p) != rank_mod(W_sector, p):
            raise RuntimeError("Paired sector: some constructed block columns lie outside W = Vq ⊕ Vq*.")
        if rank_mod(all_cols, p) != rank_mod(W_sector, p):
            raise RuntimeError("Paired sector: constructed blocks do not span the paired sector.")

        inv_data["checks_passed"].append("extract-and-remove guarantees orthogonality")
        inv_data["checks_passed"].append("each block nondegenerate")
        inv_data["checks_passed"].append("blocks span paired sector")
        inv_data["status"] = "OK"
        sector_cost = max((int(b.half_dim) for b in blocks), default=0)
        inv_data["cost_certificate"] = {
            "lower_bound": int(sector_cost),
            "attained": True,
            "complete": True,
            "sector_cost": int(sector_cost),
            "note": (
                "paired sector certified by quotient-level dual-basis construction "
                "and extract/remove orthogonalization"
            ),
        }

    except Exception as e:
        if not allow_fallback:
            raise
        inv_data["status"] = "DEGRADED"
        inv_data["note"] = f"fallback: {type(e).__name__}: {e}"

        # Single-block fallback spanning the whole paired sector.
        T_blk = darboux_basis_from_span(Ω_amb, W_sector, p)
        blocks = [AtomicBlock(mod_p(T_blk, p), int(T_blk.shape[1] // 2), key, None)]

    inv = AtomicInvariant(sector_key=key, sector_type="paired", poly_key=key, data=inv_data)
    blocks = [AtomicBlock(b.T_blk, b.half_dim, b.sector_key, inv) for b in blocks]

    # Sanity: blocks must be a direct sum (no overlap)
    if blocks:
        all_cols = np.concatenate([b.T_blk for b in blocks], axis=1)
        if rank_mod(all_cols, p) != sum(b.T_blk.shape[1] for b in blocks):
            raise RuntimeError("Paired sector: produced blocks overlap (not direct sum).")

    return blocks, inv
