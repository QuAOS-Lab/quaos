# sympleq/core/symmetries/atomic_paired.py
from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Tuple

from .modular_helpers import mod_p, independent_columns, omega_matrix, inv_mod_mat, rank_mod
from .atomic_types import AtomicBlock, AtomicInvariant
from .atomic_linear import is_nondegenerate, darboux_basis_from_span
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


def _cyclic_module_has_full_rank(
    Fp: np.ndarray,
    Np: np.ndarray,
    v: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
) -> bool:
    """
    True iff cyclic_submodule_basis produces rank exactly deg_q*L (after indep filtering).
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
    allow_fallback: bool = True,
) -> tuple[list[AtomicBlock], AtomicInvariant]:
    """
    Atomic block construction for paired sector W = V_q ⊕ V_{q*} (q != q*).

    Certified ("OK") iff:
      - for every length L, the induced chain-level pairing between the chosen generators
        is nonsingular (enforced by construction),
      - each constructed block span is nondegenerate,
      - and the blocks span the whole paired sector.
    """
    F = mod_p(F, p)

    q = primaries[key]["poly"]
    qstar = primaries[key_star]["poly"]

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
        inv_data["status"] = "DEGRADED"
        inv_data["note"] = "one side of paired sector is empty"
        inv = AtomicInvariant(sector_key=key, sector_type="paired", poly_key=key, data=inv_data)
        return [], inv

    # Restrict F to each invariant subspace (coordinates in the given bases)
    Fq = restrict_operator_invariant(F, Vq, p)
    Fqs = restrict_operator_invariant(F, Vqs, p)

    # Nilpotent parts
    Nq = q_of_F_restricted(Fq, q, p)
    Nqs = q_of_F_restricted(Fqs, qstar, p)

    max_exp = int(primaries[key]["exponent"])
    deg_q = int(primaries[key]["deg"])
    deg_qs = int(primaries[key_star]["deg"])
    if deg_q != deg_qs:
        raise RuntimeError("Paired sector: deg(q) != deg(q*) (unexpected).")

    # Jordan chain-top representatives by length on each side
    tops_left = jordan_chain_tops_nilpotent(Nq, max_exp, p)
    tops_right = jordan_chain_tops_nilpotent(Nqs, max_exp, p)

    # Semisimple fallback (exp=1) sometimes yields {}, treat as tops at L=1 spanning whole space.
    if not tops_left and max_exp >= 1:
        tops_left = {1: np.eye(Fq.shape[0], dtype=np.int64)}
    if not tops_right and max_exp >= 1:
        tops_right = {1: np.eye(Fqs.shape[0], dtype=np.int64)}

    lengths = sorted(set(tops_left.keys()) | set(tops_right.keys()))
    if not lengths and max_exp >= 1:
        lengths = [1]

    # Ambient pairing between the two primary bases
    P = _pairing_matrix_between(Vq, Vqs, p)

    n2 = F.shape[0]
    Ω_amb = omega_matrix(n2 // 2, p)

    blocks: List[AtomicBlock] = []
    built_cols: List[np.ndarray] = []

    try:
        for L in lengths:
            A_raw = tops_left.get(L,  np.zeros((Fq.shape[0], 0), dtype=np.int64))
            B_raw = tops_right.get(L, np.zeros((Fqs.shape[0], 0), dtype=np.int64))

            # Select left module generators (one per indecomposable q^L block)
            A = _select_module_generators_from_top_space(Fq, Nq, A_raw, deg_q, int(L), p)
            a_mult = int(A.shape[1])

            # Build a pool of right-side candidate generators FROM THE RIGHT TOP SPACE
            # (not forced to match left yet).
            pool = _candidate_generators_from_top_space(Fqs, Nqs, B_raw, deg_q, int(L), p)

            inv_data["length_multiplicities"][int(L)] = (a_mult, int(pool.shape[1]))

            if a_mult == 0:
                continue
            if pool.shape[1] == 0:
                raise RuntimeError(f"Paired sector: no valid right generators at length L={L}.")

            # Chain-level left transform
            Npow = _mat_pow_mod(Nq, int(L) - 1, p)
            NA = mod_p(Npow @ A, p)  # (dimVq × a_mult)

            # Choose right generators with full-rank pairing, then dualize them to get identity pairing
            W, M = _select_right_generators_with_full_pairing(NA=NA, P=P, right_pool=pool, p=p)
            rM = rank_mod(M, p)
            inv_data["pairing_rank"][int(L)] = int(rM)
            if rM != a_mult:
                raise RuntimeError(
                    f"Paired sector: chain-level top pairing singular at length {L}. "
                    f"rank={rM}, expected={a_mult}."
                )

            Minv = inv_mod_mat(M, p)
            B_dual = mod_p(W @ Minv, p)  # ensures NA^T P B_dual = I

            # Now build one atomic block per paired top
            for j in range(a_mult):
                v_top = A[:, j:j + 1]        # coords in Vq basis
                w_top = B_dual[:, j:j + 1]   # coords in Vqs basis

                # Cyclic submodule bases in sector coordinates
                C_left = cyclic_submodule_basis(Fq,  Nq,  v_top, deg_q, int(L), p)
                C_right = cyclic_submodule_basis(Fqs, Nqs, w_top, deg_q, int(L), p)

                C_left = independent_columns(mod_p(C_left, p), p)
                C_right = independent_columns(mod_p(C_right, p), p)

                # Lift to ambient
                W_left = mod_p(Vq @ C_left, p)
                W_right = mod_p(Vqs @ C_right, p)

                span = independent_columns(np.concatenate([W_left, W_right], axis=1), p)

                if not is_nondegenerate(Ω_amb, span, p):
                    raise RuntimeError(f"Paired sector: constructed atomic span is degenerate (L={L}, j={j}).")

                T_blk = darboux_basis_from_span(Ω_amb, span, p)
                blocks.append(AtomicBlock(T_blk=mod_p(T_blk, p), half_dim=int(T_blk.shape[1] // 2), sector_key=key, inv=None))
                built_cols.append(T_blk)

        # Global paired-sector span check: blocks should span W = Vq ⊕ Vqs and lie inside it.
        W_sector = independent_columns(np.concatenate([Vq, Vqs], axis=1), p)
        dim_sector = rank_mod(W_sector, p)

        if built_cols:
            all_cols = np.concatenate(built_cols, axis=1)
            if rank_mod(np.concatenate([W_sector, all_cols], axis=1), p) != dim_sector:
                raise RuntimeError("Paired sector: some constructed block columns lie outside W = Vq ⊕ Vq*.")
            dim_blocks = rank_mod(all_cols, p)
        else:
            dim_blocks = 0

        if dim_blocks != dim_sector:
            raise RuntimeError(
                f"Paired sector: constructed blocks do not span the paired sector subspace "
                f"(dim_blocks={dim_blocks}, dim_sector={dim_sector})."
            )

        inv_data["checks_passed"].append("per-length chain pairing nonsingular (constructed)")
        inv_data["checks_passed"].append("each block nondegenerate")
        inv_data["checks_passed"].append("blocks span paired sector")
        inv_data["status"] = "OK"

    except Exception as e:
        if not allow_fallback:
            # In certified mode, fail loudly.
            raise
        inv_data["status"] = "DEGRADED"
        inv_data["note"] = f"fallback: {type(e).__name__}: {e}"

        W_sector = independent_columns(np.concatenate([Vq, Vqs], axis=1), p)
        if not is_nondegenerate(Ω_amb, W_sector, p):
            raise RuntimeError("Paired sector fallback: sector span is degenerate (unexpected).") from e
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
