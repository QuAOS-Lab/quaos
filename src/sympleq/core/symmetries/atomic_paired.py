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


def atomic_blocks_in_paired_sector(
    F: np.ndarray,
    p: int,
    key: tuple[int, ...],
    key_star: tuple[int, ...],
    primaries: dict,
    allow_fallback: bool = True
) -> tuple[list[AtomicBlock], AtomicInvariant]:
    """
    Atomic block construction for paired sector W = V_q ⊕ V_{q*} (q != q*).

    Certified ("OK") iff:
      - all per-length pairings are non-singular using the correct chain-level pairing,
      - all constructed spans are nondegenerate,
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

    # Nilpotent parts: N = q(F) on V_q and N* = q*(F) on V_{q*}
    Nq = q_of_F_restricted(Fq, q, p)
    Nqs = q_of_F_restricted(Fqs, qstar, p)

    max_exp = int(primaries[key]["exponent"])
    deg_q = int(primaries[key]["deg"])
    deg_qs = int(primaries[key_star]["deg"])
    if deg_q != deg_qs:
        raise RuntimeError("Paired sector: deg(q) != deg(q*) (unexpected).")

    # Jordan chain-top representatives by length on each side
    tops_left = jordan_chain_tops_nilpotent(Nq, max_exp, p)   # Dict[L] -> (dimVq x mult_L)
    tops_right = jordan_chain_tops_nilpotent(Nqs, max_exp, p)

    if not tops_left and max_exp >= 1:
        # semisimple (or implementation returns {} for exp=1): tops at L=1 span the whole space
        tops_left = {1: np.eye(Fq.shape[0], dtype=np.int64)}
    if not tops_right and max_exp >= 1:
        tops_right = {1: np.eye(Fqs.shape[0], dtype=np.int64)}

    lengths = sorted(set(tops_left.keys()) | set(tops_right.keys()))

    # Ambient pairing between the two primary bases
    P = _pairing_matrix_between(Vq, Vqs, p)  # (dimVq x dimVqs)

    n2 = F.shape[0]
    Ω_amb = omega_matrix(n2 // 2, p)

    blocks: List[AtomicBlock] = []
    built_cols: List[np.ndarray] = []

    try:
        for L in lengths:
            A_raw = tops_left.get(L,  np.zeros((Fq.shape[0], 0), dtype=np.int64))
            B_raw = tops_right.get(L, np.zeros((Fqs.shape[0], 0), dtype=np.int64))

            A = _select_module_generators_from_top_space(Fq,  Nq,  A_raw, deg_q, int(L), p)
            B = _select_module_generators_from_top_space(Fqs, Nqs, B_raw, deg_q, int(L), p)

            a_mult = int(A.shape[1])
            b_mult = int(B.shape[1])
            inv_data["length_multiplicities"][int(L)] = (a_mult, b_mult)

            if a_mult != b_mult:
                raise RuntimeError(f"Paired sector mismatch at length {L}: {a_mult} vs {b_mult}")
            if a_mult == 0:
                continue

            # Pairing must be computed at the *chain level*:
            # M_L = (N^{L-1} A)^T P B
            Npow = _mat_pow_mod(Nq, int(L) - 1, p)
            NA = mod_p(Npow @ A, p)  # (dimVq x a_mult)

            M = mod_p(NA.T @ P @ B, p)  # (a_mult x a_mult)
            rM = rank_mod(M, p)
            inv_data["pairing_rank"][int(L)] = int(rM)
            if rM != a_mult:
                raise RuntimeError(
                    f"Paired sector: chain-level top pairing singular at length {L}. "
                    f"rank={rM}, expected={a_mult}."
                )

            Minv = inv_mod_mat(M, p)
            B_dual = mod_p(B @ Minv, p)

            # Now build one atomic block per paired top
            for j in range(a_mult):
                v_top = A[:, j:j + 1]        # coords in Vq basis
                w_top = B_dual[:, j:j + 1]   # coords in Vqs basis

                # Cyclic submodule bases in sector coordinates
                C_left = cyclic_submodule_basis(Fq, Nq, v_top, deg_q, int(L), p)    # (dimVq x deg*L)
                C_right = cyclic_submodule_basis(Fqs, Nqs, w_top, deg_q, int(L), p) # (dimVqs x deg*L)

                # Lift to ambient
                W_left = mod_p(Vq @ C_left, p)    # (2n x deg*L)
                W_right = mod_p(Vqs @ C_right, p)

                span = independent_columns(np.concatenate([W_left, W_right], axis=1), p)

                # Paired sector blocks should be nondegenerate
                if not is_nondegenerate(Ω_amb, span, p):
                    raise RuntimeError(f"Paired sector: constructed atomic span is degenerate (L={L}, j={j}).")

                T_blk = darboux_basis_from_span(Ω_amb, span, p)  # Darboux basis for the block
                blocks.append(
                    AtomicBlock(
                        T_blk=mod_p(T_blk, p),
                        half_dim=int(T_blk.shape[1] // 2),
                        sector_key=key,
                        inv=None,
                    )
                )
                built_cols.append(T_blk)
            
        # Global paired-sector span check: blocks should span W = Vq ⊕ Vqs, and lie inside it.
        W_sector = independent_columns(np.concatenate([Vq, Vqs], axis=1), p)
        dim_sector = rank_mod(W_sector, p)

        if built_cols:
            all_cols = np.concatenate(built_cols, axis=1)
            # (i) blocks lie in the sector: span(W_sector, all_cols) has same rank as W_sector
            if rank_mod(np.concatenate([W_sector, all_cols], axis=1), p) != dim_sector:
                raise RuntimeError("Paired sector: some constructed block columns lie outside W = Vq ⊕ Vq*.")
            # (ii) blocks span the sector
            dim_blocks = rank_mod(all_cols, p)
        else:
            dim_blocks = 0

        if dim_blocks != dim_sector:
            raise RuntimeError(
                f"Paired sector: constructed blocks do not span the paired sector subspace "
                f"(dim_blocks={dim_blocks}, dim_sector={dim_sector})."
            )

        # if rank_mod(W_blocks, p) != rank_mod(W_sector, p):
        #     raise RuntimeError("Paired sector: constructed blocks do not span the paired sector subspace.")

        inv_data["checks_passed"].append("per-length chain pairing nonsingular")
        inv_data["checks_passed"].append("each block nondegenerate")
        inv_data["checks_passed"].append("blocks span paired sector")
        inv_data["status"] = "OK"

    except Exception as e:
        # Non-certified fallback: return one big paired-sector block so callers don't crash.
        # (This keeps smoke tests happy but does NOT claim correctness/certification.)
        if not allow_fallback:
            raise Exception(inv_data)
        inv_data["status"] = "DEGRADED"
        inv_data["note"] = f"fallback: {type(e).__name__}: {e}"

        W_sector = independent_columns(np.concatenate([Vq, Vqs], axis=1), p)
        if not is_nondegenerate(Ω_amb, W_sector, p):
            # This *shouldn't* happen for a true paired sector. If it does, we should fail loudly.
            raise RuntimeError("Paired sector fallback: sector span is degenerate (unexpected).") from e
        T_blk = darboux_basis_from_span(Ω_amb, W_sector, p)
        blocks = [AtomicBlock(mod_p(T_blk, p), int(T_blk.shape[1] // 2), key, None)]

    inv = AtomicInvariant(sector_key=key, sector_type="paired", poly_key=key, data=inv_data)
    blocks = [AtomicBlock(b.T_blk, b.half_dim, b.sector_key, inv) for b in blocks]
    # Sanity: blocks must be a direct sum inside W = Vq ⊕ Vq*
    if blocks:
        all_cols = np.concatenate([b.T_blk for b in blocks], axis=1)
        if rank_mod(all_cols, p) != sum(b.T_blk.shape[1] for b in blocks):
            raise RuntimeError("Paired sector: produced blocks overlap (not direct sum).")

    return blocks, inv
