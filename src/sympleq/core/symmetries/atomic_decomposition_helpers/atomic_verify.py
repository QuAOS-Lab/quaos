# sympleq/core/symmetries/atomic_decomposition_helpers/atomic_verify.py
from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from ..modular_helpers import mod_p, omega_matrix, inv_mod_mat, rank_mod
from .atomic_types import AtomicBlock, AtomicInvariant, SectorCostCertificate
from .atomic_linear import split_uv


def _as_int_list(xs: Iterable[int]) -> List[int]:
    return [int(x) for x in xs]


def _omega_for_dim(dim: int, p: int) -> np.ndarray:
    if dim % 2 != 0:
        raise ValueError(f"Symplectic dimension must be even, got {dim}.")
    return omega_matrix(dim // 2, p)


def _block_frame(blocks: Sequence[AtomicBlock], n2: int, p: int) -> np.ndarray:
    """
    Reconstruct the global Darboux frame used by atomic_decomposition:
    all U halves first, then all V halves.
    """
    if not blocks:
        return np.zeros((n2, 0), dtype=np.int64)
    U_cols: List[np.ndarray] = []
    V_cols: List[np.ndarray] = []
    for block in blocks:
        U, V = split_uv(block.T_blk)
        U_cols.append(U)
        V_cols.append(V)
    return mod_p(np.concatenate(U_cols + V_cols, axis=1), p)


def _block_index_sets(blocks: Sequence[AtomicBlock]) -> List[np.ndarray]:
    """
    Indices of each atomic block inside the global frame ordering
    [U_1, ..., U_r, V_1, ..., V_r].
    """
    half_dims = [int(b.half_dim) for b in blocks]
    total_k = int(sum(half_dims))
    out: List[np.ndarray] = []
    u0 = 0
    v0 = total_k
    for k in half_dims:
        inds = list(range(u0, u0 + k)) + list(range(v0, v0 + k))
        out.append(np.array(inds, dtype=np.int64))
        u0 += k
        v0 += k
    return out


def verify_block(F: np.ndarray, T_blk: np.ndarray, p: int, *, check_invariant: bool = True) -> Dict[str, Any]:
    """
    Verify one returned atomic block basis.

    Returns a diagnostics dictionary and raises RuntimeError on failure.
    """
    F = mod_p(np.asarray(F, dtype=np.int64), p)
    T_blk = mod_p(np.asarray(T_blk, dtype=np.int64), p)
    n2 = F.shape[0]
    if F.shape != (n2, n2):
        raise RuntimeError(f"F must be square, got {F.shape}.")
    if T_blk.ndim != 2 or T_blk.shape[0] != n2:
        raise RuntimeError(f"T_blk has incompatible shape {T_blk.shape}; expected first dimension {n2}.")
    if T_blk.shape[1] % 2 != 0:
        raise RuntimeError(f"T_blk has odd number of columns {T_blk.shape[1]}.")

    k2 = int(T_blk.shape[1])
    k = k2 // 2
    Omega = _omega_for_dim(n2, p)
    Omega_k = _omega_for_dim(k2, p)

    r = rank_mod(T_blk, p)
    if r != k2:
        raise RuntimeError(f"Block columns are dependent: rank={r}, columns={k2}.")

    gram = mod_p(T_blk.T @ Omega @ T_blk, p)
    if not np.array_equal(gram % p, Omega_k % p):
        raise RuntimeError("Block is not Darboux: T_blk^T Omega T_blk != Omega_k.")

    invariant = True
    if check_invariant:
        rank_aug = rank_mod(np.concatenate([T_blk, mod_p(F @ T_blk, p)], axis=1), p)
        invariant = bool(rank_aug == k2)
        if not invariant:
            raise RuntimeError(
                f"Block span is not F-invariant: rank([T,FT])={rank_aug}, rank(T)={k2}."
            )

    return {"half_dim": int(k), "dim": int(k2), "rank": int(r), "invariant": bool(invariant)}


def verify_block_family(F: np.ndarray, blocks: Sequence[AtomicBlock], p: int) -> Dict[str, Any]:
    """Verify all atomic blocks and their pairwise symplectic orthogonality."""
    F = mod_p(np.asarray(F, dtype=np.int64), p)
    n2 = F.shape[0]
    Omega = _omega_for_dim(n2, p)

    block_infos = []
    for i, block in enumerate(blocks):
        info = verify_block(F, block.T_blk, p)
        info["index"] = int(i)
        info["sector_key"] = tuple(block.sector_key)
        block_infos.append(info)

    for i in range(len(blocks)):
        Ti = mod_p(blocks[i].T_blk, p)
        for j in range(i + 1, len(blocks)):
            Tj = mod_p(blocks[j].T_blk, p)
            gij = mod_p(Ti.T @ Omega @ Tj, p)
            if np.any(gij % p):
                nz = np.argwhere(gij % p)
                raise RuntimeError(
                    "Atomic blocks are not symplectically orthogonal: "
                    f"i={i}, j={j}, nnz={int(nz.shape[0])}, first={_as_int_list(nz[0]) if nz.size else []}."
                )

    return {
        "n_blocks": int(len(blocks)),
        "total_dim": int(sum(2 * int(b.half_dim) for b in blocks)),
        "blocks": block_infos,
    }


def verify_global_basis(F: np.ndarray, B: np.ndarray, Sigma: np.ndarray, p: int) -> Dict[str, Any]:
    """Verify B is symplectic and Sigma = B^{-1} F B."""
    F = mod_p(np.asarray(F, dtype=np.int64), p)
    B = mod_p(np.asarray(B, dtype=np.int64), p)
    Sigma = mod_p(np.asarray(Sigma, dtype=np.int64), p)
    n2 = F.shape[0]
    if F.shape != (n2, n2) or B.shape != (n2, n2) or Sigma.shape != (n2, n2):
        raise RuntimeError(f"Shape mismatch: F={F.shape}, B={B.shape}, Sigma={Sigma.shape}.")

    Omega = _omega_for_dim(n2, p)
    gram = mod_p(B.T @ Omega @ B, p)
    if not np.array_equal(gram % p, Omega % p):
        raise RuntimeError("Global basis is not symplectic: B^T Omega B != Omega.")

    Sigma_expected = mod_p(inv_mod_mat(B, p) @ F @ B, p)
    if not np.array_equal(Sigma % p, Sigma_expected % p):
        raise RuntimeError("Sigma mismatch: Sigma != B^{-1} F B.")

    return {"n2": int(n2), "rank_B": int(rank_mod(B, p)), "sigma_verified": True}


def verify_sigma_block_diagonal(Sigma: np.ndarray, blocks: Sequence[AtomicBlock], p: int) -> Dict[str, Any]:
    """
    Verify Sigma has no off-block entries in the block index sets induced by the
    frame ordering [U_1,...,U_r,V_1,...,V_r].
    """
    Sigma = mod_p(np.asarray(Sigma, dtype=np.int64), p)
    n2 = Sigma.shape[0]
    if Sigma.shape != (n2, n2):
        raise RuntimeError(f"Sigma must be square, got {Sigma.shape}.")

    idx_sets = _block_index_sets(blocks)
    covered = np.concatenate(idx_sets) if idx_sets else np.zeros(0, dtype=np.int64)
    if len(covered) != n2:
        raise RuntimeError(
            f"Blocks do not cover Sigma dimensions: covered={len(covered)}, n2={n2}."
        )

    for bi, idx in enumerate(idx_sets):
        mask = np.ones(n2, dtype=bool)
        mask[idx] = False
        off_rows = Sigma[np.ix_(idx, np.where(mask)[0])]
        off_cols = Sigma[np.ix_(np.where(mask)[0], idx)]
        if np.any(off_rows % p) or np.any(off_cols % p):
            raise RuntimeError(f"Sigma is not block diagonal with respect to atomic block {bi}.")

    return {"block_diagonal": True, "block_dims": [int(2 * b.half_dim) for b in blocks]}


def verify_cost_certificate(
    blocks: Sequence[AtomicBlock],
    invariants: Sequence[AtomicInvariant],
    p: int,
    *,
    completed: bool = False,
) -> Dict[str, Any]:
    """
    Conservative global minimal-cost certificate aggregation.

    A global minimality claim is made only if every sector invariant carries an
    explicit complete cost certificate with a finite lower_bound and attained=True.
    """
    qudit_cost = max((int(b.half_dim) for b in blocks), default=0)
    sector_certificates: List[Dict[str, Any]] = []
    lower_bounds: List[int] = []
    missing: List[Dict[str, Any]] = []
    incomplete: List[Dict[str, Any]] = []

    for inv in invariants:
        data = inv.data if isinstance(inv.data, dict) else {}
        cert_obj = data.get("sector_cost_certificate")
        cert = cert_obj.as_dict() if isinstance(cert_obj, SectorCostCertificate) else data.get("cost_certificate")
        label = {
            "sector_key": tuple(inv.sector_key),
            "sector_type": inv.sector_type,
            "poly_key": tuple(inv.poly_key),
        }
        if data.get("status") != "OK":
            item = {**label, "reason": f"status={data.get('status')}"}
            incomplete.append(item)
            sector_certificates.append({**label, "complete": False, "attained": False, "lower_bound": None})
            continue
        if not isinstance(cert, dict):
            item = {**label, "reason": "missing sector cost_certificate"}
            missing.append(item)
            sector_certificates.append({**label, "complete": False, "attained": False, "lower_bound": None})
            continue

        lb = cert.get("lower_bound")
        attained = bool(cert.get("extraction_attained", cert.get("attained", False)))
        complete = bool(cert.get("complete", False))
        try:
            lb_int: Optional[int] = None if lb is None else int(lb)
        except Exception:
            lb_int = None

        item = {
            **label,
            "lower_bound": lb_int,
            "attained": attained,
            "complete": complete,
            "sector_cost": cert.get("sector_cost"),
            "note": cert.get("note", ""),
        }
        sector_certificates.append(item)

        if not complete or not attained or lb_int is None:
            incomplete.append({**label, "reason": "incomplete/ unattained/ missing lower_bound"})
        else:
            lower_bounds.append(lb_int)

    all_complete = (not completed) and (len(missing) == 0) and (len(incomplete) == 0)
    global_lower_bound = max(lower_bounds) if all_complete and lower_bounds else (0 if all_complete else None)
    certified_minimal = bool(all_complete and global_lower_bound == qudit_cost)

    return {
        "qudit_cost": int(qudit_cost),
        "lower_bound": None if global_lower_bound is None else int(global_lower_bound),
        # The attained cost is the verified max_i k_i over constructed blocks
        # (= qudit_cost), no longer hardcoded True.
        "attained": int(qudit_cost),
        "complete": bool(all_complete),
        "certified_minimal": bool(certified_minimal),
        "completed_global_basis": bool(completed),
        "sector_certificates": sector_certificates,
        "missing": missing,
        "incomplete": incomplete,
    }


def verify_atomic_decomposition(
    F: np.ndarray,
    B: np.ndarray,
    Sigma: np.ndarray,
    blocks: Sequence[AtomicBlock],
    invariants: Sequence[AtomicInvariant],
    p: int,
    *,
    expect_full_block_cover: bool = True,
) -> Dict[str, Any]:
    """End-to-end verification oracle for returned atomic decompositions."""
    global_info = verify_global_basis(F, B, Sigma, p)
    family_info = verify_block_family(F, blocks, p)

    n2 = np.asarray(F).shape[0]
    block_frame = _block_frame(blocks, n2, p)
    block_cover = bool(block_frame.shape == B.shape and np.array_equal(block_frame % p, mod_p(B, p) % p))
    if expect_full_block_cover and not block_cover:
        raise RuntimeError(
            "Returned B is not exactly the concatenated atomic block frame; "
            "this is expected only for best_effort global completion."
        )

    sigma_info: Dict[str, Any] = {"block_diagonal_checked": False}
    if block_cover:
        sigma_info = verify_sigma_block_diagonal(Sigma, blocks, p)
        sigma_info["block_diagonal_checked"] = True

    cost_info = verify_cost_certificate(blocks, invariants, p, completed=not block_cover)
    return {
        "global": global_info,
        "blocks": family_info,
        "sigma": sigma_info,
        "cost_certificate": cost_info,
        "block_cover": bool(block_cover),
    }
