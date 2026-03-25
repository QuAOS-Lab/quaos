# sympleq/core/symmetries/atomic_decomposition.py
from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Tuple, Literal, Optional

from .atomic_decomposition_helpers.rcf_prepass import rcf_prepass, _is_x_pm_1
from .modular_helpers import (
    mod_p,
    omega_matrix,
    inv_mod_mat,
    is_symplectic,
    rank_mod,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_types import AtomicBlock, AtomicInvariant
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_linear import (
    split_uv,
    is_nondegenerate,
    darboux_basis_from_span,
    symplectic_completion_from_block,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_paired import atomic_blocks_in_paired_sector
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_self import atomic_blocks_in_self_sector_nonunipotent
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2 import (
    atomic_blocks_in_unipotent_self_sector_p2)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_linear import restrict_operator
from sympleq.core.symmetries.atomic_decomposition_helpers.module_invariants import q_of_F_restricted

Mode = Literal["auto", "certified", "best_effort"]


class CertificationError(RuntimeError):
    """
    Raised by atomic_block_decompose_certified when a certified decomposition
    cannot be produced.
    """

    def __init__(self, message: str, *, info: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.info: Dict[str, Any] = info or {}


def _concat_blocks_to_partial_basis(blocks: List[AtomicBlock], n2: int, p: int) -> np.ndarray:
    """
    Build a *partial* symplectic frame T = [U_all | V_all] (2n x 2k) from block bases.
    Does NOT require spanning the full space.
    """
    if not blocks:
        return np.zeros((n2, 0), dtype=np.int64)

    U_list: list[np.ndarray] = []
    V_list: list[np.ndarray] = []
    for b in blocks:
        U_blk, V_blk = split_uv(b.T_blk)
        U_list.append(U_blk)
        V_list.append(V_blk)

    T = np.concatenate(U_list + V_list, axis=1)
    return mod_p(T, p)


def _verify_full_symplectic_basis(B: np.ndarray, p: int) -> None:
    """Raise if B is not a full symplectic basis."""
    if B.ndim != 2 or B.shape[0] != B.shape[1]:
        raise RuntimeError(f"Global basis must be square, got shape {B.shape}.")
    n2 = B.shape[0]
    if n2 % 2 != 0:
        raise RuntimeError(f"Global basis must have even dimension, got {n2}.")
    n = n2 // 2
    Omega = omega_matrix(n, p)
    G = mod_p(B.T @ Omega @ B, p)
    if not np.array_equal(G % p, Omega % p):
        raise RuntimeError("Global basis B is not symplectic (B^T Omega B != Omega).")


def _compute_sigma(F: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """Sigma = B^{-1} F B (mod p)."""
    return mod_p(inv_mod_mat(B, p) @ F @ B, p)


def _sector_fallback_block(
    *,
    F: np.ndarray,
    p: int,
    sec: Dict[str, Any],
    reason: str,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Best-effort fallback: produce a single block spanning the whole sector subspace.

    This is NOT a claim of atomic optimality; it is only meant to keep the global
    pipeline returning a valid decomposition if the sector subspace is nondegenerate.
    """
    W = mod_p(np.asarray(sec["W_basis"], dtype=np.int64), p)
    n2 = F.shape[0]
    Omega_amb = omega_matrix(n2 // 2, p)

    if W.size == 0 or W.shape[1] == 0:
        key = tuple(sec.get("key", ()))
        inv = AtomicInvariant(
            sector_key=key,
            sector_type="self" if sec.get("type") == "self" else "paired",
            poly_key=key,
            data={"status": "DEGRADED", "note": "empty sector basis", "reason": reason},
        )
        return [], inv

    if W.shape[1] % 2 != 0:
        raise RuntimeError(
            f"Best-effort fallback: sector has odd dimension {W.shape[1]} (cannot be symplectic). "
            f"type={sec.get('type')} key={sec.get('key')}"
        )

    if not is_nondegenerate(Omega_amb, W, p):
        raise RuntimeError(
            f"Best-effort fallback: sector span is degenerate (cannot build Darboux basis). "
            f"type={sec.get('type')} key={sec.get('key')}"
        )

    T_blk = darboux_basis_from_span(Omega_amb, W, p)

    if sec.get("type") == "paired":
        key = tuple(sec["key"])
        sector_type: Literal["paired", "self"] = "paired"
    else:
        key = tuple(sec["key"])
        sector_type = "self"

    inv = AtomicInvariant(
        sector_key=key,
        sector_type=sector_type,
        poly_key=key,
        data={
            "status": "DEGRADED",
            "note": "sector fallback block (not certified)",
            "reason": reason,
        },
    )

    block = AtomicBlock(
        T_blk=mod_p(T_blk, p),
        half_dim=int(T_blk.shape[1] // 2),
        sector_key=key,
        inv=inv,
    )
    return [block], inv


def _build_sector(
    F: np.ndarray,
    p: int,
    sec: Dict[str, Any],
    prim: Dict[Tuple[int, ...], Dict[str, Any]],
    *,
    certified: bool,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:

    if sec["type"] == "paired":
        key = sec["key"]
        key_star = sec["key_star"]
        blocks, inv = atomic_blocks_in_paired_sector(F, p, key, key_star, prim, allow_fallback=not certified)
        return blocks, inv

    # self sector
    sector_key = sec["key"]
    info = prim[sector_key]

    q = info["poly"]
    if p == 2 and _is_x_pm_1(q, p):
        blocks, inv = atomic_blocks_in_unipotent_self_sector_p2(F, sector_key, prim)
        return blocks, inv

    # inside _build_sector, self sector branch, after p=2 unipotent special case:

    W = sec["W_basis"]                        # ambient basis spanning sector
    Omega_amb = omega_matrix(F.shape[0] // 2, p)  # ambient standard Omega

    # build a symplectic basis for the sector in ambient coordinates
    T_sec = darboux_basis_from_span(Omega_amb, W, p)

    # restrict F to sector coords
    F_sec = restrict_operator(F, T_sec, p)

    # N = q(F_sec) where q is the primary polynomial
    q = prim[sector_key]["poly"]
    N_sec = q_of_F_restricted(F_sec, q, p)

    deg_q = int(sec["deg"])
    max_exp = int(sec["exponent"])
    Omega_sec = omega_matrix(T_sec.shape[1] // 2, p)  # since T_sec is symplectic

    blocks, inv = atomic_blocks_in_self_sector_nonunipotent(
        F_sec=F_sec,
        T_sec=T_sec,
        Omega=Omega_sec,
        N=N_sec,
        deg_q=deg_q,
        max_exp=max_exp,
        p=p,
        sector_key=sector_key,
        poly_key=tuple(sector_key),
        allow_fallback=not certified,   # IMPORTANT
    )
    return blocks, inv


def _complete_global_basis_from_blocks(
    *,
    F: np.ndarray,
    p: int,
    blocks: List[AtomicBlock],
) -> Tuple[np.ndarray, bool]:
    """
    Return a full symplectic basis B, and a boolean `completed` indicating whether we
    had to perform a completion step beyond concatenating block bases.
    """
    n2 = F.shape[0]
    n = n2 // 2
    Omega = omega_matrix(n, p)

    T = _concat_blocks_to_partial_basis(blocks, n2, p)  # 2n × 2k
    if T.size == 0:
        # No blocks at all: fall back to identity (valid symplectic basis).
        B = np.eye(n2, dtype=np.int64)
        _verify_full_symplectic_basis(B, p)
        return mod_p(B, p), True

    # Ensure the frame is actually symplectic on its span:
    # i.e. T^T Omega T == Omega_k. If not, that's a bug upstream (block builder)
    k2 = T.shape[1]
    if k2 % 2 != 0:
        raise RuntimeError(f"Partial basis has odd number of columns {k2}, cannot be a symplectic frame.")
    k = k2 // 2
    Omega_k = omega_matrix(k, p)
    G = mod_p(T.T @ Omega @ T, p)
    if not np.array_equal(G % p, Omega_k % p):
        raise RuntimeError("Partial basis T is not a Darboux frame: T^T Omega T != Omega_k. Upstream block bug?")

    if k2 == n2:
        B = T
        _verify_full_symplectic_basis(B, p)
        return mod_p(B, p), False

    # Complete the symplectic frame to a full symplectic basis.
    B = symplectic_completion_from_block(T, p)
    _verify_full_symplectic_basis(B, p)
    return mod_p(B, p), True


def _first_nonorthogonal_pair(blocks, p):
    """Return (i, j, Gij) where T_i^T Omega T_j != 0 for any i < j."""
    if not blocks:
        return None
    n2 = blocks[0].T_blk.shape[0]
    Omega = omega_matrix(n2 // 2, p)
    for i in range(len(blocks)):
        Ti = blocks[i].T_blk
        for j in range(i + 1, len(blocks)):
            Tj = blocks[j].T_blk
            Gij = mod_p(Ti.T @ Omega @ Tj, p)
            if np.any(Gij % p):
                return (i, j, Gij)
    return None


def atomic_block_decompose_certified(F: np.ndarray, p: int) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Certified decomposition (strict):
      - Every sector must return inv.data["status"] == "OK"
      - Block bases must span the full space (no completion allowed)
      - B must be symplectic
    """
    if not is_symplectic(F, p):
        raise ValueError("Input F is not symplectic.")

    F = mod_p(F, p)
    n2 = F.shape[0]
    if n2 % 2 != 0:
        raise ValueError(f"F must be (2n)x(2n), got shape {F.shape}.")

    meta = rcf_prepass(F, p)
    prim = meta["primaries"]
    sectors = meta["sectors"]

    blocks: List[AtomicBlock] = []
    sector_invariants: List[AtomicInvariant] = []
    failures: list[dict[str, Any]] = []

    for i, sec in enumerate(sectors):
        try:
            b, inv = _build_sector(F, p, sec, prim, certified=True)
        except Exception as e:
            failures.append(
                {
                    "sector_index": i,
                    "sector_type": sec.get("type"),
                    "key": sec.get("key"),
                    "key_star": sec.get("key_star"),
                    "error": f"{type(e).__name__}: {e}",
                    "deg": sec.get("deg"),
                    "exponent": sec.get("exponent"),
                    "dim2": sec.get("dim2"),
                    "sec_note": sec.get("note", ""),
                }
            )
            continue

        blocks += b
        sector_invariants.append(inv)

        if inv.data.get("status") != "OK":
            failures.append(
                {
                    "sector_index": i,
                    "sector_type": sec.get("type"),
                    "key": sec.get("key"),
                    "key_star": sec.get("key_star"),
                    "status": inv.data.get("status"),
                    "note": inv.data.get("note", ""),
                    "deg": sec.get("deg"),
                    "exponent": sec.get("exponent"),
                    "dim2": sec.get("dim2"),
                    "sec_note": sec.get("note", ""),
                    "builder_debug": inv.data.get("last_attempts", inv.data.get("progress_lengths", None)),
                }
            )

    if failures:
        raise CertificationError(
            "Certified decomposition failed: at least one sector is uncertified or errored.",
            info={
                "p": int(p),
                "n2": int(n2),
                "failures": failures,
                "Lmin_star": int(meta.get("Lmin_star", 1)),
                "sector_signature": [(sec.get("type"), sec.get("key"), sec.get("key_star")) for sec in sectors],
            },
        )

    # In certified mode we REQUIRE spanning; no completion allowed.
    T = _concat_blocks_to_partial_basis(blocks, n2, p)
    if T.shape != (n2, n2):
        raise CertificationError(
            f"Certified decomposition failed: global basis has wrong shape {T.shape} (blocks do not span).",
            info={
                "p": int(p),
                "n2": int(n2),
                "atomic_half_dims": [b.half_dim for b in blocks],
                "rank_partial": int(rank_mod(T, p)) if T.size else 0,
            },
        )

    pair = _first_nonorthogonal_pair(blocks, p)
    if pair is not None:
        i, j, Gij = pair
        nz = np.argwhere(Gij % p)
        nz_preview = [tuple(map(int, x)) for x in nz[:10]]  # first 10 positions

        raise CertificationError(
            "Certified decomposition failed: found non-orthogonal pair of atomic blocks.",
            info={
                "p": int(p),
                "n2": int(n2),
                "atomic_half_dims": [int(b.half_dim) for b in blocks],
                "pair": {
                    "i": int(i),
                    "j": int(j),
                    "half_dims": (int(blocks[i].half_dim), int(blocks[j].half_dim)),
                    "sectors": (
                        (blocks[i].sector_key, getattr(blocks[i], "sector_type", None)),
                        (blocks[j].sector_key, getattr(blocks[j], "sector_type", None)),
                    ),
                    "Gij_nnz": int(nz.shape[0]),
                    "Gij_nz_preview": nz_preview,
                },
            },
        )

    B = T
    _verify_full_symplectic_basis(B, p)
    Sigma = _compute_sigma(F, B, p)

    info = {
        "status": "OK",
        "certified": True,
        "sector_invariants": sector_invariants,
        "atomic_half_dims": [b.half_dim for b in blocks],
        "Q_opt": max([b.half_dim for b in blocks], default=0),
        "Lmin_star": int(meta.get("Lmin_star", 1)),
        "completed": False,
    }
    return Sigma, B, info


def atomic_block_decompose_best_effort(
    F: np.ndarray, p: int, *, last_error: Optional[CertificationError] = None
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Best-effort decomposition:
      - Always attempts to return (Sigma, B) with B symplectic and Sigma = B^{-1} F B.
      - Allows sector fallbacks and ALSO global completion if blocks don't span.
    """
    if not is_symplectic(F, p):
        raise ValueError("Input F is not symplectic.")

    F = mod_p(F, p)
    n2 = F.shape[0]
    if n2 % 2 != 0:
        raise ValueError(f"F must be (2n)x(2n), got shape {F.shape}.")

    meta = rcf_prepass(F, p)
    prim = meta["primaries"]
    sectors = meta["sectors"]

    blocks: List[AtomicBlock] = []
    sector_invariants: List[AtomicInvariant] = []
    errors: list[dict[str, Any]] = []

    for i, sec in enumerate(sectors):
        try:
            b, inv = _build_sector(F, p, sec, prim, certified=False)
            blocks += b
            sector_invariants.append(inv)
        except Exception as e:
            reason = f"{type(e).__name__}: {e}"
            errors.append(
                {
                    "sector_index": i,
                    "sector_type": sec.get("type"),
                    "key": sec.get("key"),
                    "key_star": sec.get("key_star"),
                    "error": reason,
                }
            )
            b_fb, inv_fb = _sector_fallback_block(F=F, p=p, sec=sec, reason=reason)
            blocks += b_fb
            sector_invariants.append(inv_fb)

    # Global completion step (THIS is what fixes your failure)
    B, completed = _complete_global_basis_from_blocks(F=F, p=p, blocks=blocks)
    Sigma = _compute_sigma(F, B, p)

    certified = (not completed) and all(inv.data.get("status") == "OK" for inv in sector_invariants)
    status = "OK" if certified else "DEGRADED"

    info: Dict[str, Any] = {
        "status": status,
        "certified": bool(certified),
        "sector_invariants": sector_invariants,
        "atomic_half_dims": [b.half_dim for b in blocks],
        "Q_opt": max([b.half_dim for b in blocks], default=0),
        "Lmin_star": int(meta.get("Lmin_star", 1)),
        "completed": bool(completed),
    }
    if errors:
        info["errors"] = errors
    if last_error is not None:
        info["last_certification_error"] = getattr(last_error, "info", {}) or {"message": str(last_error)}

    return Sigma, B, info


def atomic_block_decompose(F: np.ndarray, p: int, mode: Mode = "auto") -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Public entry point.

    mode:
      - "certified": return only if fully certified, else raise CertificationError
      - "best_effort": always return a decomposition (may be degraded)
      - "auto": try certified first, else fall back to best-effort (default)
    """
    if mode == "certified":
        return atomic_block_decompose_certified(F, p)
    if mode == "best_effort":
        return atomic_block_decompose_best_effort(F, p)
    if mode != "auto":
        raise ValueError(f"Unknown mode={mode!r}. Expected 'auto','certified','best_effort'.")

    try:
        return atomic_block_decompose_certified(F, p)
    except CertificationError as e:
        return atomic_block_decompose_best_effort(F, p, last_error=e)
