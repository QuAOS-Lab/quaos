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
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_types import (
    AtomicBlock,
    AtomicInvariant,
    SectorContext,
    ExtractionObstruction,
    SearchBudgetExceeded,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_linear import (
    split_uv,
    is_nondegenerate,
    darboux_basis_from_span,
    symplectic_completion_from_block,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_paired import atomic_blocks_in_paired_sector
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_self import atomic_blocks_in_self_sector_nonunipotent
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_self_p2_unitary import (
    atomic_blocks_in_self_sector_p2_nonunipotent_unitary,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2 import (
    atomic_blocks_in_unipotent_self_sector_p2)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_linear import restrict_operator
from sympleq.core.symmetries.atomic_decomposition_helpers.module_invariants import q_of_F_restricted
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import (
    verify_atomic_decomposition,
    verify_global_basis,
    verify_cost_certificate,
)

Mode = Literal["auto", "certified", "best_effort"]
Convention = Literal["column", "row"]


class CertificationError(RuntimeError):
    """
    Raised by atomic_block_decompose_certified when a certified decomposition
    cannot be produced.
    """

    def __init__(self, message: str, *, info: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.info: Dict[str, Any] = info or {}


def _classify_extraction_error(e: BaseException) -> str:
    """Tag a sector failure as a search-budget limit, a mathematical
    obstruction, or something else, for CertificationError.info diagnostics."""
    if isinstance(e, SearchBudgetExceeded):
        return "budget"
    if isinstance(e, ExtractionObstruction):
        return "obstruction"
    return "other"


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


def _convert_column_result_to_row(Sigma_col: np.ndarray, B_col: np.ndarray, info: Dict[str, Any], p: int) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Convert an internally column-action result back to row/right-action convention."""
    Sigma_row = mod_p(Sigma_col.T, p)
    B_row = mod_p(B_col.T, p)
    info = dict(info)
    info["input_convention"] = "row"
    info["internal_convention"] = "column"
    info["basis_convention"] = "row"
    info["convention_note"] = (
        "Internal computation used F_col = F_row.T. Returned B has row-basis convention, "
        "so Sigma = B F_row B^{-1} = Sigma_col.T."
    )
    return Sigma_row, B_row, info


def _attach_cost_certificate(info: Dict[str, Any], blocks: List[AtomicBlock], invariants: List[AtomicInvariant], *, completed: bool, p: int) -> None:
    """Attach conservative global cost/minimality fields in-place."""
    cert = verify_cost_certificate(blocks, invariants, p, completed=completed)
    info["qudit_cost"] = int(cert["qudit_cost"])
    info["certified_lower_bound"] = cert["lower_bound"]
    info["certified_minimal"] = bool(cert["certified_minimal"])
    info["certified_minimal_qudit_cost"] = bool(cert["certified_minimal"])
    info["minimal_cost_certified"] = bool(cert["certified_minimal"])
    info["cost_certificate"] = cert
    # Backwards-compatible alias used by older notebooks.
    info["Q_opt"] = int(cert["qudit_cost"])


def _sector_contexts_from_meta(meta: Dict[str, Any]) -> List[SectorContext]:
    """
    Return typed sector contexts from rcf_prepass.

    New prepass code returns ``sector_contexts`` directly.  The small legacy
    fallback below keeps older cached/notebook prepass dictionaries usable, but
    production code should rely on the typed contexts.
    """
    ctxs = meta.get("sector_contexts")
    if ctxs is None:
        prepass_context = meta.get("prepass_context")
        ctxs = getattr(prepass_context, "sectors", None)
    if ctxs is not None:
        return list(ctxs)

    legacy_sectors = meta.get("sectors")
    primaries = meta.get("primaries")
    if legacy_sectors is None or primaries is None:
        raise RuntimeError("rcf_prepass did not return sector contexts or legacy sector data.")

    out: List[SectorContext] = []
    for index, sec in enumerate(legacy_sectors):
        key = tuple(sec["key"])
        sector_type = "paired" if sec.get("type") == "paired" else "self"
        out.append(
            SectorContext(
                sector_key=key,
                sector_type=sector_type,
                poly_key=key,
                sector_key_star=tuple(sec["key_star"]) if sec.get("key_star") is not None else None,
                p=int(meta["p"]),
                deg_q=int(sec.get("deg", primaries.get(key, {}).get("deg", 0))),
                max_exp=int(sec.get("exponent", primaries.get(key, {}).get("exponent", 0))),
                T_sec=sec.get("W_basis"),
                meta={
                    "index": int(index),
                    "legacy_sec": sec,
                    "W_basis": sec.get("W_basis"),
                    "primaries": primaries,
                    "Lmin_star": meta.get("Lmin_star", 1),
                    "coordinate_note": "legacy context constructed in atomic_decomposition.py",
                },
            )
        )
    return out


def _ctx_meta(ctx: SectorContext) -> Dict[str, Any]:
    return ctx.meta if isinstance(ctx.meta, dict) else {}


def _sector_debug(ctx: SectorContext, *, sector_index: int | None = None) -> Dict[str, Any]:
    meta = _ctx_meta(ctx)
    legacy = meta.get("legacy_sec", {}) if isinstance(meta.get("legacy_sec", {}), dict) else {}
    return {
        "sector_index": int(sector_index if sector_index is not None else meta.get("index", -1)),
        "sector_type": ctx.sector_type,
        "key": ctx.sector_key,
        "key_star": ctx.sector_key_star,
        "deg": int(ctx.deg_q),
        "exponent": int(ctx.max_exp),
        "dim2": int(ctx.T_sec.shape[1]) if isinstance(ctx.T_sec, np.ndarray) else legacy.get("dim2"),
        "sec_note": legacy.get("note", ""),
        "coordinate_note": meta.get("coordinate_note", ""),
    }


def _sector_span_basis(ctx: SectorContext, p: int) -> np.ndarray:
    """Return an ambient basis for the full sector span."""
    meta = _ctx_meta(ctx)
    W = meta.get("W_basis")
    if W is None:
        legacy = meta.get("legacy_sec", {})
        if isinstance(legacy, dict):
            W = legacy.get("W_basis")
    if W is None:
        W = ctx.T_sec
    if W is None:
        raise RuntimeError(f"Sector {ctx.sector_key} has no ambient sector basis in context.")
    return mod_p(np.asarray(W, dtype=np.int64), p)


def _sector_fallback_block(
    *,
    F: np.ndarray,
    p: int,
    ctx: SectorContext,
    reason: str,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Best-effort fallback: produce a single block spanning the whole sector subspace.

    Fallbacks are intentionally kept in this global best-effort wrapper, rather
    than hidden inside sector builders.  Certified mode never calls this routine.
    """
    W = _sector_span_basis(ctx, p)
    n2 = F.shape[0]
    Omega_amb = omega_matrix(n2 // 2, p)

    if W.size == 0 or W.shape[1] == 0:
        inv = AtomicInvariant(
            sector_key=ctx.sector_key,
            sector_type=ctx.sector_type,
            poly_key=ctx.poly_key,
            data={"status": "DEGRADED", "note": "empty sector basis", "reason": reason},
        )
        return [], inv

    if W.shape[1] % 2 != 0:
        # Best-effort must always return. An odd-dimensional sector span cannot
        # carry a symplectic form, so we cannot build a Darboux block here; emit
        # an empty DEGRADED block list and let global completion absorb the span.
        inv = AtomicInvariant(
            sector_key=ctx.sector_key,
            sector_type=ctx.sector_type,
            poly_key=ctx.poly_key,
            data={
                "status": "DEGRADED",
                "note": f"sector has odd dimension {W.shape[1]}; deferred to global completion",
                "reason": reason,
            },
        )
        return [], inv

    if not is_nondegenerate(Omega_amb, W, p):
        # Likewise, a degenerate sector span has no Darboux basis; defer to the
        # global symplectic completion rather than raising out of the loop.
        inv = AtomicInvariant(
            sector_key=ctx.sector_key,
            sector_type=ctx.sector_type,
            poly_key=ctx.poly_key,
            data={
                "status": "DEGRADED",
                "note": "sector span is degenerate; deferred to global completion",
                "reason": reason,
            },
        )
        return [], inv

    T_blk = darboux_basis_from_span(Omega_amb, W, p)

    inv = AtomicInvariant(
        sector_key=ctx.sector_key,
        sector_type=ctx.sector_type,
        poly_key=ctx.poly_key,
        data={
            "status": "DEGRADED",
            "note": "sector fallback block (not certified)",
            "reason": reason,
        },
    )

    block = AtomicBlock(
        T_blk=mod_p(T_blk, p),
        half_dim=int(T_blk.shape[1] // 2),
        sector_key=ctx.sector_key,
        inv=inv,
    )
    return [block], inv


def _require_primaries(ctx: SectorContext) -> Dict[Tuple[int, ...], Dict[str, Any]]:
    primaries = _ctx_meta(ctx).get("primaries")
    if not isinstance(primaries, dict):
        raise RuntimeError(f"Sector {ctx.sector_key} context is missing primary data.")
    return primaries


def _build_sector(
    F: np.ndarray,
    p: int,
    ctx: SectorContext,
    *,
    certified: bool,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Build one sector from a typed SectorContext.

    The ``certified`` flag is retained for diagnostics/API symmetry, but sector
    builders are always called with their internal fallbacks disabled.  Best-effort
    fallback is performed only by ``atomic_block_decompose_best_effort`` after an
    exception escapes this function.
    """
    if int(ctx.p) != int(p):
        raise ValueError(f"Sector context has p={ctx.p}, but decomposition requested p={p}.")

    prim = _require_primaries(ctx)

    if ctx.sector_type == "paired":
        if ctx.sector_key_star is None:
            raise RuntimeError(f"Paired sector {ctx.sector_key} is missing reciprocal partner key.")
        blocks, inv = atomic_blocks_in_paired_sector(
            F,
            p,
            ctx.sector_key,
            ctx.sector_key_star,
            prim,
            allow_fallback=False,
        )
        return blocks, inv

    # self sector
    sector_key = ctx.sector_key
    info = prim[sector_key]
    q = info["poly"]

    if p == 2 and _is_x_pm_1(q, p):
        blocks, inv = atomic_blocks_in_unipotent_self_sector_p2(F, sector_key, prim)
        return blocks, inv

    # Prefer precomputed sector coordinates from the context.  Recompute only if
    # the context was produced by an older prepass or coordinate construction failed.
    T_sec = ctx.T_sec
    F_sec = ctx.F_sec
    Omega_sec = ctx.Omega_sec
    N_sec = ctx.N_sec

    if T_sec is None or F_sec is None or Omega_sec is None or N_sec is None:
        W = _sector_span_basis(ctx, p)
        Omega_amb = omega_matrix(F.shape[0] // 2, p)
        T_sec = darboux_basis_from_span(Omega_amb, W, p)
        F_sec = restrict_operator(F, T_sec, p)
        Omega_sec = omega_matrix(T_sec.shape[1] // 2, p)
        N_sec = q_of_F_restricted(F_sec, q, p)

    if p == 2 and int(ctx.deg_q) > 1:
        blocks, inv = atomic_blocks_in_self_sector_p2_nonunipotent_unitary(
            F_sec=F_sec,
            T_sec=T_sec,
            Omega=Omega_sec,
            N=N_sec,
            deg_q=int(ctx.deg_q),
            max_exp=int(ctx.max_exp),
            p=p,
            sector_key=sector_key,
            poly_key=ctx.poly_key,
            allow_fallback=False,
        )
        return blocks, inv

    blocks, inv = atomic_blocks_in_self_sector_nonunipotent(
        F_sec=F_sec,
        T_sec=T_sec,
        Omega=Omega_sec,
        N=N_sec,
        deg_q=int(ctx.deg_q),
        max_exp=int(ctx.max_exp),
        p=p,
        sector_key=sector_key,
        poly_key=ctx.poly_key,
        allow_fallback=False,
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


def atomic_block_decompose_certified(
    F: np.ndarray,
    p: int,
    *,
    convention: Convention = "column",
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Certified decomposition (strict):
      - Every sector must return inv.data["status"] == "OK"
      - Block bases must span the full space (no completion allowed)
      - B must be symplectic
    """
    if convention not in ("column", "row"):
        raise ValueError(f"Unknown convention={convention!r}. Expected 'column' or 'row'.")

    F = np.asarray(F, dtype=int)

    if F.ndim != 2:
        raise ValueError(f"F must be a 2D square matrix, got ndim={F.ndim}.")

    if F.shape[0] != F.shape[1]:
        raise ValueError(f"F must be square, got shape {F.shape}.")

    n2 = F.shape[0]
    if n2 % 2 != 0:
        raise ValueError(f"F must be (2n)x(2n), got shape {F.shape}.")

    F = mod_p(F, p)

    if convention == "row":
        Sigma_col, B_col, info = atomic_block_decompose_certified(
            F.T,
            p,
            convention="column",
        )
        return _convert_column_result_to_row(Sigma_col, B_col, info, p)

    if not is_symplectic(F, p):
        raise ValueError("Input F is not symplectic in column-action convention.")

    # (F was already validated and reduced mod p above, before the convention
    # branch; the previously duplicated shape/mod checks here have been removed.)
    meta = rcf_prepass(F, p)
    sector_contexts = _sector_contexts_from_meta(meta)

    blocks: List[AtomicBlock] = []
    sector_invariants: List[AtomicInvariant] = []
    failures: list[dict[str, Any]] = []

    for i, ctx in enumerate(sector_contexts):
        try:
            b, inv = _build_sector(F, p, ctx, certified=True)
        except Exception as e:
            dbg = _sector_debug(ctx, sector_index=i)
            dbg["error"] = f"{type(e).__name__}: {e}"
            dbg["error_kind"] = _classify_extraction_error(e)
            failures.append(dbg)
            continue

        blocks += b
        sector_invariants.append(inv)

        if inv.data.get("status") != "OK":
            dbg = _sector_debug(ctx, sector_index=i)
            dbg.update(
                {
                    "status": inv.data.get("status"),
                    "note": inv.data.get("note", ""),
                    "builder_debug": inv.data.get("last_attempts", inv.data.get("progress_lengths", None)),
                }
            )
            failures.append(dbg)

    if failures:
        raise CertificationError(
            "Certified decomposition failed: at least one sector is uncertified or errored.",
            info={
                "p": int(p),
                "n2": int(n2),
                "failures": failures,
                "Lmin_star": int(meta.get("Lmin_star", 1)),
                "sector_signature": [
                    (ctx.sector_type, ctx.sector_key, ctx.sector_key_star) for ctx in sector_contexts
                ],
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

    try:
        verification = verify_atomic_decomposition(
            F, B, Sigma, blocks, sector_invariants, p, expect_full_block_cover=True
        )
    except Exception as e:
        raise CertificationError(
            "Certified decomposition failed final verification.",
            info={"p": int(p), "n2": int(n2), "error": f"{type(e).__name__}: {e}"},
        ) from e

    info = {
        "status": "OK",
        "certified": True,
        "input_convention": "column",
        "internal_convention": "column",
        "sector_invariants": sector_invariants,
        "atomic_half_dims": [int(b.half_dim) for b in blocks],
        "Lmin_star": int(meta.get("Lmin_star", 1)),
        "completed": False,
        "verification": verification,
    }
    _attach_cost_certificate(info, blocks, sector_invariants, completed=False, p=p)
    return Sigma, B, info


def atomic_block_decompose_best_effort(
    F: np.ndarray, p: int, *, last_error: Optional[CertificationError] = None, convention: Convention = "column"
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Best-effort decomposition:
      - Always attempts to return (Sigma, B) with B symplectic and Sigma = B^{-1} F B.
      - Allows sector fallbacks and ALSO global completion if blocks don't span.
    """
    if convention == "row":
        Sigma_col, B_col, info = atomic_block_decompose_best_effort(
            mod_p(F, p).T, p, last_error=last_error, convention="column"
        )
        return _convert_column_result_to_row(Sigma_col, B_col, info, p)
    if convention != "column":
        raise ValueError(f"Unknown convention={convention!r}. Expected 'column' or 'row'.")

    if not is_symplectic(F, p):
        raise ValueError("Input F is not symplectic in column-action convention.")

    F = mod_p(F, p)
    n2 = F.shape[0]
    if n2 % 2 != 0:
        raise ValueError(f"F must be (2n)x(2n), got shape {F.shape}.")

    meta = rcf_prepass(F, p)
    sector_contexts = _sector_contexts_from_meta(meta)

    blocks: List[AtomicBlock] = []
    sector_invariants: List[AtomicInvariant] = []
    errors: list[dict[str, Any]] = []

    for i, ctx in enumerate(sector_contexts):
        try:
            b, inv = _build_sector(F, p, ctx, certified=False)
            blocks += b
            sector_invariants.append(inv)
        except Exception as e:
            reason = f"{type(e).__name__}: {e}"
            dbg = _sector_debug(ctx, sector_index=i)
            dbg["error"] = reason
            errors.append(dbg)
            try:
                b_fb, inv_fb = _sector_fallback_block(F=F, p=p, ctx=ctx, reason=reason)
            except Exception as e2:
                # Best-effort must never raise; record and defer to global completion.
                dbg["fallback_error"] = f"{type(e2).__name__}: {e2}"
                b_fb, inv_fb = [], AtomicInvariant(
                    sector_key=ctx.sector_key,
                    sector_type=ctx.sector_type,
                    poly_key=ctx.poly_key,
                    data={
                        "status": "DEGRADED",
                        "note": "sector fallback failed; deferred to global completion",
                        "reason": reason,
                    },
                )
            blocks += b_fb
            sector_invariants.append(inv_fb)

    # Global completion step (THIS is what fixes your failure)
    B, completed = _complete_global_basis_from_blocks(F=F, p=p, blocks=blocks)
    Sigma = _compute_sigma(F, B, p)

    certified = (not completed) and all(inv.data.get("status") == "OK" for inv in sector_invariants)
    status = "OK" if certified else "DEGRADED"

    warnings: List[str] = []
    if completed:
        warnings.append("Global symplectic completion was used; atomic blocks do not span by themselves.")
    if errors:
        warnings.append("At least one sector used a best-effort fallback block.")
    if not certified:
        warnings.append("This result is a valid best-effort decomposition, not a certified atomic/minimal result.")

    try:
        verification = {
            "global": verify_global_basis(F, B, Sigma, p),
            "cost_certificate": verify_cost_certificate(blocks, sector_invariants, p, completed=completed),
        }
    except Exception as e:
        verification = {"error": f"{type(e).__name__}: {e}"}
        warnings.append("Final best-effort verification reported an error; inspect verification['error'].")

    info: Dict[str, Any] = {
        "status": status,
        "certified": bool(certified),
        "input_convention": "column",
        "internal_convention": "column",
        "sector_invariants": sector_invariants,
        "atomic_half_dims": [int(b.half_dim) for b in blocks],
        "Lmin_star": int(meta.get("Lmin_star", 1)),
        "completed": bool(completed),
        "warnings": warnings,
        "verification": verification,
    }
    _attach_cost_certificate(info, blocks, sector_invariants, completed=bool(completed), p=p)
    if errors:
        info["errors"] = errors
    if last_error is not None:
        info["last_certification_error"] = getattr(last_error, "info", {}) or {"message": str(last_error)}

    # Best-effort should never be confused with a certified minimal-cost proof.
    if status != "OK" or completed or errors:
        info["certified_minimal"] = False
        info["certified_minimal_qudit_cost"] = False
        info["minimal_cost_certified"] = False
        if isinstance(info.get("cost_certificate"), dict):
            info["cost_certificate"]["certified_minimal"] = False
            info["cost_certificate"]["complete"] = False

    return Sigma, B, info


def atomic_block_decompose(F: np.ndarray, p: int, mode: Mode = "auto", *, convention: Convention = "column") -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Public entry point.

    mode:
      - "certified": return only if fully certified, else raise CertificationError
      - "best_effort": always return a decomposition (may be degraded)
      - "auto": try certified first, else fall back to best-effort (default)

    convention:
      - "column": internal convention, x -> F x, returning Sigma = B^{-1} F B
      - "row": project convention, x -> x F; internally uses F.T and returns
        Sigma = B F B^{-1}
    """
    if mode == "certified":
        return atomic_block_decompose_certified(F, p, convention=convention)
    if mode == "best_effort":
        return atomic_block_decompose_best_effort(F, p, convention=convention)
    if mode != "auto":
        raise ValueError(f"Unknown mode={mode!r}. Expected 'auto','certified','best_effort'.")

    try:
        return atomic_block_decompose_certified(F, p, convention=convention)
    except CertificationError as e:
        return atomic_block_decompose_best_effort(F, p, last_error=e, convention=convention)
