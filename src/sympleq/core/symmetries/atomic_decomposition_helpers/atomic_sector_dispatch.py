from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

from ..modular_helpers import mod_p, omega_matrix
from .atomic_linear import darboux_basis_from_span, is_nondegenerate, restrict_operator
from .atomic_paired import atomic_blocks_in_paired_sector
from .atomic_self import atomic_blocks_in_self_sector_nonunipotent
from .atomic_self_p2_unitary import atomic_blocks_in_self_sector_p2_nonunipotent_unitary
from .atomic_types import (
    AtomicBlock,
    AtomicInvariant,
    ExtractionObstruction,
    SearchBudgetExceeded,
    SectorContext,
    SectorCostCertificate,
)
from .atomic_unipotent_p2 import atomic_blocks_in_unipotent_self_sector_p2
from .module_invariants import q_of_F_restricted
from .rcf_prepass import _is_x_pm_1


def classify_extraction_error(e: BaseException) -> str:
    """Classify a sector extraction failure for CertificationError diagnostics."""
    if isinstance(e, SearchBudgetExceeded):
        return "budget"
    if isinstance(e, ExtractionObstruction):
        return "obstruction"
    return "other"


def sector_contexts_from_meta(meta: Dict[str, Any]) -> List[SectorContext]:
    """
    Return typed sector contexts from rcf_prepass.

    New prepass code returns ``sector_contexts`` directly. The legacy fallback
    keeps older cached/notebook prepass dictionaries usable, but production code
    should rely on the typed contexts.
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
                    "coordinate_note": "legacy context constructed in atomic_sector_dispatch.py",
                },
            )
        )
    return out


def ctx_meta(ctx: SectorContext) -> Dict[str, Any]:
    return ctx.meta if isinstance(ctx.meta, dict) else {}


def sector_debug(ctx: SectorContext, *, sector_index: int | None = None) -> Dict[str, Any]:
    meta = ctx_meta(ctx)
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


def sector_span_basis(ctx: SectorContext, p: int) -> np.ndarray:
    """Return an ambient basis for the full sector span."""
    meta = ctx_meta(ctx)
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


def sector_fallback_block(
    *,
    F: np.ndarray,
    p: int,
    ctx: SectorContext,
    reason: str,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Best-effort fallback: produce a single block spanning the whole sector subspace.

    Fallbacks are intentionally kept in the global best-effort wrapper, rather
    than hidden inside sector builders. Certified mode never calls this routine.
    """
    W = sector_span_basis(ctx, p)
    omega_amb = omega_matrix(F.shape[0] // 2, p)

    if W.size == 0 or W.shape[1] == 0:
        inv = AtomicInvariant(
            sector_key=ctx.sector_key,
            sector_type=ctx.sector_type,
            poly_key=ctx.poly_key,
            data={"status": "DEGRADED", "note": "empty sector basis", "reason": reason},
        )
        return [], inv

    if W.shape[1] % 2 != 0:
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

    if not is_nondegenerate(omega_amb, W, p):
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

    T_blk = darboux_basis_from_span(omega_amb, W, p)
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


def require_primaries(ctx: SectorContext) -> Dict[Tuple[int, ...], Dict[str, Any]]:
    primaries = ctx_meta(ctx).get("primaries")
    if not isinstance(primaries, dict):
        raise RuntimeError(f"Sector {ctx.sector_key} context is missing primary data.")
    return primaries


def build_sector_raw(
    F: np.ndarray,
    p: int,
    ctx: SectorContext,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Build one sector from a typed SectorContext.

    Sector builders are called with their internal fallbacks disabled. A sector
    that cannot be built deterministically raises; the decomposition wrapper
    decides whether to fail strictly or use a degraded fallback.
    """
    if int(ctx.p) != int(p):
        raise ValueError(f"Sector context has p={ctx.p}, but decomposition requested p={p}.")

    prim = require_primaries(ctx)

    if ctx.sector_type == "paired":
        if ctx.sector_key_star is None:
            raise RuntimeError(f"Paired sector {ctx.sector_key} is missing reciprocal partner key.")
        return atomic_blocks_in_paired_sector(
            F,
            p,
            ctx.sector_key,
            ctx.sector_key_star,
            prim,
        )

    sector_key = ctx.sector_key
    info = prim[sector_key]
    q = info["poly"]

    if p == 2 and _is_x_pm_1(q, p):
        return atomic_blocks_in_unipotent_self_sector_p2(F, sector_key, prim)

    T_sec = ctx.T_sec
    F_sec = ctx.F_sec
    Omega_sec = ctx.Omega_sec
    N_sec = ctx.N_sec

    if T_sec is None or F_sec is None or Omega_sec is None or N_sec is None:
        W = sector_span_basis(ctx, p)
        omega_amb = omega_matrix(F.shape[0] // 2, p)
        T_sec = darboux_basis_from_span(omega_amb, W, p)
        F_sec = restrict_operator(F, T_sec, p)
        Omega_sec = omega_matrix(T_sec.shape[1] // 2, p)
        N_sec = q_of_F_restricted(F_sec, q, p)

    if p == 2 and int(ctx.deg_q) > 1:
        return atomic_blocks_in_self_sector_p2_nonunipotent_unitary(
            F_sec=F_sec,
            T_sec=T_sec,
            Omega=Omega_sec,
            N=N_sec,
            deg_q=int(ctx.deg_q),
            max_exp=int(ctx.max_exp),
            p=p,
            sector_key=sector_key,
            poly_key=ctx.poly_key,
        )

    return atomic_blocks_in_self_sector_nonunipotent(
        F_sec=F_sec,
        T_sec=T_sec,
        Omega=Omega_sec,
        N=N_sec,
        deg_q=int(ctx.deg_q),
        max_exp=int(ctx.max_exp),
        p=p,
        sector_key=sector_key,
        poly_key=ctx.poly_key,
    )


def inject_invariant_lower_bound(
    inv: AtomicInvariant, ctx: SectorContext, blocks: List[AtomicBlock]
) -> None:
    """
    Attach the invariant-derived lower bound from the prepass to one sector.

    ``sector_cost`` is the verified constructed sector cost; the certificate is
    minimal only when that attained cost equals the independent lower bound.
    """
    if not isinstance(getattr(inv, "data", None), dict):
        return
    meta = ctx.meta or {}
    lb = meta.get("cost_lower_bound")
    bound_complete = bool(meta.get("cost_lower_bound_complete", False))
    sector_cost = max((int(b.half_dim) for b in blocks), default=0)
    status_ok = inv.data.get("status", "OK") == "OK"

    legacy_cc = dict(inv.data.get("cost_certificate") or {})
    legacy_cc["lengths_present"] = list(meta.get("lengths_present", []))
    note = legacy_cc.get("note") or ""
    note = (note + " | " if note else "") + "lower_bound from invariants"

    lower_bound_complete = bool(bound_complete and lb is not None and status_ok)
    sector_cert = SectorCostCertificate(
        sector_cost=int(sector_cost),
        lower_bound=None if lb is None else int(lb),
        lower_bound_complete=lower_bound_complete,
        extraction_attained=bool(status_ok),
        certified_minimal_sector=bool(lower_bound_complete and lb is not None and sector_cost == int(lb)),
        note=note,
        extra={
            k: v
            for k, v in legacy_cc.items()
            if k
            not in {
                "sector_cost",
                "lower_bound",
                "lower_bound_complete",
                "attained",
                "extraction_attained",
                "complete",
                "certified",
                "certified_minimal_sector",
                "note",
            }
        },
    )
    inv.data["sector_cost_certificate"] = sector_cert
    inv.data["cost_certificate"] = sector_cert.as_dict()


def build_sector(
    F: np.ndarray,
    p: int,
    ctx: SectorContext,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Dispatch to the sector builders, then attach the invariant-derived cost
    lower bound to the returned certificate.
    """
    blocks, inv = build_sector_raw(F, p, ctx)
    try:
        inject_invariant_lower_bound(inv, ctx, blocks)
    except Exception:
        # Bound injection must never break extraction; leave the builder's cert.
        pass
    return blocks, inv
