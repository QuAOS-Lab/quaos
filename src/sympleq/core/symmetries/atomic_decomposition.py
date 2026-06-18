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

Convention = Literal["column", "row"]


class CertificationError(RuntimeError):
    """
    Raised by :func:`decompose_or_raise` when the decomposition is not a
    certified atomic (or certified minimal) result.
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
    info["cost_certificate"] = cert
    # Deprecated aliases of the two canonical fields above (``certified_minimal``
    # and ``qudit_cost``); retained as shims for older notebooks/callers.
    info["certified_minimal_qudit_cost"] = bool(cert["certified_minimal"])
    info["minimal_cost_certified"] = bool(cert["certified_minimal"])
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


def _build_sector_raw(
    F: np.ndarray,
    p: int,
    ctx: SectorContext,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Build one sector from a typed SectorContext.

    Sector builders are always called with their internal fallbacks disabled; a
    sector that cannot be built deterministically raises, and the single
    decomposition route absorbs it via a degraded fallback plus global
    completion (marking the result uncertified).
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
    )
    return blocks, inv


def _inject_invariant_lower_bound(
    inv: AtomicInvariant, ctx: SectorContext, blocks: List[AtomicBlock]
) -> None:
    """
    Phase 2: overwrite the sector cost certificate's ``lower_bound`` with the
    invariant-derived value computed in the prepass (``ctx.meta``), replacing the
    old circular ``lower_bound = attained``.  ``attained`` becomes the verified
    constructed sector cost; ``certified_minimal_sector`` records whether the two
    agree.  Done in this single dispatch wrapper so all four builders are covered.
    """
    if not isinstance(getattr(inv, "data", None), dict):
        return
    meta = ctx.meta or {}
    lb = meta.get("cost_lower_bound")
    bound_complete = bool(meta.get("cost_lower_bound_complete", False))
    sector_cost = max((int(b.half_dim) for b in blocks), default=0)
    status_ok = (inv.data.get("status", "OK") == "OK")

    cc = dict(inv.data.get("cost_certificate") or {})
    cc["sector_cost"] = int(sector_cost)
    cc["attained"] = True  # a decomposition was constructed; cost is sector_cost
    cc["lower_bound"] = None if lb is None else int(lb)
    cc["complete"] = bool(bound_complete and lb is not None and status_ok)
    cc["certified"] = cc["complete"]
    cc["certified_minimal_sector"] = bool(lb is not None and sector_cost == int(lb))
    cc["lengths_present"] = list(meta.get("lengths_present", []))
    note = cc.get("note") or ""
    cc["note"] = (note + " | " if note else "") + "lower_bound from invariants (Phase 2)"
    inv.data["cost_certificate"] = cc


def _build_sector(
    F: np.ndarray,
    p: int,
    ctx: SectorContext,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """Dispatch to the sector builders, then attach the invariant-derived
    (search-independent) cost lower bound to the returned certificate."""
    blocks, inv = _build_sector_raw(F, p, ctx)
    try:
        _inject_invariant_lower_bound(inv, ctx, blocks)
    except Exception:
        # Bound injection must never break extraction; leave the builder's cert.
        pass
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


def atomic_block_decompose(
    F: np.ndarray,
    p: int,
    *,
    convention: Convention = "column",
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Decompose a symplectic ``F`` into atomic invariant symplectic blocks.

    Single route: this always runs once and always returns ``(Sigma, B, info)``
    with ``B`` symplectic and ``Sigma = B^{-1} F B``.  Certification is reported
    in the output rather than selected by a mode:

      * ``info["certified"]`` -- True iff the returned basis is exactly the
        concatenated atomic block frame (no global completion was needed), every
        sector was built by an implemented family (no degraded fallback), and the
        independent verification passed.  A valid-but-not-atomic result (e.g. one
        relying on global completion) has this False.
      * ``info["cost_certificate"]["certified_minimal"]`` (mirrored at
        ``info["certified_minimal"]``) -- True iff the attained qudit cost equals
        the invariant-derived lower bound.  This minimality proof is computed from
        conjugacy invariants in the prepass, independent of the construction.

    Callers that want a hard failure on anything less than a certified minimal
    decomposition should use :func:`decompose_or_raise`.
    """
    if convention not in ("column", "row"):
        raise ValueError(f"Unknown convention={convention!r}. Expected 'column' or 'row'.")

    F = np.asarray(F, dtype=int)
    if F.ndim != 2 or F.shape[0] != F.shape[1]:
        raise ValueError(f"F must be a 2D square matrix, got shape {getattr(F, 'shape', None)}.")
    n2 = F.shape[0]
    if n2 % 2 != 0:
        raise ValueError(f"F must be (2n)x(2n), got shape {F.shape}.")
    F = mod_p(F, p)

    if convention == "row":
        Sigma_col, B_col, info = atomic_block_decompose(F.T, p, convention="column")
        return _convert_column_result_to_row(Sigma_col, B_col, info, p)

    if not is_symplectic(F, p):
        raise ValueError("Input F is not symplectic in column-action convention.")

    meta = rcf_prepass(F, p)
    sector_contexts = _sector_contexts_from_meta(meta)

    blocks: List[AtomicBlock] = []
    sector_invariants: List[AtomicInvariant] = []
    errors: list[dict[str, Any]] = []

    for i, ctx in enumerate(sector_contexts):
        try:
            b, inv = _build_sector(F, p, ctx)
            blocks += b
            sector_invariants.append(inv)
        except Exception as e:
            reason = f"{type(e).__name__}: {e}"
            dbg = _sector_debug(ctx, sector_index=i)
            dbg["error"] = reason
            dbg["error_kind"] = _classify_extraction_error(e)
            errors.append(dbg)
            # Degraded fallback: keep going and let global completion absorb the
            # span. This can never raise out of the single route.
            try:
                b_fb, inv_fb = _sector_fallback_block(F=F, p=p, ctx=ctx, reason=reason)
            except Exception as e2:
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

    # Assemble a full symplectic basis. Completion is a no-op when the atomic
    # blocks already span (the certified case); otherwise it absorbs the rest.
    B, completed = _complete_global_basis_from_blocks(F=F, p=p, blocks=blocks)
    Sigma = _compute_sigma(F, B, p)

    try:
        verification = verify_atomic_decomposition(
            F, B, Sigma, blocks, sector_invariants, p, expect_full_block_cover=False
        )
        verif_error = None
    except Exception as e:
        verification = {"error": f"{type(e).__name__}: {e}"}
        verif_error = verification["error"]

    all_sectors_ok = all(inv.data.get("status") == "OK" for inv in sector_invariants)
    block_cover = bool(verification.get("block_cover", False)) if isinstance(verification, dict) else False
    certified = bool(
        (not completed) and (not errors) and all_sectors_ok and block_cover and (verif_error is None)
    )

    warnings: List[str] = []
    if completed:
        warnings.append("Global symplectic completion was used; atomic blocks do not span by themselves.")
    if errors:
        warnings.append("At least one sector used a degraded fallback block.")
    if verif_error is not None:
        warnings.append("Independent verification reported an error; inspect verification['error'].")
    if not certified:
        warnings.append("Result is a valid decomposition but not a certified atomic one.")

    info: Dict[str, Any] = {
        "status": "OK" if certified else "DEGRADED",
        "certified": certified,
        "input_convention": "column",
        "internal_convention": "column",
        "sector_invariants": sector_invariants,
        "atomic_half_dims": [int(b.half_dim) for b in blocks],
        "Lmin_star": int(meta.get("Lmin_star", 1)),
        "completed": bool(completed),
        "warnings": warnings,
        "verification": verification,
    }
    if errors:
        info["errors"] = errors

    # Minimality lives entirely in the (invariant-derived) cost certificate.
    _attach_cost_certificate(info, blocks, sector_invariants, completed=bool(completed), p=p)
    return Sigma, B, info


def decompose_or_raise(
    F: np.ndarray,
    p: int,
    *,
    convention: Convention = "column",
    require_minimal: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Strict wrapper around :func:`atomic_block_decompose`.

    Raises :class:`CertificationError` unless the result is a certified atomic
    decomposition (``info["certified"]``) and, when ``require_minimal`` (default),
    also certified minimal (``info["certified_minimal"]``).  Library code should
    prefer :func:`atomic_block_decompose` and read the flags directly.
    """
    Sigma, B, info = atomic_block_decompose(F, p, convention=convention)
    if not info.get("certified", False):
        raise CertificationError(
            "Result is not a certified atomic decomposition.",
            info={"warnings": info.get("warnings"), "errors": info.get("errors")},
        )
    if require_minimal and not info.get("certified_minimal", False):
        raise CertificationError(
            "Result is a certified atomic decomposition but minimality is not certified.",
            info={"cost_certificate": info.get("cost_certificate")},
        )
    return Sigma, B, info
