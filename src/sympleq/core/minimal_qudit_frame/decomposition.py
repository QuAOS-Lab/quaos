from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional, Tuple

import numpy as np

from .helpers.mqf_certification import attach_cost_certificate
from .helpers.mqf_completion import (
    complete_global_basis_from_blocks,
    compute_sigma,
    convert_column_result_to_row,
)
from .helpers.mqf_sector_dispatch import (
    build_sector,
    classify_extraction_error,
    sector_contexts_from_meta,
    sector_debug,
    sector_fallback_block,
)
from .helpers.mqf_types import MQFBlock, MQFInvariant
from .helpers.mqf_verify import verify_minimal_qudit_frame
from .helpers.rcf_prepass import rcf_prepass
from sympleq.core.symmetries.modular_helpers import is_symplectic, mod_p

Convention = Literal["column", "row"]


class MinimalQuditFrameCertificationError(RuntimeError):
    """
    Raised by :func:`minimal_qudit_frame_or_raise` when the decomposition is not a
    certified MQF (or certified minimal) result.
    """

    def __init__(self, message: str, *, info: Optional[Dict[str, Any]] = None):
        super().__init__(message)
        self.info: Dict[str, Any] = info or {}


def minimal_qudit_frame(
    F: np.ndarray,
    p: int,
    *,
    convention: Convention = "column",
    allow_degraded: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Decompose a symplectic ``F`` into MQF invariant symplectic blocks.

    Strict by default: sector extraction failures raise ``MinimalQuditFrameCertificationError``
    instead of being hidden behind a global completion. Set
    ``allow_degraded=True`` to recover the previous best-effort behaviour, which
    always tries to return ``(Sigma, B, info)`` with ``B`` symplectic and reports
    certification in the output flags.

      * ``info["certified"]`` -- True iff the returned basis is exactly the
        concatenated MQF block frame (no global completion was needed), every
        sector was built by an implemented family (no degraded fallback), and the
        independent verification passed. A valid-but-not-MQF result (e.g. one
        relying on global completion) has this False.
      * ``info["cost_certificate"]["certified_minimal"]`` (mirrored at
        ``info["certified_minimal"]``) -- True iff the attained qudit cost equals
        the invariant-derived lower bound. This minimality proof is computed from
        conjugacy invariants in the prepass, independent of the construction.

    Callers that want a hard failure on anything less than a certified minimal
    decomposition should use :func:`minimal_qudit_frame_or_raise`.
    """
    if convention not in ("column", "row"):
        raise ValueError(f"Unknown convention={convention!r}. Expected 'column' or 'row'.")

    F = np.asarray(F, dtype=int)
    if F.ndim != 2 or F.shape[0] != F.shape[1]:
        raise ValueError(f"F must be a 2D square matrix, got shape {getattr(F, 'shape', None)}.")
    if F.shape[0] % 2 != 0:
        raise ValueError(f"F must be (2n)x(2n), got shape {F.shape}.")
    F = mod_p(F, p)

    if convention == "row":
        Sigma_col, B_col, info = minimal_qudit_frame(
            F.T, p, convention="column", allow_degraded=allow_degraded
        )
        return convert_column_result_to_row(Sigma_col, B_col, info, p)

    if not is_symplectic(F, p):
        raise ValueError("Input F is not symplectic in column-action convention.")

    meta = rcf_prepass(F, p)
    sector_contexts = sector_contexts_from_meta(meta)

    blocks: List[MQFBlock] = []
    sector_invariants: List[MQFInvariant] = []
    errors: list[dict[str, Any]] = []

    for i, ctx in enumerate(sector_contexts):
        try:
            sector_blocks, inv = build_sector(F, p, ctx)
            blocks += sector_blocks
            sector_invariants.append(inv)
        except Exception as e:
            reason = f"{type(e).__name__}: {e}"
            dbg = sector_debug(ctx, sector_index=i)
            dbg["error"] = reason
            dbg["error_kind"] = classify_extraction_error(e)
            errors.append(dbg)
            if not allow_degraded:
                raise MinimalQuditFrameCertificationError(
                    "Sector extraction failed in strict MQF decomposition mode.",
                    info={"errors": errors, "failures": errors, "failed_sector": dbg},
                ) from e

            try:
                fallback_blocks, inv_fb = sector_fallback_block(F=F, p=p, ctx=ctx, reason=reason)
            except Exception as e2:
                dbg["fallback_error"] = f"{type(e2).__name__}: {e2}"
                fallback_blocks, inv_fb = [], MQFInvariant(
                    sector_key=ctx.sector_key,
                    sector_type=ctx.sector_type,
                    poly_key=ctx.poly_key,
                    data={
                        "status": "DEGRADED",
                        "note": "sector fallback failed; deferred to global completion",
                        "reason": reason,
                    },
                )
            blocks += fallback_blocks
            sector_invariants.append(inv_fb)

    B, completed = complete_global_basis_from_blocks(F=F, p=p, blocks=blocks)
    Sigma = compute_sigma(F, B, p)

    try:
        verification = verify_minimal_qudit_frame(
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
        warnings.append("Global symplectic completion was used; MQF blocks do not span by themselves.")
    if errors:
        warnings.append("At least one sector used a degraded fallback block.")
    if verif_error is not None:
        warnings.append("Independent verification reported an error; inspect verification['error'].")
    if not certified:
        warnings.append("Result is a symplectic basis reduction but not a certified MQF block decomposition.")

    info: Dict[str, Any] = {
        "status": "OK" if certified else "DEGRADED",
        "certified": certified,
        "input_convention": "column",
        "internal_convention": "column",
        "sector_invariants": sector_invariants,
        "mqf_half_dims": [int(b.half_dim) for b in blocks],
        "Lmin_star": int(meta.get("Lmin_star", 1)),
        "completed": bool(completed),
        "allow_degraded": bool(allow_degraded),
        "mode": "best_effort" if allow_degraded else "certified",
        "warnings": warnings,
        "verification": verification,
    }
    if errors:
        info["errors"] = errors

    attach_cost_certificate(info, blocks, sector_invariants, completed=bool(completed), p=p)
    return Sigma, B, info


def minimal_qudit_frame_or_raise(
    F: np.ndarray,
    p: int,
    *,
    convention: Convention = "column",
    require_minimal: bool = True,
) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Strict wrapper around :func:`minimal_qudit_frame`.

    Raises :class:`MinimalQuditFrameCertificationError` unless the result is a certified MQF
    decomposition (``info["certified"]``) and, when ``require_minimal`` (default),
    also certified minimal (``info["certified_minimal"]``). Library code should
    prefer :func:`minimal_qudit_frame` and read the flags directly.
    """
    Sigma, B, info = minimal_qudit_frame(F, p, convention=convention, allow_degraded=False)
    if not info.get("certified", False):
        raise MinimalQuditFrameCertificationError(
            "Result is not a certified MQF decomposition.",
            info={"warnings": info.get("warnings"), "errors": info.get("errors")},
        )
    if require_minimal and not info.get("certified_minimal", False):
        raise MinimalQuditFrameCertificationError(
            "Result is a certified MQF decomposition but minimality is not certified.",
            info={"cost_certificate": info.get("cost_certificate")},
        )
    return Sigma, B, info
