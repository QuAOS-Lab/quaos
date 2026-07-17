from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

from .mqf_types import MQFBlock, MQFInvariant, SectorCostCertificate


def _as_int_or_none(value: Any) -> Optional[int]:
    try:
        return None if value is None else int(value)
    except Exception:
        return None


def _sector_label(inv: MQFInvariant) -> Dict[str, Any]:
    return {
        "sector_key": tuple(inv.sector_key),
        "sector_type": inv.sector_type,
        "poly_key": tuple(inv.poly_key),
    }


def _sector_certificate(inv: MQFInvariant) -> Dict[str, Any] | None:
    data = inv.data if isinstance(getattr(inv, "data", None), dict) else {}
    cert_obj = data.get("sector_cost_certificate")
    if isinstance(cert_obj, SectorCostCertificate):
        return cert_obj.as_dict()
    cert = data.get("cost_certificate")
    return cert if isinstance(cert, dict) else None


def verify_cost_certificate(
    blocks: Sequence[MQFBlock],
    invariants: Sequence[MQFInvariant],
    p: int,
    *,
    completed: bool = False,
) -> Dict[str, Any]:
    """
    Conservative global minimal-cost certificate aggregation.

    A global minimality claim is made only if every sector invariant carries an
    explicit complete cost certificate with a finite lower_bound and attained=True.
    The ``p`` argument is retained for API compatibility with the verification
    helpers; the aggregation itself is purely certificate bookkeeping.
    """
    qudit_cost = max((int(b.half_dim) for b in blocks), default=0)
    sector_certificates: List[Dict[str, Any]] = []
    lower_bounds: List[int] = []
    missing: List[Dict[str, Any]] = []
    incomplete: List[Dict[str, Any]] = []

    for inv in invariants:
        data = inv.data if isinstance(inv.data, dict) else {}
        cert = _sector_certificate(inv)
        label = _sector_label(inv)

        if data.get("status") != "OK":
            item = {**label, "reason": f"status={data.get('status')}"}
            incomplete.append(item)
            sector_certificates.append({**label, "complete": False, "attained": False, "lower_bound": None})
            continue
        if cert is None:
            item = {**label, "reason": "missing sector cost_certificate"}
            missing.append(item)
            sector_certificates.append({**label, "complete": False, "attained": False, "lower_bound": None})
            continue

        lb_int = _as_int_or_none(cert.get("lower_bound"))
        attained = bool(cert.get("extraction_attained", cert.get("attained", False)))
        complete = bool(cert.get("complete", False))

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
        "attained": int(qudit_cost),
        "complete": bool(all_complete),
        "certified_minimal": bool(certified_minimal),
        "completed_global_basis": bool(completed),
        "sector_certificates": sector_certificates,
        "missing": missing,
        "incomplete": incomplete,
    }


def build_sector_gap_report(invariants: List[MQFInvariant]) -> Dict[str, Any]:
    """
    Summarise sector-local certification gaps.

    A gap means the sector builder returned a verified block decomposition whose
    attained sector cost is larger than the invariant lower bound attached by
    the prepass.
    """
    sectors: List[Dict[str, Any]] = []
    gaps: List[Dict[str, Any]] = []
    incomplete: List[Dict[str, Any]] = []
    inconsistent: List[Dict[str, Any]] = []

    for idx, inv in enumerate(invariants):
        data = inv.data if isinstance(getattr(inv, "data", None), dict) else {}
        cert = _sector_certificate(inv) or {}
        lb_int = _as_int_or_none(cert.get("lower_bound"))
        sc_int = _as_int_or_none(cert.get("sector_cost", cert.get("qudit_cost")))
        complete = bool(cert.get("complete", cert.get("certified", False)))
        attained = bool(cert.get("extraction_attained", cert.get("attained", data.get("status") == "OK")))
        entry: Dict[str, Any] = {
            "sector_index": int(idx),
            **_sector_label(inv),
            "deg": data.get("deg"),
            "exponent": data.get("exponent"),
            "status": data.get("status"),
            "sector_cost": sc_int,
            "lower_bound": lb_int,
            "complete": complete,
            "attained": attained,
            "note": cert.get("note", cert.get("reason", data.get("note", ""))),
        }
        for key in (
            "length_summary",
            "length_multiplicities",
            "pairing_rank",
            "blocks",
            "p2_unipotent",
            "p2_self_reciprocal_nonunipotent",
            "minimality_guard",
        ):
            if key in data:
                entry[key] = data[key]
        sectors.append(entry)

        if not complete or not attained or lb_int is None or sc_int is None:
            incomplete.append(entry)
        elif sc_int > lb_int:
            gaps.append(entry)
        elif sc_int < lb_int:
            inconsistent.append(entry)

    return {
        "n_sectors": int(len(sectors)),
        "n_gaps": int(len(gaps)),
        "n_incomplete": int(len(incomplete)),
        "n_inconsistent": int(len(inconsistent)),
        "gaps": gaps,
        "incomplete": incomplete,
        "inconsistent": inconsistent,
        "sectors": sectors,
    }


def attach_cost_certificate(
    info: Dict[str, Any],
    blocks: List[MQFBlock],
    invariants: List[MQFInvariant],
    *,
    completed: bool,
    p: int,
) -> None:
    """
    Attach conservative global cost/minimality fields in-place.

    ``qudit_cost`` / ``Q_att`` is always the attained cost of the returned
    block frame. ``Q_opt`` is populated only when the invariant lower bound
    certifies that this attained cost is optimal.
    """
    cert = verify_cost_certificate(blocks, invariants, p, completed=completed)
    q_att = int(cert["qudit_cost"])
    is_minimal = bool(cert["certified_minimal"])

    info["qudit_cost"] = q_att
    info["Q_att"] = q_att
    info["attained_qudit_cost"] = q_att
    info["certified_lower_bound"] = cert["lower_bound"]
    info["certified_minimal"] = is_minimal
    info["cost_certificate"] = cert

    info["Q_opt"] = q_att if is_minimal else None
    info["optimal_qudit_cost"] = q_att if is_minimal else None

    info["certified_minimal_qudit_cost"] = is_minimal
    info["minimal_cost_certified"] = is_minimal

    gap_report = build_sector_gap_report(invariants)
    info["sector_gap_report"] = gap_report
    if gap_report["n_gaps"]:
        info.setdefault("warnings", []).append(
            "Minimality certificate gap: at least one sector attained a larger cost than its invariant lower bound."
        )
    if gap_report["n_inconsistent"]:
        info.setdefault("warnings", []).append(
            "Cost certificate inconsistency: at least one sector attained less than its invariant lower bound."
        )
