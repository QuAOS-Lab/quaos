from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple, List, Literal, TypedDict

import numpy as np
from numpy.typing import NDArray

SectorType = Literal["paired", "self"]


class CostCertificate(TypedDict, total=False):
    """Canonical (typed) shape of the global cost certificate.

    ``lower_bound`` is an invariant-derived quantity, independent of the
    constructed decomposition, so ``certified_minimal`` (``lower_bound ==
    attained``) is a genuine theorem rather than true-by-construction.
    """
    qudit_cost: int            # attained max half-dim over blocks
    lower_bound: Optional[int]  # max over invariant sector bounds
    attained: int              # verified attained cost (== qudit_cost)
    complete: bool             # every sector contributed a certified bound
    certified_minimal: bool    # lower_bound == qudit_cost and complete
    sector_certificates: List[Dict[str, Any]]
    missing: List[Dict[str, Any]]
    incomplete: List[Dict[str, Any]]


@dataclass(frozen=True, slots=True)
class SectorCostCertificate:
    """Typed sector-local cost certificate.

    ``lower_bound_complete`` records that the invariant-derived lower bound is
    available for this sector.  ``extraction_attained`` records that the sector
    extractor constructed verified blocks.  ``certified_minimal_sector`` is the
    stronger statement that the attained sector cost equals the invariant lower
    bound.

    The public ``cost_certificate`` stored on :class:`MQFInvariant` remains a
    plain dictionary for backward compatibility; use :meth:`as_dict` when
    serialising.
    """

    sector_cost: int
    lower_bound: Optional[int]
    lower_bound_complete: bool
    extraction_attained: bool
    certified_minimal_sector: bool
    note: str = ""
    extra: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> Dict[str, Any]:
        out: Dict[str, Any] = dict(self.extra)
        out.update(
            {
                "sector_cost": int(self.sector_cost),
                "lower_bound": None if self.lower_bound is None else int(self.lower_bound),
                "lower_bound_complete": bool(self.lower_bound_complete),
                "attained": bool(self.extraction_attained),
                "extraction_attained": bool(self.extraction_attained),
                "complete": bool(self.lower_bound_complete and self.extraction_attained),
                "certified_minimal_sector": bool(self.certified_minimal_sector),
                "certified": bool(self.certified_minimal_sector),
                "note": str(self.note),
            }
        )
        return out


# ---------------------------------------------------------------------------
# Extraction exception taxonomy
# ---------------------------------------------------------------------------
# Both subclass RuntimeError so existing ``except RuntimeError`` (and broader)
# handlers keep working unchanged; the distinct types let callers (e.g.
# MinimalQuditFrameCertificationError.info) tell a *search budget* failure apart from a genuine
# *mathematical obstruction*.

class ExtractionObstruction(RuntimeError):
    """A sector cannot be decomposed as required by the theory.

    Signals a genuine mathematical obstruction (e.g. no valid block can be
    extracted from the remaining invariant subspace), as opposed to merely
    running out of search budget.
    """


class SearchBudgetExceeded(RuntimeError):
    """A bounded/deterministic search gave up before exhausting possibilities.

    Signals that an iteration guard or candidate budget was hit, not that the
    decomposition is provably impossible.
    """


@dataclass(frozen=True, slots=True)
class MQFInvariant:
    """
    Deterministic sector-local conjugacy invariant summary for certification.

    - For paired sectors (q != q*): typically a multiset/partition-like summary of cyclic data.
    - For self-reciprocal sectors: includes the extra form labels needed for symplectic classification.
      In p=2 unipotent sectors this is where you store the W/V/alpha-type data.
    """
    sector_key: Tuple[int, ...]          # identifies the symplectic sector W_q (or pair)
    sector_type: SectorType              # "paired" or "self"
    poly_key: Tuple[int, ...]            # identifies q(x) itself (monic coeff tuple)
    data: Dict[str, Any]                 # structured, but kept free-form here


@dataclass(frozen=True, slots=True)
class MQFBlock:
    """
    One MQF (indecomposable) nondegenerate F-invariant symplectic summand.

    T_blk columns form a hyperbolic basis in ambient coordinates:
        T_blk^T Ω T_blk = Ω_k
    """
    T_blk: NDArray[np.int64]             # (2n)×(2k)
    half_dim: int                        # k
    sector_key: Tuple[int, ...]          # which sector this block came from
    inv: Optional[MQFInvariant] = None


@dataclass(frozen=True, slots=True)
class MinimalQuditFrameDecomposition:
    """
    Full decomposition result: F = T S T^{-1} with S block-diagonal in MQF symplectic blocks.
    """
    S: NDArray[np.int64]                 # (2n)×(2n)
    T: NDArray[np.int64]                 # (2n)×(2n), symplectic
    blocks: List[MQFBlock]            # in the order they appear in S
    invariants: List[MQFInvariant]    # per-sector certificates (optionally duplicated per block)
    qudit_cost: int                      # max_k over blocks


@dataclass(frozen=True, slots=True)
class SectorContext:
    # identity
    sector_key: Tuple[int, ...]          # q(x) coeff tuple for this sector (monic)
    sector_type: SectorType              # "paired" or "self"
    poly_key: Tuple[int, ...]            # same as sector_key for self; for paired can still be q

    # pairing partner (only for paired)
    sector_key_star: Optional[Tuple[int, ...]] = None

    # parameters
    p: int = 2
    deg_q: int = 0
    max_exp: int = 0                    # k(q)

    # sector coordinate embedding
    T_sec: np.ndarray | None = None            # (2n x dim_sec), columns basis of sector in ambient
    # optional inverse map ambient->sector coords if you want speed/stability
    T_sec_leftinv: Optional[np.ndarray] = None  # (dim_sec x 2n), s.t. leftinv @ T_sec = I

    # restricted operators in sector coords
    F_sec: np.ndarray | None = None            # (dim_sec x dim_sec)
    Omega_sec: np.ndarray | None = None        # (dim_sec x dim_sec) symplectic form restricted
    N_sec: Optional[np.ndarray] = None  # (dim_sec x dim_sec) nilpotent N = q(F)|_{V_q}

    # extra cached invariants (optional)
    meta: Dict[str, Any] | None = None


@dataclass(frozen=True, slots=True)
class PrepassContext:
    p: int
    n: int
    Omega: np.ndarray                  # ambient
    sectors: list[SectorContext]
    meta: Dict[str, Any]               # e.g. Lmin_star, signatures, etc
