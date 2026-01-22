# sympleq/core/symmetries/atomic_types.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple, List, Literal

import numpy as np
from numpy.typing import NDArray

SectorType = Literal["paired", "self"]


@dataclass(frozen=True, slots=True)
class AtomicInvariant:
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
class AtomicBlock:
    """
    One atomic (indecomposable) nondegenerate F-invariant symplectic summand.

    T_blk columns form a hyperbolic basis in ambient coordinates:
        T_blk^T Ω T_blk = Ω_k
    """
    T_blk: NDArray[np.int64]             # (2n)×(2k)
    half_dim: int                        # k
    sector_key: Tuple[int, ...]          # which sector this block came from
    inv: Optional[AtomicInvariant] = None


@dataclass(frozen=True, slots=True)
class AtomicDecomposition:
    """
    Full decomposition result: F = T S T^{-1} with S block-diagonal in atomic symplectic blocks.
    """
    S: NDArray[np.int64]                 # (2n)×(2n)
    T: NDArray[np.int64]                 # (2n)×(2n), symplectic
    blocks: List[AtomicBlock]            # in the order they appear in S
    invariants: List[AtomicInvariant]    # per-sector certificates (optionally duplicated per block)
    qudit_cost: int                      # max_k over blocks
