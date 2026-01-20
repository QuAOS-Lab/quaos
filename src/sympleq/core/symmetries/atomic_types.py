# sympleq/core/symmetries/atomic_types.py
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from typing import Any, Dict, Optional, Tuple


@dataclass(frozen=True)
class AtomicInvariant:
    """
    Stores a sector-local conjugacy invariant summary for certification.
    For paired sectors this can just be the multiset of cyclic lengths.
    For p=2 unipotent it must include the extra labels (W/V/alpha types).
    """
    sector_key: Tuple[int, ...]
    sector_type: str                  # "paired" or "self"
    poly_key: Tuple[int, ...]
    data: Dict[str, Any]              # free-form, but deterministic


@dataclass
class AtomicBlock:
    """
    One atomic symplectic indecomposable summand.
    T_blk columns are a hyperbolic basis in the ambient coordinates:
      T_blk^T Ω T_blk = Ω_k
    """
    T_blk: np.ndarray                 # (2n) x (2k)
    half_dim: int                     # k
    sector_key: Tuple[int, ...]
    inv: Optional[AtomicInvariant] = None
