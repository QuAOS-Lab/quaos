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
