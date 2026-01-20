# sympleq/core/symmetries/atomic_unipotent_p2.py
from __future__ import annotations
import numpy as np
from typing import List, Tuple

from .modular_helpers import mod_p, independent_columns
from .atomic_types import AtomicBlock, AtomicInvariant
from .symplectic_basis import symplectic_gram_schmidt_from_span
from .module_invariants import restrict_operator_invariant

def classify_unipotent_sp2(
    F_u: np.ndarray,  # restriction of F to the unipotent sector basis coords
) -> dict:
    """
    Compute the conjugacy invariants for unipotent elements in Sp(2m,2).

    MUST RETURN enough data to determine the indecomposable summands uniquely.
    Typical outputs include:
      - Jordan partition of N=F+I (over GF(2)),
      - additional 'chi'/Arf-type bits distinguishing W vs V types,
      - multiplicities of W(m), V(2k), and alpha-twisted blocks.

    This is the missing certified piece.
    """
    raise NotImplementedError("Implement Sp(2m,2) unipotent classification invariants here.")

def build_unipotent_blocks_from_invariants(
    F: np.ndarray,
    V_u: np.ndarray,   # ambient basis columns for the unipotent sector
    inv_data: dict,
    p: int = 2
) -> List[AtomicBlock]:
    """
    Construct explicit atomic block bases in the ambient coordinates, matching inv_data exactly.
    Should return a list of AtomicBlock with hyperbolic bases.
    """
    raise NotImplementedError("Construct unipotent atomic blocks from invariants.")

def atomic_blocks_in_unipotent_self_sector_p2(
    F: np.ndarray, key: Tuple[int, ...], primaries: dict
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Entry point for p=2, q=x+1 (or x-1; same in GF(2)) self sector.
    """
    p = 2
    V_u = independent_columns(primaries[key]["V_basis"], p)
    F_u = restrict_operator_invariant(F, V_u, p)

    inv_data = classify_unipotent_sp2(F_u)
    inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data=inv_data)

    blocks = build_unipotent_blocks_from_invariants(F, V_u, inv_data, p=2)
    blocks = [AtomicBlock(b.T_blk, b.half_dim, key, inv) for b in blocks]
    return blocks, inv
