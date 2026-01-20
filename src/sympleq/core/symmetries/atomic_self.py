# sympleq/core/symmetries/atomic_self.py
from __future__ import annotations
import numpy as np
from typing import List, Tuple

from .modular_helpers import mod_p, independent_columns
from .atomic_types import AtomicBlock, AtomicInvariant
from .symplectic_basis import symplectic_gram_schmidt_from_span
from .module_invariants import restrict_operator_invariant, q_of_F_restricted, jordan_chain_tops_nilpotent, cyclic_submodule_basis


def atomic_blocks_in_self_sector_nonunipotent(
    F: np.ndarray, p: int, key: Tuple[int, ...], primaries: dict
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Placeholder for self-reciprocal sector q=q* with q != x±1.
    You will ultimately:
      - build N = q(F) on V_q,
      - get chain-top data in quotient spaces,
      - compute induced form invariants on those quotients,
      - assemble indecomposable symplectic summands,
      - construct bases deterministically.
    """
    q = primaries[key]["poly"]
    V = independent_columns(primaries[key]["V_basis"], p)
    Fv = restrict_operator_invariant(F, V, p)
    N = q_of_F_restricted(Fv, q, p)
    max_exp = int(primaries[key]["exponent"])

    tops = jordan_chain_tops_nilpotent(N, max_exp, p)

    # TODO: compute self-sector indecomposable types from (tops + induced forms).
    inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data={
        "status": "TODO",
        "note": "self-reciprocal non-unipotent sector needs induced-form classification"
    })

    # Fallback: return the whole sector as one block (NOT optimal/certified)
    # Replace once classification is implemented.
    T_blk = symplectic_gram_schmidt_from_span(V, p)
    blocks = [AtomicBlock(T_blk=T_blk, half_dim=T_blk.shape[1] // 2, sector_key=key, inv=inv)]
    return blocks, inv
