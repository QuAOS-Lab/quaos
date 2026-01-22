# sympleq/core/symmetries/atomic_decomposition.py
from __future__ import annotations
import numpy as np
from typing import Dict, List, Tuple

from .rcf_prepass import rcf_prepass, _is_x_pm_1
from .modular_helpers import mod_p, omega_matrix, inv_mod_mat, is_symplectic
from .atomic_types import AtomicBlock
from .atomic_linear import split_uv
from .atomic_paired import atomic_blocks_in_paired_sector
from .atomic_self import atomic_blocks_in_self_sector_nonunipotent
from .atomic_unipotent_p2 import atomic_blocks_in_unipotent_self_sector_p2


def _concat_blocks_to_basis(blocks: List[AtomicBlock], n2: int, p: int) -> np.ndarray:
    """
    Build global symplectic basis B = [U_all | V_all] from block hyperbolic bases.
    Requires each block basis is canonical [U_blk | V_blk].
    """
    if not blocks:
        return np.eye(n2, dtype=np.int64)

    U_list = []
    V_list = []
    for b in blocks:
        U_blk, V_blk = split_uv(b.T_blk)
        U_list.append(U_blk)
        V_list.append(V_blk)

    B = np.concatenate(U_list + V_list, axis=1)
    return mod_p(B, p)


def atomic_block_decompose(F: np.ndarray, p: int) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Provably-optimal route (once unipotent p=2 is implemented):
      - CRT primary decomposition + sectors
      - sector-wise atomic indecomposable block construction
      - concatenate into global symplectic B
      - Sigma = B^{-1} F B (symplectic similarity)
    Returns: (Sigma, B, info)
    """
    assert is_symplectic(F, p)
    n2 = F.shape[0]
    n = n2 // 2
    Ω = omega_matrix(n, p)

    meta = rcf_prepass(F, p)
    prim = meta["primaries"]
    sectors = meta["sectors"]

    blocks: List[AtomicBlock] = []
    sector_invariants = []

    for sec in sectors:
        if sec["type"] == "paired":
            key = sec["key"]
            key_star = sec["key_star"]
            b, inv = atomic_blocks_in_paired_sector(F, p, key, key_star, prim)
        else:
            key = sec["key"]
            q = prim[key]["poly"]
            if p == 2 and _is_x_pm_1(q, p):
                b, inv = atomic_blocks_in_unipotent_self_sector_p2(F, key, prim)
            else:
                b, inv = atomic_blocks_in_self_sector_nonunipotent(F, p, key, prim)

        blocks += b
        sector_invariants.append(inv)

    B = _concat_blocks_to_basis(blocks, n2, p)
    if B.shape != (n2, n2):
        raise RuntimeError(f"Global basis has wrong shape {B.shape}; blocks do not span full space.")

    Ω = omega_matrix(n, p)
    if not np.array_equal(mod_p(B.T @ Ω @ B, p), Ω % p):
        raise RuntimeError("Global basis B is not symplectic.")

    Sigma = mod_p(inv_mod_mat(B, p) @ F @ B, p)

    info = {
        "sector_invariants": sector_invariants,
        "atomic_half_dims": [b.half_dim for b in blocks],
        "Q_opt": max([b.half_dim for b in blocks], default=0),
        "certified": all(inv.data.get("status") == "OK" for inv in sector_invariants),
    }
    return Sigma, B, info
