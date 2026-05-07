# sympleq/core/symmetries/atomic_decomposition_helpers/atomic_filtration.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np

from ..modular_helpers import mod_p, independent_columns, rank_mod, nullspace_mod
from .atomic_linear import mat_pow_mod, kernel_in_span


@dataclass(frozen=True, slots=True)
class NilpotentFiltration:
    """
    Cached nilpotent kernel filtration for a fixed nilpotent operator N on a
    fixed invariant working subspace.

    K[j] is an ambient-coordinate basis for ker(N^j) ∩ span(space_basis).
    denom[L] is K[L-1] + N K[L+1].
    tops[L] is a deterministic set of representatives for K[L]/denom[L].
    """

    N: np.ndarray
    p: int
    max_exp: int
    space_basis: np.ndarray
    K: Dict[int, np.ndarray]
    denom: Dict[int, np.ndarray]
    tops: Dict[int, np.ndarray]


def _basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    base = independent_columns(mod_p(base, p), p) if base.size else base
    picked = np.zeros((candidates.shape[0], 0), dtype=np.int64)
    r_base = rank_mod(base, p) if base.size else 0
    for j in range(candidates.shape[1]):
        c = mod_p(candidates[:, j:j + 1], p)
        r_try = rank_mod(np.concatenate([base, picked, c], axis=1), p)
        if r_try > r_base + picked.shape[1]:
            picked = np.concatenate([picked, c], axis=1)
            if picked.shape[1] == want:
                return picked
    raise RuntimeError("_basis_extend: could not extend by required amount.")


def build_nilpotent_filtration(N: np.ndarray, space_basis: np.ndarray, max_exp: int, p: int) -> NilpotentFiltration:
    """Build and cache K_j, quotient denominators, and top representatives."""
    N = mod_p(np.asarray(N, dtype=np.int64), p)
    space_basis = independent_columns(mod_p(space_basis, p), p)
    d = N.shape[0]
    max_exp = int(max_exp)

    K: Dict[int, np.ndarray] = {0: np.zeros((d, 0), dtype=np.int64)}
    for j in range(1, max_exp + 1):
        K[j] = kernel_in_span(mat_pow_mod(N, j, p), space_basis, p)
    K[max_exp + 1] = K[max_exp]

    denom: Dict[int, np.ndarray] = {}
    tops: Dict[int, np.ndarray] = {}
    for L in range(1, max_exp + 1):
        KL = K[L]
        if KL.shape[1] == 0:
            denom[L] = np.zeros((d, 0), dtype=np.int64)
            continue
        D = K[L - 1]
        NK_next = mod_p(N @ K[L + 1], p) if K[L + 1].shape[1] else np.zeros((d, 0), dtype=np.int64)
        if NK_next.shape[1]:
            D = np.concatenate([D, NK_next], axis=1) if D.shape[1] else NK_next
        D = independent_columns(D, p) if D.shape[1] else D
        denom[L] = D
        rD = rank_mod(D, p) if D.shape[1] else 0
        need = KL.shape[1] - rD
        if need > 0:
            tops[L] = _basis_extend(D, KL, need, p)

    return NilpotentFiltration(
        N=N,
        p=int(p),
        max_exp=max_exp,
        space_basis=space_basis,
        K=K,
        denom=denom,
        tops=tops,
    )
