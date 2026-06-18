
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple

from sympleq.core.symmetries.modular_helpers import (
    mod_p,
    rank_mod,
    nullspace_mod,
    solve_linear_many,
    independent_columns,
    basis_extend as _basis_extend,
)

from sympleq.core.symmetries.polynomials_fp import poly_eval_matrix


@dataclass(frozen=True)
class PrimaryInfo:
    key: Tuple[int, ...]
    q: np.ndarray
    deg: int
    exp: int
    V_basis: np.ndarray


def restrict_operator_invariant(F: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    B = independent_columns(mod_p(B, p), p)
    FB = mod_p(F @ B, p)
    Fs = solve_linear_many(B, FB, p)
    return mod_p(Fs, p)


def nilpotent_kernel_basis(N: np.ndarray, j: int, p: int) -> np.ndarray:
    d = N.shape[0]
    Nj = np.eye(d, dtype=np.int64)
    for _ in range(j):
        Nj = mod_p(Nj @ N, p)
    return nullspace_mod(Nj, p)


def jordan_chain_tops_nilpotent(N: np.ndarray, max_exp: int, p: int) -> Dict[int, np.ndarray]:
    """
    tops[L] columns are representatives of K_L / (K_{L-1} + N K_{L+1}),
    giving chain tops of exact length L.

    This now delegates to the shared ``build_nilpotent_filtration`` with
    ``space_basis = I`` (full ambient space), removing the previously duplicated
    filtration logic. The result is identical: with ``space_basis = I`` the
    span-restricted kernels reduce to the ordinary kernels of ``N^j``.
    """
    # Local import keeps module import order acyclic (atomic_filtration does not
    # depend on module_invariants).
    from .atomic_filtration import build_nilpotent_filtration
    d = N.shape[0]
    I = np.eye(d, dtype=np.int64)
    return build_nilpotent_filtration(N, I, int(max_exp), p).tops


def cyclic_submodule_basis(Fp: np.ndarray, Np: np.ndarray, v_top: np.ndarray, deg_q: int, L: int, p: int) -> np.ndarray:
    Fp = mod_p(Fp, p)
    Np = mod_p(Np, p)
    d = Fp.shape[0]
    v_top = mod_p(v_top.reshape(-1, 1), p)

    F_pow: List[np.ndarray] = [np.eye(d, dtype=np.int64)]
    for _ in range(1, deg_q):
        F_pow.append(mod_p(F_pow[-1] @ Fp, p))

    cols: List[np.ndarray] = []
    w = v_top.copy()
    for _t in range(L):
        for a in range(deg_q):
            cols.append(mod_p(F_pow[a] @ w, p))
        w = mod_p(Np @ w, p)

    B = np.concatenate(cols, axis=1) if cols else np.zeros((d, 0), dtype=np.int64)
    B = independent_columns(B, p)
    if B.shape[1] != deg_q * L:
        raise RuntimeError("cyclic_submodule_basis: unexpected dimension; top vector selection inconsistent.")
    return B


def q_of_F_restricted(Fp: np.ndarray, q: np.ndarray, p: int) -> np.ndarray:
    return mod_p(poly_eval_matrix(mod_p(Fp, p), q, p), p)
