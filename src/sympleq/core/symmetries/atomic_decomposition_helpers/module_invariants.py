
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Tuple

from sympleq.core.symmetries.modular_helpers import (
    mod_p,
    rank_mod,
    nullspace_mod,
    solve_linear_many,
    independent_columns,
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


def _basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    base = independent_columns(base, p) if base.size else base
    picked = np.zeros((candidates.shape[0], 0), dtype=np.int64)
    r_base = rank_mod(base, p) if base.size else 0

    for j in range(candidates.shape[1]):
        c = candidates[:, j:j + 1]
        r_try = rank_mod(np.concatenate([base, picked, c], axis=1), p)
        if r_try > r_base + picked.shape[1]:
            picked = np.concatenate([picked, c], axis=1)
            if picked.shape[1] == want:
                return picked

    raise RuntimeError("_basis_extend: could not extend by required amount.")


def jordan_chain_tops_nilpotent(N: np.ndarray, max_exp: int, p: int) -> Dict[int, np.ndarray]:
    """
    tops[L] columns are representatives of K_L / (K_{L-1} + N K_{L+1}),
    giving chain tops of exact length L.
    """
    d = N.shape[0]
    K: List[np.ndarray] = [np.zeros((d, 0), dtype=np.int64)]
    for j in range(1, max_exp + 1):
        K.append(independent_columns(nilpotent_kernel_basis(N, j, p), p))
    K.append(K[max_exp])  # K_{max+1} := K_max

    tops: Dict[int, np.ndarray] = {}
    for L in range(1, max_exp + 1):
        KL = K[L]
        if KL.shape[1] == 0:
            continue

        S = K[L - 1]
        NKLp1 = mod_p(N @ K[L + 1], p) if K[L + 1].shape[1] else np.zeros((d, 0), dtype=np.int64)
        if NKLp1.shape[1]:
            S = np.concatenate([S, NKLp1], axis=1) if S.shape[1] else NKLp1
        S = independent_columns(S, p) if S.shape[1] else S

        rS = rank_mod(S, p) if S.shape[1] else 0
        need = KL.shape[1] - rS
        if need <= 0:
            continue

        chosen = _basis_extend(S, KL, need, p)
        tops[L] = chosen

    return tops


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
