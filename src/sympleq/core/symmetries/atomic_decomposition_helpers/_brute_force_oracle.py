# sympleq/core/symmetries/atomic_decomposition_helpers/_brute_force_oracle.py
"""
Theory-free brute-force oracle for the exact qudit cost (Definition 1.1).

This deliberately contains *none* of the Chapter-5 classification: it computes
the exact minimum of ``max_i k_i`` over invariant symplectic decompositions by
exhaustive search, so it is an independent ground truth for the pipeline and the
closed-form cost formulas (review sec. 9.2).

Method.  Every F-invariant subspace is a sum of cyclic (Krylov) submodules, so
the full lattice of invariant subspaces is the sum-closure of the cyclic ones.
We enumerate that lattice (with a hard cap to stay finite on degenerate inputs
such as the identity), keep the nondegenerate ones as candidate blocks, and
recurse:

    cost(U) = min over nondegenerate invariant W <= U of max(dim W / 2, cost(W^perp ∩ U))

with ``cost(0) = 0``, branch-and-bound (smallest blocks first, prune when a
block is already as large as the best found), and memoization on the canonical
column space.  Only generic GF(p) linear algebra is used.

Intended for 2n <= 8 over small fields.
"""
from __future__ import annotations

from itertools import product
from typing import Dict, List, Tuple

import numpy as np

from ..modular_helpers import mod_p, rank_mod, nullspace_mod, omega_matrix, rref_mod


class OracleInfeasible(RuntimeError):
    """Raised when the invariant-subspace lattice exceeds the search cap."""


def _rref_rows(M: np.ndarray, p: int) -> np.ndarray:
    R, piv = rref_mod(mod_p(M, p), p)
    r = len(piv)
    return mod_p(R[:r, :], p)


def _colspace_basis(B: np.ndarray, p: int) -> np.ndarray:
    """Canonical basis of the column space, returned as columns."""
    B = mod_p(B, p)
    if B.shape[1] == 0:
        return B
    rows = _rref_rows(B.T, p)          # canonical basis of colspace, as rows
    return mod_p(rows.T, p)


def _key(B: np.ndarray, p: int) -> Tuple:
    rows = _rref_rows(B.T, p)
    return (int(B.shape[0]), int(rows.shape[0]), tuple(map(tuple, rows.tolist())))


def _krylov(F: np.ndarray, v: np.ndarray, p: int) -> np.ndarray:
    F = mod_p(F, p)
    w = mod_p(v.reshape(-1, 1), p)
    n = F.shape[0]
    cols = []
    for _ in range(n + 1):
        cols.append(w.copy())
        w = mod_p(F @ w, p)
    return _colspace_basis(np.concatenate(cols, axis=1), p)


def invariant_subspaces(F: np.ndarray, p: int, cap: int = 5000) -> List[np.ndarray]:
    """All F-invariant subspaces (as column bases), via sum-closure of cyclics."""
    F = mod_p(F, p)
    n = F.shape[0]
    cyc: Dict[Tuple, np.ndarray] = {}
    for bits in product(range(p), repeat=n):
        if not any(bits):
            continue
        C = _krylov(F, np.array(bits, dtype=np.int64), p)
        cyc[_key(C, p)] = C

    zero = np.zeros((n, 0), dtype=np.int64)
    lattice: Dict[Tuple, np.ndarray] = {_key(zero, p): zero}
    for k, C in cyc.items():
        lattice[k] = C
    cyc_list = list(cyc.values())

    # Worklist closure: every invariant subspace is reached by repeatedly adding
    # one cyclic generator, so processing each new subspace once (summed with each
    # cyclic) generates the whole lattice without the all-pairs-per-pass blowup.
    work = list(cyc.values())
    while work:
        A = work.pop()
        for C in cyc_list:
            Sb = _colspace_basis(np.concatenate([A, C], axis=1), p)
            k = _key(Sb, p)
            if k not in lattice:
                lattice[k] = Sb
                work.append(Sb)
                if len(lattice) > cap:
                    raise OracleInfeasible(
                        f"invariant-subspace lattice exceeded cap={cap} (degenerate F?)"
                    )
    return list(lattice.values())


def _is_nondeg(W: np.ndarray, Omega: np.ndarray, p: int) -> bool:
    if W.shape[1] == 0:
        return False
    G = mod_p(W.T @ Omega @ W, p)
    return rank_mod(G, p) == W.shape[1]


def _subspace_leq(W: np.ndarray, U: np.ndarray, p: int) -> bool:
    if W.shape[1] == 0:
        return True
    return rank_mod(np.concatenate([U, W], axis=1), p) == rank_mod(U, p)


def _perp_within(W: np.ndarray, U: np.ndarray, Omega: np.ndarray, p: int) -> np.ndarray:
    """Basis of { x in colspace(U) : W^T Omega x = 0 }."""
    U = mod_p(U, p)
    if U.shape[1] == 0:
        return U
    A = mod_p(W.T @ Omega @ U, p)          # constraint A c = 0, with x = U c
    Ns = nullspace_mod(A, p)               # (dimU x k)
    return _colspace_basis(mod_p(U @ Ns, p), p)


def brute_force_cost(F: np.ndarray, p: int, cap: int = 5000) -> int:
    """Exact qudit cost of F by exhaustive search. Raises OracleInfeasible if
    the invariant-subspace lattice is too large (e.g. near-identity at large n)."""
    F = mod_p(F, p)
    n2 = F.shape[0]
    if n2 == 0:
        return 0
    Omega = omega_matrix(n2 // 2, p)

    subs = invariant_subspaces(F, p, cap=cap)
    nd: List[Tuple[np.ndarray, int]] = sorted(
        ((W, W.shape[1]) for W in subs if _is_nondeg(W, Omega, p)),
        key=lambda t: t[1],
    )

    memo: Dict[Tuple, int] = {}

    def cost(U: np.ndarray) -> int:
        ku = _key(U, p)
        if ku in memo:
            return memo[ku]
        dimU = U.shape[1]
        if dimU == 0:
            memo[ku] = 0
            return 0
        best = dimU // 2  # U itself is a valid (possibly decomposable) block
        for W, dW in nd:
            k = dW // 2
            if k >= best:           # sorted ascending => nothing smaller-max remains
                break
            if not _subspace_leq(W, U, p):
                continue
            comp = _perp_within(W, U, Omega, p)
            if comp.shape[1] != dimU - dW:
                continue
            cand = max(k, cost(comp))
            if cand < best:
                best = cand
                if best == 1:
                    break
        memo[ku] = best
        return best

    whole = _colspace_basis(np.eye(n2, dtype=np.int64), p)
    return cost(whole)
