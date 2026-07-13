from __future__ import annotations

import numpy as np
import pytest

from sympleq.core.symmetries.atomic_decomposition import atomic_block_decompose
from sympleq.core.symmetries.atomic_decomposition_helpers._brute_force_oracle import (
    OracleInfeasible,
    brute_force_cost,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2_generators import (
    canonical_unipotent_p2_block,
    canonical_W_block_p2,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import inv_mod_mat, is_symplectic, mod_p, rank_mod


def _companion_multiplication_matrix(q: list[int], p: int) -> np.ndarray:
    q = [int(c) % p for c in q]
    d = len(q) - 1
    C = np.zeros((d, d), dtype=np.int64)
    for j in range(d - 1):
        C[j + 1, j] = 1
    for i in range(d):
        C[i, d - 1] = (-q[i]) % p
    return mod_p(C, p)


def _block_diag(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return np.block([
        [A, np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)],
        [np.zeros((B.shape[0], A.shape[1]), dtype=np.int64), B],
    ])


def _symplectic_scale(A: np.ndarray, p: int) -> np.ndarray:
    A = mod_p(A, p)
    return mod_p(_block_diag(A, inv_mod_mat(A, p).T), p)




def _assert_sigma_block_diagonal_by_half_dims(Sigma: np.ndarray, half_dims: list[int], p: int) -> None:
    Sigma = mod_p(Sigma, p)
    n2 = Sigma.shape[0]
    n = n2 // 2
    assert sum(int(h) for h in half_dims) == n
    idx_sets = []
    off = 0
    for h in half_dims:
        h = int(h)
        idx_sets.append(list(range(off, off + h)) + list(range(n + off, n + off + h)))
        off += h
    for i, idx_i in enumerate(idx_sets):
        for j, idx_j in enumerate(idx_sets):
            if i == j:
                continue
            assert not np.any(Sigma[np.ix_(idx_i, idx_j)] % p)


def _assert_matches_bruteforce(F: np.ndarray, p: int, *, cap: int = 30000) -> None:
    F = mod_p(F, p)
    assert is_symplectic(F, p)
    try:
        exact = brute_force_cost(F, p, cap=cap)
    except OracleInfeasible as exc:
        pytest.skip(f"brute-force oracle exceeded cap: {exc}")

    Sigma, B, info = atomic_block_decompose(F, p)
    verify_global_basis(F, B, Sigma, p)
    _assert_sigma_block_diagonal_by_half_dims(Sigma, [int(h) for h in info["atomic_half_dims"]], p)

    attained = int(info.get("Q_att", info.get("attained_qudit_cost", info["qudit_cost"])))
    assert attained == exact
    assert info["minimal_cost_certified"] is True
    assert info["Q_opt"] == exact
    assert info["cost_certificate"]["certified_minimal"] is True


@pytest.mark.parametrize("p,n", [(2, 1), (2, 2), (3, 1), (3, 2)])
def test_bruteforce_oracle_identity_cases(p: int, n: int) -> None:
    _assert_matches_bruteforce(np.eye(2 * n, dtype=np.int64), p)


@pytest.mark.parametrize(
    "F,p",
    [
        (canonical_unipotent_p2_block("V_beta", 2, beta=0), 2),
        (canonical_unipotent_p2_block("V_beta", 4, beta=1), 2),
        (canonical_W_block_p2(2), 2),
        (canonical_unipotent_p2_block("W_beta", 3, beta=1), 2),
    ],
)
def test_bruteforce_oracle_p2_canonical_blocks(F: np.ndarray, p: int) -> None:
    _assert_matches_bruteforce(F, p)


def test_bruteforce_oracle_odd_p_nonlinear_self_reciprocal_sector() -> None:
    # q=x^2+1 is irreducible and self-reciprocal over GF(3).  The companion
    # matrix has determinant 1, hence is a 2x2 symplectic block.
    p = 3
    F = _companion_multiplication_matrix([1, 0, 1], p)
    assert rank_mod(F, p) == 2
    _assert_matches_bruteforce(F, p)


def test_bruteforce_oracle_odd_p_nonlinear_paired_sector() -> None:
    # q=x^2+x+2 is irreducible over GF(3) and not self-reciprocal.  The
    # symplectic scale diag(A,A^{-T}) realizes the paired sector q + q*.
    p = 3
    A = _companion_multiplication_matrix([2, 1, 1], p)
    F = _symplectic_scale(A, p)
    _assert_matches_bruteforce(F, p)
