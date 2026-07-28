from __future__ import annotations

import numpy as np
import pytest

from sympleq.core.minimal_qudit_frame.helpers.mqf_extension import top_quotient_E_basis
from sympleq.core.minimal_qudit_frame.helpers.mqf_filtration import build_nilpotent_filtration
from sympleq.core.minimal_qudit_frame.helpers.mqf_krylov import (
    select_module_generators_from_top_quotient,
)
from sympleq.core.minimal_qudit_frame.helpers.module_invariants import cyclic_submodule_basis
from sympleq.core.symmetries.modular_helpers import mod_p, rank_mod
from sympleq.core.symmetries.polynomials_fp import poly_eval_matrix, poly_pow


def _companion_multiplication_matrix(poly: np.ndarray, p: int) -> np.ndarray:
    """Companion matrix as multiplication by x modulo poly, column convention."""
    coeff = [int(c) % p for c in np.asarray(poly, dtype=np.int64).reshape(-1)]
    d = len(coeff) - 1
    C = np.zeros((d, d), dtype=np.int64)
    for j in range(d - 1):
        C[j + 1, j] = 1
    for i in range(d):
        C[i, d - 1] = (-coeff[i]) % p
    return mod_p(C, p)


def _block_diag(mats: list[np.ndarray], p: int) -> np.ndarray:
    n = sum(M.shape[0] for M in mats)
    out = np.zeros((n, n), dtype=np.int64)
    off = 0
    for M in mats:
        d = M.shape[0]
        out[off:off + d, off:off + d] = M
        off += d
    return mod_p(out, p)


def _rand_invertible(rng: np.random.Generator, m: int, p: int) -> np.ndarray:
    for _ in range(200):
        A = rng.integers(0, p, size=(m, m), dtype=np.int64)
        if rank_mod(A, p) == m:
            return mod_p(A, p)
    raise RuntimeError("could not sample invertible matrix")


@pytest.mark.parametrize(
    "p,q,L,mult,seed",
    [
        (3, [1, 0, 1], 2, 2, 10),       # two q^2 blocks, deg(q)=2
        (5, [1, 1, 1], 2, 3, 11),       # three q^2 blocks, deg(q)=2
        (2, [1, 1, 0, 1], 2, 2, 12),    # two q^2 blocks, deg(q)=3
    ],
)
def test_top_quotient_E_basis_survives_adversarial_base_field_basis(
    p: int, q: list[int], L: int, mult: int, seed: int
) -> None:
    q_arr = np.asarray(q, dtype=np.int64)
    deg_q = len(q) - 1
    block_poly = poly_pow(q_arr, L, p)
    F_block = _companion_multiplication_matrix(block_poly, p)
    F = _block_diag([F_block.copy() for _ in range(mult)], p)
    N = mod_p(poly_eval_matrix(F, q_arr, p), p)

    filt = build_nilpotent_filtration(N, np.eye(F.shape[0], dtype=np.int64), L, p)
    top = filt.tops[L]
    denom = filt.denom[L]

    assert top.shape[1] == deg_q * mult
    assert rank_mod(np.concatenate([denom, top], axis=1), p) == rank_mod(denom, p) + deg_q * mult

    rng = np.random.default_rng(seed)
    scramble = _rand_invertible(rng, top.shape[1], p)
    top_scrambled = mod_p(top @ scramble, p)

    tq = top_quotient_E_basis(F, N, L, deg_q, p, top_reps=top_scrambled, denom=denom)
    assert tq.dim_fp == deg_q * mult
    assert tq.dim_E == mult
    assert len(tq.e_basis_reps) == mult
    assert rank_mod(np.concatenate([denom, tq.e_orbit_basis], axis=1), p) == rank_mod(denom, p) + deg_q * mult

    gens = select_module_generators_from_top_quotient(F, N, top_scrambled, deg_q, L, p, denom=denom)
    assert gens.shape[1] == mult

    cyclic_parts = []
    for j in range(gens.shape[1]):
        C = cyclic_submodule_basis(F, N, gens[:, j:j + 1], deg_q, L, p)
        assert C.shape[1] == deg_q * L
        cyclic_parts.append(C)

    total = np.concatenate(cyclic_parts, axis=1)
    assert rank_mod(total, p) == deg_q * L * mult
