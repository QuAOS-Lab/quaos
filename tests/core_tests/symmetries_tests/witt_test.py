import numpy as np
import pytest

from sympleq.core.symmetries.modular_helpers import mod_p, rank_mod
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_witt import witt_decompose_form


def rand_full_rank_matrix(m, p, rng):
    while True:
        A = rng.integers(0, p, size=(m, m), dtype=np.int64)
        if rank_mod(A, p) == m:
            return A


class TestWittDecomposition:

    @pytest.mark.parametrize("p,m,seed", [(3, 7, 1), (5, 8, 2), (7, 9, 3)])
    def test_witt_symmetric_congruence_is_valid(self, p, m, seed):
        rng = np.random.default_rng(seed)
        A = rand_full_rank_matrix(m, p, rng)
        B = mod_p(A + A.T, p)  # symmetric

        P, blocks, info = witt_decompose_form(B, p)
        Bp = mod_p(P.T @ B @ P, p)

        # P invertible
        assert rank_mod(P, p) == m

        # congruence preserved rank
        assert rank_mod(Bp, p) == rank_mod(B, p)

        # Ensure blocks list covers all dims (no missing indices)
        used = set()
        for t, i, j in blocks:
            if t == "hyp":
                used.add(i)
                used.add(j)
            else:
                used.add(i)
        assert used == set(range(m))

    @pytest.mark.parametrize("p,m,seed", [(3, 8, 1), (5, 10, 2)])
    def test_witt_alternating_decomposes_into_hyperbolic_planes(self, p, m, seed):
        rng = np.random.default_rng(seed)
        # build alternating full rank: B = A - A^T
        A = rand_full_rank_matrix(m, p, rng)
        B = mod_p(A - A.T, p)
        assert np.all(np.diag(B) % p == 0)

        P, blocks, info = witt_decompose_form(B, p)
        Bp = mod_p(P.T @ B @ P, p)

        assert rank_mod(P, p) == m
        assert rank_mod(Bp, p) == rank_mod(B, p)

        # for full-rank alternating, expect only hyperbolic blocks and no radical
        r = rank_mod(B, p)
        assert r % 2 == 0
        assert all(t in ("hyp", "rad") for (t, _, __) in blocks)
        assert sum(t == "hyp" for (t, _, __) in blocks) == r // 2
        assert sum(t == "rad" for (t, _, __) in blocks) == m - r

