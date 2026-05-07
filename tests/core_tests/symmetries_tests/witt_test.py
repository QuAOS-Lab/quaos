import numpy as np
import pytest

from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_witt import witt_decompose_form
from sympleq.core.symmetries.modular_helpers import mod_p, rank_mod


def rand_full_rank_matrix(m: int, p: int, rng: np.random.Generator) -> np.ndarray:
    for _ in range(1000):
        A = rng.integers(0, p, size=(m, m), dtype=np.int64)
        if rank_mod(A, p) == m:
            return A
    raise RuntimeError("failed to generate full-rank matrix")


class TestWittDecomposition:
    @pytest.mark.parametrize("p,m,seed", [(3, 7, 1), (5, 8, 2), (7, 9, 3)])
    def test_witt_symmetric_congruence_is_valid(self, p: int, m: int, seed: int) -> None:
        rng = np.random.default_rng(seed)
        A = rand_full_rank_matrix(m, p, rng)
        B = mod_p(A + A.T, p)
        P, blocks, info = witt_decompose_form(B, p)
        Bp = mod_p(P.T @ B @ P, p)
        assert rank_mod(P, p) == m
        assert rank_mod(Bp, p) == rank_mod(B, p)
        used = set()
        for typ, i, j in blocks:
            assert typ in {"hyp", "ani", "rad"}
            used.add(i)
            if typ == "hyp":
                used.add(j)
        assert used == set(range(m))
        assert info["rank"] == rank_mod(B, p)
        assert info["form"] == "symmetric"

    @pytest.mark.parametrize("p,m,seed", [(3, 8, 1), (5, 10, 2)])
    def test_witt_alternating_decomposes_into_hyperbolic_planes(self, p: int, m: int, seed: int) -> None:
        rng = np.random.default_rng(seed)
        A = rand_full_rank_matrix(m, p, rng)
        B = mod_p(A - A.T, p)
        assert np.all(np.diag(B) % p == 0)
        P, blocks, info = witt_decompose_form(B, p)
        Bp = mod_p(P.T @ B @ P, p)
        assert rank_mod(P, p) == m
        assert rank_mod(Bp, p) == rank_mod(B, p)
        r = rank_mod(B, p)
        assert r % 2 == 0
        assert all(typ in {"hyp", "rad"} for (typ, _, __) in blocks)
        assert sum(typ == "hyp" for (typ, _, __) in blocks) == r // 2
        assert sum(typ == "rad" for (typ, _, __) in blocks) == m - r
        assert info["form"] == "alternating"

    def test_witt_rejects_p2_bilinear_shortcut(self) -> None:
        B = np.array([[0, 1], [1, 0]], dtype=np.int64)
        with pytest.raises(NotImplementedError):
            witt_decompose_form(B, 2)

    def test_witt_rejects_general_nonsymmetric_nonalternating_form(self) -> None:
        B = np.array([[1, 1], [0, 1]], dtype=np.int64)
        with pytest.raises(ValueError):
            witt_decompose_form(B, 5)
