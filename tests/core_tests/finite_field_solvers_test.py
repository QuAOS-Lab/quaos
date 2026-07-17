import numpy as np
import pytest

from sympleq.core.finite_field_solvers import gf_inv, gf_lu, gf_rref


def _random_invertible_matrix(p: int, size: int, rng: np.random.Generator) -> np.ndarray:
    while True:
        matrix = rng.integers(0, p, size=(size, size), dtype=int)
        try:
            gf_inv(matrix, p=p)
        except ValueError:
            continue
        return matrix


def _random_matrix(p: int, shape: tuple[int, int], rng: np.random.Generator) -> np.ndarray:
    return rng.integers(0, p, size=shape, dtype=int)


def _is_permutation_matrix(P: np.ndarray) -> bool:
    return bool(
        P.ndim == 2
        and P.shape[0] == P.shape[1]
        and np.all((P == 0) | (P == 1))
        and np.all(P.sum(axis=0) == 1)
        and np.all(P.sum(axis=1) == 1)
    )


@pytest.mark.parametrize("p", [2, 3, 5, 7])
@pytest.mark.parametrize("n", [1, 2, 4])
def test_gf_inv_returns_inverse(p: int, n: int):
    rng = np.random.default_rng(seed=1000 + 10 * p + n)
    A = _random_invertible_matrix(p, n, rng)

    inv = gf_inv(A, p=p)

    assert np.array_equal((A @ inv) % p, np.eye(n, dtype=int) % p)


@pytest.mark.parametrize("p", [2, 3, 5])
def test_gf_rref_reconstructs_reduced_matrix(p: int):
    rng = np.random.default_rng(seed=2000 + p)
    A = _random_matrix(p, (4, 6), rng)

    R, M, N, rank = gf_rref(A, p=p)

    assert np.array_equal(((M @ A) % p @ N) % p, R)
    pivots = []
    for row in R:
        nz = np.nonzero(row)[0]
        if nz.size == 0:
            assert np.all(row % p == 0)
            continue
        pivot_col = int(nz[0])
        pivots.append(pivot_col)
        assert row[pivot_col] % p == 1
        assert np.all(row[:pivot_col] % p == 0)
        assert np.all(row[pivot_col + 1:] % p == 0)
    assert rank == len(pivots)


@pytest.mark.parametrize("p", [2, 5, 11])
@pytest.mark.parametrize("n", [2, 3, 5])
def test_gf_lu_reconstructs_permuted_matrix(p: int, n: int):
    rng = np.random.default_rng(seed=3000 + 10 * p + n)
    A = _random_matrix(p, (n, n), rng)

    L, U, P = gf_lu(A, p=p)

    assert _is_permutation_matrix(P)
    assert np.array_equal((P @ A) % p, (L @ U) % p)
    assert np.array_equal(np.diag(L) % p, np.ones(n, dtype=int))
    assert np.all((np.triu(L, k=1) % p) == 0)
