import numpy as np
import pytest

from sympleq.core.minimal_qudit_frame import minimal_qudit_frame, MinimalQuditFrameCertificationError
from sympleq.core.minimal_qudit_frame.helpers.mqf_linear import symplectic_left_inverse
from sympleq.core.minimal_qudit_frame.helpers.mqf_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import (
    inv_mod_mat,
    is_symplectic,
    mod_p,
    omega_matrix,
    rank_mod,
)


def _rand_mat(rng: np.random.Generator, n: int, m: int, p: int) -> np.ndarray:
    return rng.integers(0, p, size=(n, m), dtype=np.int64)


def _rand_invertible(rng: np.random.Generator, n: int, p: int, max_tries: int = 500) -> np.ndarray:
    for _ in range(max_tries):
        A = _rand_mat(rng, n, n, p)
        if rank_mod(A, p) == n:
            return mod_p(A, p)
    raise RuntimeError("could not sample invertible matrix")


def _rand_symmetric(rng: np.random.Generator, n: int, p: int) -> np.ndarray:
    M = _rand_mat(rng, n, n, p)
    return mod_p(M + M.T, p) if p == 2 else mod_p((M + M.T) * pow(2, -1, p), p)


def _block_diag(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return np.block([
        [A, np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)],
        [np.zeros((B.shape[0], A.shape[1]), dtype=np.int64), B],
    ])


def _symplectic_scale(A: np.ndarray, p: int) -> np.ndarray:
    A = mod_p(A, p)
    return mod_p(_block_diag(A, inv_mod_mat(A, p).T), p)


def _symplectic_shear_upper(B: np.ndarray, p: int) -> np.ndarray:
    n = B.shape[0]
    I = np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    return mod_p(np.block([[I, B], [Z, I]]), p)


def _symplectic_shear_lower(C: np.ndarray, p: int) -> np.ndarray:
    n = C.shape[0]
    I = np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    return mod_p(np.block([[I, Z], [C, I]]), p)


def rand_symplectic(rng: np.random.Generator, n: int, p: int, steps: int = 12) -> np.ndarray:
    F = np.eye(2 * n, dtype=np.int64)
    for _ in range(steps):
        kind = int(rng.integers(0, 3))
        if kind == 0:
            S = _symplectic_shear_upper(_rand_symmetric(rng, n, p), p)
        elif kind == 1:
            S = _symplectic_shear_lower(_rand_symmetric(rng, n, p), p)
        else:
            S = _symplectic_scale(_rand_invertible(rng, n, p), p)
        F = mod_p(S @ F, p)
    assert is_symplectic(F, p)
    return F


def _extract_MQF_blocks_from_global_basis(B: np.ndarray, half_dims: list[int], p: int) -> list[np.ndarray]:
    """Reconstruct T_blk=[U_blk|V_blk] from B=[U_all|V_all]."""
    B = mod_p(B, p)
    n2 = B.shape[0]
    n = n2 // 2
    assert sum(half_dims) == n
    U_all = B[:, :n]
    V_all = B[:, n:]
    out = []
    off = 0
    for h in half_dims:
        U = U_all[:, off:off + h]
        V = V_all[:, off:off + h]
        out.append(mod_p(np.concatenate([U, V], axis=1), p))
        off += h
    return out


def _assert_darboux_block(T: np.ndarray, p: int) -> None:
    n2, m2 = T.shape
    assert m2 % 2 == 0
    G = mod_p(T.T @ omega_matrix(n2 // 2, p) @ T, p)
    assert np.array_equal(G, omega_matrix(m2 // 2, p))


def _assert_invariant_span(F: np.ndarray, T: np.ndarray, p: int) -> None:
    T = mod_p(T, p)
    FT = mod_p(F @ T, p)
    assert rank_mod(np.concatenate([T, FT], axis=1), p) == T.shape[1]


class TestMQFBlocks:
    def test_paired_only_certified_blocks_are_invariant_and_nondegenerate(self) -> None:
        p, n = 5, 4
        # a != a^{-1}, so this is a single paired sector.
        A = 2 * np.eye(n, dtype=np.int64)
        F = _symplectic_scale(A, p)
        Sigma, B, info = minimal_qudit_frame(F, p)

        assert info["status"] == "OK"
        assert info["certified"] is True
        assert info["minimal_cost_certified"] is True
        verify_global_basis(F, B, Sigma, p)

        blocks = _extract_MQF_blocks_from_global_basis(B, info["mqf_half_dims"], p)
        assert len(blocks) >= 1
        assert sum(info["mqf_half_dims"]) == n
        for T in blocks:
            _assert_darboux_block(T, p)
            _assert_invariant_span(F, T, p)

    def test_best_effort_random_blocks_are_valid_when_MQF_cover_is_returned(self) -> None:
        rng = np.random.default_rng(123)
        for p in [2, 3, 5]:
            for n in [2, 3]:
                F = rand_symplectic(rng, n, p, steps=10)
                Sigma, B, info = minimal_qudit_frame(F, p, allow_degraded=True)
                verify_global_basis(F, B, Sigma, p)
                assert info["status"] in {"OK", "DEGRADED"}
                assert info["qudit_cost"] == max(info["mqf_half_dims"], default=0)

                # If no global completion was used, B is exactly the returned block frame.
                if not info.get("completed", False):
                    blocks = _extract_MQF_blocks_from_global_basis(B, info["mqf_half_dims"], p)
                    for T in blocks:
                        _assert_darboux_block(T, p)
                        _assert_invariant_span(F, T, p)

    def test_certified_mode_is_deterministic_on_paired_only_case(self) -> None:
        p, n = 7, 3
        F = _symplectic_scale(3 * np.eye(n, dtype=np.int64), p)
        out1 = minimal_qudit_frame(F, p)
        out2 = minimal_qudit_frame(F, p)
        Sigma1, B1, info1 = out1
        Sigma2, B2, info2 = out2
        assert np.array_equal(Sigma1, Sigma2)
        assert np.array_equal(B1, B2)
        assert info1["mqf_half_dims"] == info2["mqf_half_dims"]
        assert info1["cost_certificate"] == info2["cost_certificate"]

    def test_certified_failures_are_certification_errors(self) -> None:
        rng = np.random.default_rng(456)
        F = rand_symplectic(rng, 3, 3, steps=10)
        try:
            minimal_qudit_frame(F, 3)
        except Exception as exc:
            assert isinstance(exc, MinimalQuditFrameCertificationError)
            assert isinstance(exc.info, dict)
            assert "failures" in exc.info or "mqf_half_dims" in exc.info or "error" in exc.info
