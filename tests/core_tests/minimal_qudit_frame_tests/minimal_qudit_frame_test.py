import numpy as np
import pytest

from sympleq.core.minimal_qudit_frame import minimal_qudit_frame, MinimalQuditFrameCertificationError
from sympleq.core.minimal_qudit_frame.block_structure import (
    block_decompose,
    block_decompose_certified,
    block_decompose_optimal,
)
from sympleq.core.minimal_qudit_frame.helpers.mqf_linear import symplectic_left_inverse
from sympleq.core.minimal_qudit_frame.helpers.mqf_verify import (
    verify_global_basis,
    verify_cost_certificate,
)
from sympleq.core.minimal_qudit_frame.helpers.rcf_prepass import rcf_prepass, _is_x_pm_1
from sympleq.core.symmetries.modular_helpers import (
    independent_columns,
    inv_mod_mat,
    is_symplectic,
    mod_p,
    omega_matrix,
    rank_mod,
)
from sympleq.core.symmetries.polynomials_fp import poly_monic, poly_reciprocal


def _rand_mat(rng: np.random.Generator, n: int, m: int, p: int) -> np.ndarray:
    return rng.integers(0, p, size=(n, m), dtype=np.int64)


def _rand_invertible(rng: np.random.Generator, n: int, p: int, max_tries: int = 500) -> np.ndarray:
    for _ in range(max_tries):
        A = _rand_mat(rng, n, n, p)
        if rank_mod(A, p) == n:
            return mod_p(A, p)
    raise RuntimeError("could not sample invertible matrix")


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
    return mod_p(np.block([[I, mod_p(B, p)], [Z, I]]), p)


def _symplectic_shear_lower(C: np.ndarray, p: int) -> np.ndarray:
    n = C.shape[0]
    I = np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    return mod_p(np.block([[I, Z], [mod_p(C, p), I]]), p)


def _rand_symmetric(rng: np.random.Generator, n: int, p: int) -> np.ndarray:
    M = _rand_mat(rng, n, n, p)
    if p == 2:
        # In characteristic 2 the diagonal is allowed for symplectic shears.
        return mod_p(M + M.T + np.diag(np.diag(M)), p)
    inv2 = pow(2, -1, p)
    return mod_p((M + M.T) * inv2, p)


def rand_symplectic(rng: np.random.Generator, n: int, p: int, steps: int = 12) -> np.ndarray:
    F = np.eye(2 * n, dtype=np.int64)
    for _ in range(int(steps)):
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


def _row_is_symplectic(B: np.ndarray, p: int) -> bool:
    n = B.shape[0] // 2
    Ω = omega_matrix(n, p)
    return np.array_equal(mod_p(B @ Ω @ B.T, p), Ω)


def _attained_cost(info: dict) -> int:
    return int(info.get("Q_att", info.get("attained_qudit_cost", info["qudit_cost"])))


def _assert_cost_semantics(info: dict) -> None:
    assert _attained_cost(info) == int(info["qudit_cost"])
    if info.get("minimal_cost_certified", False):
        assert info.get("Q_opt") == _attained_cost(info)
        assert info.get("optimal_qudit_cost") == _attained_cost(info)
    else:
        assert info.get("Q_opt") is None
        assert info.get("optimal_qudit_cost") is None


def _sector_signature(meta: dict) -> list[tuple]:
    sig = []
    for ctx in meta.get("sector_contexts", []):
        if ctx.sector_type == "paired":
            sig.append(("paired", tuple(ctx.sector_key), tuple(ctx.sector_key_star), ctx.deg_q, ctx.max_exp))
        else:
            ctx_meta = ctx.meta or {}
            sig.append(("self", tuple(ctx.sector_key), ctx.deg_q, ctx.max_exp, ctx_meta.get("floor_certified")))
    return sorted(sig, key=str)


class TestRCFPrepass:
    def test_rcf_prepass_returns_context_fields(self) -> None:
        p, n = 5, 2
        F = _symplectic_scale(np.diag([2, 3]).astype(np.int64), p)
        meta = rcf_prepass(F, p)
        assert "primaries" in meta
        assert "sector_contexts" in meta
        assert "prepass_context" in meta
        assert "sectors" not in meta
        ctx = meta["sector_contexts"][0]
        assert ctx.p == p
        assert ctx.T_sec is not None
        assert ctx.sector_type in {"paired", "self"}

    def test_rcf_prepass_paired_linear_factors(self) -> None:
        p = 5
        A = np.diag([2, 3]).astype(np.int64)
        F = _symplectic_scale(A, p)
        meta = rcf_prepass(F, p)
        ctxs = meta["sector_contexts"]
        assert len(ctxs) == 1
        assert ctxs[0].sector_type == "paired"
        k2 = tuple(poly_monic(np.array([p - 2, 1], dtype=np.int64), p).tolist())
        k3 = tuple(poly_monic(np.array([p - 3, 1], dtype=np.int64), p).tolist())
        assert k2 in meta["primaries"] and k3 in meta["primaries"]
        assert tuple(poly_reciprocal(meta["primaries"][k2]["poly"], p).tolist()) == k3

    def test_rcf_prepass_p2_unipotent_floor_is_conservative(self) -> None:
        p, n = 2, 3
        F = np.eye(2 * n, dtype=np.int64)
        meta = rcf_prepass(F, p)
        assert len(meta["sector_contexts"]) == 1
        ctx = meta["sector_contexts"][0]
        assert ctx.sector_type == "self"
        assert ctx.meta["floor_certified"] is False
        assert meta["Lmin_star"] == 1

    def test_rcf_prepass_is_conjugacy_invariant_at_sector_signature_level(self) -> None:
        rng = np.random.default_rng(10)
        p, n = 5, 3
        F = _symplectic_scale(2 * np.eye(n, dtype=np.int64), p)
        S = rand_symplectic(rng, n, p, steps=8)
        F2 = mod_p(inv_mod_mat(S, p) @ F @ S, p)
        assert _sector_signature(rcf_prepass(F, p)) == _sector_signature(rcf_prepass(F2, p))


class TestMinimalQuditFrameDecompositionAPI:
    def test_certified_paired_only_has_minimal_cost_certificate(self) -> None:
        p, n = 5, 3
        F = _symplectic_scale(2 * np.eye(n, dtype=np.int64), p)
        Sigma, B, info = minimal_qudit_frame(F, p)
        verify_global_basis(F, B, Sigma, p)
        assert info["status"] == "OK"
        assert info["certified"] is True
        assert info["minimal_cost_certified"] is True
        assert info["certified_lower_bound"] == info["qudit_cost"] == info["Q_opt"]
        assert info["cost_certificate"]["complete"] is True

    def test_best_effort_random_returns_valid_global_decomposition(self) -> None:
        rng = np.random.default_rng(11)
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                F = rand_symplectic(rng, n, p, steps=8)
                Sigma, B, info = minimal_qudit_frame(F, p, allow_degraded=True)
                verify_global_basis(F, B, Sigma, p)
                assert info["status"] in {"OK", "DEGRADED"}
                _assert_cost_semantics(info)
                if info["status"] != "OK" or info.get("completed"):
                    assert info["minimal_cost_certified"] is False

    def test_auto_mode_falls_back_but_records_last_error_when_needed(self) -> None:
        rng = np.random.default_rng(12)
        F = rand_symplectic(rng, 3, 3, steps=8)
        Sigma, B, info = minimal_qudit_frame(F, 3, allow_degraded=True)
        verify_global_basis(F, B, Sigma, 3)
        assert info["status"] in {"OK", "DEGRADED"}
        if info["status"] == "DEGRADED":
            assert "warnings" in info
            assert info["minimal_cost_certified"] is False

    def test_row_convention_returns_row_action_conjugacy(self) -> None:
        p, n = 5, 2
        F_col = _symplectic_scale(2 * np.eye(n, dtype=np.int64), p)
        F_row = F_col.T
        Sigma_row, B_row, info = minimal_qudit_frame(F_row, p, convention="row")
        assert info["input_convention"] == "row"
        assert _row_is_symplectic(B_row, p)
        # Row convention relation: Sigma = B F B^{-1}.
        assert np.array_equal(mod_p(B_row @ F_row @ inv_mod_mat(B_row, p), p), Sigma_row)

    def test_column_convention_coordinate_identity(self) -> None:
        p, n = 7, 2
        F = _symplectic_scale(3 * np.eye(n, dtype=np.int64), p)
        Sigma, B, info = minimal_qudit_frame(F, p)
        assert np.array_equal(mod_p(F @ B, p), mod_p(B @ Sigma, p))

    def test_certified_failure_uses_certification_error(self) -> None:
        rng = np.random.default_rng(13)
        F = rand_symplectic(rng, 3, 3, steps=8)
        try:
            minimal_qudit_frame(F, 3)
        except Exception as exc:
            assert isinstance(exc, MinimalQuditFrameCertificationError)
            assert isinstance(exc.info, dict)

    def test_block_decomposition_wrappers_preserve_certificate_info(self) -> None:
        p, n = 5, 2
        F = _symplectic_scale(2 * np.eye(n, dtype=np.int64), p)

        Sigma, B, info = block_decompose(F, p)
        verify_global_basis(F, B, Sigma, p)
        assert "Q_att" in info
        assert "Q_opt" in info
        _assert_cost_semantics(info)

        Sigma_c, B_c, info_c = block_decompose_certified(F, p, require_minimal=True)
        verify_global_basis(F, B_c, Sigma_c, p)
        assert info_c["minimal_cost_certified"] is True
        assert info_c["Q_opt"] == info_c["Q_att"]

        Sigma_o, B_o, info_o = block_decompose_optimal(F, p)
        verify_global_basis(F, B_o, Sigma_o, p)
        assert info_o["minimal_cost_certified"] is True
        assert info_o["Q_opt"] == info_o["Q_att"]

    def test_invalid_convention_raises_value_error(self) -> None:
        F = np.eye(2, dtype=np.int64)
        with pytest.raises(ValueError):
            minimal_qudit_frame(F, 2, allow_degraded=True, convention="sideways")
