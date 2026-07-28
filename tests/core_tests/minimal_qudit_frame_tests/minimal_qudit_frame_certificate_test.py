"""
Additional decomposition-certificate tests for the MQF qudit decomposition.

These tests are intended to complement the targeted p=2 tests.  They focus on:

  - conjugation invariance of block-size/cost data,
  - mixed p=2 unipotent direct sums and random conjugation,
  - row/column convention consistency,
  - explicit Sigma block-diagonality according to returned mqf_half_dims,
  - cost-certificate consistency,
  - identity edge cases,
  - a small set of deterministic regression seeds.

Run with the focused minimal-qudit-frame test suite.
"""

from __future__ import annotations

import numpy as np
import pytest

from sympleq.core.circuits.random_symplectic import symplectic_random_transvection
from sympleq.core.minimal_qudit_frame import minimal_qudit_frame
from sympleq.core.minimal_qudit_frame.helpers.mqf_unipotent_p2 import (
    mqf_blocks_in_unipotent_self_sector_p2,
)
from sympleq.core.minimal_qudit_frame.helpers.mqf_unipotent_p2_generators import (
    direct_sum_unipotent_p2_blocks,
)
from sympleq.core.minimal_qudit_frame.helpers.rcf_prepass import (
    _is_x_pm_1,
    rcf_prepass,
)
from sympleq.core.minimal_qudit_frame.helpers.mqf_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import (
    independent_columns,
    inv_mod_mat,
    is_symplectic,
    mod_p,
    omega_matrix,
    rank_mod,
)


# ---------------------------------------------------------------------------
# Local helpers
# ---------------------------------------------------------------------------


def _as_matrix(obj):
    """Accept either F or (F, metadata)-style fixture returns."""
    if isinstance(obj, tuple):
        return obj[0]
    return obj


def _rand_symplectic(rng: np.random.Generator, n: int, p: int, steps: int) -> np.ndarray:
    """Use the project's transvection sampler with a short local alias."""
    F = symplectic_random_transvection(n, p, num_transvections=steps, rng=rng)
    F = mod_p(F, p)
    assert is_symplectic(F, p)
    return F


def _extract_p2_payload(info: dict) -> dict:
    for inv in info.get("sector_invariants", []):
        data = getattr(inv, "data", {})
        if isinstance(data, dict) and "p2_unipotent" in data:
            return data["p2_unipotent"]
    raise AssertionError("No p=2 unipotent payload found in info['sector_invariants'].")


def _find_p2_unipotent_key(meta: dict, p: int = 2):
    prim = meta["primaries"]
    keys = [k for k, d in prim.items() if _is_x_pm_1(d["poly"], p)]
    assert len(keys) == 1, f"expected exactly one p=2 x+1 primary; got {keys}"
    return keys[0]


def _assert_sigma_block_diagonal(Sigma: np.ndarray, half_dims: list[int], p: int) -> None:
    """
    Check that Sigma is block diagonal in the returned MQF block ordering.

    The global basis ordering is assumed to be grouped as [U_all | V_all], with
    per-block half dimensions h_i.  For block i, the coordinates are

        u_offset : u_offset+h_i
        n+u_offset : n+u_offset+h_i.
    """
    Sigma = mod_p(Sigma, p)
    assert Sigma.ndim == 2 and Sigma.shape[0] == Sigma.shape[1]
    n2 = Sigma.shape[0]
    assert n2 % 2 == 0
    n = n2 // 2
    assert sum(int(h) for h in half_dims) == n

    block_indices: list[list[int]] = []
    off = 0
    for h in half_dims:
        h = int(h)
        assert h > 0
        idx = list(range(off, off + h)) + list(range(n + off, n + off + h))
        block_indices.append(idx)
        off += h

    used = sorted(j for idx in block_indices for j in idx)
    assert used == list(range(n2))

    for i, idx_i in enumerate(block_indices):
        for j, idx_j in enumerate(block_indices):
            if i == j:
                continue
            off_block = Sigma[np.ix_(idx_i, idx_j)]
            assert not np.any(off_block % p), (
                f"Sigma has nonzero off-block entries between MQF blocks {i} and {j}.\n"
                f"half_dims={half_dims}\n"
                f"off_block=\n{off_block}"
            )


def _attained_cost(info: dict) -> int:
    return int(info.get("Q_att", info.get("attained_qudit_cost", info["qudit_cost"])))


def _assert_cost_semantics(info: dict) -> None:
    assert _attained_cost(info) == int(info["qudit_cost"])
    assert _attained_cost(info) == max(int(h) for h in info["mqf_half_dims"])
    if info.get("minimal_cost_certified", False):
        assert info.get("Q_opt") == _attained_cost(info)
        assert info.get("optimal_qudit_cost") == _attained_cost(info)
    else:
        assert info.get("Q_opt") is None
        assert info.get("optimal_qudit_cost") is None


def _assert_basic_certified_decomposition(F: np.ndarray, Sigma: np.ndarray, B: np.ndarray, info: dict, p: int) -> None:
    """Common global decomposition checks."""
    verify_global_basis(F, B, Sigma, p)
    assert is_symplectic(B, p)
    assert is_symplectic(Sigma, p)
    assert info["certified"] is True
    _assert_cost_semantics(info)
    assert sum(int(h) for h in info["mqf_half_dims"]) == F.shape[0] // 2
    _assert_sigma_block_diagonal(Sigma, [int(h) for h in info["mqf_half_dims"]], p)


class TestMinimalQuditFrameDecompositionCertificateCoverage:
    @pytest.mark.parametrize("p,n,seed", [(2, 5, 0), (3, 5, 1), (5, 5, 2)])
    def test_minimal_qudit_frame_conjugation_invariance(self, p: int, n: int, seed: int) -> None:
        """
        Conjugating F by a symplectic S should preserve the MQF block-size
        multiset and certified cost.
        """
        rng = np.random.default_rng(seed)

        F = _rand_symplectic(rng, n, p, steps=40)
        S = _rand_symplectic(rng, n, p, steps=40)
        F2 = mod_p(inv_mod_mat(S, p) @ F @ S, p)
        assert is_symplectic(F2, p)

        Sigma1, B1, info1 = minimal_qudit_frame(F, p)
        Sigma2, B2, info2 = minimal_qudit_frame(F2, p)

        _assert_basic_certified_decomposition(F, Sigma1, B1, info1, p)
        _assert_basic_certified_decomposition(F2, Sigma2, B2, info2, p)

        assert sorted(int(x) for x in info1["mqf_half_dims"]) == sorted(int(x) for x in info2["mqf_half_dims"])
        assert _attained_cost(info1) == _attained_cost(info2)
        assert bool(info1["minimal_cost_certified"]) == bool(info2["minimal_cost_certified"])

    def test_p2_mixed_chapter5_direct_sum_conjugated_certifies(self) -> None:
        """
        Mix W(k), V_beta(2k), and W_beta(2l+1), then hide the direct sum by
        random symplectic conjugation.  The decomposition should recover the
        same MQF half-dimension multiset.
        """
        p = 2
        F = _as_matrix(
            direct_sum_unipotent_p2_blocks(
                [
                    ("W", 3, None),
                    ("V_beta", 4, 1),
                    ("W_beta", 5, 1),
                    ("V_beta", 2, 0),
                ]
            )
        )
        F = mod_p(F, p)
        assert is_symplectic(F, p)

        n = F.shape[0] // 2
        rng = np.random.default_rng(123)
        S = _rand_symplectic(rng, n, p, steps=50)
        Fh = mod_p(inv_mod_mat(S, p) @ F @ S, p)
        assert is_symplectic(Fh, p)

        Sigma, B, info = minimal_qudit_frame(Fh, p)
        _assert_basic_certified_decomposition(Fh, Sigma, B, info, p)

        assert info["minimal_cost_certified"] is True
        assert info["certified_minimal_qudit_cost"] is True

        # Expected half dimensions:
        #   W(3)        -> 3
        #   V_beta(4)  -> 2
        #   W_beta(5)  -> 5
        #   V_beta(2)  -> 1
        assert sorted(int(x) for x in info["mqf_half_dims"]) == sorted([3, 2, 5, 1])
        assert info["minimal_cost_certified"] is True
        assert int(info["Q_opt"]) == 5

        payload = _extract_p2_payload(info)
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

    @pytest.mark.parametrize("p,n,seed", [(2, 5, 0), (3, 5, 1), (5, 5, 2)])
    def test_row_and_column_conventions_are_consistent(self, p: int, n: int, seed: int) -> None:
        """
        The row-action API should agree with the column-action API on block-size
        and cost invariants after transposing the input matrix.
        """
        rng = np.random.default_rng(seed)
        F_col = _rand_symplectic(rng, n, p, steps=35)
        F_row = F_col.T

        Sigma_col, B_col, info_col = minimal_qudit_frame(
            F_col, p, convention="column"
        )
        Sigma_row, B_row, info_row = minimal_qudit_frame(
            F_row, p, convention="row"
        )

        _assert_basic_certified_decomposition(F_col, Sigma_col, B_col, info_col, p)

        # Row convention returns row-convention objects.  The most robust API-level
        # check here is invariant equality rather than reusing the column verifier.
        assert sorted(int(x) for x in info_col["mqf_half_dims"]) == sorted(int(x) for x in info_row["mqf_half_dims"])
        assert _attained_cost(info_col) == _attained_cost(info_row)
        assert bool(info_col["minimal_cost_certified"]) == bool(info_row["minimal_cost_certified"])
        assert rank_mod(B_row, p) == 2 * n
        assert Sigma_row.shape == (2 * n, 2 * n)

    def test_p2_kernel_profile_is_conjugation_invariant_for_mixed_blocks(self) -> None:
        p = 2
        F = _as_matrix(
            direct_sum_unipotent_p2_blocks(
                [
                    ("W", 2, None),
                    ("W_beta", 3, 1),
                    ("V_beta", 4, 1),
                ]
            )
        )
        F = mod_p(F, p)
        assert is_symplectic(F, p)

        n = F.shape[0] // 2
        rng = np.random.default_rng(77)
        S = _rand_symplectic(rng, n, p, steps=50)
        Fh = mod_p(inv_mod_mat(S, p) @ F @ S, p)
        assert is_symplectic(Fh, p)

        Sigma, B, info = minimal_qudit_frame(Fh, p)
        _assert_basic_certified_decomposition(Fh, Sigma, B, info, p)

        payload = _extract_p2_payload(info)
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

    def test_p2_direct_sector_builder_certifies_mixed_unipotent_sum(self) -> None:
        """
        Sector-local version of the mixed p=2 unipotent test.  This gives a
        more direct failure if the p2 sector extractor regresses.
        """
        p = 2
        F = _as_matrix(
            direct_sum_unipotent_p2_blocks(
                [
                    ("W", 3, None),
                    ("V_beta", 4, 1),
                    ("W_beta", 5, 1),
                ]
            )
        )
        F = mod_p(F, p)
        assert is_symplectic(F, p)

        meta = rcf_prepass(F, p)
        key = _find_p2_unipotent_key(meta, p)

        blocks, inv = mqf_blocks_in_unipotent_self_sector_p2(F, key, meta["primaries"])

        assert inv.data["status"] == "OK"
        assert sum(int(b.half_dim) for b in blocks) == F.shape[0] // 2

        payload = inv.data["p2_unipotent"]
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

    @pytest.mark.parametrize("p,n,seed", [(2, 5, 0), (3, 4, 1), (5, 4, 2)])
    def test_certified_minimal_cost_matches_block_list(self, p: int, n: int, seed: int) -> None:
        rng = np.random.default_rng(seed)
        F = _rand_symplectic(rng, n, p, steps=35)

        Sigma, B, info = minimal_qudit_frame(F, p)
        _assert_basic_certified_decomposition(F, Sigma, B, info, p)

        assert _attained_cost(info) == max(int(h) for h in info["mqf_half_dims"])

        if info["minimal_cost_certified"]:
            cert = info["cost_certificate"]
            assert int(cert["attained"]) == _attained_cost(info)
            assert int(cert.get("qudit_cost", _attained_cost(info))) == _attained_cost(info)
            assert int(cert["lower_bound"]) == _attained_cost(info)
            assert cert["certified_minimal"] is True

    @pytest.mark.parametrize("p,n", [(2, 1), (2, 4), (3, 3), (5, 3)])
    def test_identity_decomposes_into_single_qudit_blocks(self, p: int, n: int) -> None:
        F = np.eye(2 * n, dtype=np.int64)

        Sigma, B, info = minimal_qudit_frame(F, p)
        _assert_basic_certified_decomposition(F, Sigma, B, info, p)

        assert sorted(int(x) for x in info["mqf_half_dims"]) == [1] * n
        assert info["minimal_cost_certified"] is True
        assert int(info["Q_opt"]) == 1
        assert info["minimal_cost_certified"] is True

    @pytest.mark.parametrize(
        "p,n,steps,seed",
        [
            (2, 5, 60, 0),
            (2, 6, 70, 1),
            (3, 5, 50, 2),
            (5, 5, 50, 3),
        ],
    )
    def test_known_regression_random_symplectics(self, p: int, n: int, steps: int, seed: int) -> None:
        rng = np.random.default_rng(seed)
        F = _rand_symplectic(rng, n, p, steps=steps)

        Sigma, B, info = minimal_qudit_frame(F, p)
        _assert_basic_certified_decomposition(F, Sigma, B, info, p)
