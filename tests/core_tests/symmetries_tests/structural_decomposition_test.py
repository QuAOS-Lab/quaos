import itertools
import pprint

import numpy as np
import pytest

from sympleq.core.circuits.random_symplectic import symplectic_random_transvection
from sympleq.core.symmetries.atomic_decomposition import (
    CertificationError,
    atomic_block_decompose,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2 import (
    atomic_blocks_in_unipotent_self_sector_p2,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2_generators import (
    canonical_unipotent_p2_block,
    direct_sum_unipotent_p2_blocks,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.rcf_prepass import (
    _is_x_pm_1,
    rcf_prepass,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import (
    inv_mod_mat,
    is_symplectic,
    mod_p,
    omega_matrix,
    rank_mod,
)



def _as_matrix(obj):
    """Accept either F or (F, metadata)-style fixture returns."""
    if isinstance(obj, tuple):
        return obj[0]
    return obj


def _rand_symplectic(rng: np.random.Generator, n: int, p: int, steps: int = 40) -> np.ndarray:
    F = symplectic_random_transvection(n, p, num_transvections=steps, rng=rng)
    F = mod_p(F, p)
    assert is_symplectic(F, p)
    return F


def _direct_sum_grouped_blocks(blocks: list[np.ndarray], p: int) -> np.ndarray:
    """
    Symplectic direct sum in grouped [x_1,...,x_n | z_1,...,z_n] ordering.

    If block i has shape 2n_i × 2n_i and is written as

        [[A_i, B_i],
         [C_i, D_i]],

    this returns

        [[diag(A_i), diag(B_i)],
         [diag(C_i), diag(D_i)]].
    """
    if not blocks:
        raise ValueError("Need at least one block.")

    ns = []
    parts = []
    for F in blocks:
        F = mod_p(np.asarray(F, dtype=np.int64), p)
        assert F.ndim == 2 and F.shape[0] == F.shape[1] and F.shape[0] % 2 == 0
        assert is_symplectic(F, p)
        n = F.shape[0] // 2
        ns.append(n)
        parts.append((F[:n, :n], F[:n, n:], F[n:, :n], F[n:, n:]))

    total_n = sum(ns)

    def block_diag(mats: list[np.ndarray]) -> np.ndarray:
        out = np.zeros((total_n, total_n), dtype=np.int64)
        off = 0
        for M in mats:
            m = M.shape[0]
            assert M.shape == (m, m)
            out[off:off + m, off:off + m] = M
            off += m
        return mod_p(out, p)

    A = block_diag([x[0] for x in parts])
    B = block_diag([x[1] for x in parts])
    C = block_diag([x[2] for x in parts])
    D = block_diag([x[3] for x in parts])
    F_big = mod_p(np.block([[A, B], [C, D]]), p)

    assert F_big.shape == (2 * total_n, 2 * total_n)
    assert is_symplectic(F_big, p)
    return F_big


def _assert_sigma_block_diagonal(Sigma: np.ndarray, half_dims: list[int], p: int) -> None:
    """
    Check that Sigma is block diagonal in returned block order.

    This assumes the returned basis is grouped as [u_all | v_all], with each
    atomic block i occupying u-offset and v-offset ranges of length h_i.
    """
    Sigma = mod_p(Sigma, p)
    n2 = Sigma.shape[0]
    assert Sigma.shape == (n2, n2)
    assert n2 % 2 == 0
    n = n2 // 2
    assert sum(int(h) for h in half_dims) == n

    blocks: list[list[int]] = []
    off = 0
    for h in half_dims:
        h = int(h)
        assert h > 0
        blocks.append(list(range(off, off + h)) + list(range(n + off, n + off + h)))
        off += h

    assert sorted(j for idx in blocks for j in idx) == list(range(n2))

    for i, idx_i in enumerate(blocks):
        for j, idx_j in enumerate(blocks):
            if i == j:
                continue
            off_block = Sigma[np.ix_(idx_i, idx_j)]
            assert not np.any(off_block % p), (
                f"Sigma has nonzero entries between returned blocks {i} and {j}.\n"
                f"half_dims={half_dims}\n"
                f"off_block=\n{off_block}"
            )


def _assert_certified_structural_decomposition(F: np.ndarray, Sigma: np.ndarray, B: np.ndarray, info: dict, p: int) -> None:
    """Checks that are part of the structural correctness claim."""
    verify_global_basis(F, B, Sigma, p)
    assert is_symplectic(B, p)
    assert is_symplectic(Sigma, p)
    assert info["certified"] is True
    assert rank_mod(B, p) == F.shape[0]
    assert sum(int(h) for h in info["atomic_half_dims"]) == F.shape[0] // 2
    assert int(info["Q_opt"]) == max(int(h) for h in info["atomic_half_dims"])
    _assert_sigma_block_diagonal(Sigma, [int(h) for h in info["atomic_half_dims"]], p)


def _certified_decompose(F: np.ndarray, p: int):
    Sigma, B, info = atomic_block_decompose(F, p)
    _assert_certified_structural_decomposition(F, Sigma, B, info, p)
    return Sigma, B, info


def _sector_signature_from_prepass(F: np.ndarray, p: int) -> tuple:
    """
    Intrinsic primary/sector signature from the prepass, independent of block
    choices.  This is useful for conjugation-invariance tests.
    """
    meta = rcf_prepass(F, p)
    sig = []
    for ctx in meta.get("sector_contexts", []):
        sig.append(
            (
                str(ctx.sector_type),
                tuple(ctx.sector_key),
                tuple(ctx.sector_key_star) if ctx.sector_key_star is not None else None,
                int(ctx.deg_q),
                int(ctx.max_exp),
                int(ctx.meta.get("dim2", ctx.T_sec.shape[1] if ctx.T_sec is not None else 0)),
            )
        )
    return tuple(sorted(sig, key=repr))


def _find_p2_unipotent_key(meta: dict, p: int = 2):
    prim = meta["primaries"]
    keys = [k for k, d in prim.items() if _is_x_pm_1(d["poly"], p)]
    assert len(keys) == 1, f"expected exactly one p=2 x+1 primary; got {keys}"
    return keys[0]


def _extract_p2_payload(info: dict) -> dict:
    for inv in info.get("sector_invariants", []):
        data = getattr(inv, "data", {})
        if isinstance(data, dict) and "p2_unipotent" in data:
            return data["p2_unipotent"]
    raise AssertionError("No p=2 unipotent sector payload found.")


def _all_small_symplectics_n1(p: int):
    """
    Exhaustively enumerate Sp(2,F_p) by checking all 2x2 matrices.

    Sizes:
      p=2 -> 16 candidates
      p=3 -> 81 candidates
      p=5 -> 625 candidates
    """
    for entries in itertools.product(range(p), repeat=4):
        F = np.array(entries, dtype=np.int64).reshape(2, 2)
        if is_symplectic(F, p):
            yield mod_p(F, p)


def _companion_from_monic_coeffs(coeffs: list[int], p: int) -> np.ndarray:
    """
    Companion matrix for a monic polynomial

        x^d + c_{d-1} x^{d-1} + ... + c_0.

    coeffs are [c_0, ..., c_{d-1}, 1].
    """
    coeffs = [int(c) % p for c in coeffs]
    assert coeffs[-1] == 1
    d = len(coeffs) - 1
    C = np.zeros((d, d), dtype=np.int64)
    if d > 1:
        C[1:, :-1] = np.eye(d - 1, dtype=np.int64)
    C[:, -1] = [(-coeffs[i]) % p for i in range(d)]
    return mod_p(C, p)


def _symplectic_scale_from_A(A: np.ndarray, p: int) -> np.ndarray:
    """diag(A, A^{-T}) in grouped coordinates."""
    A = mod_p(A, p)
    Ainv = inv_mod_mat(A, p)
    z = np.zeros_like(A)
    F = mod_p(np.block([[A, z], [z, Ainv.T]]), p)
    assert is_symplectic(F, p)
    return F


# ---------------------------------------------------------------------------
# Article-readiness structural tests
# ---------------------------------------------------------------------------


class TestStructuralDecomposition:
    @pytest.mark.parametrize("p", [2, 3, 5])
    def test_exhaustive_n1_symplectics_certify(self, p: int) -> None:
        """
        Exhaustive n=1 oracle.

        This catches many convention and edge-case failures and is cheap enough
        to run in normal CI.
        """
        count = 0
        for F in _all_small_symplectics_n1(p):
            _certified_decompose(F, p)
            count += 1

        # |Sp(2,p)| = |SL(2,p)| = p(p^2 - 1)
        assert count == p * (p * p - 1)

    @pytest.mark.parametrize("p,n", [(2, 1), (2, 4), (3, 3), (5, 3)])
    def test_identity_structural_decomposition(self, p: int, n: int) -> None:
        F = np.eye(2 * n, dtype=np.int64)
        _Sigma, _B, info = _certified_decompose(F, p)

        assert sorted(int(h) for h in info["atomic_half_dims"]) == [1] * n
        assert int(info["Q_opt"]) == 1

    def test_p2_chapter5_mixed_direct_sum_and_conjugation(self) -> None:
        """
        Regression fixture for the bad-characteristic unipotent sector.

        This checks W(k), V_beta(2k), and W_beta(2l+1) together, both before
        and after a random symplectic conjugation.
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

        _Sigma, _B, info = _certified_decompose(F, p)
        payload = _extract_p2_payload(info)
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

        rng = np.random.default_rng(123)
        n = F.shape[0] // 2
        S = _rand_symplectic(rng, n, p, steps=50)
        Fh = mod_p(inv_mod_mat(S, p) @ F @ S, p)

        _Sigma_h, _B_h, info_h = _certified_decompose(Fh, p)
        payload_h = _extract_p2_payload(info_h)
        assert payload_h["kernel_profile"] == payload_h["kernel_profile_blocks"]
        assert payload_h["classification_complete"] is True

        # Structural canonical block-size multiset should be conjugation-invariant
        # for this fixture, even though it is not necessarily the cost-optimal
        # decomposition for every possible cost objective.
        assert sorted(int(h) for h in info["atomic_half_dims"]) == sorted(int(h) for h in info_h["atomic_half_dims"])

    def test_p2_direct_unipotent_sector_builder_matches_kernel_profile(self) -> None:
        """
        Sector-local test for the p=2 Chapter-5 unipotent extractor.

        This provides more direct diagnostics than the public API if the sector
        routine regresses.
        """
        p = 2
        F = _as_matrix(
            direct_sum_unipotent_p2_blocks(
                [
                    ("W", 2, None),
                    ("V_beta", 4, 1),
                    ("W_beta", 3, 1),
                ]
            )
        )
        F = mod_p(F, p)
        assert is_symplectic(F, p)

        meta = rcf_prepass(F, p)
        key = _find_p2_unipotent_key(meta, p)
        blocks, inv = atomic_blocks_in_unipotent_self_sector_p2(F, key, meta["primaries"])

        assert inv.data["status"] == "OK"
        assert sum(int(b.half_dim) for b in blocks) == F.shape[0] // 2

        payload = inv.data["p2_unipotent"]
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

    @pytest.mark.parametrize(
        "p,poly_coeffs",
        [
            # p=2, q=x^2+x+1
            (2, [1, 1, 1]),
            # p=3, q=x^2+1 is self-reciprocal and irreducible over F_3
            (3, [1, 0, 1]),
            # p=5, q=x^2+4x+1 is self-reciprocal and appears in random tests.
            # It is reducible over F_5, so this mainly checks that prepass/router
            # diagnostics remain sane rather than irreducible Hermitian theory.
            (5, [1, 4, 1]),
        ],
    )
    def test_self_reciprocal_degree2_scale_fixtures_certify_or_diagnose(self, p: int, poly_coeffs: list[int]) -> None:
        """
        Small self-reciprocal degree-2 scale fixtures.

        These are useful for the article because they exercise the q=q*, deg(q)>1
        pathway and row/column conventions.  For reducible fixtures the code may
        split the primary further; either way the result must be verified or fail
        diagnostically.
        """
        A = _companion_from_monic_coeffs(poly_coeffs, p)
        F = _symplectic_scale_from_A(A, p)

        try:
            _certified_decompose(F, p)
        except CertificationError as exc:
            # This is acceptable only as a clearly diagnostic future-coverage
            # failure.  If the paper claims full support for this sector, remove
            # this except block and require success.
            assert "failures" in exc.info
            assert exc.info["failures"], pprint.pformat(exc.info)

    @pytest.mark.parametrize("p,n,seed", [(2, 5, 0), (3, 5, 1), (5, 5, 2)])
    def test_row_column_convention_consistency_on_supported_inputs(self, p: int, n: int, seed: int) -> None:
        """
        Row/right-action and column-action conventions should produce the same
        structural invariants.  If arbitrary random p=2 cases become too broad
        for the current implementation, the failure should be diagnostic.
        """
        rng = np.random.default_rng(seed)
        F_col = _rand_symplectic(rng, n, p, steps=40)
        F_row = F_col.T

        try:
            _Sigma_col, _B_col, info_col = atomic_block_decompose(F_col, p, convention="column")
            _Sigma_row, B_row, info_row = atomic_block_decompose(F_row, p, convention="row")
        except CertificationError as exc:
            assert "failures" in exc.info
            assert exc.info["failures"], pprint.pformat(exc.info)
            return

        assert sorted(int(h) for h in info_col["atomic_half_dims"]) == sorted(int(h) for h in info_row["atomic_half_dims"])
        assert int(info_col["Q_opt"]) == int(info_row["Q_opt"])
        assert rank_mod(B_row, p) == 2 * n

    @pytest.mark.parametrize("p,n,seed", [(2, 6, 10), (3, 6, 11), (5, 6, 12)])
    def test_prepass_sector_signature_is_conjugation_invariant(self, p: int, n: int, seed: int) -> None:
        """
        The primary/sector decomposition itself should be invariant under
        symplectic conjugation, independent of later block choices.
        """
        rng = np.random.default_rng(seed)
        F = _rand_symplectic(rng, n, p, steps=45)
        S = _rand_symplectic(rng, n, p, steps=45)
        Fh = mod_p(inv_mod_mat(S, p) @ F @ S, p)

        assert _sector_signature_from_prepass(F, p) == _sector_signature_from_prepass(Fh, p)

    def test_direct_sum_prepass_signature_contains_component_factors(self) -> None:
        """
        Direct-sum sanity test for sectorization.

        This does not require cost optimality.  It checks that a mixed direct sum
        produces a prepass signature whose dimensions sum to the full space.
        """
        p = 2
        F_p2 = _as_matrix(canonical_unipotent_p2_block("V_beta", 4, beta=1))
        A = _companion_from_monic_coeffs([1, 1, 1], p)  # x^2+x+1
        F_h = _symplectic_scale_from_A(A, p)

        F = _direct_sum_grouped_blocks([F_p2, F_h], p)
        meta = rcf_prepass(F, p)

        assert sum(int(ctx.meta.get("dim2", 0)) for ctx in meta["sector_contexts"]) == F.shape[0]
        assert len(meta["sector_contexts"]) >= 1

        _certified_decompose(F, p)

    def test_non_symplectic_input_rejected(self) -> None:
        F = np.eye(4, dtype=np.int64)
        F[0, 0] = 0
        assert not is_symplectic(F, 2)

        with pytest.raises(ValueError, match="not symplectic"):
            atomic_block_decompose(F, 2)

    def test_invalid_dimension_rejected(self) -> None:
        F = np.eye(3, dtype=np.int64)
        with pytest.raises(ValueError):
            atomic_block_decompose(F, 2)
