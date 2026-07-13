import numpy as np
import pytest

from sympleq.core.symmetries.atomic_decomposition import atomic_block_decompose
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2_generators import (
    canonical_unipotent_p2_block,
    canonical_W_block_p2,
    direct_sum_unipotent_p2_blocks,
    random_symplectic_conjugate_p2,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import inv_mod_mat, is_symplectic, mod_p


def _unipotent_jordan(n: int) -> np.ndarray:
    J = np.eye(n, dtype=np.int64)
    for i in range(n - 1):
        J[i, i + 1] = 1
    return J % 2

def _as_matrix(obj):
    """
    Accept either F or (F, metadata)-style generator returns.
    Keeps tests robust if the fixture helpers return auxiliary data.
    """
    if isinstance(obj, tuple):
        return obj[0]
    return obj

def _symplectic_blockdiag_from_A(A: np.ndarray) -> np.ndarray:
    
    A = mod_p(A, 2)
    Ainv = inv_mod_mat(A, 2)
    Z = np.zeros_like(A)
    return mod_p(np.block([[A, Z], [Z, Ainv.T]]), 2)


def _symplectic_shear_upper(A: np.ndarray) -> np.ndarray:
    n = A.shape[0]
    I = np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    return mod_p(np.block([[I, A], [Z, I]]), 2)


def _extract_p2_unipotent_payload(info: dict) -> dict:
    for inv in info.get("sector_invariants", []):
        data = getattr(inv, "data", {})
        if isinstance(data, dict) and "p2_unipotent" in data:
            return data["p2_unipotent"]
    raise AssertionError("No p=2 unipotent invariant found")


def _canon_length_invariants(payload: dict) -> tuple:
    out = []
    for L, d in payload.get("length_invariants", {}).items():
        out.append((
            int(L),
            int(d.get("top_dim", 0)),
            int(d.get("B_rank", 0)),
            int(d.get("rad_dim", 0)),
            bool(d.get("B_sym_ok", True)),
            bool(d.get("B_alt_ok", True)),
            bool(d.get("q_witness_polar_ok", True)),
            d.get("arf", None),
        ))
    return tuple(sorted(out))


class TestP2UnipotentBookkeeping:
    def test_p2_unipotent_Wn_fixture_is_certified(self) -> None:
        p, n = 2, 6
        F = _symplectic_blockdiag_from_A(_unipotent_jordan(n))
        assert is_symplectic(F, p)

        Sigma, B, info = atomic_block_decompose(F, p)

        verify_global_basis(F, B, Sigma, p)

        payload = _extract_p2_unipotent_payload(info)

        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

        assert info["certified"] is True
        assert info["minimal_cost_certified"] is True
        assert info["certified_minimal_qudit_cost"] is True

        # For F = diag(J_n, J_n^{-T}), the whole sector is one W(n)-type block.
        assert info["Q_opt"] == n
        assert sorted(info["atomic_half_dims"]) == [n]

    def test_p2_unipotent_conjugacy_invariant_payload_for_W_block(self) -> None:
        F = canonical_W_block_p2(5)
        F2, _S = random_symplectic_conjugate_p2(F, seed=4, steps=20)
        _, _, info1 = atomic_block_decompose(F, 2)
        _, _, info2 = atomic_block_decompose(F2, 2)
        p1 = _extract_p2_unipotent_payload(info1)
        p2 = _extract_p2_unipotent_payload(info2)
        assert p1["kernel_profile"] == p2["kernel_profile"]
        assert _canon_length_invariants(p1) == _canon_length_invariants(p2)

    @pytest.mark.parametrize("lengths", [[2], [3], [2, 4]])
    def test_direct_sum_W_block_fixtures_are_symplectic_and_decompose(self, lengths) -> None:
        F = direct_sum_unipotent_p2_blocks([("W", L, 0) for L in lengths])
        assert is_symplectic(F, 2)
        Sigma, B, info = atomic_block_decompose(F, 2)
        verify_global_basis(F, B, Sigma, 2)
        payload = _extract_p2_unipotent_payload(info)
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["length_dropping_used"] is False

    @pytest.mark.parametrize(
        "kind,length,beta",
        [
            ("V_beta", 2, 0),
            ("V_beta", 2, 1),
            ("V_beta", 4, 0),
            ("V_beta", 4, 1),
            ("W_beta", 3, 1),
            ("W_beta", 5, 1),
        ],
    )
    def test_chapter5_fixtures_are_implemented_and_certify(self, kind, length, beta) -> None:
        F = _as_matrix(canonical_unipotent_p2_block(kind, length, beta=beta))
        assert is_symplectic(F, 2)

        Sigma, B, info = atomic_block_decompose(F, 2)
        verify_global_basis(F, B, Sigma, 2)

        payload = _extract_p2_unipotent_payload(info)
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

        assert info["certified"] is True
        assert info["minimal_cost_certified"] is True
        assert info["certified_minimal_qudit_cost"] is True

    def test_chapter5_direct_sum_fixture_certifies(self) -> None:
        F = _as_matrix(
            direct_sum_unipotent_p2_blocks(
                [
                    ("V_beta", 2, 1),
                    ("W_beta", 3, 1),
                    ("V_beta", 4, 0),
                ]
            )
        )
        assert is_symplectic(F, 2)

        Sigma, B, info = atomic_block_decompose(F, 2)
        verify_global_basis(F, B, Sigma, 2)

        payload = _extract_p2_unipotent_payload(info)
        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

        assert info["certified"] is True
        assert info["minimal_cost_certified"] is True

    @pytest.mark.parametrize("n,seed", [(4, 0), (5, 1), (6, 2)])
    def test_p2_upper_shear_family_certifies_and_payload_is_consistent(self, n: int, seed: int) -> None:
        rng = np.random.default_rng(seed)

        # Symmetric over F_2. Random rank/diagonal data means the decomposition may
        # contain a mixture of length-2 and length-1 pieces, so do not assert exact
        # Chapter-5 labels here.
        A = rng.integers(0, 2, size=(n, n), dtype=np.int64)
        A = mod_p(A + A.T + np.diag(np.diag(A)), 2)

        F = _symplectic_shear_upper(A)
        assert is_symplectic(F, 2)

        Sigma, B, info = atomic_block_decompose(F, 2)
        verify_global_basis(F, B, Sigma, 2)

        payload = _extract_p2_unipotent_payload(info)

        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        assert payload["classification_complete"] is True

        assert info["certified"] is True
        assert info["minimal_cost_certified"] is True
        assert info["certified_minimal_qudit_cost"] is True

        # These are structural sanity checks on the diagnostic forms, not a full
        # assertion about every debug-only quadratic witness.
        for L, d in payload.get("length_invariants", {}).items():
                assert d.get("B_sym_ok", True), f"B_sym_ok failed at L={L}: {d}"
                assert int(d.get("top_dim", 0)) >= int(d.get("B_rank", 0))
                assert int(d.get("rad_dim", 0)) == int(d.get("top_dim", 0)) - int(d.get("B_rank", 0))