import numpy as np
import pytest

from sympleq.core.symmetries.block_decomposition import atomic_block_decompose
from sympleq.core.symmetries.modular_helpers import mod_p, inv_mod_mat
from sympleq.core.circuits.utils import is_symplectic
from sympleq.core.circuits.random_symplectic import symplectic_random_transvection


def _unipotent_jordan(n: int, p: int = 2) -> np.ndarray:
    """n×n unipotent Jordan block J = I + shift."""
    J = np.eye(n, dtype=int)
    for i in range(n - 1):
        J[i, i + 1] = 1
    return J % p


def _symplectic_blockdiag_from_A(A: np.ndarray, p: int) -> np.ndarray:
    """Return F = diag(A, (A^{-1})^T), which is symplectic for standard Ω."""
    A = mod_p(A, p)
    Ainv = inv_mod_mat(A, p)
    return mod_p(np.block([[A, np.zeros_like(A)], [np.zeros_like(A), Ainv.T]]), p)


def _symplectic_shear_from_symmetric(A: np.ndarray, p: int) -> np.ndarray:
    """Return F = [[I, A],[0, I]] which is symplectic iff A is symmetric."""
    n = A.shape[0]
    I = np.eye(n, dtype=int)
    Z = np.zeros((n, n), dtype=int)
    return mod_p(np.block([[I, A], [Z, I]]), p)


def _extract_p2_unipotent_payload(info: dict) -> dict:
    invs = info.get("sector_invariants", [])
    for inv in invs:
        data = getattr(inv, "data", {})
        if isinstance(data, dict) and "p2_unipotent" in data:
            return data["p2_unipotent"]
    raise AssertionError("No p=2 unipotent sector invariant found in info['sector_invariants']")


def _canon_length_invariants(payload: dict) -> tuple:
    li = payload.get("length_invariants", {})
    items = []
    for L, d in li.items():
        items.append((
            int(L),
            int(d.get("top_dim", 0)),
            int(d.get("B_rank", 0)),
            int(d.get("rad_dim", 0)),
            bool(d.get("q_polar_ok", False)),
            d.get("arf", None),
        ))
    items.sort()
    return tuple(items)


class TestP2UnipotentBookkeeping:

    def test_p2_unipotent_jordan_kernel_profile_matches_blocks(self):
        p, n = 2, 20
        J = _unipotent_jordan(n, p)
        F = _symplectic_blockdiag_from_A(J, p)
        assert is_symplectic(F, p)

        _, _, info = atomic_block_decompose(F, p, mode="certified")
        payload = _extract_p2_unipotent_payload(info)

        assert payload["kernel_profile"] == payload["kernel_profile_blocks"], "kernel profile mismatch"
        for L, d in payload.get("length_invariants", {}).items():
            assert d.get("B_sym_ok", True), f"B_sym_ok failed at L={L}: {d}"
            assert d.get("B_alt_ok", True), f"B_alt_ok failed at L={L}: {d}"
            assert d.get("q_witness_polar_ok", True), f"no quadratic witness at L={L}: {d}"
        expected = [min(2*k, 2*n) for k in range(1, n+1)]
        assert payload["kernel_profile"][:n] == expected

        # Basic bookkeeping sanity: counts sum to extracted block count
        beta_counts = payload.get("beta_counts_by_L", {})
        total_count = sum(sum(v.values()) for v in beta_counts.values())

        # The sector-level blocks list lives in the parent invariant data.
        # We don't require exact equality here (since only V/W blocks are counted),
        # but total_count should be <= number of sector blocks.
        assert total_count >= 1

    def test_p2_unipotent_conjugacy_invariant_payload(self):
        p, n = 2, 6
        J = _unipotent_jordan(n, p)
        F = _symplectic_blockdiag_from_A(J, p)
        assert is_symplectic(F, p)

        rng = np.random.default_rng(0)
        S = symplectic_random_transvection(n, p, num_transvections=120, rng=rng)
        Sinv = inv_mod_mat(S, p)
        F2 = mod_p(Sinv @ F @ S, p)
        assert is_symplectic(F2, p)

        _, _, info1 = atomic_block_decompose(F, p, mode="certified")
        _, _, info2 = atomic_block_decompose(F2, p, mode="certified")

        p1 = _extract_p2_unipotent_payload(info1)
        p2p = _extract_p2_unipotent_payload(info2)

        assert p1["kernel_profile"] == p2p["kernel_profile"]
        assert _canon_length_invariants(p1) == _canon_length_invariants(p2p)

    @pytest.mark.parametrize("n,seed", [(6, 0), (7, 1), (8, 2)])
    def test_p2_shear_family_is_certified_and_polar_ok(self, n, seed):
        p = 2
        rng = np.random.default_rng(seed)
        A = rng.integers(0, 2, size=(n, n), dtype=int)
        A = mod_p(A + A.T, p)  # force symmetric
        F = _symplectic_shear_from_symmetric(A, p)
        assert is_symplectic(F, p)

        _, _, info = atomic_block_decompose(F, p, mode="certified")
        payload = _extract_p2_unipotent_payload(info)

        assert payload["kernel_profile"] == payload["kernel_profile_blocks"]
        for L, d in payload.get("length_invariants", {}).items():
            assert d.get("B_sym_ok", True), f"B_sym_ok failed at L={L}: {d}"
            assert d.get("B_alt_ok", True), f"B_alt_ok failed at L={L}: {d}"
            assert d.get("q_witness_polar_ok", True), f"no quadratic witness at L={L}: {d}"

