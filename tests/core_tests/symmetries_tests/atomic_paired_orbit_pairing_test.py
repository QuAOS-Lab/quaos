import numpy as np

from sympleq.core.symmetries.atomic_decomposition import atomic_block_decompose
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import inv_mod_mat, is_symplectic, mod_p


def _block_diag(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    return np.block([
        [A, np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)],
        [np.zeros((B.shape[0], A.shape[1]), dtype=np.int64), B],
    ])


def _symplectic_scale(A: np.ndarray, p: int) -> np.ndarray:
    A = mod_p(A, p)
    return mod_p(_block_diag(A, inv_mod_mat(A, p).T), p)


def _paired_progress(info: dict) -> list[dict]:
    out: list[dict] = []
    for inv in info.get("sector_invariants", []):
        if getattr(inv, "sector_type", None) == "paired":
            out.extend(getattr(inv, "data", {}).get("progress", []))
    return out


def test_linear_paired_sector_uses_orbit_pairing_matrix() -> None:
    p, n = 5, 4
    F = _symplectic_scale(2 * np.eye(n, dtype=np.int64), p)
    assert is_symplectic(F, p)

    Sigma, B, info = atomic_block_decompose(F, p)
    verify_global_basis(F, B, Sigma, p)

    assert info["certified"] is True
    assert info["certified_minimal"] is True
    assert sorted(int(h) for h in info["atomic_half_dims"]) == [1, 1, 1, 1]

    progress = _paired_progress(info)
    assert progress
    assert all(step.get("attempt") == "orbit_pairing" for step in progress)
    assert all(step.get("orbit_pairing_matrix_constructed") is True for step in progress)
    assert all(step.get("dual_solve_used") is True for step in progress)
    assert all(step.get("orbit_pairing_full_rank") is True for step in progress)


def test_nonlinear_degree_two_paired_sector_uses_expanded_orbit_matrix() -> None:
    p = 3
    # A has irreducible characteristic polynomial x^2 + x + 2 over F_3.
    # The symplectic scale diag(A, A^{-T}) therefore gives one paired sector
    # with deg(q)=2 and one length-1 head representative on each side.
    A = np.array([[0, 1], [1, 2]], dtype=np.int64)
    F = _symplectic_scale(A, p)
    assert is_symplectic(F, p)

    Sigma, B, info = atomic_block_decompose(F, p)
    verify_global_basis(F, B, Sigma, p)

    assert info["certified"] is True
    assert info["certified_minimal"] is True
    assert info["atomic_half_dims"] == [2]
    assert info["Q_opt"] == 2

    progress = _paired_progress(info)
    assert len(progress) == 1
    step = progress[0]
    assert step["attempt"] == "orbit_pairing"
    assert step["deg_q"] == 2
    assert step["orbit_pairing_shape"] == (2, 2)
    assert step["orbit_pairing_rank"] == 2
    assert step["orbit_pairing_full_rank"] is True
    assert step["dual_solve_used"] is True
