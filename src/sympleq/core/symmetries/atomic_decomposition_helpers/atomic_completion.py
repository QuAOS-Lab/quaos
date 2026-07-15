from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

from ..modular_helpers import inv_mod_mat, mod_p, omega_matrix
from .atomic_linear import split_uv, symplectic_completion_from_block
from .atomic_types import AtomicBlock


def concat_blocks_to_partial_basis(blocks: List[AtomicBlock], n2: int, p: int) -> np.ndarray:
    """
    Build a partial symplectic frame T = [U_all | V_all] (2n x 2k)
    from block bases. This does not require spanning the full space.
    """
    if not blocks:
        return np.zeros((n2, 0), dtype=np.int64)

    u_list: list[np.ndarray] = []
    v_list: list[np.ndarray] = []
    for block in blocks:
        u_blk, v_blk = split_uv(block.T_blk)
        u_list.append(u_blk)
        v_list.append(v_blk)

    return mod_p(np.concatenate(u_list + v_list, axis=1), p)


def verify_full_symplectic_basis(B: np.ndarray, p: int) -> None:
    """Raise if B is not a full symplectic basis."""
    if B.ndim != 2 or B.shape[0] != B.shape[1]:
        raise RuntimeError(f"Global basis must be square, got shape {B.shape}.")
    n2 = B.shape[0]
    if n2 % 2 != 0:
        raise RuntimeError(f"Global basis must have even dimension, got {n2}.")
    omega = omega_matrix(n2 // 2, p)
    gram = mod_p(B.T @ omega @ B, p)
    if not np.array_equal(gram % p, omega % p):
        raise RuntimeError("Global basis B is not symplectic (B^T Omega B != Omega).")


def compute_sigma(F: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """Sigma = B^{-1} F B (mod p)."""
    return mod_p(inv_mod_mat(B, p) @ F @ B, p)


def convert_column_result_to_row(
    Sigma_col: np.ndarray,
    B_col: np.ndarray,
    info: Dict[str, Any],
    p: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
    """Convert an internally column-action result back to row/right-action convention."""
    Sigma_row = mod_p(Sigma_col.T, p)
    B_row = mod_p(B_col.T, p)
    info = dict(info)
    info["input_convention"] = "row"
    info["internal_convention"] = "column"
    info["basis_convention"] = "row"
    info["convention_note"] = (
        "Internal computation used F_col = F_row.T. Returned B has row-basis convention, "
        "so Sigma = B F_row B^{-1} = Sigma_col.T."
    )
    return Sigma_row, B_row, info


def complete_global_basis_from_blocks(
    *,
    F: np.ndarray,
    p: int,
    blocks: List[AtomicBlock],
) -> Tuple[np.ndarray, bool]:
    """
    Return a full symplectic basis B and whether a completion step beyond
    concatenating block bases was needed.
    """
    n2 = F.shape[0]
    omega = omega_matrix(n2 // 2, p)

    partial = concat_blocks_to_partial_basis(blocks, n2, p)
    if partial.size == 0:
        basis = np.eye(n2, dtype=np.int64)
        verify_full_symplectic_basis(basis, p)
        return mod_p(basis, p), True

    k2 = partial.shape[1]
    if k2 % 2 != 0:
        raise RuntimeError(f"Partial basis has odd number of columns {k2}, cannot be a symplectic frame.")
    omega_k = omega_matrix(k2 // 2, p)
    gram = mod_p(partial.T @ omega @ partial, p)
    if not np.array_equal(gram % p, omega_k % p):
        raise RuntimeError("Partial basis T is not a Darboux frame: T^T Omega T != Omega_k. Upstream block bug?")

    if k2 == n2:
        verify_full_symplectic_basis(partial, p)
        return mod_p(partial, p), False

    basis = symplectic_completion_from_block(partial, p)
    verify_full_symplectic_basis(basis, p)
    return mod_p(basis, p), True
