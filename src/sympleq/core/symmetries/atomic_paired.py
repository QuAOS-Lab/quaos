# sympleq/core/symmetries/atomic_paired.py
from __future__ import annotations
import numpy as np
from typing import List, Tuple

from .modular_helpers import mod_p, independent_columns, rank_mod, solve_linear_many, omega_matrix
from .atomic_types import AtomicBlock, AtomicInvariant
from .atomic_linear import restrict_operator
from .symplectic_basis import symplectic_gram_schmidt_from_span
from .module_invariants import restrict_operator_invariant, q_of_F_restricted, jordan_chain_tops_nilpotent, cyclic_submodule_basis


def _pairing_matrix_between(
    V_left: np.ndarray, V_right: np.ndarray, p: int
) -> np.ndarray:
    """
    Compute bilinear pairing matrix P = V_left^T Ω V_right (mod p),
    where columns are basis vectors.
    """
    n2 = V_left.shape[0]
    Ω = omega_matrix(n2 // 2, p)
    return mod_p(V_left.T @ Ω @ V_right, p)


def atomic_blocks_in_paired_sector(
    F: np.ndarray,
    p: int,
    key: Tuple[int, ...],
    key_star: Tuple[int, ...],
    primaries: dict,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Deterministic atomic block construction for paired sector W = V_q ⊕ V_{q*}.
    Output blocks are hyperbolic bases in ambient coordinates.
    """
    q = primaries[key]["poly"]
    qstar = primaries[key_star]["poly"]
    Vq = independent_columns(primaries[key]["V_basis"], p)
    Vqs = independent_columns(primaries[key_star]["V_basis"], p)

    # Restrict F to each invariant subspace coordinates
    Fq = restrict_operator_invariant(F, Vq, p)
    Fqs = restrict_operator_invariant(F, Vqs, p)

    # Nilpotent parts: N = q(F) on V_q, N* = q*(F) on V_{q*}
    Nq = q_of_F_restricted(Fq, q, p)
    Nqs = q_of_F_restricted(Fqs, qstar, p)

    # Extract chain tops by length in each side
    # (These are "atomic seeds": each top corresponds to a cyclic submodule.)
    max_exp = int(primaries[key]["exponent"])
    tops_left = jordan_chain_tops_nilpotent(Nq, max_exp, p)   # dict L -> (d x t_L)
    tops_right = jordan_chain_tops_nilpotent(Nqs, max_exp, p)

    # In a well-formed paired sector, multiplicities should match by length.
    lengths = sorted(set(tops_left.keys()) | set(tops_right.keys()))
    inv_data = {"length_multiplicities": {}}
    blocks: List[AtomicBlock] = []

    # Precompute ambient pairing map between Vq and Vqs
    P = _pairing_matrix_between(Vq, Vqs, p)  # shape (dim(Vq), dim(Vqs))

    for L in lengths:
        A = tops_left.get(L, np.zeros((Fq.shape[0], 0), dtype=np.int64))
        B = tops_right.get(L, np.zeros((Fqs.shape[0], 0), dtype=np.int64))
        a_mult = A.shape[1]
        b_mult = B.shape[1]
        inv_data["length_multiplicities"][int(L)] = (int(a_mult), int(b_mult))

        if a_mult != b_mult:
            raise RuntimeError(f"Paired sector mismatch at length {L}: {a_mult} vs {b_mult}")

        # For each top vector on left, build its cyclic basis and find matching partner on right.
        for j in range(a_mult):
            v_top = A[:, j:j+1]  # coordinates in Vq basis
            w_top = B[:, j:j+1]  # coordinates in Vqs basis
            # Build cyclic submodule bases (in sector coordinates)
            C_left = cyclic_submodule_basis(Fq, Nq, v_top, primaries[key]["deg"], L, p)    # (dimVq x deg*L)
            C_right = cyclic_submodule_basis(Fqs, Nqs, w_top, primaries[key_star]["deg"], L, p)

            # Lift to ambient
            W_left = mod_p(Vq @ C_left, p)     # (2n x deg*L)
            W_right = mod_p(Vqs @ C_right, p)

            # Span of the atomic block is W_left ⊕ W_right
            span = independent_columns(np.concatenate([W_left, W_right], axis=1), p)

            # Build canonical hyperbolic basis for this atomic span
            T_blk = symplectic_gram_schmidt_from_span(span, p)  # (2n x 2k)

            half_dim = T_blk.shape[1] // 2
            blocks.append(
                AtomicBlock(
                    T_blk=T_blk,
                    half_dim=half_dim,
                    sector_key=key,
                    inv=None
                )
            )

    inv = AtomicInvariant(sector_key=key, sector_type="paired", poly_key=key, data=inv_data)
    # Attach same invariant to each block for traceability
    blocks = [AtomicBlock(b.T_blk, b.half_dim, b.sector_key, inv) for b in blocks]
    return blocks, inv
