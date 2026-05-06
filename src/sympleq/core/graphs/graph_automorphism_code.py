from __future__ import annotations

import numpy as np
import galois

from sympleq.core.finite_field_solvers import gf2_inv


def check_code_automorphism(
    G: galois.FieldArray,
    basis_order: list[int],
    labels: list[int],
    pi: np.ndarray,
    G_mod2: np.ndarray | None = None,
) -> bool:
    """Linear-code test: does there exist U with U G P = G?

    This is the standard generator-matrix automorphism test for a linear code.
    Let B be the chosen basis columns (by label), and C = G[:, P(B)]. If C is
    invertible, the unique U that would satisfy the equation is U = C^{-1}.
    We then check whether U (G P) == G.

    Notes
    -----
    - Correctness is independent of the particular basis_order, but using a
      fixed basis_order makes the test deterministic.
    - The GF(2) path uses a fast XOR-based inverse.
    """
    lab_to_idx = {lab: i for i, lab in enumerate(labels)}
    B_cols = np.array([lab_to_idx[b] for b in basis_order], dtype=int)
    PBcols = pi[B_cols]

    if G_mod2 is not None:
        C = G_mod2[:, PBcols]
        try:
            C_inv = gf2_inv(C)
        except np.linalg.LinAlgError:
            return False
        Gp = G_mod2[:, pi]
        return np.array_equal((C_inv @ Gp) & 1, G_mod2)

    C = G[:, PBcols]
    try:
        U = np.linalg.inv(C)  # galois overrides numpy.linalg for FieldArray
    except np.linalg.LinAlgError:
        return False
    Gp = G[:, pi]
    return np.array_equal(U @ Gp, G)


def compute_induced_completion_matrix_gf2(
    G_mod2: np.ndarray,
    B_cols: np.ndarray,
    pi: np.ndarray,
) -> np.ndarray | None:
    """Return C = G[:, pi(B)] if invertible, else None.

    For GF(2) code automorphisms we require existence of U with U G P = G.
    Fixing the images of the basis columns B determines
        C := G[:, P(B)]   (k x k)
        U := C^{-1}.

    Once C is known and invertible, the automorphism condition implies
        g_{P(i)} = C g_i  for all i.

    This helper validates invertibility (by attempting gf2_inv) and returns C
    for use in induced-completion pruning.
    """
    PB = pi[B_cols]
    C = G_mod2[:, PB]
    try:
        _ = gf2_inv(C)
    except np.linalg.LinAlgError:
        return None
    return C


def compute_prefix_UG(
    G: galois.FieldArray,
    G_mod2: np.ndarray | None,
    B_cols: np.ndarray,
    pi: np.ndarray,
) -> np.ndarray | galois.FieldArray | None:
    """Return U G for a partial permutation once all basis columns are mapped.

    The current search stores ``pi[i] = image(i)``. Once every basis column has
    an image, the code automorphism condition determines U from
    ``U G[:, pi(B)] = G[:, B]``. The search then checks columns incrementally via
    ``(U G)[:, pi(i)] == G[:, i]``.
    """
    PB = np.asarray(pi, dtype=int)[np.asarray(B_cols, dtype=int)]
    if np.any(PB < 0):
        return None

    if G_mod2 is not None:
        C = G_mod2[:, PB]
        try:
            C_inv = gf2_inv(C)
        except np.linalg.LinAlgError:
            return None
        return (C_inv @ G_mod2) & 1

    C = G[:, PB]
    try:
        U = np.linalg.inv(C)
    except np.linalg.LinAlgError:
        return None
    return U @ G
