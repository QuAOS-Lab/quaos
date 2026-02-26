from __future__ import annotations

import numpy as np
import galois

# gf2 inverse may live either in this package (older layout) or in the shared solvers
try:  # pragma: no cover
    from sympleq.core.finite_field_solvers import gf2_inv  # type: ignore
except Exception:  # pragma: no cover
    from .graph_automorphism_gf2 import gf2_inv  # type: ignore


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
