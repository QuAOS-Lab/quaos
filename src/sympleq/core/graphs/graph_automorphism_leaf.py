from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import galois

from sympleq.core.circuits.target import find_map_to_target_pauli_sum, get_phase_vector
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Gate
from sympleq.core.symmetries.phase_correction import solve_phase_vector_h_from_residual

from .graph_automorphism_code import check_code_automorphism


@dataclass
class LeafContext:
    p: int
    two_lcm: int
    n_qudits: int
    identity_perm: np.ndarray
    S_mod: np.ndarray
    G: galois.FieldArray
    G_mod2: np.ndarray | None
    basis_order: list[int]
    labels: list[int]
    pauli_sum: PauliSum
    ref_tableau: np.ndarray
    ref_phases: np.ndarray
    ref_weights: np.ndarray
    base_tableau: np.ndarray
    base_weights: np.ndarray
    base_phases: np.ndarray
    basis_indices: np.ndarray
    basis_source_ps: PauliSum
    # Precomputed inverse of the square row-basis tableau (when rank == 2*n_qudits).
    basis_src_inv_gf2: np.ndarray | None
    basis_src_inv_gfp: galois.FieldArray | None
    row_basis_cache: dict[str, np.ndarray]


def check_leaf(pi: np.ndarray, ctx: LeafContext) -> Gate | None:
    """Run all structural and phase-correction checks for a candidate permutation.

    Returns a symmetry Gate if successful, otherwise None.
    """
    if np.array_equal(pi, ctx.identity_perm):
        return None

    # (1) Edge-colour constraint: Pi must preserve the commutation matrix.
    if not np.array_equal(ctx.S_mod[np.ix_(pi, pi)], ctx.S_mod):
        return None

    # (2) Linear code (matroid) constraint.
    if not check_code_automorphism(ctx.G, ctx.basis_order, ctx.labels, pi, ctx.G_mod2):
        return None

    # (3) Build the candidate symplectic from the permutation of the (ordered) row-basis.
    tgt_idx = pi[ctx.basis_indices]
    H_basis_tgt = PauliSum.from_tableau(
        ctx.base_tableau[tgt_idx],
        ctx.pauli_sum.dimensions,
        weights=ctx.base_weights[tgt_idx],
    )
    H_basis_tgt.set_phases(np.array(ctx.base_phases[tgt_idx], dtype=int, copy=True))
    H_basis_src = ctx.basis_source_ps

    F: np.ndarray
    h0: np.ndarray

    if ctx.p == 2 and ctx.basis_src_inv_gf2 is not None:
        T = (ctx.base_tableau[tgt_idx] & 1).astype(np.uint8, copy=False)
        F = (ctx.basis_src_inv_gf2 @ T) & 1
        F = np.asarray(F, dtype=int)
    elif ctx.p != 2 and ctx.basis_src_inv_gfp is not None:
        GF = galois.GF(int(ctx.p))
        T = GF(ctx.base_tableau[tgt_idx] % ctx.p)
        F_gf = ctx.basis_src_inv_gfp @ T
        F = (np.asarray(F_gf, dtype=int) % ctx.p).astype(int, copy=False)
    else:
        # Rank-deficient tableau: the basis alone doesn't determine F uniquely.
        F, _h0, _, _ = find_map_to_target_pauli_sum(H_basis_src, H_basis_tgt)

    # Deterministic lift from symplectic -> quadratic phase vector; linear correction solved below.
    h0 = get_phase_vector(F.T, int(ctx.pauli_sum.dimensions[0]))

    nq = ctx.n_qudits
    pauli = ctx.pauli_sum

    SG_F = Gate('Symmetry', F.T, np.asarray(h0, dtype=int))
    H_full_tg = pauli.copy()[pi]
    H_full_F = SG_F.act(pauli, tuple(range(nq)))

    delta = (H_full_tg.phases - H_full_F.phases) % (2 * int(ctx.pauli_sum.lcm))

    # For qubits, the quadratic part is not unique; if odd residuals appear, try a standard diagonal lift.
    if ctx.p == 2 and np.any(delta % 2 != 0):
        F2 = F % 2
        A, B = F2[:nq, :nq], F2[:nq, nq:]
        C, D = F2[nq:, :nq], F2[nq:, nq:]
        hx0 = np.diag((A @ B.T) % 2) % 2
        hz0 = np.diag((C @ D.T) % 2) % 2
        h0_alt = np.concatenate([hx0, hz0]).astype(int)

        SG_F_alt = Gate('Symmetry', F.T, h0_alt)
        H_full_Fa = SG_F_alt.act(pauli, tuple(range(nq)))
        delta_alt = (H_full_tg.phases - H_full_Fa.phases) % 4
        if (delta_alt % 2).sum() < (delta % 2).sum():
            h0, SG_F, H_full_F, delta = h0_alt, SG_F_alt, H_full_Fa, delta_alt

    h_lin = solve_phase_vector_h_from_residual(
        ctx.base_tableau,
        delta,
        ctx.pauli_sum.dimensions,
        debug=False,
        row_basis_cache=ctx.row_basis_cache,
    )
    if h_lin is None:
        return None

    h0_mod = np.asarray(h0, dtype=int) % ctx.two_lcm
    h_lin_mod = np.asarray(h_lin, dtype=int) % ctx.two_lcm
    h_tot = (h0_mod + h_lin_mod) % ctx.two_lcm
    SG = Gate('Symmetry', F.T, h_tot)

    # Final verification: act and compare canonical forms.
    H_out_cf = SG.act(pauli, tuple(range(nq))).to_standard_form()
    H_out_cf.weight_to_phase()

    if not np.array_equal(H_out_cf.tableau, ctx.ref_tableau):
        return None
    if not np.all((ctx.ref_phases - H_out_cf.phases) % ctx.two_lcm == 0):
        return None
    if not np.all(np.isclose(H_out_cf.weights, ctx.ref_weights, atol=1e-8, rtol=0)):
        return None

    return SG
