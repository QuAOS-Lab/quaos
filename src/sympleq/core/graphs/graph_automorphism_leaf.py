from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import galois

from sympleq.core.circuits.target import get_phase_vector
from sympleq.core.circuits.target import find_map_to_target_pauli_sum
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Gate
from sympleq.core.circuits.phase_correction import solve_phase_vector_h_from_residual
from sympleq.core.finite_field_solvers import solve_gf2
from sympleq.core.symmetries.modular_helpers import (
    _solve_linear,
    mod_p,
    nullspace_mod,
    omega_matrix,
    rank_mod,
)

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
    # Full row-basis completion of basis_source_ps.tableau. Used when the
    # Hamiltonian basis is rank-deficient so each leaf can avoid the generic
    # transvection/completion mapper.
    basis_complete_tableau: np.ndarray | None
    basis_complete_inv: np.ndarray | None
    row_basis_cache: dict[str, np.ndarray]


def _solve_qubit_phase_family_correction(base_tableau: np.ndarray, delta: np.ndarray) -> np.ndarray | None:
    """
    Complete phase-correction solve for qubits (fixed symplectic F, fixed reference h0).

    We solve over the full admissible family h = h0 + 2u (mod 4), i.e.
        P (2u) = delta (mod 4)
    which is equivalent to:
        P u = delta/2 (mod 2)
    after checking that delta is even entrywise mod 4.

    Returns
    -------
    h_lin : np.ndarray | None
        A phase correction vector h_lin = 2u (mod 4), or None if no solution exists.
    """
    delta4 = np.asarray(delta, dtype=int).reshape(-1) % 4

    # Necessary and sufficient divisibility condition for 2u-correction in Z_4
    if np.any(delta4 % 2 != 0):
        return None

    # Reduce to GF(2) system: P u = (delta/2) mod 2
    A2 = (np.asarray(base_tableau, dtype=np.uint8) & 1)
    b2 = ((delta4 // 2) % 2).astype(np.uint8)

    u = solve_gf2(A2, b2)
    if u is None:
        return None

    # Lift back to a valid linear phase correction mod 4
    h_lin = (2 * np.asarray(u, dtype=int)) % 4
    return h_lin


def _target_completion_from_source_pairings(
    target_basis: np.ndarray,
    source_complete: np.ndarray,
    p: int,
) -> np.ndarray | None:
    """
    Complete target_basis so it has the same symplectic Gram matrix as
    source_complete.

    Rows are Pauli tableau vectors and maps act on the right.  If C is the
    completed source basis and D is the completed target basis with
    C Ω C^T = D Ω D^T, then F = C^{-1}D is symplectic and maps the constrained
    source rows to the target rows.
    """
    D = mod_p(np.asarray(target_basis, dtype=int), p)
    C = mod_p(np.asarray(source_complete, dtype=int), p)
    n2 = C.shape[1]
    if C.shape != (n2, n2) or D.ndim != 2 or D.shape[1] != n2:
        return None

    r = int(D.shape[0])
    if r > n2 or not np.array_equal(
        mod_p(D @ omega_matrix(n2 // 2, p) @ D.T, p),
        mod_p(C[:r] @ omega_matrix(n2 // 2, p) @ C[:r].T, p),
    ):
        return None

    Omega = omega_matrix(n2 // 2, p)
    for k in range(r, n2):
        A = mod_p(D @ Omega, p)
        b = mod_p(C[:k] @ Omega @ C[k].reshape(-1, 1), p).reshape(-1)
        try:
            candidate = _solve_linear(A, b, p).reshape(-1)
        except RuntimeError:
            return None

        if rank_mod(np.vstack([D, candidate]), p) <= k:
            null_basis = nullspace_mod(A, p)
            found = False
            for j in range(null_basis.shape[1]):
                trial = mod_p(candidate + null_basis[:, j], p)
                if rank_mod(np.vstack([D, trial]), p) > k:
                    candidate = trial
                    found = True
                    break
            if not found:
                return None

        D = np.vstack([D, mod_p(candidate, p)])

    if rank_mod(D, p) != n2:
        return None
    if not np.array_equal(mod_p(D @ Omega @ D.T, p), mod_p(C @ Omega @ C.T, p)):
        return None
    return mod_p(D, p)


def check_leaf(pi: np.ndarray, ctx: LeafContext) -> Gate | None:
    """Run all structural and phase-correction checks for a candidate permutation.

    Returns a symmetry Gate if successful, otherwise None.
    """
    if np.array_equal(pi, ctx.identity_perm):
        return None

    # (1) Edge-colour constraint: Pi must preserve the commutation matrix.
    if not np.array_equal(ctx.S_mod[np.ix_(pi, pi)], ctx.S_mod):
        return None

    # (2) Linear code  constraint.
    if not check_code_automorphism(ctx.G, ctx.basis_order, ctx.labels, pi, ctx.G_mod2):
        # print('FAILED CODE CHECK')
        return None

    # (3) Build the candidate symplectic from the permutation of the (ordered) row-basis.
    tgt_idx = pi[ctx.basis_indices]

    def _generic_map() -> tuple[np.ndarray, np.ndarray] | None:
        H_basis_tgt = PauliSum.from_tableau(
            ctx.base_tableau[tgt_idx],
            ctx.pauli_sum.dimensions,
            weights=ctx.base_weights[tgt_idx],
        )
        H_basis_tgt.set_phases(np.array(ctx.base_phases[tgt_idx], dtype=int, copy=True))
        try:
            F_gen, h_gen, _, _ = find_map_to_target_pauli_sum(ctx.basis_source_ps, H_basis_tgt)
        except Exception:
            return None
        return np.asarray(F_gen, dtype=int), np.asarray(h_gen, dtype=int)

    def _fast_map() -> tuple[np.ndarray, np.ndarray] | None:
        F: np.ndarray

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
            if ctx.basis_complete_tableau is None or ctx.basis_complete_inv is None:
                return None
            target_complete = _target_completion_from_source_pairings(
                ctx.base_tableau[tgt_idx],
                ctx.basis_complete_tableau,
                ctx.p,
            )
            if target_complete is None:
                return None
            F = mod_p(ctx.basis_complete_inv @ target_complete, ctx.p)
            if not np.array_equal(mod_p(ctx.basis_complete_tableau @ F, ctx.p), target_complete):
                return None

        # Deterministic lift from symplectic -> quadratic phase vector; linear correction solved below.
        h0 = get_phase_vector(F.T, int(ctx.pauli_sum.dimensions[0]))
        return F, np.asarray(h0, dtype=int)

    def _verify_map(F: np.ndarray, h0: np.ndarray) -> Gate | None:
        nq = ctx.n_qudits
        pauli = ctx.pauli_sum

        SG_F = Gate('Symmetry', F.T, np.asarray(h0, dtype=int))
        H_full_tg = pauli.copy()[pi]
        H_full_F = SG_F.act(pauli, tuple(range(nq)))

        delta = (H_full_tg.phases - H_full_F.phases) % ctx.two_lcm

        # Complete qubit phase correction (fixed F): solve over the full family h = h0 + 2u (mod 4)
        if ctx.p == 2:
            h_lin = _solve_qubit_phase_family_correction(ctx.base_tableau, delta)
        else:
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

    fast = _fast_map()
    if fast is not None:
        leaf = _verify_map(*fast)
        if leaf is not None:
            return leaf

    generic = _generic_map()
    if generic is None:
        return None
    return _verify_map(*generic)
