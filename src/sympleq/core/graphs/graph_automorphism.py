from dataclasses import dataclass
import numpy as np
import galois
from numba import njit
from typing import Any, cast
from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.core.finite_field_solvers import get_linear_dependencies
from sympleq.core.circuits.target import find_map_to_target_pauli_sum, get_phase_vector
from sympleq.core.circuits.find_symplectic import map_pauli_sum_to_target_tableau
from sympleq.core.finite_field_solvers import _select_row_basis_indices
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Gate
from sympleq.core.graphs.graph_coloring import _wl_colors_from_S
from sympleq.core.symmetries.phase_correction import solve_phase_vector_h_from_residual


@njit(cache=True, fastmath=True)
def _consistent_numba(S_mod: np.ndarray, phi: np.ndarray, mapped_idx: np.ndarray, i: int, y: int) -> bool:
    """
    Require S[i, j] == S[y, phi[j]] and S[j, i] == S[phi[j], y] for all mapped j.
    """
    for t in range(mapped_idx.size):
        j = mapped_idx[t]
        yj = phi[j]
        if S_mod[i, j] != S_mod[y, yj]:
            return False
        if S_mod[j, i] != S_mod[yj, y]:
            return False
    return True


def _build_bitrows_binary(S_mod: np.ndarray) -> tuple[np.ndarray, int]:
    """
    For p=2 only. Pack each row's 0/1 into chunks of 64 bits.
    Returns (bits[n, C], chunks=C). Column j lives at chunk=j>>6, bit=(j & 63).
    """
    n = S_mod.shape[0]
    C = (n + 63) // 64
    bits = np.zeros((n, C), dtype=np.uint64)
    for i in range(n):
        for j in range(n):
            if S_mod[i, j] & 1:
                bits[i, j >> 6] |= (1 << (j & 63))
    return bits, C


@njit(cache=True, fastmath=True)
def _consistent_bitset(bits: np.ndarray, phi: np.ndarray, mapped_idx: np.ndarray, i: int, y: int) -> bool:
    """
    Same logic as _consistent_numba but reading single bits from packed rows.
    """
    for t in range(mapped_idx.size):
        j = mapped_idx[t]
        yj = phi[j]
        # read bit S[i,j]
        bi = (bits[i, j >> 6] >> (j & 63)) & 1
        # read bit S[y,yj]
        by = (bits[y, yj >> 6] >> (yj & 63)) & 1
        if bi != by:
            return False
        # read bit S[j,i] vs S[yj,y]
        bji = (bits[j, i >> 6] >> (i & 63)) & 1
        byy = (bits[yj, y >> 6] >> (y & 63)) & 1
        if bji != byy:
            return False
    return True


class _ConsistencyChecker:
    """Reusable consistency kernel; chooses bitset or direct variant once."""

    def __init__(self, S_mod: np.ndarray, p2_bitset: bool):
        self.S_mod = S_mod
        if p2_bitset:
            self.bits, _ = _build_bitrows_binary(S_mod)
            self._fn = self._bitset
        else:
            self._fn = self._direct

    def __call__(self, phi: np.ndarray, mapped_idx: np.ndarray, i: int, y: int) -> bool:
        return self._fn(phi, mapped_idx, int(i), int(y))

    def _bitset(self, phi: np.ndarray, mapped_idx: np.ndarray, i: int, y: int) -> bool:
        return _consistent_bitset(self.bits, phi, mapped_idx, i, y)

    def _direct(self, phi: np.ndarray, mapped_idx: np.ndarray, i: int, y: int) -> bool:
        return _consistent_numba(self.S_mod, phi, mapped_idx, i, y)


def _gf2_inv(M: np.ndarray) -> np.ndarray:
    """Invert a square GF(2) matrix using XOR elimination; raises LinAlgError if singular."""
    A = (M.copy().astype(np.uint8) & 1)
    n = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise np.linalg.LinAlgError("GF2 inverse requires square matrix")
    Id = np.eye(n, dtype=np.uint8)
    aug = np.hstack((A, Id))
    row = 0
    for col in range(n):
        if row >= n:
            break
        nz = np.flatnonzero(aug[row:, col])
        if nz.size == 0:
            continue
        piv = row + nz[0]
        if piv != row:
            aug[[row, piv]] = aug[[piv, row]]
        mask = aug[:, col].astype(bool)
        mask[row] = False
        aug[mask] ^= aug[row]
        row += 1
    if not np.array_equal(aug[:, :n] & 1, np.eye(n, dtype=np.uint8)):
        raise np.linalg.LinAlgError("matrix is singular over GF(2)")
    return aug[:, n:] & 1


def _check_code_automorphism(
    G: galois.FieldArray,
    basis_order: list[int],
    labels: list[int],
    pi: np.ndarray,
    G_mod2: np.ndarray | None = None,
) -> bool:
    """
    Linear-code test over GF(p): there exists U with U G P = G ?
    Let C = G[:, P(B)]; if invertible, U = C^{-1} and check U G P == G.
    """
    lab_to_idx = {lab: i for i, lab in enumerate(labels)}
    B_cols = np.array([lab_to_idx[b] for b in basis_order], dtype=int)
    PBcols = pi[B_cols]
    if G_mod2 is not None:
        C = G_mod2[:, PBcols]
        try:
            C_inv = _gf2_inv(C)
        except np.linalg.LinAlgError:
            return False
        Gp = G_mod2[:, pi]
        return np.array_equal((C_inv @ Gp) & 1, G_mod2)
    else:
        C = G[:, PBcols]
        # fall back to galois / numpy inverse; accept LinAlgError as failure
        try:
            U = np.linalg.inv(C)  # works on galois.FieldArray
        except np.linalg.LinAlgError:
            return False
        Gp = G[:, pi]
        return np.array_equal(U @ Gp, G)


@dataclass
class _LeafContext:
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
    # Coefficient-normalized copy used in leaf checks: phases are absorbed into weights
    # so we can compare full complex coefficients without any weight_to_phase ambiguity.
    pauli_coeff: PauliSum
    ref_tableau: np.ndarray
    ref_weights: np.ndarray
    base_tableau: np.ndarray
    base_weights: np.ndarray
    base_phases: np.ndarray
    basis_indices: np.ndarray
    basis_source_ps: PauliSum
    # Precomputed inverse of the square row-basis tableau (when rank == 2*n_qudits).
    # This makes the leaf mapping deterministic (no order-dependent transvections).
    basis_src_inv_gf2: np.ndarray | None
    basis_src_inv_gfp: galois.FieldArray | None
    row_basis_cache: dict[str, np.ndarray]


def _check_leaf(pi: np.ndarray,
                ctx: _LeafContext,
                known_F: np.ndarray | None = None,
                fail_loudly: bool = False) -> Gate | None:
    """
    Run all structural and phase-correction checks for a candidate permutation.
    Returns a symmetry Gate or None.
    """
    if np.array_equal(pi, ctx.identity_perm):
        return None
    if not np.array_equal(ctx.S_mod[np.ix_(pi, pi)], ctx.S_mod):
        return None
    if not _check_code_automorphism(ctx.G, ctx.basis_order, ctx.labels, pi, ctx.G_mod2):
        return None

    # For constructing candidate symplectics we must be careful about *row order*.
    #
    # - Full-rank case: we precompute an inverse for ctx.basis_indices, so we must
    #   use that exact ordered basis when constructing F.
    # - Rank-deficient case: transvection-based mapping is order-dependent, so we
    #   prefer a stable pivot-order basis (cached) over the independent-label list.
    use_precomputed_inv = (
        (ctx.p == 2 and ctx.basis_src_inv_gf2 is not None)
        or (ctx.p != 2 and ctx.basis_src_inv_gfp is not None)
    )
    basis_rows = ctx.basis_indices
    if not use_precomputed_inv and ctx.row_basis_cache is not None:
        key = "gf2" if ctx.p == 2 else "gfp"
        rb = ctx.row_basis_cache.get(key)
        if rb is not None and rb.size:
            basis_rows = rb.astype(int, copy=False)

    tgt_idx = pi[basis_rows]
    H_basis_tgt = PauliSum.from_tableau(
        ctx.base_tableau[tgt_idx],
        ctx.pauli_sum.dimensions,
        weights=ctx.base_weights[tgt_idx],
    )
    H_basis_tgt.set_phases(np.array(ctx.base_phases[tgt_idx], dtype=int, copy=True))
    H_basis_src = PauliSum.from_tableau(
        ctx.base_tableau[basis_rows],
        ctx.pauli_sum.dimensions,
        weights=ctx.base_weights[basis_rows],
    )
    H_basis_src.set_phases(np.array(ctx.base_phases[basis_rows], dtype=int, copy=True))

    # Derive the symplectic map from the row-basis permutation whenever possible.
    # This is the most direct "basis -> basis" construction and avoids order-dependent
    # transvection sequences.
    F: np.ndarray
    h0: np.ndarray
    if ctx.p == 2 and ctx.basis_src_inv_gf2 is not None:
        T = (ctx.base_tableau[tgt_idx] & 1).astype(np.uint8, copy=False)
        F = (ctx.basis_src_inv_gf2 @ T) & 1
        F = np.asarray(F, dtype=int)
        # get_phase_vector expects Gate.symplectic (not the right-action matrix F).
        h0 = get_phase_vector(F.T, int(ctx.pauli_sum.dimensions[0]))
    elif ctx.p != 2 and ctx.basis_src_inv_gfp is not None:
        GF = galois.GF(int(ctx.p))
        T = GF(ctx.base_tableau[tgt_idx] % ctx.p)
        F_gf = ctx.basis_src_inv_gfp @ T
        F = (np.asarray(F_gf, dtype=int) % ctx.p).astype(int, copy=False)
        h0 = get_phase_vector(F.T, int(ctx.pauli_sum.dimensions[0]))
    else:
        # Rank-deficient tableau: the basis alone doesn't determine F.
        # Build an F from a pivot-order basis mapping; this tends to recover a
        # "good" symplectic completion that makes phase correction solvable.
        if ctx.p == 2:
            src_basis = (H_basis_src.tableau & 1).astype(int, copy=False)
            tgt_basis = (H_basis_tgt.tableau & 1).astype(int, copy=False)
            F = map_pauli_sum_to_target_tableau(src_basis, tgt_basis)

            # Sanity-check: the candidate must implement the full row mapping implied by pi.
            # If not, fall back to a basis-first full mapping (rare; mostly for pathological
            # dependency representations).
            base = ctx.base_tableau & 1
            tgt_full = ctx.base_tableau[pi] & 1
            if not np.array_equal((base @ F) & 1, tgt_full):
                all_idx = np.arange(base.shape[0], dtype=int)
                basis_set = set(int(x) for x in basis_rows)
                rest = np.array([i for i in all_idx if i not in basis_set], dtype=int)
                order = np.concatenate([np.asarray(basis_rows, dtype=int), rest])
                F = map_pauli_sum_to_target_tableau(base[order], tgt_full[order])

            h0 = get_phase_vector(F.T, int(ctx.pauli_sum.dimensions[0]))
        else:
            # Generic fallback for non-qubit rank-deficient cases.
            F, h0, _, _ = find_map_to_target_pauli_sum(H_basis_src, H_basis_tgt)

    if fail_loudly and known_F is not None and not np.array_equal(F.T, known_F):
        print(f'[DEBUG] Found correct permutation, but incorrect F.')
        print('known_F:\n', known_F)
        print('candidate F:\n', F.T)
        print('basis rows used:\n', basis_rows)
        print('input basis:\n', H_basis_src.tableau)
        print('target indices:\n', tgt_idx)
        print('target basis:\n', H_basis_tgt.tableau)
        print('H_tableau (base order):\n', ctx.base_tableau)

    nq = ctx.n_qudits
    pauli_coeff = ctx.pauli_coeff

    SG_F = Gate('Symmetry', list(range(nq)), F.T, ctx.pauli_sum.dimensions, np.asarray(h0, dtype=int))
    H_full_F = SG_F.act(pauli_coeff)

    # We want coefficients to match under the permutation. With phases absorbed into weights,
    # the coefficient of term i after acting is w_i * omega^{phase_i}, while the target is w_{pi[i]}.
    # Solve for a per-term residual phase b_i such that omega^{b_i} = w_{pi[i]} / (w_i * omega^{phase_i}).
    two_lcm = 2 * int(ctx.pauli_sum.lcm)
    w_src = np.asarray(pauli_coeff.weights, dtype=np.complex128)
    w_tgt = w_src[pi]
    phase_out = np.asarray(H_full_F.phases, dtype=int) % two_lcm

    if two_lcm == 4:
        omega_pows = np.array([1.0 + 0.0j, 0.0 + 1.0j, -1.0 + 0.0j, 0.0 - 1.0j], dtype=np.complex128)
        denom = w_src * omega_pows[phase_out]
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(np.abs(denom) > 1e-15, w_tgt / denom, 1.0 + 0.0j)
        # Pick the closest 4th-root; reject if not close to any.
        dist = np.abs(ratio.reshape(-1, 1) - omega_pows.reshape(1, -1))
        b = np.argmin(dist, axis=1).astype(int)
        if fail_loudly:
            bad = np.where(dist[np.arange(dist.shape[0]), b] > 1e-6)[0]
            if bad.size:
                print("[DEBUG] coefficient ratio not a 4th root of unity at rows:", bad[:10])
        if np.any(dist[np.arange(dist.shape[0]), b] > 1e-6):
            return None
        delta = b % two_lcm
    else:
        # Generic (slower) path: use a small neighborhood around the nearest angle.
        omega = np.exp(2j * np.pi / two_lcm)
        denom = w_src * (omega ** phase_out)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(np.abs(denom) > 1e-15, w_tgt / denom, 1.0 + 0.0j)
        theta = np.angle(ratio)
        k0 = (np.round((two_lcm * theta) / (2.0 * np.pi)).astype(int) % two_lcm)
        candidates = np.stack([(k0 - 1) % two_lcm, k0, (k0 + 1) % two_lcm], axis=1)
        omega_pow = omega ** candidates
        dist = np.abs(ratio.reshape(-1, 1) - omega_pow)
        pick = np.argmin(dist, axis=1)
        delta = candidates[np.arange(candidates.shape[0]), pick] % two_lcm
        if np.any(dist[np.arange(dist.shape[0]), pick] > 1e-6):
            return None

    # Historical note: an older phase-correction path assumed only even corrections were possible
    # for qubits and attempted a second h0 guess. We now solve directly over Z_4, so we keep a
    # single baseline h0 and let the linear solve handle the remaining gauge freedom.

    h_lin = solve_phase_vector_h_from_residual(ctx.base_tableau, delta, ctx.pauli_sum.dimensions,
                                               debug=False, row_basis_cache=ctx.row_basis_cache)
    if h_lin is None:
        if fail_loudly:
            print('[DEBUG] Failed to find phase correction.')
        return None

    h0_mod = np.asarray(h0, dtype=int) % ctx.two_lcm
    h_lin_mod = np.asarray(h_lin, dtype=int) % ctx.two_lcm
    h_tot = (h0_mod + h_lin_mod) % ctx.two_lcm
    SG = Gate('Symmetry', list(range(nq)), F.T, ctx.pauli_sum.dimensions, h_tot)

    # Compare in coefficient form (phases absorbed into weights) for robustness.
    H_out_cf = SG.act(pauli_coeff).to_standard_form()

    if not np.array_equal(H_out_cf.tableau, ctx.ref_tableau):
        if fail_loudly:
            print('[DEBUG] tableau mismatch.')
        return None
    if not np.all(np.isclose(H_out_cf.weights, ctx.ref_weights, atol=1e-8, rtol=0)):
        if fail_loudly:
            print('[DEBUG] weight mismatch.')
        return None
    # Omega = np.zeros((2 * nq, 2 * nq), dtype=int)   # This has caused rare failures on known symmetries ... ???
    # Omega[:nq, nq:] = np.eye(nq, dtype=int)
    # Omega[nq:, :nq] = -np.eye(nq, dtype=int)
    # if np.all(((ctx.p - 1) * np.diag(F @ Omega @ F.T) + h_tot % 2) != 0):
    #     if fail_loudly:
    #         print('[DEBUG] Inconsistent phase correction.')
    #     return None

    return SG


def clifford_graph_automorphism_search(
    pauli_sum: PauliSum,
    k_wanted: int,
    dynamic_refine_every: int = 0,
    extra_column_invariants: str = "none",
    p2_bitset: str | bool = "auto",
    color_mode: str = "wl",   # "wl" | "coeffs_only" | "none"
    max_wl_rounds: int = 10,
    known_F: np.ndarray | None = None,
    debug_permutation: list[int] | None = None,
) -> list[Gate]:
    """
    Find up to k automorphisms preserving S and the vector set.
    All preprocessing based only on pauli_sum is handled internally.
    Dynamic WL (if enabled) is used only for ordering every `dynamic_refine_every` steps.
    """

    if known_F is not None and debug_permutation is None:
        # Permutation convention in this file: pi maps a source term index i to the target index pi[i]
        # such that acting by the symmetry sends row i -> row pi[i].
        #
        # If we build Tb = Ta @ F^T, then Tb[i] is the image of Ta[i], so we want pi with
        # Tb[i] == Ta[pi[i]]. Hence tableau_permutation(Tb, Ta) (not the inverse).
        debug_permutation = tableau_permutation(pauli_sum.tableau @ known_F.T % 2, pauli_sum.tableau)

    # ---- preprocessing that depends only on pauli_sum ----
    independent_labels, dependencies = get_linear_dependencies(pauli_sum.tableau, 2)
    labels = sorted(set(independent_labels) | set(dependencies.keys()))
    S_mod = pauli_sum.symplectic_product_matrix()
    G, basis_order = pauli_sum.matroid()
    # Use coefficient magnitudes only: Clifford symmetries can change phases of terms via the
    # phase-vector correction, but cannot change magnitudes.
    coeffs = np.abs(pauli_sum.weights)

    if not np.all([pauli_sum.dimensions[i] == pauli_sum.dimensions[0] for i in range(1, len(pauli_sum.dimensions))]):
        raise ValueError("All qubits must have same dimension for now. The key things to fix are: "
                         "_gf_solve_one_solution, and the symplectic_solver for F.")
    p = int(pauli_sum.lcm)
    n = len(labels)

    col_invariants = None
    if extra_column_invariants != "none":
        G_for_inv = G.copy()
        if extra_column_invariants == "hist":
            inv = np.zeros((n, min(p, 16)), dtype=np.int64)
            for j in range(n):
                col = np.array([int(x) for x in G_for_inv[:, j]])
                cnt = np.bincount(col, minlength=p)
                inv[j, :min(p, 16)] = cnt[:min(p, 16)]
            col_invariants = inv
        else:
            raise ValueError("extra_column_invariants must be 'none' or 'hist'.")

    base_colors, base_classes = _build_base_partition(
        S_mod, p,
        coeffs=coeffs if coeffs is not None else None,
        col_invariants=col_invariants if color_mode == "wl" else None,
        max_rounds=max_wl_rounds,
        color_mode=color_mode,
    )

    use_bitset = (p == 2 and (p2_bitset is True or (p2_bitset == "auto" and n <= 256)))

    results = []

    # Prepared consistency kernel
    consistency = _ConsistencyChecker(S_mod, use_bitset)

    # Static domain order from base classes (largest first)
    base_order = sorted(base_classes.keys(), key=lambda c: -len(base_classes[c]))
    domain_order = [i for c in base_order for i in base_classes[c]]

    # Track remaining candidates per color (and coeff if present) to avoid repeated scans in select_next
    rem_counts: dict[Any, int] = {}
    if coeffs is None:
        rem_counts.update({c: len(base_classes[c]) for c in base_classes})
    else:
        for idx in range(n):
            key = (int(base_colors[idx]), coeffs[idx])
            rem_counts[key] = rem_counts.get(key, 0) + 1

    # Prepare references for the leaf checks so they dont have to be computed every time.
    # We do all comparisons in coefficient form (phases absorbed into weights), which
    # avoids representation ambiguity when weights have arbitrary complex arguments.
    pauli_coeff = pauli_sum.copy()
    pauli_coeff.phase_to_weight()
    pauli_standard = pauli_coeff.to_standard_form()
    ref_tableau = pauli_standard.tableau.astype(int, copy=False)
    ref_weights = np.asarray(pauli_standard.weights)
    base_tableau = pauli_sum.tableau.astype(int, copy=False)
    base_weights = np.asarray(pauli_sum.weights)
    base_phases = np.asarray(pauli_sum.phases, dtype=int)
    basis_indices = np.asarray(independent_labels, dtype=int)
    basis_source_ps = pauli_sum[basis_indices]

    # If the tableau has full rank (rank == 2*n_qudits), the basis mapping uniquely determines
    # the symplectic matrix F via a single inversion. Precompute that inverse once.
    basis_src_inv_gf2: np.ndarray | None = None
    basis_src_inv_gfp: galois.FieldArray | None = None
    basis_src = np.asarray(basis_source_ps.tableau, dtype=int)
    if basis_src.shape[0] == basis_src.shape[1]:
        if p == 2:
            try:
                basis_src_inv_gf2 = _gf2_inv(basis_src & 1)
            except np.linalg.LinAlgError:
                basis_src_inv_gf2 = None
        else:
            GF = galois.GF(int(p))
            try:
                # galois overloads numpy.linalg for FieldArray, but type checkers often
                # don't understand that and infer a float ndarray.
                basis_src_inv_gfp = cast(galois.FieldArray, np.linalg.inv(GF(basis_src % p)))
            except np.linalg.LinAlgError:
                basis_src_inv_gfp = None

    dims_array = np.asarray(pauli_sum.dimensions, dtype=int)
    row_basis_cache: dict[str, np.ndarray] = {}
    if dims_array.size and np.all(dims_array == dims_array[0]):
        p_uni = int(dims_array[0])
        if p_uni == 2:
            row_basis_cache["gf2"] = _select_row_basis_indices(base_tableau % 2, 2, base_tableau.shape[1])
        else:
            row_basis_cache["gfp"] = _select_row_basis_indices(base_tableau % p_uni, p_uni, base_tableau.shape[1])

    phi = -np.ones(n, dtype=np.int64)
    used = np.zeros(n, dtype=bool)
    identity_perm = np.arange(pauli_sum.n_paulis(), dtype=np.int64)
    target_pi = None
    if debug_permutation is not None:
        target_pi = np.asarray(debug_permutation, dtype=np.int64)
        if target_pi.shape[0] != pauli_sum.n_paulis():
            raise ValueError("debug_permutation length must equal number of Pauli terms.")

    steps = 0
    cur_colors = base_colors.copy()

    def _dec_count(y_idx: int):
        if coeffs is None:
            rem_counts[int(base_colors[y_idx])] -= 1
        else:
            key = (int(base_colors[y_idx]), coeffs[y_idx])
            rem_counts[key] -= 1

    def _inc_count(y_idx: int):
        if coeffs is None:
            rem_counts[int(base_colors[y_idx])] += 1
        else:
            key = (int(base_colors[y_idx]), coeffs[y_idx])
            rem_counts[key] += 1

    def select_next() -> int:
        # MRV measured against remaining count in the relevant color/coeff bucket
        best_i, best_rem = -1, 10**9
        for i in domain_order:
            if phi[i] >= 0:
                continue
            if coeffs is None:
                rem = rem_counts[int(base_colors[i])]
            else:
                rem = rem_counts.get((int(base_colors[i]), coeffs[i]), 0)
            if rem < best_rem:
                best_i, best_rem = i, rem
                if rem <= 1:
                    break
        return best_i

    # For GF(2), carry an int copy to avoid FieldArray overhead in the code-automorphism test
    G_mod2 = None
    if p == 2:
        G_mod2 = (np.asarray(G, dtype=np.uint8) & 1)

    leaf_ctx = _LeafContext(
        p=p,
        two_lcm=2 * int(pauli_sum.lcm),
        n_qudits=pauli_sum.n_qudits(),
        identity_perm=identity_perm,
        S_mod=S_mod,
        G=G,
        G_mod2=G_mod2,
        basis_order=basis_order,
        labels=labels,
        pauli_sum=pauli_sum,
        pauli_coeff=pauli_coeff,
        ref_tableau=ref_tableau,
        ref_weights=ref_weights,
        base_tableau=base_tableau,
        base_weights=base_weights,
        base_phases=base_phases,
        basis_indices=basis_indices,
        basis_source_ps=basis_source_ps,
        basis_src_inv_gf2=basis_src_inv_gf2,
        basis_src_inv_gfp=basis_src_inv_gfp,
        row_basis_cache=row_basis_cache,
    )

    def dynamic_refine():
        nonlocal cur_colors
        if dynamic_refine_every <= 0:
            return
        # 1-WL just to order
        cur_colors = _wl_colors_from_S(S_mod, int(2), coeffs=coeffs, col_invariants=None, max_rounds=1)

    @dataclass
    class _DFSFrame:
        i: int
        bi: int
        mapped_idx: np.ndarray
        candidate: list[int]
        idx: int = 0  # next candidate index to try
        assigned_y: int = -1  # -1 means "unassigned"

    def _undo_assignment(frame: _DFSFrame) -> None:
        y = int(frame.assigned_y)
        if y < 0:
            return
        phi[frame.i] = -1
        used[y] = False
        _inc_count(y)
        frame.assigned_y = -1

    def _make_frame() -> _DFSFrame | None:
        nonlocal steps
        if len(results) >= k_wanted:
            return None

        if dynamic_refine_every and (steps % dynamic_refine_every == 0):
            dynamic_refine()
        steps += 1

        i = int(select_next())
        if i < 0:
            return None
        bi = int(base_colors[i])
        mapped_idx = np.where(phi >= 0)[0].astype(np.int64)

        # Order candidates by current colors (ordering heuristic only)
        if target_pi is not None:
            targ_y = int(target_pi[i])
            if used[targ_y]:
                print(f"[DEBUG] target impossible at i={i}: target {targ_y} already used")
                candidate: list[int] = []
            else:
                # bypass base-class and coeff filtering in debug mode
                candidate = [targ_y]
        else:
            candidate = [y for y in base_classes[bi] if not used[y]]
            if coeffs is not None:
                candidate = [y for y in candidate if coeffs[i] == coeffs[y]]
            candidate.sort(key=lambda y: cur_colors[y])

        return _DFSFrame(i=i, bi=bi, mapped_idx=mapped_idx, candidate=candidate)

    # Iterative DFS (explicit stack) to avoid hitting Python's recursion limit for large Pauli sums.
    stack: list[_DFSFrame] = []
    while True:
        if len(results) >= k_wanted:
            break

        # Debug pruning: current partial assignment disagrees with the target permutation.
        if target_pi is not None and np.any((phi >= 0) & (phi != target_pi)):
            print(f"[DEBUG] prune: assignment {phi} conflicts with target {target_pi}")
            # backtrack to the most recent assigned variable
            while stack and stack[-1].assigned_y < 0:
                stack.pop()
            if not stack:
                break
            _undo_assignment(stack[-1])
            continue

        # Leaf check
        if np.all(phi >= 0):
            pi = phi.copy()
            fail_loud = target_pi is not None and np.array_equal(pi, target_pi)
            leaf = _check_leaf(pi, leaf_ctx, known_F=known_F, fail_loudly=fail_loud)
            if leaf is not None:
                results.append(leaf)
                break  # match previous behavior: stop after the first found symmetry
            if fail_loud:
                print("[DEBUG] reached target permutation but leaf check failed")
            # leaf failed -> backtrack one level
            if not stack:
                break
            _undo_assignment(stack[-1])
            continue

        # Ensure there's a frame for the next variable.
        if not stack or stack[-1].assigned_y >= 0:
            fr = _make_frame()
            if fr is None:
                # No variable to assign or no candidates: fail this branch.
                if not stack:
                    break
                _undo_assignment(stack[-1])
                continue
            stack.append(fr)

        frame = stack[-1]

        # Try candidates for this frame's variable i.
        assigned = False
        while frame.idx < len(frame.candidate):
            y = int(frame.candidate[frame.idx])
            frame.idx += 1
            if used[y]:
                continue
            if not consistency(phi, frame.mapped_idx, frame.i, y):
                if target_pi is not None:
                    print(f"[DEBUG] target candidate fails consistency at i={frame.i}, y={y}")
                continue

            phi[frame.i] = y
            used[y] = True
            _dec_count(y)
            frame.assigned_y = y
            assigned = True
            break

        if assigned:
            # descend; next loop iteration will create/advance the next frame
            continue

        # No candidates left for this variable -> pop frame (no assignment) and backtrack.
        stack.pop()
        if not stack:
            break
        _undo_assignment(stack[-1])

    return results[:k_wanted]


def tableau_permutation(Ta, Tb):
    """
    Debug helper
    """
    if Ta.shape != Tb.shape:
        return None
    used = set()
    pi = [-1]*Ta.shape[0]
    for i in range(Ta.shape[0]):
        for j in range(Tb.shape[0]):
            if j in used:
                continue
            if np.array_equal(Ta[i], Tb[j]):
                pi[i] = j
                used.add(j)
                break
    return None if -1 in pi else pi
