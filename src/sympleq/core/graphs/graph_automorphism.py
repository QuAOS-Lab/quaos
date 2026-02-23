from dataclasses import dataclass
import numpy as np
import galois
from numba import njit
from typing import Any, cast
from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.core.finite_field_solvers import get_linear_dependencies
from sympleq.core.circuits.target import find_map_to_target_pauli_sum, get_phase_vector
# from sympleq.core.circuits.find_symplectic import map_pauli_sum_to_target_tableau
from sympleq.core.finite_field_solvers import _select_row_basis_indices
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Gate
from sympleq.core.graphs.graph_coloring import _wl_colors_from_S
from sympleq.core.symmetries.phase_correction import solve_phase_vector_h_from_residual
from sympleq.core.graphs.dynamic_refine_search import dynamic_refine_individualized

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

def _coeff_ids(coeffs: np.ndarray | None) -> np.ndarray | None:
    if coeffs is None:
        return None
    # map arbitrary coefficient objects -> stable small ints
    _, ids = np.unique(coeffs, return_inverse=True)
    return ids.astype(np.int64)

def _choose_base_anchors(base_classes: dict[int, list[int]], max_anchors: int) -> list[int]:
    """
    Choose anchors from *small* color classes first (most discriminative).
    """
    anchors = []
    # colors sorted by class size ascending
    for c in sorted(base_classes.keys(), key=lambda cc: len(base_classes[cc])):
        if not base_classes[c]:
            continue
        anchors.append(int(base_classes[c][0]))
        if len(anchors) >= max_anchors:
            break
    return anchors

def _compute_anchor_hash(
    S_mod: np.ndarray,
    anchors: np.ndarray,          # shape (A,)
    base_colors: np.ndarray,      # shape (M,)
    coeff_id: np.ndarray | None,  # shape (M,) or None
    seed: int = 0,
) -> np.ndarray:
    """
    Compute a 64-bit hash per vertex from (base_color, coeff_id, S[v,anchors], S[anchors,v]).
    This is for ordering only; collisions are fine.
    Cost: O(M*A).
    """
    M = S_mod.shape[0]
    A = anchors.size
    # features: [base_color, coeff_id?, S_v_to_a..., S_a_to_v...]
    width = 1 + (1 if coeff_id is not None else 0) + 2 * A
    feats = np.empty((M, width), dtype=np.int16)

    col = 0
    feats[:, col] = base_colors.astype(np.int16, copy=False); col += 1
    if coeff_id is not None:
        feats[:, col] = coeff_id.astype(np.int16, copy=False); col += 1

    feats[:, col:col + A] = S_mod[:, anchors].astype(np.int16, copy=False); col += A
    feats[:, col:col + A] = S_mod[anchors, :].T.astype(np.int16, copy=False); col += A

    rng = np.random.default_rng(seed)
    w = rng.integers(1, np.iinfo(np.uint64).max, size=width, dtype=np.uint64)
    # 64-bit dot product hash
    h = (feats.astype(np.uint64) * w).sum(axis=1, dtype=np.uint64)
    return h



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
    ref_tableau: np.ndarray
    ref_phases: np.ndarray
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
                ctx: _LeafContext) -> Gate | None:
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

    # Build the candidate symplectic from the permutation of the (ordered) row-basis.
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
        # Rank-deficient tableau: the basis alone doesn't determine F, so we fall back
        # to the generic mapper.
        F, _h0, _, _ = find_map_to_target_pauli_sum(H_basis_src, H_basis_tgt)

    # Use a deterministic lift from symplectic -> quadratic phase vector; linear correction is solved below.
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

    h_lin = solve_phase_vector_h_from_residual(ctx.base_tableau, delta, ctx.pauli_sum.dimensions,
                                               debug=False, row_basis_cache=ctx.row_basis_cache)
    if h_lin is None:
        return None

    h0_mod = np.asarray(h0, dtype=int) % ctx.two_lcm
    h_lin_mod = np.asarray(h_lin, dtype=int) % ctx.two_lcm
    h_tot = (h0_mod + h_lin_mod) % ctx.two_lcm
    SG = Gate('Symmetry', F.T, h_tot)

    H_out_cf = SG.act(pauli, tuple(range(nq))).to_standard_form()
    H_out_cf.weight_to_phase()

    if not np.array_equal(H_out_cf.tableau, ctx.ref_tableau):
        return None
    if not np.all((ctx.ref_phases - H_out_cf.phases) % ctx.two_lcm == 0):
        return None
    if not np.all(np.isclose(H_out_cf.weights, ctx.ref_weights, atol=1e-8, rtol=0)):
        return None

    # unnecessary to check this separately since it's implied by the tableau and phase checks
    # Omega = np.zeros((2 * nq, 2 * nq), dtype=int)
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
    extra_column_invariants: str = "hist",
    p2_bitset: str | bool = "auto",
    color_mode: str = "wl",   # "wl" | "coeffs_only" | "none"
    max_wl_rounds: int = 10,
    progress: bool = False,
    progress_every: int = 2048,
) -> list[Gate]:
    """
    Find up to k automorphisms preserving S and the vector set.
    All preprocessing based only on pauli_sum is handled internally.
    Dynamic WL (if enabled) is used only for ordering every `dynamic_refine_every` steps.
    """

    # ---- preprocessing that depends only on pauli_sum ----
    # Choose a canonical gauge up front so coefficients that differ only by a discrete
    # Clifford phase are represented with identical weights and differing integer phases.
    pauli = pauli_sum.copy()
    pauli.weight_to_phase()

    independent_labels, dependencies = get_linear_dependencies(pauli.tableau, 2)
    labels = sorted(set(independent_labels) | set(dependencies.keys()))
    S_mod = pauli.symplectic_product_matrix()
    G, basis_order = pauli.matroid()
    coeffs = np.asarray(pauli.weights)

    if not np.all([pauli.dimensions[i] == pauli.dimensions[0] for i in range(1, len(pauli.dimensions))]):
        raise ValueError("All qubits must have same dimension for now. The key things to fix are: "
                         "_gf_solve_one_solution, and the symplectic_solver for F.")
    p = int(pauli.lcm)
    n = len(labels)

    G_mod2: np.ndarray | None = None
    if p == 2:
        G_mod2 = (np.asarray(G, dtype=np.uint8) & 1)
    col_invariants = None
    if extra_column_invariants != "none":
        # parse tokens like "hist+lc" or "hist,lc"
        toks = {t.strip().lower() for t in extra_column_invariants.replace(",", "+").split("+") if t.strip()}
        if "none" in toks:
            toks.remove("none")

        feats: list[np.ndarray] = []

        if "hist" in toks:
            inv_hist = np.zeros((n, min(p, 16)), dtype=np.int64)
            # Use an int view for bincount
            if p == 2 and G_mod2 is not None:
                G_for_hist = G_mod2
            else:
                G_for_hist = (np.asarray(G, dtype=int) % p)
            for j in range(n):
                col = np.asarray(G_for_hist[:, j], dtype=int)
                cnt = np.bincount(col, minlength=p)
                inv_hist[j, :min(p, 16)] = cnt[:min(p, 16)]
            feats.append(inv_hist)

        if ("lc" in toks) or ("loop_coloop" in toks) or ("loops_coloops" in toks):
            # loop: column is zero
            if p == 2 and G_mod2 is not None:
                is_loop = np.all(G_mod2 == 0, axis=0)
            else:
                G_int = (np.asarray(G, dtype=int) % p)
                is_loop = np.all(G_int == 0, axis=0)

            # coloop: only defined meaningfully relative to the chosen basis_order
            # We compute it cheaply via X = C^{-1} G in basis coordinates.
            is_coloop = np.zeros(n, dtype=bool)
            try:
                lab_to_idx = {lab: i for i, lab in enumerate(labels)}
                B_cols = np.array([lab_to_idx[b] for b in basis_order], dtype=int)
                B_mask = np.zeros(n, dtype=bool)
                B_mask[B_cols] = True
                nonB = np.where(~B_mask)[0]

                if p == 2 and G_mod2 is not None:
                    C = G_mod2[:, B_cols]
                    C_inv = _gf2_inv(C)
                    X = (C_inv @ G_mod2) & 1
                    # basis element is a coloop iff its row is zero on all non-basis columns
                    basis_row_used = np.any(X[:, nonB] != 0, axis=1) if nonB.size else np.zeros(X.shape[0], dtype=bool)
                else:
                    C = G[:, B_cols]
                    U = np.linalg.inv(C)  # works for galois.FieldArray
                    X = U @ G
                    X_np = np.asarray(X, dtype=int) % p
                    basis_row_used = np.any(X_np[:, nonB] != 0, axis=1) if nonB.size else np.zeros(X_np.shape[0], dtype=bool)

                is_coloop[B_cols] = ~basis_row_used
            except Exception:
                # If anything goes wrong (should be rare), fall back to "unknown" (no pruning).
                # This preserves correctness; it just weakens early pruning.
                pass

            inv_lc = np.stack([is_loop.astype(np.int64), is_coloop.astype(np.int64)], axis=1)
            feats.append(inv_lc)

        if toks - {"hist", "lc", "loop_coloop", "loops_coloops"}:
            raise ValueError("extra_column_invariants must be 'none', 'hist', 'lc', or 'hist+lc'.")

        if feats:
            col_invariants = np.hstack(feats) if len(feats) > 1 else feats[0]

    base_colors, base_classes = _build_base_partition(
        S_mod, p,
        coeffs=coeffs if coeffs is not None else None,
        col_invariants=col_invariants if color_mode == "wl" else None,
        max_rounds=max_wl_rounds,
        color_mode=color_mode,
    )

    coeff_id = _coeff_ids(coeffs)
    base_anchors = _choose_base_anchors(base_classes, max_anchors=32)  # tune 16–64
    anchors = np.array(base_anchors, dtype=np.int64)
    key_hash = _compute_anchor_hash(S_mod, anchors, base_colors, coeff_id, seed=0)

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
    pauli_standard = pauli.to_standard_form()
    pauli_standard.weight_to_phase()
    ref_tableau = pauli_standard.tableau.astype(int, copy=False)
    ref_phases = np.asarray(pauli_standard.phases, dtype=int)
    ref_weights = np.asarray(pauli_standard.weights)
    base_tableau = pauli.tableau.astype(int, copy=False)
    base_weights = np.asarray(pauli.weights)
    base_phases = np.asarray(pauli.phases, dtype=int)
    basis_indices = np.asarray(independent_labels, dtype=int)
    basis_source_ps = pauli[basis_indices]

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

    dims_array = np.asarray(pauli.dimensions, dtype=int)
    row_basis_cache: dict[str, np.ndarray] = {}
    if dims_array.size and np.all(dims_array == dims_array[0]):
        p_uni = int(dims_array[0])
        if p_uni == 2:
            row_basis_cache["gf2"] = _select_row_basis_indices(base_tableau % 2, 2, base_tableau.shape[1])
        else:
            row_basis_cache["gfp"] = _select_row_basis_indices(base_tableau % p_uni, p_uni, base_tableau.shape[1])

    phi = -np.ones(n, dtype=np.int64)
    used = np.zeros(n, dtype=bool)
    identity_perm = np.arange(pauli.n_paulis(), dtype=np.int64)

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

    def _bucket_key(idx: int) -> Any:
        if coeffs is None:
            return int(base_colors[idx])
        return (int(base_colors[idx]), coeffs[idx])

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
        two_lcm=2 * int(pauli.lcm),
        n_qudits=pauli.n_qudits(),
        identity_perm=identity_perm,
        S_mod=S_mod,
        G=G,
        G_mod2=G_mod2,
        basis_order=basis_order,
        labels=labels,
        pauli_sum=pauli,
        ref_tableau=ref_tableau,
        ref_phases=ref_phases,
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
        nonlocal anchors, key_hash
        if dynamic_refine_every <= 0:
            return
        # augment anchors with up to K mapped domain vertices (individualization)
        mapped = np.where(phi >= 0)[0]
        K = 16  # add at most 16 individualized anchors per refresh
        extra = mapped[:K].astype(np.int64, copy=False)
        anchors_new = np.unique(np.concatenate([anchors, extra]))
        # cap total anchors
        Amax = 64
        if anchors_new.size > Amax:
            anchors_new = anchors_new[:Amax]
        anchors = anchors_new
        key_hash = _compute_anchor_hash(S_mod, anchors, base_colors, coeff_id, seed=0)

    @dataclass
    class _DFSFrame:
        i: int
        bi: int
        mapped_idx: np.ndarray
        candidate: list[int]
        branch_leaf_mass: int = 0
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
        bkey = _bucket_key(i)
        mapped_idx = np.where(phi >= 0)[0].astype(np.int64)

        candidate = [y for y in base_classes[bi] if not used[y]]
        if coeffs is not None:
            candidate = [y for y in candidate if coeffs[i] == coeffs[y]]
        cand = np.array(candidate, dtype=np.int64)
        cand = cand[np.argsort(key_hash[cand], kind="mergesort")]
        candidate = cand.tolist()

        branch_leaf_mass = 0
        if progress_enabled:
            rem_i = int(rem_counts.get(bkey, 0))
            if rem_i > 0:
                state_leaf_mass = 1
                for cnt in rem_counts.values():
                    state_leaf_mass *= factorials[int(cnt)]
                branch_leaf_mass = state_leaf_mass // rem_i

        return _DFSFrame(
            i=i,
            bi=bi,
            mapped_idx=mapped_idx,
            candidate=candidate,
            branch_leaf_mass=int(branch_leaf_mass),
        )

    # Iterative DFS
    progress_bar = None
    pending_progress = 0
    leaves_checked = 0
    progress_enabled = bool(progress)
    factorials: list[int] = []
    total_leaf_space = 0
    explored_leaf_space = 0
    if progress_enabled:
        factorials = [1] * (n + 1)
        for i in range(2, n + 1):
            factorials[i] = factorials[i - 1] * i
        total_leaf_space = 1
        for cnt in rem_counts.values():
            total_leaf_space *= factorials[int(cnt)]

    def _progress_pct_string() -> str:
        if total_leaf_space <= 0:
            return "0.00%"
        pct_x100 = (int(explored_leaf_space) * 10000) // int(total_leaf_space)
        return f"{pct_x100 / 100:.2f}%"

    if progress_enabled:
        try:
            from tqdm.auto import tqdm
            progress_bar = tqdm(
                total=None,
                desc="Clifford automorphism search",
                unit="node",
                leave=False,
                dynamic_ncols=True,
                mininterval=0.2,
            )
            progress_bar.set_postfix(found=0, leaves=0, searched="0.00%", refresh=False)
        except Exception:
            progress_bar = None

    stack: list[_DFSFrame] = []
    try:
        while True:
            if progress_bar is not None:
                pending_progress += 1
                if pending_progress >= progress_every:
                    progress_bar.update(pending_progress)
                    pending_progress = 0
                    progress_bar.set_postfix(
                        found=len(results),
                        leaves=leaves_checked,
                        searched=_progress_pct_string(),
                        refresh=False,
                    )

            if len(results) >= k_wanted:
                break

            # Leaf check
            if np.all(phi >= 0):
                leaves_checked += 1
                if progress_enabled:
                    explored_leaf_space += 1
                pi = phi.copy()
                leaf = _check_leaf(pi, leaf_ctx)
                if leaf is not None:
                    results.append(leaf)
                    break  # match previous behavior: stop after the first found symmetry
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
                    if progress_enabled and frame.branch_leaf_mass > 0:
                        explored_leaf_space += frame.branch_leaf_mass
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
    finally:
        if progress_bar is not None:
            if pending_progress:
                progress_bar.update(pending_progress)
            progress_bar.set_postfix(
                found=len(results),
                leaves=leaves_checked,
                searched=_progress_pct_string(),
                refresh=False,
            )
            progress_bar.close()

    return results[:k_wanted]
