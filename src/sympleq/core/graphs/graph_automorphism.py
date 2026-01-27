from dataclasses import dataclass
import numpy as np
import galois
from numba import njit
from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.core.finite_field_solvers import get_linear_dependencies
from sympleq.core.circuits.target import find_map_to_target_pauli_sum
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


def _check_code_automorphism(
    G: galois.FieldArray,
    basis_order: list[int],
    labels: list[int],
    pi: np.ndarray
) -> bool:
    """
    Linear-code test over GF(p): there exists U with U G P = G ?
    Let C = G[:, P(B)]; if invertible, U = C^{-1} and check U G P == G.
    """
    lab_to_idx = {lab: i for i, lab in enumerate(labels)}
    B_cols = np.array([lab_to_idx[b] for b in basis_order], dtype=int)
    PBcols = pi[B_cols]
    C = G[:, PBcols]
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
    basis_order: list[int]
    labels: list[int]
    pauli_sum: PauliSum
    pauli_weighted: PauliSum
    ref_tableau: np.ndarray
    ref_phases: np.ndarray
    ref_weights: np.ndarray
    base_tableau: np.ndarray
    base_weights: np.ndarray
    base_phases: np.ndarray
    basis_indices: np.ndarray
    basis_source_ps: PauliSum
    row_basis_cache: dict[str, np.ndarray]


def _check_leaf(pi: np.ndarray, ctx: _LeafContext) -> Gate | None:
    """
    Run all structural and phase-correction checks for a candidate permutation.
    Returns a symmetry Gate or None.
    """
    if np.array_equal(pi, ctx.identity_perm):
        return None
    if not np.array_equal(ctx.S_mod[np.ix_(pi, pi)], ctx.S_mod):
        return None
    if not _check_code_automorphism(ctx.G, ctx.basis_order, ctx.labels, pi):
        return None

    H_basis_src = ctx.basis_source_ps
    tgt_idx = pi[ctx.basis_indices]
    H_basis_tgt = PauliSum.from_tableau(ctx.base_tableau[tgt_idx], ctx.pauli_sum.dimensions,
                                        weights=ctx.base_weights[tgt_idx])
    H_basis_tgt.set_phases(np.array(ctx.base_phases[tgt_idx], dtype=int, copy=True))
    F, h0, _, _ = find_map_to_target_pauli_sum(H_basis_src, H_basis_tgt)

    nq = ctx.n_qudits
    pauli_weighted = ctx.pauli_weighted

    SG_F = Gate('Symmetry', list(range(nq)), F.T, ctx.pauli_sum.dimensions, np.asarray(h0, dtype=int))
    H_full_tg = pauli_weighted.copy()[pi]
    H_full_F = SG_F.act(pauli_weighted)

    delta = (H_full_tg.phases - H_full_F.phases) % (2 * int(ctx.pauli_sum.lcm))

    if ctx.p == 2 and np.any(delta % 2 != 0):
        F2 = F % 2
        A, B = F2[:nq, :nq], F2[:nq, nq:]
        C, D = F2[nq:, :nq], F2[nq:, nq:]
        hx0 = np.diag((A @ B.T) % 2) % 2
        hz0 = np.diag((C @ D.T) % 2) % 2
        h0_alt = np.concatenate([hx0, hz0]).astype(int)

        SG_F_alt = Gate('Symmetry', list(range(nq)), F.T, ctx.pauli_sum.dimensions, h0_alt)
        H_full_Fa = SG_F_alt.act(pauli_weighted)
        delta_alt = (H_full_tg.phases - H_full_Fa.phases) % 4
        if (delta_alt % 2).sum() < (delta % 2).sum():
            h0, SG_F, H_full_F, delta = h0_alt, SG_F_alt, H_full_Fa, delta_alt

    h_lin = solve_phase_vector_h_from_residual(ctx.base_tableau, delta, ctx.pauli_sum.dimensions,
                                               debug=True, row_basis_cache=ctx.row_basis_cache)
    if h_lin is None:
        return None

    h0_mod = np.asarray(h0, dtype=int) % ctx.two_lcm
    h_lin_mod = np.asarray(h_lin, dtype=int) % ctx.two_lcm
    h_tot = (h0_mod + h_lin_mod) % ctx.two_lcm
    SG = Gate('Symmetry', list(range(nq)), F.T, ctx.pauli_sum.dimensions, h_tot)

    H_out_cf = SG.act(pauli_weighted).to_standard_form()
    H_out_cf.weight_to_phase()

    if not np.array_equal(H_out_cf.tableau, ctx.ref_tableau):
        return None
    if not np.all((ctx.ref_phases - H_out_cf.phases) % ctx.two_lcm == 0):
        return None
    if not np.array_equal(H_out_cf.weights, ctx.ref_weights):
        return None
    Omega = np.zeros((2 * nq, 2 * nq), dtype=int)
    Omega[:nq, nq:] = np.eye(nq, dtype=int)
    Omega[nq:, :nq] = -np.eye(nq, dtype=int)
    if np.all(((ctx.p - 1) * np.diag(F @ Omega @ F.T) + h_tot % 2) != 0):
        return None

    return SG


def clifford_graph_automorphism_search(
    pauli_sum: PauliSum,
    k_wanted: int,
    dynamic_refine_every: int = 0,
    extra_column_invariants: str = "none",
    p2_bitset: str | bool = "auto",
    color_mode: str = "wl",   # "wl" | "coeffs_only" | "none"
    max_wl_rounds: int = 10,
) -> list[Gate]:
    """
    Find up to k automorphisms preserving S and the vector set.
    All preprocessing based only on pauli_sum is handled internally.
    Dynamic WL (if enabled) is used only for ordering every `dynamic_refine_every` steps.
    """

    # ---- preprocessing that depends only on pauli_sum ----
    independent_labels, dependencies = get_linear_dependencies(pauli_sum.tableau, 2)
    labels = sorted(set(independent_labels) | set(dependencies.keys()))
    S_mod = pauli_sum.symplectic_product_matrix()
    G, basis_order = pauli_sum.matroid()
    coeffs = pauli_sum.weights

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
        coeffs=coeffs,
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

    # Prepare references for the leaf checks so they dont have to be computed every time
    pauli_weighted = pauli_sum.copy()
    pauli_weighted.weight_to_phase()
    pauli_standard = pauli_weighted.to_standard_form()
    pauli_standard.weight_to_phase()
    ref_tableau = pauli_standard.tableau.astype(int, copy=False)
    ref_phases = np.asarray(pauli_standard.phases, dtype=int)
    ref_weights = np.asarray(pauli_standard.weights)
    base_tableau = pauli_sum.tableau.astype(int, copy=False)
    base_weights = np.asarray(pauli_sum.weights)
    base_phases = np.asarray(pauli_sum.phases, dtype=int)
    basis_indices = np.asarray(independent_labels, dtype=int)
    basis_source_ps = pauli_sum[basis_indices]

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

    steps = 0
    cur_colors = base_colors.copy()

    def select_next() -> int:
        # MRV measured against base feasibility
        best_i, best_rem = -1, 10**9
        for i in domain_order:
            if phi[i] >= 0:
                continue
            bi = int(base_colors[i])
            rem = sum((not used[y] and (coeffs is None or coeffs[i] == coeffs[y])) for y in base_classes[bi])
            if rem < best_rem:
                best_i, best_rem = i, rem
                if rem <= 1:
                    break
        return best_i

    leaf_ctx = _LeafContext(
        p=p,
        two_lcm=2 * int(pauli_sum.lcm),
        n_qudits=pauli_sum.n_qudits(),
        identity_perm=identity_perm,
        S_mod=S_mod,
        G=G,
        basis_order=basis_order,
        labels=labels,
        pauli_sum=pauli_sum,
        pauli_weighted=pauli_weighted,
        ref_tableau=ref_tableau,
        ref_phases=ref_phases,
        ref_weights=ref_weights,
        base_tableau=base_tableau,
        base_weights=base_weights,
        base_phases=base_phases,
        basis_indices=basis_indices,
        basis_source_ps=basis_source_ps,
        row_basis_cache=row_basis_cache,
    )

    def dynamic_refine():
        nonlocal cur_colors
        if dynamic_refine_every <= 0:
            return
        # 1-WL just to order
        cur_colors = _wl_colors_from_S(S_mod, int(2), coeffs=coeffs, col_invariants=None, max_rounds=1)

    def dfs() -> bool:
        nonlocal steps
        if len(results) >= k_wanted:
            return True
        if np.all(phi >= 0):
            pi = phi.copy()
            leaf = _check_leaf(pi, leaf_ctx)
            if leaf is not None:
                results.append(leaf)
                return True
            return False

        if dynamic_refine_every and (steps % dynamic_refine_every == 0):
            dynamic_refine()
        steps += 1

        i = select_next()
        bi = int(base_colors[i])
        mapped_idx = np.where(phi >= 0)[0].astype(np.int64)

        # Order candidates by current colors (ordering heuristic only)
        candidate = [y for y in base_classes[bi] if not used[y]]
        if coeffs is not None:
            candidate = [y for y in candidate if coeffs[i] == coeffs[y]]
        candidate.sort(key=lambda y: cur_colors[y])

        for y in candidate:
            if not consistency(phi, mapped_idx, i, y):
                continue
            phi[i] = y
            used[y] = True
            if dfs():
                return True
            phi[i] = -1
            used[y] = False
        return False

    dfs()
    return results[:k_wanted]
