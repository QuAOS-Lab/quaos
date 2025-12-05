from __future__ import annotations
from typing import Optional
import numpy as np
import galois
from numba import njit
from sympleq.core.circuits.target import find_map_to_target_pauli_sum
from sympleq.core.finite_field_solvers import solve_linear_system_over_gf
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Gate
from sympleq.core.graphs.graph_coloring import _wl_colors_from_S

Label = int
DepPairs = dict[Label, list[tuple[Label, int]]]


# =============================================================================
# Incremental S-consistency (p-agnostic and p=2 bitset variant)
# =============================================================================


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


# ---- Optional p=2 bitset kernel --------------------------------------------

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


# =============================================================================
# Full checks at a leaf
# =============================================================================


def _check_code_automorphism(
    G: galois.FieldArray,
    basis_order: list[int],
    labels: list[int],
    pi: np.ndarray
) -> bool:
    """
    Linear-code test over GF(p): ∃ U with U G P = G ?
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


# =============================================================================
# Fallback full DFS (feasible by base partition; complete). Kept simple/serial.
# =============================================================================


def _full_dfs_complete(
    pauli_sum: PauliSum,
    independent_labels: list[int],
    S_mod: np.ndarray,
    coeffs: Optional[np.ndarray],
    base_colors: np.ndarray,
    base_classes: dict[int, list[int]],
    G: galois.FieldArray,
    basis_order: list[int],
    labels: list[int],
    k_wanted: int,
    p2_bitset: bool,
    dynamic_refine_every: int = 0,
    F_known_debug: Optional[np.ndarray] = None,
) -> list[Gate]:
    """
    Complete interleaved DFS that maps all labels with feasibility constrained ONLY by base partition.
    Dynamic WL (if enabled) is used only for ordering every `dynamic_refine_every` steps.
    """
    n = S_mod.shape[0]
    results = []

    p = int(pauli_sum.lcm)

    # Prepared consistency kernel
    if p2_bitset:
        bits, _ = _build_bitrows_binary(S_mod)

        def consistent(phi, mapped, i, y):
            return _consistent_bitset(bits, phi, mapped, int(i), int(y))
    else:
        def consistent(phi, mapped, i, y):
            return _consistent_numba(S_mod, phi, mapped, int(i), int(y))

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
            row_basis_cache["gf2"] = _row_basis_indices(base_tableau % 2, 2, base_tableau.shape[1])
        else:
            row_basis_cache["gfp"] = _row_basis_indices(base_tableau % p_uni, p_uni, base_tableau.shape[1])

    phi = -np.ones(n, dtype=np.int64)
    used = np.zeros(n, dtype=bool)

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

    def at_leaf(pi: np.ndarray) -> Gate | None:
        """
        This is where the heavy checks happen - if we are here, we are at a leaf
        (i.e. all mapped). If it returns None, then the permutation is not correct, even though it satisfies al colour
        constraints. This is suitably rare that we can leave all of the expensive checks to this stage.

        :param pi: current permutation
        :return: symmetry circuit or None
        """
        # ---- structural checks ----
        if np.array_equal(pi, np.arange(pauli_sum.n_paulis(), dtype=pi.dtype)):
            return None
        if not np.array_equal(S_mod[np.ix_(pi, pi)], S_mod):
            return None
        if not _check_code_automorphism(G, basis_order, labels, pi):
            return None

        # Phase correction part

        # ---- 1) find F, h0 from the permutation
        H_basis_src = basis_source_ps
        tgt_idx = pi[basis_indices]
        H_basis_tgt = PauliSum.from_tableau(base_tableau[tgt_idx], pauli_sum.dimensions,
                                            weights=base_weights[tgt_idx])
        H_basis_tgt.set_phases(np.array(base_phases[tgt_idx], dtype=int, copy=True))
        F, h0, _, _ = find_map_to_target_pauli_sum(H_basis_src, H_basis_tgt)
        if F_known_debug is not None and np.array_equal(F.T, F_known_debug):
            print('DEBUG: found known F failure means it cant correct h')

        nq = pauli_sum.n_qudits()

        # ---- 2) apply that lift ON THE FULL HAMILTONIAN ----
        SG_F = Gate('Symmetry', list(range(nq)), F.T, pauli_sum.dimensions, np.asarray(h0, dtype=int))
        H_full_tg = pauli_weighted.copy()[pi]
        H_full_F = SG_F.act(pauli_weighted)

        # ---- 3) build residual phases Δ on the FULL SET (mod 2p) ----
        delta = (H_full_tg.phases - H_full_F.phases) % (2 * int(pauli_sum.lcm))

        # Optional qubit lift tweak: if any residual is odd, try canonical h0 from F
        if p == 2 and np.any(delta % 2 != 0):
            # print('Trying alternate h0 for F to reduce odd residuals')
            F2 = F % 2
            A, B = F2[:nq, :nq], F2[:nq, nq:]
            C, D = F2[nq:, :nq], F2[nq:, nq:]
            hx0 = np.diag((A @ B.T) % 2) % 2
            hz0 = np.diag((C @ D.T) % 2) % 2
            h0_alt = np.concatenate([hx0, hz0]).astype(int)

            SG_F_alt = Gate('Symmetry', list(range(nq)), F.T, pauli_sum.dimensions, h0_alt)
            H_full_Fa = SG_F_alt.act(pauli_weighted)
            delta_alt = (H_full_tg.phases - H_full_Fa.phases) % 4
            if (delta_alt % 2).sum() < (delta % 2).sum():
                h0, SG_F, H_full_F, delta = h0_alt, SG_F_alt, H_full_Fa, delta_alt

        # ---- 4) solve for h on the FULL post-F tableau ----
        h_lin = solve_phase_vector_h_from_residual(base_tableau, delta, pauli_sum.dimensions,
                                                   debug=True, row_basis_cache=row_basis_cache)
        if h_lin is None:
            return None

        # ---- 5) compose final symmetry and verify against H (perm-insensitive via standard_form) ----
        two_lcm = 2 * int(pauli_sum.lcm)
        h0_mod = np.asarray(h0, dtype=int) % two_lcm
        h_lin_mod = np.asarray(h_lin, dtype=int) % two_lcm
        h_tot = (h0_mod + h_lin_mod) % two_lcm
        SG = Gate('Symmetry', list(range(nq)), F.T, pauli_sum.dimensions, h_tot)

        H_out_cf = SG.act(pauli_weighted).to_standard_form()
        H_out_cf.weight_to_phase()

        if not np.array_equal(H_out_cf.tableau, ref_tableau):
            return None
        if not np.all((ref_phases - H_out_cf.phases) % two_lcm == 0):
            return None
        if not np.array_equal(H_out_cf.weights, ref_weights):
            return None
        Omega = np.zeros((2 * nq, 2 * nq), dtype=int)
        Omega[:nq, nq:] = np.eye(nq, dtype=int)
        Omega[nq:, :nq] = -np.eye(nq, dtype=int)
        if np.all(((p - 1) * np.diag(F @ Omega @ F.T) + h_tot % 2) != 0):
            # print('Warning: found F,h that dont satisfy Clifford constraint!')
            return None

        return SG

    def maybe_dynamic_refine():
        nonlocal cur_colors
        if dynamic_refine_every <= 0:
            return
        # very light 1-WL just to order (we do not change feasibility!)
        cur_colors = _wl_colors_from_S(S_mod, int(2), coeffs=coeffs, col_invariants=None, max_rounds=1)

    def dfs() -> bool:
        nonlocal steps
        if len(results) >= k_wanted:
            return True
        if np.all(phi >= 0):
            pi = phi.copy()
            leaf = at_leaf(pi)
            if leaf is not None:
                results.append(leaf)
                return True
            return False

        if dynamic_refine_every and (steps % dynamic_refine_every == 0):
            maybe_dynamic_refine()
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
            if not consistent(phi, mapped_idx, i, y):
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

# =============================================================================
# Public API
# =============================================================================


def _gf_solve_one_solution(A_int: np.ndarray, b_int: np.ndarray, p: int) -> Optional[np.ndarray]:
    """Gauss–Jordan elimination over GF(p). Return one solution (free vars=0) or None if inconsistent."""
    GF = galois.GF(p)
    A = GF(A_int % p)
    b = GF(b_int % p)
    N, M = A.shape
    Ab = np.concatenate([A, b.reshape(-1, 1)], axis=1)

    row = 0
    pivots = []
    for col in range(M):
        pivot = None
        for r in range(row, N):
            if Ab[r, col] != GF(0):
                pivot = r
                break
        if pivot is None:
            continue
        if pivot != row:
            tmp = Ab[row, :].copy()
            Ab[row, :], Ab[pivot, :] = Ab[pivot, :], tmp
        inv = GF(1) / Ab[row, col]
        Ab[row, :] = Ab[row, :] * inv
        for r in range(N):
            if r == row:
                continue
            if Ab[r, col] != GF(0):
                Ab[r, :] = Ab[r, :] - Ab[r, col] * Ab[row, :]
        pivots.append(col)
        row += 1
        if row == N:
            break

    # inconsistency: [0 ... 0 | nonzero]
    for r in range(N):
        if np.all(Ab[r, :M] == GF(0)) and Ab[r, M] != GF(0):
            return None

    x = GF.Zeros(M)
    for r, c in enumerate(pivots):
        if r < N:
            x[c] = Ab[r, M]
    return np.array([int(v) for v in x], dtype=int)

# ---------- main solver with diagnostics ----------


def solve_phase_vector_h_from_residual(
    tableau_in: np.ndarray,   # (N, 2n) ints; rows of input Paulis
    delta_2L: np.ndarray,     # (N,) ints mod 2L; desired phase corrections
    dimensions: np.ndarray | list[int],
    debug: bool = False,
    row_basis_cache: dict[str, np.ndarray] | None = None,
) -> Optional[np.ndarray]:
    """
    Solve tableau_in @ h ≡ delta_2L (mod 2L), where h is the additional phase vector.

    Returns h or None if inconsistent.
    """
    A = np.asarray(tableau_in, dtype=int)
    b = np.asarray(delta_2L, dtype=int)
    dims = np.asarray(dimensions, dtype=int)
    if A.ndim != 2:
        raise ValueError("tableau_in must be a 2D array")
    N, M = A.shape
    if M != 2 * dims.size:
        raise ValueError("tableau_in columns must equal 2 * number of qudits")

    L = int(np.lcm.reduce(dims))
    modulus = 2 * L

    if dims.size and np.all(dims == dims[0]):
        p_uni = int(dims[0])
        if p_uni == 2:
            if np.any(b % 2 != 0):
                if debug:
                    print("[phase] qubit fallback: residual has odd entries, cannot fix")
                return None
            rows = None
            if row_basis_cache is not None:
                rows = row_basis_cache.get("gf2")
            if rows is None or rows.size == 0:
                rows = _row_basis_indices(A % 2, 2, A.shape[1])
                if row_basis_cache is not None:
                    row_basis_cache["gf2"] = rows
            if rows.size == 0:
                return None
            A2 = (A[rows] % 2).astype(int, copy=False)
            b2 = ((b[rows] // 2) % 2).astype(int, copy=False)
            sol2 = _gf_solve_one_solution(A2, b2, 2)
            if sol2 is not None:
                return (2 * sol2.astype(int)) % modulus
        else:
            rows = None
            if row_basis_cache is not None:
                rows = row_basis_cache.get("gfp")
            if rows is None or rows.size == 0:
                rows = _row_basis_indices(A % p_uni, p_uni, A.shape[1])
                if row_basis_cache is not None:
                    row_basis_cache["gfp"] = rows
            if rows.size == 0:
                return None
            A_p = (A[rows] % p_uni).astype(int, copy=False)
            b_p = (b[rows] % p_uni).astype(int, copy=False)
            sol_p = _gf_solve_one_solution(A_p, b_p, p_uni)
            if sol_p is not None:
                if debug:
                    print(f"[phase] GF({p_uni}) fallback succeeded")
                return sol_p.astype(int) % modulus

    sol = solve_linear_system_over_gf(A % modulus, b % modulus, modulus)
    if sol is not None:
        if debug:
            print("[phase] direct solve mod", modulus, "succeeded")
        return sol % modulus

    if debug:
        print("[phase] direct mod", modulus, "solve failed")

    return None


def _row_basis_indices(A_int: np.ndarray, p: int, want_cols: int) -> np.ndarray:
    """Select indices of rows giving up to `want_cols` independent equations over GF(p)."""
    GF = galois.GF(p)
    A = GF(A_int % p).copy()
    N, M = A.shape
    used = np.zeros(N, dtype=bool)
    basis = []
    col = 0
    for _ in range(N):
        if col >= M:
            break
        pick = None
        for r in range(N):
            if used[r]:
                continue
            if A[r, col] != GF(0):
                pick = r
                break
        if pick is None:
            col += 1
            continue
        basis.append(pick)
        used[pick] = True
        inv = GF(1) / A[pick, col]
        A[pick, :] = A[pick, :] * inv
        for r in range(N):
            if r == pick or used[r]:
                continue
            if A[r, col] != GF(0):
                A[r, :] = A[r, :] - A[r, col] * A[pick, :]
        col += 1
        if len(basis) >= want_cols:
            break
    return np.array(basis, dtype=int)
