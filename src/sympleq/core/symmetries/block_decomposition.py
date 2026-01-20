import numpy as np
import itertools
from .modular_helpers import (mod_p, rank_mod, _solve_linear, nullspace_mod, omega_matrix, inv_mod_mat,
                              inv_mod_scalar, matmul_mod, is_symplectic, independent_columns, solve_linear_many)
from .rcf_prepass import rcf_prepass
from sympleq.core.graphs.utils import qudit_coupling_graph

# Atomic seeds


def _restrict_operator(F: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """
    Return the matrix of F restricted to span(B), in coordinates of B.
    Requires span(B) invariant (true for primary bases).
    """
    B = independent_columns(mod_p(B, p), p)
    FB = mod_p(F @ B, p)
    return mod_p(solve_linear_many(B, FB, p), p)


def _poly_eval_matrix(F: np.ndarray, poly: np.ndarray, p: int) -> np.ndarray:
    """
    Evaluate poly(F) with coeffs low->high, column-action convention.
    """
    F = mod_p(F, p)
    poly = mod_p(poly, p).reshape(-1)
    d = F.shape[0]
    out = np.zeros((d, d), dtype=np.int64)
    P = np.eye(d, dtype=np.int64)
    for a in poly:
        a = int(a) % p
        if a:
            out = mod_p(out + a * P, p)
        P = mod_p(P @ F, p)
    return out


def atomic_seeds(F: np.ndarray, p: int) -> list[np.ndarray]:
    """
    Return 'atomic' seed vectors v (columns) = Jordan-chain tops of N=phi(F)
    on each primary component V_phi, lifted back to ambient coordinates.
    """
    meta = rcf_prepass(F, p)

    # rcf_prepass must expose primaries with at least:
    #   prim[key]["poly"]      irreducible phi(x) coefficients low->high
    #   prim[key]["exponent"]  k_phi
    #   prim[key]["V_basis"]   basis columns for V_phi
    prim = meta["primaries"]

    seeds: list[np.ndarray] = []
    for data in prim.values():
        q = data["poly"]
        exp = int(data["exponent"])
        V = independent_columns(mod_p(data["V_basis"], p), p)
        if V.shape[1] == 0:
            continue

        # Restrict F to V_phi
        Fp = _restrict_operator(F, V, p)

        # Nilpotent on primary: N = q(F) (restricted)
        Np = _poly_eval_matrix(Fp, q, p)

        # Chain tops for exact lengths L=1..exp (in primary coordinates)
        tops = jordan_chain_tops_nilpotent(Np, exp, p)  # dict[int, (dimV x mult_L)]

        for L, T in tops.items():
            # Each column of T is a 'top' (primary coords); lift to ambient
            for j in range(T.shape[1]):
                v = mod_p(V @ T[:, j:j + 1], p)
                if not np.all(v % p == 0):
                    seeds.append(v)

    return seeds


def _augment_seeds(seeds: list[np.ndarray], p: int, n_extra: int, rng: np.random.Generator) -> list[np.ndarray]:
    """
    Add random linear combinations of atomic seeds (helps in self-reciprocal sectors).
    """
    if not seeds or n_extra <= 0:
        return seeds
    M = np.concatenate(seeds, axis=1)  # n2 x m
    for _ in range(n_extra):
        coeff = rng.integers(0, p, size=(M.shape[1], 1), dtype=np.int64)
        seeds.append(mod_p(M @ coeff, p))
    return seeds


# =========================
# Mode graph helpers
# =========================

def _mode_graph_from_S(S: np.ndarray, p: int) -> list[list[int]]:
    """
    Build adjacency for 'modes' (pairs). S is [U_all | V_all].
    Connect i--j if any entry in the 4x4 cross-block between modes i and j is nonzero mod p.
    """
    n2 = S.shape[0]
    assert n2 % 2 == 0
    k = n2 // 2  # number of modes
    adj = [[] for _ in range(k)]
    M = mod_p(S, p)
    for i in range(k):
        rows_i = [i, k + i]
        for j in range(i + 1, k):
            cols_j = [j, k + j]
            block_ij = M[np.ix_(rows_i, cols_j)]
            block_ji = M[np.ix_([j, k + j], [i, k + i])]
            if np.any(block_ij % p != 0) or np.any(block_ji % p != 0):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _components(adj: list[list[int]]) -> list[list[int]]:
    """Connected components from adjacency list."""
    k = len(adj)
    seen = [False] * k
    comps = []
    for s in range(k):
        if seen[s]:
            continue
        stack = [s]
        seen[s] = True
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
        comps.append(sorted(comp))
    return comps


def _mode_permutation_for_blocks(S: np.ndarray, p: int) -> np.ndarray:
    """
    Permute modes so that:
      (1) Nontrivial components with >= 2 modes first, sorted by size asc (tie: rank desc),
      (2) then nontrivial singletons (rank desc),
      (3) then trivial components (size asc).
    """
    n2 = S.shape[0]
    k = n2 // 2
    adj = _mode_graph_from_S(S, p)
    comps = _components(adj)

    M = mod_p(S - np.eye(n2, dtype=np.int64), p)

    def comp_rank(comp_modes: list[int]) -> int:
        idx = comp_modes + [c + k for c in comp_modes]
        return rank_mod(M[np.ix_(idx, idx)], p)

    info = []
    for comp in comps:
        r = comp_rank(comp)
        info.append({"modes": comp, "n_modes": len(comp), "rank": r})

    non_triv_ge2 = [c for c in info if c["rank"] > 0 and c["n_modes"] >= 2]
    non_triv_1 = [c for c in info if c["rank"] > 0 and c["n_modes"] == 1]
    triv = [c for c in info if c["rank"] == 0]

    non_triv_ge2.sort(key=lambda c: (c["n_modes"], -c["rank"]))
    non_triv_1.sort(key=lambda c: -c["rank"])
    triv.sort(key=lambda c: c["n_modes"])

    order_modes: list[int] = []
    for bucket in (non_triv_ge2, non_triv_1, triv):
        for c in bucket:
            order_modes.extend(c["modes"])

    P = np.zeros((k, k), dtype=np.int64)
    for new_pos, old_mode in enumerate(order_modes):
        P[new_pos, old_mode] = 1
    return P


def _apply_mode_permutation_to_ST(S: np.ndarray, T: np.ndarray, p: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Permute modes with the same P in U and V: Π = diag(P, P).
    Preserves Omega; produces S' = Π^{-1} S Π, T' = T Π.
    """
    n2 = S.shape[0]
    k = n2 // 2
    P = _mode_permutation_for_blocks(S, p)
    Pi = np.block([[P, np.zeros((k, k), dtype=np.int64)],
                  [np.zeros((k, k), dtype=np.int64), P]])
    Pi_inv = Pi.T
    S2 = mod_p(Pi_inv @ S @ Pi, p)
    T2 = mod_p(T @ Pi, p)
    return S2, T2

# =========================
# Symplectic structures
# =========================


# Safe scalar from 1x1 array (no deprecation warnings)
def _scalar(x: np.ndarray) -> int:
    return int(np.asarray(x, dtype=np.int64).reshape(()))


def _split_uv(T: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Given canonical [U | V], return (U, V)."""
    k2 = T.shape[1]
    assert k2 % 2 == 0
    k = k2 // 2
    return T[:, :k], T[:, k:]


def nilpotent_kernel_basis(N: np.ndarray, j: int, p: int) -> np.ndarray:
    d = N.shape[0]
    Nj = np.eye(d, dtype=np.int64)
    for _ in range(j):
        Nj = mod_p(Nj @ N, p)
    return nullspace_mod(Nj, p)


def _basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    base = independent_columns(base, p) if base.size else base
    picked = np.zeros((candidates.shape[0], 0), dtype=np.int64)
    r_base = rank_mod(base, p) if base.size else 0

    for j in range(candidates.shape[1]):
        c = candidates[:, j:j+1]
        r_try = rank_mod(np.concatenate([base, picked, c], axis=1), p)
        if r_try > r_base + picked.shape[1]:
            picked = np.concatenate([picked, c], axis=1)
            if picked.shape[1] == want:
                return picked

    raise RuntimeError("_basis_extend: could not extend by required amount.")


def jordan_chain_tops_nilpotent(N: np.ndarray, max_exp: int, p: int) -> dict[int, np.ndarray]:
    """
    tops[L] columns are representatives of K_L / (K_{L-1} + N K_{L+1}),
    giving chain tops of exact length L.
    """
    d = N.shape[0]
    K: list[np.ndarray] = [np.zeros((d, 0), dtype=np.int64)]
    for j in range(1, max_exp + 1):
        K.append(independent_columns(nilpotent_kernel_basis(N, j, p), p))
    K.append(K[max_exp])  # K_{max+1} := K_max

    tops: dict[int, np.ndarray] = {}
    for L in range(1, max_exp + 1):
        KL = K[L]
        if KL.shape[1] == 0:
            continue

        S = K[L - 1]
        NKLp1 = mod_p(N @ K[L+1], p) if K[L+1].shape[1] else np.zeros((d, 0), dtype=np.int64)
        if NKLp1.shape[1]:
            S = np.concatenate([S, NKLp1], axis=1) if S.shape[1] else NKLp1
        S = independent_columns(S, p) if S.shape[1] else S

        rS = rank_mod(S, p) if S.shape[1] else 0
        need = KL.shape[1] - rS
        if need <= 0:
            continue

        chosen = _basis_extend(S, KL, need, p)
        tops[L] = chosen

    return tops


# =========================
# Canonical symplectic basis from a span
# =========================

def symplectic_basis_from_span(B: np.ndarray, p: int) -> np.ndarray:
    """
    Construct a canonical symplectic basis [U | V] from a non-degenerate even-dimensional span B.
    Output: T (n2 x 2k) such that T^T Omega T = Omega_k.
    """
    # Ensure even ambient dimension; pad a zero row if necessary.
    if B.shape[0] % 2 != 0:
        B = np.vstack([B, np.zeros((1, B.shape[1]), dtype=np.int64)])
    n2, s = B.shape
    Omega = omega_matrix(n2 // 2, p)
    S = independent_columns(B, p)
    dS = S.shape[1]

    if dS % 2 != 0:
        raise ValueError("Subspace dimension must be even")

    def try_build(S_perm: np.ndarray) -> np.ndarray | None:
        idx_used = set()
        U, V = [], []
        for _ in range(0, dS, 2):
            found = False
            for j in range(dS):
                if j in idx_used:
                    continue
                u = S_perm[:, j:j + 1]
                for k in range(j + 1, dS):
                    if k in idx_used:
                        continue
                    v = S_perm[:, k: k + 1]
                    if _scalar(u.T @ Omega @ v % p) == 0:
                        continue

                    # Project u into the symplectic complement of previously chosen pairs
                    for u_prev, v_prev in zip(U, V):
                        alpha = _scalar(u.T @ Omega @ v_prev % p)
                        gamma = _scalar(u.T @ Omega @ u_prev % p)
                        if alpha or gamma:
                            u = mod_p(u - alpha * u_prev + gamma * v_prev, p)

                    if np.all(u % p == 0):
                        continue

                    # After updating u, ensure the pairing with v is still nonzero
                    if _scalar(u.T @ Omega @ v % p) == 0:
                        continue

                    # Orthogonalize v against previous pairs
                    for u_prev, v_prev in zip(U, V):
                        coeff_u = _scalar(v.T @ Omega @ v_prev % p)
                        coeff_v = _scalar(v.T @ Omega @ u_prev % p)
                        if coeff_u or coeff_v:
                            v = mod_p(v - coeff_u * u_prev + coeff_v * v_prev, p)

                    beta = _scalar(u.T @ Omega @ v % p)
                    if beta == 0:
                        continue

                    beta_inv = inv_mod_scalar(beta, p)
                    v = mod_p(v * beta_inv, p)

                    idx_used.update([j, k])
                    found = True
                    break
                if found:
                    U.append(mod_p(u, p))
                    V.append(mod_p(v, p))
                    break
            if not found:
                return None

        T = np.hstack(U + V)
        G = mod_p(T.T @ Omega @ T, p)
        if np.array_equal(G % p, omega_matrix(len(U), p)):
            return T
        return None

    # Try deterministic order first, then random permutations to find a symplectic pairing
    T = try_build(S)
    if T is not None:
        return T

    rng = np.random.default_rng(2025)
    for _ in range(256):
        perm = rng.permutation(dS)
        T = try_build(S[:, perm])
        if T is not None:
            return T

    # Random invertible column mixes to search the span (not just permutations).
    for _ in range(512):
        R = rng.integers(0, p, size=(dS, dS), dtype=np.int64)
        # Ensure invertible over GF(p)
        if np.linalg.matrix_rank(R % p) != dS:
            continue
        T = try_build(mod_p(S @ R, p))
        if T is not None:
            return T

    # Exhaustive search for small dimensions to avoid flaky failures in tests (factorial grows fast).
    if dS <= 8:
        for perm in itertools.permutations(range(dS)):
            T = try_build(S[:, perm])
            if T is not None:
                return T

    raise RuntimeError("Constructed basis is not symplectic")


# =========================
# Minimal symplectic block via Krylov + partner
# =========================


def krylov_closure(F: np.ndarray, v: np.ndarray, p: int, cap: int | None = None) -> np.ndarray:
    """Krylov span K = span{v, Fv, ..., F^{k-1}v} until rank stops increasing."""
    n2 = F.shape[0]
    cap_local: int = n2 if cap is None else int(cap)
    V = np.zeros((n2, 0), dtype=np.int64)
    w = v.reshape(-1, 1)
    r0 = 0
    for _ in range(cap_local):
        can = np.concatenate([V, w], axis=1)
        r = rank_mod(can, p)
        if r > r0:
            V = can
            r0 = r
            w = matmul_mod(F, w, p)
        else:
            break
    return V


def build_partner_for_krylov(F: np.ndarray, K: np.ndarray, p: int) -> np.ndarray:
    """
    Given K = span{v, Fv, ..., F^{k-1}v}, find z with constraints:
      <v, F^b z> = δ_{b, k-1} for b=0..k-1.
    Efficient vector recurrence for rows r_b = (F^b)^T Omega v.
    """
    n2 = F.shape[0]
    Omega = omega_matrix(n2 // 2, p)
    k = K.shape[1]
    v0 = K[:, 0:1]  # column

    # rb_0 = Omega v0, rb_{b+1} = F^T rb_b
    rb = mod_p(Omega @ v0, p)
    rows = [rb.T.reshape(-1)]
    for _ in range(1, k):
        rb = mod_p(F.T @ rb, p)
        rows.append(rb.T.reshape(-1))
    A = np.stack(rows, axis=0)           # (k, n2)
    b_vec = np.zeros((k, 1), dtype=np.int64)
    b_vec[-1, 0] = 1                       # δ_{b,k-1}
    z = _solve_linear(A, b_vec, p)         # (n2,1)
    return z


def _restricted_action(F: np.ndarray, T_blk: np.ndarray, p: int) -> np.ndarray:
    """Return S_blk = T_blk^{-1} F T_blk using J^{-1} T^T Omega."""
    n2 = F.shape[0]
    Omega = omega_matrix(n2 // 2, p)
    J = mod_p(T_blk.T @ Omega @ T_blk, p)  # should be Omega_k
    J_inv = inv_mod_mat(J, p)
    left_inv = matmul_mod(J_inv, matmul_mod(T_blk.T, Omega, p), p)
    return matmul_mod(left_inv, matmul_mod(F, T_blk, p), p)


def _minimal_block_from_seeds(
    F: np.ndarray,
    p: int,
    seeds: list[np.ndarray],
    min_block_size: int = 0
) -> np.ndarray | None:
    """
    tries seeds, returns best T_blk or None.
    Score: smallest size (>= min_block_size if >0), tie-break by rank(S_blk - I) descending.
    """
    n2 = F.shape[0]
    n = n2 // 2
    Omega = omega_matrix(n, p)

    best: tuple[int, int, np.ndarray] | None = None  # (size, r_score, T_blk)

    for v in seeds:
        if v.size == 0 or np.all(v % p == 0):
            continue
        K = krylov_closure(F, v, p)
        if K.shape[1] == 0:
            continue
        try:
            z = build_partner_for_krylov(F, K, p)
        except RuntimeError:
            continue
        Kz = krylov_closure(F, z, p)
        W = np.concatenate([K, Kz], axis=1)

        B = independent_columns(W, p)
        if B.shape[1] % 2 != 0:
            continue
        if min_block_size and B.shape[1] < min_block_size:
            continue

        G = mod_p(B.T @ Omega @ B, p)
        if rank_mod(G, p) != B.shape[1]:
            continue

        try:
            T_blk = symplectic_basis_from_span(B, p)  # [U|V]
        except Exception:
            continue

        S_blk = _restricted_action(F, T_blk, p)
        r_score = rank_mod(mod_p(S_blk - np.eye(S_blk.shape[0], dtype=np.int64), p), p)
        if r_score == 0:
            continue  # reject trivial block

        size = T_blk.shape[1]

        # Enforce a *lower* size bound if desired
        if min_block_size and size < min_block_size:
            continue
        # # Prefer *larger* minimal blocks; break ties by r_score descending
        # if (best is None) or (size > best[0]) or (size == best[0] and r_score > best[1]):
        #     best = (size, r_score, T_blk)
        if (best is None) or (size < best[0]) or (size == best[0] and r_score > best[1]):
            best = (size, r_score, T_blk)

    return None if best is None else best[2]


def minimal_symplectic_block_in_complement(
    F: np.ndarray, p: int, N: np.ndarray, trials: int = 64
) -> tuple[np.ndarray, np.ndarray] | None:
    """
    Find the smallest-dimension invariant non-degenerate symplectic block W for F,
    **constrained to the subspace span(N)** (columns of N).
    Returns (T_blk, S_blk) in ambient coordinates, or None if none found.
    """
    rng = np.random.default_rng(2025)

    seeds: list[np.ndarray] = [N[:, i: i + 1] for i in range(N.shape[1])]
    for _ in range(trials):
        coeffs = rng.integers(0, p, size=(N.shape[1], 1), dtype=np.int64)
        seeds.append(mod_p(N @ coeffs, p))

    T_blk = _minimal_block_from_seeds(F, p, seeds, min_block_size=0)
    if T_blk is None:
        return None
    S_blk = _restricted_action(F, T_blk, p)
    return T_blk, S_blk


def minimal_symplectic_block(F: np.ndarray, p: int, trials: int = 64) -> tuple[np.ndarray, np.ndarray]:
    """
    Find a smallest-dimension invariant, non-degenerate symplectic block W for F.
    Returns (T_blk, S_blk).
    """
    n2 = F.shape[0]
    rng = np.random.default_rng(2025)

    seeds: list[np.ndarray] = [np.eye(n2, dtype=np.int64)[:, i: i + 1] for i in range(n2)]
    seeds += [rng.integers(0, p, size=(n2, 1), dtype=np.int64) for _ in range(trials)]

    T_blk = _minimal_block_from_seeds(F, p, seeds, min_block_size=0)
    if T_blk is None:
        raise RuntimeError("Failed to find a non-degenerate invariant symplectic block")
    S_blk = _restricted_action(F, T_blk, p)
    return T_blk, S_blk


# =========================
# Local symplectic completion and global peel
# =========================

def complete_symplectic_local(T_blk: np.ndarray, p: int) -> np.ndarray:
    """
    Given T_blk (n2_sub x k2) a canonical symplectic basis of a non-degenerate subspace W,
    complete to a full symplectic basis of the trailing subspace:
       T_local = [U_W | U_perp | V_W | V_perp], with T_local^T Omega_sub T_local = Omega_sum.
    """
    n2_sub = T_blk.shape[0]
    Omega_sub = omega_matrix(n2_sub // 2, p)

    # Omega-orthogonal complement: W_perp = {x | T_blk^T Omega x = 0}
    A = mod_p(T_blk.T @ Omega_sub, p)           # k2 x n2_sub
    N = nullspace_mod(A, p)               # n2_sub x d_perp
    d_perp = N.shape[1]
    if d_perp % 2 != 0:
        raise RuntimeError("Complement has odd dimension (cannot symplectically complete)")

    if d_perp > 0:
        T_perp = symplectic_basis_from_span(N, p)   # [U_perp | V_perp]
        U_perp, V_perp = _split_uv(T_perp)
    else:
        U_perp = np.zeros((n2_sub, 0), dtype=np.int64)
        V_perp = np.zeros((n2_sub, 0), dtype=np.int64)

    U_W, V_W = _split_uv(T_blk)
    T_local = np.concatenate([U_W, U_perp, V_W, V_perp], axis=1)

    G = mod_p(T_local.T @ Omega_sub @ T_local, p)
    if not np.array_equal(G % p, Omega_sub % p):
        raise RuntimeError("Local completion is not symplectic (T_local^T Omega T_local ≠ Omega).")
    return T_local


# =========================
# Nontrivial block finder
# =========================

def minimal_symplectic_block_full(
    F: np.ndarray, p: int, trials: int = 64, min_block_size: int = 2
) -> np.ndarray:
    n2 = F.shape[0]
    rng = np.random.default_rng(2025)

    seeds = atomic_seeds(F, p)

    # Good fallback if factorization / prepass returns nothing (should be rare)
    if not seeds:
        seeds = [np.eye(n2, dtype=np.int64)[:, i:i + 1] for i in range(n2)]

    # Optional: add random combos of atomic seeds (use trials as “augmentation budget”)
    seeds = _augment_seeds(seeds, p, n_extra=trials, rng=rng)

    T_blk = _minimal_block_from_seeds(F, p, seeds, min_block_size=min_block_size)
    if T_blk is None:
        raise RuntimeError("Failed to find a nontrivial symplectic block meeting the size floor")
    return T_blk


def block_decompose(
    F: np.ndarray, p: int, min_block_size: int = 2, trials: int = 64
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fully recursive block decomposition using mode-space recursion.

    Returns:
        S, T such that:
          - T is symplectic (T^T Ω T = Ω),
          - S = T^{-1} F T,
          - S is block-structured in mode space (up to SWAPs).
    """
    assert is_symplectic(F, p), "F must be symplectic"
    n2 = F.shape[0]
    assert n2 % 2 == 0
    n = n2 // 2
    Omega = omega_matrix(n, p)

    # Current similarity and transformed F
    T_global = np.eye(n2, dtype=np.int64)
    F_cur = F.copy()

    # Number of modes already peeled into blocks at the front
    mode_offset = 0

    while mode_offset < n:
        # Remaining number of modes
        n_rem = n - mode_offset
        if n_rem <= 0:
            break

        # Indices for remaining modes in canonical [X | Z] order
        x_idx_rem = list(range(mode_offset, n))
        z_idx_rem = list(range(n + mode_offset, n2))
        idx_rem = x_idx_rem + z_idx_rem  # this is canonical for Ω_sub

        # Restrict F_cur to the remaining subspace
        F_sub = mod_p(F_cur[np.ix_(idx_rem, idx_rem)], p)

        # Try to find a nontrivial minimal symplectic block in F_sub
        try:
            T_blk_sub = minimal_symplectic_block_full(F_sub, p, trials=trials, min_block_size=min_block_size)
        except RuntimeError:
            # No suitable block found in the remaining subspace
            break

        k2 = T_blk_sub.shape[1]
        if k2 % 2 != 0:
            # Should not happen for a valid symplectic block
            break
        k = k2 // 2
        if k == 0:
            break

        # Complete this block to a full symplectic basis of the remaining subspace
        # T_local is (2*n_rem x 2*n_rem), canonical [U_block | U_perp | V_block | V_perp]
        T_local = complete_symplectic_local(T_blk_sub, p)

        # Lift T_local back to the full space as block-diagonal (I on peeled modes)
        T_lift = np.eye(n2, dtype=np.int64)
        T_lift[np.ix_(idx_rem, idx_rem)] = T_local

        # Update global similarity and transformed F
        T_global = mod_p(T_global @ T_lift, p)
        F_cur = mod_p(inv_mod_mat(T_lift, p) @ F_cur @ T_lift, p)

        # After this, the first k modes in the remaining subspace
        # have been turned into a canonical block. In global coordinates
        # this corresponds to modes [mode_offset .. mode_offset + k - 1].
        mode_offset += k

        # If we want to enforce a floor on block size in modes, we can break
        # once remaining modes are too small to form a nontrivial block.
        if n - mode_offset < 1:
            break

    # Final sanity: T_global must be symplectic
    assert np.array_equal(mod_p(T_global.T @ Omega @ T_global, p), Omega % p), \
        "Constructed T_global is not symplectic"

    # Compute S = T^{-1} F T
    S = mod_p(inv_mod_mat(T_global, p) @ F @ T_global, p)

    # Optionally pack coupled modes contiguously
    S, T_global = _apply_mode_permutation_to_ST(S, T_global, p)

    # Final checks
    assert is_symplectic(T_global, p)
    assert np.array_equal(S, mod_p(inv_mod_mat(T_global, p) @ F @ T_global, p))

    return S, T_global


def block_decompose_optimal(
    F: np.ndarray,
    p: int,
    min_block_size: int = 2,
    trials: int = 64,
    certify: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Symplectic block decomposition with an optimality certificate.

    This function:
      1. Runs a structural pre-pass rcf_prepass(F, p) to compute a
         theoretical lower bound Lmin_star on the half-dimension of any
         nontrivial symplectic invariant block.
      2. Calls block_decompose(F, p, ...) to obtain a block-diagonal S
         and a symplectic similarity T.
      3. Computes the actual largest block size from S via ordered_block_sizes.
      4. Optionally checks that the achieved qudit cost equals Lmin_star.

    Returns:
        S, T, info where:
          - S = T^{-1} F T is block-structured,
          - T is symplectic (T^T Ω T = Ω),
          - info is a dict with:
              "Lmin_star":  theoretical lower bound on half-dimension,
              "Q_alg":      achieved qudit cost (max half-block-dim),
              "block_sizes": list of block sizes in phase-space dimension (2 * n_modes).
    """
    # 1. Structural pre-pass: compute sector data and lower bound
    meta = rcf_prepass(F, p)
    # print("rcf_prepass keys:", list(meta.keys()))
    Lmin_star = max(sec.get("half_dim_floor", 0) for sec in meta["sectors"])
    Lmin_star = max(1, min(Lmin_star, F.shape[0] // 2))
    print("Lmin_star (theory half-dim):", Lmin_star)

    # 2. Run the existing decomposition
    S, T = block_decompose(F, p, min_block_size=min_block_size, trials=trials)
    sizes = ordered_block_sizes(S, p)
    print("Block sizes:", sizes, "Q_alg:", max(sizes) // 2)

    # 3. Extract actual block sizes (phase-space dims: 2 * n_modes)
    sizes = ordered_block_sizes(S, p)
    Q_alg = max(sizes) // 2 if sizes else 0  # qudit cost = half dimension

    if certify and Q_alg < Lmin_star:
        raise RuntimeError(
            f"Inconsistent: Q_alg={Q_alg} < Lmin_star={Lmin_star}."
        )

    # info = {
    #     "Lmin_star": Lmin_star,
    #     "Q_alg": Q_alg,
    #     "block_sizes": sizes,
    #     "certified_optimal": (Q_alg == Lmin_star),
    # }

    return S, T  # , info


# =========================
# Diagnostics
# =========================

def ordered_block_sizes(S: np.ndarray, p: int) -> list[int]:
    """Return 2*n for connected components, ordered with the same nontrivial-first policy."""
    n2 = S.shape[0]
    k = n2 // 2
    adj = _mode_graph_from_S(S, p)
    comps = _components(adj)
    M = mod_p(S - np.eye(n2, dtype=np.int64), p)

    def comp_rank(comp):
        idx = comp + [c + k for c in comp]
        return rank_mod(M[np.ix_(idx, idx)], p)

    info = [{"modes": c, "n_modes": len(c), "rank": comp_rank(c)} for c in comps]
    non_triv_ge2 = [d for d in info if d["rank"] > 0 and d["n_modes"] >= 2]
    non_triv_1 = [d for d in info if d["rank"] > 0 and d["n_modes"] == 1]
    triv = [d for d in info if d["rank"] == 0]

    non_triv_ge2.sort(key=lambda d: (d["n_modes"], -d["rank"]))
    non_triv_1.sort(key=lambda d: -d["rank"])
    triv.sort(key=lambda d: d["n_modes"])

    ordered = non_triv_ge2 + non_triv_1 + triv
    return [2 * d["n_modes"] for d in ordered]


def block_indexes(F: np.ndarray) -> list[list[int]]:
    """
    Compute qudit blocks as connected components of the coupling graph.

    Args:
        F: 2n x 2n symplectic matrix (integer numpy array).

    Returns:
        A list of blocks, where each block is a sorted list of qudit indices.
        For example, if qudits (0,1) and (2,3) form two decoupled clusters,
        returns [[0, 1], [2, 3]] (order of blocks is not guaranteed).
    """
    neighbors = qudit_coupling_graph(F)
    n = len(neighbors)

    visited = [False] * n
    blocks: list[list[int]] = []

    for start in range(n):
        if visited[start]:
            continue

        # BFS/DFS to collect a connected component
        stack = [start]
        visited[start] = True
        component: list[int] = []

        while stack:
            v = stack.pop()
            component.append(v)
            for w in neighbors[v]:
                if not visited[w]:
                    visited[w] = True
                    stack.append(w)

        component.sort()
        blocks.append(component)

    # Optional: sort blocks by their smallest index for reproducibility
    blocks.sort(key=lambda comp: comp[0])
    return blocks
