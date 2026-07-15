# sympleq/core/symmetries/atomic_decomposition_helpers/atomic_paired.py
from __future__ import annotations

import numpy as np
from typing import Any, Dict, List, Tuple

from ..modular_helpers import (
    mod_p,
    independent_columns,
    omega_matrix,
    rank_mod,
    nullspace_mod,
    mat_pow_mod,
    solve_linear,
)
from .atomic_types import (AtomicBlock, AtomicInvariant,
                           ExtractionObstruction, SearchBudgetExceeded)
from .atomic_linear import (
    is_nondegenerate,
    darboux_basis_from_span,
    symplectic_orthogonal_complement_in_span,
)
from .module_invariants import (
    restrict_operator_invariant,
    q_of_F_restricted,
    cyclic_submodule_basis,
)
from .atomic_krylov import select_module_generators_from_top_quotient
from .atomic_filtration import build_nilpotent_filtration

# Single shared matrix-power implementation (B7); private name kept for callers.
_mat_pow_mod = mat_pow_mod


def _pairing_matrix_between(V_left: np.ndarray, V_right: np.ndarray, p: int) -> np.ndarray:
    """P = V_left^T Ω V_right in GF(p)."""
    n2 = V_left.shape[0]
    Omega = omega_matrix(n2 // 2, p)
    return mod_p(V_left.T @ Omega @ V_right, p)


def _intersection_basis(A: np.ndarray, B: np.ndarray, p: int) -> np.ndarray:
    """Return an independent column basis for span(A) ∩ span(B) over GF(p)."""
    A = independent_columns(mod_p(A, p), p)
    B = independent_columns(mod_p(B, p), p)
    if A.size == 0 or B.size == 0:
        return np.zeros((A.shape[0], 0), dtype=np.int64)

    # Solve A x = B y  <=>  [A | -B] [x;y] = 0
    M = np.concatenate([A, mod_p(-B, p)], axis=1)
    N = nullspace_mod(M, p)  # (a+b)×k
    if N.size == 0:
        return np.zeros((A.shape[0], 0), dtype=np.int64)
    X = N[: A.shape[1], :]
    Id = mod_p(A @ X, p)
    return independent_columns(Id, p)


def _cyclic_module_has_full_rank(
    Fp: np.ndarray,
    Np: np.ndarray,
    v: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
) -> bool:
    """
    True iff cyclic_submodule_basis produces rank exactly deg_q*L (after independent filtering).
    """
    v = mod_p(v.reshape(-1, 1), p)
    target = int(deg_q) * int(L)
    try:
        C = cyclic_submodule_basis(Fp, Np, v, int(deg_q), int(L), p)
    except RuntimeError:
        return False
    C = independent_columns(mod_p(C, p), p)
    return int(C.shape[1]) == target


def _candidate_generators_from_top_space(
    Fp: np.ndarray,
    Np: np.ndarray,
    top_basis: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
    *,
    denom: np.ndarray | None = None,
) -> np.ndarray:
    """
    Deterministic module-generator pool drawn from a top quotient.

    This uses quotient F-orbits rather than individual GF(p) top-basis columns,
    so a deg(q)-dimensional base-field orbit contributes one generator.
    """
    return select_module_generators_from_top_quotient(
        Fp, Np, top_basis, int(deg_q), int(L), int(p), denom=denom
    )


def _select_right_generators_with_full_pairing(
    *,
    NA: np.ndarray,                # (dimVq × m)
    P: np.ndarray,                 # (dimVq × dimVqs)
    right_pool: np.ndarray,        # (dimVqs × t)
    p: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Deterministically pick m columns W from right_pool such that M = NA^T P W is invertible.
    Returns (W, M).
    """
    NA = mod_p(NA, p)
    P = mod_p(P, p)
    right_pool = independent_columns(mod_p(right_pool, p), p)

    m = NA.shape[1]
    if m == 0:
        return np.zeros((right_pool.shape[0], 0), dtype=np.int64), np.zeros((0, 0), dtype=np.int64)
    if right_pool.shape[1] == 0:
        raise RuntimeError("Paired sector: empty right generator pool; cannot form dual pairing.")

    W_list: List[np.ndarray] = []
    M_cols: List[np.ndarray] = []
    rM = 0

    for j in range(right_pool.shape[1]):
        w = right_pool[:, j:j + 1]
        col = mod_p(NA.T @ (P @ w), p)  # (m×1)

        # rank test on the growing pairing matrix
        if not M_cols:
            M_try = col
        else:
            M_try = np.concatenate(M_cols + [col], axis=1)

        r_try = rank_mod(M_try, p)
        if r_try > rM:
            W_list.append(w)
            M_cols.append(col)
            rM = r_try
            if len(W_list) == m:
                break

    if len(W_list) != m:
        raise RuntimeError(
            f"Paired sector: could not select {m} right generators with full-rank pairing; "
            f"got {len(W_list)} (rank={rM})."
        )

    W = mod_p(np.concatenate(W_list, axis=1), p)          # (dimVqs × m)
    M = mod_p(np.concatenate(M_cols, axis=1), p)          # (m × m)
    if rank_mod(M, p) != m:
        raise RuntimeError("Paired sector: internal error, selected pairing matrix is not full rank.")
    return W, M




def _ordered_orbit_columns(Fp: np.ndarray, reps: np.ndarray, deg_q: int, p: int) -> np.ndarray:
    """
    Return the ordered orbit matrix

        [v_0, F v_0, ..., F^{r-1} v_0 | v_1, ..., F^{r-1} v_1 | ...]

    for chain-head representatives stored as columns of ``reps``.  The order is
    intentionally not reduced by ``independent_columns`` because callers use the
    row/column labels (representative index, orbit component) to build dual
    pairings.
    """
    Fp = mod_p(Fp, p)
    reps = mod_p(reps, p)
    deg_q = int(deg_q)
    if reps.shape[1] == 0:
        return np.zeros((Fp.shape[0], 0), dtype=np.int64)

    F_pows = [np.eye(Fp.shape[0], dtype=np.int64)]
    for _ in range(1, deg_q):
        F_pows.append(mod_p(F_pows[-1] @ Fp, p))

    cols: list[np.ndarray] = []
    for j in range(reps.shape[1]):
        v = reps[:, j:j + 1]
        for Ppow in F_pows:
            cols.append(mod_p(Ppow @ v, p))
    return mod_p(np.concatenate(cols, axis=1), p)


def _orbit_pairing_matrix(
    *,
    F_left: np.ndarray,
    N_left: np.ndarray,
    F_right: np.ndarray,
    P: np.ndarray,
    left_reps: np.ndarray,
    right_reps: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the expanded orbit-pairing matrix over GF(p).

    If ``left_reps`` and ``right_reps`` contain m_l and m_r chain-head
    representatives, this returns

        M[(i,a),(j,b)] = <N^{L-1} F_left^a u_i, F_right^b v_j>

    in the current left/right sector coordinates, together with the ordered left
    head orbit columns and right head orbit columns.  It is the explicit
    base-field representation of the full pairing data between the two
    F-orbits; for deg(q)=1 it reduces to the usual scalar head-tail matrix.
    """
    deg_q = int(deg_q)
    L = int(L)
    N_pow = _mat_pow_mod(N_left, L - 1, p)

    left_heads = _ordered_orbit_columns(F_left, left_reps, deg_q, p)
    right_heads = _ordered_orbit_columns(F_right, right_reps, deg_q, p)
    if left_heads.shape[1] == 0 or right_heads.shape[1] == 0:
        return (
            np.zeros((left_heads.shape[1], right_heads.shape[1]), dtype=np.int64),
            left_heads,
            right_heads,
        )

    left_tails = mod_p(N_pow @ left_heads, p)
    M = mod_p(left_tails.T @ (P @ right_heads), p)
    return M, left_heads, right_heads


def _verified_block_from_orbit_pairing(
    *,
    Fq: np.ndarray,
    Nq: np.ndarray,
    Fqs: np.ndarray,
    Nqs: np.ndarray,
    Vq_r: np.ndarray,
    Vqs_r: np.ndarray,
    P: np.ndarray,
    left_reps: np.ndarray,
    right_reps: np.ndarray,
    deg_q: int,
    L: int,
    p: int,
    Omega_amb: np.ndarray,
    rem: np.ndarray,
    rem_dim: int,
) -> tuple[np.ndarray, dict[str, Any]] | tuple[None, dict[str, Any]]:
    """
    Construct one paired block using the full expanded orbit-pairing matrix.

    The routine first forms the ordinary GF(p) matrix of pairings between all
    F-orbit components of the selected length-L heads.  If this matrix is full
    rank, Gaussian elimination gives dual right-hand orbit combinations.  Each
    dual combination is then tried as a right chain head and the lifted cyclic
    block is verified directly in the ambient symplectic space.
    """
    M, left_heads, right_heads = _orbit_pairing_matrix(
        F_left=Fq,
        N_left=Nq,
        F_right=Fqs,
        P=P,
        left_reps=left_reps,
        right_reps=right_reps,
        deg_q=deg_q,
        L=L,
        p=p,
    )

    orbit_dim_left = int(left_heads.shape[1])
    orbit_dim_right = int(right_heads.shape[1])
    rM = rank_mod(M, p)
    diag: dict[str, Any] = {
        "L": int(L),
        "deg_q": int(deg_q),
        "n_left_reps": int(left_reps.shape[1]),
        "n_right_reps": int(right_reps.shape[1]),
        "orbit_pairing_shape": (int(M.shape[0]), int(M.shape[1])),
        "orbit_pairing_rank": int(rM),
        "orbit_pairing_full_rank": bool(
            M.shape[0] == M.shape[1] and rM == M.shape[0]
        ),
        "orbit_pairing_matrix_constructed": True,
        "dual_solve_used": False,
    }

    if orbit_dim_left == 0 or orbit_dim_right == 0:
        diag["reason"] = "empty orbit basis"
        return None, diag
    if orbit_dim_left != orbit_dim_right:
        diag["reason"] = "left/right orbit dimensions differ"
        return None, diag
    if rM != orbit_dim_left:
        diag["reason"] = "orbit-pairing matrix is rank deficient"
        return None, diag

    # M is square and invertible.  For a chosen left orbit component e_k, solve
    # M x = e_k.  The right head z = right_heads @ x has canonical pairings with
    # the selected left orbit components.  We then verify the cyclic block built
    # from the corresponding left orbit component and z.
    for target_row in range(orbit_dim_left):
        rhs = np.zeros((orbit_dim_left, 1), dtype=np.int64)
        rhs[target_row, 0] = 1
        try:
            coeff = solve_linear(M, rhs, p)
        except RuntimeError:
            continue

        v_top = left_heads[:, target_row:target_row + 1]
        w_top = mod_p(right_heads @ coeff, p)
        diag["dual_solve_used"] = True
        diag["chosen_left_orbit_component"] = int(target_row % int(deg_q))
        diag["chosen_left_representative"] = int(target_row // int(deg_q))

        try:
            C_left = cyclic_submodule_basis(Fq, Nq, v_top, deg_q, int(L), p)
            C_right = cyclic_submodule_basis(Fqs, Nqs, w_top, deg_q, int(L), p)
        except RuntimeError as exc:
            diag["last_cyclic_error"] = f"{type(exc).__name__}: {exc}"
            continue

        C_left = independent_columns(mod_p(C_left, p), p)
        C_right = independent_columns(mod_p(C_right, p), p)

        W_left = mod_p(Vq_r @ C_left, p)
        W_right = mod_p(Vqs_r @ C_right, p)
        span = independent_columns(np.concatenate([W_left, W_right], axis=1), p)

        expected_dim = 2 * int(deg_q) * int(L)
        if span.shape[1] != expected_dim:
            diag["last_span_error"] = (
                f"wrong span dimension at L={L}: expected {expected_dim}, got {span.shape[1]}"
            )
            continue
        if not is_nondegenerate(Omega_amb, span, p):
            diag["last_span_error"] = f"constructed span degenerate at L={L}"
            continue
        if rank_mod(np.concatenate([rem, span], axis=1), p) != rem_dim:
            diag["last_span_error"] = f"constructed span not contained in remaining subspace at L={L}"
            continue

        diag["verified"] = True
        diag["block_dim"] = int(span.shape[1])
        return span, diag

    diag["reason"] = "all dual orbit representatives failed direct block verification"
    return None, diag


def atomic_blocks_in_paired_sector(
    F: np.ndarray,
    p: int,
    key: tuple[int, ...],
    key_star: tuple[int, ...],
    primaries: dict,
) -> tuple[list[AtomicBlock], AtomicInvariant]:
    """
    Atomic block construction for paired sector W = V_q ⊕ V_{q*} (q != q*).

    Certified ("OK") iff:
      - for every length L, the induced chain-level pairing between the chosen generators
        is non-singular (enforced by construction),
      - each constructed block span is nondegenerate,
      - and the blocks span the whole paired sector.
    """
    F = mod_p(F, p)

    q = primaries[key]["poly"]
    q_star = primaries[key_star]["poly"]

    Vq = independent_columns(mod_p(primaries[key]["V_basis"], p), p)
    Vqs = independent_columns(mod_p(primaries[key_star]["V_basis"], p), p)

    inv_data: Dict[str, Any] = {
        "status": "PENDING",
        "deg": int(primaries[key]["deg"]),
        "exponent": int(primaries[key]["exponent"]),
        "length_multiplicities": {},
        "pairing_rank": {},
        "checks_passed": [],
        "note": "",
    }

    if Vq.shape[1] == 0 or Vqs.shape[1] == 0:
        raise RuntimeError("Paired sector: one side of paired sector is empty.")

    # --- Certified-by-construction algorithm (extract-and-remove) ---
    # The previous "bulk" builder can produce invariant nondegenerate summands that are not
    # symplectically orthogonal to each other, which later makes the global basis non-symplectic.
    # Here we instead extract one atomic block at a time and remove its symplectic orthogonal
    # complement inside the paired sector. This guarantees pairwise orthogonality of blocks.

    n2 = F.shape[0]
    Omega_amb = omega_matrix(n2 // 2, p)

    max_exp = int(primaries[key]["exponent"])
    deg_q = int(primaries[key]["deg"])
    deg_qs = int(primaries[key_star]["deg"])
    if deg_q != deg_qs:
        raise RuntimeError("Paired sector: deg(q) != deg(q*) (unexpected).")

    # The full paired sector span W = Vq ⊕ Vq*
    W_sector = independent_columns(np.concatenate([Vq, Vqs], axis=1), p)
    if not is_nondegenerate(Omega_amb, W_sector, p):
        raise RuntimeError("Paired sector: sector span is degenerate (unexpected).")

    inv_data.setdefault("progress", [])  # list of dicts per extracted block

    blocks: List[AtomicBlock] = []

    try:
        rem = W_sector
        rem_dim = rank_mod(rem, p)

        # Safety: prevent infinite loops if something goes wrong.
        max_iters = rem_dim // (2 * deg_q) + 5
        it = 0

        while rem_dim > 0:
            it += 1
            if it > max_iters:
                raise SearchBudgetExceeded("Paired sector: extraction stuck (too many iterations).")

            # Recompute the current left/right parts inside the remaining invariant subspace.
            Vq_r = _intersection_basis(Vq, rem, p)
            Vqs_r = _intersection_basis(Vqs, rem, p)
            if Vq_r.shape[1] == 0 or Vqs_r.shape[1] == 0:
                raise RuntimeError("Paired sector: remaining subspace lost one side (unexpected).")

            # Restrict to each side in its own coordinates.
            Fq = restrict_operator_invariant(F, Vq_r, p)
            Fqs = restrict_operator_invariant(F, Vqs_r, p)
            Nq = q_of_F_restricted(Fq, q, p)
            Nqs = q_of_F_restricted(Fqs, q_star, p)

            filt_left = build_nilpotent_filtration(Nq, np.eye(Fq.shape[0], dtype=np.int64), max_exp, p)
            filt_right = build_nilpotent_filtration(Nqs, np.eye(Fqs.shape[0], dtype=np.int64), max_exp, p)
            tops_left = filt_left.tops
            tops_right = filt_right.tops
            if not tops_left and max_exp >= 1:
                tops_left = {1: np.eye(Fq.shape[0], dtype=np.int64)}
            if not tops_right and max_exp >= 1:
                tops_right = {1: np.eye(Fqs.shape[0], dtype=np.int64)}

            lengths = sorted(set(tops_left.keys()) | set(tops_right.keys()), reverse=True)
            if not lengths and max_exp >= 1:
                lengths = [1]

            # Pairing between current left/right bases.
            P = _pairing_matrix_between(Vq_r, Vqs_r, p)

            extracted = False
            last_err: str | None = None
            for L in lengths:
                A_raw = tops_left.get(L, np.zeros((Fq.shape[0], 0), dtype=np.int64))
                B_raw = tops_right.get(L, np.zeros((Fqs.shape[0], 0), dtype=np.int64))

                A = select_module_generators_from_top_quotient(
                    Fq, Nq, A_raw, deg_q, int(L), p,
                    denom=filt_left.denom.get(int(L), np.zeros((Fq.shape[0], 0), dtype=np.int64)),
                )
                pool = _candidate_generators_from_top_space(
                    Fqs, Nqs, B_raw, deg_q, int(L), p,
                    denom=filt_right.denom.get(int(L), np.zeros((Fqs.shape[0], 0), dtype=np.int64)),
                )

                # Bookkeeping only (not used by the algorithm)
                inv_data["length_multiplicities"][int(L)] = (int(A.shape[1]), int(pool.shape[1]))

                if A.shape[1] == 0 or pool.shape[1] == 0:
                    continue

                # Build the full expanded orbit-pairing matrix for the selected
                # length-L head representatives and solve its dual system over GF(p).
                # This is the explicit linear-algebraic version of the reciprocal
                # orbit pairing; the older scalar-component scan is no longer used.
                span, pair_diag = _verified_block_from_orbit_pairing(
                    Fq=Fq,
                    Nq=Nq,
                    Fqs=Fqs,
                    Nqs=Nqs,
                    Vq_r=Vq_r,
                    Vqs_r=Vqs_r,
                    P=P,
                    left_reps=A,
                    right_reps=pool,
                    deg_q=deg_q,
                    L=int(L),
                    p=p,
                    Omega_amb=Omega_amb,
                    rem=rem,
                    rem_dim=rem_dim,
                )
                inv_data["pairing_rank"][int(L)] = int(pair_diag.get("orbit_pairing_rank", 0))

                if span is None:
                    last_err = f"orbit-pairing dual solve failed at L={L}: {pair_diag.get('reason', pair_diag)}"
                    inv_data["progress"].append({"L": int(L), "attempt": "orbit_pairing", **pair_diag})
                    continue

                T_blk = darboux_basis_from_span(Omega_amb, span, p)
                blocks.append(
                    AtomicBlock(
                        T_blk=mod_p(T_blk, p),
                        half_dim=int(T_blk.shape[1] // 2),
                        sector_key=key,
                        inv=None,
                    )
                )

                # Remove its symplectic orthogonal complement within the remaining subspace.
                rem2 = symplectic_orthogonal_complement_in_span(Omega_amb, span, rem, p)
                rem2 = independent_columns(mod_p(rem2, p), p)
                rem2_dim = rank_mod(rem2, p)
                drop = rem_dim - rem2_dim
                if drop != span.shape[1]:
                    raise RuntimeError(
                        f"Paired sector: rank-drop mismatch when removing block (expected {span.shape[1]}, got {drop})."
                    )
                progress_entry = {
                    "L": int(L),
                    "block_dim": int(span.shape[1]),
                    "rem_dim_before": int(rem_dim),
                    "rem_dim_after": int(rem2_dim),
                    "attempt": "orbit_pairing",
                }
                progress_entry.update(pair_diag)
                inv_data["progress"].append(progress_entry)

                rem = rem2
                rem_dim = rem2_dim
                extracted = True
                break

            if not extracted:
                raise ExtractionObstruction(
                    "Paired sector: could not extract a valid block from remaining subspace"
                    + (f" (last_err={last_err})" if last_err else "")
                )

        # Final span check
        all_cols = np.concatenate([b.T_blk for b in blocks], axis=1) if blocks else np.zeros((n2, 0), dtype=np.int64)
        if rank_mod(np.concatenate([W_sector, all_cols], axis=1), p) != rank_mod(W_sector, p):
            raise RuntimeError("Paired sector: some constructed block columns lie outside W = Vq ⊕ Vq*.")
        if rank_mod(all_cols, p) != rank_mod(W_sector, p):
            raise RuntimeError("Paired sector: constructed blocks do not span the paired sector.")

        inv_data["checks_passed"].append("extract-and-remove guarantees orthogonality")
        inv_data["checks_passed"].append("each block nondegenerate")
        inv_data["checks_passed"].append("blocks span paired sector")
        inv_data["status"] = "OK"
        sector_cost = max((int(b.half_dim) for b in blocks), default=0)
        inv_data["cost_certificate"] = {
            "lower_bound": int(sector_cost),
            "attained": True,
            "complete": True,
            "sector_cost": int(sector_cost),
            "note": (
                "paired sector certified by expanded orbit-pairing matrix dual solve "
                "and extract/remove orthogonalization"
            ),
        }

    except Exception:
        # Single certified route: a paired sector that cannot be built
        # deterministically propagates here and is absorbed by global completion
        # upstream (the result is then reported uncertified). No best-effort
        # single-block fallback is constructed.
        raise

    inv = AtomicInvariant(sector_key=key, sector_type="paired", poly_key=key, data=inv_data)
    blocks = [AtomicBlock(b.T_blk, b.half_dim, b.sector_key, inv) for b in blocks]

    # Sanity: blocks must be a direct sum (no overlap)
    if blocks:
        all_cols = np.concatenate([b.T_blk for b in blocks], axis=1)
        if rank_mod(all_cols, p) != sum(b.T_blk.shape[1] for b in blocks):
            raise RuntimeError("Paired sector: produced blocks overlap (not direct sum).")

    return blocks, inv
