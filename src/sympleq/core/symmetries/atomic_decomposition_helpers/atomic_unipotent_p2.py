# sympleq/core/symmetries/atomic_decomposition_helpers/atomic_unipotent_p2.py
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

from ..modular_helpers import mod_p, independent_columns, rank_mod, omega_matrix
from .atomic_types import (AtomicBlock, AtomicInvariant,
                           ExtractionObstruction, SearchBudgetExceeded)
from .atomic_linear import (
    restrict_operator,
    is_nondegenerate,
    darboux_basis_from_span,
    mat_pow_mod,
    symplectic_orthogonal_complement_in_span,
)
from .module_invariants import (
    cyclic_submodule_basis,
    q_of_F_restricted,
)
from .atomic_krylov import _select_module_generators_from_top_space
from .atomic_filtration import build_nilpotent_filtration


# ---------------------------------------------------------------------------
# GF(2) scalar helpers
# ---------------------------------------------------------------------------

def _scalar_mod2(x) -> int:
    """Safely extract a GF(2) scalar from a scalar-like NumPy expression."""
    arr = np.asarray(x, dtype=np.int64) % 2
    if arr.size != 1:
        raise ValueError(f"Expected scalar-like array, got shape={arr.shape}")
    return int(arr.reshape(-1)[0])



def jordan_chain_tops_nilpotent_in_span(
    N: np.ndarray, space_basis: np.ndarray, max_exp: int, p: int
) -> Dict[int, np.ndarray]:
    """
    Restricted length-top quotient representatives, now backed by the shared
    NilpotentFiltration cache.  Kept under the old name for compatibility with
    existing callers/tests.
    """
    return build_nilpotent_filtration(N, space_basis, max_exp, p).tops


# ---------------------------------------------------------------------------
# p=2 top bilinear/quadratic data
# ---------------------------------------------------------------------------

def _top_pairing_matrix(A_top: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int) -> np.ndarray:
    """Matrix of b_L([v],[w]) = <v, N^{L-1} w> on top representatives."""
    p = 2
    A_top = independent_columns(mod_p(A_top, p), p)
    Nr = np.eye(N.shape[0], dtype=np.int64) if L <= 1 else mat_pow_mod(N, L - 1, p)
    return mod_p(A_top.T @ Omega @ (Nr @ A_top), p)


def _mid_exponent_p2(L: int) -> int:
    """
    Exponent used by the Chapter-5-style quadratic refinement.

    For L=2k use k-1; for L=2l+1 use l.
    """
    L = int(L)
    if L % 2 == 0:
        return max(L // 2 - 1, 0)
    return (L - 1) // 2


def _q_values_on_top(A_top: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int) -> np.ndarray:
    """Return q_L(a_i)=<a_i, N^e a_i> for the selected top representatives."""
    p = 2
    A_top = independent_columns(mod_p(A_top, p), p)
    if A_top.shape[1] == 0:
        return np.zeros(0, dtype=np.int64)
    e = _mid_exponent_p2(int(L))
    Ne = np.eye(N.shape[0], dtype=np.int64) if e == 0 else mat_pow_mod(N, e, p)
    Q = mod_p(A_top.T @ Omega @ (Ne @ A_top), p)
    return np.diag(Q).astype(np.int64) % 2


def _q_value_from_coeff(coeff: np.ndarray, A_top: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int) -> int:
    """Evaluate q_L on a top vector represented by coefficient vector ``coeff``."""
    p = 2
    coeff = mod_p(np.asarray(coeff, dtype=np.int64).reshape(-1, 1), p)
    v = mod_p(A_top @ coeff, p)
    e = _mid_exponent_p2(int(L))
    Ne = np.eye(N.shape[0], dtype=np.int64) if e == 0 else mat_pow_mod(N, e, p)
    return int(mod_p(v.T @ Omega @ (Ne @ v), p).reshape(())) & 1


def _beta_from_top_generators_p2(
    gens: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int
) -> int:
    """
    Chapter-5-style beta label of an extracted p=2 indecomposable.

    For a single top generator this is q_L(v).  For a two-generator W-plane
    normalized by b_L(v,w)=1, this returns the Arf contribution q_L(v) q_L(w).
    """
    p = 2
    gens = mod_p(gens, p)
    if gens.shape[1] == 0:
        return 0

    e = _mid_exponent_p2(int(L))
    Ne = np.eye(N.shape[0], dtype=np.int64) if e == 0 else mat_pow_mod(N, e, p)

    def q(vcol: np.ndarray) -> int:
        vcol = mod_p(vcol.reshape(-1, 1), p)
        return int(mod_p(vcol.T @ Omega @ (Ne @ vcol), p).reshape(())) & 1

    if gens.shape[1] == 1:
        return q(gens[:, 0])
    if gens.shape[1] == 2:
        return q(gens[:, 0]) & q(gens[:, 1])
    return 0


def _p2_hyperbolic_pairs_from_alternating_form(
    B: np.ndarray, p: int = 2
) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Deterministic Gram-Schmidt for an alternating form over GF(2).

    Returns coefficient vectors in the original top basis: hyperbolic pairs
    ``(e_i,f_i)`` and radical lines.  This is a quotient-level construction;
    actual cyclic modules are only accepted after direct nondegeneracy checks.
    """
    if p != 2:
        raise ValueError("_p2_hyperbolic_pairs_from_alternating_form is p=2 only")
    B = mod_p(np.asarray(B, dtype=np.int64), 2)
    if B.ndim != 2 or B.shape[0] != B.shape[1]:
        raise ValueError(f"B must be square, got {B.shape}.")
    if np.any(np.diag(B) % 2):
        raise RuntimeError("p=2 top pairing is not alternating: nonzero diagonal.")

    m = B.shape[0]
    R = [np.eye(m, dtype=np.int64)[:, i] for i in range(m)]

    e_list: List[np.ndarray] = []
    f_list: List[np.ndarray] = []
    rad_list: List[np.ndarray] = []

    def pair(u: np.ndarray, v: np.ndarray) -> int:
        return _scalar_mod2(u.reshape(1, -1) @ B @ v.reshape(-1, 1))

    while R:
        a = R.pop(0)
        found = None
        for idx, cand in enumerate(R):
            if pair(a, cand) == 1:
                found = idx
                break

        if found is None:
            rad_list.append(a.copy())
            continue

        b = R.pop(found)

        newR = []
        for x in R:
            ax = pair(x, b)
            bx = pair(x, a)
            if ax:
                x = (x ^ a)
            if bx:
                x = (x ^ b)
            newR.append(x)
        R = newR

        e_list.append(a.copy())
        f_list.append(b.copy())

    return e_list, f_list, rad_list


def _radical_candidates_adapted_to_q(
    rad_list: List[np.ndarray], A_top: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int
) -> List[np.ndarray]:
    """
    Return deterministic radical-line candidates, lightly adapted to q_L.

    On the radical of b_L the polar form of q_L vanishes, so q_L is additive.
    If a q=1 line exists, the first combination found below exposes it; q=0
    lines are retained as well because V_0-type blocks can occur.
    """
    if not rad_list:
        return []

    out: List[np.ndarray] = []
    seen: set[Tuple[int, ...]] = set()

    def push(c: np.ndarray) -> None:
        cc = mod_p(np.asarray(c, dtype=np.int64).reshape(-1), 2)
        if not np.any(cc):
            return
        key = tuple(int(x) for x in cc)
        if key not in seen:
            seen.add(key)
            out.append(cc)

    for c in rad_list:
        push(c)

    # If the basis lines all have q=0 but a q=1 combination exists, a short
    # deterministic combination exposes it in the common cases.  For modest
    # radical dimension enumerate all combinations; otherwise use pair sums.
    r = len(rad_list)
    if r <= 12:
        R = np.stack(rad_list, axis=1)
        for mask in range(1, 1 << r):
            coeff = np.zeros(r, dtype=np.int64)
            for i in range(r):
                if (mask >> i) & 1:
                    coeff[i] = 1
            push(R @ coeff)
    else:
        for i in range(r):
            for j in range(i + 1, r):
                push(rad_list[i] ^ rad_list[j])

    # Stable order: prefer q=1 witnesses first, then q=0, preserving discovery
    # order inside each bucket.
    ones = [c for c in out if _q_value_from_coeff(c, A_top, Omega, N, L) == 1]
    zeros = [c for c in out if _q_value_from_coeff(c, A_top, Omega, N, L) == 0]
    return ones + zeros



def _top_coeff_candidates(m: int, *, exhaustive_limit: int = 10) -> List[np.ndarray]:
    """
    Deterministic nonzero coefficient vectors in GF(2)^m.

    For modest top dimension this is exhaustive, which makes the implemented
    W(k) and V_beta(2k) extraction independent of the arbitrary top basis.  For
    larger dimensions we use a deterministic spanning stress set; the direct
    quotient decompositions normally make the exhaustive branch sufficient for
    atomic top pieces encountered after repeated extraction.
    """
    m = int(m)
    if m <= 0:
        return []
    out: List[np.ndarray] = []
    seen: set[Tuple[int, ...]] = set()

    def push(c: np.ndarray) -> None:
        cc = mod_p(np.asarray(c, dtype=np.int64).reshape(-1), 2)
        if cc.shape[0] != m or not np.any(cc):
            return
        key = tuple(int(x) for x in cc)
        if key not in seen:
            seen.add(key)
            out.append(cc)

    # Basis vectors first for reproducibility and readable diagnostics.
    for i in range(m):
        e = np.zeros(m, dtype=np.int64)
        e[i] = 1
        push(e)

    if m <= exhaustive_limit:
        for mask in range(1, 1 << m):
            c = np.zeros(m, dtype=np.int64)
            for i in range(m):
                if (mask >> i) & 1:
                    c[i] = 1
            push(c)
    else:
        # Deterministic but bounded supplement for large quotient dimensions.
        for i in range(m):
            for j in range(i + 1, m):
                c = np.zeros(m, dtype=np.int64)
                c[i] = 1
                c[j] = 1
                push(c)
        push(np.ones(m, dtype=np.int64))

    return out


def _hyperbolic_pair_candidates(Btop: np.ndarray) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Deterministic quotient-normal-form W-pair extraction.

    This certified helper does **not** enumerate arbitrary coefficient vectors.
    For alternating nonsingular top forms it uses symplectic Gram-Schmidt.  For
    nonsingular non-alternating symmetric GF(2) forms, it first diagonalizes the
    form into orthogonal anisotropic lines ``u_i`` and then pairs adjacent lines
    as ``(u_i, u_i + u_j)``.  The resulting pair has

        b(u_i, u_i + u_j) = 1,

    with a nonsingular two-dimensional top Gram, which is the W/W_beta lifting
    certificate used by the p=2 unipotent extractor.
    """
    Btop = mod_p(np.asarray(Btop, dtype=np.int64), 2)
    if Btop.ndim != 2 or Btop.shape[0] != Btop.shape[1]:
        raise ValueError(f"Btop must be square, got {Btop.shape}.")
    m = int(Btop.shape[0])
    if m == 0:
        return []
    if not np.array_equal(Btop, Btop.T):
        raise RuntimeError("p=2 top form is not symmetric; cannot extract certified W-pairs.")

    if np.all(np.diag(Btop) % 2 == 0):
        e_list, f_list, rad_list = _p2_hyperbolic_pairs_from_alternating_form(Btop, p=2)
        if rad_list:
            raise RuntimeError("alternating top form has a radical in certified W-pair extraction.")
        return [(mod_p(e, 2), mod_p(f, 2)) for e, f in zip(e_list, f_list)]

    anis = _anisotropic_orthogonal_line_candidates(Btop)
    if len(anis) != m:
        raise RuntimeError(
            f"non-alternating top form did not diagonalize to a full anisotropic basis: {len(anis)} != {m}"
        )
    if len(anis) % 2 != 0:
        raise RuntimeError("non-alternating W-pair extraction needs an even number of anisotropic lines.")

    out: List[Tuple[np.ndarray, np.ndarray]] = []
    for i in range(0, len(anis), 2):
        u = mod_p(np.asarray(anis[i], dtype=np.int64).reshape(-1), 2)
        v = mod_p(np.asarray(anis[i + 1], dtype=np.int64).reshape(-1), 2)
        e = u
        f = mod_p(u + v, 2)
        if _scalar_mod2(e.reshape(1, -1) @ Btop @ f.reshape(-1, 1)) != 1:
            raise RuntimeError("internal error: constructed non-alternating W-pair has zero top pairing.")
        out.append((e, f))
    return out


def _anisotropic_orthogonal_line_candidates(Btop: np.ndarray) -> List[np.ndarray]:
    """
    Deterministically diagonalize the non-alternating part of a symmetric
    GF(2) top pairing.

    Returns coefficient vectors c in the original top basis such that

        c_i^T Btop c_i = 1,
        c_i^T Btop c_j = 0  for i != j

    for the extracted anisotropic lines.  In the length-2 unipotent case,
    these are precisely the quotient lines that should lift to one-qudit
    V_beta(2)-type blocks when the direct cyclic-module check succeeds.

    If the remaining form becomes alternating, extraction stops and the W-pair
    path handles that alternating remainder.
    """
    B = mod_p(np.asarray(Btop, dtype=np.int64), 2).copy()
    if B.ndim != 2 or B.shape[0] != B.shape[1]:
        raise ValueError(f"Btop must be square, got {B.shape}.")
    m = int(B.shape[0])
    if m == 0:
        return []

    # P columns are current basis vectors in the original top coordinates.
    P = np.eye(m, dtype=np.int64)
    out: List[np.ndarray] = []
    k = 0

    def swap(i: int, j: int) -> None:
        nonlocal B, P
        if i == j:
            return
        B[[i, j], :] = B[[j, i], :]
        B[:, [i, j]] = B[:, [j, i]]
        P[:, [i, j]] = P[:, [j, i]]

    def shear_col(i: int, j: int) -> None:
        """Basis update e_i <- e_i + e_j, with congruence update."""
        nonlocal B, P
        B[:, i] = mod_p(B[:, i] + B[:, j], 2)
        B[i, :] = mod_p(B[i, :] + B[j, :], 2)
        P[:, i] = mod_p(P[:, i] + P[:, j], 2)

    while k < m:
        piv = None
        for i in range(k, m):
            if int(B[i, i] % 2) == 1:
                piv = i
                break
        if piv is None:
            break
        swap(piv, k)

        # Orthogonalize the remaining basis vectors against e_k.  Since
        # <e_k,e_k>=1, replacing e_t by e_t + <e_k>e_t e_k kills <e_k,e_t>.
        for t in range(k + 1, m):
            if int(B[k, t] % 2) == 1:
                shear_col(t, k)

        c = mod_p(P[:, k].reshape(-1), 2)
        if np.any(c):
            out.append(c.copy())
        k += 1

    # If a greedy anisotropic pivot leaves an alternating hyperbolic remainder,
    # do not return the bad decomposition [anisotropic] + H.  Over GF(2),
    # an anisotropic line u orthogonal to a hyperbolic pair (e,f) can be
    # replaced by three mutually orthogonal anisotropic lines
    #
    #     u+e,  u+f,  u+e+f.
    #
    # This is the key normalization needed for length-2 p=2 unipotent sectors:
    # a full-rank non-alternating top form should be diagonalized into V-type
    # one-qudit lines whenever the lifted cyclic-module checks succeed, rather
    # than leaving an alternating remainder to be extracted as W blocks.
    if k < m and out:
        Brem = mod_p(B[k:, k:], 2)
        if Brem.size and np.all(np.diag(Brem) % 2 == 0):
            try:
                e_rem, f_rem, _rad_rem = _p2_hyperbolic_pairs_from_alternating_form(Brem, p=2)
            except RuntimeError:
                e_rem, f_rem = [], []

            if e_rem:
                Prem = mod_p(P[:, k:], 2)

                def lift_rem(c_rem: np.ndarray) -> np.ndarray:
                    return mod_p(Prem @ np.asarray(c_rem, dtype=np.int64).reshape(-1, 1), 2).reshape(-1)

                carrier = out.pop()
                for er, fr in zip(e_rem, f_rem):
                    e = lift_rem(er)
                    f = lift_rem(fr)
                    a = mod_p(carrier + e, 2)
                    b = mod_p(carrier + f, 2)
                    c = mod_p(carrier + e + f, 2)
                    out.append(a.reshape(-1).copy())
                    out.append(b.reshape(-1).copy())
                    carrier = c.reshape(-1)
                out.append(carrier.reshape(-1).copy())

    # Remove accidental duplicates while preserving order.
    unique: List[np.ndarray] = []
    seen: set[Tuple[int, ...]] = set()
    for c in out:
        cc = mod_p(np.asarray(c, dtype=np.int64).reshape(-1), 2)
        key = tuple(int(x) for x in cc)
        if np.any(cc) and key not in seen:
            seen.add(key)
            unique.append(cc)

    return unique

def _p2_length_form_invariants(A_top: np.ndarray, Omega: np.ndarray, N: np.ndarray, L: int) -> Dict[str, Any]:
    p = 2
    A_top = independent_columns(mod_p(A_top, p), p)
    m = int(A_top.shape[1])
    if m == 0:
        return {
            "top_dim": 0,
            "B_rank": 0,
            "rad_dim": 0,
            "B_sym_ok": True,
            "B_alt_ok": True,
            "q_values": [],
            "q1_count": 0,
            "n_hyp": 0,
            "arf": None,
        }

    B = _top_pairing_matrix(A_top, Omega, N, int(L))
    qvals = _q_values_on_top(A_top, Omega, N, int(L))
    B_rank = int(rank_mod(B, p))
    rad_dim = int(m - B_rank)
    B_sym_ok = bool(np.array_equal(B, B.T))
    B_alt_ok = bool(np.all(np.diag(B) % 2 == 0))

    arf = None
    n_hyp = 0
    if B_alt_ok:
        e_list, f_list, _ = _p2_hyperbolic_pairs_from_alternating_form(B, p=2)
        n_hyp = len(e_list)
        arf_val = 0
        for e, f in zip(e_list, f_list):
            arf_val ^= (_q_value_from_coeff(e, A_top, Omega, N, L) & _q_value_from_coeff(f, A_top, Omega, N, L))
        if rad_dim == 0:
            arf = int(arf_val)

    return {
        "top_dim": int(m),
        "B_rank": int(B_rank),
        "rad_dim": int(rad_dim),
        "B_sym_ok": bool(B_sym_ok),
        "B_alt_ok": bool(B_alt_ok),
        "q_values": [int(x) for x in qvals.reshape(-1)],
        "q1_count": int(np.sum(qvals % 2)),
        "n_hyp": int(n_hyp),
        "arf": arf,
    }


# ---------------------------------------------------------------------------
# Main builder (p=2 unipotent self sector)
# ---------------------------------------------------------------------------

def _p2_full_sector_invariants(N, Omega, F_sec, deg_q, max_exp0, m, p):
    """Step 2 (sec. 4.9): invariants first.

    Compute, before any extraction, the sector kernel profile and the per-length
    top-form invariants b_L (via :func:`_p2_length_form_invariants`) plus a
    generator-count diagnostic. These determine the V/W type at each length
    (Lemma 4.5) and feed the cost certificate.
    """
    # Full-sector invariant summary.
    tops_full = jordan_chain_tops_nilpotent_in_span(N, np.eye(2 * m, dtype=np.int64), max_exp0, p)
    kernel_profile = [
        int((2 * m) - rank_mod(mat_pow_mod(N, k, p), p))
        for k in range(1, int(max_exp0) + 1)
    ]

    length_summary: Dict[int, Dict[str, Any]] = {}
    length_invariants: Dict[int, Dict[str, Any]] = {}
    for L, Araw in sorted(tops_full.items()):
        Araw = independent_columns(mod_p(Araw, p), p)
        invL = _p2_length_form_invariants(Araw, Omega, N, int(L))
        length_invariants[int(L)] = invL
        try:
            A_gen = _select_module_generators_from_top_space(F_sec, N, Araw, deg_q, int(L), p)
            gen_dim = int(A_gen.shape[1])
        except RuntimeError:
            gen_dim = -1
        length_summary[int(L)] = {
            "mult": int(invL["top_dim"]),
            "gen_dim": int(gen_dim),
            "gen_dim_ok": bool(gen_dim == int(invL["top_dim"])),
            **invL,
        }
    return kernel_profile, length_summary, length_invariants


def _p2_extract_blocks(F_sec, N, Omega, T_sec, m, max_exp0, deg_q, key, p):
    """Steps 3-4 (sec. 4.9): type-by-invariant extraction engine.

    Peels the longest remaining length one block at a time. ``accept_block``
    records a lifted, verified summand (Cor. 4.12 one-scalar certificate, with a
    direct nondegeneracy/invariance assertion); ``try_invariant_candidates``
    performs the mandatory b_L dispatch (V-only / W-only, Prop. 4.7, Cor. 4.12).
    Returns the blocks, their sector-coordinate bases, per-block metadata, and
    the extraction log.
    """
    blocks: List[AtomicBlock] = []
    blocks_meta: List[Dict[str, Any]] = []
    built_cols_sec: List[np.ndarray] = []
    extraction_log: List[Dict[str, Any]] = []

    space_basis = np.eye(2 * m, dtype=np.int64)
    iteration = 0

    def accept_block(span: np.ndarray, top_gens: np.ndarray, L: int, typ_base: str) -> bool:
        nonlocal space_basis
        span = independent_columns(mod_p(span, p), p)
        if span.shape[1] == 0 or span.shape[1] % 2 != 0:
            return False
        if not is_nondegenerate(Omega, span, p):
            return False

        # Invariance guardrail: the span must be genuinely F_sec- and N-invariant.
        # NOTE: restrict_operator is a symplectic left-inverse *projection* and does
        # not raise on a non-invariant span, so it cannot be used as a guard. Test
        # invariance directly via rank([span | op @ span]) == rank(span).
        def _is_op_invariant(op: np.ndarray, sub: np.ndarray) -> bool:
            r = rank_mod(sub, p)
            aug = np.concatenate([sub, mod_p(op @ sub, p)], axis=1)
            return rank_mod(aug, p) == r

        if not _is_op_invariant(F_sec, span) or not _is_op_invariant(N, span):
            return False

        beta = _beta_from_top_generators_p2(top_gens, Omega, N, int(L))
        if typ_base == "V":
            # Implemented here: single cyclic even-length V_beta(2k) blocks.
            # A q=0 single-chain block is accepted only if the direct
            # nondegeneracy check above succeeds; it is still an implemented
            # V-type summand for the purposes of this sector decomposition.
            if int(L) % 2 != 0:
                return False
            typ = "V_beta" if beta else "V"
        else:
            # Implemented here: W(k) hyperbolic paired-chain blocks and
            # W_beta(2l+1) odd quadratic-pair blocks.  The same lifted
            # two-chain span is verified directly; the beta label records the
            # Chapter-5 odd quadratic refinement when L is odd.
            if int(L) % 2 == 1 and int(L) >= 3 and beta:
                typ = "W_beta"
            else:
                typ = "W"

        T_blk_sec = darboux_basis_from_span(Omega, span, p)
        rem = symplectic_orthogonal_complement_in_span(Omega, T_blk_sec, space_basis, p)
        T_blk_amb = mod_p(T_sec @ T_blk_sec, p)

        blocks.append(AtomicBlock(T_blk=T_blk_amb, half_dim=int(T_blk_amb.shape[1] // 2), sector_key=key, inv=None))
        built_cols_sec.append(T_blk_sec)
        blocks_meta.append(
            {
                "type": typ,
                "L": int(L),
                "half_dim": int(T_blk_amb.shape[1] // 2),
                "beta": int(beta),
                "top_dim": int(top_gens.shape[1]),
            }
        )
        space_basis = rem
        return True

    def try_invariant_candidates(A_top: np.ndarray, L: int) -> bool:
        """
        Type-by-invariant extraction (Phase 3 / sec. 4.9, Step 3).

        Dispatch is *mandatory* on the top form b_L:
          * even L, b_L non-alternating (diag != 0): extract V-blocks only, by
            orthonormalising the top form (Prop. 4.7, the three-line replacement
            in ``_anisotropic_orthogonal_line_candidates``) and lifting each
            anisotropic line via Cor. 4.12(a).  There is NO fall-through to
            W-pairs: a W here would have half-dimension L instead of L/2 and is
            therefore cost-pessimal.
          * even L with b_L alternating, or odd L: extract W-pairs only, via
            symplectic Gram-Schmidt on the alternating form
            (``_p2_hyperbolic_pairs_from_alternating_form``), lifting each
            hyperbolic top pair via Cor. 4.12(b).

        Acceptance carries the one-scalar Cor. 4.12 certificate (b_L(v,v)=1 for a
        V-line; b_L(v,w)=1 with the top Gram nonsingular for a W-pair); the direct
        nondegeneracy/invariance test in accept_block is retained as a defensive
        runtime assertion.  No arbitrary coefficient enumeration is used on the
        certified path.
        """
        L = int(L)
        Btop = _top_pairing_matrix(A_top, Omega, N, L)
        invL = _p2_length_form_invariants(A_top, Omega, N, L)

        # Prop. 4.3 (runtime assertion): on a nondegenerate sector the top form is
        # nondegenerate, i.e. rad b_L = 0.  A nonzero radical means an invariant
        # was computed wrong; surface it rather than extracting a wrong block.
        if int(invL["rad_dim"]) != 0:
            raise ExtractionObstruction(
                f"p=2 unipotent: top form b_L has nonzero radical at L={L} "
                f"(rad_dim={invL['rad_dim']}); Prop. 4.3 violated."
            )

        b_alt = bool(invL["B_alt_ok"])

        def _bv(c1: np.ndarray, c2: np.ndarray) -> int:
            return _scalar_mod2(c1.reshape(1, -1) @ Btop @ c2.reshape(-1, 1))

        if L % 2 == 0 and not b_alt:
            # --- V-type: non-alternating symmetric top form -> single chains ---
            for c in _anisotropic_orthogonal_line_candidates(Btop):
                c = mod_p(np.asarray(c, dtype=np.int64).reshape(-1), 2)
                if _bv(c, c) != 1:            # Cor. 4.12(a): require b_L(v,v) = 1
                    continue
                v_top = mod_p(A_top @ c.reshape(-1, 1), p)
                try:
                    Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, L, p)
                except RuntimeError:
                    continue
                Cv = independent_columns(Cv, p)
                if accept_block(Cv, v_top, L, "V"):
                    return True
            # A non-alternating top form always has an orthonormal basis
            # (Prop. 4.7(B)); if none lifted, that is an obstruction to surface,
            # never a licence to fall back to a cost-pessimal W-pair.
            return False

        # --- W-type: alternating top form (even L, diag 0) or odd L ---
        for ce, cf in _hyperbolic_pair_candidates(Btop):
            ce = mod_p(np.asarray(ce, dtype=np.int64).reshape(-1), 2)
            cf = mod_p(np.asarray(cf, dtype=np.int64).reshape(-1), 2)
            # Cor. 4.12(b): need b_L(v,w)=1 and a nonsingular 2x2 top Gram, i.e.
            # NOT both lines anisotropic (that would be two V-lines, not a W-pair).
            if _bv(ce, cf) != 1 or (_bv(ce, ce) == 1 and _bv(cf, cf) == 1):
                continue
            v_top = mod_p(A_top @ ce.reshape(-1, 1), p)
            w_top = mod_p(A_top @ cf.reshape(-1, 1), p)
            try:
                Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, L, p)
                Cw = cyclic_submodule_basis(F_sec, N, w_top, deg_q, L, p)
            except RuntimeError:
                continue
            span = independent_columns(np.concatenate([Cv, Cw], axis=1), p)
            if accept_block(span, np.concatenate([v_top, w_top], axis=1), L, "W"):
                return True

        return False

    while space_basis.shape[1] > 0:
        iteration += 1
        if iteration > 2 * m + 5:
            raise SearchBudgetExceeded("Unipotent p=2 self sector: extraction exceeded iteration guard.")

        tops = jordan_chain_tops_nilpotent_in_span(N, space_basis, max_exp0, p)
        if not tops:
            raise RuntimeError(
                "Unipotent p=2 self sector: remaining invariant subspace has no nilpotent top spaces."
            )

        L = max(tops.keys())
        A_top = independent_columns(mod_p(tops[L], p), p)
        invL = _p2_length_form_invariants(A_top, Omega, N, int(L))
        log_entry: Dict[str, Any] = {
            "iteration": int(iteration),
            "remaining_dim": int(space_basis.shape[1]),
            "L": int(L),
            **invL,
        }

        progressed = try_invariant_candidates(A_top, int(L))
        log_entry["progressed"] = bool(progressed)
        extraction_log.append(log_entry)

        if not progressed:
            raise ExtractionObstruction(
                "Unipotent p=2 self sector: certified quotient-driven extraction made no progress "
                f"at length L={int(L)}. Diagnostics={log_entry}"
            )
    return blocks, built_cols_sec, blocks_meta, extraction_log


def _p2_verify_kernel_profile(built_cols_sec, N, kernel_profile, max_exp0, m, p, key):
    """Step 4 (dagger) (sec. 4.9): consistency assertions.

    Assert the extracted blocks span the sector with non-overlapping bases, and
    that the per-power kernel profile is additive over the blocks (block sum ==
    sector). Returns the block kernel profile. Raises on any mismatch.
    """
    dim_sector = 2 * m
    all_cols_sec = np.concatenate(built_cols_sec, axis=1) if built_cols_sec else np.zeros((dim_sector, 0), dtype=np.int64)
    dim_blocks = rank_mod(all_cols_sec, p) if all_cols_sec.size else 0
    if dim_blocks != dim_sector:
        raise RuntimeError(
            f"Unipotent p=2 self sector: blocks do not span sector "
            f"(dim_blocks={dim_blocks}, dim_sector={dim_sector}). key={key}, exp0={max_exp0}"
        )
    if rank_mod(all_cols_sec, p) != all_cols_sec.shape[1]:
        raise RuntimeError("Unipotent p=2 self sector: extracted block bases overlap.")

    N_blks = [restrict_operator(N, T_blk, p) for T_blk in built_cols_sec]
    kernel_profile_blocks: List[int] = []
    for k in range(1, int(max_exp0) + 1):
        tot = 0
        for Nb in N_blks:
            Nk = mat_pow_mod(Nb, k, p)
            tot += int(Nb.shape[0] - rank_mod(Nk, p))
        kernel_profile_blocks.append(int(tot))

    if kernel_profile_blocks != kernel_profile:
        raise RuntimeError(
            "Unipotent p=2 self sector: kernel profile mismatch (block sum != sector). "
            f"key={key}, profile={kernel_profile}, block_profile={kernel_profile_blocks}"
        )
    return kernel_profile_blocks


def _p2_build_certificate(blocks_meta, kernel_profile, kernel_profile_blocks,
                          length_summary, length_invariants, extraction_log,
                          deg_q, max_exp0):
    """Steps 5-6 (sec. 4.9): cost certificate and invariant assembly.

    Assemble the sector invariant payload (beta/Arf counts per length, the
    implemented-family classification, and the sector cost certificate). The
    final ``lower_bound`` is overwritten downstream with the invariant-derived
    bound (Thm. 4.13) in _build_sector (Phase 2).
    """
    beta_counts_by_L: Dict[int, Dict[str, int]] = {}
    for bm in blocks_meta:
        L = int(bm["L"])
        beta = int(bm.get("beta", 0)) & 1
        typ = str(bm.get("type", ""))
        d = beta_counts_by_L.setdefault(L, {"V0": 0, "V1": 0, "W0": 0, "W1": 0})
        if typ.startswith("V"):
            d["V1" if beta else "V0"] += 1
        else:
            d["W1" if (typ.endswith("beta") or beta) else "W0"] += 1

    sector_cost = max((int(b["half_dim"]) for b in blocks_meta), default=0)
    used_best_effort = False  # the best-effort sweep was removed (deterministic path only)
    implemented_types = {"V", "V_beta", "W", "W_beta"}
    unimplemented = sorted({str(b.get("type", "")) for b in blocks_meta if str(b.get("type", "")) not in implemented_types})
    implemented_complete = (not used_best_effort) and (len(unimplemented) == 0)

    inv_data: Dict[str, Any] = {
        "status": "OK",
        "deg": int(deg_q),
        "exponent": int(max_exp0),
        "length_summary": length_summary,
        "blocks": blocks_meta,
        "extraction_log": extraction_log,
        "p2_unipotent": {
            "kernel_profile": kernel_profile,
            "kernel_profile_blocks": kernel_profile_blocks,
            "length_invariants": length_invariants,
            "beta_counts_by_L": beta_counts_by_L,
            "classification_complete": bool(implemented_complete),
            "classification_status": (
                "implemented Chapter-5 extraction for W(k), V_beta(2k), and W_beta(2l+1) succeeded"
                if implemented_complete else
                "p=2 unipotent extraction used a fallback component"
            ),
            "implemented_block_families": ["W(k)", "V_beta(2k)", "W_beta(2l+1)"],
            "unimplemented_block_families": [],
            "unimplemented_blocks_seen": unimplemented,
            "length_dropping_used": False,
            "used_best_effort_sweep": used_best_effort,
            "candidate_enumeration_used": False,
            "certified_extraction_policy": "quotient_top_normal_form_no_enumeration",
        },
        "cost_certificate": {
            "lower_bound": int(sector_cost) if implemented_complete else None,
            "attained": True,
            "complete": bool(implemented_complete),
            "sector_cost": int(sector_cost),
            "forced_block_types": blocks_meta,
            "note": (
                "p=2 unipotent sector certified using implemented W(k), V_beta(2k), and W_beta(2l+1) quotient-top extraction."
                if implemented_complete else
                "p=2 unipotent sector decomposed, but minimality is not certified because a fallback path was involved."
            ),
        },
    }
    return inv_data


def atomic_blocks_in_unipotent_self_sector_p2(
    F: np.ndarray,
    key: Tuple[int, ...],
    primaries: dict,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    p=2, self sector with q(x)=x±1.

    Deterministic type-by-invariant policy (no search/fallback):
      * at each top length L, dispatch on the top form b_L -- V-blocks when it is
        non-alternating (Cor. 4.12(a)), W-pairs when alternating or L is odd
        (Cor. 4.12(b));
      * lift quotient pieces to cyclic modules and accept with the one-scalar
        certificate (the direct nondegeneracy/invariance test is a defensive
        assertion);
      * never drop chain length.
    A sector that cannot be built this way raises; the single decomposition route
    absorbs it as a degraded fallback (marking the result uncertified).
    """
    p = 2
    F = mod_p(F, p)

    q = primaries[key]["poly"]
    deg_q = int(primaries[key]["deg"])
    max_exp0 = int(primaries[key]["exponent"])
    if deg_q != 1:
        raise RuntimeError(f"p=2 unipotent self-sector expected deg(q)=1, got {deg_q}.")

    V = independent_columns(mod_p(primaries[key]["V_basis"], p), p)
    if V.shape[1] == 0:
        inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data={"status": "OK", "empty": True})
        return [], inv

    n2 = F.shape[0]
    Omega_amb = omega_matrix(n2 // 2, p)
    if not is_nondegenerate(Omega_amb, V, p):
        raise RuntimeError("Unipotent p=2 self sector basis is degenerate; cannot proceed.")

    T_sec = darboux_basis_from_span(Omega_amb, V, p)
    m2 = T_sec.shape[1]
    if m2 % 2 != 0:
        raise RuntimeError("Unipotent p=2 self sector dimension must be even.")
    m = m2 // 2

    F_sec = restrict_operator(F, T_sec, p)
    Omega = omega_matrix(m, p)
    N = q_of_F_restricted(F_sec, q, p)

    # Step 2: invariants first (kernel profile + per-length top-form invariants).
    kernel_profile, length_summary, length_invariants = _p2_full_sector_invariants(
        N, Omega, F_sec, deg_q, max_exp0, m, p
    )

    # Steps 3-4: type-by-invariant extraction.
    blocks, built_cols_sec, blocks_meta, extraction_log = _p2_extract_blocks(
        F_sec, N, Omega, T_sec, m, max_exp0, deg_q, key, p
    )

    # Step 4 (dagger): kernel-profile / span consistency.
    kernel_profile_blocks = _p2_verify_kernel_profile(
        built_cols_sec, N, kernel_profile, max_exp0, m, p, key
    )

    # Steps 5-6: cost certificate + assemble.
    inv_data = _p2_build_certificate(
        blocks_meta, kernel_profile, kernel_profile_blocks,
        length_summary, length_invariants, extraction_log, deg_q, max_exp0
    )
    inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data=inv_data)
    blocks = [AtomicBlock(b.T_blk, b.half_dim, key, inv) for b in blocks]
    return blocks, inv
