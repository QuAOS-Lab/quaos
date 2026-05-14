# sympleq/core/symmetries/atomic_unipotent_p2.py
from __future__ import annotations

import itertools
from typing import Any, Dict, List, Tuple

import numpy as np

from ..modular_helpers import mod_p, independent_columns, rank_mod, omega_matrix
from .atomic_types import AtomicBlock, AtomicInvariant
from .atomic_linear import (
    restrict_operator,
    is_nondegenerate,
    darboux_basis_from_span,
    mat_pow_mod,
    kernel_in_span,
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


# ---------------------------------------------------------------------------
# Restricted nilpotent top spaces
# ---------------------------------------------------------------------------

def _basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    """Deterministically pick ``want`` columns from ``candidates`` extending ``span(base)``."""
    base = independent_columns(mod_p(base, p), p) if base.size else base
    picked = np.zeros((candidates.shape[0], 0), dtype=np.int64)
    r_base = rank_mod(base, p) if base.size else 0

    for j in range(candidates.shape[1]):
        c = mod_p(candidates[:, j:j + 1], p)
        r_try = rank_mod(np.concatenate([base, picked, c], axis=1), p)
        if r_try > r_base + picked.shape[1]:
            picked = np.concatenate([picked, c], axis=1)
            if picked.shape[1] == want:
                return picked

    raise RuntimeError("_basis_extend: could not extend by required amount.")


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
    Deterministic candidate top pairs (c,d) with c^T Btop d = 1.

    When Btop is alternating, try quotient-level Gram-Schmidt pairs first.
    For non-alternating top forms (which occur for V_beta(2k) single-chain
    blocks), skip Gram-Schmidt and use exhaustive coefficient pairs.  Direct
    cyclic-module nondegeneracy checks decide whether a candidate is accepted.
    """
    Btop = mod_p(np.asarray(Btop, dtype=np.int64), 2)
    m = int(Btop.shape[0])

    e_list: List[np.ndarray] = []
    f_list: List[np.ndarray] = []
    if m > 0 and np.all(np.diag(Btop) % 2 == 0):
        try:
            e_list, f_list, _ = _p2_hyperbolic_pairs_from_alternating_form(Btop, p=2)
        except RuntimeError:
            # Fall through to exhaustive coefficient-pair enumeration.  The
            # acceptor verifies every lifted block, so this is safe.
            e_list, f_list = [], []

    out: List[Tuple[np.ndarray, np.ndarray]] = []
    seen: set[Tuple[Tuple[int, ...], Tuple[int, ...]]] = set()

    def pair_value(c: np.ndarray, d: np.ndarray) -> int:
        return _scalar_mod2(c.reshape(1, -1) @ Btop @ d.reshape(-1, 1))

    def push(c: np.ndarray, d: np.ndarray) -> None:
        cc = mod_p(np.asarray(c, dtype=np.int64).reshape(-1), 2)
        dd = mod_p(np.asarray(d, dtype=np.int64).reshape(-1), 2)
        if cc.shape[0] != m or dd.shape[0] != m:
            return
        if not np.any(cc) or not np.any(dd):
            return
        if pair_value(cc, dd) != 1:
            return
        key = (tuple(int(x) for x in cc), tuple(int(x) for x in dd))
        if key not in seen:
            seen.add(key)
            out.append((cc, dd))

    for c, d in zip(e_list, f_list):
        push(c, d)

    coeffs = _top_coeff_candidates(m, exhaustive_limit=8)
    for c in coeffs:
        for d in coeffs:
            push(c, d)

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

def atomic_blocks_in_unipotent_self_sector_p2(
    F: np.ndarray,
    key: Tuple[int, ...],
    primaries: dict,
    *,
    allow_fallback: bool = False,
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    p=2, self sector with q(x)=x±1.

    Certified policy:
      * decompose the current length-top quotient by its p=2 top bilinear
        pairing b_L and mid-chain quadratic labels q_L;
      * lift quotient pieces to cyclic modules;
      * accept a block only after direct nondegeneracy/invariance checks;
      * never drop chain length in certified mode.

    ``allow_fallback=True`` permits one last deterministic candidate sweep at the
    current length before failing.  The global certified wrapper calls this with
    the default ``False``.
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
        except Exception:
            gen_dim = -1
        length_summary[int(L)] = {
            "mult": int(invL["top_dim"]),
            "gen_dim": int(gen_dim),
            "gen_dim_ok": bool(gen_dim == int(invL["top_dim"])),
            **invL,
        }

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

        # Invariance guardrail: F and N restrictions must be well-defined by restrict_operator.
        try:
            _ = restrict_operator(F_sec, span, p)
            _ = restrict_operator(N, span, p)
        except Exception:
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
        Btop = _top_pairing_matrix(A_top, Omega, N, int(L))

        # Only alternating top pairings have a quotient-level radical for the
        # bilinear Gram-Schmidt helper.  V_beta(2k) one-line tops can have
        # a nonzero diagonal in Btop; in that case the whole point is to try
        # the top line as an anisotropic single-chain candidate below, not to
        # reject the length before extraction starts.
        rad_list: List[np.ndarray] = []
        if Btop.shape[0] > 0 and np.all(np.diag(Btop) % 2 == 0):
            try:
                _e_list, _f_list, rad_list = _p2_hyperbolic_pairs_from_alternating_form(Btop, p=2)
            except RuntimeError:
                rad_list = []

        # V_beta(2k): a single cyclic even-length block may appear as an
        # anisotropic one-dimensional top line, so it is NOT always in the
        # radical of b_L.  Try all quotient-top lines deterministically,
        # preferring q=1 witnesses, and let direct cyclic-module
        # nondegeneracy/invariance checks decide acceptance.
        if int(L) % 2 == 0:
            # For non-alternating symmetric top forms, first diagonalize the
            # anisotropic part.  This prevents the greedy extractor from
            # prematurely grouping diagonalizable V_beta(2k) lines into W pairs.
            preferred = _anisotropic_orthogonal_line_candidates(Btop)
            coeffs_all = _top_coeff_candidates(A_top.shape[1], exhaustive_limit=10)

            seen_coeffs: set[Tuple[int, ...]] = set()
            coeffs: List[np.ndarray] = []
            for c0 in preferred + coeffs_all:
                cc = mod_p(np.asarray(c0, dtype=np.int64).reshape(-1), 2)
                keyc = tuple(int(x) for x in cc)
                if keyc not in seen_coeffs:
                    seen_coeffs.add(keyc)
                    coeffs.append(cc)

            def _top_norm(c: np.ndarray) -> int:
                return _scalar_mod2(c.reshape(1, -1) @ Btop @ c.reshape(-1, 1))

            coeffs = sorted(
                coeffs,
                # Prefer anisotropic top lines for V_beta extraction; q-value is
                # only a secondary Chapter-5 label and may vanish for all lines.
                key=lambda c: (1 - _top_norm(c), 1 - _q_value_from_coeff(c, A_top, Omega, N, int(L)), tuple(int(x) for x in c)),
            )
            for c in coeffs:
                v_top = mod_p(A_top @ c.reshape(-1, 1), p)
                try:
                    Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, int(L), p)
                except Exception:
                    continue
                Cv = independent_columns(Cv, p)
                if accept_block(Cv, v_top, int(L), "V"):
                    return True

            # Keep radical-adapted candidates as a redundant fallback for large
            # top spaces where exhaustive coefficient enumeration is capped.
            for c in _radical_candidates_adapted_to_q(rad_list, A_top, Omega, N, int(L)):
                v_top = mod_p(A_top @ c.reshape(-1, 1), p)
                try:
                    Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, int(L), p)
                except Exception:
                    continue
                Cv = independent_columns(Cv, p)
                if accept_block(Cv, v_top, int(L), "V"):
                    return True

        # W(k): choose dual/hyperbolic top vectors at quotient level.  We try
        # the canonical Gram-Schmidt pairs first and, for small quotient
        # dimensions, all coefficient pairs with c^T B d = 1.
        for ce, cf in _hyperbolic_pair_candidates(Btop):
            v_top = mod_p(A_top @ ce.reshape(-1, 1), p)
            w_top = mod_p(A_top @ cf.reshape(-1, 1), p)
            try:
                Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, int(L), p)
                Cw = cyclic_submodule_basis(F_sec, N, w_top, deg_q, int(L), p)
            except Exception:
                continue
            span = independent_columns(np.concatenate([Cv, Cw], axis=1), p)
            if accept_block(span, np.concatenate([v_top, w_top], axis=1), int(L), "W"):
                return True

        return False

    def try_best_effort_candidates(A_top: np.ndarray, L: int) -> bool:
        """Last-resort deterministic sweep used only outside certified mode."""
        if not allow_fallback:
            return False
        A_gen = _select_module_generators_from_top_space(F_sec, N, A_top, deg_q, int(L), p)
        cols = [A_gen[:, j:j + 1] for j in range(A_gen.shape[1])]
        if int(L) % 2 == 0:
            for v_top in cols:
                try:
                    Cv = cyclic_submodule_basis(F_sec, N, v_top, deg_q, int(L), p)
                except Exception:
                    continue
                if accept_block(Cv, v_top, int(L), "V"):
                    return True
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                try:
                    Cv = cyclic_submodule_basis(F_sec, N, cols[i], deg_q, int(L), p)
                    Cw = cyclic_submodule_basis(F_sec, N, cols[j], deg_q, int(L), p)
                except Exception:
                    continue
                span = independent_columns(np.concatenate([Cv, Cw], axis=1), p)
                if accept_block(span, np.concatenate([cols[i], cols[j]], axis=1), int(L), "W"):
                    return True
        return False

    while space_basis.shape[1] > 0:
        iteration += 1
        if iteration > 2 * m + 5:
            raise RuntimeError("Unipotent p=2 self sector: extraction exceeded iteration guard.")

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
        if not progressed:
            progressed = try_best_effort_candidates(A_top, int(L))
            log_entry["used_best_effort_sweep"] = bool(progressed)

        log_entry["progressed"] = bool(progressed)
        extraction_log.append(log_entry)

        if not progressed:
            raise RuntimeError(
                "Unipotent p=2 self sector: certified quotient-driven extraction made no progress "
                f"at length L={int(L)}. Diagnostics={log_entry}"
            )

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
    used_best_effort = bool(any(x.get("used_best_effort_sweep", False) for x in extraction_log))
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

    inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data=inv_data)
    blocks = [AtomicBlock(b.T_blk, b.half_dim, key, inv) for b in blocks]
    return blocks, inv
