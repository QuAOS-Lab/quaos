# sympleq/core/symmetries/atomic_unipotent_p2.py
from __future__ import annotations
import numpy as np
from typing import List, Tuple, Dict

from sympleq.core.symmetries.modular_helpers import (
    mod_p, independent_columns, rank_mod, _solve_linear, omega_matrix
)
from sympleq.core.symmetries.atomic_types import AtomicBlock, AtomicInvariant
from sympleq.core.symmetries.atomic_linear import (
    restrict_operator,
    is_nondegenerate, darboux_basis_from_span,
    mat_pow_mod, kernel_in_span,
    symplectic_orthogonal_complement_in_span
)
from sympleq.core.symmetries.atomic_krylov import _col, krylov_chain


# ---------------------------
# Public entry points
# ---------------------------

def classify_unipotent_sp2(F_u: np.ndarray) -> dict:
    """
    Deterministic, sector-internal invariant summary for unipotent elements in Sp(2m,2).

    IMPORTANT:
      - This assumes F_u is expressed in a *canonical* symplectic basis, i.e.
            F_u^T Ω F_u = Ω
        with Ω = omega_matrix(m,2).

    What we return:
      - Jordan profile of N = F_u - I.
      - “Even-layer anisotropy” bits: for each even m that occurs, whether ∃ v in ker N^m \\ ker N^(m-1)
        with  <v, N^(m-1) v> = 1  (a useful structural distinguisher in char 2).
      - A canonical block list produced by the same deterministic extraction that we use for construction.

    This is strong enough to:
      - certify the extracted atomic block multiset for your decomposition pipeline,
      - and serve as a reproducible “fingerprint” for regression tests.
    """
    p = 2
    F_u = mod_p(F_u, p)
    n2 = F_u.shape[0]
    if n2 % 2 != 0 or F_u.shape[1] != n2:
        raise ValueError("classify_unipotent_sp2: F_u must be square with even dimension.")
    Omega = omega_matrix(n2 // 2, p)

    N = unipotent_N(F_u, p)
    prof = jordan_profile_nilpotent(N, p)

    # Even-layer anisotropy bits (helpful in char 2; cheap-ish to compute)
    anisotropic_even_layers: Dict[int, int] = {}
    Id = np.eye(n2, dtype=np.int64)
    max_size = int(prof.get("max_size", 0))
    for m in range(2, max_size + 1):
        if m % 2 != 0:
            continue
        # Only bother if blocks of size >= m exist
        blocks_ge = prof.get("blocks_ge", [])
        if len(blocks_ge) < m or int(blocks_ge[m - 1]) == 0:
            continue
        anisotropic_even_layers[m] = int(_exists_anisotropic_in_layer(N, Omega, m, Id, p))

    # Canonical block list from deterministic extraction in coordinate space
    blocks_meta: List[dict] = []
    space_basis = Id
    while space_basis.shape[1] > 0:
        _, space_basis, meta = extract_one_unipotent_block_char2(F_u, Omega, space_basis, p=2)
        blocks_meta.append(meta)

    return {
        "p": 2,
        "n2": int(n2),
        "jordan_profile": prof,
        "anisotropic_even_layers": anisotropic_even_layers,
        "blocks": blocks_meta,
    }


def build_unipotent_blocks_from_invariants(
    F: np.ndarray,
    V_u: np.ndarray,
    inv_data: dict,
    p: int = 2
) -> List[AtomicBlock]:
    """
    Construct explicit atomic unipotent blocks in AMBIENT coordinates.

    Inputs:
      F       : ambient symplectic matrix (2nx2n)
      V_u     : ambient columns spanning the unipotent sector subspace
      inv_data: output of classify_unipotent_sp2() computed on the restricted operator
               in a canonical symplectic basis of this sector.

    Output:
      list[AtomicBlock] where each block basis T_blk is ambient and Darboux:
          T_blk^T Omega_amb T_blk = Omega_block
    """
    if p != 2:
        raise ValueError("build_unipotent_blocks_from_invariants: this file is p=2 only.")
    F = mod_p(F, 2)
    V_u = independent_columns(mod_p(V_u, 2), 2)
    if V_u.shape[1] == 0:
        return []

    n2 = F.shape[0]
    Omega_amb = omega_matrix(n2 // 2, 2)

    # Canonicalize the unipotent sector basis so the restricted Ω is standard.
    # This is CRUCIAL: CRT projector bases are not guaranteed Darboux.
    if not is_nondegenerate(Omega_amb, V_u, 2):
        raise RuntimeError("Unipotent sector basis is degenerate; cannot build symplectic blocks.")
    T_u = darboux_basis_from_span(Omega_amb, V_u, 2)   # ambient (2n × 2m)
    m2 = T_u.shape[1]
    if m2 % 2 != 0:
        raise RuntimeError("Unipotent sector dimension must be even.")
    m = m2 // 2

    # Restrict operator into canonical sector coordinates
    F_u = restrict_operator(F, T_u, 2)
    Omega_u = omega_matrix(m, 2)

    # Extract blocks in coordinate space, then lift them to ambient via T_u
    blocks: List[AtomicBlock] = []
    metas: List[dict] = []

    space_basis = np.eye(2 * m, dtype=np.int64)
    while space_basis.shape[1] > 0:
        T_blk_u, space_basis, meta = extract_one_unipotent_block_char2(F_u, Omega_u, space_basis, p=2)
        metas.append(meta)

        # Lift to ambient
        T_blk_amb = mod_p(T_u @ T_blk_u, 2)
        half_dim = T_blk_amb.shape[1] // 2
        blocks.append(AtomicBlock(T_blk=T_blk_amb, half_dim=int(half_dim), sector_key=(), inv=None))

    # Optional consistency check: extracted meta list matches inv_data “blocks”
    if isinstance(inv_data, dict) and "blocks" in inv_data:
        want = inv_data["blocks"]
        if want != metas:
            raise RuntimeError(
                "Unipotent reconstruction mismatch: extracted block meta differs from inv_data['blocks'].\n"
                f"want={want}\n"
                f"got ={metas}"
            )

    return blocks


def atomic_blocks_in_unipotent_self_sector_p2(
    F: np.ndarray, key: Tuple[int, ...], primaries: dict
) -> Tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Entry point for p=2, q=x+1 self sector.

    This routine:
      1) takes the ambient sector span primaries[key]["V_basis"],
      2) canonicalizes it into a Darboux basis T_u,
      3) restricts F into that canonical basis,
      4) classifies (invariants),
      5) reconstructs explicit blocks in ambient coordinates.
    """
    p = 2
    F = mod_p(F, p)
    V_u = independent_columns(mod_p(primaries[key]["V_basis"], p), p)
    if V_u.shape[1] == 0:
        inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data={"p": 2, "n2": 0, "blocks": []})
        return [], inv

    # Build invariants in canonical sector coords
    n2 = F.shape[0]
    Ω_amb = omega_matrix(n2 // 2, p)
    if not is_nondegenerate(Ω_amb, V_u, p):
        raise RuntimeError("Unipotent primary component is degenerate; cannot proceed.")
    T_u = darboux_basis_from_span(Ω_amb, V_u, p)         # ambient Darboux basis for this sector
    F_u = restrict_operator(F, T_u, p)                   # sector coords with standard Ω
    inv_data = classify_unipotent_sp2(F_u)

    inv = AtomicInvariant(sector_key=key, sector_type="self", poly_key=key, data=inv_data)

    # Rebuild blocks in ambient coordinates; pass T_u as the sector span (also fine)
    blocks = build_unipotent_blocks_from_invariants(F, T_u, inv_data, p=2)
    # Fill in metadata expected by your pipeline
    out = [AtomicBlock(b.T_blk, b.half_dim, key, inv) for b in blocks]
    return out, inv


# ---------------------------
# Core helpers (unchanged logic + small robustness tweaks)
# ---------------------------

def unipotent_N(F: np.ndarray, p: int) -> np.ndarray:
    """Return N = F - I mod p. (For p=2, this equals F + I.)"""
    n = F.shape[0]
    return mod_p(F - np.eye(n, dtype=np.int64), p)


def jordan_profile_nilpotent(N: np.ndarray, p: int) -> dict:
    """
    Jordan profile from kernel dimensions of N^r.

    Returns dict with:
      - ker_dims[r] = dim ker N^(r+1)  (0-indexed list)
      - blocks_ge[r] = # Jordan blocks of size >= (r+1)
      - blocks_exact[r] = # blocks of size exactly (r+1)
      - max_size = largest Jordan size
    """
    N = mod_p(N, p)
    n = N.shape[0]

    ker_dims: List[int] = []
    P = np.eye(n, dtype=np.int64)
    for r in range(1, n + 1):
        P = mod_p(P @ N, p)   # P = N^r
        rk = rank_mod(P, p)
        ker_dims.append(n - rk)
        if ker_dims[-1] == n:
            break

    blocks_ge: List[int] = []
    prev = 0
    for d in ker_dims:
        blocks_ge.append(d - prev)
        prev = d

    blocks_exact: List[int] = []
    for i in range(len(blocks_ge)):
        nxt = blocks_ge[i + 1] if (i + 1) < len(blocks_ge) else 0
        blocks_exact.append(blocks_ge[i] - nxt)

    max_size = 0
    for i, c in enumerate(blocks_exact):
        if c > 0:
            max_size = i + 1

    return {
        "ker_dims": ker_dims,
        "blocks_ge": blocks_ge,
        "blocks_exact": blocks_exact,
        "max_size": max_size,
    }


def jordan_profile_nilpotent_in_span(N: np.ndarray, span_basis: np.ndarray, p: int) -> dict:
    """
    Jordan profile for N restricted to span(span_basis), computed via kernel_in_span(N^r,...).
    """
    N = mod_p(N, p)
    span_basis = independent_columns(mod_p(span_basis, p), p)
    d = span_basis.shape[1]
    if d == 0:
        return {"ker_dims": [], "blocks_ge": [], "blocks_exact": [], "max_size": 0, "dim": 0}

    ker_dims: List[int] = []
    Nr = np.eye(N.shape[0], dtype=np.int64)
    for r in range(1, N.shape[0] + 1):
        Nr = mod_p(Nr @ N, p)
        K = kernel_in_span(Nr, span_basis, p)
        ker_dims.append(K.shape[1])
        if ker_dims[-1] == d:
            break

    blocks_ge: List[int] = []
    prev = 0
    for dd in ker_dims:
        blocks_ge.append(dd - prev)
        prev = dd

    blocks_exact: List[int] = []
    for i in range(len(blocks_ge)):
        nxt = blocks_ge[i + 1] if (i + 1) < len(blocks_ge) else 0
        blocks_exact.append(blocks_ge[i] - nxt)

    max_size = 0
    for i, c in enumerate(blocks_exact):
        if c > 0:
            max_size = i + 1

    return {
        "ker_dims": ker_dims,
        "blocks_ge": blocks_ge,
        "blocks_exact": blocks_exact,
        "max_size": max_size,
        "dim": d,
    }


def _exists_anisotropic_in_layer(N: np.ndarray, Omega: np.ndarray, m: int,
                                 space_basis: np.ndarray, p: int = 2) -> bool:
    """
    Test if there exists v in ker(N^m) \\ ker(N^(m-1)) within span(space_basis)
    such that q(v) := <v, N^(m-1) v> = 1.

    This is a robust “layer” diagnostic that your extraction uses implicitly for even m.
    """
    if p != 2:
        raise ValueError("_exists_anisotropic_in_layer: p=2 only.")
    m = int(m)
    if m <= 0:
        return False

    space_basis = independent_columns(mod_p(space_basis, p), p)
    if space_basis.shape[1] == 0:
        return False

    Nm = mat_pow_mod(N, m, p)
    Nm_1 = mat_pow_mod(N, m - 1, p)

    KerNm = kernel_in_span(Nm, space_basis, p)
    if KerNm.shape[1] == 0:
        return False

    def q(v: np.ndarray) -> int:
        return int(mod_p(v.T @ Omega @ (Nm_1 @ v), p).reshape(()))

    # scan a spanning set for ker(N^m)\ker(N^(m-1))
    for j in range(KerNm.shape[1]):
        v = KerNm[:, j:j + 1]
        if np.all(mod_p(Nm_1 @ v, p) == 0):
            continue
        if q(v) == 1:
            return True

    # pairwise sums can reveal anisotropy when basis vectors don't
    for i in range(KerNm.shape[1]):
        vi = KerNm[:, i:i + 1]
        if np.all(mod_p(Nm_1 @ vi, p) == 0):
            continue
        for j in range(i + 1, KerNm.shape[1]):
            vj = KerNm[:, j:j + 1]
            if np.all(mod_p(Nm_1 @ vj, p) == 0):
                continue
            v = mod_p(vi + vj, p)
            if np.all(mod_p(Nm_1 @ v, p) == 0):
                continue
            if q(v) == 1:
                return True

    return False


# ---------------------------
# Your constructive block extraction (kept, but now used end-to-end)
# ---------------------------

def find_v_for_V_block_char2(
    N: np.ndarray,
    Omega: np.ndarray,
    k: int,
    space_basis: np.ndarray,
    p: int = 2,
) -> np.ndarray | None:
    """
    Find v such that:
      - N^(2k) v = 0 but N^(2k-1) v != 0,
      - <v, N^(2k-1) v> = 1.
    """
    if p != 2:
        raise ValueError("find_v_for_V_block_char2 is intended for p=2 only.")
    k = int(k)
    if k <= 0:
        return None

    space_basis = independent_columns(mod_p(space_basis, p), p)
    if space_basis.shape[1] == 0:
        return None

    N2k = mat_pow_mod(N, 2 * k, p)
    N2k_1 = mat_pow_mod(N, 2 * k - 1, p)

    Ker2k = kernel_in_span(N2k, space_basis, p)
    if Ker2k.shape[1] == 0:
        return None

    candidates: List[np.ndarray] = []
    for j in range(Ker2k.shape[1]):
        v = Ker2k[:, j:j + 1]
        if np.all(mod_p(N2k_1 @ v, p) == 0):
            continue
        candidates.append(v)

    if not candidates:
        return None

    def q(v: np.ndarray) -> int:
        return int(mod_p(v.T @ Omega @ (N2k_1 @ v), p).reshape(()))

    for v in candidates:
        if q(v) == 1:
            return mod_p(v, p)

    for i in range(len(candidates)):
        for j in range(i + 1, len(candidates)):
            v = mod_p(candidates[i] + candidates[j], p)
            if np.all(mod_p(N2k_1 @ v, p) == 0):
                continue
            if q(v) == 1:
                return v

    return None


def build_V_block_subspace_char2(
    N: np.ndarray,
    Omega: np.ndarray,
    v: np.ndarray,
    k: int,
    p: int = 2,
) -> np.ndarray:
    """
    Build a V(2k)-type invariant subspace from v via chain:
        span{ v, Nv, ..., N^(2k-1) v }.
    Returns Darboux basis T_blk (n×2k).
    """
    if p != 2:
        raise ValueError("build_V_block_subspace_char2 is intended for p=2 only.")
    v = _col(v)
    k = int(k)

    chain = krylov_chain(N, v, 2 * k, p)
    X = independent_columns(chain, p)
    if X.shape[1] != 2 * k:
        raise RuntimeError("build_V_block_subspace_char2: chain did not have expected dimension.")
    if not is_nondegenerate(Omega, X, p):
        raise RuntimeError("build_V_block_subspace_char2: constructed span is degenerate (not V-type).")

    return darboux_basis_from_span(Omega, X, p)


def find_wx_for_W_block_char2(
    N: np.ndarray,
    Omega: np.ndarray,
    m: int,
    space_basis: np.ndarray,
    p: int = 2,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find w, x such that:
      - N^m w = 0 but N^(m-1) w != 0,
      - x ∈ ker(N^m), and <N^b w, x> = δ_{b,m-1},
      - and x has length m (enforced by a witness constraint).
    """
    if p != 2:
        raise ValueError("find_wx_for_W_block_char2 is intended for p=2 only.")
    m = int(m)
    space_basis = independent_columns(mod_p(space_basis, p), p)
    if space_basis.shape[1] == 0:
        raise RuntimeError("find_wx_for_W_block_char2: empty space_basis.")

    Nm = mat_pow_mod(N, m, p)
    Nm_1 = mat_pow_mod(N, m - 1, p)

    KerNm = kernel_in_span(Nm, space_basis, p)
    if KerNm.shape[1] == 0:
        raise RuntimeError("find_wx_for_W_block_char2: ker(N^m) is empty in this space.")

    w: np.ndarray | None = None
    for j in range(KerNm.shape[1]):
        cand = KerNm[:, j:j + 1]
        if np.all(mod_p(Nm_1 @ cand, p) == 0):
            continue
        w = cand
        break
    if w is None:
        raise RuntimeError("find_wx_for_W_block_char2: could not find w of chain length m.")

    K = krylov_chain(N, w, m, p)

    A = mod_p(K.T @ Omega @ KerNm, p)  # m×d
    b = np.zeros((m, 1), dtype=np.int64)
    b[-1, 0] = 1

    P = mod_p(Nm_1 @ KerNm, p)  # n×d

    for coord in range(P.shape[0]):
        row = P[coord: coord + 1, :]
        if np.all(row % p == 0):
            continue
        A2 = np.concatenate([A, row], axis=0)
        b2 = np.concatenate([b, np.array([[1]], dtype=np.int64)], axis=0)
        try:
            coeff = _solve_linear(A2, b2, p)
        except RuntimeError:
            continue
        x = mod_p(KerNm @ coeff, p)
        if np.all(mod_p(Nm_1 @ x, p) == 0):
            continue
        return mod_p(w, p), x

    coeff0 = _solve_linear(A, b, p)
    x0 = mod_p(KerNm @ coeff0, p)
    return mod_p(w, p), x0


def build_W_block_subspace_char2(
    N: np.ndarray,
    Omega: np.ndarray,
    w: np.ndarray,
    x: np.ndarray,
    m: int,
    p: int = 2,
) -> np.ndarray:
    """
    Build W(m)-type invariant subspace:
      span{ N^i w, N^i x : i=0..m-1 }.
    Returns Darboux basis T_blk (n×2m).
    """
    if p != 2:
        raise ValueError("build_W_block_subspace_char2 is intended for p=2 only.")
    w = _col(w)
    x = _col(x)
    m = int(m)

    W1 = krylov_chain(N, w, m, p)
    W2 = krylov_chain(N, x, m, p)
    X = independent_columns(np.concatenate([W1, W2], axis=1), p)
    if X.shape[1] != 2 * m:
        raise RuntimeError("build_W_block_subspace_char2: span did not have expected dimension 2m.")
    if not is_nondegenerate(Omega, X, p):
        raise RuntimeError("build_W_block_subspace_char2: constructed span is degenerate.")
    return darboux_basis_from_span(Omega, X, p)


def extract_one_unipotent_block_char2(
    F: np.ndarray,
    Omega: np.ndarray,
    space_basis: np.ndarray,
    p: int = 2,
) -> Tuple[np.ndarray, np.ndarray, dict]:
    """
    Extract one invariant nondegenerate unipotent block inside current space (p=2).

    Returns:
      T_blk    : Darboux basis for extracted block (ambient coords)
      rem      : basis for (current_space ∩ span(T_blk)^⊥)
      meta     : dict describing the block
    """
    if p != 2:
        raise ValueError("extract_one_unipotent_block_char2 is intended for p=2 only.")
    space_basis = independent_columns(mod_p(space_basis, p), p)
    if space_basis.shape[1] == 0:
        raise RuntimeError("extract_one_unipotent_block_char2: empty space_basis.")

    N = unipotent_N(F, p)
    prof = jordan_profile_nilpotent_in_span(N, space_basis, p)
    m_max = int(prof["max_size"])
    if m_max <= 0:
        raise RuntimeError("extract_one_unipotent_block_char2: no nilpotent action detected in space.")

    # Try “V-type” when m_max is even
    if (m_max % 2) == 0:
        k = m_max // 2
        v = find_v_for_V_block_char2(N, Omega, k, space_basis, p)
        if v is not None:
            T_blk = build_V_block_subspace_char2(N, Omega, v, k, p)
            rem = symplectic_orthogonal_complement_in_span(Omega, T_blk, space_basis, p)
            return T_blk, rem, {"type": "V", "m": m_max, "k": k, "jordan_max": m_max}

    # Otherwise build W-block of length m_max
    w, x = find_wx_for_W_block_char2(N, Omega, m_max, space_basis, p)
    T_blk = build_W_block_subspace_char2(N, Omega, w, x, m_max, p)
    rem = symplectic_orthogonal_complement_in_span(Omega, T_blk, space_basis, p)
    return T_blk, rem, {"type": "W", "m": m_max, "jordan_max": m_max}
