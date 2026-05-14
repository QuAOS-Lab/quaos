# sympleq/core/symmetries/atomic_decomposition_helpers/atomic_self_p2_unitary.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..modular_helpers import mod_p, independent_columns, inv_mod_mat, omega_matrix, rank_mod, nullspace_mod
from .atomic_types import AtomicBlock, AtomicInvariant
from .atomic_linear import (
    darboux_basis_from_span,
    is_nondegenerate,
    restrict_operator,
    symplectic_orthogonal_complement_in_span,
)
from .module_invariants import cyclic_submodule_basis


# ---------------------------------------------------------------------------
# Small GF(2^d) implementation
# ---------------------------------------------------------------------------


class GF2Extension:
    """
    GF(2)[x]/(modulus), represented by integers whose binary expansion stores
    coefficients in the power basis 1, alpha, ..., alpha^{d-1}.

    This class is deliberately small and local to the p=2 self-reciprocal
    non-unipotent sector.  It provides exactly the arithmetic needed for the
    Hermitian top-space decomposition.
    """

    def __init__(self, modulus_coeffs: Sequence[int]):
        q = [int(c) & 1 for c in modulus_coeffs]
        while len(q) > 1 and q[-1] == 0:
            q.pop()
        if len(q) < 3 or q[-1] != 1:
            raise ValueError(f"Expected a monic irreducible polynomial of degree >1 over GF(2), got {modulus_coeffs!r}")
        self.modulus_coeffs = tuple(q)
        self.d = len(q) - 1
        self.mask = (1 << self.d) - 1

        # Integer with bits for q(x) without leading x^d term.  During reduction,
        # x^d == q_low(x) because subtraction is addition in characteristic two.
        self._q_low = 0
        for i, c in enumerate(q[:-1]):
            if c & 1:
                self._q_low |= 1 << i

        self.alpha = 2 if self.d > 1 else 1
        self._trace_matrix_inv: Optional[np.ndarray] = None

    def add(self, a: int, b: int) -> int:
        return (int(a) ^ int(b)) & self.mask

    sub = add

    def reduce(self, x: int) -> int:
        x = int(x)
        # Reduce from high degree down to d.
        while x.bit_length() > self.d:
            k = x.bit_length() - 1
            shift = k - self.d
            x ^= 1 << k
            x ^= self._q_low << shift
        return x & self.mask

    def mul(self, a: int, b: int) -> int:
        a = int(a) & self.mask
        b = int(b) & self.mask
        out = 0
        aa = a
        bb = b
        while bb:
            if bb & 1:
                out ^= aa
            aa <<= 1
            bb >>= 1
        return self.reduce(out)

    def pow(self, a: int, e: int) -> int:
        a = int(a) & self.mask
        e = int(e)
        out = 1
        base = a
        while e > 0:
            if e & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            e >>= 1
        return out

    def inv(self, a: int) -> int:
        a = int(a) & self.mask
        if a == 0:
            raise ZeroDivisionError("inverse of zero in GF(2^d)")
        # Multiplicative group has order 2^d - 1.
        return self.pow(a, (1 << self.d) - 2)

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.inv(b))

    def conj(self, a: int) -> int:
        """
        The involution alpha -> alpha^{-1}.  For an irreducible self-reciprocal
        polynomial of even degree d over GF(2), this equals Frobenius 2^(d/2).
        """
        if self.d % 2 != 0:
            raise ValueError("self-reciprocal non-linear irreducibles over GF(2) should have even degree")
        return self.pow(a, 1 << (self.d // 2))

    def trace_to_F2(self, a: int) -> int:
        a = int(a) & self.mask
        acc = 0
        x = a
        for _ in range(self.d):
            acc ^= x
            x = self.mul(x, x)
        # The trace lies in GF(2), represented as 0 or 1.
        return acc & 1

    def coeff_vector(self, a: int) -> np.ndarray:
        return np.array([(int(a) >> i) & 1 for i in range(self.d)], dtype=np.int64).reshape(self.d, 1)

    def from_coeff_vector(self, c: np.ndarray) -> int:
        cc = np.asarray(c, dtype=np.int64).reshape(-1) % 2
        if cc.size != self.d:
            raise ValueError(f"expected coefficient vector of length {self.d}, got {cc.size}")
        out = 0
        for i, bit in enumerate(cc):
            if int(bit) & 1:
                out |= 1 << i
        return out & self.mask

    def trace_pairing_inverse(self) -> np.ndarray:
        """
        Inverse of the matrix M_{a,k} = Tr(alpha^a alpha^k), 0<=a,k<d.
        It maps trace coordinates to power-basis coefficients.
        """
        if self._trace_matrix_inv is not None:
            return self._trace_matrix_inv
        M = np.zeros((self.d, self.d), dtype=np.int64)
        powers = [1]
        for _ in range(1, 2 * self.d):
            powers.append(self.mul(powers[-1], self.alpha))
        for a in range(self.d):
            for k in range(self.d):
                M[a, k] = self.trace_to_F2(powers[a + k])
        if rank_mod(M, 2) != self.d:
            raise RuntimeError("trace pairing matrix is singular; polynomial may not define a field")
        self._trace_matrix_inv = inv_mod_mat(M, 2)
        return self._trace_matrix_inv

    def from_trace_values(self, values: Sequence[int]) -> int:
        """
        Given b_a = Tr(alpha^a h), recover h in the power basis.
        """
        b = np.asarray(values, dtype=np.int64).reshape(self.d, 1) % 2
        c = mod_p(self.trace_pairing_inverse() @ b, 2)
        return self.from_coeff_vector(c)


# ---------------------------------------------------------------------------
# Basic base-field linear algebra helpers
# ---------------------------------------------------------------------------


def _mat_pow_mod(A: np.ndarray, e: int, p: int) -> np.ndarray:
    A = mod_p(np.asarray(A, dtype=np.int64), p)
    n = A.shape[0]
    R = np.eye(n, dtype=np.int64)
    B = A.copy()
    ee = int(e)
    while ee > 0:
        if ee & 1:
            R = mod_p(R @ B, p)
        B = mod_p(B @ B, p)
        ee >>= 1
    return R


def _kernel(A: np.ndarray, p: int) -> np.ndarray:
    return independent_columns(nullspace_mod(mod_p(A, p), p), p)


def _basis_extend(base: np.ndarray, candidates: np.ndarray, want: int, p: int) -> np.ndarray:
    base = independent_columns(mod_p(base, p), p) if base.size else base
    candidates = independent_columns(mod_p(candidates, p), p)
    picked = np.zeros((candidates.shape[0], 0), dtype=np.int64)
    r_base = rank_mod(base, p) if base.size else 0
    for j in range(candidates.shape[1]):
        c = candidates[:, j:j + 1]
        trial = np.concatenate([base, picked, c], axis=1) if base.size or picked.size else c
        if rank_mod(trial, p) > r_base + picked.shape[1]:
            picked = np.concatenate([picked, c], axis=1) if picked.size else c
            if picked.shape[1] == want:
                return mod_p(picked, p)
    raise RuntimeError("_basis_extend: could not extend by required amount")


@dataclass
class TopQuotient:
    L: int
    K_L: np.ndarray
    denom: np.ndarray
    top_reps_f2: np.ndarray
    e_basis_reps: List[np.ndarray]
    e_orbit_basis: np.ndarray
    dim_f2: int
    dim_E: int


def _top_quotient_E_basis(F: np.ndarray, N: np.ndarray, L: int, deg_q: int, p: int) -> TopQuotient:
    """
    Construct an E-basis of T_L = K_L/(K_{L-1}+N K_{L+1}).

    The output e_basis_reps contains base-field representatives v_i.  Their
    F-orbits modulo the quotient denominator form an E-basis.
    """
    if p != 2:
        raise ValueError("_top_quotient_E_basis is currently used only over GF(2)")
    d0 = F.shape[0]
    L = int(L)
    deg_q = int(deg_q)

    K_L = _kernel(_mat_pow_mod(N, L, p), p)
    K_Lm1 = np.zeros((d0, 0), dtype=np.int64) if L <= 1 else _kernel(_mat_pow_mod(N, L - 1, p), p)
    K_Lp1 = _kernel(_mat_pow_mod(N, L + 1, p), p)
    NK_Lp1 = mod_p(N @ K_Lp1, p) if K_Lp1.shape[1] else np.zeros((d0, 0), dtype=np.int64)

    denom = K_Lm1
    if NK_Lp1.shape[1]:
        denom = np.concatenate([denom, NK_Lp1], axis=1) if denom.shape[1] else NK_Lp1
    denom = independent_columns(denom, p) if denom.shape[1] else denom

    rank_K = rank_mod(K_L, p) if K_L.shape[1] else 0
    rank_D = rank_mod(denom, p) if denom.shape[1] else 0
    dim_f2 = int(rank_K - rank_D)
    if dim_f2 < 0:
        raise RuntimeError("top quotient has negative computed dimension")
    if dim_f2 == 0:
        return TopQuotient(
            L=L,
            K_L=K_L,
            denom=denom,
            top_reps_f2=np.zeros((d0, 0), dtype=np.int64),
            e_basis_reps=[],
            e_orbit_basis=np.zeros((d0, 0), dtype=np.int64),
            dim_f2=0,
            dim_E=0,
        )
    if dim_f2 % deg_q != 0:
        raise RuntimeError(f"top quotient dimension {dim_f2} is not divisible by deg_q={deg_q}")

    top_reps = _basis_extend(denom, K_L, dim_f2, p)
    dim_E = dim_f2 // deg_q

    selected: List[np.ndarray] = []
    selected_orbits = np.zeros((d0, 0), dtype=np.int64)
    base = denom
    F_pows = [np.eye(d0, dtype=np.int64)]
    for _ in range(1, deg_q):
        F_pows.append(mod_p(F_pows[-1] @ F, p))

    # Deterministically scan the quotient representatives.  Since q is irreducible,
    # every nonzero quotient vector has a full deg(q)-dimensional E-orbit.
    for j in range(top_reps.shape[1]):
        v = top_reps[:, j:j + 1]
        orbit = np.concatenate([mod_p(P @ v, p) for P in F_pows], axis=1)
        orbit = independent_columns(orbit, p)
        trial = np.concatenate([base, selected_orbits, orbit], axis=1) if (base.shape[1] or selected_orbits.shape[1]) else orbit
        old_rank = (rank_mod(base, p) if base.shape[1] else 0) + selected_orbits.shape[1]
        new_rank = rank_mod(trial, p)
        if new_rank - old_rank == deg_q:
            selected.append(v)
            selected_orbits = np.concatenate([selected_orbits, orbit], axis=1) if selected_orbits.shape[1] else orbit
            if len(selected) == dim_E:
                break

    if len(selected) != dim_E:
        raise RuntimeError(
            f"could not construct E-basis of top quotient: needed {dim_E}, got {len(selected)}; "
            f"dim_f2={dim_f2}, deg_q={deg_q}"
        )

    return TopQuotient(
        L=L,
        K_L=K_L,
        denom=denom,
        top_reps_f2=top_reps,
        e_basis_reps=selected,
        e_orbit_basis=selected_orbits,
        dim_f2=dim_f2,
        dim_E=dim_E,
    )


# ---------------------------------------------------------------------------
# Hermitian top form and decomposition
# ---------------------------------------------------------------------------


def _field_scalar_apply_to_vector(F: np.ndarray, field: GF2Extension, scalar: int, v: np.ndarray, p: int) -> np.ndarray:
    """Return scalar(F) v for scalar in GF(2)[alpha]/q(alpha)."""
    coeff = field.coeff_vector(scalar).reshape(-1)
    out = np.zeros_like(v.reshape(-1, 1), dtype=np.int64)
    w = v.reshape(-1, 1)
    P = np.eye(F.shape[0], dtype=np.int64)
    for a in range(field.d):
        if int(coeff[a]) & 1:
            out = mod_p(out + P @ w, p)
        P = mod_p(P @ F, p)
    return mod_p(out, p)


def _top_vector_from_E_coords(F: np.ndarray, field: GF2Extension, e_reps: Sequence[np.ndarray], coords: Sequence[int], p: int) -> np.ndarray:
    if len(e_reps) != len(coords):
        raise ValueError("coordinate length mismatch")
    out = np.zeros_like(e_reps[0].reshape(-1, 1), dtype=np.int64)
    for v, c in zip(e_reps, coords):
        c = int(c) & field.mask
        if c:
            # The Hermitian top form is built with the first argument sampled
            # along the F-orbit, so the coordinate vector returned by the
            # Hermitian Gram-Schmidt routine is in the conjugate/right-module
            # convention relative to our base-field lift.  Lifting with c
            # directly can select a Hermitian-anisotropic coordinate line whose
            # GF(2)[F]-cyclic module is nevertheless degenerate in the base
            # symplectic form for the standard q=x^2+x+1 fixture.  Applying the
            # involution here aligns the E-coordinate convention with the
            # base-field scalar action.
            out = mod_p(out + _field_scalar_apply_to_vector(F, field, field.conj(c), v, p), p)
    return mod_p(out, p)


def _symp_scalar(Omega: np.ndarray, u: np.ndarray, v: np.ndarray, p: int) -> int:
    val = mod_p(u.reshape(1, -1) @ Omega @ v.reshape(-1, 1), p)
    return int(val.reshape(-1)[0])


def _hermitian_top_matrix(
    F: np.ndarray,
    N: np.ndarray,
    Omega: np.ndarray,
    field: GF2Extension,
    e_reps: Sequence[np.ndarray],
    L: int,
    p: int,
) -> np.ndarray:
    """
    Build H_ij in E from base pairings

        b_a = <F^a e_i, N^{L-1} e_j>,  a=0,...,d-1,

    by solving b_a = Tr(alpha^a H_ij).

    The trace pairing is nondegenerate, so this recovers the unique E element
    representing the induced Hermitian form in the chosen convention.  Direct
    base-field verification of lifted blocks remains the final certificate.
    """
    m = len(e_reps)
    H = np.zeros((m, m), dtype=np.int64)
    Npow = _mat_pow_mod(N, int(L) - 1, p)
    F_pows = [np.eye(F.shape[0], dtype=np.int64)]
    for _ in range(1, field.d):
        F_pows.append(mod_p(F_pows[-1] @ F, p))

    for i, ui in enumerate(e_reps):
        for j, vj in enumerate(e_reps):
            Nv = mod_p(Npow @ vj, p)
            values = [_symp_scalar(Omega, mod_p(P @ ui, p), Nv, p) for P in F_pows]
            H[i, j] = field.from_trace_values(values)

    return H


def _h_pair(field: GF2Extension, H: np.ndarray, x: Sequence[int], y: Sequence[int]) -> int:
    """Hermitian product x^* H y, conjugate-linear in x and linear in y."""
    m = H.shape[0]
    acc = 0
    for i in range(m):
        xi = int(x[i]) & field.mask
        if xi == 0:
            continue
        cxi = field.conj(xi)
        for j in range(m):
            yj = int(y[j]) & field.mask
            if yj == 0 or int(H[i, j]) == 0:
                continue
            acc = field.add(acc, field.mul(field.mul(cxi, int(H[i, j])), yj))
    return acc


def _vec_add_scaled(field: GF2Extension, x: List[int], y: Sequence[int], a: int) -> List[int]:
    """x + y*a, where vectors are right E-coordinate vectors."""
    a = int(a) & field.mask
    if a == 0:
        return [int(z) & field.mask for z in x]
    return [field.add(int(xi), field.mul(int(yi), a)) for xi, yi in zip(x, y)]


def _canonical_E_basis(m: int) -> List[List[int]]:
    out: List[List[int]] = []
    for i in range(m):
        v = [0] * m
        v[i] = 1
        out.append(v)
    return out


def _hermitian_decompose(field: GF2Extension, H: np.ndarray) -> Tuple[List[Tuple[str, List[int], Optional[List[int]]]], Dict[str, Any]]:
    """
    Deterministic Hermitian Gram-Schmidt over E/E0.

    Returns blocks:
      - ("line", u, None) for anisotropic/self-dual lines;
      - ("pair", u, v) for hyperbolic pairs.
    """
    H = np.asarray(H, dtype=np.int64)
    if H.ndim != 2 or H.shape[0] != H.shape[1]:
        raise ValueError(f"H must be square, got {H.shape}")
    m = H.shape[0]

    remaining = _canonical_E_basis(m)
    blocks: List[Tuple[str, List[int], Optional[List[int]]]] = []
    diagnostics: Dict[str, Any] = {
        "dim_E": int(m),
        "n_lines": 0,
        "n_pairs": 0,
        "steps": [],
    }

    def is_zero_vec(v: Sequence[int]) -> bool:
        return all((int(x) & field.mask) == 0 for x in v)

    while remaining:
        # Prefer anisotropic lines because one Hermitian line lifts to the
        # smallest self-dual cyclic GF(2)[F]-module.  It is not enough to scan
        # only the current basis vectors: a hyperbolic Hermitian plane can have
        # an isotropic basis (e,f), while e + c f is anisotropic for a suitable c.
        line_idx: Optional[int] = None
        line_norm = 0
        line_vec: Optional[List[int]] = None
        line_remaining: Optional[List[List[int]]] = None

        # First try the current basis vectors.
        for idx, u in enumerate(remaining):
            nrm = _h_pair(field, H, u, u)
            if nrm != 0:
                line_idx = idx
                line_norm = nrm
                line_vec = u
                line_remaining = [w for k, w in enumerate(remaining) if k != idx]
                break

        # If the current basis is totally isotropic, manufacture an anisotropic
        # line from a coupled pair.  Earlier versions used the closed-form choice
        # c = a^{-1} alpha for w = u + v*c.  That depends on the exact left/right
        # Hermitian coordinate convention used when H is reconstructed from the
        # base symplectic trace pairings.  To avoid convention mistakes, search
        # deterministically over all field coefficients.  This is only O(|E|) per
        # coupled pair and is polynomial in deg(q), unlike the old exponential
        # search over the full GF(2)-top quotient.
        if line_vec is None:
            pair_for_line: Optional[Tuple[int, int, int]] = None
            for i, u in enumerate(remaining):
                for j in range(i + 1, len(remaining)):
                    a = _h_pair(field, H, u, remaining[j])
                    if a != 0:
                        pair_for_line = (i, j, a)
                        break
                if pair_for_line is not None:
                    break

            if pair_for_line is not None:
                i, j, a = pair_for_line
                u = remaining[i]
                v = remaining[j]

                w_line: Optional[List[int]] = None
                nrm = 0

                # Try w = u + v*c first in a canonical coefficient order.
                for c in range(1, field.mask + 1):
                    cand = _vec_add_scaled(field, u, v, c)
                    cand_norm = _h_pair(field, H, cand, cand)
                    if cand_norm != 0:
                        w_line = cand
                        nrm = cand_norm
                        break

                # Defensive fallback for possible opposite coordinate convention:
                # w = u*c + v.  This should rarely be needed, but it keeps the
                # decomposition robust against changes in Hermitian coordinate
                # orientation while preserving deterministic behavior.
                if w_line is None:
                    for c in range(1, field.mask + 1):
                        cand = _vec_add_scaled(field, v, u, c)
                        cand_norm = _h_pair(field, H, cand, cand)
                        if cand_norm != 0:
                            w_line = cand
                            nrm = cand_norm
                            break

                if w_line is None or nrm == 0:
                    raise RuntimeError(
                        "failed to manufacture anisotropic Hermitian line from isotropic pair after exhaustive field-coefficient search; "
                        f"a={int(a)}, field_size={int(field.mask + 1)}, alpha={int(field.alpha)}"
                    )

                line_vec = w_line
                line_norm = nrm
                # Replacing one vector of a coupled pair by the anisotropic
                # combination preserves the span with the other vector still in
                # the remaining list.  Remove only index i.
                line_remaining = [remaining[k] for k in range(len(remaining)) if k != i]

        if line_vec is not None:
            u = line_vec
            inv_norm = field.inv(line_norm)
            new_remaining: List[List[int]] = []
            assert line_remaining is not None
            for w in line_remaining:
                # w <- w + u * (h(u,u)^-1 h(u,w)) to kill h(u,w).
                coeff = field.mul(inv_norm, _h_pair(field, H, u, w))
                ww = _vec_add_scaled(field, w, u, coeff)
                if not is_zero_vec(ww):
                    new_remaining.append(ww)
            remaining = new_remaining
            blocks.append(("line", u, None))
            diagnostics["n_lines"] += 1
            diagnostics["steps"].append({"type": "line", "norm": int(line_norm), "remaining": len(remaining)})
            continue

        # If no anisotropic line can be manufactured, fall back to an explicit
        # hyperbolic pair.  This should only occur for genuinely degenerate or
        # convention-mismatched Hermitian data, but keeping the branch makes the
        # diagnostic clearer.
        pair: Optional[Tuple[int, int, int]] = None
        for i, u in enumerate(remaining):
            for j in range(i + 1, len(remaining)):
                hv = _h_pair(field, H, u, remaining[j])
                if hv != 0:
                    pair = (i, j, hv)
                    break
            if pair is not None:
                break
        if pair is None:
            raise RuntimeError("Hermitian top form appears degenerate: no anisotropic line or hyperbolic pair found")

        i, j, huv = pair
        u = remaining[i]
        v = remaining[j]
        inv_huv = field.inv(huv)
        v = [field.mul(x, inv_huv) for x in v]

        old_remaining = [w for k, w in enumerate(remaining) if k not in (i, j)]
        new_remaining = []
        for w in old_remaining:
            a = _h_pair(field, H, v, w)
            b = _h_pair(field, H, u, w)
            ww = _vec_add_scaled(field, w, u, a)
            ww = _vec_add_scaled(field, ww, v, b)
            if not is_zero_vec(ww):
                new_remaining.append(ww)
        remaining = new_remaining
        blocks.append(("pair", u, v))
        diagnostics["n_pairs"] += 1
        diagnostics["steps"].append({"type": "pair", "h_uv": int(huv), "remaining": len(remaining)})

    return blocks, diagnostics


# ---------------------------------------------------------------------------
# Block lifting and public sector extractor
# ---------------------------------------------------------------------------


def _cyclic_basis_or_error(F: np.ndarray, N: np.ndarray, v: np.ndarray, deg_q: int, L: int, p: int) -> np.ndarray:
    C = cyclic_submodule_basis(F, N, v.reshape(-1, 1), int(deg_q), int(L), p)
    C = independent_columns(mod_p(C, p), p)
    expected = int(deg_q) * int(L)
    if C.shape[1] != expected:
        raise RuntimeError(f"cyclic module has dimension {C.shape[1]}, expected {expected}")
    return C


def _lift_hermitian_blocks(
    F: np.ndarray,
    N: np.ndarray,
    Omega: np.ndarray,
    field: GF2Extension,
    top: TopQuotient,
    herm_blocks: Sequence[Tuple[str, List[int], Optional[List[int]]]],
    deg_q: int,
    p: int,
) -> Tuple[List[np.ndarray], List[Dict[str, Any]]]:
    spans: List[np.ndarray] = []
    records: List[Dict[str, Any]] = []
    for block_type, u_coords, v_coords in herm_blocks:
        u_top = _top_vector_from_E_coords(F, field, top.e_basis_reps, u_coords, p)
        Cu = _cyclic_basis_or_error(F, N, u_top, deg_q, top.L, p)

        if block_type == "line":
            if not is_nondegenerate(Omega, Cu, p):
                raise RuntimeError(
                    "Hermitian anisotropic line lifted to a degenerate cyclic module; "
                    "check trace/Hermitian convention"
                )
            spans.append(Cu)
            records.append(
                {
                    "type": "unitary_line",
                    "L": int(top.L),
                    "deg_q": int(deg_q),
                    "dim": int(Cu.shape[1]),
                    "half_dim": int(Cu.shape[1] // 2),
                }
            )
        elif block_type == "pair":
            if v_coords is None:
                raise RuntimeError("Hermitian pair missing second coordinate vector")
            v_top = _top_vector_from_E_coords(F, field, top.e_basis_reps, v_coords, p)
            Cv = _cyclic_basis_or_error(F, N, v_top, deg_q, top.L, p)
            Cuv = independent_columns(np.concatenate([Cu, Cv], axis=1), p)
            expected = 2 * int(deg_q) * int(top.L)
            if Cuv.shape[1] != expected:
                raise RuntimeError(f"Hermitian pair lift dimension {Cuv.shape[1]}, expected {expected}")
            if not is_nondegenerate(Omega, Cuv, p):
                raise RuntimeError(
                    "Hermitian hyperbolic pair lifted to a degenerate module; check trace/Hermitian convention"
                )
            spans.append(Cuv)
            records.append(
                {
                    "type": "unitary_pair",
                    "L": int(top.L),
                    "deg_q": int(deg_q),
                    "dim": int(Cuv.shape[1]),
                    "half_dim": int(Cuv.shape[1] // 2),
                }
            )
        else:
            raise ValueError(f"unknown Hermitian block type {block_type!r}")
    return spans, records


def atomic_blocks_in_self_sector_p2_nonunipotent_unitary(
    F_sec: np.ndarray,
    T_sec: np.ndarray,
    Omega: np.ndarray,
    N: np.ndarray,
    deg_q: int,
    max_exp: int,
    p: int,
    sector_key: Tuple[int, ...],
    poly_key: Tuple[int, ...],
    allow_fallback: bool = False,
    *,
    max_top_dim: int = 4096,  # retained for API compatibility; no exponential search is used.
) -> tuple[List[AtomicBlock], AtomicInvariant]:
    """
    Certified extractor for p=2, q=q*, deg(q)>1, q != x+1 self sectors.

    The sector is treated as a Hermitian/unitary primary sector over

        E = GF(2)[x]/(q),   involution alpha -> alpha^{-1}.

    For every nilpotent chain length L, the exact top quotient

        T_L = ker(N^L)/(ker(N^{L-1}) + N ker(N^{L+1}))

    is converted into an E-vector space.  The induced Hermitian form is recovered
    from base-field symplectic pairings via the trace pairing on E.  A deterministic
    Hermitian Gram-Schmidt decomposition gives E-lines and hyperbolic E-pairs;
    these are lifted to GF(2)[F]-cyclic modules and certified by direct base-field
    rank/nondegeneracy checks.

    No exponential enumeration of top vectors is used.
    """
    p = int(p)
    if p != 2:
        raise ValueError("atomic_blocks_in_self_sector_p2_nonunipotent_unitary is p=2 only")
    if int(deg_q) <= 1:
        raise ValueError("expected non-linear self-reciprocal primary polynomial")
    if int(deg_q) % 2 != 0:
        raise ValueError("self-reciprocal non-linear irreducibles over GF(2) should have even degree")

    F_sec = mod_p(np.asarray(F_sec, dtype=np.int64), p)
    N = mod_p(np.asarray(N, dtype=np.int64), p)
    Omega = mod_p(np.asarray(Omega, dtype=np.int64), p)
    T_sec = mod_p(np.asarray(T_sec, dtype=np.int64), p)
    q = tuple(int(c) & 1 for c in poly_key)
    field = GF2Extension(q)
    if field.d != int(deg_q):
        raise ValueError(f"poly degree {field.d} does not match deg_q={deg_q}")

    m0 = int(F_sec.shape[0])
    if F_sec.shape != (m0, m0) or N.shape != (m0, m0) or Omega.shape != (m0, m0):
        raise ValueError("shape mismatch in p=2 non-unipotent self-sector inputs")
    if m0 % 2 != 0:
        raise ValueError("self-sector dimension must be even")

    inv_data: Dict[str, Any] = {
        "status": "PENDING",
        "deg": int(deg_q),
        "exponent": int(max_exp),
        "p2_self_reciprocal_nonunipotent": True,
        "algorithm": "extension_field_hermitian_top_decomposition",
        "field_modulus": q,
        "field_degree": int(field.d),
        "involution": "a -> a^(2^(deg/2)) = a^{-1}",
        "note": "",
        "checks_passed": [],
        "top_forms": [],
        "blocks": [],
        "max_top_dim_parameter_ignored": int(max_top_dim),
    }

    # Normalize sector form to the standard symplectic matrix if necessary.
    Omega_std = omega_matrix(m0 // 2, p)
    T_map = T_sec
    if not np.array_equal(mod_p(Omega - Omega_std, p), np.zeros_like(Omega)):
        S = darboux_basis_from_span(Omega, np.eye(m0, dtype=np.int64), p)
        Sinv = inv_mod_mat(S, p)
        F = mod_p(Sinv @ F_sec @ S, p)
        Nstd = mod_p(Sinv @ N @ S, p)
        T_map = mod_p(T_sec @ S, p)
        inv_data["checks_passed"].append("darboux_normalized_sector_form")
    else:
        F = F_sec
        Nstd = N

    Omega_full = omega_matrix(m0 // 2, p)
    R = np.eye(m0, dtype=np.int64)
    blocks: List[AtomicBlock] = []
    block_records: List[Dict[str, Any]] = []

    while R.shape[1] > 0:
        dim_r = int(R.shape[1])
        if dim_r % 2 != 0:
            raise RuntimeError(f"remaining p=2 unitary self-sector dimension is odd ({dim_r})")

        # Keep R Darboux in the standard sector coordinates.
        if R.shape[1] and not np.array_equal(
            mod_p(R.T @ Omega_full @ R - omega_matrix(dim_r // 2, p), p),
            np.zeros((dim_r, dim_r), dtype=np.int64),
        ):
            R = darboux_basis_from_span(Omega_full, R, p)

        F_r = restrict_operator(F, R, p)
        N_r = restrict_operator(Nstd, R, p)
        Omega_r = omega_matrix(dim_r // 2, p)

        extracted = False
        for L in range(int(max_exp), 0, -1):
            top = _top_quotient_E_basis(F_r, N_r, int(L), int(deg_q), p)
            if top.dim_E == 0:
                continue

            H = _hermitian_top_matrix(F_r, N_r, Omega_r, field, top.e_basis_reps, int(L), p)
            herm_blocks, herm_diag = _hermitian_decompose(field, H)
            spans_r, records_r = _lift_hermitian_blocks(F_r, N_r, Omega_r, field, top, herm_blocks, int(deg_q), p)

            inv_data["top_forms"].append(
                {
                    "remaining_dim": int(dim_r),
                    "L": int(L),
                    "dim_f2": int(top.dim_f2),
                    "dim_E": int(top.dim_E),
                    "H": [[int(x) for x in row] for row in H.tolist()],
                    "hermitian_decomposition": herm_diag,
                    "n_lifted_blocks": int(len(spans_r)),
                }
            )

            if not spans_r:
                continue

            # Extract one block, then recompute the top quotient on the new
            # symplectic complement.  This avoids mixing coordinates from the old
            # remaining space with coordinates from the updated complement.
            span_r, rec = spans_r[0], records_r[0]
            if span_r.shape[1] % 2 != 0 or not is_nondegenerate(Omega_r, span_r, p):
                raise RuntimeError(f"lifted Hermitian block failed nondegeneracy check: {rec}")
            T_blk_r = darboux_basis_from_span(Omega_r, span_r, p)
            T_blk_std = mod_p(R @ T_blk_r, p)
            T_blk_amb = mod_p(T_map @ T_blk_std, p)
            half_dim = int(T_blk_amb.shape[1] // 2)
            rec = dict(rec)
            rec["half_dim"] = int(half_dim)
            block_records.append(rec)
            blocks.append(
                AtomicBlock(
                    T_blk=T_blk_amb,
                    half_dim=half_dim,
                    sector_key=sector_key,
                    inv=None,
                )
            )

            R_new = symplectic_orthogonal_complement_in_span(Omega_full, T_blk_std, R, p)
            R_new = independent_columns(mod_p(R_new, p), p)
            if R_new.shape[1] == 0:
                R = R_new
            else:
                if R_new.shape[1] % 2 != 0:
                    raise RuntimeError(f"remaining complement dimension is odd ({int(R_new.shape[1])})")
                R = darboux_basis_from_span(Omega_full, R_new, p)

            extracted = True
            break

        if not extracted:
            if allow_fallback:
                inv_data["status"] = "DEGRADED"
                inv_data["note"] = "p=2 Hermitian self-sector extraction stuck"
                break
            raise RuntimeError(
                "p=2 Hermitian self-sector extraction stuck; "
                f"remaining_dim={dim_r}, top_forms={inv_data['top_forms'][-3:]}"
            )

    if inv_data["status"] == "PENDING":
        inv_data["status"] = "OK"
        inv_data["checks_passed"].append("p2_self_reciprocal_nonunipotent_hermitian_blocks_extracted")

    half_dims = [int(b.half_dim) for b in blocks]
    inv_data["blocks"] = block_records
    inv_data["classification_complete"] = bool(inv_data["status"] == "OK")
    inv_data["implemented_block_families"] = [
        "Hermitian anisotropic E-lines lifted to self-dual cyclic GF(2)[F]-modules",
        "Hermitian hyperbolic E-pairs lifted to paired cyclic GF(2)[F]-modules",
    ]
    inv_data["cost_certificate"] = {
        "attained": bool(inv_data["status"] == "OK"),
        "complete": bool(inv_data["status"] == "OK"),
        "certified_minimal": bool(inv_data["status"] == "OK"),
        "lower_bound": max(half_dims, default=0),
        "qudit_cost": max(half_dims, default=0),
        "block_half_dims": half_dims,
        "reason": (
            "p=2 self-reciprocal non-unipotent sector certified by extension-field Hermitian "
            "top-space decomposition and direct base-field symplectic verification."
        ),
    }

    inv = AtomicInvariant(sector_key=sector_key, sector_type="self", poly_key=poly_key, data=inv_data)
    return blocks, inv
