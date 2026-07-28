from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any, List, Sequence, Tuple

import numpy as np

from sympleq.core.symmetries.modular_helpers import mod_p, inv_mod_mat, rank_mod, independent_columns, mat_pow_mod
from .module_invariants import cyclic_submodule_basis


Element = Tuple[int, ...]


class GFpExtension:
    """
    Small arithmetic layer for E = GF(p)[x]/(q), with q monic irreducible.

    Elements are tuples of length d containing power-basis coefficients
    a_0 + a_1 alpha + ... + a_{d-1} alpha^{d-1}.  The class is deliberately
    minimal and is meant for quotient-top linear algebra, not as a general CAS.
    """

    def __init__(self, modulus_coeffs: Sequence[int], p: int):
        self.p = int(p)
        q = [int(c) % self.p for c in modulus_coeffs]
        while len(q) > 1 and q[-1] == 0:
            q.pop()
        if len(q) < 2 or q[-1] % self.p != 1:
            raise ValueError(f"Expected a monic polynomial over GF({p}), got {modulus_coeffs!r}")
        self.modulus_coeffs = tuple(q)
        self.d = len(q) - 1
        self.zero: Element = tuple([0] * self.d)
        self.one: Element = tuple([1] + [0] * (self.d - 1))
        if self.d == 1:
            self.alpha: Element = (q[0] * 0 + 0,)  # not used for d=1
        else:
            self.alpha = tuple([0, 1] + [0] * (self.d - 2))
        self._alpha_inv: Element | None = None
        self._trace_matrix_inv: np.ndarray | None = None

    def _coerce(self, a: Sequence[int] | np.ndarray | int) -> Element:
        if isinstance(a, (int, np.integer)):
            out = [0] * self.d
            out[0] = int(a) % self.p
            return tuple(out)
        aa = np.asarray(a, dtype=np.int64).reshape(-1)
        if aa.size != self.d:
            raise ValueError(f"expected element length {self.d}, got {aa.size}")
        return tuple(int(x) % self.p for x in aa)

    def add(self, a: Sequence[int] | int, b: Sequence[int] | int) -> Element:
        aa = self._coerce(a)
        bb = self._coerce(b)
        return tuple((aa[i] + bb[i]) % self.p for i in range(self.d))

    def neg(self, a: Sequence[int] | int) -> Element:
        aa = self._coerce(a)
        return tuple((-aa[i]) % self.p for i in range(self.d))

    def sub(self, a: Sequence[int] | int, b: Sequence[int] | int) -> Element:
        return self.add(a, self.neg(b))

    def reduce_coeffs(self, coeffs: Sequence[int]) -> Element:
        c = [int(x) % self.p for x in coeffs]
        if not c:
            return self.zero
        # q(x) = x^d + q_{d-1}x^{d-1}+...+q_0, so x^d = -sum q_i x^i.
        while len(c) > self.d:
            lead = c[-1] % self.p
            if lead:
                shift = len(c) - self.d - 1
                for i in range(self.d):
                    c[shift + i] = (c[shift + i] - lead * self.modulus_coeffs[i]) % self.p
            c.pop()
        c += [0] * (self.d - len(c))
        return tuple(x % self.p for x in c)

    def mul(self, a: Sequence[int] | int, b: Sequence[int] | int) -> Element:
        aa = self._coerce(a)
        bb = self._coerce(b)
        tmp = [0] * (2 * self.d - 1)
        for i, ai in enumerate(aa):
            if ai:
                for j, bj in enumerate(bb):
                    if bj:
                        tmp[i + j] = (tmp[i + j] + ai * bj) % self.p
        return self.reduce_coeffs(tmp)

    def pow(self, a: Sequence[int] | int, e: int) -> Element:
        ee = int(e)
        if ee < 0:
            return self.pow(self.inv(a), -ee)
        out = self.one
        base = self._coerce(a)
        while ee:
            if ee & 1:
                out = self.mul(out, base)
            base = self.mul(base, base)
            ee >>= 1
        return out

    def inv(self, a: Sequence[int] | int) -> Element:
        aa = self._coerce(a)
        if aa == self.zero:
            raise ZeroDivisionError("inverse of zero in finite-field extension")
        return self.pow(aa, self.p ** self.d - 2)

    def div(self, a: Sequence[int] | int, b: Sequence[int] | int) -> Element:
        return self.mul(a, self.inv(b))

    def is_zero(self, a: Sequence[int] | int) -> bool:
        return self._coerce(a) == self.zero

    def coeff_vector(self, a: Sequence[int] | int) -> np.ndarray:
        return np.asarray(self._coerce(a), dtype=np.int64).reshape(self.d, 1)

    def from_coeff_vector(self, c: np.ndarray) -> Element:
        return self._coerce(np.asarray(c, dtype=np.int64).reshape(-1))

    def elements(self) -> List[Element]:
        return [tuple(int(x) for x in coeffs) for coeffs in product(range(self.p), repeat=self.d)]

    def nonzero_elements(self) -> List[Element]:
        return [a for a in self.elements() if a != self.zero]

    def trace_to_base(self, a: Sequence[int] | int) -> int:
        """Trace E -> GF(p)."""
        aa = self._coerce(a)
        acc = self.zero
        x = aa
        for _ in range(self.d):
            acc = self.add(acc, x)
            x = self.pow(x, self.p)
        # The trace lies in the prime field, represented by a constant element.
        # Round small numerical/representation noise by checking nonconstant entries.
        if any(acc[i] % self.p for i in range(1, self.d)):
            raise RuntimeError(f"trace did not land in the base field: {acc}")
        return int(acc[0]) % self.p

    def alpha_power(self, k: int) -> Element:
        if k == 0:
            return self.one
        if self.d == 1:
            return self.one
        return self.pow(self.alpha, int(k))

    def alpha_inv(self) -> Element:
        if self._alpha_inv is None:
            if self.d == 1:
                self._alpha_inv = self.one
            else:
                self._alpha_inv = self.inv(self.alpha)
        return self._alpha_inv

    def conj(self, a: Sequence[int] | int) -> Element:
        """Reciprocal involution alpha -> alpha^{-1}."""
        aa = self._coerce(a)
        if self.d == 1:
            return aa
        out = self.zero
        pow_ai = self.one
        alpha_inv = self.alpha_inv()
        for ci in aa:
            if ci:
                out = self.add(out, self.mul(ci, pow_ai))
            pow_ai = self.mul(pow_ai, alpha_inv)
        return out

    def trace_pairing_inverse(self) -> np.ndarray:
        """
        Inverse of M_{a,k}=Tr(alpha^a alpha^k), 0<=a,k<d.
        It maps trace coordinates to power-basis coefficients.
        """
        if self._trace_matrix_inv is not None:
            return self._trace_matrix_inv
        M = np.zeros((self.d, self.d), dtype=np.int64)
        powers = [self.alpha_power(k) for k in range(2 * self.d)]
        for a in range(self.d):
            for k in range(self.d):
                M[a, k] = self.trace_to_base(self.mul(powers[a], powers[k]))
        if rank_mod(M, self.p) != self.d:
            raise RuntimeError("trace pairing matrix is singular; modulus may not be irreducible")
        self._trace_matrix_inv = inv_mod_mat(M, self.p)
        return self._trace_matrix_inv

    def from_trace_values(self, values: Sequence[int]) -> Element:
        b = np.asarray(values, dtype=np.int64).reshape(self.d, 1) % self.p
        c = mod_p(self.trace_pairing_inverse() @ b, self.p)
        return self.from_coeff_vector(c)


class GF2Extension:
    """
    Packed GF(2)[x]/(q) arithmetic for high-degree binary extension fields.

    Generic :class:`GFpExtension` represents elements as coefficient tuples,
    which is useful for arbitrary p but too slow for the large p=2
    self-reciprocal sectors that appear in random symplectic decompositions.
    This class provides the arithmetic needed by the p=2 Hermitian sector while
    representing elements by Python integers whose binary expansion stores the
    power-basis coefficients.
    """

    def __init__(self, modulus_coeffs: Sequence[int]):
        q = [int(c) & 1 for c in modulus_coeffs]
        while len(q) > 1 and q[-1] == 0:
            q.pop()
        if len(q) < 3 or q[-1] != 1:
            raise ValueError(
                f"Expected a monic irreducible polynomial of degree >1 over GF(2), got {modulus_coeffs!r}"
            )
        self.modulus_coeffs = tuple(q)
        self.d = len(q) - 1
        self.mask = (1 << self.d) - 1

        self._q_low = 0
        for i, c in enumerate(q[:-1]):
            if c & 1:
                self._q_low |= 1 << i

        self.alpha = 2 if self.d > 1 else 1
        self._trace_matrix_inv: np.ndarray | None = None

    def add(self, a: int, b: int) -> int:
        return (int(a) ^ int(b)) & self.mask

    sub = add

    def reduce(self, x: int) -> int:
        x = int(x)
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
        return self.pow(a, (1 << self.d) - 2)

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.inv(b))

    def conj(self, a: int) -> int:
        """
        The reciprocal involution alpha -> alpha^{-1}; for self-reciprocal
        irreducibles over GF(2), this is Frobenius 2^(d/2).
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
        Inverse of M_{a,k} = Tr(alpha^a alpha^k), 0<=a,k<d.
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
        b = np.asarray(values, dtype=np.int64).reshape(self.d, 1) % 2
        c = mod_p(self.trace_pairing_inverse() @ b, 2)
        return self.from_coeff_vector(c)


@dataclass(slots=True)
class TopQuotientEBasis:
    L: int
    denom: np.ndarray
    top_reps_fp: np.ndarray
    e_basis_reps: List[np.ndarray]
    e_orbit_basis: np.ndarray
    dim_fp: int
    dim_E: int


def _orbit_columns(F: np.ndarray, v: np.ndarray, deg_q: int, p: int) -> np.ndarray:
    cols: List[np.ndarray] = []
    P = np.eye(F.shape[0], dtype=np.int64)
    vv = mod_p(v.reshape(-1, 1), p)
    for _ in range(int(deg_q)):
        cols.append(mod_p(P @ vv, p))
        P = mod_p(P @ F, p)
    return independent_columns(np.concatenate(cols, axis=1), p)


def top_quotient_E_basis(
    F: np.ndarray,
    N: np.ndarray,
    L: int,
    deg_q: int,
    p: int,
    *,
    top_reps: np.ndarray,
    denom: np.ndarray | None = None,
) -> TopQuotientEBasis:
    """
    Convert GF(p)-top representatives for K_L/(K_{L-1}+NK_{L+1}) into an
    E-basis, E=GF(p)[x]/(q), using only ordinary rank tests modulo the quotient
    denominator.  One returned representative corresponds to one q^L cyclic
    module, not to one base-field vector.
    """
    F = mod_p(F, p)
    N = mod_p(N, p)
    top_reps = independent_columns(mod_p(top_reps, p), p)
    d0 = F.shape[0]
    if denom is None:
        denom = np.zeros((d0, 0), dtype=np.int64)
    denom = independent_columns(mod_p(denom, p), p) if denom.size else np.zeros((d0, 0), dtype=np.int64)
    deg_q = int(deg_q)
    L = int(L)
    dim_fp = int(top_reps.shape[1])
    if dim_fp == 0:
        return TopQuotientEBasis(L, denom, top_reps, [], np.zeros((d0, 0), dtype=np.int64), 0, 0)
    if deg_q <= 0 or dim_fp % deg_q != 0:
        raise RuntimeError(f"top quotient dimension {dim_fp} is not divisible by deg_q={deg_q}")
    dim_E = dim_fp // deg_q

    selected: List[np.ndarray] = []
    selected_orbits = np.zeros((d0, 0), dtype=np.int64)
    base_rank = rank_mod(denom, p) if denom.shape[1] else 0

    for j in range(top_reps.shape[1]):
        v = top_reps[:, j:j + 1]
        orbit = _orbit_columns(F, v, deg_q, p)
        if orbit.shape[1] != deg_q:
            continue
        old = np.concatenate([denom, selected_orbits], axis=1) if selected_orbits.shape[1] else denom
        old_rank = base_rank + selected_orbits.shape[1]
        trial = np.concatenate([old, orbit], axis=1) if old.shape[1] else orbit
        new_rank = rank_mod(trial, p)
        if new_rank - old_rank != deg_q:
            continue
        # Guardrail: the selected top vector must generate a full q^L cyclic module.
        try:
            C = cyclic_submodule_basis(F, N, v, deg_q, L, p)
        except RuntimeError:
            continue
        if C.shape[1] != deg_q * L:
            continue
        selected.append(v)
        selected_orbits = np.concatenate([selected_orbits, orbit], axis=1) if selected_orbits.shape[1] else orbit
        if len(selected) == dim_E:
            break

    if len(selected) != dim_E:
        raise RuntimeError(f"could not construct E-basis of top quotient: needed {dim_E}, got {len(selected)}")

    return TopQuotientEBasis(
        L=L,
        denom=denom,
        top_reps_fp=top_reps,
        e_basis_reps=selected,
        e_orbit_basis=selected_orbits,
        dim_fp=dim_fp,
        dim_E=dim_E,
    )


def field_scalar_apply_to_vector(
    F: np.ndarray,
    field: GFpExtension,
    scalar: Sequence[int] | int,
    v: np.ndarray,
    p: int,
) -> np.ndarray:
    coeff = field.coeff_vector(scalar).reshape(-1)
    out = np.zeros_like(v.reshape(-1, 1), dtype=np.int64)
    P = np.eye(F.shape[0], dtype=np.int64)
    vv = v.reshape(-1, 1)
    for a in range(field.d):
        ca = int(coeff[a]) % p
        if ca:
            out = mod_p(out + ca * (P @ vv), p)
        P = mod_p(P @ F, p)
    return mod_p(out, p)


def top_vector_from_E_coords(
    F: np.ndarray,
    field: GFpExtension,
    e_reps: Sequence[np.ndarray],
    coords: Sequence[Sequence[int] | int],
    p: int,
    *,
    conjugate_coeffs: bool = False,
) -> np.ndarray:
    if len(e_reps) != len(coords):
        raise ValueError("coordinate length mismatch")
    if not e_reps:
        raise ValueError("no E representatives supplied")
    out = np.zeros_like(e_reps[0].reshape(-1, 1), dtype=np.int64)
    for v, c in zip(e_reps, coords):
        cc = field.conj(c) if conjugate_coeffs else field._coerce(c)
        if cc != field.zero:
            out = mod_p(out + field_scalar_apply_to_vector(F, field, cc, v, p), p)
    return mod_p(out, p)


def symp_scalar(Omega: np.ndarray, u: np.ndarray, v: np.ndarray, p: int) -> int:
    return int(mod_p(u.reshape(1, -1) @ Omega @ v.reshape(-1, 1), p).reshape(-1)[0]) % p


def hermitian_top_matrix(
    F: np.ndarray,
    N: np.ndarray,
    Omega: np.ndarray,
    field: GFpExtension,
    e_reps: Sequence[np.ndarray],
    L: int,
    p: int,
) -> np.ndarray:
    """
    Recover the E-valued top form H_ij from base-field pairings
    b_a = <F^a e_i, N^{L-1} e_j> by trace-duality in the power basis.
    """
    m = len(e_reps)
    H = np.empty((m, m), dtype=object)
    Npow = mat_pow_mod(N, int(L) - 1, p)
    F_pows = [np.eye(F.shape[0], dtype=np.int64)]
    for _ in range(1, field.d):
        F_pows.append(mod_p(F_pows[-1] @ F, p))
    for i, ui in enumerate(e_reps):
        for j, vj in enumerate(e_reps):
            Nv = mod_p(Npow @ vj, p)
            values = [symp_scalar(Omega, mod_p(P @ ui, p), Nv, p) for P in F_pows]
            H[i, j] = field.from_trace_values(values)
    return H


def h_pair(
    field: GFpExtension,
    H: np.ndarray,
    x: Sequence[Sequence[int] | int],
    y: Sequence[Sequence[int] | int],
) -> Element:
    m = H.shape[0]
    acc = field.zero
    for i in range(m):
        xi = field._coerce(x[i])
        if xi == field.zero:
            continue
        cxi = field.conj(xi)
        for j in range(m):
            yj = field._coerce(y[j])
            hij = H[i, j]
            if yj == field.zero or field.is_zero(hij):
                continue
            acc = field.add(acc, field.mul(field.mul(cxi, hij), yj))
    return acc


def _vec_add_scaled(
    field: GFpExtension,
    x: List[Element],
    y: Sequence[Element],
    a: Sequence[int] | int,
) -> List[Element]:
    aa = field._coerce(a)
    return [field.add(xi, field.mul(yi, aa)) for xi, yi in zip(x, y)]


def _vec_scale(field: GFpExtension, x: Sequence[Element], a: Sequence[int] | int) -> List[Element]:
    aa = field._coerce(a)
    return [field.mul(xi, aa) for xi in x]


def canonical_E_basis(field: GFpExtension, m: int) -> List[List[Element]]:
    out: List[List[Element]] = []
    for i in range(m):
        v = [field.zero] * m
        v[i] = field.one
        out.append(v)
    return out


def hermitian_orthogonal_lines(field: GFpExtension, H: np.ndarray) -> Tuple[List[List[Element]], dict[str, Any]]:
    """
    Deterministic Hermitian Gram-Schmidt.  Returns E-coordinate lines with
    nonzero norm.  If the current basis is isotropic, it scans one field
    coefficient to manufacture an anisotropic vector from a coupled pair.
    """
    m = int(H.shape[0])
    remaining = canonical_E_basis(field, m)
    lines: List[List[Element]] = []
    info: dict[str, Any] = {"dim_E": m, "n_lines": 0, "manufactured": 0, "steps": []}

    def is_zero_vec(v: Sequence[Element]) -> bool:
        return all(field.is_zero(x) for x in v)

    while remaining:
        line_vec = None
        line_norm = field.zero
        line_remaining: List[List[Element]] | None = None

        for idx, u in enumerate(remaining):
            nrm = h_pair(field, H, u, u)
            if not field.is_zero(nrm):
                line_vec = u
                line_norm = nrm
                line_remaining = [w for k, w in enumerate(remaining) if k != idx]
                break

        if line_vec is None:
            pair = None
            for i, u in enumerate(remaining):
                for j in range(i + 1, len(remaining)):
                    hv = h_pair(field, H, u, remaining[j])
                    if not field.is_zero(hv):
                        pair = (i, j)
                        break
                if pair is not None:
                    break
            if pair is None:
                raise RuntimeError("Hermitian top form appears degenerate")
            i, j = pair
            u = remaining[i]
            v = remaining[j]
            for c in field.nonzero_elements():
                cand = _vec_add_scaled(field, u, v, c)
                nrm = h_pair(field, H, cand, cand)
                if not field.is_zero(nrm):
                    line_vec = cand
                    line_norm = nrm
                    line_remaining = [remaining[k] for k in range(len(remaining)) if k != i]
                    info["manufactured"] += 1
                    break
            if line_vec is None:
                raise RuntimeError("failed to manufacture anisotropic Hermitian line from coupled pair")

        assert line_vec is not None and line_remaining is not None
        inv_norm = field.inv(line_norm)
        new_remaining: List[List[Element]] = []
        for w in line_remaining:
            coeff = field.mul(inv_norm, h_pair(field, H, line_vec, w))
            ww = _vec_add_scaled(field, w, line_vec, field.neg(coeff))
            if not is_zero_vec(ww):
                new_remaining.append(ww)
        lines.append(line_vec)
        remaining = new_remaining
        info["n_lines"] += 1
        info["steps"].append({"norm": tuple(line_norm), "remaining": len(remaining)})

    return lines, info
