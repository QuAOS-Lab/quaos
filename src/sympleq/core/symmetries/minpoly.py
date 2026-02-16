from __future__ import annotations
import numpy as np
from typing import List, Tuple, Optional

from sympleq.core.symmetries.modular_helpers import mod_p, matmul_mod, inv_mod_scalar
from sympleq.core.symmetries.polynomials_fp import (poly_monic, poly_lcm, poly_trim, poly_is_zero, poly_divmod,
                                                    poly_gcd, poly_sub, poly_mul, poly_add)


class _SpanBasisSimple:
    """
    Span basis over GF(p) for vectors in GF(p)^n, no coefficient tracking.
    Supports membership test without mutating the basis.
    """
    def __init__(self, n: int, p: int):
        self.n = int(n)
        self.p = int(p)
        self.pivots: List[int] = []
        self.vecs: List[np.ndarray] = []  # each 1D length n, pivot normalized to 1 (if p odd)

    def _reduce_inplace(self, v: np.ndarray) -> None:
        """Reduce v by current basis in-place."""
        p = self.p
        if p == 2:
            for i, piv in enumerate(self.pivots):
                if v[piv] & 1:
                    v ^= self.vecs[i]
        else:
            for i, piv in enumerate(self.pivots):
                a = int(v[piv] % p)
                if a:
                    v[:] = mod_p(v - a * self.vecs[i], p)

    def in_span(self, v: np.ndarray) -> bool:
        """Return True iff v lies in the current span (does not mutate the basis)."""
        v = mod_p(np.asarray(v, dtype=np.int64).reshape(-1), self.p).copy()
        if v.size != self.n:
            raise ValueError("wrong vector length")
        self._reduce_inplace(v)
        return not np.any(v)  # v == 0

    def add(self, v: np.ndarray) -> bool:
        """Add v if independent; return True iff span increased."""
        p = self.p
        v = mod_p(np.asarray(v, dtype=np.int64).reshape(-1), p).copy()
        if v.size != self.n:
            raise ValueError("wrong vector length")

        self._reduce_inplace(v)
        nz = np.flatnonzero(v)
        if nz.size == 0:
            return False
        piv = int(nz[0])

        # Normalize pivot (odd p)
        if p != 2:
            inv = inv_mod_scalar(int(v[piv]), p)
            v = mod_p(v * inv, p)

        # Eliminate pivot from existing basis
        if p == 2:
            for i in range(len(self.vecs)):
                if self.vecs[i][piv] & 1:
                    self.vecs[i] ^= v
        else:
            for i in range(len(self.vecs)):
                a = int(self.vecs[i][piv] % p)
                if a:
                    self.vecs[i] = mod_p(self.vecs[i] - a * v, p)

        # Insert by pivot order
        ins = 0
        while ins < len(self.pivots) and self.pivots[ins] < piv:
            ins += 1
        self.pivots.insert(ins, piv)
        self.vecs.insert(ins, v)
        return True


class _SpanBasisWithCombo:
    """
    Span basis over GF(p) with *fixed-length* combo vectors.

    Invariant: each stored basis vector vecs[i] equals a linear combination of
    generator columns g_0,...,g_{gen_dim-1} given by combos[i, :gen_dim].

    Crucial: combos are stored at fixed length (max_gen), so we never resize/copy old ones.
    """
    def __init__(self, n: int, p: int, max_gen: int):
        self.n = int(n)
        self.p = int(p)
        self.max_gen = int(max_gen)
        self.gen_dim = 0

        self.pivots: List[int] = []
        self.vecs: List[np.ndarray] = []       # 1D length n
        self.combos: List[np.ndarray] = []     # 1D length max_gen (fixed)

    def set_gen_dim(self, t: int) -> None:
        t = int(t)
        if t < self.gen_dim:
            raise ValueError("gen_dim cannot shrink")
        if t > self.max_gen:
            raise ValueError("gen_dim exceeds max_gen")
        self.gen_dim = t

    def _first_nz(self, v: np.ndarray) -> int | None:
        nz = np.flatnonzero(v)
        return int(nz[0]) if nz.size else None

    def add_or_relation(self, v: np.ndarray, combo_seed_index: int) -> Tuple[bool, np.ndarray]:
        """
        Attempt to add vector v as a new independent direction.
        Start with combo = e_{combo_seed_index} (in generator coords).

        Returns:
          (added, relation_combo_prefix)

        If added=True: relation_combo_prefix is undefined (returned as zeros).
        If added=False: v reduced to zero, and we return combo[:gen_dim], which encodes
                        a relation sum_j combo[j] * g_j = 0 with combo[combo_seed_index] = 1.
        """
        p = self.p
        v = mod_p(np.asarray(v, dtype=np.int64).reshape(-1), p).copy()
        if v.size != self.n:
            raise ValueError("wrong vector length")
        if not (0 <= combo_seed_index < self.gen_dim):
            raise ValueError("combo_seed_index out of range for current gen_dim")

        combo = np.zeros(self.max_gen, dtype=np.int64)
        combo[combo_seed_index] = 1

        # Reduce v, and apply identical row-ops to combo
        if p == 2:
            for i, piv in enumerate(self.pivots):
                if v[piv] & 1:
                    v ^= self.vecs[i]
                    combo[: self.gen_dim] ^= self.combos[i][: self.gen_dim]
        else:
            for i, piv in enumerate(self.pivots):
                a = int(v[piv] % p)
                if a:
                    v = mod_p(v - a * self.vecs[i], p)
                    combo[: self.gen_dim] = mod_p(combo[: self.gen_dim] - a * self.combos[i][: self.gen_dim], p)

        piv = self._first_nz(v)
        if piv is None:
            # dependent: 0 = sum_j combo[j] g_j, with combo[seed]=1 untouched
            return False, mod_p(combo[: self.gen_dim].copy(), p)

        # independent: normalize pivot
        if p != 2:
            inv = inv_mod_scalar(int(v[piv]), p)
            v = mod_p(v * inv, p)
            combo[: self.gen_dim] = mod_p(combo[: self.gen_dim] * inv, p)

        # eliminate new pivot from existing basis
        if p == 2:
            for i in range(len(self.vecs)):
                if self.vecs[i][piv] & 1:
                    self.vecs[i] ^= v
                    self.combos[i][: self.gen_dim] ^= combo[: self.gen_dim]
        else:
            for i in range(len(self.vecs)):
                a = int(self.vecs[i][piv] % p)
                if a:
                    self.vecs[i] = mod_p(self.vecs[i] - a * v, p)
                    self.combos[i][: self.gen_dim] = mod_p(self.combos[i][: self.gen_dim] - a * combo[: self.gen_dim], p)

        # insert by pivot order
        ins = 0
        while ins < len(self.pivots) and self.pivots[ins] < piv:
            ins += 1
        self.pivots.insert(ins, piv)
        self.vecs.insert(ins, v)

        combo_store = np.zeros(self.max_gen, dtype=np.int64)
        combo_store[: self.gen_dim] = mod_p(combo[: self.gen_dim], p)
        self.combos.insert(ins, combo_store)

        return True, np.zeros(self.gen_dim, dtype=np.int64)


def minimal_poly_for_vector(F: np.ndarray, v: np.ndarray, p: int) -> np.ndarray:
    """
    Compute minimal polynomial m_v such that m_v(F) v = 0.
    Deterministic, and fast in p=2 because we never resize combos.
    """
    F = mod_p(F, p)
    v = mod_p(np.asarray(v, dtype=np.int64).reshape(-1), p)
    n2 = F.shape[0]
    if v.size != n2:
        raise ValueError("v has wrong length")

    B = _SpanBasisWithCombo(n2, p, max_gen=n2 + 1)

    w = v.copy()
    for t in range(0, n2 + 1):
        B.set_gen_dim(t + 1)

        added, rel = B.add_or_relation(w, combo_seed_index=t)
        if not added:
            # relation: sum_{j=0}^t rel[j] F^j v = 0 and rel[t]=1
            return poly_monic(rel, p)

        # advance Krylov
        w = matmul_mod(F, w.reshape(-1, 1), p).reshape(-1)

    raise RuntimeError("minimal_poly_for_vector: exceeded ambient dimension (unexpected)")


def minimal_polynomial(F: np.ndarray, p: int) -> np.ndarray:
    """
    Deterministic minimal polynomial of the matrix F over GF(p).

    Provably correct: m_F is the lcm of minimal polynomials m_v over any
    generating set; the standard basis generates GF(p)^{n2}.
    """
    F = mod_p(F, p)
    n2 = F.shape[0]
    m = np.array([1], dtype=np.int64)

    for i in range(n2):
        ei = np.zeros(n2, dtype=np.int64)
        ei[i] = 1
        mi = minimal_poly_for_vector(F, ei, p)
        m = poly_lcm(m, mi, p)

        # optional early stop: degree can't exceed n2
        if len(m) - 1 == n2:
            break

    return poly_monic(m, p)


def _poly_deg(f: np.ndarray) -> int:
    f = poly_trim(f)
    if poly_is_zero(f):
        return -1
    return len(f) - 1


def _poly_is_one(f: np.ndarray, p: int) -> bool:
    f = poly_monic(f, p)
    return (len(f) == 1 and int(f[0]) % p == 1)


def _poly_mod(a: np.ndarray, m: np.ndarray, p: int) -> np.ndarray:
    a = mod_p(poly_trim(a), p)
    m = poly_monic(m, p)  # modulus should be monic
    _, r = poly_divmod(a, m, p)
    return mod_p(poly_trim(r), p)


def _poly_pow_mod(a: np.ndarray, e: int, m: np.ndarray, p: int) -> np.ndarray:
    """
    Compute a(x)^e mod m(x) over GF(p) in the quotient ring.
    Does not normalize to monic here.
    """
    a = _poly_mod(a, m, p)
    res = np.array([1], dtype=np.int64)
    base = a.copy()
    ee = int(e)
    while ee > 0:
        if ee & 1:
            res = _poly_mod(poly_mul(res, base, p), m, p)
        base = _poly_mod(poly_mul(base, base, p), m, p)
        ee >>= 1
    return mod_p(poly_trim(res), p)


def _derivative_poly(f: np.ndarray, p: int) -> np.ndarray:
    f = mod_p(poly_trim(f), p)
    if len(f) <= 1:
        return np.array([0], dtype=np.int64)
    df = np.array([(i * int(f[i])) % p for i in range(1, len(f))], dtype=np.int64)
    return poly_trim(df)


def _pth_root(f: np.ndarray, p: int) -> np.ndarray:
    """
    If f(x) = g(x)^p (i.e. only coefficients at degrees multiple of p),
    return g. Used when derivative is identically zero.
    """
    f = poly_trim(f)
    out = np.zeros((len(f) + p - 1) // p, dtype=np.int64)
    for i in range(0, len(f), p):
        out[i // p] = int(f[i])
    return poly_trim(out)


def _squarefree_decomposition(f: np.ndarray, p: int) -> List[Tuple[np.ndarray, int]]:
    """
    Return [(f1,1),(f2,2),...] such that f = Π fi^i, fi squarefree and pairwise coprime.
    """
    f = poly_monic(f, p)
    if _poly_deg(f) <= 0:
        return [(f, 1)]

    df = _derivative_poly(f, p)
    if poly_is_zero(df):
        # f is a p-th power: f(x) = g(x^p)
        g = _pth_root(f, p)
        sub = _squarefree_decomposition(g, p)
        return [(h, e * p) for (h, e) in sub]

    g = poly_gcd(f, df, p)
    w, r = poly_divmod(f, g, p)
    if not poly_is_zero(r):
        raise RuntimeError("squarefree_decomposition: division remainder nonzero (unexpected)")

    out: List[Tuple[np.ndarray, int]] = []
    i = 1
    while not _poly_is_one(w, p):
        y = poly_gcd(w, g, p)
        fi, r = poly_divmod(w, y, p)
        if not poly_is_zero(r):
            raise RuntimeError("squarefree_decomposition: division remainder nonzero (unexpected)")
        if not _poly_is_one(fi, p):
            out.append((poly_monic(fi, p), i))
        w = y
        g, r = poly_divmod(g, y, p)
        if not poly_is_zero(r):
            # g should divide y
            g = poly_monic(g, p)
        i += 1
        if _poly_deg(w) <= 0:
            break

    if not _poly_is_one(g, p) and not poly_is_zero(g):
        # remaining part is a p-th power
        g_root = _pth_root(g, p)
        sub = _squarefree_decomposition(g_root, p)
        out.extend([(h, e * p) for (h, e) in sub])

    if not out:
        out = [(f, 1)]
    return out


def _distinct_degree_factorization(f: np.ndarray, p: int) -> List[Tuple[np.ndarray, int]]:
    """
    Factor squarefree monic f into products of irreducibles of each degree.
    Returns [(g1,1),(g2,2),...] where gi is product of all irreducibles of degree i.
    """
    f = poly_monic(f, p)
    n = _poly_deg(f)
    if n <= 1:
        return [(f, n)]

    x = np.array([0, 1], dtype=np.int64)
    h = x.copy()
    res: List[Tuple[np.ndarray, int]] = []

    # h <- x^{p^d} mod f iteratively using Frobenius: h <- h^p mod f
    for d in range(1, n + 1):
        h = _poly_pow_mod(h, p, f, p)  # x^{p^d} mod f
        g = poly_gcd(poly_sub(h, x, p), f, p)
        if not _poly_is_one(g, p) and _poly_deg(g) >= 1 and _poly_deg(g) < _poly_deg(f):
            res.append((poly_monic(g, p), d))
            q, r = poly_divmod(f, g, p)
            if not poly_is_zero(r):
                raise RuntimeError("DDF: division remainder nonzero (unexpected)")
            f = poly_monic(q, p)
            h = _poly_mod(h, f, p)
        if _poly_deg(f) <= 1:
            break

    if _poly_deg(f) >= 1 and not _poly_is_one(f, p):
        res.append((poly_monic(f, p), _poly_deg(f)))
    return res


def _random_poly(deg_bound: int, p: int, rng: np.random.Generator) -> np.ndarray:
    """
    Random polynomial of degree < deg_bound (coeffs in 0..p-1), not identically zero.
    """
    if deg_bound <= 1:
        return np.array([rng.integers(0, p)], dtype=np.int64)
    a = rng.integers(0, p, size=deg_bound, dtype=np.int64)
    if np.all(a % p == 0):
        a[0] = 1
    return poly_trim(a)


def _equal_degree_factorization(f: np.ndarray, d: int, p: int, rng: np.random.Generator) -> List[np.ndarray]:
    """
    Split squarefree f where all irreducible factors have degree d.
    Returns monic irreducible factors.

    Uses:
      - odd p: standard Cantor–Zassenhaus (quadratic character) split
      - p=2  : Cantor–Zassenhaus trace split (required in characteristic 2)
    """
    f = poly_monic(f, p)
    n = _poly_deg(f)
    if n == d:
        return [f]

    if p == 2:
        return _equal_degree_factorization_char2(f, d, rng)

    # ---- odd characteristic splitter (your existing logic, with a try cap) ----
    exp = (p**d - 1) // 2
    max_tries = 2000
    for _ in range(max_tries):
        a = _random_poly(n, p, rng)
        a = _poly_mod(a, f, p)
        if poly_is_zero(a):
            continue

        a_pow = _poly_pow_mod(a, exp, f, p)

        g = poly_gcd(poly_sub(a_pow, np.array([1], dtype=np.int64), p), f, p)
        if _poly_is_one(g, p) or np.array_equal(g, f):
            # optional: also try gcd(a_pow + 1, f)
            g2 = poly_gcd(poly_add(a_pow, np.array([1], dtype=np.int64), p), f, p)
            if _poly_is_one(g2, p) or np.array_equal(g2, f):
                continue
            g = g2

        q, r = poly_divmod(f, g, p)
        if not poly_is_zero(r):
            continue

        return _equal_degree_factorization(poly_monic(g, p), d, p, rng) + \
               _equal_degree_factorization(poly_monic(q, p), d, p, rng)

    raise RuntimeError("equal_degree_factorization: exceeded max_tries (odd characteristic splitter)")


def _equal_degree_factorization_char2(f: np.ndarray, d: int, rng: np.random.Generator) -> List[np.ndarray]:
    """
    Cantor–Zassenhaus equal-degree factorization for p=2.

    Split using the trace map:
      Tr_{GF(2^d)/GF(2)}(a) = a + a^{2} + a^{2^2} + ... + a^{2^{d-1}}
    Compute h = Tr(a) mod f, then g = gcd(h, f) gives a nontrivial split with good probability.
    """
    p = 2
    f = poly_monic(f, p)
    n = _poly_deg(f)
    if n == d:
        return [f]

    max_tries = 4000
    zero = np.array([0], dtype=np.int64)

    for _ in range(max_tries):
        a = _random_poly(n, p, rng)
        a = _poly_mod(a, f, p)
        if poly_is_zero(a):
            continue

        # h = trace(a) = a + a^2 + ... + a^{2^{d-1}}  (mod f)
        h = zero.copy()
        t = a.copy()
        for __ in range(d):
            h = _poly_mod(poly_add(h, t, p), f, p)
            # t <- t^2 mod f
            t = _poly_mod(poly_mul(t, t, p), f, p)

        g = poly_gcd(h, f, p)
        if _poly_is_one(g, p) or np.array_equal(g, f):
            continue

        q, r = poly_divmod(f, g, p)
        if not poly_is_zero(r):
            continue

        return _equal_degree_factorization_char2(poly_monic(g, p), d, rng) + \
            _equal_degree_factorization_char2(poly_monic(q, p), d, rng)

    raise RuntimeError("equal_degree_factorization_char2: exceeded max_tries (p=2 trace splitter)")


def factor_poly_over_fp(f: np.ndarray, p: int, rng: Optional[np.random.Generator] = None) -> List[np.ndarray]:
    """
    Factor monic polynomial f over GF(p) into monic irreducibles.
    Returns a list of irreducible factors, repeated by multiplicity.
    """
    if rng is None:
        rng = np.random.default_rng()

    f0 = poly_monic(f, p)          # keep original for sanity check
    f = f0.copy()

    lin_factors: List[np.ndarray] = []

    # Optional linear-root peeling in odd characteristic
    if p > 2 and _poly_deg(f) >= 1:
        for r in range(p):
            while _poly_deg(f) >= 1 and _poly_eval_at(f, r, p) == 0:
                lin = np.array([(-r) % p, 1], dtype=np.int64)  # (x - r)
                q, rem = poly_divmod(f, lin, p)
                if not poly_is_zero(rem):
                    break
                lin_factors.append(poly_monic(lin, p))
                f = poly_monic(q, p)

    n = _poly_deg(f)
    if n <= 0:
        # f is constant (for monic inputs, typically 1). Return the peeled linear factors.
        out = lin_factors if lin_factors else [f]
        return out

    if n == 1:
        out = lin_factors + [poly_monic(f, p)]
        # sanity
        prod = np.array([1], dtype=np.int64)
        for h in out:
            prod = poly_mul(prod, h, p)
        prod = poly_monic(prod, p)
        if not np.array_equal(prod, f0):
            raise RuntimeError("factor_poly_over_fp sanity check failed: product(factors) != f0")
        return out

    # Squarefree decomposition
    sq = _squarefree_decomposition(f, p)

    out: List[np.ndarray] = []
    for sf, e in sq:
        sf = poly_monic(sf, p)
        if _poly_deg(sf) <= 0:
            continue

        # Distinct degree factorization
        ddf = _distinct_degree_factorization(sf, p)
        for g, d in ddf:
            g = poly_monic(g, p)
            if _poly_deg(g) <= 0:
                continue
            if _poly_deg(g) == d:
                factors_d = [g]
            else:
                factors_d = _equal_degree_factorization(g, d, p, rng)

            for _ in range(e):
                out.extend([poly_monic(h, p) for h in factors_d])

    out = lin_factors + out

    # Final sanity: multiply all factors = original f0
    prod = np.array([1], dtype=np.int64)
    for h in out:
        prod = poly_mul(prod, h, p)
    prod = poly_monic(prod, p)
    if not np.array_equal(prod, f0):
        raise RuntimeError("factor_poly_over_fp sanity check failed: product(factors) != f0")

    return out



def _poly_eval_at(f: np.ndarray, x: int, p: int) -> int:
    """Evaluate f(x) mod p, coeffs low->high."""
    x = int(x) % p
    acc = 0
    for a in reversed(poly_trim(f)):
        acc = (acc * x + int(a)) % p
    return acc


if __name__ == "__main__":
    # Basic tests

    p = 2
    # f = (x+1)^3 * (x^2+x+1)
    x1 = np.array([1, 1], dtype=np.int64)              # x+1
    q2 = np.array([1, 1, 1], dtype=np.int64)           # x^2+x+1 irreducible over GF(2)
    f = poly_mul(poly_mul(poly_mul(x1, x1, p), x1, p), q2, p)
    facs = factor_poly_over_fp(f, p)
    # check multiplicities
    keys = [tuple(poly_monic(g, p).tolist()) for g in facs]
    assert keys.count(tuple(x1.tolist())) == 3
    assert keys.count(tuple(poly_monic(q2, p).tolist())) == 1
    print("factor_poly_over_fp tests passed")

    p = 2
    F = np.eye(4, dtype=np.int64)
    mF = minimal_polynomial(F, p)
    # identity has minimal polynomial (x-1) = x+1 in p=2 => coeffs [1,1]
    assert np.array_equal(mF, np.array([1, 1], dtype=np.int64))
    print("minpoly.py tests passed")
