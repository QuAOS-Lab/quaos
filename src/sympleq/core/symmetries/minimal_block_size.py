import numpy as np
from typing import Dict, List, Tuple
from .modular_helpers import mod_p, rank_mod, _solve_linear, nullspace_mod

# --------------------------
# Polynomial arithmetic over GF(p): coeffs low->high, monic when expected
# --------------------------


def poly_trim(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=np.int64)
    i = len(a) - 1
    while i > 0 and a[i] == 0:
        i -= 1
    return a[:i + 1]


def poly_monic(a: np.ndarray, p: int) -> np.ndarray:
    a = mod_p(poly_trim(a), p)
    if a.size == 0 or (a.size == 1 and a[0] == 0):
        return a
    inv = pow(int(a[-1]), p - 2, p)
    return mod_p(a * inv, p)


def poly_add(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    n = max(len(a), len(b))
    c = np.zeros(n, dtype=np.int64)
    c[:len(a)] += a
    c[:len(b)] += b
    return mod_p(poly_trim(c), p)


def poly_sub(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    n = max(len(a), len(b))
    c = np.zeros(n, dtype=np.int64)
    c[:len(a)] += a
    c[:len(b)] -= b
    return mod_p(poly_trim(c), p)


def poly_mul(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    if len(a) == 0 or len(b) == 0:
        return np.zeros(1, dtype=np.int64)
    c = np.zeros(len(a) + len(b) - 1, dtype=np.int64)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        c[i:i + len(b)] += ai * b
    return mod_p(poly_trim(c), p)


def poly_divmod(a: np.ndarray, b: np.ndarray, p: int) -> Tuple[np.ndarray, np.ndarray]:
    a = mod_p(poly_trim(a.copy()), p)
    b = mod_p(poly_trim(b.copy()), p)
    if len(b) == 1 and b[0] == 0:
        raise ZeroDivisionError("poly div by 0")
    if len(a) < len(b):
        return np.zeros(1, dtype=np.int64), a
    inv_lead = pow(int(b[-1]), p - 2, p)
    q = np.zeros(len(a) - len(b) + 1, dtype=np.int64)
    r = a.copy()
    while len(r) >= len(b) and not (len(r) == 1 and r[0] == 0):
        k = len(r) - len(b)
        c = mod_p(r[-1] * inv_lead, p)
        q[k] = c
        r[k: k + len(b)] -= c * b
        r = mod_p(poly_trim(r), p)
        if len(r) == 0:
            r = np.zeros(1, dtype=np.int64)
            break
        if len(r) < len(b):
            break
    return mod_p(poly_trim(q), p), mod_p(poly_trim(r), p)


def poly_mod(a: np.ndarray, m: np.ndarray, p: int) -> np.ndarray:
    _, r = poly_divmod(a, m, p)
    return r


def poly_gcd(a: np.ndarray, b: np.ndarray, p: int) -> np.ndarray:
    a = mod_p(poly_trim(a), p)
    b = mod_p(poly_trim(b), p)
    while not (len(b) == 1 and b[0] == 0):
        _, r = poly_divmod(a, b, p)
        a, b = b, r
    return poly_monic(a, p)


def poly_pow_mod(a: np.ndarray, e: int, m: np.ndarray, p: int) -> np.ndarray:
    res = np.array([1], dtype=np.int64)
    base = a.copy()
    while e > 0:
        if e & 1:
            res = poly_mod(poly_mul(res, base, p), m, p)
        base = poly_mod(poly_mul(base, base, p), m, p)
        e >>= 1
    return poly_monic(res, p)


def poly_eval_matrix(F: np.ndarray, poly: np.ndarray, p: int) -> np.ndarray:
    n2 = F.shape[0]
    M = np.zeros_like(F)
    P_pow = np.eye(n2, dtype=np.int64)
    for a in poly:
        if a != 0:
            M = mod_p(M + a * P_pow, p)
        P_pow = mod_p(P_pow @ F, p)
    return M


def poly_reciprocal(poly: np.ndarray, p: int) -> np.ndarray:
    poly = poly_trim(poly)
    rev = poly[::-1].copy()
    return poly_monic(rev, p)


# --------------------------
# Distinct-degree + equal-degree factorization (Cantor–Zassenhaus)
# --------------------------

def squarefree_decomposition(f: np.ndarray, p: int) -> List[Tuple[np.ndarray, int]]:
    f = poly_monic(f, p)
    df = derivative_poly(f, p)
    if np.all(df == 0):
        g = poly_pth_root(f, p)
        sub = squarefree_decomposition(g, p)
        return [(h, e * p) for (h, e) in sub]
    g = poly_gcd(f, df, p)
    w, _ = poly_divmod(f, g, p)
    i, res = 1, []
    while not (len(w) == 1 and w[0] == 0):
        y = poly_gcd(w, g, p)
        fi, _ = poly_divmod(w, y, p)
        if not (len(fi) == 1 and fi[0] == 1):
            res.append((fi, i))
        w = y
        i += 1
        if len(g) == 1 and g[0] == 1:
            break
        g = poly_mod(g, np.array([0, 1], dtype=np.int64), p)
    if len(g) != 1 or g[0] != 1:
        g_root = poly_pth_root(g, p)
        sub = squarefree_decomposition(g_root, p)
        res += [(h, e * p) for (h, e) in sub]
    return res


def derivative_poly(f: np.ndarray, p: int) -> np.ndarray:
    if len(f) <= 1:
        return np.zeros(1, dtype=np.int64)
    df = np.array([(i * f[i]) % p for i in range(1, len(f))], dtype=np.int64)
    return poly_trim(df)


def poly_pth_root(f: np.ndarray, p: int) -> np.ndarray:
    g = np.zeros((len(f) + p - 1) // p, dtype=np.int64)
    for i in range(0, len(f), p):
        g[i // p] = f[i]
    return poly_trim(g)


def distinct_degree_factorization(f: np.ndarray, p: int) -> List[Tuple[np.ndarray, int]]:
    f = poly_monic(f, p)
    res = []
    x = np.array([0, 1], dtype=np.int64)
    h = x.copy()
    n = len(f) - 1
    for d in range(1, n + 1):
        h = poly_pow_mod(h, p, f, p)
        g = poly_gcd(poly_sub(h, x, p), f, p)
        if len(g) > 1 and not np.array_equal(g, np.array([1], dtype=np.int64)) and not np.array_equal(g, f):
            res.append((poly_monic(g, p), d))
            q, r = poly_divmod(f, g, p)
            assert np.all(r == 0)
            f = poly_monic(q, p)
            h = poly_mod(h, f, p)
        if len(f) <= 1 or np.array_equal(f, np.array([1], dtype=np.int64)):
            break
    if len(f) > 1 and not np.array_equal(f, np.array([1], dtype=np.int64)):
        res.append((f, len(f) - 1))
    return res


def equal_degree_factorization(f: np.ndarray, d: int, p: int, rng=np.random.default_rng(2025)) -> List[np.ndarray]:
    if len(f) - 1 == d:
        return [f]
    while True:
        a = mod_p(rng.integers(0, p, size=len(f) - 1, dtype=np.int64), p)
        a = poly_trim(a)
        if len(a) == 0:
            a = np.array([1], dtype=np.int64)
        exp = (p**d - 1) // 2
        a_pow = poly_pow_mod(a, exp, f, p)
        g = poly_gcd(poly_sub(a_pow, np.array([1], dtype=np.int64), p), f, p)
        if len(g) > 1 and not np.array_equal(g, f):
            return (equal_degree_factorization(g, d, p, rng) +
                    equal_degree_factorization(poly_monic(poly_divmod(f, g, p)[0], p), d, p, rng))


def factor_poly_over_fp(f: np.ndarray, p: int) -> List[np.ndarray]:
    f = poly_monic(f, p)
    if len(f) <= 1:
        return [f]
    sq = squarefree_decomposition(f, p)
    out = []
    for fi, ei in sq:
        dd = distinct_degree_factorization(fi, p)
        for hd, d in dd:
            parts = equal_degree_factorization(hd, d, p)
            for _ in range(ei):
                out.extend(parts)
    return [poly_monic(h, p) for h in out]


# --------------------------
# Krylov cyclic decomposition: companion polynomials of cyclic blocks
# --------------------------

def minimal_poly_for_vector(F: np.ndarray, v: np.ndarray, p: int) -> np.ndarray:
    v = mod_p(v.reshape(-1, 1), p)
    V = v.copy()
    powers = [v]
    while True:
        r = rank_mod(V, p)
        w = mod_p(F @ powers[-1], p)
        candidate = np.concatenate([V, w], axis=1)
        r2 = rank_mod(candidate, p)
        if r2 == r:
            c = _solve_linear(mod_p(V, p), mod_p(-w, p), p)
            coeffs = np.concatenate([c.reshape(-1), np.array([1], dtype=np.int64)])
            return poly_monic(coeffs, p)
        V = candidate
        powers.append(w)


def cyclic_decomposition(F: np.ndarray, p: int) -> List[Tuple[np.ndarray, np.ndarray]]:
    n2 = F.shape[0]
    basis_blocks: List[np.ndarray] = []
    polys: List[np.ndarray] = []
    B = np.zeros((n2, 0), dtype=np.int64)

    seeds = [np.eye(n2, dtype=np.int64)[:, i: i + 1] for i in range(n2)]
    rng = np.random.default_rng(2025)
    seeds += [rng.integers(0, p, size=(n2, 1), dtype=np.int64) for _ in range(min(8, n2))]

    s_idx = 0
    while B.shape[1] < n2:
        if s_idx < len(seeds):
            v = mod_p(seeds[s_idx], p)
            s_idx += 1
        else:
            v = rng.integers(0, p, size=(n2, 1), dtype=np.int64)

        if rank_mod(np.concatenate([B, v], axis=1), p) == B.shape[1]:
            continue

        m_v = minimal_poly_for_vector(F, v, p)
        deg = len(m_v) - 1

        block_cols = []
        C = B.copy()
        w = v.copy()
        steps_cap = deg + n2
        steps = 0
        while steps < steps_cap and len(block_cols) < deg:
            steps += 1
            if rank_mod(np.concatenate([C, w], axis=1), p) > C.shape[1]:
                block_cols.append(w.copy())
                C = np.concatenate([C, w], axis=1)
            w = mod_p(F @ w, p)

        if not block_cols:
            continue

        block = np.concatenate(block_cols, axis=1)
        B = np.concatenate([B, block], axis=1)
        basis_blocks.append(block)
        polys.append(m_v)

    assert B.shape[1] == n2, "Cyclic decomposition failed to span the space"
    return list(zip(basis_blocks, polys))


# --------------------------
# Symplectic-aware ±1 floor helpers
# --------------------------

def is_t_minus_lambda(q: np.ndarray, p: int):
    """
    Detect q(t)=t-1 or q(t)=t+1 (monic, low->high).
    Returns lam=+1 for t-1, lam=-1 for t+1, else None.
    Over p=2 we only get lam=+1 since +1=-1.
    """
    q = poly_monic(q, p)
    if len(q) != 2:
        return None
    a0, a1 = int(q[0] % p), int(q[1] % p)
    if a1 != 1:
        return None
    if a0 == (p - 1) % p:   # t - 1
        return 1
    if a0 == 1 % p and p != 2:  # t + 1 (distinct only if p!=2)
        return -1
    return None


def symplectic_half_floor_pm1(F: np.ndarray, p: int, lam: int, m_lam: int) -> int:
    """
    Symplectic-aware half-dimension floor for q=t-lam (lam=+1 or -1).

    Let K_j = ker((F-lam I)^j). The smallest j for which
    rank(K_j^T Ω K_j) >= 2 gives the earliest appearance of a hyperbolic pair.
    Any nondegenerate invariant symplectic subspace in this primary must have
    half-dimension >= j.
    """
    n2 = F.shape[0]
    n = n2 // 2
    # standard Omega
    Omega = np.block([
        [np.zeros((n, n), dtype=np.int64), np.eye(n, dtype=np.int64)],
        [-np.eye(n, dtype=np.int64), np.zeros((n, n), dtype=np.int64)]
    ])
    Omega = mod_p(Omega, p)

    A = mod_p(F - lam * np.eye(n2, dtype=np.int64), p)
    Aj = np.eye(n2, dtype=np.int64)

    best_j = m_lam
    for j in range(1, m_lam + 1):
        Aj = mod_p(Aj @ A, p)          # Aj = (F - lam I)^j
        Kj = nullspace_mod(Aj, p)      # basis for generalized kernel
        if Kj.shape[1] < 2:
            continue
        Gj = mod_p(Kj.T @ Omega @ Kj, p)
        if rank_mod(Gj, p) >= 2:       # at least one hyperbolic pair
            best_j = j
            break

    return int(best_j)


# --------------------------
# Primary components via irreducible factors and exponents
# --------------------------

def primary_components_from_cyclic(F: np.ndarray, p: int) -> Dict:
    """
    Build primary components V^(q) for irreducible q(t) and sector grouping.
    Symplectic-aware half-dimension floors are applied for q=t±1.
    """
    n2 = F.shape[0]
    cyclic_blocks = cyclic_decomposition(F, p)

    mult: Dict[Tuple[int, ...], int] = {}
    irreducible_polys: Dict[Tuple[int, ...], np.ndarray] = {}

    for _, m_i in cyclic_blocks:
        factors = factor_poly_over_fp(m_i, p)
        uniq = {}
        for q in factors:
            key = tuple(q.tolist())
            irreducible_polys[key] = q
            uniq.setdefault(key, 0)
            uniq[key] += 1
        for key, e in uniq.items():
            mult[key] = max(mult.get(key, 0), e)

    primaries: Dict[Tuple[int, ...], Dict] = {}
    for key, e_max in mult.items():
        q = irreducible_polys[key]
        Mq = poly_eval_matrix(F, q, p)

        # compute q(F)^{e_max}
        P = np.eye(n2, dtype=np.int64)
        base = Mq.copy()
        e = e_max
        while e > 0:
            if e & 1:
                P = mod_p(P @ base, p)
            base = mod_p(base @ base, p)
            e >>= 1

        Vq = nullspace_mod(P, p)
        primaries[key] = {
            "poly": q,
            "deg": len(q) - 1,
            "exponent": e_max,
            "V_basis": Vq
        }

    keys = list(primaries.keys())
    for key in keys:
        q = primaries[key]["poly"]
        q_star = poly_reciprocal(q, p)
        k_star = tuple(q_star.tolist())
        primaries[key]["reciprocal_key"] = k_star
        primaries[key]["self_reciprocal"] = (k_star == key)

    used = set()
    sectors = []
    for key in primaries.keys():
        if key in used:
            continue
        data = primaries[key]

        # default RCF-scale floor
        half_floor = data["deg"] * data["exponent"]

        # symplectic refinement for q=t±1
        lam = is_t_minus_lambda(data["poly"], p)
        if lam is not None:
            half_floor = symplectic_half_floor_pm1(F, p, lam, data["exponent"])

        if data["self_reciprocal"]:
            W = data["V_basis"]
            sectors.append({
                "type": "self",
                "q_key": key,
                "W_basis": W,
                "half_dim_floor": int(half_floor)
            })
            used.add(key)
        else:
            k_star = data["reciprocal_key"]
            if k_star not in primaries:
                W = data["V_basis"]
                sectors.append({
                    "type": "self",  # fallback robustness
                    "q_key": key,
                    "W_basis": W,
                    "half_dim_floor": int(half_floor)
                })
                used.add(key)
            else:
                W = np.concatenate([data["V_basis"], primaries[k_star]["V_basis"]], axis=1)
                sectors.append({
                    "type": "paired",
                    "q_key": key,
                    "qstar_key": k_star,
                    "W_basis": W,
                    "half_dim_floor": int(half_floor)
                })
                used.add(key)
                used.add(k_star)

    L_min_star = max([sec["half_dim_floor"] for sec in sectors] + [0])

    return {
        "irreducibles": primaries,
        "sectors": sectors,
        "Lmin_star": int(L_min_star)
    }


# --------------------------
# Public API
# --------------------------

def rcf_prepass(F: np.ndarray, p: int) -> Dict:
    """
    Deterministic structural pre-pass for symplectic F over GF(p).
    - Computes cyclic decomposition via Krylov subspaces.
    - Factors block polynomials over GF(p).
    - Assembles primary components V^(q) = ker(q(F)^{m_q}).
    - Pairs reciprocals q and q* into sectors W_q (or self-recognizes q=q*).
    - Applies symplectic-aware half-dimension floors for q=t±1 primaries.
    Returns per-irreducible metadata, sector bases, and a half-dimension lower bound Lmin_star.
    """
    assert F.shape[0] == F.shape[1], "F must be square"
    return primary_components_from_cyclic(F, p)
