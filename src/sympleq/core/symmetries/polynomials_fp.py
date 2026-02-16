from __future__ import annotations
import numpy as np

from sympleq.core.symmetries.modular_helpers import mod_p


def poly_coerce_1d(a) -> np.ndarray:
    arr = np.asarray(a)
    if arr.ndim == 2 and 1 in arr.shape:
        arr = arr.reshape(-1)
    if arr.ndim != 1:
        raise TypeError(f"Polynomial must be 1D coeff vector (low->high). Got shape={arr.shape}")
    return arr.astype(np.int64, copy=False)


def poly_trim(a) -> np.ndarray:
    a = poly_coerce_1d(a)
    i = len(a) - 1
    while i > 0 and int(a[i]) == 0:
        i -= 1
    return a[: i + 1]


def poly_is_zero(a) -> bool:
    a = poly_trim(a)
    return len(a) == 1 and int(a[0]) == 0


def poly_monic(a, p: int) -> np.ndarray:
    a = mod_p(poly_trim(a), p)
    if poly_is_zero(a):
        return np.array([0], dtype=np.int64)
    inv = pow(int(a[-1]) % p, -1, p)
    return mod_p(a * inv, p)


def poly_add(a, b, p: int) -> np.ndarray:
    a = poly_coerce_1d(a)
    b = poly_coerce_1d(b)
    n = max(len(a), len(b))
    c = np.zeros(n, dtype=np.int64)
    c[: len(a)] += a
    c[: len(b)] += b
    return mod_p(poly_trim(c), p)


def poly_sub(a, b, p: int) -> np.ndarray:
    a = poly_coerce_1d(a)
    b = poly_coerce_1d(b)
    n = max(len(a), len(b))
    c = np.zeros(n, dtype=np.int64)
    c[: len(a)] += a
    c[: len(b)] -= b
    return mod_p(poly_trim(c), p)


def poly_mul(a, b, p: int) -> np.ndarray:
    a = poly_coerce_1d(a)
    b = poly_coerce_1d(b)
    if poly_is_zero(a) or poly_is_zero(b):
        return np.array([0], dtype=np.int64)
    c = np.zeros(len(a) + len(b) - 1, dtype=np.int64)
    for i, ai in enumerate(a):
        ai = int(ai) % p
        if ai == 0:
            continue
        c[i:i + len(b)] += ai * b
    return mod_p(poly_trim(c), p)


def poly_divmod(a, b, p: int) -> tuple[np.ndarray, np.ndarray]:
    a = mod_p(poly_trim(a.copy()), p)
    b = mod_p(poly_trim(b.copy()), p)
    if poly_is_zero(b):
        raise ZeroDivisionError("poly_divmod: division by zero")
    if len(a) < len(b):
        return np.array([0], dtype=np.int64), a

    inv_lead = pow(int(b[-1]) % p, -1, p)
    q = np.zeros(len(a) - len(b) + 1, dtype=np.int64)
    r = a.copy()

    while len(r) >= len(b) and not poly_is_zero(r):
        k = len(r) - len(b)
        c = (int(r[-1]) * inv_lead) % p
        if c:
            q[k] = c
            r[k:k + len(b)] -= c * b
        r = poly_trim(mod_p(r, p))

    q = poly_trim(mod_p(q, p))
    r = poly_trim(mod_p(r, p))
    return q, r


def poly_gcd(a, b, p: int) -> np.ndarray:
    a = poly_monic(a, p)
    b = poly_monic(b, p)
    while not poly_is_zero(b):
        _, r = poly_divmod(a, b, p)
        a, b = b, r
    return poly_monic(a, p)


def poly_lcm(a, b, p: int) -> np.ndarray:
    a = poly_monic(a, p)
    b = poly_monic(b, p)
    g = poly_gcd(a, b, p)
    q, r = poly_divmod(poly_mul(a, b, p), g, p)
    if not poly_is_zero(r):
        raise RuntimeError("poly_lcm: unexpected nonzero remainder")
    return poly_monic(q, p)


def poly_xgcd(a, b, p: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extended gcd: returns (s,t,g) with s*a + t*b = g, g monic.
    """
    a = mod_p(poly_trim(a), p)
    b = mod_p(poly_trim(b), p)

    r0, r1 = a, b
    s0, s1 = np.array([1], dtype=np.int64), np.array([0], dtype=np.int64)
    t0, t1 = np.array([0], dtype=np.int64), np.array([1], dtype=np.int64)

    while not poly_is_zero(r1):
        q, r = poly_divmod(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, poly_sub(s0, poly_mul(q, s1, p), p)
        t0, t1 = t1, poly_sub(t0, poly_mul(q, t1, p), p)

    if poly_is_zero(r0):
        return np.array([0], dtype=np.int64), np.array([0], dtype=np.int64), np.array([0], dtype=np.int64)

    lead = int(r0[-1]) % p
    inv = pow(lead, -1, p)
    g = mod_p(r0 * inv, p)
    s = mod_p(s0 * inv, p)
    t = mod_p(t0 * inv, p)
    return poly_trim(s), poly_trim(t), poly_trim(g)


def poly_pow(a, e: int, p: int) -> np.ndarray:
    a = poly_monic(a, p)
    res = np.array([1], dtype=np.int64)
    base = a.copy()
    ee = int(e)
    while ee > 0:
        if ee & 1:
            res = poly_mul(res, base, p)
        base = poly_mul(base, base, p)
        ee >>= 1
    return poly_monic(res, p)


def poly_reciprocal(q, p: int) -> np.ndarray:
    q = poly_monic(q, p)
    return poly_monic(q[::-1].copy(), p)


def poly_eval_matrix(F: np.ndarray, poly: np.ndarray, p: int) -> np.ndarray:
    """
    Evaluate poly(F) for column-action convention, coeffs low->high.
    """
    F = mod_p(F, p)
    poly = mod_p(poly, p)
    n2 = F.shape[0]
    M = np.zeros((n2, n2), dtype=np.int64)
    P = np.eye(n2, dtype=np.int64)
    for a in poly:
        aa = int(a) % p
        if aa:
            M = mod_p(M + aa * P, p)
        P = mod_p(P @ F, p)
    return M


if __name__ == "__main__":
    p = 2
    a = np.array([1, 1, 1], dtype=np.int64)  # 1 + x + x^2
    b = np.array([1, 1], dtype=np.int64)     # 1 + x
    q, r = poly_divmod(a, b, p)
    assert np.array_equal(poly_add(poly_mul(q, b, p), r, p), poly_monic(a, p))
    g = poly_gcd(a, b, p)
    assert np.array_equal(g, np.array([1], dtype=np.int64))
    print("polynomials_fp.py tests passed")
