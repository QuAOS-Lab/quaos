from __future__ import annotations
import numpy as np
from .modular_helpers import mod_p, omega_matrix, inv_mod_scalar, independent_columns


def symplectic_gram_schmidt_from_span(B: np.ndarray, p: int) -> np.ndarray:
    """
    Deterministically build a hyperbolic basis [U|V] from a nondegenerate even-dim span.
    Input: B (2n x m) whose columns span W.
    Output: T (2n x 2k) with T^T Ω T = Ω_k and col(T)=W.
    """
    B = independent_columns(B, p)
    n2, m = B.shape
    if m == 0:
        return np.zeros((n2, 0), dtype=np.int64)
    if m % 2 != 0:
        raise ValueError("Span dimension must be even for symplectic Gram-Schmidt.")
    Ω = omega_matrix(n2 // 2, p)

    # Work list of vectors that span W
    V = [mod_p(B[:, i: i + 1], p) for i in range(m)]

    U_list: list[np.ndarray] = []
    V_list: list[np.ndarray] = []

    def pair(u: np.ndarray, v: np.ndarray) -> int:
        return int(mod_p((u.T @ Ω @ v).reshape(()), p))

    # Deterministic: iterate in given order, find first nonzero pairing
    while V:
        u = V.pop(0)
        if np.all(u % p == 0):
            continue

        # find v with <u,v> != 0
        j = None
        beta = 0
        for idx, cand in enumerate(V):
            beta = pair(u, cand)
            if beta != 0:
                j = idx
                break
        if j is None:
            # If W is nondegenerate this should not happen
            raise RuntimeError("Failed to find symplectic partner; span appears degenerate.")

        v = V.pop(j)

        # Normalize v so that <u,v>=1
        beta_inv = inv_mod_scalar(beta, p)
        v = mod_p(v * beta_inv, p)

        # Orthogonalize remaining vectors against (u,v)
        newV = []
        for w in V:
            au = pair(w, v)          # <w,v>
            av = pair(u, w)          # <u,w> = -<w,u> but sign irrelevant mod p=2
            # w <- w - au*u + av*v  (standard hyperbolic elimination)
            w2 = mod_p(w - au * u + av * v, p)
            newV.append(w2)
        V = newV

        U_list.append(u)
        V_list.append(v)

    T = np.concatenate(U_list + V_list, axis=1)  # [U|V]
    # Optional sanity check could go here; keep deterministic and cheap
    return mod_p(T, p)
