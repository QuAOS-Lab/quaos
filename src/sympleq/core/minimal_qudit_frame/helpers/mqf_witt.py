import numpy as np
from typing import List, Tuple, Dict, Any
from sympleq.core.symmetries.modular_helpers import mod_p, rank_mod


def witt_decompose_form(B: np.ndarray, p: int) -> Tuple[np.ndarray, List[Tuple[str, int, int]], Dict[str, Any]]:
    """
    Deterministic Witt-style decomposition for a bilinear form matrix B over GF(p), p odd.

    Returns
    -------
    P : (m,m) ndarray over GF(p)
        Invertible change-of-basis matrix such that B' = P^T B P is in a block-diagonal Witt normal form
        (up to the level needed by the library/tests).
    blocks : list of (type, i, j)
        Describes how basis coordinates 0..m-1 of the transformed basis are grouped.
        - ("hyp", i, j): a 2D hyperbolic plane spanning indices i and j (always j=i+1 here)
        - ("ani", i, i): a 1D anisotropic line at index i (symmetric forms only)
        - ("rad", i, i): a 1D radical vector at index i
        This *matches your tests*: they add both i and j only when type=="hyp", else just i.
    info : dict
        Diagnostics (rank, form type, counts).
    """
    B = mod_p(np.asarray(B, dtype=np.int64), p)
    m = B.shape[0]
    if B.ndim != 2 or B.shape[1] != m:
        raise ValueError(f"B must be square, got {B.shape}.")
    if p <= 1:
        raise ValueError("p must be a prime >= 2.")
    if p == 2:
        raise NotImplementedError(
            "witt_decompose_form for p=2 is not a bilinear Witt decomposition problem "
            "(quadratic refinements matter). Use the dedicated p=2 routines."
        )

    Z = np.zeros_like(B)
    is_alt = np.array_equal(mod_p(B + B.T, p), Z) and np.all(np.diag(B) % p == 0)
    is_sym = np.array_equal(mod_p(B - B.T, p), Z)
    if not (is_alt or is_sym):
        raise ValueError("B must be symmetric (B=B^T) or alternating (B=-B^T with zero diagonal).")

    # Working congruence reduction
    P = np.eye(m, dtype=np.int64)
    Bm = B.copy()

    def swap_basis(i: int, j: int) -> None:
        if i == j:
            return
        Bm[[i, j], :] = Bm[[j, i], :]
        Bm[:, [i, j]] = Bm[:, [j, i]]
        P[:, [i, j]] = P[:, [j, i]]

    def scale_basis(i: int, s: int) -> None:
        s %= p
        if s == 1:
            return
        Bm[:, i] = mod_p(Bm[:, i] * s, p)
        Bm[i, :] = mod_p(Bm[i, :] * s, p)
        P[:, i] = mod_p(P[:, i] * s, p)

    def shear_basis(i: int, j: int, a: int) -> None:
        """
        Congruence shear corresponding to basis update e_i <- e_i + a e_j.
        """
        a %= p
        if a == 0:
            return
        Bm[:, i] = mod_p(Bm[:, i] + a * Bm[:, j], p)
        Bm[i, :] = mod_p(Bm[i, :] + a * Bm[j, :], p)
        P[:, i] = mod_p(P[:, i] + a * P[:, j], p)

    blocks: List[Tuple[str, int, int]] = []
    k = 0

    # For symmetric hyperbolic normalization we need 1/2 mod p
    inv2 = pow(2, p - 2, p)  # p odd so OK

    while k < m:
        # Find a pivot row in the remaining subspace with some nonzero coupling.
        piv = None
        for i in range(k, m):
            if np.any(Bm[i, k:] % p):
                piv = i
                break

        if piv is None:
            # Remaining part is radical.
            for idx in range(k, m):
                blocks.append(("rad", idx, idx))
            break

        if piv != k:
            swap_basis(piv, k)

        if is_sym and int(Bm[k, k] % p) != 0:
            # 1D anisotropic line: eliminate couplings to the rest.
            a = int(Bm[k, k] % p)
            inva = pow(a, p - 2, p)
            for t in range(k + 1, m):
                bkt = int(Bm[k, t] % p)
                if bkt != 0:
                    # t <- t - (b(u,t)/b(u,u)) u
                    shear_basis(t, k, (-bkt * inva) % p)
            blocks.append(("ani", k, k))
            k += 1
            continue

        # Otherwise we try to form a hyperbolic plane with u := e_k isotropic.
        j = None
        for jj in range(k + 1, m):
            if int(Bm[k, jj] % p) != 0:
                j = jj
                break

        if j is None:
            # u is orthogonal to all remaining vectors -> radical direction.
            blocks.append(("rad", k, k))
            k += 1
            continue

        if j != k + 1:
            swap_basis(j, k + 1)

        # Normalize <u,v>=1 by scaling v := e_{k+1}
        a = int(Bm[k, k + 1] % p)
        inva = pow(a, p - 2, p)
        scale_basis(k + 1, inva)  # now Bm[k,k+1]=1

        if is_sym:
            # Make v isotropic (so the 2x2 block is exactly hyperbolic):
            # v <- v - (B(v,v)/2) u, since <u,v>=1 and <u,u>=0.
            bvv = int(Bm[k + 1, k + 1] % p)
            if bvv != 0:
                half = (bvv * inv2) % p
                shear_basis(k + 1, k, (-half) % p)

        # Kill couplings of u and v to the remaining basis vectors.
        for t in range(k + 2, m):
            bu = int(Bm[k, t] % p)       # <u, w>
            bv = int(Bm[k + 1, t] % p)   # <v, w>
            if bu == 0 and bv == 0:
                continue
            if is_alt:
                # alternating case: w <- w + <v,w> u - <u,w> v
                shear_basis(t, k, bv)
                shear_basis(t, k + 1, (-bu) % p)
            else:
                # symmetric case:   w <- w - <v,w> u - <u,w> v
                shear_basis(t, k, (-bv) % p)
                shear_basis(t, k + 1, (-bu) % p)

        blocks.append(("hyp", k, k + 1))
        k += 2

    # Diagnostics
    rankB = int(rank_mod(B, p))
    rad_dim = m - rankB
    info = {
        "form": "alternating" if is_alt else "symmetric",
        "rank": rankB,
        "radical_dim": int(rad_dim),
        "n_hyp": int(sum(t == "hyp" for (t, _, __) in blocks)),
        "n_ani": int(sum(t == "ani" for (t, _, __) in blocks)),
        "n_rad": int(sum(t == "rad" for (t, _, __) in blocks)),
    }
    return mod_p(P, p), blocks, info
