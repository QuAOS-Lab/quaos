# sympleq/core/symmetries/atomic_decomposition_helpers/atomic_unipotent_p2_generators.py
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

from ..modular_helpers import mod_p, inv_mod_mat, is_symplectic, rank_mod, nullspace_mod
from .atomic_linear import darboux_basis_from_span


BlockSpec = Tuple[str, int, int | None]


def _jordan_unipotent(k: int) -> np.ndarray:
    if k <= 0:
        raise ValueError("Jordan block size k must be positive.")
    U = np.eye(k, dtype=np.int64)
    for i in range(k - 1):
        U[i, i + 1] = 1
    return U % 2



def _invariant_alternating_form_for_jordan_p2(m: int) -> np.ndarray:
    """
    Find a nondegenerate alternating form G preserved by one unipotent
    Jordan block of size m over GF(2).  Used to build V_beta(2k) fixtures.
    """
    m = int(m)
    if m <= 0 or (m % 2):
        raise ValueError("V_beta block size must be positive and even.")

    U = _jordan_unipotent(m)
    pairs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    nvar = len(pairs)
    mats = []
    for (i, j) in pairs:
        G = np.zeros((m, m), dtype=np.int64)
        G[i, j] = 1
        G[j, i] = 1
        mats.append(G)

    rows = []
    for a in range(m):
        for b in range(m):
            coeff = [int((U.T @ G @ U + G)[a, b] % 2) for G in mats]
            if any(coeff):
                rows.append(coeff)
    A = np.array(rows, dtype=np.int64) if rows else np.zeros((0, nvar), dtype=np.int64)
    Z = nullspace_mod(A, 2)  # nvar x d
    if Z.shape[1] == 0:
        raise RuntimeError("could not find invariant alternating form for Jordan block")

    def form_from_coeff(c: np.ndarray) -> np.ndarray:
        G = np.zeros((m, m), dtype=np.int64)
        for bit, M in zip(c.reshape(-1) % 2, mats):
            if int(bit):
                G ^= M
        return G % 2

    # Exhaustive for normal regression sizes; deterministic sparse search otherwise.
    d = int(Z.shape[1])
    if d <= 20:
        masks = range(1, 1 << d)
    else:
        basis_masks = [1 << i for i in range(d)]
        pair_masks = [(1 << i) | (1 << j) for i in range(d) for j in range(i + 1, min(d, i + 8))]
        masks = basis_masks + pair_masks + [(1 << d) - 1]

    for mask in masks:
        coeff = np.zeros((d, 1), dtype=np.int64)
        for i in range(d):
            if (int(mask) >> i) & 1:
                coeff[i, 0] = 1
        c = mod_p(Z @ coeff, 2)
        G = form_from_coeff(c)
        if rank_mod(G, 2) == m and np.all(np.diag(G) % 2 == 0):
            if np.array_equal(mod_p(U.T @ G @ U, 2), G % 2):
                return G

    raise RuntimeError("invariant alternating forms were all degenerate for this Jordan block")


def canonical_V_beta_block_p2(length: int, beta: int = 1) -> np.ndarray:
    """
    Canonical V_beta(2k)-style p=2 unipotent block in standard coordinates.

    ``length`` is the even single Jordan-chain length 2k and also the total
    matrix dimension.  Over GF(2), beta is 0 or 1.  The construction follows
    the Chapter-5 representative in its natural basis and conjugates the
    preserved alternating form to the standard symplectic form.
    """
    L = int(length)
    b = int(beta) & 1
    if L <= 0 or (L % 2):
        raise ValueError("V_beta(2k) requires an even positive chain length.")

    k = L // 2
    idxs = list(range(-(2 * k - 1), 2 * k, 2))
    dim = L

    def pos(i: int) -> int:
        return idxs.index(int(i))

    G = np.zeros((dim, dim), dtype=np.int64)
    for i in idxs:
        a = pos(i)
        c = pos(-i)
        if a != c:
            G[a, c] = 1
            G[c, a] = 1

    U = np.zeros((dim, dim), dtype=np.int64)

    def add(src_i: int, dst_i: int, coeff: int = 1) -> None:
        if coeff & 1:
            U[pos(dst_i), pos(src_i)] ^= 1

    if k == 1:
        add(-1, -1)
        add(-1, 1)
        add(1, 1)
    else:
        for i in idxs:
            if i == 2 * k - 1:
                add(i, i)
            elif i <= -3:
                for j in range(i, 2, 2):  # i, i+2, ..., 1
                    if j in idxs:
                        add(i, j)
                if 3 in idxs:
                    add(i, 3, b)
            else:
                # -1 <= i < 2k-1: v_i -> v_i + v_{i+2}
                add(i, i)
                if i + 2 in idxs:
                    add(i, i + 2)

    U %= 2
    if not np.array_equal(mod_p(U.T @ G @ U, 2), G % 2):
        raise RuntimeError("internal error: Chapter-5 V_beta block does not preserve its defining symplectic form")

    S = darboux_basis_from_span(G, np.eye(dim, dtype=np.int64), 2)
    Fstd = mod_p(inv_mod_mat(S, 2) @ U @ S, 2)
    if not is_symplectic(Fstd, 2):
        raise RuntimeError("internal error: constructed V_beta block is not symplectic in standard coordinates")
    return Fstd


def _symplectic_direct_sum_standard(F1: np.ndarray, F2: np.ndarray, p: int = 2) -> np.ndarray:
    """Symplectic direct sum for matrices in grouped [x...|z...] ordering."""
    F1 = mod_p(np.asarray(F1, dtype=np.int64), p)
    F2 = mod_p(np.asarray(F2, dtype=np.int64), p)
    if F1.size == 0:
        return F2
    if F2.size == 0:
        return F1
    n1 = F1.shape[0] // 2
    n2 = F2.shape[0] // 2
    A1, B1, C1, D1 = F1[:n1, :n1], F1[:n1, n1:], F1[n1:, :n1], F1[n1:, n1:]
    A2, B2, C2, D2 = F2[:n2, :n2], F2[:n2, n2:], F2[n2:, :n2], F2[n2:, n2:]
    Z12 = np.zeros((n1, n2), dtype=np.int64)
    Z21 = np.zeros((n2, n1), dtype=np.int64)
    A = np.block([[A1, Z12], [Z21, A2]])
    B = np.block([[B1, Z12], [Z21, B2]])
    C = np.block([[C1, Z12], [Z21, C2]])
    D = np.block([[D1, Z12], [Z21, D2]])
    F = mod_p(np.block([[A, B], [C, D]]), p)
    if not is_symplectic(F, p):
        raise RuntimeError("internal error: symplectic direct sum is not symplectic")
    return F

def canonical_W_block_p2(k: int) -> np.ndarray:
    """
    Canonical W(k)-style p=2 unipotent symplectic block in standard coordinates.

    The block has dimension 2k and consists of reciprocal cyclic halves
    U and U^{-T}: F = diag(U, U^{-T}).  This gives two Jordan chains of
    length k and is useful as a deterministic hard paired-chain regression case.
    """
    U = _jordan_unipotent(int(k))
    D = inv_mod_mat(U, 2).T % 2
    F = np.block([[U, np.zeros((k, k), dtype=np.int64)], [np.zeros((k, k), dtype=np.int64), D]]) % 2
    if not is_symplectic(F, 2):
        raise RuntimeError("internal error: constructed W(k) block is not symplectic")
    return F


def canonical_W_beta_block_p2(length: int, beta: int = 1) -> np.ndarray:
    """
    Canonical W_beta(2*l+1)-style p=2 unipotent block in standard symplectic coordinates.

    ``length`` is the odd Jordan-chain length 2*l+1 of each of the two cyclic
    halves, so the total matrix dimension is 2*length = 4*l+2.  For GF(2) the
    nonzero Chapter-5 parameter is beta=1.  beta=0 is the ordinary W(length)
    block.
    """
    L = int(length)
    b = int(beta) & 1
    if L <= 0 or L % 2 == 0:
        raise ValueError("W_beta(2*l+1) requires an odd positive length.")
    if L == 1 or b == 0:
        return canonical_W_block_p2(L)

    ell = (L - 1) // 2
    idxs = list(range(-2 * ell, 2 * ell + 1, 2))
    dim = 2 * L

    def pos(kind: str, i: int) -> int:
        t = idxs.index(int(i))
        return 2 * t + (0 if kind == "w" else 1)

    # Symplectic form in the Chapter-5 basis: (w_i, x_{-i}) = 1.
    G = np.zeros((dim, dim), dtype=np.int64)
    for i in idxs:
        a = pos("w", i)
        c = pos("x", -i)
        G[a, c] = 1
        G[c, a] = 1

    U = np.zeros((dim, dim), dtype=np.int64)

    def add_image(src_kind: str, src_i: int, dst_kind: str, dst_i: int, coeff: int = 1) -> None:
        if coeff & 1:
            U[pos(dst_kind, dst_i), pos(src_kind, src_i)] ^= 1

    # w_i images, using the corrected Chapter-5 W_beta formula.
    for i in idxs:
        if i == 2 * ell:
            add_image("w", i, "w", i)
        elif i == -4 and i in idxs:
            add_image("w", i, "w", -4)
            add_image("w", i, "w", -2)
            if 2 in idxs:
                add_image("w", i, "x", 2, b)
        elif i == -2 and i in idxs:
            add_image("w", i, "w", -2)
            add_image("w", i, "w", 0)
            if 2 in idxs:
                add_image("w", i, "x", 2, b)
        else:
            add_image("w", i, "w", i)
            if i + 2 in idxs:
                add_image("w", i, "w", i + 2)

    # x_i images.
    for i in idxs:
        if i == 2 * ell:
            add_image("x", i, "x", i)
            continue
        for j in range(i, 2 * ell + 1, 2):
            if j in idxs:
                add_image("x", i, "x", j)
        if i <= -2 and 2 in idxs:
            add_image("x", i, "w", 2, b)

    U %= 2
    if not np.array_equal(mod_p(U.T @ G @ U, 2), G % 2):
        raise RuntimeError("internal error: Chapter-5 W_beta block does not preserve its defining symplectic form")

    # Columns of S are a Darboux basis for G: S^T G S = Omega.  Convert the
    # Chapter-5 representative into standard grouped symplectic coordinates.
    S = darboux_basis_from_span(G, np.eye(dim, dtype=np.int64), 2)
    Fstd = mod_p(inv_mod_mat(S, 2) @ U @ S, 2)
    if not is_symplectic(Fstd, 2):
        raise RuntimeError("internal error: constructed W_beta block is not symplectic in standard coordinates")
    return Fstd


def canonical_unipotent_p2_block(block_type: str, length: int, beta: int = 0) -> np.ndarray:
    """
    Return a deterministic canonical p=2 unipotent test block where implemented.

    Implemented:
      - block_type == "W": returns a W(length)-style paired-chain block.
      - block_type in {"V_beta", "VBETA"}: returns a nonzero-beta V_beta(length) block,
        where length must be even.
      - block_type in {"W_beta", "WBETA"}: returns W_beta(length), where length
        must be odd; beta=0 reduces to W(length).
    """
    typ = str(block_type).strip().upper()
    if typ == "W":
        return canonical_W_block_p2(int(length))
    if typ in {"V_BETA", "VBETA", "V"}:
        return canonical_V_beta_block_p2(int(length), int(beta))
    if typ in {"W_BETA", "WBETA"}:
        return canonical_W_beta_block_p2(int(length), int(beta))
    raise ValueError(f"unknown p=2 unipotent block_type={block_type!r}")


def direct_sum_W_blocks_p2(lengths: Sequence[int]) -> np.ndarray:
    """
    Direct sum of W(k) blocks, returned in global standard symplectic ordering
    [all U halves | all V halves].
    """
    lengths = [int(k) for k in lengths]
    if not lengths:
        return np.zeros((0, 0), dtype=np.int64)
    U_blocks = [_jordan_unipotent(k) for k in lengths]
    D_blocks = [inv_mod_mat(U, 2).T % 2 for U in U_blocks]
    n = sum(lengths)
    U = np.zeros((n, n), dtype=np.int64)
    D = np.zeros((n, n), dtype=np.int64)
    off = 0
    for k, Ub, Db in zip(lengths, U_blocks, D_blocks):
        U[off:off + k, off:off + k] = Ub
        D[off:off + k, off:off + k] = Db
        off += k
    F = np.block([[U, np.zeros((n, n), dtype=np.int64)], [np.zeros((n, n), dtype=np.int64), D]]) % 2
    if not is_symplectic(F, 2):
        raise RuntimeError("internal error: direct sum W-block matrix is not symplectic")
    return F


def direct_sum_unipotent_p2_blocks(specs: Sequence[BlockSpec]) -> np.ndarray:
    """
    Build a direct sum of implemented p=2 unipotent test blocks.

    Each spec is ``(block_type, length, beta)``.  The beta entry is ignored for
    ``W`` blocks, so specs such as ``("W", 3, None)`` are accepted.  For
    beta-labelled families, ``None`` defaults to ``0`` to match
    ``canonical_unipotent_p2_block``.
    """
    if all(str(t).strip().upper() == "W" for (t, _L, _b) in specs):
        return direct_sum_W_blocks_p2([int(L) for (_t, L, _b) in specs])

    blocks = []
    for t, L, beta in specs:
        typ = str(t).strip().upper()
        beta_i = 0 if beta is None else int(beta)
        # canonical_unipotent_p2_block ignores beta_i for W blocks, but passing
        # a normalized integer keeps mixed direct sums simple and robust.
        blocks.append(canonical_unipotent_p2_block(typ, int(L), beta_i))

    F = np.zeros((0, 0), dtype=np.int64)
    for b in blocks:
        F = _symplectic_direct_sum_standard(F, b, 2)
    if not is_symplectic(F, 2):
        raise RuntimeError("direct-sum block ordering is not standard symplectic for these block types")
    return F % 2


def random_symplectic_matrix_p2(n: int, *, seed: int = 0, steps: int = 32) -> np.ndarray:
    """Generate a reproducible symplectic matrix over GF(2) by elementary factors."""
    rng = np.random.default_rng(seed)
    n = int(n)
    if n <= 0:
        return np.zeros((0, 0), dtype=np.int64)
    S = np.eye(2 * n, dtype=np.int64)

    for _ in range(int(steps)):
        kind = int(rng.integers(0, 3))
        if kind == 0:
            # GL elementary transvection A = I + E_ij, i != j.
            i, j = rng.choice(n, size=2, replace=False)
            A = np.eye(n, dtype=np.int64)
            A[int(i), int(j)] ^= 1
            AinvT = inv_mod_mat(A, 2).T % 2
            E = np.block([[A, np.zeros((n, n), dtype=np.int64)], [np.zeros((n, n), dtype=np.int64), AinvT]]) % 2
        elif kind == 1:
            R = rng.integers(0, 2, size=(n, n), dtype=np.int64)
            Sym = (R + R.T) % 2
            # In characteristic 2 the symmetric shear block may carry a nonzero
            # diagonal: symplectic transvections x -> x + <x,v> v require shear
            # blocks B = e_i e_i^T with nonzero diagonal. (R + R.T) always zeroes
            # the diagonal, so randomize it back to reach the full subgroup.
            np.fill_diagonal(Sym, rng.integers(0, 2, size=n, dtype=np.int64))
            E = np.block([[np.eye(n, dtype=np.int64), Sym], [np.zeros((n, n), dtype=np.int64), np.eye(n, dtype=np.int64)]]) % 2
        else:
            R = rng.integers(0, 2, size=(n, n), dtype=np.int64)
            Sym = (R + R.T) % 2
            # Same characteristic-2 diagonal correction for the lower shear block.
            np.fill_diagonal(Sym, rng.integers(0, 2, size=n, dtype=np.int64))
            E = np.block([[np.eye(n, dtype=np.int64), np.zeros((n, n), dtype=np.int64)], [Sym, np.eye(n, dtype=np.int64)]]) % 2
        S = (S @ E) % 2

    if not is_symplectic(S, 2):
        raise RuntimeError("internal error: generated matrix is not symplectic")
    return S


def random_symplectic_conjugate_p2(F: np.ndarray, *, seed: int = 0, steps: int = 32) -> Tuple[np.ndarray, np.ndarray]:
    """Return (S^{-1} F S, S) for a reproducible random symplectic S over GF(2)."""
    F = mod_p(np.asarray(F, dtype=np.int64), 2)
    if F.shape[0] != F.shape[1] or F.shape[0] % 2 != 0:
        raise ValueError(f"F must be square of even dimension, got {F.shape}")
    n = F.shape[0] // 2
    S = random_symplectic_matrix_p2(n, seed=seed, steps=steps)
    Fc = mod_p(inv_mod_mat(S, 2) @ F @ S, 2)
    if not is_symplectic(Fc, 2):
        raise RuntimeError("internal error: conjugated matrix is not symplectic")
    return Fc, S
