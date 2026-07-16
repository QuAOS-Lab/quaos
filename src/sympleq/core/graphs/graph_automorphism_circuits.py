# Auto-generated helper for circuit augmentation (true circuits from nullspace)

from typing import Optional, Iterable

import numpy as np


def _gf2_nullspace_basis(A: np.ndarray) -> list[np.ndarray]:
    """
    Return a basis for the nullspace of A over GF(2).

    Parameters
    ----------
    A : np.ndarray
        m x n matrix with entries in {0,1} (will be reduced mod 2).

    Returns
    -------
    basis : list[np.ndarray]
        List of length-r vectors x (shape (n,)) over GF(2) such that A x = 0.
        If nullity r=0, returns [].
    """
    A = (np.asarray(A, dtype=np.uint8) & 1).copy()
    m, n = A.shape
    row = 0
    pivots: list[int] = []
    pivot_rows: list[int] = [-1] * n  # pivot_rows[c] = r if c is pivot

    # Gauss-Jordan to RREF
    for col in range(n):
        if row >= m:
            break
        # find pivot row with 1 in this column
        piv = -1
        for r in range(row, m):
            if A[r, col] & 1:
                piv = r
                break
        if piv < 0:
            continue
        # swap into position
        if piv != row:
            A[[row, piv]] = A[[piv, row]]
        pivots.append(col)
        pivot_rows[col] = row

        # eliminate other rows
        for r in range(m):
            if r != row and (A[r, col] & 1):
                A[r, :] ^= A[row, :]

        row += 1

    pivot_set = set(pivots)
    free_cols = [c for c in range(n) if c not in pivot_set]
    if not free_cols:
        return []

    basis: list[np.ndarray] = []
    for f in free_cols:
        x = np.zeros(n, dtype=np.uint8)
        x[f] = 1
        # For each pivot column p at row r: x[p] = A[r, f]
        for p in pivots:
            r = pivot_rows[p]
            if r >= 0 and (A[r, f] & 1):
                x[p] = 1
        basis.append(x)
    return basis


def _gf2_vec_to_bitmask(x: np.ndarray) -> int:
    """Convert a GF(2) vector x (len M) to a Python int bitmask."""
    x = (np.asarray(x, dtype=np.uint8) & 1).reshape(-1)
    bm = 0
    # Using Python int shifts; OK for M up to a few thousand.
    for i, b in enumerate(x.tolist()):
        if b:
            bm |= (1 << i)
    return bm


def _bitcount(bm: int) -> int:
    return int(bm.bit_count())


def extract_circuits_from_nullspace_gf2(
    P: np.ndarray,
    *,
    max_nullity: int = 12,
    max_circuits: int = 5000,
) -> list[list[int]]:
    """
    Compute true matroid circuits (minimal dependent supports) from the nullspace of P^T.

    Here P is the M x d matrix of vectors (rows). Dependencies are a in ker(P^T) \ {0}
    and circuits are inclusion-minimal supports of such a.

    Only enabled for GF(2) (qubit labels); P is reduced mod 2 internally.

    Safety caps:
      - if nullity r > max_nullity, returns [] (skip augmentation)
      - if number of circuits exceeds max_circuits, returns the smallest max_circuits circuits by size
    """
    P2 = (np.asarray(P, dtype=np.uint8) & 1)
    M, d = P2.shape
    A = P2.T  # d x M
    basis = _gf2_nullspace_basis(A)
    r = len(basis)
    if r == 0:
        return []
    if r > int(max_nullity):
        return []

    basis_bm = [_gf2_vec_to_bitmask(v) for v in basis]

    # Enumerate all nonzero combinations (2^r - 1) via XOR
    deps: set[int] = set()
    for mask in range(1, 1 << r):
        bm = 0
        mm = mask
        bit = 0
        while mm:
            if mm & 1:
                bm ^= basis_bm[bit]
            mm >>= 1
            bit += 1
        if bm != 0:
            deps.add(bm)

    # Convert supports and keep inclusion-minimal ones (circuits)
    dep_list = sorted(deps, key=_bitcount)
    circuits: list[int] = []
    for bm in dep_list:
        # discard if any existing circuit is a subset
        is_min = True
        for c in circuits:
            if (c & bm) == c:
                is_min = False
                break
        if is_min:
            circuits.append(bm)
            if len(circuits) >= int(max_circuits):
                break

    # Convert to index lists
    out: list[list[int]] = []
    for bm in circuits:
        idxs = [i for i in range(M) if (bm >> i) & 1]
        out.append(idxs)
    return out


def augment_S_with_circuits(
    S_mod: np.ndarray,
    circuits: list[list[int]],
    *,
    incidence_label: int = 2,
) -> np.ndarray:
    """
    Build an augmented edge-colour matrix S_aug of size (M+C)x(M+C),
    where C = len(circuits). Term-term block is S_mod, and each circuit node is
    connected to its member terms with edge label 'incidence_label'. Circuit nodes
    are otherwise disconnected (0 labels).
    """
    S = np.asarray(S_mod, dtype=np.int64)
    M = S.shape[0]
    C = len(circuits)
    if C == 0:
        return S

    S_aug = np.zeros((M + C, M + C), dtype=np.int64)
    S_aug[:M, :M] = S

    for ci, members in enumerate(circuits):
        node = M + ci
        for i in members:
            S_aug[i, node] = incidence_label
            S_aug[node, i] = incidence_label
    return S_aug

__all__ = ['extract_circuits_from_nullspace_gf2', 'augment_S_with_circuits']
