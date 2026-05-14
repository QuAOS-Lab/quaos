from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import galois


@dataclass(frozen=True, slots=True)
class FundamentalRelation:
    """A normalized fundamental relation r_j in F_p^M.

    The relation is stored sparsely as (index, coefficient) pairs.  The
    non-basis index has coefficient +1; basis indices have coefficients
    -lambda.  For p=2 all nonzero coefficients are 1, so this reduces to an
    ordinary circuit support.
    """

    dependent_index: int
    indices: tuple[int, ...]
    coefficients: tuple[int, ...]

    @property
    def support(self) -> tuple[int, ...]:
        return self.indices

    @property
    def size(self) -> int:
        return len(self.indices)


def _gf2_inv(A: np.ndarray) -> np.ndarray:
    """Invert a square binary matrix over GF(2)."""
    A = (np.asarray(A, dtype=np.uint8) & 1).copy()
    n, m = A.shape
    if n != m:
        raise np.linalg.LinAlgError("matrix must be square")

    aug = np.concatenate([A, np.eye(n, dtype=np.uint8)], axis=1)
    row = 0
    for col in range(n):
        piv = -1
        for r in range(row, n):
            if aug[r, col] & 1:
                piv = r
                break
        if piv < 0:
            raise np.linalg.LinAlgError("matrix is singular over GF(2)")
        if piv != row:
            aug[[row, piv]] = aug[[piv, row]]
        for r in range(n):
            if r != row and (aug[r, col] & 1):
                aug[r, :] ^= aug[row, :]
        row += 1
    return aug[:, n:] & 1


def extract_fundamental_relations_from_matroid(
    G: np.ndarray | galois.FieldArray,
    basis_cols: np.ndarray,
    *,
    p: int,
) -> list[FundamentalRelation]:
    """Construct the fundamental relation generators from a matroid matrix.

    Parameters
    ----------
    G:
        n_b x M coordinate/generator matrix whose columns are Pauli labels in
        some coordinate basis.
    basis_cols:
        Column indices that form the chosen basis.  They need not be the first
        n_b columns.
    p:
        Prime field size.

    Returns
    -------
    list[FundamentalRelation]
        One normalized relation for each non-basis column j:

            e_j - sum_t X[t, j] e_{basis_cols[t]},

        where X = G[:, basis_cols]^{-1} G, so the chosen basis columns become
        the identity.
    """
    p = int(p)
    basis_cols = np.asarray(basis_cols, dtype=np.int64).reshape(-1)
    n_b = int(basis_cols.size)

    G_int = np.asarray(G, dtype=int) % p
    if G_int.ndim != 2:
        raise ValueError("G must be a 2D matrix")
    if G_int.shape[0] != n_b:
        raise ValueError("number of basis columns must equal rank / number of rows of G")

    M = int(G_int.shape[1])
    basis_mask = np.zeros(M, dtype=bool)
    basis_mask[basis_cols] = True
    nonbasis_cols = [j for j in range(M) if not basis_mask[j]]

    C = G_int[:, basis_cols] % p
    if p == 2:
        C_inv = _gf2_inv(C)
        X = (C_inv @ (G_int & 1)) & 1
    else:
        GF = galois.GF(p)
        C_inv = np.linalg.inv(GF(C))
        X = np.asarray(C_inv @ GF(G_int), dtype=int) % p

    relations: list[FundamentalRelation] = []
    for j in nonbasis_cols:
        idxs: list[int] = []
        coeffs: list[int] = []

        # Basis part: coefficient is -lambda.
        for t, bcol in enumerate(basis_cols.tolist()):
            lam = int(X[t, j]) % p
            if lam:
                idxs.append(int(bcol))
                coeffs.append((-lam) % p)

        # Non-basis/dependent part: normalized coefficient +1.
        idxs.append(int(j))
        coeffs.append(1 % p)

        # Keep a deterministic order by Pauli index.  Coefficients follow the
        # same permutation.
        order = np.argsort(np.asarray(idxs, dtype=np.int64), kind="stable")
        idxs_t = tuple(int(idxs[k]) for k in order.tolist())
        coeffs_t = tuple(int(coeffs[k]) for k in order.tolist())
        relations.append(
            FundamentalRelation(
                dependent_index=int(j),
                indices=idxs_t,
                coefficients=coeffs_t,
            )
        )

    return relations


def relation_edge_label(coefficient: int, *, p: int) -> int:
    """Map a nonzero relation coefficient to an edge colour disjoint from S_mod.

    Term-term symplectic colours occupy 0, ..., p-1.  Relation-incidence
    colours are encoded as p + coefficient, with coefficient in {1, ..., p-1}.
    """
    c = int(coefficient) % int(p)
    if c == 0:
        raise ValueError("zero coefficients should not be encoded as incidence edges")
    return int(p) + c


def augment_S_with_fundamental_relations(
    S_mod: np.ndarray,
    relations: list[FundamentalRelation],
    *,
    p: int,
    coefficient_labels: bool = True,
) -> tuple[np.ndarray, int]:
    """Build an augmented edge-colour matrix using fundamental relation nodes.

    Returns
    -------
    S_aug, p_for_wl:
        S_aug has one auxiliary node for each fundamental relation.  Edge
        colours between term nodes are unchanged.  Incidence edges from a
        relation node to its Pauli-label members encode the relation
        coefficient when coefficient_labels=True.  The returned p_for_wl is one
        plus the largest edge colour and can be passed to _build_base_partition.
    """
    S = np.asarray(S_mod, dtype=np.int64)
    M = int(S.shape[0])
    R = len(relations)
    if R == 0:
        return S, max(int(np.max(S)) + 1 if S.size else 1, int(p))

    S_aug = np.zeros((M + R, M + R), dtype=np.int64)
    S_aug[:M, :M] = S

    for ri, rel in enumerate(relations):
        node = M + ri
        for idx, coeff in zip(rel.indices, rel.coefficients):
            if coefficient_labels:
                label = relation_edge_label(coeff, p=p)
            else:
                # Useful for binary support-only experiments.  Keep disjoint
                # from term-term colours.
                label = int(p) + 1
            S_aug[int(idx), node] = label
            S_aug[node, int(idx)] = label

    p_for_wl = int(np.max(S_aug)) + 1
    return S_aug, p_for_wl


# ---------------------------------------------------------------------------
# Backwards-compatible helpers for the old all-circuits GF(2) augmentation.
# These are retained so existing imports/tests do not break, but the preferred
# helper for the paper's construction is extract_fundamental_relations_from_matroid.
# ---------------------------------------------------------------------------


def _gf2_nullspace_basis(A: np.ndarray) -> list[np.ndarray]:
    A = (np.asarray(A, dtype=np.uint8) & 1).copy()
    m, n = A.shape
    row = 0
    pivots: list[int] = []
    pivot_rows: list[int] = [-1] * n

    for col in range(n):
        if row >= m:
            break
        piv = -1
        for r in range(row, m):
            if A[r, col] & 1:
                piv = r
                break
        if piv < 0:
            continue
        if piv != row:
            A[[row, piv]] = A[[piv, row]]
        pivots.append(col)
        pivot_rows[col] = row
        for r in range(m):
            if r != row and (A[r, col] & 1):
                A[r, :] ^= A[row, :]
        row += 1

    pivot_set = set(pivots)
    free_cols = [c for c in range(n) if c not in pivot_set]
    basis: list[np.ndarray] = []
    for f in free_cols:
        x = np.zeros(n, dtype=np.uint8)
        x[f] = 1
        for pcol in pivots:
            r = pivot_rows[pcol]
            if r >= 0 and (A[r, f] & 1):
                x[pcol] = 1
        basis.append(x)
    return basis


def _gf2_vec_to_bitmask(x: np.ndarray) -> int:
    x = (np.asarray(x, dtype=np.uint8) & 1).reshape(-1)
    bm = 0
    for i, b in enumerate(x.tolist()):
        if b:
            bm |= 1 << i
    return bm


def _bitcount(bm: int) -> int:
    return int(bm.bit_count())


def extract_circuits_from_nullspace_gf2(
    P: np.ndarray,
    *,
    max_nullity: int = 12,
    max_circuits: int = 5000,
) -> list[list[int]]:
    P2 = (np.asarray(P, dtype=np.uint8) & 1)
    M, _ = P2.shape
    basis = _gf2_nullspace_basis(P2.T)
    r = len(basis)
    if r == 0 or r > int(max_nullity):
        return []

    basis_bm = [_gf2_vec_to_bitmask(v) for v in basis]
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
        if bm:
            deps.add(bm)

    circuits: list[int] = []
    for bm in sorted(deps, key=_bitcount):
        if all((c & bm) != c for c in circuits):
            circuits.append(bm)
            if len(circuits) >= int(max_circuits):
                break

    return [[i for i in range(M) if (bm >> i) & 1] for bm in circuits]


def augment_S_with_circuits(
    S_mod: np.ndarray,
    circuits: list[list[int]],
    *,
    incidence_label: int = 2,
) -> np.ndarray:
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


__all__ = [
    "FundamentalRelation",
    "extract_fundamental_relations_from_matroid",
    "augment_S_with_fundamental_relations",
    "relation_edge_label",
    "extract_circuits_from_nullspace_gf2",
    "augment_S_with_circuits",
]
