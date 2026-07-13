from .atomic_decomposition import (
    atomic_block_decompose,
    decompose_or_raise,
)
import numpy as np
from .modular_helpers import mod_p, rank_mod
from sympleq.core.graphs.utils import qudit_coupling_graph


def block_decompose(F: np.ndarray, p: int, **kwargs):
    """
    General decomposition entry point.

    Returns ``(Sigma, B, info)``.  The returned ``info`` dictionary distinguishes
    the attained qudit cost (``Q_att``/``qudit_cost``) from the optimal cost
    (``Q_opt``), which is only populated when minimality is certified.
    """
    return atomic_block_decompose(F, p, **kwargs)


def block_decompose_certified(
    F: np.ndarray,
    p: int,
    *,
    require_minimal: bool = True,
    **kwargs,
):
    """
    Certified decomposition entry point.

    Returns ``(Sigma, B, info)`` and raises ``CertificationError`` unless the
    result is a certified atomic decomposition.  By default it also requires the
    qudit cost to be certified minimal; set ``require_minimal=False`` to accept
    certified atomic decompositions whose minimality is not certified.
    """
    min_block_size = kwargs.pop("min_block_size", None)
    Sigma, B, info = decompose_or_raise(F, p, require_minimal=require_minimal, **kwargs)
    if min_block_size is not None:
        info = dict(info)
        info["requested_min_block_size"] = int(min_block_size)
        info.setdefault("warnings", []).append(
            "min_block_size is a legacy compatibility option and is not used by the atomic decomposer."
        )
    return Sigma, B, info


def block_decompose_optimal(F: np.ndarray, p: int, **kwargs):
    """
    Return a certified-minimal decomposition ``(Sigma, B, info)``.

    This function intentionally raises rather than returning an uncertified
    result.  It also returns the certificate information instead of discarding it,
    so downstream numerics cannot confuse an attained cost with a proven optimum.
    """
    if "require_minimal" in kwargs:
        raise TypeError("block_decompose_optimal always requires certified minimality; use block_decompose_certified for require_minimal=False.")
    min_block_size = kwargs.pop("min_block_size", None)
    Sigma, B, info = decompose_or_raise(F, p, require_minimal=True, **kwargs)
    if min_block_size is not None:
        info = dict(info)
        info["requested_min_block_size"] = int(min_block_size)
        info.setdefault("warnings", []).append(
            "min_block_size is a legacy compatibility option and is not used by the atomic decomposer."
        )
    return Sigma, B, info


def _mode_graph_from_S(S: np.ndarray, p: int) -> list[list[int]]:
    """
    Build adjacency for 'modes' (pairs). S is [U_all | V_all].
    Connect i--j if any entry in the 4x4 cross-block between modes i and j is nonzero mod p.
    """
    n2 = S.shape[0]
    assert n2 % 2 == 0
    k = n2 // 2  # number of modes
    adj = [[] for _ in range(k)]
    M = mod_p(S, p)
    for i in range(k):
        rows_i = [i, k + i]
        for j in range(i + 1, k):
            cols_j = [j, k + j]
            block_ij = M[np.ix_(rows_i, cols_j)]
            block_ji = M[np.ix_([j, k + j], [i, k + i])]
            if np.any(block_ij % p != 0) or np.any(block_ji % p != 0):
                adj[i].append(j)
                adj[j].append(i)
    return adj


def _components(adj: list[list[int]]) -> list[list[int]]:
    """Connected components from adjacency list."""
    k = len(adj)
    seen = [False] * k
    comps = []
    for s in range(k):
        if seen[s]:
            continue
        stack = [s]
        seen[s] = True
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True
                    stack.append(v)
        comps.append(sorted(comp))
    return comps


def ordered_block_sizes(S: np.ndarray, p: int) -> list[int]:
    """Return 2*n for connected components, ordered with the same nontrivial-first policy."""
    n2 = S.shape[0]
    k = n2 // 2
    adj = _mode_graph_from_S(S, p)
    comps = _components(adj)
    M = mod_p(S - np.eye(n2, dtype=np.int64), p)

    def comp_rank(comp):
        idx = comp + [c + k for c in comp]
        return rank_mod(M[np.ix_(idx, idx)], p)

    info = [{"modes": c, "n_modes": len(c), "rank": comp_rank(c)} for c in comps]
    non_triv_ge2 = [d for d in info if d["rank"] > 0 and d["n_modes"] >= 2]
    non_triv_1 = [d for d in info if d["rank"] > 0 and d["n_modes"] == 1]
    triv = [d for d in info if d["rank"] == 0]

    non_triv_ge2.sort(key=lambda d: (d["n_modes"], -d["rank"]))
    non_triv_1.sort(key=lambda d: -d["rank"])
    triv.sort(key=lambda d: d["n_modes"])

    ordered = non_triv_ge2 + non_triv_1 + triv
    return [2 * d["n_modes"] for d in ordered]


def block_indexes(F: np.ndarray) -> list[list[int]]:
    """
    Compute qudit blocks as connected components of the coupling graph.

    Args:
        F: 2n x 2n symplectic matrix (integer numpy array).

    Returns:
        A list of blocks, where each block is a sorted list of qudit indices.
        For example, if qudits (0,1) and (2,3) form two decoupled clusters,
        returns [[0, 1], [2, 3]] (order of blocks is not guaranteed).
    """
    neighbors = qudit_coupling_graph(F)
    n = len(neighbors)

    visited = [False] * n
    blocks: list[list[int]] = []

    for start in range(n):
        if visited[start]:
            continue

        # BFS/DFS to collect a connected component
        stack = [start]
        visited[start] = True
        component: list[int] = []

        while stack:
            v = stack.pop()
            component.append(v)
            for w in neighbors[v]:
                if not visited[w]:
                    visited[w] = True
                    stack.append(w)

        component.sort()
        blocks.append(component)

    # Optional: sort blocks by their smallest index for reproducibility
    blocks.sort(key=lambda comp: comp[0])
    return blocks
