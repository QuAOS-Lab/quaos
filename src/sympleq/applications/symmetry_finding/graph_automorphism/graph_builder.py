from __future__ import annotations

import numpy as np

try:
    import igraph as ig
except ImportError as exc:  # pragma: no cover - exercised only without optional extra
    raise ImportError(
        "Graph automorphism symmetry finding requires the optional "
        "`symmetry-finding` extra: install sympleq[symmetry-finding]."
    ) from exc


def _compress_colors(colors: np.ndarray | list[int]) -> list[int]:
    values = np.asarray(colors, dtype=object)
    unique = sorted(set(values.tolist()), key=repr)
    color_map = {value: idx for idx, value in enumerate(unique)}
    return [color_map[value] for value in values.tolist()]


def build_subdivision_graph_from_s_mod(
    S_mod: np.ndarray,
    vertex_colors: np.ndarray | list[int],
) -> tuple[ig.Graph, list[int]]:
    """
    Build an igraph graph for a complete edge-coloured symplectic product matrix.

    For binary matrices we use the minority edge value to keep the graph sparse.
    For larger fields we subdivide non-background coloured edges into auxiliary
    vertices because igraph/bliss supports vertex colours, not edge colours.
    """
    S = np.asarray(S_mod, dtype=int)
    if S.ndim != 2 or S.shape[0] != S.shape[1]:
        raise ValueError("S_mod must be a square matrix.")

    n = S.shape[0]
    h_colors = _compress_colors(vertex_colors)

    if np.all((S == 0) | (S == 1)):
        upper_rows, upper_cols = np.triu_indices(n, k=1)
        upper_values = S[upper_rows, upper_cols]
        one_count = int(np.count_nonzero(upper_values))
        zero_count = int(upper_values.size - one_count)
        edge_value = int(one_count <= zero_count)
        selected = upper_values == edge_value
        edges = list(zip(upper_rows[selected].tolist(), upper_cols[selected].tolist()))
        graph = ig.Graph(n=n, edges=edges, directed=False)
        graph.vs["vertex_color"] = h_colors
        graph.vs["original_id"] = list(range(n))
        return graph, h_colors

    upper_rows, upper_cols = np.triu_indices(n, k=1)
    upper_values = S[upper_rows, upper_cols]
    values, counts = np.unique(upper_values, return_counts=True)
    background_value = int(values[np.argmax(counts)])
    selected = upper_values != background_value
    base_edges = list(zip(upper_rows[selected].tolist(), upper_cols[selected].tolist()))
    edge_colors = upper_values[selected].astype(int).tolist()

    edge_color_offset = max(h_colors, default=-1) + 1
    unique_edge_colors = sorted(set(edge_colors), key=repr)
    edge_color_map = {value: edge_color_offset + idx for idx, value in enumerate(unique_edge_colors)}

    graph = ig.Graph(n=n + len(edge_colors), directed=False)
    subdivided_edges = []
    for edge_idx, (u, v) in enumerate(base_edges):
        w = n + edge_idx
        subdivided_edges.append((u, w))
        subdivided_edges.append((w, v))
        h_colors.append(edge_color_map[edge_colors[edge_idx]])

    graph.add_edges(subdivided_edges)
    graph.vs["vertex_color"] = h_colors
    graph.vs["original_id"] = list(range(n)) + [None] * len(edge_colors)
    return graph, h_colors
