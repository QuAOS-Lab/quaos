from __future__ import annotations

from typing import Any

import numpy as np


def compress_vertex_colors(coeffs: Any) -> list[int]:
    """Compress coefficient labels to contiguous integer vertex colors."""
    _, inverse = np.unique(np.asarray(coeffs), return_inverse=True)
    return inverse.astype(int).tolist()


def build_igraph_graph_from_s_mod(S_mod: np.ndarray, vertex_colors: Any) -> tuple[Any, list[int]]:
    """Build an igraph-compatible colored graph from a symplectic product matrix.

    The original Hamiltonian term vertices occupy indices ``0..n-1`` in the
    returned graph. For binary edge labels, this builds the sparser of the graph
    or its complement. For higher-cardinality edge labels, the most common edge
    value is treated as the implicit no-edge value, and the remaining edge values
    are represented by colored subdivision vertices.
    """
    import igraph as ig

    S = np.asarray(S_mod, dtype=int)
    n = int(S.shape[0])
    if S.shape != (n, n):
        raise ValueError(f"S_mod must be square, got shape {S.shape}.")

    vertex_colors_arr = np.asarray(vertex_colors)
    if vertex_colors_arr.shape[0] != n:
        raise ValueError(
            f"vertex_colors length must match S_mod size {n}, got {vertex_colors_arr.shape[0]}."
        )

    unique_vertex_colors = sorted(set(vertex_colors_arr.tolist()), key=str)
    vertex_color_map = {c: i for i, c in enumerate(unique_vertex_colors)}
    h_colors = [vertex_color_map[c] for c in vertex_colors_arr.tolist()]

    upper_rows, upper_cols = np.triu_indices(n, k=1)
    upper_values = S[upper_rows, upper_cols]

    if np.all((S == 0) | (S == 1)):
        one_count = int(np.count_nonzero(upper_values))
        zero_count = int(upper_values.size - one_count)
        edge_value = int(one_count <= zero_count)
        selected = upper_values == edge_value
        edges = list(zip(upper_rows[selected].tolist(), upper_cols[selected].tolist()))

        graph = ig.Graph(n=n, edges=edges, directed=False)
        graph.vs["vertex_color"] = h_colors
        graph.vs["original_id"] = list(range(n))
        graph.vs["edge_color"] = [None] * n
        graph["implicit_edge_value"] = 1 - edge_value
        return graph, h_colors

    values, counts = np.unique(upper_values.astype(int), return_counts=True)
    no_edge_value = int(values[np.argmax(counts)])
    selected = upper_values != no_edge_value
    edges = list(zip(upper_rows[selected].tolist(), upper_cols[selected].tolist()))
    edge_colors = upper_values[selected].astype(int).tolist()

    edge_color_offset = max(h_colors, default=-1) + 1
    unique_edge_colors = sorted(set(edge_colors), key=str)
    edge_color_map = {c: edge_color_offset + i for i, c in enumerate(unique_edge_colors)}

    graph = ig.Graph(n=n + len(edge_colors), directed=False)
    subdivided_edges = []
    for edge_idx, (u, v) in enumerate(edges):
        w = n + edge_idx
        subdivided_edges.append((u, w))
        subdivided_edges.append((w, v))
        h_colors.append(edge_color_map[edge_colors[edge_idx]])

    graph.add_edges(subdivided_edges)
    graph.vs["vertex_color"] = h_colors
    graph.vs["original_id"] = list(range(n)) + [None] * len(edges)
    graph.vs["edge_color"] = [None] * n + edge_colors
    graph["implicit_edge_value"] = no_edge_value
    return graph, h_colors


__all__ = ["compress_vertex_colors", "build_igraph_graph_from_s_mod"]
