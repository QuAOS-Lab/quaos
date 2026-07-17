import multiprocessing as mp
from pathlib import Path
from queue import Empty
from time import perf_counter

try:
    import igraph as ig
except ImportError as exc:  # pragma: no cover - exercised only without optional extra
    raise ImportError(
        "Graph automorphism symmetry finding requires the optional "
        "`symmetry-finding` extra: install sympleq[symmetry-finding]."
    ) from exc
import matplotlib
import numpy as np

from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.models.Ising import ising_2d_hamiltonian
from sympleq.models.heisenberg import all_to_all_heisenberg_hamiltonian
from sympleq.models.toric_code import ToricCode


def compress_vertex_colors(coeffs):
    # coeffs: coefficients of the Pauli Strings in the Hamiltonian;
    # returns a list of integers that can be used as vertex colors.
    _, inverse = np.unique(np.asarray(coeffs), return_inverse=True)
    return inverse.astype(int).tolist()


def build_simple_graph_from_s_mod(S_mod, vertex_colors):
    # S_mod : Symplectic product matrix corr. to the Hamiltonian;
    # Vertex colors: Coefficients of the Pauli Strings in the Hamiltonian;
    # returns the simple igraph-graph corr to the Hamiltonians;
    # Works only for qubits.
    assert np.all((S_mod == 0) | (S_mod == 1))
    n_vertices = S_mod.shape[0]
    edges = []
    for i in range(n_vertices):
        for j in range(i + 1, n_vertices):
            if int(S_mod[i, j]) == 1:
                edges.append((i, j))

    graph = ig.Graph(n=n_vertices, edges=edges, directed=False)
    graph.vs["vertex_color"] = vertex_colors
    graph.vs["original_id"] = list(range(n_vertices))
    return graph, vertex_colors


# For QUDITS

def permutation_to_tuple(permutation):
    if hasattr(permutation, "mapping"):
        return tuple(permutation.mapping)
    return tuple(permutation)

def build_subdivision_graph_from_s_mod(S_mod, vertex_colors):
    # S_mod[i, j] == 0 means no edge.
    # S_mod[i, j] != 0 means an edge with color S_mod[i, j].
    # Since igraph/bliss consumes vertex colors, each colored edge is
    # represented by a subdivision vertex whose color is the edge color.
    # Works for qubits and qubits
    n = S_mod.shape[0]
    edges = []
    edge_colors = []

    if np.all((S_mod == 0) | (S_mod == 1)):
        upper_rows, upper_cols = np.triu_indices(n, k=1)
        upper_values = S_mod[upper_rows, upper_cols]
        one_count = int(np.count_nonzero(upper_values))
        zero_count = int(upper_values.size - one_count)
        edge_value = int(one_count <= zero_count)
        print(f"qubit S_mod counts: 0={zero_count}, 1={one_count}, edge_value={edge_value}")
        selected = upper_values == edge_value
        edges = list(zip(upper_rows[selected].tolist(), upper_cols[selected].tolist()))
        unique_vertex_colors = sorted(set(vertex_colors), key=str)
        vertex_color_map = {c: i for i, c in enumerate(unique_vertex_colors)}
        h_colors = [vertex_color_map[int(c)] for c in vertex_colors]
        h = ig.Graph(n=n, edges=edges, directed=False)
        h.vs["vertex_color"] = h_colors
        h.vs["original_id"] = list(range(n))
        h.vs["edge_color"] = [None] * n
        return h, h_colors

    upper_rows, upper_cols = np.triu_indices(n, k=1)
    upper_values = S_mod[upper_rows, upper_cols].astype(int)
    values, counts = np.unique(upper_values, return_counts=True)
    no_edge_value = int(values[np.argmax(counts)])
    print(
        f"non-qubit S_mod counts: {dict(zip(values.tolist(), counts.tolist()))}, "
        f"no_edge_value={no_edge_value}"
    )
    selected = upper_values != no_edge_value
    edges = list(zip(upper_rows[selected].tolist(), upper_cols[selected].tolist()))
    edge_colors = upper_values[selected].tolist()

    unique_vertex_colors = sorted(set(vertex_colors), key=str)
    vertex_color_map = {c: i for i, c in enumerate(unique_vertex_colors)}
    h_colors = [vertex_color_map[int(c)] for c in vertex_colors]

    edge_color_offset = max(h_colors, default=-1) + 1
    unique_edge_colors = sorted(set(edge_colors), key=str)
    edge_color_map = {c: edge_color_offset + i for i, c in enumerate(unique_edge_colors)}
    print(edge_color_map)

    h = ig.Graph(n=n + len(edge_colors), directed=False)
    subdivided_edges = []
    for edge_idx, (u, v) in enumerate(edges):
        w = n + edge_idx
        subdivided_edges.append((u, w))
        subdivided_edges.append((w, v))
        h_colors.append(edge_color_map[edge_colors[edge_idx]])
    h.add_edges(subdivided_edges)
    h.vs["vertex_color"] = h_colors
    h.vs["original_id"] = list(range(n)) + [None] * len(edges)
    h.vs["edge_color"] = [None] * n + edge_colors
    return h, h_colors


############### Visualization ###############
def print_input_and_bliss_optimized_nodes(graph, colors):
    original_vertices = sum(original_id is not None for original_id in graph.vs["original_id"])
    subdivision_vertices = graph.vcount() - original_vertices
    print(f"Original Hamiltonian vertices: {original_vertices}")
    print(f"Edge-color subdivision vertices: {subdivision_vertices}")
    print(f"Converted graph vertices: {graph.vcount()}")
    try:
        canonical_perm = permutation_to_tuple(graph.canonical_permutation(color=colors))
        optimized_graph = graph.permute_vertices(canonical_perm)
        print(f"Bliss-canonical graph vertices: {optimized_graph.vcount()}")
        final_vertex_mapping = {}
        for new_idx, old_idx in enumerate(canonical_perm):
            edge_color = graph.vs[old_idx]["edge_color"]
            original_id = graph.vs[old_idx]["original_id"]
            if edge_color is None:
                label = f"term:{original_id}"
            else:
                endpoints = tuple(
                    sorted(
                        graph.vs[neighbor]["original_id"]
                        for neighbor in graph.neighbors(old_idx)
                    )
                )
                label = f"edge:{endpoints}:color:{edge_color}"
            final_vertex_mapping[label] = new_idx
        print(f"Final vertex mapping: {final_vertex_mapping}")
    except Exception as exc:
        print(f"Could not build Bliss-optimized graph view: {type(exc).__name__}: {exc}")


def graph_filename(model_name, nx, ny, suffix):
    return SCRIPT_DIR / f"igraph_{model_name}_{nx}_{ny}_graph.{suffix}"


def save_graph_view(graph, filename):
    layout = np.asarray(graph.layout("fr").coords, dtype=float)
    colors = graph.vs["vertex_color"]
    unique_colors = sorted(set(colors))
    color_map = {
        color: plt.cm.tab20(idx % 20)
        for idx, color in enumerate(unique_colors)
    }

    fig, ax = plt.subplots(figsize=(8, 8))
    for u, v in graph.get_edgelist():
        ax.plot(
            [layout[u, 0], layout[v, 0]],
            [layout[u, 1], layout[v, 1]],
            color="0.75",
            linewidth=0.8,
            zorder=1,
        )

    ax.scatter(
        layout[:, 0],
        layout[:, 1],
        c=[color_map[color] for color in colors],
        s=45,
        edgecolors="black",
        linewidths=0.5,
        zorder=2,
    )
    for idx, (x_pos, y_pos) in enumerate(layout):
        ax.text(x_pos, y_pos, str(idx), ha="center", va="center", fontsize=6, zorder=3)

    ax.set_title(
        f"Simple S_mod Graph: {graph.vcount()} vertices, {graph.ecount()} edges"
    )
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(filename, dpi=180)
    plt.close(fig)
    print(f"Saved graph PNG: {filename}")


def save_graph_files(graph, model_name, nx, ny):
    # graphml_filename = graph_filename(model_name, nx, ny, "graphml")
    png_filename = graph_filename(model_name, nx, ny, "png")
    save_graph_view(graph, png_filename)
