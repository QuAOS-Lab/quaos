import igraph as ig
from collections import Counter
from time import perf_counter


# ============================================================
# USER SETTINGS
# ============================================================

n = 5

edges = [
    (0, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (4, 0),
]

vertex_colors = [0, 0, 1, 1, 2]

edge_colors = [0, 1, 0, 1, 0]

# ============================================================


def canonical_edge(u, v):
    return (u, v) if u < v else (v, u)


def validate_inputs(n, edges, vertex_colors, edge_colors):
    if n <= 0:
        raise ValueError("n must be positive.")

    if len(vertex_colors) != n:
        raise ValueError("vertex_colors must have length n.")

    if len(edge_colors) != len(edges):
        raise ValueError("edge_colors must have the same length as edges.")

    seen_edges = set()

    for e in edges:
        if len(e) != 2:
            raise ValueError(f"Invalid edge {e}. Each edge must have two vertices.")

        u, v = e

        if not (0 <= u < n and 0 <= v < n):
            raise ValueError(f"Edge {e} contains a vertex outside 0,...,{n-1}.")

        if u == v:
            raise ValueError("This script assumes no self-loops.")

        ce = canonical_edge(u, v)

        if ce in seen_edges:
            raise ValueError(f"Duplicate undirected edge found: {ce}")

        seen_edges.add(ce)


def build_general_colored_graph(n, edges, vertex_colors, edge_colors):
    validate_inputs(n, edges, vertex_colors, edge_colors)

    g = ig.Graph(n=n, edges=edges, directed=False)
    g.vs["vertex_color"] = vertex_colors
    g.es["edge_color"] = edge_colors

    return g


def edge_colored_to_vertex_colored(g, edge_colors, vertex_colors):
    if g.is_directed():
        raise ValueError("This conversion is for undirected graphs only.")

    n = g.vcount()
    m = g.ecount()

    unique_vertex_colors = sorted(set(vertex_colors), key=str)
    vertex_color_map = {c: i for i, c in enumerate(unique_vertex_colors)}

    h_vertex_colors = [vertex_color_map[c] for c in vertex_colors]

    edge_color_offset = max(h_vertex_colors) + 1

    unique_edge_colors = sorted(set(edge_colors), key=str)
    edge_color_map = {
        c: edge_color_offset + i
        for i, c in enumerate(unique_edge_colors)
    }

    h = ig.Graph(n=n + m, directed=False)

    new_edges = []
    edge_vertex_map = {}

    for e_idx, e in enumerate(g.es):
        u, v = e.tuple
        w = n + e_idx

        edge_vertex_map[e_idx] = w

        new_edges.append((u, w))
        new_edges.append((w, v))

        h_vertex_colors.append(edge_color_map[edge_colors[e_idx]])

    h.add_edges(new_edges)
    h.vs["vertex_color"] = h_vertex_colors

    return h, h_vertex_colors, edge_vertex_map


def map_labels_to_palette(labels):
    palette = [
        "red", "blue", "green", "orange",
        "purple", "cyan", "pink", "yellow",
        "brown", "gray", "magenta", "lime",
    ]

    unique_labels = sorted(set(labels), key=str)
    label_to_color = {
        label: palette[i % len(palette)]
        for i, label in enumerate(unique_labels)
    }

    return [label_to_color[label] for label in labels]


def draw_original_graph(g, vertex_colors, edge_colors, filename="original_graph.png"):
    plot_vertex_colors = map_labels_to_palette(vertex_colors)
    plot_edge_colors = map_labels_to_palette(edge_colors)

    layout = g.layout("fruchterman_reingold")

    ig.plot(
        g,
        target=filename,
        layout=layout,
        vertex_label=list(range(g.vcount())),
        vertex_color=plot_vertex_colors,
        edge_color=plot_edge_colors,
        vertex_size=35,
        edge_width=3,
        bbox=(700, 700),
        margin=70,
    )


def draw_vertex_colored_graph(h, h_vertex_colors, filename="converted_graph.png"):
    plot_vertex_colors = map_labels_to_palette(h_vertex_colors)

    layout = h.layout("fruchterman_reingold")

    ig.plot(
        h,
        target=filename,
        layout=layout,
        vertex_label=list(range(h.vcount())),
        vertex_color=plot_vertex_colors,
        vertex_size=30,
        edge_width=1,
        bbox=(800, 800),
        margin=80,
    )


def analyze_colored_graph(n, edges, vertex_colors, edge_colors):
    g = build_general_colored_graph(
        n=n,
        edges=edges,
        vertex_colors=vertex_colors,
        edge_colors=edge_colors,
    )

    print("\nOriginal graph")
    print("Number of vertices:", g.vcount())
    print("Number of edges:", g.ecount())

    print("\nVertex color counts:")
    print(Counter(vertex_colors))

    print("\nEdge color counts:")
    print(Counter(edge_colors))

    h, h_vertex_colors, edge_vertex_map = edge_colored_to_vertex_colored(
        g,
        edge_colors=edge_colors,
        vertex_colors=vertex_colors,
    )

    print("\nConverted vertex-colored subdivision graph")
    print("Number of vertices:", h.vcount())
    print("Number of edges:", h.ecount())

    print("\nConverted vertex color counts:")
    print(Counter(h_vertex_colors))

    print("\nColor-preserving automorphism count:")
    aut_count = h.count_automorphisms(color=h_vertex_colors)
    print(aut_count)

    print("\nGenerators of the automorphism group:")
    start_time = perf_counter()
    gens = h.automorphism_group(color=h_vertex_colors)
    elapsed_time = perf_counter() - start_time
    print(len(gens))
    print(f"Automorphism generator time: {elapsed_time:.6f} seconds")

    return g, h, aut_count, gens


if __name__ == "__main__":
    g, h, aut_count, gens = analyze_colored_graph(
        n=n,
        edges=edges,
        vertex_colors=vertex_colors,
        edge_colors=edge_colors,
    )
