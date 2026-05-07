"""Generate original-graph automorphism generators using the augmented graph."""

from __future__ import annotations

from time import perf_counter

from circuit_augmented_graph_builder import build_graphs, build_model


def permutation_to_list(permutation) -> list[int]:
    if hasattr(permutation, "mapping"):
        return list(permutation.mapping)
    return list(permutation)


def preserves_s_matrix(S_mod, vertex_colors, permutation) -> bool:
    n = S_mod.shape[0]
    permutation = list(map(int, permutation))

    if len(permutation) != n:
        return False
    if sorted(permutation) != list(range(n)):
        return False

    for i in range(n):
        if vertex_colors[i] != vertex_colors[permutation[i]]:
            return False

    for i in range(n):
        for j in range(n):
            if S_mod[i, j] != S_mod[permutation[i], permutation[j]]:
                return False

    return True


def main() -> None:
    pauli_sum = build_model()
    result = build_graphs(pauli_sum)

    base_graph = result["base_graph"]
    base_graph_colors = result["base_graph_colors"]
    augmented_graph = result["augmented_graph"]
    augmented_graph_colors = result["augmented_graph_colors"]
    S_mod = result["S_mod"]
    base_vertex_colors = result["base_vertex_colors"]
    n_original_vertices = S_mod.shape[0]

    print("=== Original graph automorphisms from augmented graph ===")
    print(f"original_vertices={n_original_vertices}")
    print(f"base_igraph_vertices={base_graph.vcount()} base_igraph_edges={base_graph.ecount()}")
    print(
        f"augmented_igraph_vertices={augmented_graph.vcount()} "
        f"augmented_igraph_edges={augmented_graph.ecount()}"
    )
    print(f"dependency_sets_found={len(result['dependency_sets'])}")

    base_count_start = perf_counter()
    base_aut_count = base_graph.count_automorphisms(color=base_graph_colors)
    base_count_seconds = perf_counter() - base_count_start
    print(f"original_igraph_automorphism_count={base_aut_count}")
    print(f"original_count_seconds={base_count_seconds:.6f}")

    base_group_start = perf_counter()
    base_aut_group = base_graph.automorphism_group(color=base_graph_colors)
    base_group_seconds = perf_counter() - base_group_start
    base_generators = (
        base_aut_group.generators
        if hasattr(base_aut_group, "generators")
        else base_aut_group
    )
    print(f"original_num_generators={len(base_generators)}")
    print(f"original_generator_seconds={base_group_seconds:.6f}")

    print("\nOriginal igraph generators:")
    base_valid_generators = []
    for idx, generator in enumerate(base_generators):
        full_permutation = permutation_to_list(generator)
        original_vertex_permutation = full_permutation[:n_original_vertices]
        is_valid = preserves_s_matrix(
            S_mod,
            base_vertex_colors,
            original_vertex_permutation,
        )
        base_valid_generators.append(is_valid)

        print(f"original_generator_{idx}_valid_on_original_graph={is_valid}")
        print(f"original_generator_{idx}_permutation={original_vertex_permutation}")

    print(f"all_original_generators_valid={all(base_valid_generators)}")

    augmented_count_start = perf_counter()
    augmented_aut_count = augmented_graph.count_automorphisms(color=augmented_graph_colors)
    augmented_count_seconds = perf_counter() - augmented_count_start
    print(f"\naugmented_graph_automorphism_count={augmented_aut_count}")
    print(f"augmented_count_seconds={augmented_count_seconds:.6f}")

    group_start = perf_counter()
    aut_group = augmented_graph.automorphism_group(color=augmented_graph_colors)
    group_seconds = perf_counter() - group_start
    generators = aut_group.generators if hasattr(aut_group, "generators") else aut_group
    print(f"augmented_num_generators={len(generators)}")
    print(f"augmented_generator_seconds={group_seconds:.6f}")

    print("\nAugmented igraph generators projected to original vertices:")
    valid_generators = []
    for idx, generator in enumerate(generators):
        full_permutation = permutation_to_list(generator)
        original_vertex_permutation = full_permutation[:n_original_vertices]
        is_valid = preserves_s_matrix(
            S_mod,
            base_vertex_colors,
            original_vertex_permutation,
        )
        valid_generators.append(is_valid)

        print(f"generator_{idx}_valid_on_original_graph={is_valid}")
        print(f"generator_{idx}_original_vertex_permutation={original_vertex_permutation}")

    print(f"all_generators_valid_on_original_graph={all(valid_generators)}")


if __name__ == "__main__":
    main()
