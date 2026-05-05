from sympleq.models import ToricCode, ising_2d_hamiltonian
from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.core.graphs.graph_automorphism import clifford_from_pauli_permutation
from igraph_aut_gen import analyze_colored_graph, draw_original_graph, draw_vertex_colored_graph


def print_graph_before_automorphism(pauli_sum, color_mode='wl', max_wl_rounds=0):
    # Match the preprocessing used inside clifford_graph_automorphism_search
    pauli = pauli_sum.copy()
    # print(pauli.shape)
    pauli.weight_to_phase()

    S_mod = pauli.symplectic_product_matrix()
    p = int(pauli.lcm)
    coeffs = pauli.weights

    base_colors, base_classes = _build_base_partition(
        S_mod,
        p,
        coeffs=coeffs,
        col_invariants=None,
        max_rounds=max_wl_rounds,
        color_mode=color_mode,
    )

    print("Graph adjacency / edge-colour matrix S_mod")
    print(S_mod)

    print("Vertex colours")
    for i, color in enumerate(base_colors):
        print(f"vertex {i}: weight={coeffs[i]}, color={color}")

    print("Colour classes")
    for color, vertices in base_classes.items():
        print(f"color {color}: {vertices}")

    print("Non-zero graph edges")
    n = S_mod.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            if S_mod[i, j] != 0:
                print(f"{i} -- {j}, edge_colour={S_mod[i, j]}")

    return S_mod, base_colors, base_classes


def dense_edge_colored_graph_from_s_mod(S_mod, vertex_colors):
    edges = []
    edge_colors = []

    n = S_mod.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            edges.append((i, j))
            edge_colors.append(int(S_mod[i, j]))

    return n, edges, list(map(int, vertex_colors)), edge_colors


def is_valid_s_mod_automorphism(S_mod, vertex_colors, permutation):
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


def permutation_to_list(permutation):
    if hasattr(permutation, "mapping"):
        return list(permutation.mapping)
    return list(permutation)


def check_igraph_generators_on_s_mod(h, S_mod, vertex_colors, group=None):
    if group is None:
        group = h.automorphism_group(color=h.vs["vertex_color"])
    generators = group.generators if hasattr(group, "generators") else group
    valid = []

    for generator in generators:
        generator = permutation_to_list(generator)
        original_vertex_permutation = generator[:S_mod.shape[0]]
        valid.append(
            is_valid_s_mod_automorphism(
                S_mod,
                vertex_colors,
                original_vertex_permutation,
            )
        )

    print("\nValid igraph generator checks on original S_mod graph:")
    print(valid)
    print("All generators valid:", all(valid))

    return all(valid)


def build_cliffords_from_igraph_generators(pauli_sum, h, S_mod, vertex_colors, group=None):
    if group is None:
        group = h.automorphism_group(color=h.vs["vertex_color"])
    generators = group.generators if hasattr(group, "generators") else group

    symmetries = []
    for idx, generator in enumerate(generators):
        generator = permutation_to_list(generator)
        pauli_term_permutation = generator[:S_mod.shape[0]]
        if idx in (0, 1):
            print(f"igraph automorphism generator {idx}:")
            print(pauli_term_permutation)

        if not is_valid_s_mod_automorphism(S_mod, vertex_colors, pauli_term_permutation):
            print(f"generator {idx}: skipped, not valid on original S_mod graph")
            continue

        symmetry, reason = clifford_from_pauli_permutation(
            pauli_sum,
            pauli_term_permutation,
            return_reason=True,
        )
        if symmetry is None:
            print(f"generator {idx}: valid graph automorphism, failed Clifford leaf check: {reason}")
            continue

        print(f"generator {idx}: Clifford symmetry found")
        print(symmetry)
        symmetries.append(symmetry)

    print(f"\nClifford symmetries from igraph generators: {len(symmetries)}")
    return symmetries


def analyze_pauli_sum_with_igraph(pauli_sum, label, output_prefix):
    print(f"\n=== {label} ===")
    print(f"{label} PauliSum tableau shape:", pauli_sum.tableau.shape)

    S_mod, colors, classes = print_graph_before_automorphism(pauli_sum)

    print(f"\nigraph automorphism result for {label} graph")
    n, edges, vertex_colors, edge_colors = dense_edge_colored_graph_from_s_mod(S_mod, colors)
    g, h, aut_count, gens = analyze_colored_graph(
        n=n,
        edges=edges,
        vertex_colors=vertex_colors,
        edge_colors=edge_colors,
    )
    check_igraph_generators_on_s_mod(h, S_mod, vertex_colors, gens)
    igraph_cliffords = build_cliffords_from_igraph_generators(pauli_sum, h, S_mod, vertex_colors, gens)

    print("\nSymmetries from igraph generators:")
    for idx, symmetry in enumerate(igraph_cliffords):
        print(f"symmetry {idx}:")
        print("Stored symplectic F:")
        print(symmetry.symplectic)
        print("Applied tableau map F.T:")
        print(symmetry.symplectic.T)
        print("Phase vector h:")
        print(symmetry.phase_vector(pauli_sum.lcm))

    original_filename = f"{output_prefix}_original_edge_colored_graph.png"
    converted_filename = f"{output_prefix}_converted_vertex_colored_graph.png"
    draw_original_graph(
        g,
        vertex_colors=vertex_colors,
        edge_colors=edge_colors,
        filename=original_filename,
    )
    draw_vertex_colored_graph(
        h,
        h_vertex_colors=h.vs["vertex_color"],
        filename=converted_filename,
    )

    print(f"\nSaved {label} graph drawings:")
    print(original_filename)
    print(converted_filename)

    return g, h, S_mod, colors, gens, igraph_cliffords


def compose_permutations(left, right):
    return tuple(left[i] for i in right)


def enumerate_group_from_generators(generators, n_vertices):
    generators = [tuple(permutation_to_list(generator)) for generator in generators]
    identity = tuple(range(n_vertices))
    seen = {identity}
    queue = [(identity, ())]

    while queue:
        current, word = queue.pop(0)
        for generator_idx, generator in enumerate(generators):
            for side, composed in (
                ("L", compose_permutations(generator, current)),
                ("R", compose_permutations(current, generator)),
            ):
                if composed not in seen:
                    seen.add(composed)
                    composed_word = word + ((generator_idx, side),)
                    queue.append((composed, composed_word))
                    yield composed, composed_word


def first_clifford_from_igraph_generated_group(pauli_sum, h, S_mod, vertex_colors, group):
    generators = group.generators if hasattr(group, "generators") else group
    identity_pauli_perm = tuple(range(pauli_sum.n_paulis()))

    print("\nSearching Clifford lifts inside the group generated by igraph generators")
    print("igraph generator Pauli-term permutations:")
    for generator_idx, generator in enumerate(generators):
        pauli_term_generator = tuple(permutation_to_list(generator)[:S_mod.shape[0]])
        print(f"generator {generator_idx}:")
        print(pauli_term_generator)

    for idx, (full_permutation, generator_word) in enumerate(enumerate_group_from_generators(generators, h.vcount())):
        pauli_term_permutation = tuple(full_permutation[:S_mod.shape[0]])
        if pauli_term_permutation == identity_pauli_perm:
            continue

        if not is_valid_s_mod_automorphism(S_mod, vertex_colors, pauli_term_permutation):
            continue

        symmetry, reason = clifford_from_pauli_permutation(
            pauli_sum,
            pauli_term_permutation,
            return_reason=True,
        )
        if symmetry is None:
            continue

        print("\nSuccessful generated-group permutation:")
        print("First Clifford symmetry found from igraph-generated group")
        print("Generated-group permutation index:", idx)
        print("Generator composition word:")
        print(generator_word)
        print("Pauli-term permutation:")
        print(pauli_term_permutation)
        print("Stored symplectic F:")
        print(symmetry.symplectic)
        print("Applied tableau map F.T:")
        print(symmetry.symplectic.T)
        print("Phase vector h:")
        print(symmetry.phase_vector(pauli_sum.lcm))
        return symmetry, pauli_term_permutation

    print("No non-identity Clifford symmetry found from the igraph-generated group.")
    return None, None


# Load an Ising ladder model: a 2 x L transverse-field Ising strip.
ising_ladder = ising_2d_hamiltonian(
    n_x=2,
    n_y=200,
    J_zz=1.0,
    h_x=0.7,
    periodic=False,
)
ising_ladder_cliffords = analyze_pauli_sum_with_igraph(
    ising_ladder,
    label="Ising ladder",
    output_prefix="ising_ladder",
)

# toric = ToricCode(
#     Nx=2,
#     Ny=2,
#     c_x=1.0,
#     c_z=2.0,
#     c_g=3.0,
#     periodic=True,
# )
# toric_hamiltonian = toric.hamiltonian()
# toric_g, toric_h, toric_S_mod, toric_colors, toric_gens, toric_cliffords = analyze_pauli_sum_with_igraph(
#     toric_hamiltonian,
#     label="Toric code",
#     output_prefix="toric",
# )
# if not toric_cliffords:
#     first_clifford_from_igraph_generated_group(
#         toric_hamiltonian,
#         toric_h,
#         toric_S_mod,
#         toric_colors,
#         toric_gens,
#     )

# # # Now run the automorphism search
# symmetries = find_clifford_symmetries(H, num_symmetries=1)

# print(f"Found {len(symmetries)} symmetry/symmetries")
# for sym in symmetries:
#     print(sym)
