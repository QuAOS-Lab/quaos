import numpy as np
from sympleq.core.circuits import Gate  # , Circuit
from sympleq.core.circuits.phase_correction import clifford_phase_decomposition
from sympleq.core.paulis import PauliSum
from sympleq.core.graphs.graph_automorphism_search import clifford_graph_automorphism_search
from sympleq.core.symmetries.block_decomposition import block_decompose_optimal, ordered_block_sizes


def min_qudit_clifford_symmetry(pauli_sum: PauliSum, talk: bool = False, progress: bool = False,
                                heuristic: bool = False, circuit_augmented_graph: bool | str = False,
                                max_nullity_for_circuits: int = 12, max_circuits: int = 5000,
                                ) -> tuple[Gate, Gate, Gate]:
    """
    Find a single Clifford symmetry g of the given PauliSum, and decompose it into blocks with a minimal qudit cost
    via a symplectic similarity transform: g = T S T^{-1}.

    :param pauli_sum: Input Hamiltonian
    :type pauli_sum: PauliSum
    :return: Description
    :rtype: tuple[Gate, Gate, Gate]
    """

    if heuristic:
        extra_column_invariants = "hist"
    else:
        extra_column_invariants = "none"

    G = find_clifford_symmetries(pauli_sum, num_symmetries=1,
                                 dynamic_refine_every=0, progress=progress,
                                 extra_column_invariants=extra_column_invariants,
                                 circuit_augmented_graph=circuit_augmented_graph,
                                 max_nullity_for_circuits=max_nullity_for_circuits,
                                 max_circuits=max_circuits)
    if len(G) == 0:
        # save pauli_sum to file for debugging, tableau, weights, phases
        raise RuntimeError("No non-trivial Clifford symmetry found for the given PauliSum.")
    g = G[0]

    if talk:
        print('Got symmetry - decomposing')

    S, T, _info = block_decompose_optimal(g.symplectic, int(pauli_sum.lcm), min_block_size=4)
    h_S, h_T = clifford_phase_decomposition(g.symplectic, g.phase_vector(), S, T, int(pauli_sum.lcm))
    S_gate = Gate('S', S, h_S)
    T_gate = Gate('T', T, h_T)

    return g, S_gate, T_gate


def multiple_min_qudit_clifford_symmetries(pauli_sum: PauliSum,
                                           n_symmetries: int = 1,
                                           ) -> tuple[list[Gate], list[Gate], list[Gate]]:
    """
    Find multiple Clifford symmetries of the given PauliSum, and decompose each into blocks via a symplectic similarity

    """
    G = find_clifford_symmetries(pauli_sum, num_symmetries=n_symmetries,
                                 dynamic_refine_every=0)

    Ss = []
    Ts = []
    for i, g in enumerate(G):
        S, T, _info = block_decompose_optimal(g.symplectic, pauli_sum.lcm)
        h_S, h_T = clifford_phase_decomposition(g.symplectic, g.phase_vector(), S, T, int(pauli_sum.lcm))
        S_gate = Gate(f'S{i}', S, h_S)
        T_gate = Gate(f'T{i}', T, h_T)
        Ss.append(S_gate)
        Ts.append(T_gate)

    return G, Ss, Ts


def block_structure(gate: Gate, lcm: int):
    symp = gate.symplectic
    sizes = np.asarray(ordered_block_sizes(symp, lcm), dtype=int) / 2
    return sizes


def qudit_cost(gate: Gate, lcm: int):
    return int(max(block_structure(gate, lcm)))


def find_clifford_symmetries(
    pauli_sum: PauliSum,
    num_symmetries: int = 1,
    # Strategy
    dynamic_refine_every: int = 0,
    use_code_induced_completion: bool = True,
    extra_column_invariants: str = "none",
    p2_bitset: str = "auto",
    color_mode: str = "wl",
    max_wl_rounds: int = 10,
    circuit_augmented_graph: bool | str = False,
    max_nullity_for_circuits: int = 12,
    max_circuits: int = 5000,
    shuffle_domain_order: bool = False,
    progress: bool = False,
    progress_every: int = 2048,
) -> list[Gate]:
    """
    Return up to k automorphisms preserving S and the vector set. See flags above.
    """
    return clifford_graph_automorphism_search(
        pauli_sum,
        k_wanted=num_symmetries,
        use_code_induced_completion=bool(use_code_induced_completion),
        extra_column_invariants=extra_column_invariants,
        p2_bitset=p2_bitset,
        color_mode=color_mode,
        max_wl_rounds=max_wl_rounds,
        circuit_augmented_graph=circuit_augmented_graph,
        max_nullity_for_circuits=int(max_nullity_for_circuits),
        max_circuits=int(max_circuits),
        shuffle_domain_order=bool(shuffle_domain_order),
        dynamic_refine_every=int(dynamic_refine_every),
        progress=progress,
        progress_every=int(progress_every),
    )
