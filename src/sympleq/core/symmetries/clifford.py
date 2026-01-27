import numpy as np
from sympleq.core.circuits import Gate, Circuit
from sympleq.core.paulis import PauliSum
from sympleq.core.graphs.graph_automorphism import clifford_graph_automorphism_search
from sympleq.core.symmetries.phase_correction import clifford_phase_decomposition
from sympleq.core.symmetries.block_decomposition import block_decompose, ordered_block_sizes


def min_qudit_clifford_symmetry(pauli_sum: PauliSum,
                                ) -> tuple[Gate, Gate, Gate]:
    """
    Find a single Clifford symmetry g of the given PauliSum, and decompose it into blocks with a minimal qudit cost
    via a symplectic similarity transform: g = T S T^{-1}.

    :param pauli_sum: Input Hamiltonian
    :type pauli_sum: PauliSum
    :return: Description
    :rtype: tuple[Gate, Gate, Gate]
    """

    G = find_clifford_symmetries(pauli_sum, num_symmetries=1,
                                 dynamic_refine_every=0)
    g = G[0]

    S, T = block_decompose(g.symplectic, int(pauli_sum.lcm), min_block_size=4)
    h_S, h_T = clifford_phase_decomposition(g.symplectic, g.phase_vector, S, T, int(pauli_sum.lcm))
    S_gate = Gate('S', g.qudit_indices, S, g.dimensions, h_S)
    T_gate = Gate('T', g.qudit_indices, T, g.dimensions, h_T)

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
        S, T = block_decompose(g.symplectic, pauli_sum.lcm)
        h_S, h_T = clifford_phase_decomposition(g.symplectic, g.phase_vector, S, T, int(pauli_sum.lcm))
        S_gate = Gate(f'S{i}', g.qudit_indices, S, g.dimensions, h_S)
        T_gate = Gate(f'T{i}', g.qudit_indices, T, g.dimensions, h_T)
        Ss.append(S_gate)
        Ts.append(T_gate)

    return G, Ss, Ts


def block_structure(gate: Gate):
    symp = gate.symplectic
    sizes = np.asarray(ordered_block_sizes(symp, int(gate.lcm)), dtype=int) / 2
    return sizes


def qudit_cost(gate: Gate):
    return int(max(block_structure(gate)))


def _labels_union(independent: list[int], dependencies: dict[int, list[tuple[int, int]]]) -> list[int]:
    return sorted(set(independent) | set(dependencies.keys()))


def find_clifford_symmetries(
    pauli_sum: PauliSum,
    num_symmetries: int = 1,
    # Strategy
    dynamic_refine_every: int = 0,
    extra_column_invariants: str = "none",
    p2_bitset: str = "auto",
    color_mode: str = "wl",
    max_wl_rounds: int = 10,
) -> list[Gate]:
    """
    Return up to k automorphisms preserving S and the vector set. See flags above.
    """
    return clifford_graph_automorphism_search(
        pauli_sum,
        k_wanted=num_symmetries,
        extra_column_invariants=extra_column_invariants,
        p2_bitset=p2_bitset,
        color_mode=color_mode,
        max_wl_rounds=max_wl_rounds,
        dynamic_refine_every=int(dynamic_refine_every),
    )
