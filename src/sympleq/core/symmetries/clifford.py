# import numpy as np
from sympleq.core.circuits import Gate
from sympleq.core.paulis import PauliSum
from sympleq.core.graphs.graph_automorphism import clifford_graph_automorphism_search


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
