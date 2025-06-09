import numpy as np
from paulis import PauliSum
print('Warning: This module quaos.pauli_utils is deprecated and will be removed in a future version. Use quaos.paulis instead.')


def check_mappable_via_clifford(pauli_sum: PauliSum, target_pauli_sum: PauliSum) -> bool:
    return bool(np.all(pauli_sum.symplectic_product_matrix() == target_pauli_sum.symplectic_product_matrix()))


def concatenate_pauli_sums(pauli_sums: list[PauliSum]) -> PauliSum:
    """
    Concatenate a list of Pauli sums into a single Pauli sum.
    """
    if len(pauli_sums) == 0:
        raise ValueError("List of Pauli sums is empty")
    if not all(isinstance(p, PauliSum) for p in pauli_sums):
        raise ValueError("All elements of the list must be Pauli sums")
    
    new_pauli_strings = pauli_sums[0].pauli_strings.copy()
    new_dimensions = pauli_sums[0].dimensions.copy()
    new_weights = pauli_sums[0].weights.copy()
    new_phases = pauli_sums[0].phases.copy()
    for p in pauli_sums[1:]:
        new_dimensions = np.concatenate((new_dimensions, p.dimensions))
        new_weights *= p.weights
        new_phases += p.phases
        for i in range(len(new_pauli_strings)):
            new_pauli_strings[i] = new_pauli_strings[i] @ p.pauli_strings[i]

    concatenated = PauliSum(new_pauli_strings, weights=new_weights, phases=new_phases, dimensions=new_dimensions,
                            standardise=False)
    return concatenated


def are_subsets_equal(pauli_sum_1: PauliSum, pauli_sum_2: PauliSum,
                      subset_1: list[tuple[int, int]], subset_2: list[tuple[int, int]] | None = None):
    """
    Check if two subsets of Pauli sums are equal.
    """
    if subset_2 is None:
        subset_2 = subset_1
    else:
        if len(subset_1) != len(subset_2):
            raise ValueError("Subsets must be of the same length")
        if not all(isinstance(i, tuple) and len(i) == 2 for i in subset_1):
            raise ValueError("Subsets must be lists of tuples of length 2")
        if not all(isinstance(i, tuple) and len(i) == 2 for i in subset_2):
            raise ValueError("Subsets must be lists of tuples of length 2")
        
    for i in range(len(subset_1)):
        if pauli_sum_1[subset_1[i]] != pauli_sum_2[subset_2[i]]:
            return False
    return True
