"""Codes for finding target Paulis and gates which map a given Pauli to a target Pauli."""
from sympleq.core.paulis import PauliString, PauliSum
import numpy as np
from collections import defaultdict
from itertools import product
from sympleq.core.circuits.find_symplectic import map_pauli_sum_to_target_tableau


# FIXME: improve function name and docstring
def find_map_to_target_pauli_sum(input_pauli: PauliSum, target_pauli: PauliSum) -> tuple[np.ndarray, np.ndarray]:
    """
    TODO: For efficiency improvement act only on target qudits

    Find a gate that maps Pauli P to target Pauli.

    Args:
        input_pauli: PauliSum
            The Pauli to be mapped.
        target_pauli: : PauliSum
            The target Pauli.

    Returns:
        images (list[np.ndarray]): The images of the gate.
        h (np.ndarray): The phase vector of the gate.
        """
    if not np.array_equal(input_pauli.dimensions(), target_pauli.dimensions()):
        raise ValueError(
            f"Input and output PauliSums must have the same dimensions\
            (got {input_pauli.dimensions()} and {target_pauli.dimensions()}).")

    if not np.array_equal(input_pauli.symplectic_product_matrix(), target_pauli.symplectic_product_matrix()):
        raise ValueError("Input and target PauliSum must be symplectically equivalent.")

    input_symplectic = input_pauli.tableau()  # [:, qudit_indices]
    target_symplectic = target_pauli.tableau()  # [:, qudit_indices]

    F = map_pauli_sum_to_target_tableau(input_symplectic, target_symplectic)
    # FIXME: should take input_pauli phases into consideration
    h = get_phase_vector(F, input_pauli.lcm())

    return F, h


def find_allowed_target(pauli_sum, target_pauli_list):
    pauli_list = [target_pauli_list[i][0] for i in range(len(target_pauli_list))]
    pauli_list = str_to_int(pauli_list)
    string_indices = [target_pauli_list[i][1] for i in range(len(target_pauli_list))]
    qudit_indices = [target_pauli_list[i][2] for i in range(len(target_pauli_list))]

    dims = pauli_sum.dimensions()
    combined_indices = list(zip(string_indices, qudit_indices))
    index_dict = defaultdict(list)

    for idx, tup in enumerate(combined_indices):
        index_dict[tup].append(idx)

    underdetermined_pauli_indices = [combined_indices[indexes[0]]
                                     for indexes in index_dict.values() if len(indexes) > 1]
    underdetermined_pauli_options = []
    for indices in underdetermined_pauli_indices:
        options = []
        for i in range(len(target_pauli_list)):
            if target_pauli_list[i][1] == indices[0] and target_pauli_list[i][2] == indices[1]:
                if target_pauli_list[i][0] not in options:
                    options.append(target_pauli_list[i][0])
        underdetermined_pauli_options.append(options)
    underdetermined_pauli_options = [str_to_int(upo) for upo in underdetermined_pauli_options]

    determined_pauli_indices = [combined_indices[indexes[0]] for indexes in index_dict.values() if len(indexes) == 1]

    options_matrix = np.empty([pauli_sum.n_paulis(), pauli_sum.n_qudits()], dtype=object)
    for i in range(pauli_sum.n_paulis()):
        for j in range(pauli_sum.n_qudits()):
            if (i, j) in underdetermined_pauli_indices:
                options_matrix[i, j] = underdetermined_pauli_options[underdetermined_pauli_indices.index((i, j))]
            elif (i, j) in determined_pauli_indices:
                options_matrix[i, j] = [pauli_list[determined_pauli_indices.index((i, j))]]
            else:
                options_matrix[i, j] = [0, 1, 2, 3]

    # flag matrix to track which indices are determined
    flag_matrix = np.zeros([pauli_sum.n_paulis(), pauli_sum.n_qudits()], dtype=int)
    # 0 for not determined, 1 for underdetermined, -1 for determined
    for idx, tup in enumerate(combined_indices):
        if tup in underdetermined_pauli_indices:
            flag_matrix[tup[0], tup[1]] = -1  # -1 for underdetermined
        else:
            flag_matrix[tup[0], tup[1]] = 1  # 1 for determined

    spm = pauli_sum.symplectic_product_matrix()

    possible_targets = []
    for combo in product(*options_matrix.flatten()):
        combo = np.reshape(combo, (pauli_sum.n_paulis(), pauli_sum.n_qudits()))
        pauli_sum_candidate = matrix_to_pauli_sum(combo, pauli_sum.weights, pauli_sum.phases, dims)
        candidate_spm = pauli_sum_candidate.symplectic_product_matrix()

        if np.all(spm == candidate_spm):
            possible_targets.append(pauli_sum_candidate)
    return possible_targets


def get_phase_vector(gate_symplectic: np.ndarray, lcm: int) -> np.ndarray:
    """
    Calculate the phase vector for a gate given its symplectic matrix.

    See PRA 71, 042315 (2005) Eq. (10).
    Solves for h

    Args:
        gate_symplectic (np.ndarray): The symplectic matrix of the gate.
        dimension (int): The dimension of the qudit.

    Returns:
        np.ndarray: The phase vector of the gate.
    """
    n_qudits = gate_symplectic.shape[0] // 2

    U = np.zeros((2 * n_qudits, 2 * n_qudits), dtype=int)
    U[n_qudits:, :n_qudits] = np.eye(n_qudits, dtype=int)
    return (lcm - 1) * np.diag(gate_symplectic.T @ U @ gate_symplectic) % 2


def str_to_int(string):
    output = []

    for s in string:
        if s == 'I':
            output.append(0)
        elif s == 'X':
            output.append(1)
        elif s == 'Y':
            output.append(2)
        elif s == 'Z':
            output.append(3)
    return output


def int_to_pauli(integer):
    if integer == 0:
        return [0, 0]
    elif integer == 1:
        return [1, 0]
    elif integer == 2:
        return [1, 1]
    elif integer == 3:
        return [0, 1]
    else:
        raise ValueError(f"Invalid integer for Pauli representation: {integer}. Must be in [0, 3].")


def matrix_to_pauli_sum(matrix, weights, phases, dimensions):

    n_paulis = matrix.shape[0]
    n_qudits = matrix.shape[1]

    pauli_strings = []
    for i in range(n_paulis):
        x_exp = np.zeros(n_qudits, dtype=int)
        z_exp = np.zeros(n_qudits, dtype=int)
        for j in range(n_qudits):
            x_exp[j] = int_to_pauli(matrix[i, j])[0]
            z_exp[j] = int_to_pauli(matrix[i, j])[1]
        ps = PauliString.from_exponents(x_exp, z_exp, dimensions)
        pauli_strings.append(ps)

    return PauliSum.from_pauli_strings(pauli_strings, weights, phases)
