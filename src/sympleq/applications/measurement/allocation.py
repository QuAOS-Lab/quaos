import numpy as np
import random
import itertools
from sympleq.core.paulis import PauliSum
from sympleq.applications.measurement.covariance_graph import graph
from sympleq.core.circuits import Circuit
from sympleq.core.circuits.gates import GATES

# Convenience aliases
H = GATES.H
CX = GATES.CX
S = GATES.S


def mcmc_number_initial_samples(shots: int, n_0: int = 500, scaling_factor: float = 1 / 10000) -> int:
    """
    Standard method to calculate the number of initial samples for the MCMC algorithm.

    Parameters:
        shots (int): The number of shots to be taken in the experiment.
        n_0 (int): The base number of initial samples to be taken.
        scaling_factor (float): A factor to scale the number of initial samples by.

    Returns:
        int: The total number of initial samples to be taken.
    """
    return int(shots * scaling_factor) + n_0


def mcmc_number_max_samples(shots: int, n_0: int = 2001, scaling_factor: float = 1 / 10000) -> int:
    """
    Standard method to calculate the maximum number of samples for the MCMC algorithm.

    Parameters:
        shots (int): The number of shots to be taken in the experiment.
        n_0 (int): The base number of maximum samples to be taken.
        scaling_factor (float): A factor to scale the number of maximum samples by.

    Returns:
        int: The total number of maximum samples to be taken.
    """
    return 4 * int(shots * scaling_factor) + n_0


def sort_hamiltonian(P: PauliSum) -> tuple[PauliSum, np.ndarray]:
    """
    Sorts the Hamiltonian's Pauli operators based on hermiticity, with hermitian ones first and then pairs of
    Paulis and their hermitian conjugate. !!! Also removes identity !!!

    Parameters:
        P (pauli): A set of Pauli operators.
        coefficients (list): Corresponding coefficients of the Pauli operators.

    Returns:
        tuple: Sorted Pauli operators, coefficients, and the size of Pauli blocks.
    """
    pauli_count = P.n_paulis()
    indices = list(range(pauli_count))

    hermitian_indices = []
    non_hermitian_indices = []
    paired_conjugates = []

    while indices:
        i = indices.pop(0)
        P0 = P[i, :]
        P0_str = str(P0)
        P0_conjugate = P0.hermitian_conjugate()
        P0_conj_str = str(P0_conjugate)

        if P0_str == P0_conj_str:
            if P0.x_exp.any() or P0.z_exp.any():
                hermitian_indices.append(i)
            continue
        else:
            non_hermitian_indices.append(i)

        for j in indices:
            P1 = P[j, :]
            P1_str = str(P1)
            if P0_conj_str == P1_str:
                paired_conjugates.append(j)
                indices.remove(j)
                break

    # Rebuild Pauli set and coefficients
    sorted_indices = []
    pauli_block_sizes = []

    # Handle hermitian indices
    for i in hermitian_indices:
        sorted_indices.append(i)
        pauli_block_sizes.append(1)

    # Handle non-hermitian indices and their conjugates
    for i, j in zip(non_hermitian_indices, paired_conjugates):
        sorted_indices.extend([i, j])
        pauli_block_sizes.append(2)

    # Extract and reorder Pauli strings and coefficients
    P1 = P.copy()
    P1.reorder(sorted_indices)

    # Remove identity
    for i in range(P1.n_paulis()):
        if P1.weights[i] != 0 and P1[i, :].is_identity():
            P1._delete_paulis([i])

    return P1, np.array(pauli_block_sizes)


def choose_measurement(S: np.ndarray, V: np.ndarray, aaa, allocation_mode: str) -> list[int]:
    """
    Choose a set of Pauli operators (based on the covariance graph) to measure in order to reduce the estimation error
    for the observable.

    Parameters
    ----------
    S : np.ndarray
        Covariance matrix of the Pauli operators.
    V : np.ndarray
        The covariance graph.
    aaa : list
        A list of sets of Pauli operators.
    allocation_mode : str
        The allocation mode. Can be either 'set' or 'rand'.

    Returns
    -------
    list
        A list of indices of the chosen Pauli operators.

    Raises
    ------
    ValueError
        If the allocation mode is not implemented.
    """
    p = S.shape[0]
    Ones = [np.ones((i, i), dtype=int) for i in range(p + 1)]
    index_set = set(range(p))
    S1 = S + Ones[p]
    s = 1 / (S.diagonal() | (S.diagonal() == 0))
    s1 = 1 / S1.diagonal()
    factor = p - np.count_nonzero(S.diagonal())
    S1[range(p), range(p)] = [a if a != 1 else -factor for a in S1.diagonal()]
    V1 = V * (S * s * s[:, None] - S1 * s1 * s1[:, None])
    V2 = 2 * V * (S * s * s[:, None] - S * s * s1[:, None])
    aaa, aaa1 = itertools.tee(aaa, 2)
    if allocation_mode == 'set':
        aa = sorted(max(aaa1, key=lambda xx: np.abs(
            V1[xx][:, xx].sum() + V2[xx][:, list(index_set.difference(xx))].sum())))
    elif allocation_mode == 'rand':
        aa = sorted(random.sample(list(set(frozenset(aa1) for aa1 in aaa1)), 1)[0])
    else:
        raise NotImplementedError('Allocation mode not implemented')
    return aa


def construct_circuit_list(P: PauliSum, xxx: list[list[int]], D: dict) -> tuple[list[Circuit], dict]:
    """
    Constructs a list of circuits which diagonalize the cliques of the given Hamiltonian.

    Parameters:
        P (PauliSum): The Hamiltonian to be diagonalized.
        xxx (list): A list of sets of Pauli operators to be diagonalized.
        D (dict): A dictionary mapping the string representation of a clique to a diagonalization circuit.

    Returns:
        tuple: A list of circuits which diagonalize the given Hamiltonian and the updated dictionary.

    Raises:
        ValueError: If a circuit is not diagonalizing.
    """
    circuit_list = []
    for aa in xxx:
        C, D = construct_diagonalization_circuit(P, aa, D=D)
        if not is_diagonalizing_circuit(P, C, aa):
            raise ValueError('Circuit is not diagonalizing')
        circuit_list.append(C)
    return circuit_list, D


def construct_diagonalization_circuit(P: PauliSum, aa: list[int], D={}) -> tuple[Circuit, dict]:
    """
    Constructs a circuit which diagonalizes the singular given clique of a Hamiltonian.

    Parameters:
        P (PauliSum): The Hamiltonian to be diagonalized.
        aa (list): A list of indices of Pauli operators to be diagonalized.
        D (dict): A dictionary mapping the string representation of a clique to a diagonalization circuit.

    Returns:
        tuple: A circuit which diagonalizes the given clique and the updated dictionary.

    Raises:
        ValueError: If a circuit is not diagonalizing.
    """
    if str(aa) in D:
        P1, C, k_dict = D[str(aa)]
    else:
        P1 = P.copy()
        P1._delete_paulis([i for i in range(P.n_paulis()) if i not in aa])
        # add products
        k_dict = {str(j0): [(a0, a0, P1.phases[j0])] for j0, a0 in enumerate(aa)}
        for j0, a0 in enumerate(aa):
            for j1, a1 in enumerate(aa):
                if j0 != j1:
                    # isolate a pair of paulis in the clique
                    P_a0 = P1[[j0]]
                    P_a0c = P_a0.H()
                    P_a1 = P1[[j1]]
                    # compute their product pauli
                    P2 = P_a0c * P_a1
                    P2.weights[0] = 1
                    # check if the product is in the original pauli list
                    P1_pauli_string_list = [str(P1[k]) for k in range(P1.n_paulis())]
                    if str(P2[0]) not in P1_pauli_string_list:
                        k_dict[str(P1.n_paulis())] = [(a0, a1, P2.phases[0])]
                        # add the product but make sure to account for possibly different phases
                        P2.phases[0] = 0
                        P1 = P1 + P2
                    else:
                        k_dict[str(P1_pauli_string_list.index(str(P2[0])))].append((a0, a1, P2.phases[0]))

        C = diagonalize(P1)
        P1.set_phases(np.zeros(P1.n_paulis()))
        P1 = C.act(P1)
        D[str(aa)] = (P1, C, k_dict)
    return C, D


# TODO: Make this part of the core function
def diagonalize(P: PauliSum) -> Circuit:
    """
    Diagonalize a PauliSum.

    This function takes a PauliSum P and returns a Circuit C such that C.act(P) is diagonal.
    If the PauliSum is already diagonal, it simply returns the input Circuit.

    If the PauliSum is not diagonal but is quditwise commuting, it diagonalizes the PauliSum by
    calling diagonalize_iter_quditwise_ on the qudits of the same dimension.

    If the PauliSum is not diagonal and is not quditwise commuting, it diagonalizes the PauliSum by
    calling diagonalize_iter_ on the qudits of the same dimension.

    Finally, it applies H gates to make any X qudits Z.

    Parameters:
        P (PauliSum): The PauliSum to be diagonalized.

    Returns:
        Circuit: The diagonalizing circuit.
    """
    dims = P.dimensions
    q = P.n_qudits()

    if len(set(P.dimensions)) == 1:
        # Currently commutation check is not working for mixed species
        if not P.is_commuting():
            raise Exception("Paulis must be pairwise commuting to be diagonalized")
    P1 = P.copy()
    C = Circuit.empty(dims)

    if P.is_quditwise_commuting():
        # for each dimension, call diagonalize_iter_quditwise_ on the qudits of the same dimension
        for d in sorted(set(dims)):
            aa = [i for i in range(q) if dims[i] == d]
            while aa:
                C = diagonalize_iter_quditwise_(P1, C, aa)
        P1 = C.act(P1)
    else:
        # for each dimension, call diagonalize_iter_ on the qudits of same dimension
        for d in sorted(set(dims)):
            aa = [i for i in range(q) if dims[i] == d]
            while aa:
                C = diagonalize_iter_(P1, C, aa)
        P1 = C.act(P1)

    # if any qudits are X rather than Z, apply H to make them Z
    if [i for i in range(q) if any(P1.x_exp[:, i])]:
        C1 = Circuit.empty(dims)
        for i in range(q):
            if any(P1.x_exp[:, i]):
                C.add_gate(H, i)
                C1.add_gate(H, i)
        P1 = C1.act(P1)
    return C


def diagonalize_iter_(P: PauliSum, C: Circuit, aa: list[int]) -> Circuit:
    """
    Diagonalize the given PauliSum on the given qudit indices.

    Parameters
    ----------
    P : PauliSum
        The PauliSum to be diagonalized.
    C : Circuit
        The circuit to be modified in place.
    aa : list[int]
        The qudit indices to be diagonalized.

    Returns
    -------
    Circuit
        The modified circuit.
    """
    p = P.n_paulis()
    P = C.act(P)
    a = aa.pop(0)

    # if all Paulis have no X-part on qudit a, return C
    if not any(P.x_exp[:, a]):
        return C

    # set a1 to be the index of the minimum Pauli with non-zero X-part of qudit a
    a1 = min(i for i in range(p) if P.x_exp[i, a])

    # add CNOT gates to cancel out all non-zero X-parts on Pauli a1, qudits in aa
    while any(P.x_exp[a1, i] for i in aa):
        C1 = Circuit.empty(P.dimensions)
        for i in aa:
            if P.x_exp[a1, i]:
                C.add_gate(CX, a, i)
                C1.add_gate(CX, a, i)
        P = C1.act(P)

    # check whether there are any non-zero Z-parts on Pauli a1, qudits in aa
    while any(P.z_exp[a1, i] for i in aa):

        # if Pauli a1, qudit a is X, apply S gate to make it Y
        if not P.z_exp[a1, a]:
            C.add_gate(S, a)
            P = S.act(P, a)

        # add backwards CNOT gates to cancel out all non-zero Z-parts on Pauli a1, qudits in aa
        for i in aa:
            if P.z_exp[a1, i]:
                C.add_gate(CX, i, a)
                P = CX.act(P, (i, a))

    # if Pauli a1, qudit a is Y, add S gate to make it X
    while P.z_exp[a1, a]:
        C.add_gate(S, a)
        P = S.act(P, a)
    return C


def diagonalize_iter_quditwise_(P: PauliSum, C: Circuit, aa: list[int]) -> Circuit:
    """
    Diagonalize the given PauliSum on the given qudit indices in a quditwise manner.

    Parameters
    ----------
    P : PauliSum
        The PauliSum to be diagonalized.
    C : Circuit
        The circuit to be modified in place.
    aa : list[int]
        The qudit indices to be diagonalized.

    Returns
    -------
    Circuit
        The modified circuit.
    """
    p = P.n_paulis()
    P = C.act(P)
    a = aa.pop(0)

    # if all Paulis have no X-part on qudit a, return C
    if not any(P.x_exp[:, a]):
        return C

    # set a1 to be the index of the minimum Pauli with non-zero X-part of qudit a
    a1 = min(i for i in range(p) if P.x_exp[i, a])

    # if Pauli a1, qudit a is Y, add S gate to make it X
    while P.z_exp[a1, a]:
        C.add_gate(S, a)
        P = S.act(P, a)
    return C


def is_diagonalizing_circuit(P: PauliSum, C: Circuit, aa: list[int]) -> bool:
    """
    Checks whether the given circuit diagonalizes the given PauliSum on the given qudit indices.

    Parameters
    ----------
    P : PauliSum
        The PauliSum to be checked.
    C : Circuit
        The circuit to be checked.
    aa : list[int]
        The qudit indices to be checked.

    Returns
    -------
    bool
        Whether the circuit diagonalizes the PauliSum on the given qudit indices.
    """
    P1 = P.copy()
    P1._delete_paulis([i for i in range(P.n_paulis()) if i not in aa])
    P1 = C.act(P1)
    return P1.is_z()


def update_data(xxx: list[list[int]], rr: list[np.ndarray], X: np.ndarray, D: dict) -> np.ndarray:
    """
    Updates the given data matrix with the given measurement results, PauliSums and dictionaries.

    Parameters
    ----------
    xxx : list[str]
        The list of PauliSum keys to be updated.
    rr : list[np.ndarray]
        The list of measurement results to be updated.
    X : np.ndarray
        The data matrix to be updated.
    D : dict
        The dictionary of PauliSums and dictionaries to be updated.

    Returns
    -------
    np.ndarray
        The updated data matrix.
    """
    d = len(X[0, 0])
    for i, aa in enumerate(xxx):
        (P1, _, k_dict) = D[str(aa)]
        p1, q1, phases1 = P1.n_paulis(), P1.n_qudits(), P1.phases
        bases_a1 = rr[i]
        ss = [(sum((bases_a1[i1] * P1.z_exp[i0, i1] * P1.lcm) // P1.dimensions[i1]
               for i1 in range(q1))) % P1.lcm for i0 in range(p1)]
        for j0, s0 in enumerate(ss):
            for a0, a1, s1 in k_dict[str(j0)]:
                if (phases1[j0] + s1) % 2 == 1:
                    raise Exception('Odd phase detected for sorting into data matrix')
                X[a0, a1, int(s0 + (phases1[j0] + s1) / 2) % d] += 1
    return X


def update_diagnostic_data(cliques: list[list[int]], diagnostic_results: list[np.ndarray],
                           diagnostic_data: np.ndarray, mode: str = 'zero') -> np.ndarray:
    """
    Updates the given diagnostic data matrix with the given measurement results and PauliSums.

    Parameters
    ----------
    cliques : list[int]
        The list of PauliSum keys to be updated.
    diagnostic_results : list[np.ndarray]
        The list of measurement results to be updated.
    diagnostic_data : np.ndarray
        The diagnostic data matrix to be updated.
    mode : str, default='zero'
        The mode of updating the diagnostic data matrix. Options are 'zero' and 'random'.

    Returns
    -------
    np.ndarray
        The updated diagnostic data matrix.
    """
    if mode == 'zero':
        for i in range(len(diagnostic_results)):
            if np.any(diagnostic_results[i]):
                diagnostic_data[cliques[i], 0] += 1
            else:
                diagnostic_data[cliques[i], 1] += 1
    elif mode == 'random':
        raise Exception('Random diagnostic states not yet implemented')
    else:
        raise Exception('Diagnostic state mode not recognized')
    return diagnostic_data


def scale_variances(A: graph, S: np.ndarray) -> graph:
    # Inputs:
    #     A - (graph)       - variance matrix
    #     S - (numpy.array) - array for tracking number of measurements
    """
    Scales the variance matrix by the number of measurements.

    Parameters
    ----------
    A : graph
        The variance matrix to be scaled.
    S : numpy.array
        The array for tracking the number of measurements.

    Returns
    -------
    graph
        The scaled variance matrix.
    """
    p = A.ord()
    S1 = S.copy()
    S1[range(p), range(p)] = [a if a != 0 else 1 for a in S1.diagonal()]
    s1 = 1 / S1.diagonal()
    return graph(S1 * A.adj * s1 * s1[:, None])


def construct_diagnostic_circuits(circuit_list: list[Circuit]) -> list[Circuit]:
    """
    Constructs a list of diagnostic circuits (for hardware error measurements) from a given list of circuits.

    Parameters
    ----------
    circuit_list : list[Circuit]
        The list of circuits to be converted into diagnostic circuits.

    Returns
    -------
    list[Circuit]
        The list of constructed diagnostic circuits.
    """
    diagnostic_circuits = []
    for circ in circuit_list:
        C_diag = Circuit.empty(dimensions=circ.dimensions)
        for gate, qudits in zip(circ.gates, circ.qudit_indices):
            C_diag.add_gate(gate, *qudits)
            if gate.name == 'H':
                C_diag.add_gate(gate, *qudits)
                C_diag.add_gate(gate, *qudits)
                C_diag.add_gate(gate, *qudits)
        diagnostic_circuits.append(C_diag)
    return diagnostic_circuits


def construct_diagnostic_states(diagnostic_circuits: list[Circuit],
                                mode: str = 'zero') -> tuple[list[np.ndarray], list[Circuit]]:
    """
    Constructs a list of diagnostic states (for hardware error measurements) from a given list of diagnostic circuits.

    Parameters
    ----------
    diagnostic_circuits : list[Circuit]
        The list of diagnostic circuits to be used for constructing the diagnostic states.
    mode : str, optional
        The mode of constructing the diagnostic states. Can be either 'zero' or 'random'. If 'zero', the diagnostic
        states are the computational basis states. If 'random', the diagnostic states are random states
        (not yet implemented). Defaults to 'zero'.

    Returns
    -------
    tuple
        A tuple containing the list of constructed diagnostic states and the list of constructed state preparation
        circuits.
    """
    if mode == 'zero':
        state = [0] * np.prod(diagnostic_circuits[0].dimensions)
        state[0] = 1
        state = np.array(state)
        state_preparation_circuits = [Circuit.empty(diagnostic_circuits[0].dimensions)] * len(diagnostic_circuits)
        return ([state] * len(diagnostic_circuits), state_preparation_circuits)
    elif mode == 'random':
        raise Exception('Random diagnostic states not yet implemented')
    else:
        raise Exception('Diagnostic state mode not recognized')


def standard_noise_probability_function(circuit: Circuit, p_entangling=0.03,
                                        p_local=0.001, p_measurement=0.001) -> float:
    """
    Standard example function to calculate the probability of a noise event occurring in a given circuit.

    Parameters
    ----------
    circuit : Circuit
        The circuit for which the noise probability is to be calculated.
    p_entangling : float, optional
        The probability of an entangling gate causing a noise event. Defaults to 0.03.
    p_local : float, optional
        The probability of a local gate causing a noise event. Defaults to 0.001.
    p_measurement : float, optional
        The probability of a measurement gate causing a noise event. Defaults to 0.001.

    Returns
    -------
    float
        The calculated probability of a noise event occurring in the given circuit.
    """
    n_local = 0
    n_entangling = 0
    for gate in circuit.gates:
        if gate.name == 'CX':
            n_entangling += 1
        else:
            n_local += 1
    noise_prob = 1 - ((1 - p_measurement) * (1 - p_entangling)**n_entangling * (1 - p_local)**n_local)
    return noise_prob


def standard_error_function(result: np.ndarray, dimensions: list[int]) -> np.ndarray:
    """
    Standard example function to simulate the statistical error of a given result.

    Parameters
    ----------
    result : np.ndarray
        The result for which the statistical error is to be simulated.
    dimensions : list
        The dimensions of the result.

    Returns
    -------
    np.ndarray
        The calculated statistical error of the given result.
    """
    return np.array([np.random.randint(dimensions[j]) for j in range(len(dimensions))])


def extract_phase(weight: complex, dimension: int) -> tuple[int, float]:
    """
    Extracts the phase and remainder of a given weight with respect to a given dimension.

    Parameters
    ----------
    weight : complex
        The weight for which the phase and remainder are to be extracted.
    dimension : int
        The dimension with respect to which the phase and remainder are to be extracted.

    Returns
    -------
    tuple[int, float]
        A tuple containing the extracted phase and remainder.
    """
    phase = np.floor(dimension * np.angle(weight) / (2 * np.pi))
    remainder = np.angle(weight) - phase * 2 * np.pi / dimension
    return phase, remainder
