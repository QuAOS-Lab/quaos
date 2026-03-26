import numpy as np
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits import Circuit, SWAP
from collections import defaultdict
from sympleq.utils import int_to_bases, multi_kron
from sympleq.core.circuits.gate_decomposition_to_circuit import gate_to_circuit
from sympleq.core.symmetries.block_decomposition import block_indexes
from sympleq.core.circuits.gates import Gate
from scipy.linalg import block_diag
from itertools import combinations
import itertools


def group_indices(lst):
    """
    Groups indices of the same value in a list into sublists.
    For example, if the input list is [1, 2, 1, 3, 2], the output will be [[0, 2], [1, 4], [3]].
    """
    index_dict = defaultdict(list)
    for idx, value in enumerate(lst):
        index_dict[value].append(idx)

    return [indices for indices in index_dict.values()]


def prime_factors(num):
    factors = []
    factor = 2
    while (num >= 2):
        if (num % factor == 0):
            factors.append(factor)
            num = num / factor
        else:
            factor += 1
    return factors


def flatten(list_of_lists):
    L = []
    for lst in list_of_lists:
        if isinstance(lst, list):
            L += flatten(lst)
        else:
            L.append(lst)
    return L


def permutation_swaps(perm):
    n = len(perm)
    visited = [False] * n
    swaps = []

    for i in range(n):
        if visited[i] or perm[i] == i:
            continue

        cycle = []
        j = i
        while not visited[j]:
            visited[j] = True
            cycle.append(j)
            j = perm[j]

        for k in range(len(cycle) - 1, 0, -1):
            swaps.append((cycle[0], cycle[k]))

    return swaps


def test_symmetry(symmetry, hamiltonian):
    if isinstance(symmetry, Gate):
        return symmetry.act(hamiltonian, tuple(np.arange(hamiltonian.n_qudits()))).is_close(hamiltonian, literal=False)
    elif isinstance(symmetry, Circuit):
        return symmetry.act(hamiltonian).is_close(hamiltonian, literal=False)
    else:
        raise Exception('Symmetry must be a Gate or Circuit')


class ConditionalHamiltonian:
    def __init__(self, hamiltonian: PauliSum, symmetry: Gate, transformation_gate: Gate):
        # original matrices
        self.original_hamiltonian = hamiltonian
        self.symmetry_gate = symmetry
        self.transformation_gate = transformation_gate
        self.symmetrised_hamiltonian = self.transformation_gate.inverse().act(self.original_hamiltonian,
                                                                              tuple(np.arange(hamiltonian.n_qudits())))
        initial_symmetry_circuit = gate_to_circuit(self.symmetry_gate, dimensions=[
                                                   2 for i in range(hamiltonian.n_qudits())])
        initial_T_circuit = gate_to_circuit(self.transformation_gate, dimensions=[
                                            2 for i in range(hamiltonian.n_qudits())])

        # check symmetry
        assert test_symmetry(self.symmetry_gate, self.symmetrised_hamiltonian)

        initial_blocks = block_indexes(self.symmetry_gate.symplectic)
        initial_block_variables = {}
        initial_block_variables['circuits'] = []
        initial_block_variables['unitaries'] = []
        initial_block_variables['eigenvalues'] = []
        initial_block_variables['eigenvectors'] = []

        commuting_blocks = []
        non_commuting_blocks = []

        # find commuting and non commuting blocks as well as their circuits and eigenvalues
        for block_index, block in enumerate(initial_blocks):
            # circuits
            block_circuit = initial_symmetry_circuit.local_circuit(block)
            initial_block_variables['circuits'].append(block_circuit)

            # unitaries
            block_unitary = np.around(block_circuit.unitary().toarray(), 14)
            initial_block_variables['unitaries'].append(block_unitary)

            # eigenvalues and eigenvectors
            if np.isclose(block_unitary, np.eye(len(block_unitary))).all():
                eigenvalues, eigenvectors = np.array([1, 1]), np.array([np.array([1, 0]), np.array([0, 1])])
            else:
                eigenvalues, eigenvectors = np.linalg.eig(block_unitary)

            assert np.all(np.around(eigenvectors.transpose().conj() @ block_unitary @
                          eigenvectors - np.diag(eigenvalues), 10) == 0)

            # sorting eigenvalues and eigenvectors
            idx = np.argsort(eigenvalues.real)
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]

            # saving eigenvalues and eigenvectors
            initial_block_variables['eigenvalues'].append(eigenvalues)
            initial_block_variables['eigenvectors'].append(eigenvectors)

            # find commuting blocks
            embedded_block_C = self.embedded_block_circuit([block_index], initial_blocks, initial_block_variables)
            if test_symmetry(embedded_block_C, self.symmetrised_hamiltonian):
                commuting_blocks.append(block_index)
            else:
                non_commuting_blocks.append(block_index)
            pass

        # find composite commuting blocks
        while len(non_commuting_blocks) > 0:
            commposite_blocks = self.search_commuting_blocks(
                non_commuting_blocks, initial_blocks, initial_block_variables)
            commuting_blocks.insert(0, commposite_blocks)
            for b in commposite_blocks:
                non_commuting_blocks.remove(b)

        # turn list of commuting block indices into list of blocks (with qudit incices)
        source_blocks = []
        for block_index in commuting_blocks:
            if isinstance(block_index, int):
                block = initial_blocks[block_index]
                source_blocks.append(block)
            elif isinstance(block_index, list):
                block = []
                for b in block_index:
                    block += initial_blocks[b]
                source_blocks.append(block)
            else:
                raise Exception('Block index must be an int or list of ints')
        # TODO: Test that all blocks are in-fact commuting

        # reorder such that qudits in commuting blocks are next to each other
        new_qudit_order = flatten(source_blocks.copy())
        swap_list = permutation_swaps(new_qudit_order)
        swap_circuit = Circuit(dimensions=np.array([2 for _ in range(hamiltonian.n_qudits())]),
                               gates=[SWAP() for _ in swap_list],
                               qudit_indices=[(s[0], s[1]) for s in swap_list])
        inv_swap_circuit = swap_circuit.inverse()
        self.ordered_symmetry_circuit = inv_swap_circuit + initial_symmetry_circuit + swap_circuit
        self.ordered_symmetry_gate = self.ordered_symmetry_circuit.composite_gate()
        self.ordered_T_circuit = inv_swap_circuit + initial_T_circuit
        self.ordered_T_gate = self.ordered_T_circuit.composite_gate()
        self.ordered_symmetrised_hamiltonian = swap_circuit.act(self.symmetrised_hamiltonian)
        assert test_symmetry(self.ordered_symmetry_gate, self.ordered_symmetrised_hamiltonian)
        assert self.ordered_T_gate.inverse().act(self.original_hamiltonian,
                                                 qudits=tuple(np.arange(hamiltonian.n_qudits()))
                                                 ).is_close(self.ordered_symmetrised_hamiltonian, literal=False)
        assert test_symmetry(self.ordered_symmetry_gate, self.ordered_T_gate.inverse().act(
            self.original_hamiltonian, qudits=tuple(np.arange(hamiltonian.n_qudits()))))

        # relabel blocks
        self.blocks = []
        qudit_counter = 0
        for block in source_blocks:
            temp_block = []
            for b in block:
                temp_block.append(qudit_counter)
                qudit_counter += 1
            self.blocks.append(temp_block)
        assert len(self.blocks) == len(commuting_blocks)

        # assign circuits, eigenvalues and eigenvectors to each block
        self.block_variables = {}
        self.block_variables['unitaries'] = []
        self.block_variables['eigenvalues'] = []
        self.block_variables['eigenvectors'] = []
        self.block_variables['degeneracies'] = {}
        self.block_variables['degeneracy_dimensions'] = {}
        self.block_variables['conditional_dimensions'] = []

        for block_index, block in enumerate(self.blocks):
            # block parameters
            block_size = len(block)
            hilbert_space_size = 2**block_size
            if isinstance(commuting_blocks[block_index], int):
                block_origin = [commuting_blocks[block_index]]
            elif isinstance(commuting_blocks[block_index], list):
                block_origin = commuting_blocks[block_index]
            else:
                raise Exception('Block origin must be an int or list of ints')

            # collect unitaries
            unitary_list = [initial_block_variables['unitaries'][origin] for origin in block_origin]
            block_unitary = multi_kron(unitary_list)
            self.block_variables['unitaries'].append(block_unitary)

            # collect eigenvalues for component blocks
            list_of_indices = []
            list_of_eigenvalues = []
            list_of_eigenvectors = []
            for origin in block_origin:
                list_of_eigenvalues.append(initial_block_variables['eigenvalues'][origin])
                list_of_eigenvectors.append(initial_block_variables['eigenvectors'][origin])
                list_of_indices.append(np.arange(len(initial_block_variables['eigenvalues'][origin])))

            index_combinations = list(itertools.product(*list_of_indices))
            eigenvalues = np.zeros(len(index_combinations), dtype=complex)
            eigenvectors = np.zeros((hilbert_space_size, hilbert_space_size), dtype=complex)
            for comb_index, combination in enumerate(index_combinations):
                eigenvalues[comb_index] = np.prod([list_of_eigenvalues[ic][comb]
                                                  for ic, comb in enumerate(combination)])
                eigenvectors[:, comb_index] = multi_kron([list_of_eigenvectors[ic][:, comb]
                                                          for ic, comb in enumerate(combination)])

            # sort eigenvalues and eigenvectors
            idx = np.argsort(eigenvalues.real)
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]
            eigenvalues_rounded = np.around(eigenvalues, 10)
            assert np.all(np.around(eigenvectors.T.conj() @ eigenvectors - np.eye(len(eigenvectors)), 10) == 0)
            assert np.all(np.around(eigenvectors.transpose().conj() @ block_unitary @
                          eigenvectors - np.diag(eigenvalues), 10) == 0)

            # save eigenvalues and eigenvectors
            self.block_variables['eigenvalues'].append(eigenvalues)
            self.block_variables['eigenvectors'].append(eigenvectors)
            self.block_variables['degeneracies'][block_index] = group_indices(eigenvalues_rounded)
            d_dimension = [len(g) for g in group_indices(eigenvalues_rounded)]
            self.block_variables['degeneracy_dimensions'][block_index] = d_dimension
            self.block_variables['conditional_dimensions'].append(len(d_dimension))

    def embedded_block_circuit(self, blocks, block_list, block_variables):
        new_dimension = self.original_hamiltonian.dimensions
        embedded_block_circuit = Circuit.empty(dimensions=new_dimension)
        for b in blocks:
            block_circuit = block_variables['circuits'][b]
            block_indexes = block_list[b]
            new_indices = [tuple([block_indexes[i] for i in t]) for t in block_circuit.qudit_indices]
            embedded_block_circuit += Circuit(dimensions=new_dimension,
                                              gates=block_circuit.gates,
                                              qudit_indices=new_indices)
        return embedded_block_circuit

    def search_commuting_blocks(self, non_commuting_blocks, block_list, block_variables):
        initial_index = non_commuting_blocks[0]
        for r in range(1, len(non_commuting_blocks[1:]) + 1):  # subset size
            for subset in combinations(non_commuting_blocks[1:], r):
                embedded_block_C = self.embedded_block_circuit([initial_index] + list(subset),
                                                               block_list, block_variables)
                if test_symmetry(embedded_block_C, self.symmetrised_hamiltonian):
                    return [initial_index] + list(subset)

        raise Exception('No commuting blocks found')

    def selection_to_qudit(self, selection: PauliString, qudit_dimension: int, eigenvectors: list[np.ndarray]):
        partial_hamiltonian_matrix = selection.to_hilbert_space()

        # rewritten relevant pauli string in terms of eigenvector projectors
        G = np.zeros((qudit_dimension, qudit_dimension), dtype=complex)
        for i in range(qudit_dimension):
            for j in range(qudit_dimension):
                G[i, j] = eigenvectors[i].conj() @ partial_hamiltonian_matrix @ eigenvectors[j].T

        # Turn projector matrix into pauli sum
        if qudit_dimension == 1:
            return G[0]
        else:
            P_qudit = PauliSum.from_hilbert_space(G, prime_factors(qudit_dimension)[::-1])
            return P_qudit

    def select_hamiltonian(self, selection: list[int]) -> PauliSum:
        n_p = self.original_hamiltonian.n_paulis()
        decomposed_paulistrings = []
        for i in range(n_p):
            # string
            ps = self.ordered_symmetrised_hamiltonian[i].copy()
            ps.weights = np.array([1])
            ps.phases = np.array([0])
            # weight component
            weight = self.ordered_symmetrised_hamiltonian.weights[i]
            # reconstruct phase component later
            phase = self.ordered_symmetrised_hamiltonian.phases[i]
            previous_lcm = self.ordered_symmetrised_hamiltonian.lcm

            # pauli strings of symmetry blocks that are later tensored
            conditional_pauli_string_blocks = []
            # for each symmetry block decompose the pauli string into the eigenbasis
            for block, selected_eigenvalue in enumerate(selection):
                # get the dimension of the degegeneracy
                d = self.block_variables['degeneracy_dimensions'][block][selected_eigenvalue]
                relevant_qudits = self.blocks[block]
                relevant_indices = np.array(self.block_variables['degeneracies'][block][selected_eigenvalue])
                eigenvectors = [self.block_variables['eigenvectors'][block][:, r] for r in relevant_indices]

                # pauli string of relevant qudits
                ps_selection = ps.copy()
                qubits_to_delete = np.array([q for q in range(ps.n_qudits()) if q not in relevant_qudits])
                ps_selection = ps_selection._delete_qudits(qubits_to_delete)

                P_qudit = self.selection_to_qudit(ps_selection, d, eigenvectors)
                if d == 1:
                    weight *= P_qudit
                else:
                    conditional_pauli_string_blocks.append(P_qudit)

            # tensor the conditional pauli strings together
            if len(conditional_pauli_string_blocks) == 0:
                new_ps = weight * np.exp(phase * 2 * np.pi * 1j / (2 * previous_lcm))
            else:
                new_ps = conditional_pauli_string_blocks[0]
                for i in range(1, len(conditional_pauli_string_blocks)):
                    new_ps = new_ps @ conditional_pauli_string_blocks[i]

                new_ps.weights *= weight
                new_ps.phases += int(phase * new_ps.lcm / previous_lcm)

            decomposed_paulistrings.append(new_ps)

        new_P = decomposed_paulistrings[0]
        for i in range(1, len(decomposed_paulistrings)):
            new_P = new_P + decomposed_paulistrings[i]

        if isinstance(new_P, PauliSum):
            new_P.combine_equivalent_paulis()
            new_P.remove_zero_weight_paulis()
        return new_P

    def diagonalising_unitary(self, inverse=False):
        U_blocks = []
        for block in range(len(self.blocks)):
            eigenvectors = self.block_variables['eigenvectors'][block]
            if inverse:
                U_blocks.append(np.linalg.inv(eigenvectors))
            else:
                U_blocks.append(eigenvectors)

        U = multi_kron(U_blocks)
        return U

    def reconstruct_full_hamiltonian(self):
        list_of_blocks = []
        for block_hamiltonian in self:
            if isinstance(block_hamiltonian, PauliSum):
                list_of_blocks.append(block_hamiltonian.to_hilbert_space().toarray())
            elif isinstance(block_hamiltonian, complex) or isinstance(block_hamiltonian, float):
                list_of_blocks.append([[block_hamiltonian]])
            elif isinstance(block_hamiltonian, np.ndarray):
                list_of_blocks.append([block_hamiltonian])
            else:
                raise ValueError('Unexpected type of block hamiltonian')

        full_hamiltonian_matrix = block_diag(*list_of_blocks)
        U = self.diagonalising_unitary(inverse=True)
        reconstructed_H = U.T.conj() @ full_hamiltonian_matrix @ U
        return reconstructed_H

    def test_conditional_hamiltonian(self):
        comparison_hamiltonian = self.ordered_symmetrised_hamiltonian.to_hilbert_space().toarray()
        reconstructed_hamiltonian = self.reconstruct_full_hamiltonian()
        return bool(np.max(np.abs(comparison_hamiltonian - reconstructed_hamiltonian)) < 1e-10)

    def __iter__(self):
        n = np.prod(self.block_variables['conditional_dimensions'])
        for i in range(n):
            selection = list(int_to_bases(i, self.block_variables['conditional_dimensions']))
            block_hamiltonian = self.select_hamiltonian(selection)
            yield block_hamiltonian
