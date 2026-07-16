import numpy as np
import scipy.sparse as sp
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits import Circuit, SWAP
from sympleq.utils import int_to_bases, bases_to_int, multi_kron
from sympleq.core.circuits.gate_decomposition_to_circuit import gate_to_circuit
from sympleq.core.symmetries.block_decomposition import block_indexes
from sympleq.core.circuits.gates import Gate
from scipy.linalg import block_diag
from itertools import combinations
import itertools


def group_indices_by_tolerance(values: np.ndarray, tol: float) -> list[list[int]]:
    """
    Group indices of complex eigenvalues by absolute-distance tolerance.
    """
    vals = np.asarray(values, dtype=complex).reshape(-1)
    used = np.zeros(vals.shape[0], dtype=bool)
    groups: list[list[int]] = []
    for i, v in enumerate(vals):
        if used[i]:
            continue
        group = [i]
        used[i] = True
        for j in range(i + 1, vals.shape[0]):
            if used[j]:
                continue
            if np.abs(vals[j] - v) <= float(tol):
                group.append(j)
                used[j] = True
        groups.append(group)
    return groups


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
    def __init__(
        self,
        hamiltonian: PauliSum,
        symmetry: Gate | Circuit,
        transformation_gate: Gate | Circuit,
        verify: bool = True,
        degeneracy_tol: float = 1e-8,
        zero_tol: float = 1e-12,
    ):
        # original matrices
        self.original_hamiltonian = hamiltonian
        self.symmetry_gate = symmetry
        self.transformation_gate = transformation_gate
        self.verify = bool(verify)
        self.degeneracy_tol = float(degeneracy_tol)
        self.zero_tol = float(zero_tol)
        self.input_dimensions = np.asarray(hamiltonian.dimensions, dtype=int)

        self.symmetry_gate, self.symmetry_circuit = self._coerce_gate_and_circuit(symmetry)
        self.transformation_gate, self.transformation_circuit = self._coerce_gate_and_circuit(transformation_gate)

        self.symmetrised_hamiltonian = self.transformation_circuit.inverse().act(self.original_hamiltonian)
        initial_symmetry_circuit = self.symmetry_circuit
        initial_T_circuit = self.transformation_circuit

        # check symmetry
        if self.verify:
            assert test_symmetry(self.symmetry_circuit, self.symmetrised_hamiltonian)

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
            if np.isclose(block_unitary, np.eye(block_unitary.shape[0], dtype=complex)).all():
                block_hilbert_space_size = int(block_unitary.shape[0])
                eigenvalues = np.ones(block_hilbert_space_size, dtype=complex)
                eigenvectors = np.eye(block_hilbert_space_size, dtype=complex)
            else:
                eigenvalues, eigenvectors = np.linalg.eig(block_unitary)

            if self.verify:
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
            composite_blocks = self.search_commuting_blocks(
                non_commuting_blocks, initial_blocks, initial_block_variables)
            commuting_blocks.insert(0, composite_blocks)
            for b in composite_blocks:
                non_commuting_blocks.remove(b)

        # turn list of commuting block indices into list of blocks (with qudit indices)
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
        new_qudit_order = [qudit for block in source_blocks for qudit in block]
        swap_list = permutation_swaps(new_qudit_order)
        swap_circuit = Circuit(dimensions=self.input_dimensions,
                               gates=[SWAP() for _ in swap_list],
                               qudit_indices=[(s[0], s[1]) for s in swap_list])
        inv_swap_circuit = swap_circuit.inverse()
        self.ordered_symmetry_circuit = inv_swap_circuit + initial_symmetry_circuit + swap_circuit
        self.ordered_symmetry_gate = self.ordered_symmetry_circuit.composite_gate()
        self.ordered_T_circuit = inv_swap_circuit + initial_T_circuit
        self.ordered_T_gate = self.ordered_T_circuit.composite_gate()
        self.ordered_symmetrised_hamiltonian = swap_circuit.act(self.symmetrised_hamiltonian)
        if self.verify:
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
        if self.verify:
            assert len(self.blocks) == len(commuting_blocks)

        # assign circuits, eigenvalues and eigenvectors to each block
        self.block_variables = {}
        self.block_variables['unitaries'] = []
        self.block_variables['eigenvalues'] = []
        self.block_variables['eigenvectors'] = []
        self.block_variables['degeneracies'] = {}
        self.block_variables['degeneracy_dimensions'] = {}
        self.block_variables['conditional_dimensions'] = []
        self.P_qudit_dict = {}
        self._block_qudits_to_delete = [
            np.array(
                [q for q in range(self.ordered_symmetrised_hamiltonian.n_qudits()) if q not in block],
                dtype=int,
            )
            for block in self.blocks
        ]
        self._ordered_single_terms_cache = None
        self._term_block_selections_cache = None
        self._term_block_selection_keys_cache = None

        for block_index, block in enumerate(self.blocks):
            # block parameters
            block_dimensions = np.asarray(self.ordered_symmetrised_hamiltonian.dimensions[block], dtype=int)
            hilbert_space_size = int(np.prod(block_dimensions))
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
            degeneracy_groups = group_indices_by_tolerance(eigenvalues, self.degeneracy_tol)
            if self.verify:
                assert np.all(np.around(eigenvectors.T.conj() @ eigenvectors - np.eye(len(eigenvectors)), 10) == 0)
                assert np.all(np.around(eigenvectors.transpose().conj() @ block_unitary @
                              eigenvectors - np.diag(eigenvalues), 10) == 0)

            # save eigenvalues and eigenvectors
            self.block_variables['eigenvalues'].append(eigenvalues)
            self.block_variables['eigenvectors'].append(eigenvectors)
            self.block_variables['degeneracies'][block_index] = [list(g) for g in degeneracy_groups]
            d_dimension = [len(g) for g in degeneracy_groups]
            self.block_variables['degeneracy_dimensions'][block_index] = d_dimension
            self.block_variables['conditional_dimensions'].append(len(d_dimension))

        self._sector_basis_indices = None
        self._sector_permutation = None
        self._inverse_sector_permutation = None

    def _coerce_gate_and_circuit(self, operator: Gate | Circuit) -> tuple[Gate, Circuit]:
        if isinstance(operator, Circuit):
            return operator.composite_gate(), operator
        return operator, gate_to_circuit(operator, dimensions=self.input_dimensions)

    def _initialize_sector_index_maps(self) -> None:
        block_hilbert_dimensions = [
            int(len(self.block_variables['eigenvalues'][block]))
            for block in range(len(self.blocks))
        ]
        conditional_dimensions = [
            int(x) for x in self.block_variables.get('conditional_dimensions', [])
        ]
        n_sectors = int(np.prod(conditional_dimensions)) if len(conditional_dimensions) > 0 else 1

        sector_basis_indices: list[np.ndarray] = []
        for sector_idx in range(n_sectors):
            selection = [int(x) for x in int_to_bases(sector_idx, conditional_dimensions)]
            selected_groups = [
                [int(i) for i in self.block_variables['degeneracies'][block][selection[block]]]
                for block in range(len(self.blocks))
            ]

            local_indices: list[int] = []
            for basis_tuple in itertools.product(*selected_groups):
                local_indices.append(int(bases_to_int(np.asarray(basis_tuple, dtype=int), block_hilbert_dimensions)))
            sector_basis_indices.append(np.asarray(local_indices, dtype=int))

        if len(sector_basis_indices) == 0:
            sector_permutation = np.array([], dtype=int)
        else:
            sector_permutation = np.concatenate(sector_basis_indices).astype(int)

        hilbert_dim = int(np.prod(block_hilbert_dimensions))
        if sector_permutation.size != hilbert_dim:
            raise RuntimeError(
                "Internal sector-index map size mismatch: "
                f"got {sector_permutation.size}, expected {hilbert_dim}."
            )
        if hilbert_dim > 0 and not np.array_equal(np.sort(sector_permutation), np.arange(hilbert_dim, dtype=int)):
            raise RuntimeError("Internal sector-index map is not a valid permutation.")

        self._sector_basis_indices = sector_basis_indices
        # sector_permutation[k] gives the kron-order basis index at sector-order position k.
        self._sector_permutation = sector_permutation
        # inverse map: kron-order -> sector-order.
        self._inverse_sector_permutation = np.argsort(sector_permutation).astype(int)

    def sector_basis_indices(self) -> list[np.ndarray]:
        if self._sector_basis_indices is None:
            self._initialize_sector_index_maps()
        return [np.asarray(idx, dtype=int).copy() for idx in self._sector_basis_indices]

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
        if hasattr(partial_hamiltonian_matrix, "toarray"):
            partial_hamiltonian_matrix = partial_hamiltonian_matrix.toarray()
        partial_hamiltonian_matrix = np.asarray(partial_hamiltonian_matrix, dtype=complex)

        V = np.column_stack(eigenvectors).astype(complex)

        # Rewritten relevant Pauli string in terms of selected-eigenspace projectors.
        G = V.conj().T @ partial_hamiltonian_matrix @ V

        # Turn projector matrix into pauli sum
        if qudit_dimension == 1:
            return G[0, 0]
        else:
            # Guard: a numerically zero matrix cannot be expanded into a non-empty PauliSum.
            if float(np.max(np.abs(G))) < self.zero_tol:
                return 0.0 + 0.0j
            P_qudit = PauliSum.from_hilbert_space(G, prime_factors(qudit_dimension)[::-1])
            return P_qudit

    def _selection_key(self, selection: list[int] | np.ndarray | tuple[int, ...]) -> tuple[int, ...]:
        return tuple(int(x) for x in selection)

    def _selection_output_dimensions(self, selection_key: tuple[int, ...]) -> tuple[list[int], bool]:
        output_dimensions: list[int] = []
        for block, selected_eigenvalue in enumerate(selection_key):
            d = int(self.block_variables['degeneracy_dimensions'][block][selected_eigenvalue])
            if d > 1:
                output_dimensions += [int(x) for x in prime_factors(d)[::-1]]
        return output_dimensions, len(output_dimensions) == 0

    def _selection_block_data(
        self,
        selection_key: tuple[int, ...],
    ) -> list[tuple[int, list[np.ndarray], tuple[int, int]]]:
        block_data: list[tuple[int, list[np.ndarray], tuple[int, int]]] = []
        for block, selected_eigenvalue in enumerate(selection_key):
            relevant_indices = np.array(self.block_variables['degeneracies'][block][selected_eigenvalue])
            eigenvectors = [self.block_variables['eigenvectors'][block][:, r] for r in relevant_indices]
            block_data.append(
                (
                    int(self.block_variables['degeneracy_dimensions'][block][selected_eigenvalue]),
                    eigenvectors,
                    (int(block), int(selected_eigenvalue)),
                )
            )
        return block_data

    def _pauli_string_cache_key(self, pauli_string: PauliString) -> tuple[tuple[int, ...], tuple[int, ...]]:
        return (
            tuple(int(x) for x in np.asarray(pauli_string.tableau, dtype=int).reshape(-1)),
            tuple(int(x) for x in np.asarray(pauli_string.dimensions, dtype=int).reshape(-1)),
        )

    def _ordered_single_terms(self) -> list[PauliString]:
        if self._ordered_single_terms_cache is None:
            ordered_single_terms: list[PauliString] = []
            for i in range(self.ordered_symmetrised_hamiltonian.n_paulis()):
                ps = self.ordered_symmetrised_hamiltonian[i].copy()
                ps.weights = np.array([1])
                ps.phases = np.array([0])
                ordered_single_terms.append(ps)
            self._ordered_single_terms_cache = ordered_single_terms
        return self._ordered_single_terms_cache

    def _term_block_selections(
        self,
    ) -> tuple[list[list[PauliString]], list[list[tuple[tuple[int, ...], tuple[int, ...]]]]]:
        if self._term_block_selections_cache is None or self._term_block_selection_keys_cache is None:
            ordered_single_terms = self._ordered_single_terms()
            term_block_selections: list[list[PauliString]] = []
            term_block_selection_keys: list[list[tuple[tuple[int, ...], tuple[int, ...]]]] = []
            for qudits_to_delete in self._block_qudits_to_delete:
                block_terms: list[PauliString] = []
                block_keys: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
                for ps in ordered_single_terms:
                    ps_selection = ps._delete_qudits(qudits_to_delete)
                    block_terms.append(ps_selection)
                    block_keys.append(self._pauli_string_cache_key(ps_selection))
                term_block_selections.append(block_terms)
                term_block_selection_keys.append(block_keys)
            self._term_block_selections_cache = term_block_selections
            self._term_block_selection_keys_cache = term_block_selection_keys
        return self._term_block_selections_cache, self._term_block_selection_keys_cache

    def select_hamiltonian(self, selection: list[int]) -> PauliSum | complex:
        selection_key = self._selection_key(selection)

        def _zero_paulisum(dimensions: list[int]) -> PauliSum:
            x0 = np.zeros(len(dimensions), dtype=int)
            z0 = np.zeros(len(dimensions), dtype=int)
            identity_ps = PauliString.from_exponents(x0, z0, dimensions)
            return PauliSum.from_pauli_strings(identity_ps, weights=[0.0 + 0.0j], phases=[0])

        output_dimensions, output_is_scalar = self._selection_output_dimensions(selection_key)
        block_data = self._selection_block_data(selection_key)
        ordered_single_terms = self._ordered_single_terms()
        term_block_selections, term_block_selection_keys = self._term_block_selections()

        n_p = len(ordered_single_terms)
        decomposed_pauli_strings = []
        scalar_offset = 0.0 + 0.0j
        ordered_weights = self.ordered_symmetrised_hamiltonian.weights
        ordered_phases = self.ordered_symmetrised_hamiltonian.phases
        previous_lcm = self.ordered_symmetrised_hamiltonian.lcm
        for i in range(n_p):
            # weight component
            weight = ordered_weights[i]
            # reconstruct phase component later
            phase = ordered_phases[i]

            # pauli strings of symmetry blocks that are later tensored
            conditional_pauli_string_blocks = []
            term_is_zero = False
            # for each symmetry block decompose the pauli string into the eigenbasis
            for block_idx, (d, eigenvectors, block_cache_key) in enumerate(block_data):
                ps_selection = term_block_selections[block_idx][i]
                ps_selection_key = term_block_selection_keys[block_idx][i]

                block_cache = self.P_qudit_dict.setdefault(block_cache_key, {})
                if ps_selection_key in block_cache:
                    P_qudit = block_cache[ps_selection_key]
                else:
                    P_qudit = self.selection_to_qudit(ps_selection, d, eigenvectors)
                    block_cache[ps_selection_key] = P_qudit

                # selection_to_qudit may legitimately return a scalar (e.g., d==1 or zero block).
                if np.isscalar(P_qudit):
                    p_scalar = P_qudit
                    if np.abs(p_scalar) < 1e-12:
                        term_is_zero = True
                        break
                    weight *= p_scalar
                else:
                    conditional_pauli_string_blocks.append(P_qudit)
            if term_is_zero:
                continue

            # tensor the conditional pauli strings together
            if len(conditional_pauli_string_blocks) == 0:
                scalar_offset += weight * np.exp(phase * 2 * np.pi * 1j / (2 * previous_lcm))
            else:
                new_ps = conditional_pauli_string_blocks[0].copy()
                for i in range(1, len(conditional_pauli_string_blocks)):
                    new_ps = new_ps @ conditional_pauli_string_blocks[i]

                new_ps.weights *= weight
                new_ps.phases += int(phase * new_ps.lcm / previous_lcm)
                decomposed_pauli_strings.append(new_ps)

        if output_is_scalar:
            return scalar_offset

        new_P = None
        if len(decomposed_pauli_strings) > 0:
            new_P = decomposed_pauli_strings[0]
            for i in range(1, len(decomposed_pauli_strings)):
                new_P = new_P + decomposed_pauli_strings[i]

        if np.abs(scalar_offset) > 1e-12:
            x0 = np.zeros(len(output_dimensions), dtype=int)
            z0 = np.zeros(len(output_dimensions), dtype=int)
            identity_ps = PauliString.from_exponents(x0, z0, output_dimensions)
            scalar_ps = PauliSum.from_pauli_strings(identity_ps, weights=[scalar_offset], phases=[0])
            new_P = scalar_ps if new_P is None else new_P + scalar_ps

        if new_P is None:
            return _zero_paulisum(output_dimensions)

        new_P.combine_equivalent_paulis()
        new_P.remove_zero_weight_paulis()
        if new_P.n_paulis() == 0:
            return _zero_paulisum(output_dimensions)
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
        if self._inverse_sector_permutation is None:
            self._initialize_sector_index_maps()
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
        # full_hamiltonian_matrix is sector-ordered; map it into kron-order so it aligns
        # with U from diagonalising_unitary(inverse=True) before reconstructing.
        if full_hamiltonian_matrix.shape[0] > 0:
            inv_perm = self._inverse_sector_permutation
            full_hamiltonian_matrix = full_hamiltonian_matrix[np.ix_(inv_perm, inv_perm)]
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


class SwapConditionalHamiltonian(ConditionalHamiltonian):
    def __init__(
        self,
        hamiltonian: PauliSum,
        symmetry: Circuit,
        transformation_gate: Gate | Circuit,
        verify: bool = True,
        degeneracy_tol: float = 1e-8,
        zero_tol: float = 1e-12,
    ):
        if not isinstance(symmetry, Circuit):
            raise TypeError("SwapConditionalHamiltonian requires a Circuit of SWAP gates.")

        self.original_hamiltonian = hamiltonian
        self.verify = bool(verify)
        self.degeneracy_tol = float(degeneracy_tol)
        self.zero_tol = float(zero_tol)
        self.input_dimensions = np.asarray(hamiltonian.dimensions, dtype=int)
        self.symmetry_circuit = symmetry
        self.symmetry_gate = symmetry.composite_gate()
        self.transformation_gate, self.transformation_circuit = self._coerce_gate_and_circuit(transformation_gate)
        self.symmetrised_hamiltonian = self.transformation_circuit.inverse().act(self.original_hamiltonian)

        if self.verify:
            assert test_symmetry(self.symmetry_circuit, self.symmetrised_hamiltonian)

        source_blocks = self._swap_source_blocks()
        new_qudit_order = [qudit for block in source_blocks for qudit in block]
        swap_list = permutation_swaps(new_qudit_order)
        swap_circuit = Circuit(
            dimensions=self.input_dimensions,
            gates=[SWAP() for _ in swap_list],
            qudit_indices=[(s[0], s[1]) for s in swap_list],
        )
        inv_swap_circuit = swap_circuit.inverse()
        self.ordered_symmetry_circuit = inv_swap_circuit + self.symmetry_circuit + swap_circuit
        self.ordered_symmetry_gate = self.ordered_symmetry_circuit.composite_gate()
        self.ordered_T_circuit = inv_swap_circuit + self.transformation_circuit
        self.ordered_T_gate = self.ordered_T_circuit.composite_gate()
        self.ordered_symmetrised_hamiltonian = swap_circuit.act(self.symmetrised_hamiltonian)

        if self.verify:
            assert test_symmetry(self.ordered_symmetry_circuit, self.ordered_symmetrised_hamiltonian)

        self.blocks = []
        qudit_counter = 0
        for block in source_blocks:
            relabelled_block = []
            for _ in block:
                relabelled_block.append(qudit_counter)
                qudit_counter += 1
            self.blocks.append(relabelled_block)

        self.block_variables = {}
        self.block_variables['unitaries'] = []
        self.block_variables['eigenvalues'] = []
        self.block_variables['eigenvectors'] = []
        self.block_variables['degeneracies'] = {}
        self.block_variables['degeneracy_dimensions'] = {}
        self.block_variables['conditional_dimensions'] = []
        self.P_qudit_dict = {}
        self._block_qudits_to_delete = [
            np.array(
                [q for q in range(self.ordered_symmetrised_hamiltonian.n_qudits()) if q not in block],
                dtype=int,
            )
            for block in self.blocks
        ]
        self._ordered_single_terms_cache = None
        self._term_block_selections_cache = None
        self._term_block_selection_keys_cache = None

        for block_index, block in enumerate(self.blocks):
            block_dimensions = np.asarray(self.ordered_symmetrised_hamiltonian.dimensions[block], dtype=int)
            if len(block) == 2:
                eigenvalues, eigenvectors, unitary, degeneracies = self._swap_block_eigendecomposition(
                    int(block_dimensions[0])
                )
            elif len(block) == 1:
                d = int(block_dimensions[0])
                eigenvalues = np.ones(d, dtype=complex)
                eigenvectors = np.eye(d, dtype=complex)
                unitary = np.eye(d, dtype=complex)
                degeneracies = [list(range(d))]
            else:
                raise ValueError(f"Unexpected SWAP source block size {len(block)}.")

            self.block_variables['unitaries'].append(unitary)
            self.block_variables['eigenvalues'].append(eigenvalues)
            self.block_variables['eigenvectors'].append(eigenvectors)
            self.block_variables['degeneracies'][block_index] = degeneracies
            self.block_variables['degeneracy_dimensions'][block_index] = [len(g) for g in degeneracies]
            self.block_variables['conditional_dimensions'].append(len(degeneracies))

        self._sector_basis_indices = None
        self._sector_permutation = None
        self._inverse_sector_permutation = None

    def _swap_source_blocks(self) -> list[list[int]]:
        swap_blocks: list[list[int]] = []
        occupied: set[int] = set()
        for gate, qudits in zip(self.symmetry_circuit.gates, self.symmetry_circuit.qudit_indices):
            if getattr(gate, "name", None) != "SWAP":
                raise ValueError("SwapConditionalHamiltonian only supports SWAP gates.")
            if len(qudits) != 2:
                raise ValueError(f"SWAP gate must have two qudit indices, got {qudits}.")
            a, b = int(qudits[0]), int(qudits[1])
            if a == b:
                raise ValueError("SWAP gate cannot act on the same qudit twice.")
            if a in occupied or b in occupied:
                raise ValueError("SwapConditionalHamiltonian requires independent SWAP gates.")
            if int(self.input_dimensions[a]) != int(self.input_dimensions[b]):
                raise ValueError("SWAP fast path requires equal dimensions within each pair.")
            occupied.add(a)
            occupied.add(b)
            swap_blocks.append(sorted([a, b]))

        blocks = [list(block) for block in swap_blocks]
        for idx in range(int(self.original_hamiltonian.n_qudits())):
            if idx not in occupied:
                blocks.append([idx])
        blocks.sort(key=lambda block: block[0])
        return blocks

    @staticmethod
    def _swap_block_eigendecomposition(dimension: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[list[int]]]:
        d = int(dimension)
        hilbert_dim = d * d
        columns: list[np.ndarray] = []
        eigenvalues: list[complex] = []

        for j in range(d):
            for k in range(j + 1, d):
                vec = np.zeros(hilbert_dim, dtype=complex)
                vec[j * d + k] = 1.0 / np.sqrt(2.0)
                vec[k * d + j] = -1.0 / np.sqrt(2.0)
                columns.append(vec)
                eigenvalues.append(-1.0 + 0.0j)

        for j in range(d):
            vec = np.zeros(hilbert_dim, dtype=complex)
            vec[j * d + j] = 1.0 + 0.0j
            columns.append(vec)
            eigenvalues.append(1.0 + 0.0j)

        for j in range(d):
            for k in range(j + 1, d):
                vec = np.zeros(hilbert_dim, dtype=complex)
                vec[j * d + k] = 1.0 / np.sqrt(2.0)
                vec[k * d + j] = 1.0 / np.sqrt(2.0)
                columns.append(vec)
                eigenvalues.append(1.0 + 0.0j)

        eigenvectors = np.column_stack(columns).astype(complex)
        eigenvalues_arr = np.asarray(eigenvalues, dtype=complex)
        unitary = np.zeros((hilbert_dim, hilbert_dim), dtype=complex)
        for j in range(d):
            for k in range(d):
                unitary[k * d + j, j * d + k] = 1.0 + 0.0j

        antisymmetric_dim = d * (d - 1) // 2
        degeneracies: list[list[int]] = []
        if antisymmetric_dim > 0:
            degeneracies.append(list(range(antisymmetric_dim)))
        degeneracies.append(list(range(antisymmetric_dim, hilbert_dim)))
        return eigenvalues_arr, eigenvectors, unitary, degeneracies

    def _ensure_sector_index_maps(self) -> None:
        if self._sector_basis_indices is None:
            self._initialize_sector_index_maps()

    def sector_basis_indices(self) -> list[np.ndarray]:
        self._ensure_sector_index_maps()
        return [np.asarray(idx, dtype=int).copy() for idx in self._sector_basis_indices]

    def reconstruct_full_hamiltonian(self):
        self._ensure_sector_index_maps()
        return super().reconstruct_full_hamiltonian()


class ConditionalHamiltonian2:
    def __init__(
        self,
        hamiltonian: PauliSum,
        symmetry: Gate | Circuit,
        transformation_gate: Gate | Circuit,
        verify: bool = True,
        degeneracy_tol: float = 1e-8,
        zero_tol: float = 1e-12,
    ):
        # original matrices
        self.original_hamiltonian = hamiltonian
        self.input_dimensions = np.asarray(hamiltonian.dimensions, dtype=int)
        self.input_qudits = hamiltonian.n_qudits()

        self.symmetry_gate, self.symmetry_circuit = self._coerce_gate_and_circuit(symmetry)
        self.transformation_gate, self.transformation_circuit = self._coerce_gate_and_circuit(transformation_gate)
        self.symmetrised_hamiltonian = self.transformation_circuit.inverse().act(self.original_hamiltonian)

        # general flags and tolerances
        self.verify = bool(verify)
        self.degeneracy_tol = float(degeneracy_tol)
        self.zero_tol = float(zero_tol)

        # check symmetry
        if self.verify:
            assert test_symmetry(self.symmetry_circuit, self.symmetrised_hamiltonian)

        self.blocks = block_indexes(self.symmetry_gate.symplectic)
        self.block_vars = {}
        self.block_vars['eigenvalues'] = []
        self.block_vars['eigenvectors'] = []
        self.block_vars['circuits'] = []
        self.block_vars['unitaries'] = []

        commuting_blocks, non_commuting_blocks = self._find_commuting_blocks()

        # find composite commuting blocks
        commuting_blocks = self._find_composite_commuting_blocks(non_commuting_blocks, commuting_blocks)

        # turn list of commuting block indices into list of blocks (with qudit indices)
        source_blocks = self._construct_source_blocks(commuting_blocks)

        # reorder such that qudits in commuting blocks are next to each other
        self._reorder_operators(source_blocks)

        # relabel blocks
        self._construct_symmetry_group(source_blocks, commuting_blocks)

        # assign circuits, eigenvalues and eigenvectors to each block
        self.block_origins = []
        self.symmetry_group_vars = {}
        self.symmetry_group_vars['eigenvalues'] = []
        self.symmetry_group_vars['eigenvectors'] = []
        self.symmetry_group_vars['degeneracies'] = {}
        self.symmetry_group_vars['degeneracy_dimensions'] = {}
        self.symmetry_group_vars['conditional_dimensions'] = []
        self.symmetry_group_vars['eigenvector_order'] = []

        self.block_projectors = {}
        self.block_projector_bases = {}
        self.sector_projectors = {}
        self._block_qudits_to_delete = []
        self._ordered_single_terms_cache = None
        self._term_block_selections_cache = None
        self._term_block_selection_keys_cache = None

        for block_index, block in enumerate(self.symmetry_groups):

            # block parameters
            block_dimensions = np.asarray(self.ordered_symmetrised_hamiltonian.dimensions[block], dtype=int)
            hilbert_space_size = int(np.prod(block_dimensions))

            # list where entries are indices of original blocks
            block_origin = self._find_block_origin(block_index, commuting_blocks)
            self.block_origins.append(block_origin)

            # collect unitaries
            unitary_list = [self.block_vars['unitaries'][origin] for origin in block_origin]
            block_unitary = multi_kron(unitary_list)

            # collect eigenvalues for component blocks
            ev, evec, icm, dg = self._find_symmetry_group_eigen_data(block_origin, hilbert_space_size, block_unitary)

            # save eigenvalues and eigenvectors
            degeneracy_indeces = [list(g) for g in dg]
            d_dimension = [len(g) for g in dg]
            cond_dimension = len(d_dimension)
            self.symmetry_group_vars['eigenvector_order'].append(icm)
            self.symmetry_group_vars['eigenvalues'].append(ev)
            self.symmetry_group_vars['eigenvectors'].append(evec)
            self.symmetry_group_vars['degeneracies'][block_index] = degeneracy_indeces
            self.symmetry_group_vars['degeneracy_dimensions'][block_index] = d_dimension
            self.symmetry_group_vars['conditional_dimensions'].append(cond_dimension)

            self._precompute_block_dictionaries(degeneracy_indeces, block_origin, block_index, evec)

        self._block_qudits_to_delete = [
            np.array(
                [q for q in range(self.ordered_symmetrised_hamiltonian.n_qudits()) if q not in block],
                dtype=int,
            )
            for block in self.symmetry_groups
        ]
        self._precompute_actual_block_projectors()

        self._sector_basis_indices = None
        self._sector_permutation = None
        self._inverse_sector_permutation = None

    def _coerce_gate_and_circuit(self, operator: Gate | Circuit) -> tuple[Gate, Circuit]:
        if isinstance(operator, Circuit):
            return operator.composite_gate(), operator
        return operator, gate_to_circuit(operator, dimensions=self.input_dimensions)

    def _find_commuting_blocks(self):
        commuting_blocks = []
        non_commuting_blocks = []

        # find commuting and non commuting blocks as well as their circuits and eigenvalues
        for block_index, block in enumerate(self.blocks):
            # circuits
            block_circuit = self.symmetry_circuit.local_circuit(block)
            self.block_vars['circuits'].append(block_circuit)

            # unitaries
            block_unitary = np.around(block_circuit.unitary().toarray(), 14)
            self.block_vars['unitaries'].append(block_unitary)

            # eigenvalues and eigenvectors
            if np.isclose(block_unitary, np.eye(block_unitary.shape[0], dtype=complex)).all():
                block_hilbert_space_size = int(block_unitary.shape[0])
                eigenvalues = np.ones(block_hilbert_space_size, dtype=complex)
                eigenvectors = np.eye(block_hilbert_space_size, dtype=complex)
            else:
                eigenvalues, eigenvectors = np.linalg.eig(block_unitary)

            if self.verify:
                assert np.all(np.around(eigenvectors.transpose().conj() @ block_unitary @
                              eigenvectors - np.diag(eigenvalues), 10) == 0)

            # sorting eigenvalues and eigenvectors
            idx = np.argsort(eigenvalues.real)
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]

            # saving eigenvalues and eigenvectors
            self.block_vars['eigenvalues'].append(eigenvalues)
            self.block_vars['eigenvectors'].append(eigenvectors)

            # find commuting blocks
            embedded_block_C = self.embedded_block_circuit([block_index], self.blocks, self.block_vars)
            if test_symmetry(embedded_block_C, self.symmetrised_hamiltonian):
                commuting_blocks.append(block_index)
            else:
                non_commuting_blocks.append(block_index)
        return commuting_blocks, non_commuting_blocks

    def _find_composite_commuting_blocks(self, non_commuting_blocks, commuting_blocks):
        while len(non_commuting_blocks) > 0:
            composite_blocks = self._find_composite_partners(non_commuting_blocks)
            commuting_blocks.insert(0, composite_blocks)
            for b in composite_blocks:
                non_commuting_blocks.remove(b)
        return commuting_blocks

    def _find_composite_partners(self, non_commuting_blocks):
        initial_index = non_commuting_blocks[0]
        for r in range(1, len(non_commuting_blocks[1:]) + 1):  # subset size
            for subset in combinations(non_commuting_blocks[1:], r):
                embedded_block_C = self.embedded_block_circuit([initial_index] + list(subset),
                                                               self.blocks, self.block_vars)
                if test_symmetry(embedded_block_C, self.symmetrised_hamiltonian):
                    return [initial_index] + list(subset)

        raise Exception('No commuting blocks found')

    def _construct_source_blocks(self, commuting_blocks):
        source_blocks = []
        for block_index in commuting_blocks:
            if isinstance(block_index, int):
                block = self.blocks[block_index]
                source_blocks.append(block)
            elif isinstance(block_index, list):
                block = []
                for b in block_index:
                    block += self.blocks[b]
                source_blocks.append(block)
            else:
                raise Exception('Block index must be an int or list of ints')
        return source_blocks

    def _reorder_operators(self, source_blocks):
        new_qudit_order = [qudit for block in source_blocks for qudit in block]
        swap_list = permutation_swaps(new_qudit_order)
        swap_circuit = Circuit(dimensions=self.input_dimensions,
                               gates=[SWAP() for _ in swap_list],
                               qudit_indices=[(s[0], s[1]) for s in swap_list])
        inv_swap_circuit = swap_circuit.inverse()
        self.ordered_symmetry_circuit = inv_swap_circuit + self.symmetry_circuit.copy() + swap_circuit
        self.ordered_symmetry_gate = self.ordered_symmetry_circuit.composite_gate()
        self.ordered_T_circuit = inv_swap_circuit + self.transformation_circuit.copy()
        self.ordered_T_gate = self.ordered_T_circuit.composite_gate()
        self.ordered_symmetrised_hamiltonian = swap_circuit.act(self.symmetrised_hamiltonian)
        if self.verify:
            assert test_symmetry(self.ordered_symmetry_gate, self.ordered_symmetrised_hamiltonian)
            assert self.ordered_T_gate.inverse().act(self.original_hamiltonian,
                                                     qudits=tuple(np.arange(self.input_qudits))
                                                     ).is_close(self.ordered_symmetrised_hamiltonian, literal=False)
            assert test_symmetry(self.ordered_symmetry_gate, self.ordered_T_gate.inverse().act(
                self.original_hamiltonian, qudits=tuple(np.arange(self.input_qudits))))

    def _construct_symmetry_group(self, source_blocks, commuting_blocks):
        self.symmetry_groups = []
        qudit_counter = 0
        for block in source_blocks:
            temp_block = []
            for b in block:
                temp_block.append(qudit_counter)
                qudit_counter += 1
            self.symmetry_groups.append(temp_block)
        if self.verify:
            assert len(self.symmetry_groups) == len(commuting_blocks)

    def _find_block_origin(self, block_index, commuting_blocks):
        if isinstance(commuting_blocks[block_index], int):
            block_origin = [commuting_blocks[block_index]]
        elif isinstance(commuting_blocks[block_index], list):
            block_origin = commuting_blocks[block_index]
        else:
            raise Exception('Block origin must be an int or list of ints')
        return block_origin

    def _find_symmetry_group_eigen_data(self, block_origin, hilbert_space_size, block_unitary):
        list_of_indices = []
        list_of_eigenvalues = []
        list_of_eigenvectors = []
        for origin in block_origin:
            list_of_eigenvalues.append(self.block_vars['eigenvalues'][origin])
            list_of_eigenvectors.append(self.block_vars['eigenvectors'][origin])
            list_of_indices.append(np.arange(len(self.block_vars['eigenvalues'][origin])))

        index_combinations = list(itertools.product(*list_of_indices))
        index_combination_matrix = np.zeros((len(index_combinations[0]), len(index_combinations)), dtype=int)
        eigenvalues = np.zeros(len(index_combinations), dtype=complex)
        eigenvectors = np.zeros((hilbert_space_size, hilbert_space_size), dtype=complex)
        for comb_index, combination in enumerate(index_combinations):
            eigenvalues[comb_index] = np.prod([list_of_eigenvalues[ic][comb]
                                               for ic, comb in enumerate(combination)])
            eigenvectors[:, comb_index] = multi_kron([list_of_eigenvectors[ic][:, comb]
                                                      for ic, comb in enumerate(combination)])
            index_combination_matrix[:, comb_index] = combination

        # sort eigenvalues and eigenvectors
        idx = np.argsort(eigenvalues.real)
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        index_combination_matrix = index_combination_matrix[:, idx]

        degeneracy_groups = group_indices_by_tolerance(eigenvalues, self.degeneracy_tol)
        if self.verify:
            assert np.all(np.around(eigenvectors.T.conj() @ eigenvectors - np.eye(len(eigenvectors)), 10) == 0)
            assert np.all(np.around(eigenvectors.transpose().conj() @ block_unitary @
                                    eigenvectors - np.diag(eigenvalues), 10) == 0)
        return eigenvalues, eigenvectors, index_combination_matrix, degeneracy_groups

    def _precompute_block_dictionaries(self, degeneracy_indeces, block_origin, block_index, eigenvectors):
        for i, inds in enumerate(degeneracy_indeces):
            self.sector_projectors[(int(block_index), int(i))] = {}
            for o, origin in enumerate(block_origin):
                evs = self.block_vars['eigenvectors'][origin]
                id_comb_m = self.symmetry_group_vars['eigenvector_order'][block_index]
                block_eigenvectors = [evs[:, id_comb_m[o, ii]] for ii in inds]
                V = np.column_stack(block_eigenvectors).astype(complex)
                dims = self.original_hamiltonian.dimensions[self.blocks[origin]]
                key = (int(block_index), int(i), int(origin))
                self.block_projectors[key] = {}
                self.block_projector_bases[key] = (np.asarray(dims, dtype=int), V)

    def _project_local_pauli(self, pauli_string: PauliString, block_projector_key: tuple[int, int, int]) -> np.ndarray:
        dims, eigenvectors = self.block_projector_bases[block_projector_key]
        partial_hamiltonian_matrix = pauli_string.to_hilbert_space()
        if hasattr(partial_hamiltonian_matrix, "toarray"):
            partial_hamiltonian_matrix = partial_hamiltonian_matrix.toarray()
        partial_hamiltonian_matrix = np.asarray(partial_hamiltonian_matrix, dtype=complex)
        return eigenvectors.conj().T @ partial_hamiltonian_matrix @ eigenvectors

    def _split_selection_by_origin(self, selection: PauliString, block: int) -> list[tuple[int, PauliString, str]]:
        pieces: list[tuple[int, PauliString, str]] = []
        start = 0
        for origin in self.block_origins[block]:
            n_origin_qudits = len(self.blocks[origin])
            stop = start + n_origin_qudits
            dims = np.asarray(selection.dimensions[start:stop], dtype=int)
            x_exp = np.asarray(selection.x_exp[start:stop], dtype=int)
            z_exp = np.asarray(selection.z_exp[start:stop], dtype=int)
            local_ps = PauliString.from_exponents(x_exp, z_exp, dims)
            pieces.append((int(origin), local_ps, str(local_ps).strip()))
            start = stop
        return pieces

    def _precompute_actual_block_projectors(self) -> None:
        term_block_selections, _ = self._term_block_selections()
        for block in range(len(self.symmetry_groups)):
            for ps_selection in term_block_selections[block]:
                local_pieces = self._split_selection_by_origin(ps_selection, block)
                for selected_eigenvalue in range(int(self.symmetry_group_vars['conditional_dimensions'][block])):
                    for origin, local_ps, local_key in local_pieces:
                        projector_key = (int(block), int(selected_eigenvalue), int(origin))
                        projector_cache = self.block_projectors[projector_key]
                        if local_key not in projector_cache:
                            projector_cache[local_key] = self._project_local_pauli(local_ps, projector_key)

    def _initialize_sector_index_maps(self) -> None:
        block_hilbert_dimensions = [
            int(len(self.symmetry_group_vars['eigenvalues'][block]))
            for block in range(len(self.symmetry_groups))
        ]
        conditional_dimensions = [
            int(x) for x in self.symmetry_group_vars.get('conditional_dimensions', [])
        ]
        n_sectors = int(np.prod(conditional_dimensions)) if len(conditional_dimensions) > 0 else 1

        sector_basis_indices: list[np.ndarray] = []
        for sector_idx in range(n_sectors):
            selection = [int(x) for x in int_to_bases(sector_idx, conditional_dimensions)]
            selected_groups = [
                [int(i) for i in self.symmetry_group_vars['degeneracies'][block][selection[block]]]
                for block in range(len(self.symmetry_groups))
            ]

            local_indices: list[int] = []
            for basis_tuple in itertools.product(*selected_groups):
                local_indices.append(int(bases_to_int(np.asarray(basis_tuple, dtype=int), block_hilbert_dimensions)))
            sector_basis_indices.append(np.asarray(local_indices, dtype=int))

        if len(sector_basis_indices) == 0:
            sector_permutation = np.array([], dtype=int)
        else:
            sector_permutation = np.concatenate(sector_basis_indices).astype(int)

        hilbert_dim = int(np.prod(block_hilbert_dimensions))
        if sector_permutation.size != hilbert_dim:
            raise RuntimeError(
                "Internal sector-index map size mismatch: "
                f"got {sector_permutation.size}, expected {hilbert_dim}."
            )
        if hilbert_dim > 0 and not np.array_equal(np.sort(sector_permutation), np.arange(hilbert_dim, dtype=int)):
            raise RuntimeError("Internal sector-index map is not a valid permutation.")

        self._sector_basis_indices = sector_basis_indices
        # sector_permutation[k] gives the kron-order basis index at sector-order position k.
        self._sector_permutation = sector_permutation
        # inverse map: kron-order -> sector-order.
        self._inverse_sector_permutation = np.argsort(sector_permutation).astype(int)

    def sector_basis_indices(self) -> list[np.ndarray]:
        if self._sector_basis_indices is None:
            self._initialize_sector_index_maps()
        return [np.asarray(idx, dtype=int).copy() for idx in self._sector_basis_indices]

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

    def selection_to_qudit(self, selection: PauliString, qudit_dimension: int,
                           block: int, selected_eigenvalue: int):
        block_G = []

        for origin, local_ps, local_key in self._split_selection_by_origin(selection, block):
            projector_key = (int(block), int(selected_eigenvalue), int(origin))
            projector_cache = self.block_projectors[projector_key]
            if local_key not in projector_cache:
                projector_cache[local_key] = self._project_local_pauli(local_ps, projector_key)
            G = projector_cache[local_key]
            block_G.append(G)

        G = np.prod(block_G, axis=0)

        # Turn projector matrix into pauli sum
        if qudit_dimension == 1:
            return G[0, 0]
        else:
            # Guard: a numerically zero matrix cannot be expanded into a non-empty PauliSum.
            if float(np.max(np.abs(G))) < self.zero_tol:
                return 0.0 + 0.0j
            P_qudit = PauliSum.from_hilbert_space(G, prime_factors(qudit_dimension)[::-1])
            return P_qudit

    def selection_to_matrix(self, selection: PauliString, qudit_dimension: int,
                            block: int, selected_eigenvalue: int):
        block_G = []
        for origin, local_ps, local_key in self._split_selection_by_origin(selection, block):
            projector_key = (int(block), int(selected_eigenvalue), int(origin))
            projector_cache = self.block_projectors[projector_key]
            if local_key not in projector_cache:
                projector_cache[local_key] = self._project_local_pauli(local_ps, projector_key)
            block_G.append(projector_cache[local_key])

        G = np.prod(block_G, axis=0)
        if qudit_dimension == 1:
            return G[0, 0]
        if float(np.max(np.abs(G))) < self.zero_tol:
            return 0.0 + 0.0j
        return np.asarray(G, dtype=complex)

    def _selection_key(self, selection: list[int] | np.ndarray | tuple[int, ...]) -> tuple[int, ...]:
        return tuple(int(x) for x in selection)

    def _selection_output_dimensions(self, selection_key: tuple[int, ...]) -> tuple[list[int], bool]:
        output_dimensions: list[int] = []
        for block, selected_eigenvalue in enumerate(selection_key):
            d = int(self.symmetry_group_vars['degeneracy_dimensions'][block][selected_eigenvalue])
            if d > 1:
                output_dimensions += [int(x) for x in prime_factors(d)[::-1]]
        return output_dimensions, len(output_dimensions) == 0

    def _selection_block_data(self, selection_key: tuple[int, ...]):
        block_data = []
        for block, selected_eigenvalue in enumerate(selection_key):
            block_data.append(
                (
                    int(block),
                    int(selected_eigenvalue),
                    int(self.symmetry_group_vars['degeneracy_dimensions'][block][selected_eigenvalue]),
                    (int(block), int(selected_eigenvalue)),
                )
            )

        return block_data

    def _pauli_string_cache_key(self, pauli_string: PauliString) -> tuple[tuple[int, ...], tuple[int, ...]]:
        return (
            tuple(int(x) for x in np.asarray(pauli_string.tableau, dtype=int).reshape(-1)),
            tuple(int(x) for x in np.asarray(pauli_string.dimensions, dtype=int).reshape(-1)),
        )

    def _ordered_single_terms(self) -> list[PauliString]:
        if self._ordered_single_terms_cache is None:
            ordered_single_terms: list[PauliString] = []
            for i in range(self.ordered_symmetrised_hamiltonian.n_paulis()):
                ps = self.ordered_symmetrised_hamiltonian[i].copy()
                ps.weights = np.array([1])
                ps.phases = np.array([0])
                ordered_single_terms.append(ps)
            self._ordered_single_terms_cache = ordered_single_terms
        return self._ordered_single_terms_cache

    def _term_block_selections(
        self,
    ) -> tuple[list[list[PauliString]], list[list[tuple[tuple[int, ...], tuple[int, ...]]]]]:
        if self._term_block_selections_cache is None or self._term_block_selection_keys_cache is None:
            ordered_single_terms = self._ordered_single_terms()
            term_block_selections: list[list[PauliString]] = []
            term_block_selection_keys: list[list[tuple[tuple[int, ...], tuple[int, ...]]]] = []
            for qudits_to_delete in self._block_qudits_to_delete:
                block_terms: list[PauliString] = []
                block_keys: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
                for ps in ordered_single_terms:
                    ps_selection = ps._delete_qudits(qudits_to_delete)
                    block_terms.append(ps_selection)
                    block_keys.append(self._pauli_string_cache_key(ps_selection))
                term_block_selections.append(block_terms)
                term_block_selection_keys.append(block_keys)
            self._term_block_selections_cache = term_block_selections
            self._term_block_selection_keys_cache = term_block_selection_keys
        return self._term_block_selections_cache, self._term_block_selection_keys_cache

    def select_hamiltonian(self, selection: list[int]) -> PauliSum | complex:
        selection_key = self._selection_key(selection)

        def _zero_paulisum(dimensions: list[int]) -> PauliSum:
            x0 = np.zeros(len(dimensions), dtype=int)
            z0 = np.zeros(len(dimensions), dtype=int)
            identity_ps = PauliString.from_exponents(x0, z0, dimensions)
            return PauliSum.from_pauli_strings(identity_ps, weights=[0.0 + 0.0j], phases=[0])

        output_dimensions, output_is_scalar = self._selection_output_dimensions(selection_key)
        block_data = self._selection_block_data(selection_key)
        ordered_single_terms = self._ordered_single_terms()
        term_block_selections, term_block_selection_keys = self._term_block_selections()

        n_p = len(ordered_single_terms)
        decomposed_pauli_strings = []
        scalar_offset = 0.0 + 0.0j
        ordered_weights = self.ordered_symmetrised_hamiltonian.weights
        ordered_phases = self.ordered_symmetrised_hamiltonian.phases
        previous_lcm = self.ordered_symmetrised_hamiltonian.lcm
        for i in range(n_p):
            # weight component
            weight = ordered_weights[i]
            # reconstruct phase component later
            phase = ordered_phases[i]

            # pauli strings of symmetry blocks that are later tensored
            conditional_pauli_string_blocks = []
            term_is_zero = False
            # for each symmetry block decompose the pauli string into the eigenbasis
            for block_idx, (block, selected_eigenvalue, d, block_cache_key) in enumerate(block_data):
                ps_selection = term_block_selections[block_idx][i]
                ps_selection_key = term_block_selection_keys[block_idx][i]
                block_cache = self.sector_projectors.setdefault(block_cache_key, {})
                if ps_selection_key in block_cache:
                    P_qudit = block_cache[ps_selection_key]
                else:
                    P_qudit = self.selection_to_qudit(ps_selection, d, block, selected_eigenvalue)
                    block_cache[ps_selection_key] = P_qudit

                # selection_to_qudit may legitimately return a scalar (e.g., d==1 or zero block).
                if np.isscalar(P_qudit):
                    p_scalar = P_qudit
                    if np.abs(p_scalar) < 1e-12:
                        term_is_zero = True
                        break
                    weight *= p_scalar
                else:
                    conditional_pauli_string_blocks.append(P_qudit)
            if term_is_zero:
                continue
            # tensor the conditional pauli strings together
            if len(conditional_pauli_string_blocks) == 0:
                scalar_offset += weight * np.exp(phase * 2 * np.pi * 1j / (2 * previous_lcm))
            else:
                new_ps = conditional_pauli_string_blocks[0].copy()
                for i in range(1, len(conditional_pauli_string_blocks)):
                    new_ps = new_ps @ conditional_pauli_string_blocks[i]

                new_ps.weights *= weight
                new_ps.phases += int(phase * new_ps.lcm / previous_lcm)
                decomposed_pauli_strings.append(new_ps)

        if output_is_scalar:
            return scalar_offset

        new_P = None
        if len(decomposed_pauli_strings) > 0:
            new_P = decomposed_pauli_strings[0]
            for i in range(1, len(decomposed_pauli_strings)):
                new_P = new_P + decomposed_pauli_strings[i]

        if np.abs(scalar_offset) > 1e-12:
            x0 = np.zeros(len(output_dimensions), dtype=int)
            z0 = np.zeros(len(output_dimensions), dtype=int)
            identity_ps = PauliString.from_exponents(x0, z0, output_dimensions)
            scalar_ps = PauliSum.from_pauli_strings(identity_ps, weights=[scalar_offset], phases=[0])
            new_P = scalar_ps if new_P is None else new_P + scalar_ps

        if new_P is None:
            return _zero_paulisum(output_dimensions)

        new_P.combine_equivalent_paulis()
        new_P.remove_zero_weight_paulis()
        if new_P.n_paulis() == 0:
            return _zero_paulisum(output_dimensions)
        return new_P

    def select_hamiltonian_matrix(self, selection: list[int], *, sparse: bool = False):
        selection_key = self._selection_key(selection)
        output_dimensions, output_is_scalar = self._selection_output_dimensions(selection_key)
        block_data = self._selection_block_data(selection_key)
        ordered_single_terms = self._ordered_single_terms()
        term_block_selections, _ = self._term_block_selections()

        output_dim = int(np.prod(output_dimensions)) if len(output_dimensions) > 0 else 1
        H = (
            sp.csr_matrix((output_dim, output_dim), dtype=complex)
            if sparse
            else np.zeros((output_dim, output_dim), dtype=complex)
        )
        scalar_offset = 0.0 + 0.0j
        ordered_weights = self.ordered_symmetrised_hamiltonian.weights
        ordered_phases = self.ordered_symmetrised_hamiltonian.phases
        previous_lcm = self.ordered_symmetrised_hamiltonian.lcm

        for term_idx in range(len(ordered_single_terms)):
            weight = ordered_weights[term_idx]
            phase = ordered_phases[term_idx]
            term_blocks = []
            term_is_zero = False

            for block_idx, (block, selected_eigenvalue, d, _) in enumerate(block_data):
                ps_selection = term_block_selections[block_idx][term_idx]
                block_matrix = self.selection_to_matrix(ps_selection, d, block, selected_eigenvalue)
                if np.isscalar(block_matrix):
                    block_scalar = complex(block_matrix)
                    if np.abs(block_scalar) < self.zero_tol:
                        term_is_zero = True
                        break
                    weight *= block_scalar
                else:
                    term_blocks.append(block_matrix)

            if term_is_zero:
                continue

            coeff = weight * np.exp(phase * 2 * np.pi * 1j / (2 * previous_lcm))
            if len(term_blocks) == 0:
                scalar_offset += coeff
                if not output_is_scalar:
                    identity = sp.eye(output_dim, format="csr", dtype=complex) if sparse else np.eye(output_dim, dtype=complex)
                    H = H + coeff * identity
                continue

            if sparse:
                term_matrix = sp.csr_matrix(term_blocks[0])
                for block_matrix in term_blocks[1:]:
                    term_matrix = sp.kron(term_matrix, sp.csr_matrix(block_matrix), format="csr")
                H = H + coeff * term_matrix
            else:
                term_matrix = np.asarray(term_blocks[0], dtype=complex)
                for block_matrix in term_blocks[1:]:
                    term_matrix = np.kron(term_matrix, np.asarray(block_matrix, dtype=complex))
                H = H + coeff * term_matrix

        if output_is_scalar:
            return scalar_offset
        return H

    def diagonalising_unitary(self, inverse=False):
        U_blocks = []
        for block in range(len(self.symmetry_groups)):
            eigenvectors = self.symmetry_group_vars['eigenvectors'][block]
            if inverse:
                U_blocks.append(np.linalg.inv(eigenvectors))
            else:
                U_blocks.append(eigenvectors)

        U = multi_kron(U_blocks)
        return U

    def reconstruct_full_hamiltonian(self):
        if self._inverse_sector_permutation is None:
            self._initialize_sector_index_maps()
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
        # full_hamiltonian_matrix is sector-ordered; map it into kron-order so it aligns
        # with U from diagonalising_unitary(inverse=True) before reconstructing.
        if full_hamiltonian_matrix.shape[0] > 0:
            inv_perm = self._inverse_sector_permutation
            full_hamiltonian_matrix = full_hamiltonian_matrix[np.ix_(inv_perm, inv_perm)]
        U = self.diagonalising_unitary(inverse=True)
        reconstructed_H = U.T.conj() @ full_hamiltonian_matrix @ U
        return reconstructed_H

    def test_conditional_hamiltonian(self):
        comparison_hamiltonian = self.ordered_symmetrised_hamiltonian.to_hilbert_space().toarray()
        reconstructed_hamiltonian = self.reconstruct_full_hamiltonian()
        return bool(np.max(np.abs(comparison_hamiltonian - reconstructed_hamiltonian)) < 1e-10)

    def __iter__(self):
        n = np.prod(self.symmetry_group_vars['conditional_dimensions'])
        for i in range(n):
            selection = list(int_to_bases(i, self.symmetry_group_vars['conditional_dimensions']))
            block_hamiltonian = self.select_hamiltonian(selection)
            yield block_hamiltonian
