import numpy as np
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

        # In the Kronecker eigenbasis U = kron(U_block_0, ..., U_block_n), sector subspaces
        # are generally interleaved. Precompute index maps from sector-order <-> kron-order.
        self._initialize_sector_index_maps()

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
    ) -> list[tuple[int, list[np.ndarray], np.ndarray, np.ndarray]]:
        block_data: list[tuple[int, list[np.ndarray], np.ndarray, np.ndarray]] = []
        for block, selected_eigenvalue in enumerate(selection_key):
            relevant_qudits = self.blocks[block]
            relevant_indices = np.array(self.block_variables['degeneracies'][block][selected_eigenvalue])
            eigenvectors = [self.block_variables['eigenvectors'][block][:, r] for r in relevant_indices]
            block_data.append(
                (
                    int(self.block_variables['degeneracy_dimensions'][block][selected_eigenvalue]),
                    eigenvectors,
                    np.column_stack(eigenvectors).astype(complex),
                    np.array(
                        [q for q in range(self.ordered_symmetrised_hamiltonian.n_qudits()) if q not in relevant_qudits],
                        dtype=int,
                    ),
                )
            )
        return block_data

    def select_hamiltonian(self, selection: list[int]) -> PauliSum | complex:
        selection_key = self._selection_key(selection)

        def _zero_paulisum(dimensions: list[int]) -> PauliSum:
            x0 = np.zeros(len(dimensions), dtype=int)
            z0 = np.zeros(len(dimensions), dtype=int)
            identity_ps = PauliString.from_exponents(x0, z0, dimensions)
            return PauliSum.from_pauli_strings(identity_ps, weights=[0.0 + 0.0j], phases=[0])

        output_dimensions, output_is_scalar = self._selection_output_dimensions(selection_key)
        block_data = self._selection_block_data(selection_key)

        n_p = self.original_hamiltonian.n_paulis()
        decomposed_pauli_strings = []
        scalar_offset = 0.0 + 0.0j
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
            term_is_zero = False
            # for each symmetry block decompose the pauli string into the eigenbasis
            for d, eigenvectors, _, qudits_to_delete in block_data:
                # get the dimension of the degeneracy
                ps_selection = ps.copy()
                ps_selection = ps_selection._delete_qudits(qudits_to_delete)

                key = str(np.around(np.array(eigenvectors), 6))
                if key in self.P_qudit_dict.keys():
                    if str(ps_selection) in self.P_qudit_dict[key].keys():
                        P_qudit = self.P_qudit_dict[key][str(ps_selection)]
                    else:
                        P_qudit = self.selection_to_qudit(ps_selection, d, eigenvectors)
                        self.P_qudit_dict[key][str(ps_selection)] = P_qudit
                else:
                    self.P_qudit_dict[key] = {}
                    self.P_qudit_dict[key][str(ps_selection)] = self.selection_to_qudit(ps_selection, d, eigenvectors)

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
                new_ps = conditional_pauli_string_blocks[0]
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
