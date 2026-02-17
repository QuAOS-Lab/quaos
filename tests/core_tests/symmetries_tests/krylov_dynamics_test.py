# import numpy as np
# import pytest

# from sympleq.core.paulis import PauliSum, PauliString
# from sympleq.core.circuits.utils import pauli_unitary_qudit, I_mat
# from sympleq.core.states.state import State
# from sympleq.core.circuits import Circuit
# from sympleq.core.circuits.gates import Gate, PauliGate
# from sympleq.core.circuits.gate_decomposition_to_circuit import gate_to_circuit
# from sympleq.models.Ising import ising_chain_hamiltonian, heuristic_clifford_symmetry, product_state_ising
# from sympleq.core.symmetries.block_decomposition import block_decompose_optimal, block_indexes
# from sympleq.core.symmetries.clifford import clifford_phase_decomposition
# from sympleq.core.paulis.utils import apply_paulisum_to_state_dense
# # Adjust this import to wherever you put the Krylov helpers:
# from sympleq.core.symmetries.dynamics.krylov_time_evolve import (
#     _apply_single_qudit_unitary,
#     apply_pauli_row_to_state,
#     apply_paulisum_to_state,
#     project_to_symmetry_sector,
#     krylov_observable_dynamics_symmetry, krylov_observable_dynamics_plain,
#     project_to_symmetry_sector_statevector
# )
# from sympleq.core.symmetries.dynamics.time_evolution_examples import observable_dynamics_full_hilbert
# from sympleq.core.symmetries.dynamics.symmetry_reduced_krylov import (
#     build_symmetry_sector_basis,
# )


# class TestKrylovDynamics:
#     def test_apply_single_qudit_unitary_matches_tensor(self):
#         """
#         _apply_single_qudit_unitary should match explicit tensor-product application
#         for a simple 2-qubit example.
#         """
#         dims = np.array([2, 2], dtype=int)
#         D = int(np.prod(dims))

#         # Initial basis state |10> in ordering |00>,|01>,|10>,|11> -> index 2
#         psi = np.zeros(D, dtype=np.complex128)
#         psi[2] = 1.0

#         # Local X on qudit 0 (sparse)
#         U_loc_sparse = pauli_unitary_qudit(2, 1, 0)

#         # Helper application (X on qudit 0)
#         psi_helper = _apply_single_qudit_unitary(psi, dims, q=0, U_local=U_loc_sparse)

#         # Explicit tensor X ⊗ I in dense form
#         X_dense = U_loc_sparse.toarray()
#         I_dense = I_mat(2).toarray()
#         H_full = np.kron(X_dense, I_dense)  # 4x4 dense

#         psi_ref = H_full @ psi

#         assert psi_helper.shape == psi_ref.shape
#         assert np.allclose(psi_helper, psi_ref, atol=1e-12)

#     def test_apply_pauli_row_to_state_matches_full_operator(self):
#         """
#         apply_pauli_row_to_state should match the full operator built from a PauliSum
#         with a single term, for random state.
#         """
#         dims = np.array([2, 2], dtype=int)
#         D = int(np.prod(dims))

#         # Pauli: X on qudit 0, Z on qudit 1
#         x_exp = np.array([1, 0], dtype=int)
#         z_exp = np.array([0, 1], dtype=int)
#         tableau = np.concatenate([x_exp, z_exp])[None, :]  # (1, 4)
#         weights = np.array([1.0], dtype=np.complex128)

#         P_sum = PauliSum.from_tableau(tableau, weights=weights, dimensions=dims.tolist())

#         # Random state
#         rng = np.random.default_rng(123)
#         psi = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi /= np.linalg.norm(psi)

#         # Helper
#         psi_helper = apply_pauli_row_to_state(x_exp, z_exp, dims, psi)

#         # Full operator
#         P_full = P_sum.to_hilbert_space().toarray()
#         psi_ref = P_full @ psi

#         assert psi_helper.shape == psi_ref.shape
#         assert np.allclose(psi_helper, psi_ref, atol=1e-12)

#     def test_apply_paulisum_to_state_matches_hilbert(self):
#         """
#         apply_paulisum_to_state should match H.to_hilbert_space() @ psi for
#         small random PauliSum and state.
#         """
#         dims = [2, 2]
#         N = 2
#         D = 4

#         # Build a PauliSum with two terms: X0, Z1
#         tableau = np.zeros((2, 2 * N), dtype=int)
#         # Term 0: X0
#         tableau[0, 0] = 1
#         # Term 1: Z1
#         tableau[1, N + 1] = 1

#         weights = np.array([0.7 + 0.1j, -1.3], dtype=np.complex128)

#         H = PauliSum.from_tableau(tableau, weights=weights, dimensions=dims)

#         rng = np.random.default_rng(321)
#         psi = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi /= np.linalg.norm(psi)

#         # Helper
#         psi_helper = apply_paulisum_to_state(H, psi)

#         # Full Hilbert-space operator
#         H_full = H.to_hilbert_space().toarray()
#         psi_ref = H_full @ psi

#         assert psi_helper.shape == psi_ref.shape
#         assert np.allclose(psi_helper, psi_ref, atol=1e-12)

#     def test_apply_paulisum_identity(self):
#         """
#         If H is a single all-identity Pauli term with weight 1, H|psi> = |psi>.
#         """
#         dims = [2, 2]
#         N = 2
#         D = 4

#         # Tableau row of all zeros => identity on both qubits
#         tableau = np.zeros((1, 2 * N), dtype=int)
#         weights = np.array([1.0], dtype=np.complex128)

#         H = PauliSum.from_tableau(tableau, weights=weights, dimensions=dims)

#         rng = np.random.default_rng(777)
#         psi = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi /= np.linalg.norm(psi)

#         psi_out = apply_paulisum_to_state(H, psi)

#         assert psi_out.shape == psi.shape
#         assert np.allclose(psi_out, psi, atol=1e-12)

#     def test_project_to_symmetry_sector_yields_eigenstate(self):
#         """
#         project_to_symmetry_sector should return a state that is an eigenvector of S
#         with the requested eigenvalue, for a simple S = Z0 Z1 symmetry.
#         """
#         dims = np.array([2, 2], dtype=int)
#         D = 4

#         # Define S = Z0 Z1 via a PauliGate
#         x_exp = np.array([0, 0], dtype=int)
#         z_exp = np.array([1, 1], dtype=int)
#         ps = PauliString.from_exponents(x_exp, z_exp, dims)
#         S_gate_pauli = PauliGate(ps, name="Z0Z1")

#         # Circuit that just applies S once
#         C_S = Circuit(dimensions=dims.tolist(), gates=[S_gate_pauli])

#         # Random state
#         rng = np.random.default_rng(999)
#         psi = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi /= np.linalg.norm(psi)

#         # Project onto eigenvalue +1 sector
#         psi_proj = project_to_symmetry_sector(psi, dims, C_S, eigval=+1.0)

#         # Check norm
#         assert np.allclose(np.linalg.norm(psi_proj), 1.0, atol=1e-12)

#         # Apply S and check S|psi_proj> = +1 * |psi_proj>
#         state = State(psi_proj, dims)
#         state_S = C_S.act_on_state(state)
#         psi_S = state_S.as_array()

#         assert psi_S.shape == psi_proj.shape
#         assert np.allclose(psi_S, psi_proj, atol=1e-10)

#     def test_project_to_symmetry_sector_orthogonal_components_removed(self):
#         """
#         For S = Z0 Z1, basis states |01> and |10> have eigenvalue -1.
#         Projecting a pure |-1> eigenstate into the +1 sector should give zero norm
#         (caught as a ValueError in our implementation).
#         """
#         dims = np.array([2, 2], dtype=int)

#         x_exp = np.array([0, 0], dtype=int)
#         z_exp = np.array([1, 1], dtype=int)
#         ps = PauliString.from_exponents(x_exp, z_exp, dims)
#         S_gate_pauli = PauliGate(ps, name="Z0Z1")
#         C_S = Circuit(dimensions=dims.tolist(), gates=[S_gate_pauli])

#         # |01> basis state (pure eigenvector with eigenvalue -1)
#         psi = np.zeros(4, dtype=np.complex128)
#         psi[1] = 1.0

#         with pytest.raises(ValueError):
#             _ = project_to_symmetry_sector(psi, dims, C_S, eigval=+1.0)

#     def test_krylov_observable_expectation_is_real_for_hermitian(self):
#         """
#         For Hermitian H and O, <O(t)> should be (numerically) real;
#         imag part should be small.
#         """
#         N = 4
#         J = 0.5
#         h = 0.3

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         tableau = np.zeros((1, 2 * N), dtype=int)
#         tableau[0, 0] = 1
#         O_x = PauliSum.from_tableau(tableau, weights=[1.0], dimensions=[2] * N)

#         psi0, _ = product_state_ising(N, kind="x_plus")

#         times = np.linspace(0.0, 1.0, 21)

#         exp_krylov = krylov_observable_dynamics_symmetry(
#             H=H,
#             O=O_x,
#             T_gate=T_gate,
#             S_gate=S_gate,
#             psi0=psi0,
#             times=times,
#             m_max=32,
#             symmetry_eigval=None,
#         )

#         # Imag part should be tiny for Hermitian H and O
#         max_imag = np.max(np.abs(exp_krylov.imag))
#         assert max_imag < 1e-8

#     def test_T_conjugation_matches_hilbert(self):
#         """
#         Check that H' = T^{-1} H T computed at PauliSum level matches
#         U_T^dagger H_hilbert U_T at Hilbert level for a small system.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)

#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         # PauliSum-level conjugation
#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)

#         # Hilbert-space U_T
#         from sympleq.core.circuits import gate_to_circuit
#         C_T = gate_to_circuit(T_gate)
#         U_T = C_T.unitary().toarray()

#         # Compare matrices H_prime_hilbert vs U_T^\dagger H_hilbert U_T
#         H_h = H.to_hilbert_space().toarray()
#         H_prime_h = H_prime.to_hilbert_space().toarray()

#         H_prime_h_ref = U_T.conj().T @ H_h @ U_T

#         assert H_prime_h.shape == H_prime_h_ref.shape
#         assert np.allclose(H_prime_h, H_prime_h_ref, atol=1e-10)

#     def test_krylov_vs_full_hilbert_small_ising(self):
#         """
#         For a small Ising chain (N=4), verify that plain Krylov dynamics of <X_0(t)>
#         matches direct full Hilbert-space evolution to good accuracy.

#         This isolates the Lanczos/Krylov machinery from any symmetry transforms.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         # Hamiltonian
#         H = ising_chain_hamiltonian(N, J, h, periodic=True)

#         # Observable X_0
#         tableau = np.zeros((1, 2 * N), dtype=int)
#         tableau[0, 0] = 1
#         O_x = PauliSum.from_tableau(tableau, weights=[1.0], dimensions=[2] * N)

#         # Initial |+...+>
#         psi0, dims = product_state_ising(N, kind="x_plus")

#         # Time grid
#         times = np.linspace(0.0, 2.0, 51)

#         # Plain Krylov (no T, no S)
#         exp_krylov = krylov_observable_dynamics_plain(
#             H=H,
#             O=O_x,
#             psi0=psi0,
#             times=times,
#             m_max=16,   # dim(Hilbert) = 16, so 16 should be enough for exact recovery
#         )

#         # Full Hilbert reference
#         exp_full = observable_dynamics_full_hilbert(
#             H=H,
#             O=O_x,
#             psi0=psi0,
#             times=times,
#         )

#         assert exp_krylov.shape == exp_full.shape
#         max_diff = np.max(np.abs(exp_krylov - exp_full))
#         # With m_max=16 and N=4 this should basically be machine precision
#         assert max_diff < 1e-8

#     def test_krylov_symmetry_vs_full_hilbert_small_ising(self):
#         """
#         Sanity check: symmetry-based Krylov dynamics should still agree with full
#         Hilbert-space evolution, but this test also depends on T/S consistency.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)

#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         tableau = np.zeros((1, 2 * N), dtype=int)
#         tableau[0, 0] = 1
#         O_x = PauliSum.from_tableau(tableau, weights=[1.0], dimensions=[2] * N)

#         psi0, _ = product_state_ising(N, kind="x_plus")
#         times = np.linspace(0.0, 2.0, 51)

#         exp_krylov_sym = krylov_observable_dynamics_symmetry(
#             H=H,
#             O=O_x,
#             T_gate=T_gate,
#             S_gate=S_gate,
#             psi0=psi0,
#             times=times,
#             m_max=16,
#             symmetry_eigval=None,
#         )

#         exp_full = observable_dynamics_full_hilbert(
#             H=H,
#             O=O_x,
#             psi0=psi0,
#             times=times,
#         )

#         assert exp_krylov_sym.shape == exp_full.shape
#         max_diff = np.max(np.abs(exp_krylov_sym - exp_full))
#         # Depending on how confident you are in T/S, you can tune this:
#         assert max_diff < 1e-6

#     def test_krylov_symmetry_sector_vs_full_sector(self):
#         """
#         Explicitly exploit a simple Z2 symmetry S = Z0 Z1 on a 2-qubit system.

#         We:
#           - build H that commutes with S (diagonal in Z),
#           - build O = X0,
#           - choose a random initial state |psi0>,
#           - restrict to the +1 symmetry sector of S in full Hilbert space,
#           - compute <O(t)> in that sector exactly,
#           - compare with krylov_observable_dynamics_symmetry using symmetry_eigval=+1.

#         This tests the *symmetry-exploiting* branch of
#         krylov_observable_dynamics_symmetry (i.e. the path that uses
#         T_gate, S_gate, project_to_symmetry_sector, and H' / O').
#         """
#         N = 2
#         dims = [2, 2]
#         D = 4

#         # --- Symmetry S = Z0 Z1 as a PauliGate ---
#         x_exp_S = np.array([0, 0], dtype=int)
#         z_exp_S = np.array([1, 1], dtype=int)
#         ps_S = PauliString.from_exponents(x_exp_S, z_exp_S, dims)
#         S_gate = PauliGate(ps_S, name="Z0Z1")

#         # --- Identity T_gate (so H' = H, O' = O, psi' = psi) ---
#         symp_id = np.eye(2 * N, dtype=int)
#         phase_zero = np.zeros(2 * N, dtype=int)
#         T_gate = Gate(
#             "Id",
#             list(range(N)),
#             symp_id,
#             dimensions=dims,
#             phase_vector=phase_zero,
#         )

#         # --- Hamiltonian H = Z0 + 0.5 Z1 (commutes with S) ---
#         tableau_H = np.zeros((2, 2 * N), dtype=int)
#         # Term 0: Z0 -> [x0,x1 | z0,z1] = [0,0 | 1,0]
#         tableau_H[0, N + 0] = 1
#         # Term 1: Z1 -> [0,0 | 0,1]
#         tableau_H[1, N + 1] = 1

#         weights_H = np.array([1.0, 0.5], dtype=np.complex128)
#         H = PauliSum.from_tableau(tableau_H, weights=weights_H, dimensions=dims)

#         # --- Observable O = X0 ---
#         tableau_O = np.zeros((1, 2 * N), dtype=int)
#         tableau_O[0, 0] = 1  # X on qudit 0
#         Op = PauliSum.from_tableau(tableau_O, weights=[1.0], dimensions=dims)

#         # --- Random initial state |psi0> with components in both sectors ---
#         rng = np.random.default_rng(2025)
#         psi0 = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi0 /= np.linalg.norm(psi0)

#         times = np.linspace(0.0, 3.0, 61)

#         # --- Full Hilbert-space reference in +1 sector of S ---

#         # Full matrix
#         S_h = S_gate.unitary(dims=dims).toarray()

#         # Projector onto +1 sector: P_+ = (I + S)/2
#         I_full = np.eye(D, dtype=np.complex128)
#         P_plus = 0.5 * (I_full + S_h)

#         psi_plus = P_plus @ psi0
#         norm_plus = np.linalg.norm(psi_plus)
#         # Make sure we didn't project to (near) zero
#         assert norm_plus > 1e-10
#         psi_plus /= norm_plus

#         # Exact sector dynamics: evolve psi_plus under H, compute <O(t)>
#         exp_sector_full = observable_dynamics_full_hilbert(
#             H=H,
#             O=Op,
#             psi0=psi_plus,
#             times=times,
#         )

#         # --- Krylov symmetry-based dynamics ---
#         exp_krylov_sym = krylov_observable_dynamics_symmetry(
#             H=H,
#             O=Op,
#             T_gate=T_gate,
#             S_gate=S_gate,
#             psi0=psi0,
#             times=times,
#             m_max=4,         # Hilbert dim is 4; 4 is enough for exact recovery
#             symmetry_eigval=+1.0,
#         )

#         assert exp_krylov_sym.shape == exp_sector_full.shape
#         max_diff = np.max(np.abs(exp_krylov_sym - exp_sector_full))
#         # With such a tiny system and m_max=4, we expect near machine precision
#         assert max_diff < 1e-8

#     def test_symmetry_projection_matches_hilbert_small_ising(self):
#         """
#         For small N, check that project_to_symmetry_sector using S_gate
#         matches the explicit +1 projector built from U_S in the computational basis.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)

#         dims = np.asarray([2] * N, dtype=int)
#         D = int(np.prod(dims))

#         rng = np.random.default_rng(4242)
#         psi0 = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi0 /= np.linalg.norm(psi0)

#         # S in Hilbert space (computational basis)
#         from sympleq.core.circuits import gate_to_circuit
#         C_S = gate_to_circuit(S_gate)
#         U_S = C_S.unitary().toarray()

#         # Explicit +1 projector: P_+ = (I + U_S)/2
#         I_full = np.eye(D, dtype=np.complex128)
#         P_plus = 0.5 * (I_full + U_S)

#         psi_proj_ref = P_plus @ psi0
#         norm_ref = np.linalg.norm(psi_proj_ref)
#         assert norm_ref > 1e-10
#         psi_proj_ref /= norm_ref

#         # Code path: project_to_symmetry_sector(psi0, S_gate)
#         psi_proj_code = project_to_symmetry_sector(
#             psi0, dims, C_S, eigval=+1.0
#         )

#         # Compare up to a global phase
#         overlap = np.vdot(psi_proj_ref, psi_proj_code)
#         assert np.abs(overlap) > 1 - 1e-8

#     def test_Hprime_commutes_with_S(self):
#         """
#         Check [H', S] ≈ 0 in Hilbert space for small N, where
#         H' = T^{-1} H T and S is the symmetry from the heuristic Clifford.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)

#         H_prime_h = H_prime.to_hilbert_space().toarray()

#         from sympleq.core.circuits import gate_to_circuit
#         C_S = gate_to_circuit(S_gate)
#         U_S = C_S.unitary().toarray()

#         comm = H_prime_h @ U_S - U_S @ H_prime_h
#         norm_comm = np.linalg.norm(comm)
#         assert norm_comm < 1e-10

#     def test_symmetry_sector_dynamics_matches_full_sector_small_ising(self):
#         """
#         Small-N check: symmetry-aware Krylov dynamics in the +1 sector
#         (in the H' = T^-1 H T frame) matches explicit sector dynamics
#         computed in Hilbert space.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         # Original Hamiltonian
#         H = ising_chain_hamiltonian(N, J, h, periodic=True)

#         # Symmetry via heuristic Clifford
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         # Observable X_0
#         tableau = np.zeros((1, 2 * N), dtype=int)
#         tableau[0, 0] = 1
#         Op = PauliSum.from_tableau(tableau, weights=[1.0], dimensions=[2] * N)

#         dims = np.asarray([2] * N, dtype=int)
#         D = int(np.prod(dims))

#         rng = np.random.default_rng(1234)
#         psi0 = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi0 /= np.linalg.norm(psi0)

#         times = np.linspace(0.0, 2.0, 41)

#         # ----- Reference: full Hilbert-space sector dynamics in H' frame -----
#         # H' and O' in Pauli and Hilbert form
#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)
#         O_prime = T_inv_gate.act(Op)

#         H_prime_h = H_prime.to_hilbert_space().toarray()
#         O_prime_h = O_prime.to_hilbert_space().toarray()

#         # Unitary for T and S
#         C_T = gate_to_circuit(T_gate)
#         U_T = C_T.unitary().toarray()

#         C_S = gate_to_circuit(S_gate)
#         U_S = C_S.unitary().toarray()

#         # Initial state in H' frame: |psi0'> = U_T^\dagger |psi0>
#         psi0_prime = U_T.conj().T @ psi0

#         # Project |psi0'> onto +1 eigenspace of S (in H' basis)
#         I_full = np.eye(D, dtype=np.complex128)
#         P_plus = 0.5 * (I_full + U_S)
#         psi_plus_prime = P_plus @ psi0_prime
#         norm_plus = np.linalg.norm(psi_plus_prime)
#         assert norm_plus > 1e-10
#         psi_plus_prime /= norm_plus

#         # Exact spectral evolution under H'
#         evals, V = np.linalg.eigh(H_prime_h)
#         V_dag = V.conj().T
#         coeff0 = V_dag @ psi_plus_prime

#         exp_exact = np.empty_like(times, dtype=np.complex128)
#         for idx, t in enumerate(times):
#             phase = np.exp(-1j * evals * t)
#             psi_t = V @ (phase * coeff0)
#             exp_exact[idx] = np.vdot(psi_t, O_prime_h @ psi_t)

#         # ----- Code path: symmetry-aware Krylov in +1 sector (no Hilbert matrices) -----
#         exp_krylov_sym = krylov_observable_dynamics_symmetry(
#             H=H,
#             O=Op,
#             T_gate=T_gate,
#             S_gate=S_gate,
#             psi0=psi0,
#             times=times,
#             m_max=16,          # D=16, so this is enough to be essentially exact
#             symmetry_eigval=+1.0,
#         )

#         assert exp_krylov_sym.shape == exp_exact.shape
#         max_diff = np.max(np.abs(exp_krylov_sym - exp_exact))
#         # Tune tolerance if needed, but with D=16 and m_max=16 this should be tiny
#         assert max_diff < 1e-7

#     def test_plain_krylov_on_Hprime_matches_full_Hprime_small_ising(self):
#         """
#         Check that plain Krylov dynamics using H' (T^-1 H T) matches exact
#         Hilbert-space dynamics under H' for a small system and random state.
#         This does NOT use any symmetry projection.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         # Conjugate H and O by T^-1
#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)

#         # Take a simple observable, say X_0 transformed as well
#         tableau = np.zeros((1, 2 * N), dtype=int)
#         tableau[0, 0] = 1
#         Op = PauliSum.from_tableau(tableau, weights=[1.0], dimensions=[2] * N)
#         O_prime = T_inv_gate.act(Op)

#         dims = np.asarray([2] * N, dtype=int)
#         D = int(np.prod(dims))
#         rng = np.random.default_rng(999)
#         psi0 = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi0 /= np.linalg.norm(psi0)

#         times = np.linspace(0.0, 2.0, 41)

#         # Exact Hilbert-space dynamics under H'
#         H_prime_h = H_prime.to_hilbert_space().toarray()
#         O_prime_h = O_prime.to_hilbert_space().toarray()

#         evals, V = np.linalg.eigh(H_prime_h)
#         V_dag = V.conj().T
#         coeff0 = V_dag @ psi0

#         exp_exact = np.empty_like(times, dtype=np.complex128)
#         for idx, t in enumerate(times):
#             phase = np.exp(-1j * evals * t)
#             psi_t = V @ (phase * coeff0)
#             exp_exact[idx] = np.vdot(psi_t, O_prime_h @ psi_t)

#         # Plain Krylov on H' (no symmetry)
#         exp_krylov = krylov_observable_dynamics_plain(
#             H=H_prime,
#             O=O_prime,
#             psi0=psi0,
#             times=times,
#             m_max=16,  # H' acts on 16-dim Hilbert space
#         )

#         assert exp_krylov.shape == exp_exact.shape
#         max_diff = np.max(np.abs(exp_krylov - exp_exact))
#         assert max_diff < 1e-7

#     def test_paulisum_matvec_matches_hilbert_for_Hprime(self):
#         """
#         For H' = T^-1 H T, check that apply_paulisum_to_state(H', v)
#         matches H_prime.to_hilbert_space() @ v for random v.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)
#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)

#         H_prime_h = H_prime.to_hilbert_space().toarray()

#         dims = np.asarray(H_prime.dimensions, dtype=int)
#         D = int(np.prod(dims))
#         rng = np.random.default_rng(444)
#         v = rng.normal(size=D) + 1j * rng.normal(size=D)
#         v /= np.linalg.norm(v)

#         # Hilbert-space matvec
#         y_ref = H_prime_h @ v

#         # PauliSum-based matvec
#         y_code = apply_paulisum_to_state_dense(H_prime, v)

#         assert y_code.shape == y_ref.shape
#         max_diff = np.max(np.abs(y_code - y_ref))
#         assert max_diff < 1e-10

#     def test_single_term_matches_hilbert_for_Hprime(self):
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)
#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)

#         dims = np.asarray(H_prime.dimensions, dtype=int)

#         tableau = H_prime.tableau        # (n_terms, 2N)
#         weights = np.asarray(H_prime.weights, dtype=np.complex128).reshape(-1)
#         phases = np.asarray(H_prime.phases, dtype=int).reshape(-1)
#         n_terms, twoN = tableau.shape
#         n_qudits = twoN // 2

#         for k in range(n_terms):
#             w_k = weights[k]
#             if np.abs(w_k) == 0:
#                 continue

#             row = tableau[k]
#             x_exp = row[:n_qudits]
#             z_exp = row[n_qudits:]

#             # PauliString for this term
#             ps = PauliString.from_exponents(x_exp, z_exp, dims)
#             gate_k = PauliGate(ps)
#             Pk_unitary = gate_k.unitary(dims=dims).toarray()  # (D, D)

#             # Build a PauliSum with *only* this one Pauli term,
#             # and use *exactly* the phase/weight data from H_prime:
#             tableau_k = row.reshape(1, -1)
#             weights_k = np.array([w_k], dtype=np.complex128)
#             phases_k = np.array([phases[k]], dtype=int)

#             Hk = PauliSum(
#                 tableau=tableau_k,
#                 weights=weights_k,
#                 phases=phases_k,
#                 dimensions=dims,
#             )

#             Hk_h = Hk.to_hilbert_space().toarray()  # (D, D)

#             # Now we want to know: what complex scalar c_k satisfies
#             #    Hk_h  ≈  c_k * Pk_unitary ?
#             #
#             # Use Hilbert-Schmidt inner product to estimate c_k:
#             num = np.trace(Pk_unitary.conj().T @ Hk_h)
#             den = np.trace(Pk_unitary.conj().T @ Pk_unitary)  # = D for unitary Pauli
#             c_k = num / den

#             # Compare Hk_h with c_k * Pk_unitary
#             diff = np.max(np.abs(Hk_h - c_k * Pk_unitary))
#             assert diff < 1e-12, f"Term {k} does not match single-Pauli decomposition"

#     def test_paulisum_matvec_matches_hilbert_for_original_H(self):
#         N = 4
#         J = 1.0
#         h = 0.7
#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         H_h = H.to_hilbert_space().toarray()

#         dims = np.asarray(H.dimensions, dtype=int)
#         D = int(np.prod(dims))
#         rng = np.random.default_rng(555)
#         v = rng.normal(size=D) + 1j * rng.normal(size=D)
#         v /= np.linalg.norm(v)

#         y_ref = H_h @ v
#         y_code = apply_paulisum_to_state_dense(H, v)
#         max_diff = np.max(np.abs(y_code - y_ref))
#         assert max_diff < 1e-10

#     def test_symmetry_projection_statevector_matches_full_hilbert(self):
#         """
#         Projection via project_to_symmetry_sector_statevector in the basis where
#         S_gate acts should match the explicit Hilbert projector
#             P_+ = (I + U_S)/2
#         constructed from the unitary U_S of S_gate.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)

#         dims = np.asarray([2] * N, dtype=int)
#         D = int(np.prod(dims))

#         rng = np.random.default_rng(4242)
#         psi = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi /= np.linalg.norm(psi)

#         # Hilbert-space U_S and projector P_+
#         C_S = gate_to_circuit(S_gate)
#         U_S = C_S.unitary().toarray()
#         I_full = np.eye(D, dtype=np.complex128)
#         P_plus = 0.5 * (I_full + U_S)

#         psi_ref = P_plus @ psi
#         norm_ref = np.linalg.norm(psi_ref)
#         assert norm_ref > 1e-10
#         psi_ref /= norm_ref

#         # Code path: projector via state-level action of S_gate
#         psi_code = project_to_symmetry_sector_statevector(
#             psi, dims, S_gate, eigval=+1.0
#         )

#         assert psi_code.shape == psi_ref.shape
#         # Compare up to a global phase
#         overlap = np.vdot(psi_ref, psi_code)
#         assert np.abs(overlap) > 1.0 - 1e-8


# class TestSymmetryReducedKrylov:

#     def test_sector_basis_is_eigenspace_of_S(self):
#         """
#         Check that build_symmetry_sector_basis returns vectors that satisfy
#         S |phi> = target_eigval |phi>, using the Hilbert-space unitary of S.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)

#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)

#         dims = [2] * N
#         blocks = block_indexes(S_gate.symplectic)

#         B_plus, lam_vec = build_symmetry_sector_basis(
#             S_gate=S_gate,
#             blocks=blocks,
#             dims=dims,
#             target_eigval=+1.0,
#         )

#         # Hilbert-space unitary for S
#         C_S = gate_to_circuit(S_gate)
#         U_S = C_S.unitary().toarray()

#         # Check that each column v satisfies U_S v ≈ v (since eigenvalue +1)
#         for j in range(B_plus.shape[1]):
#             v = B_plus[:, j]
#             lhs = U_S @ v
#             assert np.allclose(lhs, v, atol=1e-8)

#         # Check orthonormality of B_plus
#         G = B_plus.conj().T @ B_plus
#         assert np.allclose(G, np.eye(G.shape[0]), atol=1e-8)

#     def test_reduced_H_matches_sector_restriction(self):
#         """
#         For small N, check that the reduced Hamiltonian H_sec constructed from
#         B_sector is equivalent (up to a unitary change of basis) to the
#         restriction of H' to the +1 eigenspace of S, computed directly from
#         the eigen-decomposition of U_S.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         dims = [2] * N

#         # H' = T^{-1} H T
#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)
#         H_prime_h = H_prime.to_hilbert_space().toarray()

#         # --- Reference: restrict H' to +1 eigenspace of S via global eigendecomp ---

#         # Global U_S
#         C_S = gate_to_circuit(S_gate)
#         U_S = C_S.unitary().toarray()

#         # Diagonalise U_S
#         evals_S, vecs_S = np.linalg.eig(U_S)

#         # Pick eigenvectors with eigenvalue close to +1
#         tol = 1e-8
#         idx_plus = [i for i, lam in enumerate(evals_S) if np.abs(lam - 1.0) < tol]
#         assert len(idx_plus) > 0, "No +1 eigenvectors found for S."

#         B_ref = vecs_S[:, idx_plus]  # (D, D_plus)

#         # Orthonormalise
#         Q_ref, _ = np.linalg.qr(B_ref)
#         H_sec_ref = Q_ref.conj().T @ H_prime_h @ Q_ref

#         # --- Our reduced H_sec via block-wise sector basis ---

#         from sympleq.core.symmetries.dynamics.symmetry_reduced_krylov import (
#             build_symmetry_sector_basis,
#             reduce_operator_to_sector,
#         )
#         from sympleq.core.symmetries.block_decomposition import block_indexes

#         blocks = block_indexes(S_gate.symplectic)
#         B_plus, _ = build_symmetry_sector_basis(
#             S_gate=S_gate, blocks=blocks, dims=dims, target_eigval=+1.0
#         )
#         H_sec = reduce_operator_to_sector(H_prime_h, B_plus)

#         # Compare spectra (spaces are the same up to a unitary within the sector)
#         evals_ref = np.linalg.eigvalsh(H_sec_ref)
#         evals = np.linalg.eigvalsh(H_sec)

#         assert len(evals_ref) == len(evals)
#         assert np.allclose(np.sort(evals_ref), np.sort(evals), atol=1e-8)

#     def test_symmetry_reduced_krylov_matches_full_sector_dynamics(self):
#         """
#         For small N, check that krylov_observable_dynamics_symmetry_reduced
#         matches explicit spectral evolution in the +1 eigenspace of S for H'.

#         The +1 sector is defined using the global eigen-decomposition of U_S,
#         NOT via (I + U_S)/2, which is only a projector when S^2 = I and
#         spec(S) ⊂ {±1}.
#         """
#         N = 4
#         J = 1.0
#         h = 0.7

#         H = ising_chain_hamiltonian(N, J, h, periodic=True)
#         F = heuristic_clifford_symmetry(N)
#         S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)
#         h_S, h_T = clifford_phase_decomposition(
#             F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
#         )
#         S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
#         T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

#         # Observable X_0
#         tableau = np.zeros((1, 2 * N), dtype=int)
#         tableau[0, 0] = 1
#         Op = PauliSum.from_tableau(tableau, weights=[1.0], dimensions=[2] * N)

#         D = 2**N

#         rng = np.random.default_rng(1234)
#         psi0 = rng.normal(size=D) + 1j * rng.normal(size=D)
#         psi0 /= np.linalg.norm(psi0)

#         times = np.linspace(0.0, 2.0, 41)

#         # ----- Reference: full Hilbert-space sector dynamics in H' frame -----

#         # H', O' in Pauli and Hilbert form
#         T_inv_gate = T_gate.inv()
#         H_prime = T_inv_gate.act(H)
#         O_prime = T_inv_gate.act(Op)
#         H_prime_h = H_prime.to_hilbert_space().toarray()
#         O_prime_h = O_prime.to_hilbert_space().toarray()

#         # Unitary for T and S
#         C_T = gate_to_circuit(T_gate)
#         U_T = C_T.unitary().toarray()

#         C_S = gate_to_circuit(S_gate)
#         U_S = C_S.unitary().toarray()

#         # Initial state in H' frame: |psi0'> = U_T^\dagger |psi0>
#         psi0_prime = U_T.conj().T @ psi0

#         # Diagonalise U_S and build projector onto +1 eigenspace
#         evals_S, vecs_S = np.linalg.eig(U_S)
#         tol = 1e-8
#         idx_plus = [i for i, lam in enumerate(evals_S) if np.abs(lam - 1.0) < tol]
#         assert len(idx_plus) > 0, "No +1 eigenvectors found for S."

#         B_plus = vecs_S[:, idx_plus]
#         Q_plus, _ = np.linalg.qr(B_plus)  # orthonormal basis of +1 sector

#         # Project |psi0'> into +1 sector
#         psi_plus_prime = Q_plus @ (Q_plus.conj().T @ psi0_prime)
#         norm_plus = np.linalg.norm(psi_plus_prime)
#         assert norm_plus > 1e-10
#         psi_plus_prime /= norm_plus

#         # Exact spectral evolution under H' (state remains in +1 sector
#         # because [H', S] = 0)
#         evals_H, V_H = np.linalg.eigh(H_prime_h)
#         V_H_dag = V_H.conj().T
#         coeff0 = V_H_dag @ psi_plus_prime

#         exp_exact = np.empty_like(times, dtype=np.complex128)
#         for idx_t, t in enumerate(times):
#             phase = np.exp(-1j * evals_H * t)
#             psi_t = V_H @ (phase * coeff0)
#             exp_exact[idx_t] = np.vdot(psi_t, O_prime_h @ psi_t)

#         # ----- Code path: symmetry-reduced Krylov in +1 sector -----

#         from sympleq.core.symmetries.dynamics.symmetry_reduced_krylov import (
#             krylov_observable_dynamics_symmetry_reduced,
#         )

#         exp_krylov = krylov_observable_dynamics_symmetry_reduced(
#             H=H,
#             O=Op,
#             T_gate=T_gate,
#             S_gate=S_gate,
#             psi0=psi0,
#             times=times,
#             symmetry_eigval=+1.0,
#             m_max=16,   # D_sec is small; 16 is plenty for exact recovery
#         )

#         assert exp_krylov.shape == exp_exact.shape
#         max_diff = np.max(np.abs(exp_krylov - exp_exact))
#         # With the reduced-space Lanczos, this should be very tight
#         assert max_diff < 1e-7

