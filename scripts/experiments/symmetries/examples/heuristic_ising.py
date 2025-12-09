import numpy as np
from sympleq.core.symmetries.clifford import clifford_phase_decomposition
from sympleq.core.symmetries.block_decomposition import block_decompose_optimal, block_indexes
from sympleq.models.Ising import ising_chain_hamiltonian, heuristic_clifford_symmetry
from sympleq.core.circuits import Gate, Circuit, gate_to_circuit
from sympleq.core.paulis import PauliSum
from typing import List, Tuple


def _local_U_eig(U_block: np.ndarray, tol: float = 1e-10) -> List[Tuple[complex, np.ndarray]]:
    """Return eigenpairs of U_block (unitary), keeping degeneracy info."""
    eigvals_u, eigvecs_u = np.linalg.eig(U_block)

    # Group indices by eigenvalue up to tolerance
    groups: List[List[int]] = []
    used = np.zeros(len(eigvals_u), dtype=bool)
    for i, lam in enumerate(eigvals_u):
        if used[i]:
            continue
        group = [i]
        used[i] = True
        for j in range(i + 1, len(eigvals_u)):
            if used[j]:
                continue
            if np.abs(lam - eigvals_u[j]) < tol:
                group.append(j)
                used[j] = True
        groups.append(group)

    out = []
    for g in groups:
        for idx in g:
            vec = eigvecs_u[:, idx]
            vec = vec / np.linalg.norm(vec)
            out.append((eigvals_u[idx], vec))
    return out


def blockwise_joint_eigensystem(H_prime: PauliSum,
                                C_S: Circuit,
                                blocks: List[List[int]],
                                tol: float = 1e-6) -> Tuple[List[complex], List[float], List[np.ndarray]]:
    """
    Build a joint eigenbasis of U_S and H' by:
      1) Tensoring eigenvectors of each U_block to form a basis of each global U eigen-subspace.
      2) Diagonalizing H' restricted to each U eigen-subspace.
    Returns (eigvals_U, eigvals_H, eigenvectors).
    """
    # Step 1: local eigenpairs of U blocks
    local_eig = []
    for block in blocks:
        U_block = C_S.local_circuit(block).unitary().toarray()
        local_eig.append(_local_U_eig(U_block, tol=tol))

    # Build global basis from tensor products of local eigenvectors
    global_u = []
    global_vecs = []
    for choice in np.array(np.meshgrid(*[range(len(le)) for le in local_eig])).T.reshape(-1, len(local_eig)):
        lam_u = 1.0 + 0j
        vec = np.array([1.0], dtype=complex)
        for block_idx, sel in enumerate(choice):
            u_val, v = local_eig[block_idx][sel]
            lam_u *= u_val
            vec = np.kron(vec, v)
        vec = vec / np.linalg.norm(vec)
        global_u.append(lam_u)
        global_vecs.append(vec)

    # Step 2: diagonalize H' within each degenerate U eigen-subspace
    H_full = H_prime.to_hilbert_space().toarray()
    eigvals_U_out: List[complex] = []
    eigvals_H_out: List[float] = []
    eigvecs_out: List[np.ndarray] = []

    used = np.zeros(len(global_u), dtype=bool)
    for i, lam in enumerate(global_u):
        if used[i]:
            continue
        indices = [i]
        used[i] = True
        for j in range(i + 1, len(global_u)):
            if used[j]:
                continue
            if np.abs(lam - global_u[j]) < tol:
                indices.append(j)
                used[j] = True

        V = np.stack([global_vecs[k] for k in indices], axis=1)
        # Orthonormalize the basis for the degenerate subspace
        V, _ = np.linalg.qr(V)
        H_sub = V.conj().T @ H_full @ V
        eigvals_h, vecs_h = np.linalg.eigh(H_sub)
        for idx, lam_h in enumerate(eigvals_h):
            vec_global = V @ vecs_h[:, idx]
            vec_global = vec_global / np.linalg.norm(vec_global)
            eigvals_U_out.append(lam)
            eigvals_H_out.append(float(lam_h))
            eigvecs_out.append(vec_global)

    return eigvals_U_out, eigvals_H_out, eigvecs_out


def local_joint_eigensystem_for_block(H_prime: PauliSum,
                                      C_S: Circuit,
                                      block: List[int],
                                      tol: float = 1e-8) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute simultaneous eigenpairs (U_block, H'_block) for a single block.
    Returns (eigvals_U, eigvals_H, eigvecs) with eigvecs as columns.
    """
    U_block = C_S.local_circuit(block).unitary().toarray()
    pauli_indices = list(range(H_prime.n_paulis()))
    H_block_ps = H_prime.get_subspace(qudit_indices=block, pauli_indices=pauli_indices)
    H_block = H_block_ps.to_hilbert_space().toarray()

    eigvals_u, eigvecs_u = np.linalg.eig(U_block)

    # Group by U eigenvalue (tolerance to handle numerical noise)
    groups: List[List[int]] = []
    used = np.zeros(len(eigvals_u), dtype=bool)
    for i, lam in enumerate(eigvals_u):
        if used[i]:
            continue
        group = [i]
        used[i] = True
        for j in range(i + 1, len(eigvals_u)):
            if used[j]:
                continue
            if np.abs(lam - eigvals_u[j]) < tol:
                group.append(j)
                used[j] = True
        groups.append(group)

    eigvals_u_out: List[complex] = []
    eigvals_h_out: List[float] = []
    eigvecs_out: List[np.ndarray] = []

    for g in groups:
        Vg = eigvecs_u[:, g]
        # Orthonormalize within the degenerate subspace
        Vg, _ = np.linalg.qr(Vg)
        H_sub = Vg.conj().T @ H_block @ Vg
        evals_h, evecs_h = np.linalg.eigh(H_sub)
        for idx, lam_h in enumerate(evals_h):
            vec = Vg @ evecs_h[:, idx]
            vec = vec / np.linalg.norm(vec)
            eigvals_u_out.append(eigvals_u[g[0]])
            eigvals_h_out.append(float(lam_h))
            eigvecs_out.append(vec)

    return np.array(eigvals_u_out), np.array(eigvals_h_out), np.column_stack(eigvecs_out)


def evolve_state_in_block(state: np.ndarray, eigvals_h: np.ndarray, eigvecs: np.ndarray, t: float) -> np.ndarray:
    """
    Evolve a state |psi> under H_block with eigen-decomposition (eigvals_h, eigvecs).
    state: column vector in the block Hilbert space.
    """
    coeffs = eigvecs.conj().T @ state
    phases = np.exp(-1j * eigvals_h * t)
    return eigvecs @ (coeffs * phases)


def run_example():
    """This example shows we can create joint eigenstates
      but requires full Hilbert space representations in the check """
    N = 10
    J = 1
    h = 0.5
    H = ising_chain_hamiltonian(N, J, h, periodic=True)

    F = heuristic_clifford_symmetry(N)
    S, T = block_decompose_optimal(F.symplectic, 2)

    h_S, h_T = clifford_phase_decomposition(F.symplectic, F.phase_vector, S, T, int(H.lcm))
    S_gate = Gate('S', F.qudit_indices, S, F.dimensions, h_S)
    T_gate = Gate('T', F.qudit_indices, T, F.dimensions, h_T)

    # Conjugate H by T to get H' commuting with S
    H_prime = T_gate.inv().act(H)

    C_S = gate_to_circuit(S_gate)
    U_S = C_S.unitary().toarray()

    blocks = block_indexes(S_gate.symplectic)
    print("Blocks:", blocks)

    eigvals_U, eigvals_H, eig_vecs = blockwise_joint_eigensystem(H_prime, C_S, blocks)

    # Verification against full matrices (small N so still affordable)
    H_prime_hilbert = H_prime.to_hilbert_space().toarray()
    for lam_u, lam_h, v in zip(eigvals_U, eigvals_H, eig_vecs):
        lhs_u = U_S @ v
        lhs_h = H_prime_hilbert @ v
        assert np.allclose(lhs_u, lam_u * v, atol=1e-8)
        assert np.allclose(lhs_h, lam_h * v, atol=1e-8)
    print("Verified joint eigenvectors for U_S and H'.")

    # Demonstrate block-local evolution without building full H'
    block_systems = []
    for block in blocks:
        eig_u_b, eig_h_b, eigvecs_b = local_joint_eigensystem_for_block(H_prime, C_S, block)
        block_systems.append((eig_h_b, eigvecs_b))

    t = 0.123
    evolved_blocks = []
    for (eig_h_b, eigvecs_b), block in zip(block_systems, blocks):
        dim_block = int(np.prod([H_prime.dimensions[i] for i in block]))
        psi0 = np.zeros(dim_block, dtype=complex)
        psi0[0] = 1.0  # |00...> on the block
        psi_t = evolve_state_in_block(psi0, eig_h_b, eigvecs_b, t)
        evolved_blocks.append(psi_t)

    print("Evolved each block locally without forming the full Hamiltonian.")


if __name__ == "__main__":
    run_example()
