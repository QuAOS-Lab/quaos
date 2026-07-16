import numpy as np
from typing import Sequence, Tuple, List

from sympleq.models.Ising import ising_chain_hamiltonian, heuristic_clifford_symmetry, product_state_ising
from sympleq.core.symmetries.block_decomposition import block_decompose_optimal, block_indexes
from sympleq.core.symmetries.clifford import clifford_phase_decomposition
from sympleq.core.circuits import Gate, Circuit, gate_to_circuit
from sympleq.core.paulis import PauliSum
from sympleq.core.symmetries.dynamics.krylov_time_evolve import krylov_observable_dynamics_symmetry
import matplotlib.pyplot as plt
from time import time


def _local_U_eig(U_block: np.ndarray, tol: float = 1e-10) -> List[Tuple[complex, np.ndarray]]:
    """Return eigen-pairs of U_block (unitary), keeping degeneracy info."""
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
    # Step 1: local eigen-pairs of U blocks
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


def observable_dynamics_from_joint_eigenbasis(
    H: PauliSum,
    O: PauliSum,
    S_gate: Gate,
    T_gate: Gate,
    blocks: List[List[int]],
    psi0: np.ndarray,
    times: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """
    Inefficient proof of concept: Evaluate ⟨O(t)⟩ using the simultaneous eigenbasis of H' and S.

    Compute ⟨O(t)⟩ using the simultaneous eigenbasis of H' and S.

    Parameters
    ----------
    H : PauliSum
        Original Hamiltonian (in original basis).
    O : PauliSum
        Observable of interest (in original basis).
    S_gate : Gate
        Gate implementing symmetry S (in the T-basis picture).
    T_gate : Gate
        Gate implementing T, with H' = T^{-1} H T commuting with S.
    blocks : list of blocks (from block_indexes(S_gate.symplectic))
    psi0 : (D,) complex ndarray
        Initial statevector in the original basis.
    times : 1D array-like
        Times at which to evaluate ⟨O(t)⟩.

    Returns
    -------
    exp_t : (len(times),) complex ndarray
        Expectation values ⟨O(t)⟩.
    """
    times = np.asarray(times, dtype=float)
    dims = np.asarray(H.dimensions, dtype=int)
    D = int(np.prod(dims))
    assert psi0.shape == (D,)

    # 1) Conjugate H and O into the T-basis:
    #    H' = T^{-1} H T,  O' = T^{-1} O T
    T_inv_gate = T_gate.inv()
    H_prime = T_inv_gate.act(H)
    O_prime = T_inv_gate.act(O)

    # 2) Build circuit for S and get joint eigenbasis of H' and S
    C_S = gate_to_circuit(S_gate)

    eigvals_U, eigvals_H, eig_vecs = blockwise_joint_eigensystem(
        H_prime, C_S, blocks
    )

    # Sanity: pack eigenvectors into a unitary-like matrix V
    V = np.column_stack(eig_vecs)  # shape (D, D_expected)
    if V.shape[0] != D:
        raise ValueError(f"Eigenvector dimension {V.shape} inconsistent with Hilbert dimension {D}.")
    if V.shape[1] != D:
        print("Warning: number of joint eigenvectors != full Hilbert dimension; "
              "are you restricting to a subspace?")

    # 3) Build O' in the Hilbert space and then in the joint eigenbasis
    O_prime_hilbert = O_prime.to_hilbert_space().toarray()   # (D, D)
    O_eig = V.conj().T @ O_prime_hilbert @ V                 # O in eigenbasis of H'

    # 4) Transform initial state ψ0 into the T-basis:
    #    |ψ'⟩ = U_T |ψ⟩
    C_T = gate_to_circuit(T_gate)
    U_T = C_T.unitary().toarray()           # (D, D)
    psi_prime0 = U_T @ psi0                 # in H'/S basis

    # 5) Expand ψ' in the joint eigenbasis: c_k = v_k† ψ'
    c = V.conj().T @ psi_prime0             # (D,)

    # 6) Time evolution in the eigenbasis:
    #    |ψ'(t)⟩ = ∑_k c_k e^{-i E_k t} |v_k⟩
    #    ⟨O(t)⟩ = ⟨ψ'(t)| O' |ψ'(t)⟩ = c(t)† O_eig c(t)
    E = np.asarray(eigvals_H, dtype=float)
    if E.shape[0] != c.shape[0]:
        raise ValueError("Mismatch between number of eigenvalues and expansion coefficients.")

    exp_t = np.empty(times.shape[0], dtype=np.complex128)
    for idx, t in enumerate(times):
        phase = np.exp(-1j * E * t)   # (D,)
        c_t = phase * c               # state coefficients at time t
        exp_t[idx] = np.vdot(c_t, O_eig @ c_t)

    return exp_t


def observable_dynamics_full_hilbert(
    H: PauliSum,
    O: PauliSum,
    psi0: np.ndarray,
    times: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """
    Direct spectral time evolution with the full Hilbert-space Hamiltonian H.

    Parameters
    ----------
    H : PauliSum
        Hamiltonian (in original basis).
    O : PauliSum
        Observable (in original basis).
    psi0 : (D,) complex ndarray
        Initial statevector in original basis.
    times : array-like
        Times at which to evaluate ⟨O(t)⟩.

    Returns
    -------
    exp_t : (len(times),) complex ndarray
        Expectation values ⟨O(t)⟩.
    """
    times = np.asarray(times, dtype=float)
    dims = np.asarray(H.dimensions, dtype=int)
    D = int(np.prod(dims))
    assert psi0.shape == (D,)

    # Build full Hilbert-space matrices
    H_h = H.to_hilbert_space().toarray()  # (D, D)
    O_h = O.to_hilbert_space().toarray()  # (D, D)

    # Diagonalise H
    e_vals, V = np.linalg.eigh(H_h)        # H = V diag(evals) V†

    # Expand initial state in eigenbasis
    c0 = V.conj().T @ psi0

    exp_t = np.empty(times.shape[0], dtype=np.complex128)
    for idx, t in enumerate(times):
        phase = np.exp(-1j * e_vals * t)
        c_t = phase * c0
        psi_t = V @ c_t
        exp_t[idx] = np.vdot(psi_t, O_h @ psi_t)

    return exp_t


def run_joint_eigen_dynamics_example():
    # Model parameters
    N = 10
    J = 1.0
    h = 0.5

    # 1) Build Ising Hamiltonian H
    H = ising_chain_hamiltonian(N, J, h, periodic=True)

    # 2) Find heuristic Clifford symmetry F and its block decomposition S, T
    F = heuristic_clifford_symmetry(N)  # some Gate-like object with symplectic + phase
    S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)

    # Phases for S and T
    h_S, h_T = clifford_phase_decomposition(
        F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
    )

    S_gate = Gate('S', F.qudit_indices, S_symp, F.dimensions, h_S)
    T_gate = Gate('T', F.qudit_indices, T_symp, F.dimensions, h_T)

    # 3) Build a simple observable: X on qudit 0
    obs_tableau = np.zeros((1, 2 * N), dtype=int)
    obs_tableau[0, 0] = 1  # X on qudit 0
    O_x = PauliSum.from_tableau(obs_tableau, weights=[1.0], dimensions=[2] * N)

    # 4) Get block structure from S (in symplectic space)
    blocks = block_indexes(S_gate.symplectic)
    print("Blocks from S:", blocks)

    # 5) Choose initial product state |ψ0>
    psi0, dims = product_state_ising(N, kind="x_plus")  # |+...+>

    # 6) Time grid
    times = np.linspace(0.0, 10.0, 201)

    # 7) Compute ⟨X_0(t)⟩ using joint eigenbasis of H' and S
    exp_t = observable_dynamics_from_joint_eigenbasis(
        H=H,
        O=O_x,
        S_gate=S_gate,
        T_gate=T_gate,
        blocks=blocks,
        psi0=psi0,
        times=times,
    )

    # 8) Print / inspect results
    print("⟨X_0(t)⟩ at a few sample times:")
    for t, val in zip(times[::50], exp_t[::50]):
        print(f" t = {t:.3f},  <X_0(t)> = {val.real:.6f} + {val.imag:.2e}i")


def run_krylov_symmetry_example():
    """
    Example: exploit a Clifford symmetry S of an Ising chain in a Krylov
    time-evolution scheme.

    Steps:
      1) Build Ising Hamiltonian H.
      2) Find a Clifford symmetry F and its block decomposition S, T.
      3) Build Clifford gates S_gate, T_gate.
      4) Choose observable O_x = X_0 and initial product state |+...+>.
      5) Project the initial state into the +1 symmetry sector of S in Hilbert space
         (for reference).
      6) Compute exact sector dynamics <X_0(t)> in +1 sector via full Hilbert matrices.
      7) Compute the same sector dynamics using symmetry-aware Krylov evolution:
           krylov_observable_dynamics_symmetry(..., symmetry_eigval=+1.0).
      8) Plot and compare.
    """
    # --- 1) Ising chain Hamiltonian H ---
    Ns = [10, 11, 12, 13, 14]  # , 15, 16
    fig, ax = plt.subplots(figsize=(7, 4))
    times_symmetry = []
    times_no_symmetry = []
    for N in Ns:
        J = 1.0
        h = 0.5
        H = ising_chain_hamiltonian(N, J, h, periodic=True)

        # --- 2) Symmetry via heuristic Clifford F, block decompose into S, T ---
        F = heuristic_clifford_symmetry(N)
        S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)

        # --- 3) Build Clifford gates S and T with correct phases ---
        h_S, h_T = clifford_phase_decomposition(
            F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
        )

        S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
        T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

        # --- 4) Observable O_x = X on qubit 0 ---
        obs_tableau = np.zeros((1, 2 * N), dtype=int)
        obs_tableau[0, 0] = 1  # X on qudit 0
        O_x = PauliSum.from_tableau(obs_tableau, weights=[1.0], dimensions=[2] * N)

        # Initial product state |ψ0> = |+...+>
        psi0, dims = product_state_ising(N, kind="x_plus")
        dims = np.asarray(dims, dtype=int)
        D = int(np.prod(dims))

        # Time grid
        times = np.linspace(0.0, 5.0, 201)

        # -------------------------------------------------------------------------
        # 5) Reference: project |ψ0> into +1 sector of S and compute exact dynamics
        # -------------------------------------------------------------------------

        # Build S as a full Hilbert-space unitary
        C_S = gate_to_circuit(S_gate)
        U_S = C_S.unitary().toarray()   # dimension 2^N x 2^N (for qubits)

        # Projector onto +1 sector: P_+ = (I + U_S)/2
        I_full = np.eye(D, dtype=np.complex128)
        P_plus = 0.5 * (I_full + U_S)

        psi_plus = P_plus @ psi0
        norm_plus = np.linalg.norm(psi_plus)
        if norm_plus < 1e-12:
            raise RuntimeError(
                "Projection onto +1 sector of S is (numerically) zero. "
                "Choose a different initial state or symmetry sector."
            )
        psi_plus /= norm_plus
        t0 = time()
        # Exact sector dynamics using full Hilbert matrices
        exp_exact_sector = observable_dynamics_full_hilbert(
            H=H,
            O=O_x,
            psi0=psi_plus,
            times=times)
        t1 = time()
        times_no_symmetry.append(t1 - t0)
        print(f"Exact +1 sector dynamics computed in {t1 - t0:.3f} seconds.")

        # -------------------------------------------------------------------------
        # 6) Symmetry-aware Krylov dynamics in the +1 sector
        # -------------------------------------------------------------------------

        # This will:
        #   - conjugate H, O by T: H' = T^{-1} H T, O' = T^{-1} O T
        #   - transform the state to the T-basis
        #   - project the T-basis state into the chosen symmetry sector of S
        #   - run Lanczos/Krylov in that sector using PauliSum matvecs only.
        t0 = time()
        exp_krylov_sym = krylov_observable_dynamics_symmetry(
            H=H,
            O=O_x,
            T_gate=T_gate,
            S_gate=S_gate,
            psi0=psi0,
            times=times,
            m_max=128,       # Krylov dimension; 128 is generous for N=10
            symmetry_eigval=+1.0,
        )
        t1 = time()
        times_symmetry.append(t1 - t0)
        print(f"Symmetry aware Krylov +1 sector dynamics computed in {t1 - t0:.3f} seconds.")

        # -------------------------------------------------------------------------
        # 7) Plot comparison
        # -------------------------------------------------------------------------

        ax.plot(
            times,
            exp_exact_sector.real,
            label="Exact (+1 sector, full Hilbert)",
            linestyle="-",
        )
        ax.plot(
            times,
            exp_krylov_sym.real,
            label="Krylov (+1 sector, symmetry-exploiting)",
            linestyle="--",
        )
        ax.set_xlabel("$t$")
        ax.set_ylabel("Re$\\langle X_0(t) \\rangle$")
        ax.legend()
        plt.tight_layout()

        # Optional: print max difference as a sanity check
        max_diff = np.max(np.abs(exp_krylov_sym - exp_exact_sector))
        print(f"Max |Krylov_sym - exact_sector| over times = {max_diff:.3e}")

    fig, ax2 = plt.subplots(figsize=(7, 4))
    ax2.plot(Ns, times_no_symmetry, marker='o', label='No symmetry')
    ax2.plot(Ns, times_symmetry, marker='o', label='With symmetry')
    ax2.set_xlabel('Number of spins N')
    ax2.set_ylabel('Time (s)')

    plt.show()


def symmetry_aware_testing():
    fig, ax = plt.subplots(figsize=(7, 4))
    Ns = [25]
    for N in Ns:
        print(f"Running symmetry-aware Krylov example for N={N}...")
        t0 = time()
        J = 1.0
        h = 0.5
        H = ising_chain_hamiltonian(N, J, h, periodic=True)

        # --- 2) Symmetry via heuristic Clifford F, block decompose into S, T ---
        F = heuristic_clifford_symmetry(N)
        S_symp, T_symp = block_decompose_optimal(F.symplectic, 2)

        # --- 3) Build Clifford gates S and T with correct phases ---
        h_S, h_T = clifford_phase_decomposition(
            F.symplectic, F.phase_vector, S_symp, T_symp, int(H.lcm)
        )

        S_gate = Gate("S", F.qudit_indices, S_symp, F.dimensions, h_S)
        T_gate = Gate("T", F.qudit_indices, T_symp, F.dimensions, h_T)

        # --- 4) Observable O_x = X on qubit 0 ---
        obs_tableau = np.zeros((1, 2 * N), dtype=int)
        obs_tableau[0, 0] = 1  # X on qudit 0
        O_x = PauliSum.from_tableau(obs_tableau, weights=[1.0], dimensions=[2] * N)

        # Initial product state |ψ0> = |+...+>
        psi0, dims = product_state_ising(N, kind="x_plus")
        dims = np.asarray(dims, dtype=int)

        # Time grid
        times = np.linspace(0.0, 5.0, 101)

        # -------------------------------------------------------------------------
        # 6) Symmetry-aware Krylov dynamics in the +1 sector
        # -------------------------------------------------------------------------

        # This will:
        #   - conjugate H, O by T: H' = T^{-1} H T, O' = T^{-1} O T
        #   - transform the state to the T-basis
        #   - project the T-basis state into the chosen symmetry sector of S
        #   - run Lanczos/Krylov in that sector using PauliSum matvecs only.
        exp_krylov_sym = krylov_observable_dynamics_symmetry(
            H=H,
            O=O_x,
            T_gate=T_gate,
            S_gate=S_gate,
            psi0=psi0,
            times=times,
            m_max=128,       # Krylov dimension; 128 is generous for N=10
            symmetry_eigval=+1.0,
        )
        t1 = time()
        print(f"Symmetry aware Krylov +1 sector dynamics N = {N} computed in {t1 - t0:.3f} seconds.")
        ax.plot(times, exp_krylov_sym.real, label=f"N={N}")
    ax.set_xlabel("t")
    ax.set_ylabel("$Re \\langle X_0 (t)\\rangle$")
    ax.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # run_joint_eigen_dynamics_example()
    # run_compare_joint_vs_full()
    run_krylov_symmetry_example()
    # symmetry_aware_testing()
