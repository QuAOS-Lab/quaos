import numpy as np
import scipy.sparse as sp
from typing import List, Tuple, Sequence
from sympleq.core.paulis.utils import apply_paulisum_to_state_dense
from sympleq.core.circuits import Circuit, gate_to_circuit
from sympleq.core.circuits.gates import Gate
from sympleq.core.paulis import PauliSum
from sympleq.core.states.state import State
from sympleq.core.symmetries.block_decomposition import block_indexes
from dataclasses import dataclass
from itertools import product




@dataclass
class DenseSectorBasis:
    """
    Dense basis for a symmetry sector.

    basis: shape (D, D_sec), columns are an orthonormal basis of the sector.
    """
    basis: np.ndarray  # (D, D_sec)

    @property
    def dim_full(self) -> int:
        return self.basis.shape[0]

    @property
    def dim_sector(self) -> int:
        return self.basis.shape[1]

    def project_state(self, psi_full: np.ndarray) -> np.ndarray:
        """
        Project a full Hilbert state |psi> onto sector coordinates:
            phi = B^\dagger |psi>.
        """
        psi_full = np.asarray(psi_full, dtype=np.complex128).reshape(self.dim_full)
        return self.basis.conj().T @ psi_full

    def lift_state(self, phi_sec: np.ndarray) -> np.ndarray:
        """
        Lift sector coordinates |phi> back to full Hilbert space:
            |psi> = B |phi>.
        """
        phi_sec = np.asarray(phi_sec, dtype=np.complex128).reshape(self.dim_sector)
        return self.basis @ phi_sec


def _local_symmetry_block_eigendecomposition(
    S_gate: Gate,
    block: List[int],
    dims: Sequence[int] | np.ndarray,
    tol: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Diagonalise S restricted to a single block of qudits.

    Parameters
    ----------
    S_gate : Gate
        Clifford gate implementing the symmetry S on some set of qudits.
    block : list[int]
        Qudit indices belonging to this block.
    dims : sequence[int]
        Local dimensions of all qudits (global ordering).
    tol : float
        Numerical tolerance for sanity checks.

    Returns
    -------
    eigvals_block : (d_block,) complex ndarray
        Eigenvalues of S|_block.
    eigvecs_block : (d_block, d_block) complex ndarray
        Columns are eigenvectors in the block Hilbert space.
    """
    # Build a circuit from S_gate and restrict to the block
    C_S = gate_to_circuit(S_gate)
    C_S_block = C_S.local_circuit(block)

    # Unitary for S restricted to this block
    U_block = C_S_block.unitary().toarray()
    eigvals, eigvecs = np.linalg.eig(U_block)

    # Optional: normalise eigenvectors explicitly
    for j in range(eigvecs.shape[1]):
        v = eigvecs[:, j]
        nrm = np.linalg.norm(v)
        if nrm < tol:
            raise RuntimeError("Found near-zero eigenvector in S block eigendecomposition.")
        eigvecs[:, j] = v / nrm

    return eigvals, eigvecs


def build_symmetry_sector_basis(
    S_gate: Gate,
    blocks: List[List[int]],
    dims: Sequence[int] | np.ndarray,
    target_eigval: complex = +1.0,
    tol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Construct an explicit orthonormal basis for the symmetry sector
        H_{target_eigval} = { |phi> : S |phi> = target_eigval |phi> }.

    We use the fact that S decomposes into independent blocks. For each block
    we diagonalise S_block, then build global eigenvectors as tensor products
    of block eigenvectors, and retain only those with global eigenvalue
    matching target_eigval up to `tol`.

    Parameters
    ----------
    S_gate : Gate
        Clifford implementing S in the computational basis.
    blocks : list[list[int]]
        List of blocks of qudit indices, e.g. [[0,1], [2,3], ...].
    dims : sequence[int]
        Local dimensions for all qudits in global order.
    target_eigval : complex
        Desired global eigenvalue (e.g. +1 or -1 for Z2 symmetry).
    tol : float
        Tolerance for matching eigenvalues.

    Returns
    -------
    B_sector : (D, D_sector) complex ndarray
        Columns form an orthonormal basis of the target symmetry sector,
        expressed in the global computational basis.
    sector_eigvals : (D_sector,) complex ndarray
        The corresponding S eigenvalue for each basis vector (should all
        equal target_eigval up to `tol`).
    """
    dims = np.asarray(dims, dtype=int)
    n = len(dims)
    D = int(np.prod(dims))

    # 1) Local eigenpairs for each block
    block_eigvals: List[np.ndarray] = []
    block_eigvecs: List[np.ndarray] = []
    block_dims: List[int] = []

    for block in blocks:
        ev, evec = _local_symmetry_block_eigendecomposition(S_gate, block, dims, tol=tol)
        block_eigvals.append(ev)
        block_eigvecs.append(evec)
        # local Hilbert dim for this block
        d_block = int(np.prod([dims[i] for i in block]))
        if d_block != evec.shape[0]:
            raise ValueError("Block dimension mismatch in S block eigendecomposition.")
        block_dims.append(d_block)

    # 2) Build all tensor-product eigenvectors and select those with global eigenvalue = target_eigval
    # For B blocks, index tuples index_j choose an eigenvector in each block.
    num_blocks = len(blocks)
    index_ranges = [range(len(block_eigvals[b])) for b in range(num_blocks)]

    # meshgrid over eigenvector choices per block
    all_choices = np.array(
        np.meshgrid(*index_ranges, indexing="ij")
    ).reshape(num_blocks, -1).T  # shape: (n_choices, num_blocks)

    basis_vectors = []
    sector_eigvals = []

    for choice in all_choices:
        lam_global = 1.0 + 0.0j
        # Start from a 1-dim vector and tensor up
        v_global = np.array([1.0 + 0.0j])

        for b, idx in enumerate(choice):
            lam_b = block_eigvals[b][idx]
            lam_global *= lam_b
            v_block = block_eigvecs[b][:, idx]  # shape (d_block,)

            # Kronecker with current v_global
            v_global = np.kron(v_global, v_block)

        # Check if this global eigenvector belongs to the target sector
        if np.abs(lam_global - target_eigval) < tol:
            # Normalise and store
            nrm = np.linalg.norm(v_global)
            if nrm < tol:
                continue
            basis_vectors.append(v_global / nrm)
            sector_eigvals.append(lam_global)

    if len(basis_vectors) == 0:
        raise RuntimeError(
            f"No eigenvectors found in S sector with eigenvalue {target_eigval}."
        )

    B_sector = np.column_stack(basis_vectors)  # (D, D_sector)
    sector_eigvals = np.asarray(sector_eigvals, dtype=complex)

    # Orthonormalise columns of B_sector (just in case)
    Q, R = np.linalg.qr(B_sector)
    # Ensure phases are consistent (optional)
    diagR = np.diag(R)
    phases = np.where(diagR == 0, 1.0, diagR / np.abs(diagR))
    Q = Q * phases  # Fix arbitrary column phases

    # Recompute eigenvalues using Q (they should still be target_eigval)
    # but we keep the original `sector_eigvals` which are all ~= target.

    return Q, sector_eigvals


def reduce_operator_to_sector(
    A_h: np.ndarray,
    B_sector: np.ndarray,
) -> np.ndarray:
    """
    Reduce a full Hilbert-space operator A_h to the symmetry sector
    spanned by columns of B_sector, via

        A_sector = B_sector^\dagger A_h B_sector.

    Parameters
    ----------
    A_h : (D, D) ndarray
        Operator in the global computational basis.
    B_sector : (D, D_sec) ndarray
        Orthonormal basis of the sector.

    Returns
    -------
    A_sector : (D_sec, D_sec) ndarray
        Operator restricted to the sector.
    """
    return B_sector.conj().T @ A_h @ B_sector


def project_state_to_sector(
    psi: np.ndarray,
    B_sector: np.ndarray,
) -> np.ndarray:
    """
    Represent a global state |psi> in the symmetry sector basis:

        |psi_sec> = B_sector^\dagger |psi>

    Parameters
    ----------
    psi : (D,) ndarray
        Statevector in the global computational basis.
    B_sector : (D, D_sec) ndarray
        Orthonormal basis of the symmetry sector.

    Returns
    -------
    psi_sec : (D_sec,) ndarray
        Coordinates of |psi> in the sector basis.
    """
    return B_sector.conj().T @ psi

def krylov_observable_dynamics_symmetry_reduced(
    H: PauliSum,
    O: PauliSum,
    T_gate: Gate,
    S_gate: Gate,
    psi0: np.ndarray,
    times: Sequence[float] | np.ndarray,
    symmetry_eigval: complex = 1.0,
    m_max: int = 64,
) -> np.ndarray:
    """
    Symmetry-reduced Krylov evolution of <O(t)> using H, O and a Clifford symmetry S
    with block structure, without ever forming full Hilbert matrices for H' or O'.

    Steps:
      1. Conjugate H, O by T^{-1} to get H', O' that commute with S.
      2. Build a dense basis B_λ for the symmetry sector eigval=symmetry_eigval
         using the block structure of S.
      3. Transform the initial state |psi0> to the H' frame, |psi0'> = U_T^\dagger |psi0>,
         and project it into the sector, obtaining |psi0_sec>.
      4. Reduce H' and O' to the sector: H_sec = B_λ^\dagger H' B_λ, O_sec likewise.
      5. Run dense-matrix Lanczos/Krylov on (H_sec, O_sec, |psi0_sec>) to get <O(t)>.
    """
    times = np.asarray(times, dtype=float)
    dims = np.asarray(H.dimensions, dtype=int).reshape(-1)
    D = int(np.prod(dims))
    psi0 = np.asarray(psi0, dtype=np.complex128).reshape(D)

    # 1) H' = T^{-1} H T, O' = T^{-1} O T
    T_inv_gate = T_gate.inv()
    H_prime = T_inv_gate.act(H)
    O_prime = T_inv_gate.act(O)

    # 2) Build sector basis B_λ from block structure of S
    blocks = block_indexes(S_gate.symplectic)
    sector_basis = build_dense_sector_basis(
        S_gate=S_gate,
        blocks=blocks,
        dims=dims,
        target_eigval=symmetry_eigval,
    )

    # 3) Transform initial state to H' frame: |psi0'> = U_T^\dagger |psi0>
    state0 = State(psi0, dims)
    C_T_inv = gate_to_circuit(T_inv_gate)
    state0_prime = C_T_inv.act_on_state(state0)
    psi0_prime = state0_prime.as_array()

    # Project |psi0'> into the symmetry sector:
    psi0_sec = sector_basis.project_state(psi0_prime)
    norm_sec = np.linalg.norm(psi0_sec)
    if norm_sec < 1e-14:
        raise ValueError(
            "Initial state has (near) zero overlap with the requested symmetry sector."
        )
    psi0_sec = psi0_sec / norm_sec

    # 4) Reduced operators in the sector H_sec, O_sec
    H_sec = reduce_paulisum_to_sector_dense(H_prime, sector_basis)
    O_sec = reduce_observable_to_sector_dense(O_prime, sector_basis)

    # 5) Krylov in the reduced sector
    return _krylov_observable_from_dense_matrices(
        H_mat=H_sec,
        O_mat=O_sec,
        psi0=psi0_sec,
        times=times,
        m_max=m_max,
    )



# ---------- 1. Local eigen-decomposition of S on a single block ----------

def _local_block_eigendecomp(
    S_gate: Gate,
    block: Sequence[int],
    tol: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Diagonalise the restriction of S_gate to a given block.

    Parameters
    ----------
    S_gate : Gate
        Symmetry gate S acting on all qudits.
    block : list[int]
        Indices of qudits in this block.
    dims : Sequence[int]
        Local dimensions of all qudits (len = N).
    tol : float
        Numerical tolerance (unused here but kept for future use).

    Returns
    -------
    evals_block : (d_block,) complex ndarray
        Eigenvalues of U_S restricted to the block.
    evecs_block : (d_block, d_block) complex ndarray
        Eigenvectors as columns, in the local Hilbert space for the block.
    """
    # Build circuit for S, then restrict to block
    C_S = gate_to_circuit(S_gate)
    C_block: Circuit = C_S.local_circuit(list(block))

    # Local unitary U_S_block
    U_S_block = C_block.unitary().toarray()
    evals_block, evecs_block = np.linalg.eig(U_S_block)

    # Optional: normalise columns (should already be unitary)
    for j in range(evecs_block.shape[1]):
        norm = np.linalg.norm(evecs_block[:, j])
        if norm > tol:
            evecs_block[:, j] /= norm

    return evals_block, evecs_block


# ---------- 2. Tensorised sector basis over all blocks ----------

@dataclass
class TensorSectorBasis:
    """
    Tensor-product symmetry sector basis for a given eigenvalue of S.

    The idea:
      - For each block b, we diagonalise U_S on that block -> (λ_{b,r}, |φ_{b,r}>).
      - Global eigenvectors of S are tensor products of chosen local eigenvectors.
      - We keep only those global products with total eigenvalue = target_eigval.

    This class does NOT store a full D x D_sector dense matrix. Instead, it stores:
      - per-block eigenvalues/eigenvectors,
      - a list of global multi-indices labeling sector basis vectors.
    """
    dims: np.ndarray                      # (N,)
    blocks: List[List[int]]               # list of blocks
    local_evals: List[np.ndarray]         # per block: (d_block,)
    local_evecs: List[np.ndarray]         # per block: (d_block, d_block) columns
    sector_eigval: complex                # target eigenvalue λ
    sector_multi_indices: List[Tuple[int, ...]]  # list of r = (r_b) for sector basis

    @property
    def n_blocks(self) -> int:
        return len(self.blocks)

    @property
    def dim_sector(self) -> int:
        return len(self.sector_multi_indices)

    @property
    def dim_total(self) -> int:
        return int(np.prod(self.dims))

    @classmethod
    def from_gate_and_blocks(
        cls,
        S_gate: Gate,
        blocks: List[List[int]],
        dims: Sequence[int] | np.ndarray,
        target_eigval: complex = 1.0,
        tol: float = 1e-8,
    ) -> "TensorSectorBasis":
        """
        Build a tensorised sector basis for S_gate at eigenvalue target_eigval.

        Steps:
          1. For each block, diagonalise U_S_block -> (evals_b, evecs_b).
          2. Enumerate all multi-indices (r_b) over blocks.
          3. Keep those for which prod_b evals_b[r_b] ≈ target_eigval.
        """
        dims = np.asarray(dims, dtype=int).reshape(-1)
        local_evals: List[np.ndarray] = []
        local_evecs: List[np.ndarray] = []

        for block in blocks:
            evals_b, evecs_b = _local_block_eigendecomp(S_gate, block)
            local_evals.append(evals_b)
            local_evecs.append(evecs_b)

        # Enumerate all combinations of local eigenvector indices
        sector_multi_indices: List[Tuple[int, ...]] = []
        for idx_tuple in product(
            *[range(len(evals_b)) for evals_b in local_evals]
        ):
            lam = 1.0 + 0j
            for b, r_b in enumerate(idx_tuple):
                lam *= local_evals[b][r_b]
            if abs(lam - target_eigval) < tol:
                sector_multi_indices.append(tuple(idx_tuple))

        if not sector_multi_indices:
            raise ValueError(
                f"No global eigenvectors found for eigenvalue {target_eigval} "
                f"with the given blocks."
            )

        return cls(
            dims=dims,
            blocks=[list(b) for b in blocks],
            local_evals=local_evals,
            local_evecs=local_evecs,
            sector_eigval=target_eigval,
            sector_multi_indices=sector_multi_indices,
        )

    # ----- utilities for constructing full Hilbert vectors from sector basis -----

    def _block_dim(self, block_index: int) -> int:
        block = self.blocks[block_index]
        return int(np.prod(self.dims[block]))

    def basis_vector_full(self, k: int) -> np.ndarray:
        """
        Construct the k-th sector basis vector as a full Hilbert-space vector
        in computational basis.

        This is used for:
          - small-N tests,
          - building reduced operators H_sec = B† H B,
          - projecting arbitrary states into the sector.

        Shape: (D,)
        """
        if k < 0 or k >= self.dim_sector:
            raise IndexError("Sector basis index out of range.")

        # For each block, take the chosen local eigenvector in that block's Hilbert space,
        # then tensor across blocks.
        idx_tuple = self.sector_multi_indices[k]

        # local_evecs[b]: (d_block_b, d_block_b)
        # pick column idx_tuple[b] for block b
        psi_block_list: List[np.ndarray] = []
        for b, r_b in enumerate(idx_tuple):
            v_b = self.local_evecs[b][:, r_b]   # (d_block_b,)
            psi_block_list.append(v_b)

        psi_full = psi_block_list[0]
        for v_b in psi_block_list[1:]:
            psi_full = np.kron(psi_full, v_b)

        # Hilbert space dimension check
        if psi_full.shape[0] != self.dim_total:
            raise ValueError(
                f"Constructed sector basis vector has dim {psi_full.shape[0]}, "
                f"expected {self.dim_total}."
            )

        return psi_full

    def project_state_to_sector(self, psi: np.ndarray) -> np.ndarray:
        """
        Project a full Hilbert-space state |psi> into the sector and return
        coordinates in the sector basis (length = dim_sector).

        c_k = <φ_k | psi>, where |φ_k> is the k-th sector basis vector.
        """
        psi = np.asarray(psi, dtype=np.complex128).reshape(self.dim_total)

        coeffs = np.empty(self.dim_sector, dtype=np.complex128)
        for k in range(self.dim_sector):
            phi_k = self.basis_vector_full(k)
            coeffs[k] = np.vdot(phi_k, psi)

        # We do not re-normalise here; the caller can decide.
        return coeffs

    def lift_sector_vector(self, coeffs: np.ndarray) -> np.ndarray:
        """
        Embed a sector-coordinate vector c (length = dim_sector) into the full
        Hilbert space:
            |psi> = sum_k c_k |φ_k>.
        """
        coeffs = np.asarray(coeffs, dtype=np.complex128).reshape(self.dim_sector)
        psi_full = np.zeros(self.dim_total, dtype=np.complex128)

        for k, c_k in enumerate(coeffs):
            if abs(c_k) < 1e-15:
                continue
            phi_k = self.basis_vector_full(k)
            psi_full += c_k * phi_k

        return psi_full


def reduce_paulisum_to_sector_tensor(
    H_prime: PauliSum,
    sector_basis: TensorSectorBasis,
) -> np.ndarray:
    """
    Build the reduced Hamiltonian H_sec = B† H' B in the sector basis defined
    by `sector_basis`, where columns of B are |φ_k> (sector basis vectors).
    
    This uses PauliSum matvec and TensorSectorBasis.basis_vector_full(k),
    so no full H' matrix is needed. Complexity is O(D * dim_sector^2), so
    this is still only for small/intermediate-N, but sets the structure
    for later tensor/MPS implementations.
    """
    D_sec = sector_basis.dim_sector
    H_sec = np.zeros((D_sec, D_sec), dtype=np.complex128)

    for j in range(D_sec):
        # column j of H_sec is B† H' |φ_j>
        phi_j = sector_basis.basis_vector_full(j)
        Hphi_j = apply_paulisum_to_state_dense(H_prime, phi_j)
        for i in range(D_sec):
            phi_i = sector_basis.basis_vector_full(i)
            H_sec[i, j] = np.vdot(phi_i, Hphi_j)

    return H_sec


def reduce_observable_to_sector_tensor(
    O_prime: PauliSum,
    sector_basis: TensorSectorBasis,
) -> np.ndarray:
    """
    Build the reduced observable O_sec = B† O' B in the same sector basis.
    """
    D_sec = sector_basis.dim_sector
    O_sec = np.zeros((D_sec, D_sec), dtype=np.complex128)

    for j in range(D_sec):
        phi_j = sector_basis.basis_vector_full(j)
        Ophi_j = apply_paulisum_to_state_dense(O_prime, phi_j)
        for i in range(D_sec):
            phi_i = sector_basis.basis_vector_full(i)
            O_sec[i, j] = np.vdot(phi_i, Ophi_j)

    return O_sec


def _krylov_observable_from_dense_matrices(
    H_mat: np.ndarray,
    O_mat: np.ndarray,
    psi0: np.ndarray,
    times: np.ndarray,
    m_max: int,
) -> np.ndarray:
    """
    Lanczos/Krylov evolution of <O(t)> using dense H_mat and O_mat,
    starting from psi0 in the same space.
    """
    times = np.asarray(times, dtype=float)
    H_mat = np.asarray(H_mat, dtype=np.complex128)
    O_mat = np.asarray(O_mat, dtype=np.complex128)

    n = H_mat.shape[0]
    assert H_mat.shape == (n, n)
    assert O_mat.shape == (n, n)

    psi0 = np.asarray(psi0, dtype=np.complex128).reshape(n)
    norm0 = np.linalg.norm(psi0)
    if norm0 < 1e-14:
        raise ValueError("Initial sector state has (near) zero norm.")
    v0 = psi0 / norm0

    m_max = int(m_max)
    if m_max <= 0:
        raise ValueError("m_max must be positive")

    # Lanczos
    V = np.zeros((n, m_max), dtype=np.complex128)
    alpha = np.zeros(m_max, dtype=float)
    beta = np.zeros(m_max - 1, dtype=float)

    V[:, 0] = v0
    w = H_mat @ v0
    alpha[0] = np.vdot(v0, w).real
    w = w - alpha[0] * v0

    k = 1
    for j in range(1, m_max):
        beta_jm1 = np.linalg.norm(w)
        if beta_jm1 < 1e-14:
            k = j
            break
        beta[j - 1] = beta_jm1
        vj = w / beta_jm1
        V[:, j] = vj

        w = H_mat @ vj
        alpha[j] = np.vdot(vj, w).real
        w = w - alpha[j] * vj - beta[j - 1] * V[:, j - 1]
        k = j + 1

    alpha = alpha[:k]
    beta = beta[: max(0, k - 1)]
    V = V[:, :k]

    # Tridiagonal T_k
    T_k = np.diag(alpha)
    if k > 1:
        T_k += np.diag(beta, 1) + np.diag(beta, -1)

    # Diagonalise T_k
    evals, U = np.linalg.eigh(T_k)
    evals = evals.astype(float)
    U_dag = U.conj().T

    # Observable in Krylov basis K_k
    O_K = np.zeros((k, k), dtype=np.complex128)
    for j in range(k):
        Oj = O_mat @ V[:, j]
        for i in range(k):
            O_K[i, j] = np.vdot(V[:, i], Oj)

    # Transform O_K to eigenbasis of T_k
    O_K_eig = U_dag @ O_K @ U

    # Initial Krylov vector is e1
    e1 = np.zeros(k, dtype=np.complex128)
    e1[0] = 1.0
    e1_eig = U_dag @ e1

    # Time evolution in eigenbasis
    expectations = np.empty(times.shape[0], dtype=np.complex128)
    for idx, t in enumerate(times):
        phase = np.exp(-1j * evals * t)  # exp(-i E t)
        # O(t) in eigenbasis: O_ij e^{i(λ_i - λ_j)t}
        O_t_eig = O_K_eig * (phase[None, :] * phase.conj()[:, None])
        expectations[idx] = np.vdot(e1_eig, O_t_eig @ e1_eig)

    return expectations


def reduce_observable_to_sector_dense(
    O: PauliSum,
    sector_basis: DenseSectorBasis,
) -> np.ndarray:
    """
    Reduced observable O_sec = B^\dagger O B, using only PauliSum matvecs.
    """
    D_sec = sector_basis.dim_sector
    O_sec = np.zeros((D_sec, D_sec), dtype=np.complex128)

    for j in range(D_sec):
        e_j = np.zeros(D_sec, dtype=np.complex128)
        e_j[j] = 1.0
        psi_j = sector_basis.lift_state(e_j)
        y_full = apply_paulisum_to_state_dense(O, psi_j)
        O_sec[:, j] = sector_basis.project_state(y_full)

    return O_sec

def reduce_paulisum_to_sector_dense(
    H: PauliSum,
    sector_basis: DenseSectorBasis,
) -> np.ndarray:
    """
    Build the reduced Hamiltonian H_sec = B^\dagger H B in the symmetry sector
    defined by `sector_basis`, using only PauliSum matvecs:
        H_sec[:, j] = B^\dagger ( H (B e_j) ).
    """
    D_sec = sector_basis.dim_sector
    H_sec = np.zeros((D_sec, D_sec), dtype=np.complex128)

    for j in range(D_sec):
        e_j = np.zeros(D_sec, dtype=np.complex128)
        e_j[j] = 1.0
        psi_j = sector_basis.lift_state(e_j)               # full vector
        y_full = apply_paulisum_to_state_dense(H, psi_j)   # full vector
        H_sec[:, j] = sector_basis.project_state(y_full)   # sector coords

    return H_sec


def build_dense_sector_basis(
    S_gate: Gate,
    blocks: list[list[int]],
    dims: Sequence[int] | np.ndarray,
    target_eigval: complex = 1.0,
    tol: float = 1e-8,
) -> DenseSectorBasis:
    """
    Build a DenseSectorBasis for the eigenspace of S_gate with eigenvalue target_eigval.

    Uses the block decomposition of S_gate (given by `blocks`), diagonalises each
    block unitary U_block, and forms global eigenvectors as tensor products
    of local eigenvectors. We keep only those global eigenvectors whose
    eigenvalue product matches target_eigval up to tolerance.
    """
    dims = np.asarray(dims, dtype=int).reshape(-1)
    N = len(dims)
    D = int(np.prod(dims))

    from sympleq.core.circuits import gate_to_circuit

    # Local block eigenpairs
    local_eigs: list[tuple[np.ndarray, np.ndarray]] = []
    C_S = gate_to_circuit(S_gate)

    for block in blocks:
        # U_block acts only on this block; represent in Hilbert space
        C_block = C_S.local_circuit(block)
        U_block = C_block.unitary().toarray()
        evals, vecs = np.linalg.eig(U_block)

        # Normalise eigenvectors
        for j in range(vecs.shape[1]):
            v = vecs[:, j]
            nrm = np.linalg.norm(v)
            if nrm > 0:
                vecs[:, j] = v / nrm
        local_eigs.append((evals, vecs))

    # Build global eigenvectors as tensor products of block eigenvectors
    basis_vecs: list[np.ndarray] = []
    # Indices over local eigenvectors in each block
    index_grids = np.array(
        np.meshgrid(*[np.arange(len(evals_b)) for (evals_b, _) in local_eigs]),
        dtype=int,
    ).reshape(len(local_eigs), -1).T

    for choice in index_grids:
        lam_global = 1.0 + 0.0j
        vec_global = np.array([1.0 + 0.0j])
        total_dim = 1

        for (evals_b, vecs_b), block_idx in zip(local_eigs, range(len(local_eigs))):
            j = choice[block_idx]
            lam_b = evals_b[j]
            v_b = vecs_b[:, j]  # shape (d_block,)

            lam_global *= lam_b
            # Kronecker product to build full vector
            vec_global = np.kron(vec_global, v_b)
            total_dim *= v_b.shape[0]

        # Filter by eigenvalue
        if np.abs(lam_global - target_eigval) < tol:
            # Normalise (should already be norm 1, but just in case)
            nrm = np.linalg.norm(vec_global)
            if nrm > 0:
                vec_global = vec_global / nrm
            # Embed in full Hilbert space (the Kronecker we used already has full dimension)
            assert total_dim == D, "Global vector dimension mismatch."
            basis_vecs.append(vec_global)

    if not basis_vecs:
        raise ValueError(f"No eigenvectors found for eigenvalue {target_eigval}.")

    B = np.column_stack(basis_vecs)  # (D, D_sec)
    # Orthonormalise to get a clean basis
    Q, _ = np.linalg.qr(B)
    return DenseSectorBasis(Q)
