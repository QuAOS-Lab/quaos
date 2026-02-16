from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import List, Sequence, Tuple

import numpy as np

from sympleq.core.circuits import Circuit, gate_to_circuit
from sympleq.core.circuits.gates import Gate
from sympleq.core.paulis import PauliSum
from sympleq.core.paulis.utils import apply_paulisum_to_state_dense
from sympleq.core.states.state import State
from sympleq.core.symmetries.block_decomposition import block_indexes


# ======================================================================
# 0. Core Krylov / Lanczos engine on a *dense* reduced space
# ======================================================================


def _krylov_observable_from_dense_matrices(
    H_mat: np.ndarray,
    O_mat: np.ndarray,
    psi0: np.ndarray,
    times: np.ndarray,
    m_max: int,
) -> np.ndarray:
    """
    Lanczos/Krylov evolution of <O(t)> using dense H_mat and O_mat,
    starting from psi0 in the same (reduced) space.

    This is the core engine; it does not know about symmetries or Paulis.
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


# ======================================================================
# 1. Dense sector basis (current production representation)
# ======================================================================


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


def _block_unitary_eig(
    S_gate: Gate,
    block: Sequence[int],
    tol: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Diagonalise the restriction of S_gate to a given block of qudits.

    Returns
    -------
    evals_block : (d_block,) complex
    evecs_block : (d_block, d_block) complex, columns are eigenvectors.
    """
    C_S = gate_to_circuit(S_gate)
    C_block: Circuit = C_S.local_circuit(list(block))
    U_block = C_block.unitary().toarray()
    evals, evecs = np.linalg.eig(U_block)

    # Normalise columns
    for j in range(evecs.shape[1]):
        v = evecs[:, j]
        nrm = np.linalg.norm(v)
        if nrm > tol:
            evecs[:, j] = v / nrm

    return evals, evecs


def build_dense_sector_basis(
    S_gate: Gate,
    blocks: List[List[int]],
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
    D = int(np.prod(dims))

    # Local block eigenpairs
    local_eigs: List[Tuple[np.ndarray, np.ndarray]] = []
    for block in blocks:
        evals_b, evecs_b = _block_unitary_eig(S_gate, block, tol=tol)
        local_eigs.append((evals_b, evecs_b))

    # Build global eigenvectors as tensor products of block eigenvectors
    basis_vecs: List[np.ndarray] = []

    index_grids = np.array(
        np.meshgrid(
            *[np.arange(len(evals_b)) for (evals_b, _) in local_eigs],
            indexing="ij",
        ),
        dtype=int,
    ).reshape(len(local_eigs), -1).T

    for choice in index_grids:
        lam_global = 1.0 + 0.0j
        vec_global = np.array([1.0 + 0.0j])
        total_dim = 1

        for (evals_b, vecs_b), block_idx in zip(local_eigs, range(len(local_eigs))):
            j = choice[block_idx]
            lam_b = evals_b[j]
            v_b = vecs_b[:, j]  # (d_block,)

            lam_global *= lam_b
            vec_global = np.kron(vec_global, v_b)
            total_dim *= v_b.shape[0]

        if np.abs(lam_global - target_eigval) < tol:
            nrm = np.linalg.norm(vec_global)
            if nrm > 0:
                vec_global = vec_global / nrm
            assert total_dim == D, "Global vector dimension mismatch."
            basis_vecs.append(vec_global)

    if not basis_vecs:
        raise ValueError(f"No eigenvectors found for eigenvalue {target_eigval}.")

    B = np.column_stack(basis_vecs)  # (D, D_sec)
    # Orthonormalise to get a clean basis
    Q, _ = np.linalg.qr(B)
    return DenseSectorBasis(Q)


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


# ======================================================================
# 2. Tensor-sector basis (MPS/MPO-ready abstraction)
# ======================================================================


@dataclass
class TensorSectorBasis:
    """
    Tensor-product symmetry sector basis for a given eigenvalue of S.

    The idea:
      - For each block b, diagonalise U_S on that block -> (λ_{b,r}, |φ_{b,r}>).
      - Global eigenvectors of S are tensor products of chosen local eigenvectors.
      - We keep only those global products with total eigenvalue = target_eigval.

    This class currently constructs full Hilbert vectors on demand via
    `basis_vector_full`, but its interface is compatible with future
    MPS/MPO implementations that avoid explicit length-D vectors.
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
            evals_b, evecs_b = _block_unitary_eig(S_gate, block, tol=tol)
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

        idx_tuple = self.sector_multi_indices[k]

        psi_block_list: List[np.ndarray] = []
        for b, r_b in enumerate(idx_tuple):
            v_b = self.local_evecs[b][:, r_b]   # (d_block_b,)
            psi_block_list.append(v_b)

        psi_full = psi_block_list[0]
        for v_b in psi_block_list[1:]:
            psi_full = np.kron(psi_full, v_b)

        if psi_full.shape[0] != self.dim_total:
            raise ValueError(
                f"Constructed sector basis vector has dim {psi_full.shape[0]}, "
                f"expected {self.dim_total}."
            )

        return psi_full

    def project_state_to_sector(self, psi: np.ndarray) -> np.ndarray:
        """
        Project a full Hilbert-space state |psi> into the sector and return
        coordinates in the sector basis (length = dim_sector):

            c_k = <φ_k | psi>, where |φ_k> is the k-th sector basis vector.
        """
        psi = np.asarray(psi, dtype=np.complex128).reshape(self.dim_total)

        coeffs = np.empty(self.dim_sector, dtype=np.complex128)
        for k in range(self.dim_sector):
            phi_k = self.basis_vector_full(k)
            coeffs[k] = np.vdot(phi_k, psi)

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


# ======================================================================
# 3. High-level symmetry-reduced Krylov driver
# ======================================================================


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

    This function is the "production" entry point. The DenseSectorBasis it uses
    could later be replaced by a TensorSectorBasis + MPS/MPO backend, without
    changing its signature.
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


# ======================================================================
# 4. Thin compatibility wrappers for older tests / APIs
# ======================================================================


def _local_symmetry_block_eigendecomposition(
    S_gate: Gate,
    block: List[int],
    dims: Sequence[int] | np.ndarray,
    tol: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Legacy helper used in some tests.

    Now just a thin wrapper around `_block_unitary_eig`. The `dims` argument
    is unused but kept for backwards compatibility.
    """
    _ = dims  # unused; kept for signature compatibility
    return _block_unitary_eig(S_gate, block, tol=tol)


def build_symmetry_sector_basis(
    S_gate: Gate,
    blocks: List[List[int]],
    dims: Sequence[int] | np.ndarray,
    target_eigval: complex = +1.0,
    tol: float = 1e-8,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Legacy API returning (B_sector, sector_eigvals).

    Internally uses `build_dense_sector_basis`. The eigenvalues returned are
    computed via the full U_S for diagnostic purposes and should all be
    ~ target_eigval within tolerance.
    """
    dims = np.asarray(dims, dtype=int).reshape(-1)
    D = int(np.prod(dims))

    dense_basis = build_dense_sector_basis(
        S_gate=S_gate,
        blocks=blocks,
        dims=dims,
        target_eigval=target_eigval,
        tol=tol,
    )
    B_sector = dense_basis.basis  # (D, D_sec)

    # Compute eigenvalues of S for each basis vector (mainly for tests)
    C_S = gate_to_circuit(S_gate)
    U_S = C_S.unitary().toarray()
    sector_eigvals = np.empty(B_sector.shape[1], dtype=complex)
    for j in range(B_sector.shape[1]):
        v = B_sector[:, j]
        sector_eigvals[j] = np.vdot(v, U_S @ v)

    assert B_sector.shape[0] == D
    return B_sector, sector_eigvals


def reduce_operator_to_sector(
    A_h: np.ndarray,
    B_sector: np.ndarray,
) -> np.ndarray:
    """
    Legacy helper: reduce a full Hilbert-space operator A_h to the symmetry sector
    spanned by columns of B_sector, via

        A_sector = B_sector^\dagger A_h B_sector.
    """
    return B_sector.conj().T @ A_h @ B_sector


def project_state_to_sector(
    psi: np.ndarray,
    B_sector: np.ndarray,
) -> np.ndarray:
    """
    Legacy helper: represent a global state |psi> in the symmetry sector basis:

        |psi_sec> = B_sector^\dagger |psi>
    """
    return B_sector.conj().T @ psi
