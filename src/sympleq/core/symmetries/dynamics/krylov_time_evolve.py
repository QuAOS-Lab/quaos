import numpy as np
from typing import Sequence, Optional
from sympleq.core.symmetries.dynamics.symmetry_reduced_krylov import (
    krylov_observable_dynamics_symmetry_reduced,
)
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits.utils import pauli_unitary_qudit
from sympleq.core.circuits import Circuit, gate_to_circuit
from sympleq.core.states.state import State
from sympleq.core.circuits.gates import Gate
from sympleq.core.paulis.utils import apply_paulisum_to_state_dense


def _apply_single_qudit_unitary(
    psi: np.ndarray,
    dims: np.ndarray,
    q: int,
    U_local,
) -> np.ndarray:
    """
    Apply a single-qudit unitary U_local on qudit q to the flat statevector psi.

    psi shape: (D,), D = prod(dims)
    dims: (n_qudits,)
    q: index of qudit to transform
    U_local: (d_q, d_q) array or sparse matrix
    """
    psi = np.asarray(psi, dtype=np.complex128)
    dims = np.asarray(dims, dtype=int).reshape(-1)
    n_qudits = dims.size
    if not (0 <= q < n_qudits):
        raise ValueError(f"Qudit index {q} out of range for dims {dims}.")

    D = int(np.prod(dims))
    if psi.shape != (D,):
        raise ValueError(f"Incompatible psi shape {psi.shape} for dims {dims}, D={D}.")

    # Reshape to tensor and move axis q to front: shape (d_q, rest)
    psi_tensor = psi.reshape(dims)
    psi_perm = np.moveaxis(psi_tensor, q, 0)
    d_q = dims[q]
    rest_dim = int(D // d_q)
    psi_perm_flat = psi_perm.reshape(d_q, rest_dim)

    # Apply U_local on the leading index
    psi_perm_flat_out = U_local @ psi_perm_flat

    # Reshape back and invert axis move
    psi_perm_out = psi_perm_flat_out.reshape((d_q,) + psi_perm.shape[1:])
    psi_out_tensor = np.moveaxis(psi_perm_out, 0, q)
    return psi_out_tensor.reshape(D)


def apply_pauli_row_to_state(
    x_row: np.ndarray,
    z_row: np.ndarray,
    dims: np.ndarray,
    psi: np.ndarray,
) -> np.ndarray:
    """
    Apply a single Pauli term (specified by exponents x_row, z_row over all qudits)
    to the flat statevector psi, without building the full operator.

    Parameters
    ----------
    x_row, z_row : (n_qudits,) int arrays
        Exponents of X and Z on each qudit.
    dims : (n_qudits,) int array
        Local dimensions.
    psi : (D,) complex ndarray
        Statevector in computational basis with D = prod(dims).

    Returns
    -------
    P_psi : (D,) complex ndarray
        Result of applying the Pauli operator to psi.
    """
    psi_out = np.asarray(psi, dtype=np.complex128)
    dims = np.asarray(dims, dtype=int).reshape(-1)
    n_qudits = dims.size

    for q in range(n_qudits):
        x_q = int(x_row[q])
        z_q = int(z_row[q])
        if x_q == 0 and z_q == 0:
            continue  # identity on this site

        d_q = int(dims[q])
        U_q = pauli_unitary_qudit(d_q, x_q, z_q)  # (d_q, d_q) operator on that qudit
        psi_out = _apply_single_qudit_unitary(psi_out, dims, q, U_q)

    return psi_out


def apply_paulisum_to_state(H: PauliSum, psi: np.ndarray) -> np.ndarray:
    """
    Compute (H psi) where H is a PauliSum and psi is a flat statevector.

      H = ∑_k w_k P_k  ⇒  H|ψ⟩ = ∑_k w_k P_k |ψ⟩

    Uses tableau rows of H, no full Hamiltonian matrix.
    """
    psi = np.asarray(psi, dtype=np.complex128)
    dims = np.asarray(H.dimensions, dtype=int).reshape(-1)
    D = int(np.prod(dims))
    if psi.shape != (D,):
        raise ValueError(f"Incompatible psi shape {psi.shape} for dims {dims}, D={D}.")

    tableau = H.tableau  # shape (n_terms, 2*n_qudits)
    weights = np.asarray(H.weights, dtype=np.complex128)
    n_terms, twon = tableau.shape
    n_qudits = twon // 2

    if weights.shape[0] != n_terms:
        raise ValueError(
            f"PauliSum.weights has length {weights.shape[0]}, "
            f"but tableau has {n_terms} rows."
        )
    if n_qudits != dims.size:
        raise ValueError("Mismatch between tableau qudit count and dims.")

    x_all = tableau[:, :n_qudits]
    z_all = tableau[:, n_qudits:]

    out = np.zeros_like(psi)
    for k in range(n_terms):
        x_row = x_all[k]
        z_row = z_all[k]
        w_k = weights[k]

        P_psi = apply_pauli_row_to_state(x_row, z_row, dims, psi)
        out += w_k * P_psi

    return out

def project_to_symmetry_sector(
    psi: np.ndarray,
    dims: np.ndarray,
    C_S: Circuit,
    eigval: complex = +1.0,
) -> np.ndarray:
    """
    Project psi onto the eigen-sector of S with eigenvalue `eigval`,
    assuming S^2 = I and eigenvalues in {+1, -1} or phases on the unit circle.

    Uses C_S.act_on_state rather than building U_S explicitly.

    P_eigval ~ (I + eigval* S)/2  (for Z2-like symmetry).
    """
    psi = np.asarray(psi, dtype=np.complex128)
    dims = np.asarray(dims, dtype=int).reshape(-1)
    D = int(np.prod(dims))
    if psi.shape != (D,):
        raise ValueError(f"Incompatible psi shape {psi.shape} for dims {dims}, D={D}.")

    state = State(psi, dims)
    state_S = C_S.act_on_state(state)
    psi_S = state_S.as_array()

    # Simple projector P = (I + eigval * S)/2
    psi_proj = 0.5 * (psi + eigval * psi_S)
    norm = np.linalg.norm(psi_proj)
    if norm < 1e-14:
        raise ValueError(
            "Projection onto symmetry sector nearly zero; "
            "initial state has negligible weight in this sector."
        )
    return psi_proj / norm


def krylov_observable_dynamics_symmetry(
    H: PauliSum,
    O: PauliSum,
    T_gate: Gate,
    S_gate: Gate,
    psi0: np.ndarray,
    times,
    m_max: int = 64,
    symmetry_eigval: complex | None = None,
) -> np.ndarray:
    """
    Symmetry-aware Krylov evolution <O(t)> without explicit Hilbert matrices.
    """
    times = np.asarray(times, dtype=float)
    dims = np.asarray(H.dimensions, dtype=int).reshape(-1)
    D = int(np.prod(dims))
    psi0 = np.asarray(psi0, dtype=np.complex128).reshape(D)

    # 1) Conjugate H and O by T^{-1} at PauliSum level
    T_inv_gate = T_gate.inv()
    H_prime = T_inv_gate.act(H)
    O_prime = T_inv_gate.act(O)

    # 2) Map |psi0> -> |psi0'> = U_T^{-1} |psi0> using *circuit* decomposition of T^{-1}
    state0 = State(psi0, dims)
    state0_prime = _apply_clifford_gate_via_circuit(T_inv_gate, state0)
    psi0_prime = state0_prime.as_array()

    # 3) Optional projection into a symmetry sector of S in the H' basis
    if symmetry_eigval is not None:
        psi0_prime = project_to_symmetry_sector_statevector(
            psi0_prime, dims, S_gate, eigval=symmetry_eigval
        )

    # 4) Krylov dynamics in H' basis using PauliSum matvecs only
    expvals_prime = krylov_observable_dynamics_plain(
        H=H_prime,
        O=O_prime,
        psi0=psi0_prime,
        times=times,
        m_max=m_max,
    )

    # 5) No need to transform back; expectations are invariant under conjugation
    return expvals_prime


def krylov_observable_dynamics_plain(
    H: PauliSum,
    O: PauliSum,
    psi0: np.ndarray,
    times: Sequence[float] | np.ndarray,
    m_max: int = 64,
) -> np.ndarray:
    """
    Krylov/Lanczos time evolution of ⟨O(t)⟩ under a Hamiltonian H given
    as a PauliSum, without using any symmetry information.

    The algorithm:
      1. Run Lanczos on H from |psi0⟩ using apply_paulisum_to_state_dense(H, ·)
         to build an orthonormal Krylov basis V and tridiagonal T_k.
      2. Diagonalise T_k = U diag(λ) U†.
      3. Represent O in the Krylov basis: O_K = V† O V.
      4. Propagate the initial Krylov vector e₁ in the eigenbasis of T_k:
         e₁ → e^{-i λ t} e₁, and compute
             ⟨O(t)⟩ = e₁† O_K(t) e₁
         where O_K(t) picks up phases e^{i(λ_i - λ_j)t}.

    Parameters
    ----------
    H : PauliSum
        Hamiltonian.
    O : PauliSum
        Observable.
    psi0 : (D,) complex ndarray
        Initial statevector.
    times : array-like
        Times at which to compute ⟨O(t)⟩.
    m_max : int
        Maximum Krylov dimension.

    Returns
    -------
    expectations : (len(times),) complex ndarray
        Approximate ⟨O(t)⟩.
    """
    times = np.asarray(times, dtype=float)
    dims = np.asarray(H.dimensions, dtype=int).reshape(-1)
    D = int(np.prod(dims))
    psi0 = np.asarray(psi0, dtype=np.complex128).reshape(D)

    # 0) Normalise initial state
    norm0 = np.linalg.norm(psi0)
    if norm0 < 1e-14:
        raise ValueError("Initial state has (near) zero norm.")
    v0 = psi0 / norm0

    # 1) Lanczos: build V, alpha, beta
    m_max = int(m_max)
    if m_max <= 0:
        raise ValueError("m_max must be positive")

    V = np.zeros((D, m_max), dtype=np.complex128)
    alpha = np.zeros(m_max, dtype=float)
    beta = np.zeros(m_max - 1, dtype=float)

    V[:, 0] = v0
    w = apply_paulisum_to_state_dense(H, v0)
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

        w = apply_paulisum_to_state_dense(H, vj)
        alpha[j] = np.vdot(vj, w).real
        w = w - alpha[j] * vj - beta[j - 1] * V[:, j - 1]
        k = j + 1

    alpha = alpha[:k]
    beta = beta[: max(0, k - 1)]
    V = V[:, :k]

    # 2) Tridiagonal T_k
    T_k = np.diag(alpha)
    if k > 1:
        T_k += np.diag(beta, 1) + np.diag(beta, -1)

    # 3) Diagonalise T_k
    evals, U = np.linalg.eigh(T_k)
    evals = evals.astype(float)
    U_dag = U.conj().T

    # 4) Observable in Krylov basis
    def apply_O(vec: np.ndarray) -> np.ndarray:
        return apply_paulisum_to_state_dense(O, vec)

    O_K = np.zeros((k, k), dtype=np.complex128)
    for j in range(k):
        Oj = apply_O(V[:, j])
        for i in range(k):
            O_K[i, j] = np.vdot(V[:, i], Oj)

    # Transform to eigenbasis of T_k
    O_K_eig = U_dag @ O_K @ U

    # Initial Krylov vector is e1
    e1 = np.zeros(k, dtype=np.complex128)
    e1[0] = 1.0
    e1_eig = U_dag @ e1

    # 5) Time evolution inside Krylov space
    expectations = np.empty(times.shape[0], dtype=np.complex128)
    for idx, t in enumerate(times):
        phase = np.exp(-1j * evals * t)  # exp(-i E t)
        # O(t) in eigenbasis: O_ij e^{i(λ_i - λ_j)t}
        O_t_eig = O_K_eig * (phase[None, :] * phase.conj()[:, None])
        expectations[idx] = np.vdot(e1_eig, O_t_eig @ e1_eig)

    return expectations


def project_to_symmetry_sector_statevector(
    psi: np.ndarray,
    dims: Sequence[int] | np.ndarray,
    S_gate: Gate,
    eigval: complex = +1.0,
) -> np.ndarray:
    """
    Project a statevector psi onto the eigenspace of the symmetry S with
    eigenvalue `eigval` (typically ±1), without ever forming U_S as a matrix.

    We:
      1. Decompose the generic S_gate into a Circuit of known gates
         using gate_to_circuit(S_gate).
      2. Apply that circuit gate-by-gate to the State, using each
         primitive gate's act_on_state implementation.
      3. Form the projected state
             |psi_proj> ∝ |psi> + eigval * S|psi>.
    """
    dims = np.asarray(dims, dtype=int).reshape(-1)
    D = int(np.prod(dims))
    psi = np.asarray(psi, dtype=np.complex128).reshape(D)

    # Build State from psi
    state = State(psi, dims)

    # Decompose S into a circuit of primitive gates and apply it
    C_S = gate_to_circuit(S_gate)
    state_S = state
    for g in C_S.gates:
        # Each primitive gate (Hadamard, PHASE, SUM, SWAP, etc.)
        # should have an efficient act_on_state implementation.
        state_S = g.act_on_state(state_S)

    psi_S = state_S.as_array()

    # Project onto eigval sector: (I + eigval * U_S) |psi>
    psi_proj = psi + eigval * psi_S

    norm = np.linalg.norm(psi_proj)
    if norm < 1e-12:
        raise ValueError(
            "Projection annihilated the state; the given state has negligible "
            f"overlap with the symmetry eigenspace eigval={eigval}."
        )

    psi_proj /= norm
    return psi_proj


def _apply_clifford_gate_via_circuit(gate: Gate, state: State) -> State:
    """
    Apply a (possibly generic/composite) Clifford Gate to a State by first
    decomposing it into a Circuit of primitive gates via gate_to_circuit
    and then applying those primitives gate-by-gate.

    This avoids calling Gate.unitary() for generic symplectic Gates, and
    relies only on the act_on_state implementations of the known gate
    subclasses (Hadamard, PHASE, SUM, SWAP, CNOT, PauliGate, ...).
    """
    circ = gate_to_circuit(gate)
    out = state
    for g in circ.gates:
        out = g.act_on_state(out)
    return out

def observable_dynamics_krylov(
    H: PauliSum,
    O: PauliSum,
    psi0: np.ndarray,
    times: Sequence[float] | np.ndarray,
    *,
    T_gate: Optional[Gate] = None,
    S_gate: Optional[Gate] = None,
    symmetry_eigval: complex | None = None,
    use_reduced_basis: bool = False,
    m_max: int = 64,
) -> np.ndarray:
    """
    High-level, production-ready entry point for Krylov-based dynamics
    of ⟨O(t)⟩ under a Hamiltonian H.

    This function chooses between:
      - plain Krylov (no symmetry),
      - symmetry-aware Krylov in the full Hilbert space,
      - symmetry-reduced Krylov in a symmetry sector,

    depending on the presence of (T_gate, S_gate) and use_reduced_basis.

    Parameters
    ----------
    H : PauliSum
        Hamiltonian.
    O : PauliSum
        Observable.
    psi0 : np.ndarray
        Initial statevector in the *original* basis of H.
        Shape (D,), where D = ∏ dims.
    times : array-like
        Times at which to compute ⟨O(t)⟩.
    T_gate : Gate | None, optional
        Clifford gate implementing the similarity transform T, such that
        H' = T^{-1} H T commutes with S_gate. If None, no symmetry is used.
    S_gate : Gate | None, optional
        Clifford symmetry gate S that commutes with H' and admits a block
        structure (for symmetry reduction). If None, no symmetry is used.
    symmetry_eigval : complex | None, optional
        Desired eigenvalue of S used to select a symmetry sector. Only
        relevant if S_gate is not None. If None, the full spectrum of S
        is implicitly included (no sector projection).
    use_reduced_basis : bool, optional
        If True and (T_gate, S_gate) are provided, use the symmetry-reduced
        Krylov algorithm that works in a symmetry sector's reduced basis
        (no explicit H' or O' Hilbert matrices).
    m_max : int, optional
        Maximum Krylov dimension.

    Returns
    -------
    expectations : np.ndarray
        Array of shape (len(times),) containing ⟨O(t)⟩.
    """
    times = np.asarray(times, dtype=float)

    # Case 1: no symmetry information – plain Krylov
    if T_gate is None or S_gate is None:
        return krylov_observable_dynamics_plain(
            H=H,
            O=O,
            psi0=psi0,
            times=times,
            m_max=m_max,
        )

    # Case 2: symmetry-aware but *no* reduced basis: full Hilbert sector
    if not use_reduced_basis:
        return krylov_observable_dynamics_symmetry(
            H=H,
            O=O,
            T_gate=T_gate,
            S_gate=S_gate,
            psi0=psi0,
            times=times,
            m_max=m_max,
            symmetry_eigval=symmetry_eigval,
        )

    # Case 3: symmetry-reduced Krylov in a symmetry sector
    if symmetry_eigval is None:
        raise ValueError(
            "use_reduced_basis=True requires an explicit symmetry_eigval, "
            "e.g. +1.0 or a phase on the unit circle."
        )

    return krylov_observable_dynamics_symmetry_reduced(
        H=H,
        O=O,
        T_gate=T_gate,
        S_gate=S_gate,
        psi0=psi0,
        times=times,
        symmetry_eigval=symmetry_eigval,
        m_max=m_max,
    )
