from __future__ import annotations
from abc import ABC
import numpy as np
from typing import TypeVar, Self
import scipy.sparse as sp

from sympleq.core.paulis import PauliObject
from sympleq.core.circuits.utils import embed_symplectic, transvection_matrix
from sympleq.core.circuits.random_symplectic import symplectic_random_transvection
from sympleq.core.circuits.find_symplectic import map_pauli_sum_to_target_tableau
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.core.circuits.target import get_phase_vector


# We define a type using TypeVar to let the type checker know that
# the input and output of the `act` function share the same type.
P = TypeVar("P", bound="PauliObject")


class Gate(ABC):
    """
    Abstract base class for dimension-independent Clifford gates.

    Gates are defined purely by their symplectic matrix and phase vector(s).
    They do not store qudit indices or dimensions - these are passed at act-time.

    The symplectic matrix defines how Pauli operators transform under conjugation.
    The phase vector handles the phase acquired during the transformation.
    Some gates (like PHASE) have dimension-dependent phase vectors, which are
    stored in `_exceptional_phase_vectors`.
    """

    def __init__(self, name: str, symplectic: np.ndarray,
                 phase_vector: np.ndarray | None = None,
                 exceptional_phase_vectors: dict[int, np.ndarray] | None = None):
        self._name = name
        self._n_qudits = symplectic.shape[0] // 2
        self._symplectic = symplectic.astype(int)
        self._symplectic.setflags(write=False)

        if phase_vector is None:
            phase_vector = np.zeros(2 * self._n_qudits, dtype=int)
        self._phase_vector = np.asarray(phase_vector, dtype=int)
        self._phase_vector.setflags(write=False)

        # Store dimension-specific phase vectors (e.g., qubits are special for PHASE gate)
        self._exceptional_phase_vectors: dict[int, np.ndarray] = {}
        if exceptional_phase_vectors is not None:
            for dim, pv in exceptional_phase_vectors.items():
                pv_arr = np.asarray(pv, dtype=int)
                pv_arr.setflags(write=False)
                self._exceptional_phase_vectors[dim] = pv_arr

        # Precompute matrices for phase calculation (see Eq.[7] in PHYSICAL REVIEW A 71, 042315 (2005))
        # U = [[0_n, 0_n],
        #      [I_n, 0_n]]
        U = np.zeros((2 * self._n_qudits, 2 * self._n_qudits), dtype=int)
        U[self._n_qudits:, :self._n_qudits] = np.eye(self._n_qudits, dtype=int)
        self._U_symplectic_conjugated = self._symplectic.T @ U @ self._symplectic

        self._V_diag = np.diag(self._U_symplectic_conjugated)
        # This is the part associated with the quadratic form.
        # Remove diagonal part to match definition in Eq.[7].
        self._p_part = 2 * np.triu(self._U_symplectic_conjugated) - np.diag(self._V_diag)

        # Will be set by _Gates to link to inverse gate
        self._inverse: Self | None = None

    @classmethod
    def from_random(cls, n_qudits: int, dimension: int, num_transvections: int | None = None) -> "Gate":
        """
        Generate a random Clifford gate by composing random transvections.

        Parameters
        ----------
        n_qudits : int
            Number of qudits the gate acts on.
        dimension : int
            Local Hilbert space dimension (e.g., 2 for qubits).
        num_transvections : int | None
            Number of transvections to compose. If None, defaults to 4*n_qudits.

        Returns
        -------
        Gate
            A random Clifford gate with the generated symplectic matrix.
        """

        symplectic = symplectic_random_transvection(n_qudits, dimension, num_transvections)
        # For random gates, we use zero phase vector (phases depend on specific gate sequence)
        phase_vector = get_phase_vector(symplectic, dimension)

        return _GenericGate("random", symplectic, phase_vector)

    @classmethod
    def solve_from_target(cls, input_tableau: np.ndarray, target_tableau: np.ndarray) -> "Gate":
        """
        Find a Clifford gate that maps the input Pauli tableau to the target tableau.

        Uses symplectic transvections to find a symplectic matrix F such that
        input_tableau @ F = target_tableau (mod 2).

        Parameters
        ----------
        input_tableau : np.ndarray
            Input Pauli tableau of shape (m, 2n) where m is the number of Paulis
            and n is the number of qudits.
        target_tableau : np.ndarray
            Target Pauli tableau of the same shape.

        Returns
        -------
        Gate
            A Clifford gate whose symplectic matrix performs the mapping.

        Raises
        ------
        ValueError
            If the tableaus have different shapes or are not mappable via Clifford.

        Notes
        -----
        Currently only works for GF(2) (qubits). The input and target must have
        matching symplectic product matrices for a Clifford mapping to exist.
        """

        input_tableau = np.asarray(input_tableau, dtype=int)
        target_tableau = np.asarray(target_tableau, dtype=int)

        if input_tableau.shape != target_tableau.shape:
            raise ValueError(
                f"Tableau shapes must match: {input_tableau.shape} vs {target_tableau.shape}"
            )

        if input_tableau.ndim == 1:
            input_tableau = input_tableau.reshape(1, -1)
            target_tableau = target_tableau.reshape(1, -1)

        n_qudits = input_tableau.shape[1] // 2

        symplectic = map_pauli_sum_to_target_tableau(input_tableau, target_tableau)
        phase_vector = np.zeros(2 * n_qudits, dtype=int)

        return _GenericGate("target", symplectic, phase_vector)

    @property
    def name(self) -> str:
        return self._name

    @property
    def n_qudits(self) -> int:
        return self._n_qudits

    @property
    def symplectic(self) -> np.ndarray:
        return self._symplectic

    def phase_vector(self, dimension: int = 0) -> np.ndarray:
        """
        Get the phase vector for a given dimension.

        Some gates have dimension-specific phase vectors (e.g., PHASE gate for qubits).
        If no exceptional phase vector exists for the given dimension, returns the default.
        """
        if dimension in self._exceptional_phase_vectors:
            return self._exceptional_phase_vectors[dimension]
        return self._phase_vector

    def __repr__(self) -> str:
        return f"Gate(name={self._name}, n_qudits={self._n_qudits})"

    def act(self, pauli: P, qudits: int | tuple[int, ...]) -> P:
        """
        Apply this gate to a Pauli object at the specified qudit indices.

        Returns the updated Pauli with transformed tableau and phases.
        See Eq.[7] in PHYSICAL REVIEW A 71, 042315 (2005).

        Parameters
        ----------
        pauli : Pauli | PauliString | PauliSum
            The Pauli object to transform.
        qudits : int | tuple[int, ...]
            The qudit index (for single-qudit gates) or tuple of indices
            (for multi-qudit gates) on which the gate acts.

        Returns
        -------
        The transformed Pauli object of the same type as the input.
        """
        if isinstance(qudits, int):
            qudits = (qudits,)

        affected_qudits = np.asarray(qudits, dtype=int)

        if len(affected_qudits) != self._n_qudits:
            raise ValueError(f"Gate acts on {self._n_qudits} qudits, but {len(affected_qudits)} indices provided.")

        T = pauli.tableau

        # Precompute tableau mask to select affected columns
        tableau_mask = np.concatenate([affected_qudits, affected_qudits + pauli.n_qudits()])
        T_affected = T[:, tableau_mask]

        pauli_dimensions = pauli.dimensions[affected_qudits]
        relevant_dimensions = np.tile(pauli_dimensions, 2)
        relevant_lcm = int(np.lcm.reduce(pauli_dimensions))

        # Apply symplectic transformation with modulo reduction
        updated_tableau = np.mod(T_affected @ self._symplectic.T, relevant_dimensions)
        new_tableau = T.copy()
        new_tableau[:, tableau_mask] = updated_tableau

        # Compute phase contribution
        phase_vec = self.phase_vector(relevant_lcm) % (2 * pauli.lcm)
        modified_phase_vector = phase_vec - self._V_diag
        linear_terms = T_affected @ modified_phase_vector
        quadratic_terms = np.sum(T_affected * (T_affected @ self._p_part), axis=1)

        dimensional_factor = pauli.lcm // relevant_lcm
        acquired_phases = (linear_terms + quadratic_terms) * dimensional_factor

        new_phases = (pauli.phases + acquired_phases) % (2 * pauli.lcm)

        return pauli.__class__(tableau=new_tableau, dimensions=pauli.dimensions,
                               weights=pauli.weights, phases=new_phases)

    def local_unitary(self, dimension: int | None = None) -> sp.csr_matrix:
        """
        Compute the unitary matrix representation of this gate.

        Parameters
        ----------
        dimension : int
            The local Hilbert space dimension (e.g., 2 for qubits, 3 for qutrits).

        Returns
        -------
        np.ndarray
            The d^n x d^n unitary matrix, where n is the number of qudits the gate acts on.
        """
        # Default implementation for generic gates - subclasses should override
        raise NotImplementedError(
            f"unitary() not implemented for {self.__class__.__name__}. "
            "Override in subclass or use a specific gate class."
        )

    def inverse(self) -> Self:
        """Return the inverse of this gate.

        For singleton gates (from GATES), returns the pre-linked inverse.
        Otherwise computes the inverse symplectic matrix.
        """
        if self._inverse is not None:
            return self._inverse

        # Compute inverse symplectic: C^{-1} = -Omega^T @ C^T @ Omega
        n = self._n_qudits
        zero_block = np.zeros((n, n), dtype=int)
        identity_block = np.eye(n, dtype=int)
        Omega = np.block([[zero_block, identity_block], [-identity_block, zero_block]])

        C_inv = -Omega.T @ self._symplectic.T @ Omega

        U = np.zeros((2 * n, 2 * n), dtype=int)
        U[n:, :n] = np.eye(n, dtype=int)

        U_C = (self._symplectic.T @ U @ self._symplectic)
        U_C_inv = (C_inv.T @ U @ C_inv)

        P_C = (2 * np.triu(U_C) - np.diag(np.diag(U_C)))
        P_C_inv = (2 * np.triu(U_C_inv) - np.diag(np.diag(U_C_inv)))

        # FIXME: handle exceptional phase vector correctly

        term1 = (-self.phase_vector() @ C_inv)
        term2 = (np.diag(U_C.T)) @ C_inv
        term3 = np.diag((C_inv.T @ P_C @ C_inv))
        term4 = np.diag(U_C_inv)
        term5 = np.diag(P_C_inv)

        phase_vector_inv = (term1 + term2 - term3 + term4 - term5)

        inv_name = self._name + "_inv" if not self._name.endswith("_inv") else self._name[:-4]

        return self.__class__(inv_name, C_inv, phase_vector_inv)

    def transvection(self, transvection_vector: np.ndarray | list, transvection_weight: int = 1) -> Gate:
        """
        Returns a new gate that is the transvection of this gate by the given vector.

        The transvection vector should be a 2n-dimensional vector where n is the number of qudits.
        Note: This returns a generic Gate, not a subclass instance.
        """
        if not isinstance(transvection_weight, int) and not isinstance(transvection_weight, np.int64):
            raise TypeError("Transvection weight must be an integer.")

        if isinstance(transvection_vector, list):
            transvection_vector = np.asarray(transvection_vector)

        T = transvection_matrix(transvection_vector, multiplier=transvection_weight)
        new_name = self._name if self._name.startswith("T-") else "T-" + self._name

        return _GenericGate(new_name, self._symplectic @ T, self._phase_vector.copy())

    def full_symplectic(self, qudits: tuple[int, ...] | int, n_qudits: int, p: int) -> np.ndarray:
        """
        Get the full 2n x 2n symplectic matrix for a gate acting on specific qudits.

        Parameters
        ----------
        qudits : tuple[int, ...] | int
            The qudit index(es) the gate acts on
        n_qudits : int
            Total number of qudits in the system
        p : int
            The prime dimension for modular arithmetic

        Returns
        -------
        np.ndarray
            The full 2n x 2n symplectic matrix mod p
        """
        if isinstance(qudits, int):
            qudits = (qudits,)
        F, _ = embed_symplectic(self.symplectic, self.phase_vector(p), qudits, n_qudits)
        return F % p


class _GenericGate(Gate):
    """
    A generic gate for computed gates (e.g., from transvection, composite gates).
    These are not singletons and don't have pre-computed inverses.
    Uses parent's inverse() which computes the inverse dynamically.
    """
    pass


class _HADAMARD(Gate):
    """Hadamard gate: X -> -Z, Z -> X (or inverse: X -> Z, Z -> -X)"""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse

        if is_inverse:
            symplectic = np.array([
                [0, 1],    # image of X:  X -> Z
                [-1, 0]    # image of Z:  Z -> -X
            ], dtype=int)
            name = "H_inv"
        else:
            symplectic = np.array([
                [0, -1],   # image of X:  X -> -Z
                [1, 0]     # image of Z:  Z -> X
            ], dtype=int)
            name = "H"

        super().__init__(name, symplectic)

    def local_unitary(self, dimension: int | None = None) -> sp.csr_matrix:
        from sympleq.core.circuits.utils import H_mat
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        U = H_mat(dimension)
        if self._is_inverse:
            return U.conj().T
        return U

    to_local_hilbert_space = local_unitary


class _PHASE(Gate):
    """Phase gate (S): X -> XZ, Z -> Z. Has special phase vector for qubits."""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse

        if is_inverse:
            symplectic = np.array([
                [1, -1],  # image of X:  X -> XZ^{-1}
                [0, 1]    # image of Z:  Z -> Z
            ], dtype=int).T
            name = "S_inv"
        else:
            symplectic = np.array([
                [1, 1],   # image of X:  X -> XZ
                [0, 1]    # image of Z:  Z -> Z
            ], dtype=int).T
            name = "S"

        # Qubits (dimension=2) have a special phase vector
        if is_inverse:
            exceptional = {2: np.array([-1, 0], dtype=int)}
        else:
            exceptional = {2: np.array([1, 0], dtype=int)}

        super().__init__(name, symplectic, exceptional_phase_vectors=exceptional)

    def local_unitary(self, dimension: int | None = None) -> sp.csr_matrix:
        from sympleq.core.circuits.utils import S_mat
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        U = S_mat(dimension)
        if self._is_inverse:
            return U.conj().T
        return U

    to_local_hilbert_space = local_unitary


class _CX(Gate):
    """CX gate: X0 -> X0 X1, X1 -> X1, Z0 -> Z0, Z1 -> Z0^{-1} Z1"""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse

        # CX is self-inverse for qubits, but not in general
        if is_inverse:
            symplectic = np.array([
                [1, -1, 0, 0],   # image of X0:  X0 -> X0 X1^{-1}
                [0, 1, 0, 0],    # image of X1:  X1 -> X1
                [0, 0, 1, 0],    # image of Z0:  Z0 -> Z0
                [0, 0, 1, 1]     # image of Z1:  Z1 -> Z0 Z1
            ], dtype=int).T
            name = "CX_inv"
        else:
            symplectic = np.array([
                [1, 1, 0, 0],    # image of X0:  X0 -> X0 X1
                [0, 1, 0, 0],    # image of X1:  X1 -> X1
                [0, 0, 1, 0],    # image of Z0:  Z0 -> Z0
                [0, 0, -1, 1]    # image of Z1:  Z1 -> Z0^{-1} Z1
            ], dtype=int).T
            name = "CX"

        super().__init__(name, symplectic)

    def local_unitary(self, dimension: int | None = None) -> sp.csr_matrix:
        """CX acts as |j,k⟩ -> |j, (j+k) mod d⟩ (or |j, (k-j) mod d⟩ for inverse)."""

        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        D = d * d  # total dimension
        U = np.zeros((D, D), dtype=complex)
        for j in range(d):
            for k in range(d):
                in_idx = j * d + k  # |j,k⟩
                if self._is_inverse:
                    out_k = (k - j) % d
                else:
                    out_k = (j + k) % d
                out_idx = j * d + out_k  # |j, out_k⟩
                U[out_idx, in_idx] = 1.0
        return sp.csr_matrix(U)

    to_local_hilbert_space = local_unitary


class _SWAP(Gate):
    """SWAP gate: X0 <-> X1, Z0 <-> Z1. Self-inverse."""

    def __init__(self):
        symplectic = np.array([
            [0, 1, 0, 0],  # image of X0:  X0 -> X1
            [1, 0, 0, 0],  # image of X1:  X1 -> X0
            [0, 0, 0, 1],  # image of Z0:  Z0 -> Z1
            [0, 0, 1, 0]   # image of Z1:  Z1 -> Z0
        ], dtype=int).T

        super().__init__("SWAP", symplectic)

    def local_unitary(self, dimension: int | None = None) -> sp.csr_matrix:
        """SWAP acts as |j,k⟩ -> |k,j⟩."""
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        D = d * d
        U = np.zeros((D, D), dtype=complex)
        for j in range(d):
            for k in range(d):
                in_idx = j * d + k   # |j,k⟩
                out_idx = k * d + j  # |k,j⟩
                U[out_idx, in_idx] = 1.0
        return sp.csr_matrix(U)

    to_local_hilbert_space = local_unitary

    def inverse(self) -> _SWAP:
        # SWAP is self-inverse
        return self


class _CZ(Gate):
    """Controlled-Z gate. Self-inverse."""

    def __init__(self):
        symplectic = np.array([
            [1, 0, 0, 1],  # image of X0:  X0 -> X0 Z1
            [0, 1, 1, 0],  # image of X1:  X1 -> Z0 X1
            [0, 0, 1, 0],  # image of Z0:  Z0 -> Z0
            [0, 0, 0, 1]   # image of Z1:  Z1 -> Z1
        ], dtype=int).T

        super().__init__("CZ", symplectic)

    def local_unitary(self, dimension: int | None = None) -> sp.csr_matrix:
        """CZ adds phase ω^{jk} to |j,k⟩ where ω = exp(2πi/d)."""
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        D = d * d
        omega = np.exp(2j * np.pi / d)
        U = np.zeros((D, D), dtype=complex)
        for j in range(d):
            for k in range(d):
                idx = j * d + k
                U[idx, idx] = omega ** (j * k)
        return sp.csr_matrix(U)

    to_local_hilbert_space = local_unitary

    def inverse(self) -> _CZ:
        # CZ is self-inverse
        return self


class _Gates:
    """
    Singleton container for pre-instantiated gates.

    Gates are accessed as properties, e.g., GATES.H, GATES.S, GATES.CX.
    Inverse gates are linked so that gate.inverse() returns the inverse singleton.
    """

    def __init__(self):
        # Single-qudit gates
        self._H = _HADAMARD(is_inverse=False)
        self._H_inv = _HADAMARD(is_inverse=True)
        self._H._inverse = self._H_inv
        self._H_inv._inverse = self._H

        self._S = _PHASE(is_inverse=False)
        self._S_inv = _PHASE(is_inverse=True)
        self._S._inverse = self._S_inv
        self._S_inv._inverse = self._S

        # Two-qudit gates
        self._CX = _CX(is_inverse=False)
        self._CX_inv = _CX(is_inverse=True)
        self._CX._inverse = self._CX_inv
        self._CX_inv._inverse = self._CX

        self._SWAP = _SWAP()
        # SWAP is self-inverse, already handled in the class

        self._CZ = _CZ()
        # CZ is self-inverse, already handled in the class

    # Hadamard
    @property
    def H(self) -> _HADAMARD:
        return self._H

    @property
    def H_inv(self) -> _HADAMARD:
        return self._H_inv

    # Phase (S)
    @property
    def S(self) -> _PHASE:
        return self._S

    @property
    def S_inv(self) -> _PHASE:
        return self._S_inv

    # CX (controlled-X / CNOT)
    @property
    def CX(self) -> _CX:
        return self._CX

    @property
    def CX_inv(self) -> _CX:
        return self._CX_inv

    # SWAP
    @property
    def SWAP(self) -> _SWAP:
        return self._SWAP

    # CZ
    @property
    def CZ(self) -> _CZ:
        return self._CZ


# Global singleton instance
GATES = _Gates()


class PauliGate(Gate):
    """
    A gate constructed from a PauliString.

    Unlike the singleton gates, PauliGate is dynamically created based on the input PauliString.
    The symplectic matrix is identity, and the phase vector encodes the Pauli conjugation effect.
    """

    def __init__(self, pauli):
        # Import here to avoid circular imports
        from sympleq.core.paulis import PauliString
        from sympleq.core.circuits.utils import symplectic_form

        if not isinstance(pauli, PauliString):
            raise TypeError("PauliGate requires a PauliString")

        self.pauli_string = pauli
        n = pauli.n_qudits()
        lcm = int(pauli.lcm)

        symplectic = np.eye(2 * n, dtype=int)
        phase_vector = (2 * symplectic_form(n, lcm) @ np.concatenate([pauli.x_exp, pauli.z_exp])) % (2 * lcm)

        # Store dimensions for this gate (needed for act method compatibility)
        self._dimensions = np.asarray(pauli.dimensions, dtype=int)

        super().__init__("Pauli", symplectic, phase_vector)

    @property
    def dimensions(self) -> np.ndarray:
        return self._dimensions

    def inverse(self) -> PauliGate:
        # Pauli gates are self-inverse (up to phase)
        return self

    def local_unitary(self, dimension: int | None = None) -> sp.csr_matrix:
        """
        Compute the unitary for this PauliGate.

        For PauliGate, dimension is optional since it's determined by the stored PauliString.
        """
        from sympleq.core.circuits.utils import pauli_unitary_from_tableau
        # Use the dimension from the PauliString
        d = int(self._dimensions[0])  # assumes uniform dimensions
        x = self.pauli_string.x_exp
        z = self.pauli_string.z_exp
        return pauli_unitary_from_tableau(d, x, z, convention="bare")

    to_local_hilbert_space = local_unitary

    def act(self, pauli: P, qudits: int | tuple[int, ...] | None = None) -> P:
        """
        Apply this PauliGate to a Pauli object.

        For PauliGate, qudits defaults to all qudits in order (0, 1, 2, ..., n-1)
        since the gate was constructed for a specific number of qudits.
        """
        if qudits is None:
            qudits = tuple(range(self._n_qudits))
        return super().act(pauli, qudits)


# Convenience aliases for backward compatibility
# These are the gate classes, not instances
HADAMARD = _HADAMARD
PHASE = _PHASE
CX = _CX
SWAP = _SWAP
CZ = _CZ
