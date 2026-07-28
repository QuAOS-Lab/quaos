from __future__ import annotations
from abc import ABC
import numpy as np
from typing import Self, overload

from sympleq._typing import IntArrayLike
from sympleq.core.paulis import PauliObject
from sympleq.core.paulis._typing import (
    TableauType, TableauLike, PhasesType, DimensionsType, HilbertOperator
)
from sympleq.core.circuits.utils import embed_symplectic, embed_unitary, transvection_matrix
from sympleq.core.circuits.random_symplectic import (
    symplectic_random_koenig_smolin_gf2,
    symplectic_random_transvection,
)
from sympleq.core.circuits.find_symplectic import map_paulisum_to_target_tableau
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.core.circuits.target import get_phase_vector
from sympleq.core.paulis.pauli_string import PauliString
from sympleq.core.paulis.pauli_sum import PauliSum


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

    def __init__(self, name: str, symplectic: TableauType,
                 phase_vector: PhasesType | None = None,
                 exceptional_phase_vectors: dict[int, PhasesType] | None = None):
        """
        Initialize a Clifford gate.

        Parameters
        ----------
        name : str
            Human-readable name for the gate.
        symplectic : TableauType
            The 2n x 2n symplectic matrix defining the Pauli transformation.
        phase_vector : PhasesType | None
            Phase vector for the gate. If None, defaults to zeros.
        exceptional_phase_vectors : dict[int, PhasesType] | None
            Dimension-specific phase vectors (e.g., PHASE gate for qubits).
        """
        self._name = name
        self._n_qudits = symplectic.shape[0] // 2
        self._symplectic = symplectic.astype(int)
        self._symplectic.setflags(write=False)

        if phase_vector is None:
            phase_vector = np.zeros(2 * self._n_qudits, dtype=int)
        self._phase_vector = np.asarray(phase_vector, dtype=int)
        self._phase_vector.setflags(write=False)

        # Store dimension-specific phase vectors (e.g., qubits are special for PHASE gate)
        self._exceptional_phase_vectors: dict[int, PhasesType] = {}
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
    def from_random(
        cls,
        n_qudits: int,
        dimension: int,
        num_transvections: int | None = None,
        *,
        sampler: str = "transvection",
        rng: np.random.Generator | None = None,
    ) -> Gate:
        """
        Generate a random Clifford gate.

        Parameters
        ----------
        n_qudits : int
            Number of qudits the gate acts on.
        dimension : int
            Local Hilbert space dimension (e.g., 2 for qubits).
        num_transvections : int | None
            Number of transvections to compose for ``sampler="transvection"``.
            If None, defaults to 4*n_qudits.
        sampler : str
            Random symplectic sampler to use. ``"transvection"`` preserves the
            historical behavior. ``"koenig-smolin"`` uses the Koenig-Smolin
            uniform index sampler and is available only for qubits
            (``dimension == 2``).
        rng : np.random.Generator | None
            Optional random generator for the selected sampler.

        Returns
        -------
        Gate
            A random Clifford gate with the generated symplectic matrix.
        """

        sampler_key = str(sampler).strip().lower().replace("_", "-")
        if sampler_key == "transvection":
            symplectic = symplectic_random_transvection(
                n_qudits,
                dimension,
                num_transvections,
                rng=rng,
            )
        elif sampler_key == "koenig-smolin":
            if dimension != 2:
                raise ValueError("sampler='koenig-smolin' is only implemented for dimension=2.")
            if num_transvections is not None:
                raise ValueError("num_transvections is not used with sampler='koenig-smolin'.")
            symplectic = symplectic_random_koenig_smolin_gf2(n_qudits, rng=rng)
        else:
            raise ValueError(
                "Unknown random Clifford sampler "
                f"{sampler!r}. Expected 'transvection' or 'koenig-smolin'."
            )

        # For random gates, we use zero phase vector (phases depend on specific gate sequence)
        phase_vector = get_phase_vector(symplectic, dimension)

        return _GenericGate("random", symplectic, phase_vector)

    # TODO: the following function should work for mixed qudits. the gate method should actually take two PauliSums
    #       (inclusive of phases) and should return a gate that maps the first to the second. This is the mthod that
    #       will use the functions in the new file that will contain a polished version of the functions in, e.g.,
    #       find_symplectic.py
    @classmethod
    def solve_from_target(cls, input_tableau: TableauLike, target_tableau: TableauLike,
                          dimension: int = DEFAULT_QUDIT_DIMENSION) -> Gate:
        """
        Find a Clifford gate that maps the input Pauli tableau to the target tableau.

        Uses symplectic transvections to find a symplectic matrix F such that
        input_tableau @ F = target_tableau (mod p), with p=`dimension`.

        Parameters
        ----------
        input_tableau : TableauLike
            Input Pauli tableau of shape (m, 2n) where m is the number of Paulis
            and n is the number of qudits.
        target_tableau : TableauLike
            Target Pauli tableau of the same shape.
        dimension : int
            Local Hilbert space dimension (e.g., 2 for qubits).

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
        Supports GF(p) for prime `dimension` via the compatibility layer in
        `find_symplectic.py`. The input and target must have matching symplectic
        product matrices for a Clifford mapping to exist.
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

        symplectic = map_paulisum_to_target_tableau(
            input_tableau,
            target_tableau,
            p=int(dimension),
            method="auto",
        )
        phase_vector = np.zeros(2 * n_qudits, dtype=int)

        return _GenericGate("target", symplectic, phase_vector)

    @property
    def name(self) -> str:
        """str : Human-readable name of the gate."""
        return self._name

    @property
    def n_qudits(self) -> int:
        """int : Number of qudits the gate acts on."""
        return self._n_qudits

    @property
    def symplectic(self) -> TableauType:
        """TableauType : The 2n x 2n symplectic matrix."""
        return self._symplectic

    def phase_vector(self, dimension: int | None = None) -> PhasesType:
        """
        Get the phase vector for a given dimension.

        Some gates have dimension-specific phase vectors (e.g., PHASE gate for qubits).
        If no exceptional phase vector exists for the given dimension, returns the default.

        Parameters
        ----------
        dimension : int
            The local Hilbert space dimension. Used to look up exceptional phase vectors.

        Returns
        -------
        PhasesType
            The phase vector for the given dimension.
        """
        if dimension in self._exceptional_phase_vectors:
            return self._exceptional_phase_vectors[dimension]
        return self._phase_vector

    def __repr__(self) -> str:
        return f"Gate(name={self._name}, n_qudits={self._n_qudits})"

    def __hash__(self):
        return hash((
            self.name,
            tuple(self.symplectic.flatten().tobytes()),
            tuple(self._phase_vector.tobytes()),
        ))

    def __eq__(self, other):
        if not isinstance(other, Gate):
            return False
        return np.all(self.symplectic == other.symplectic) and \
            np.all(self._phase_vector == other._phase_vector) and \
            np.all(self._exceptional_phase_vectors == other._exceptional_phase_vectors)

    @overload
    def act(self, pauli: PauliSum, qudits: int | tuple[int, ...]) -> PauliSum:
        ...

    @overload
    def act(self, pauli: PauliString, qudits: int | tuple[int, ...]) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliObject, qudits: int | tuple[int, ...]) -> PauliObject:
        ...

    def act(self, pauli, qudits):
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
        PauliObject
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

        return pauli.__class__(
            tableau=new_tableau, dimensions=pauli.dimensions,
            weights=pauli.weights, phases=new_phases
        )

    def act_in_hilbert_space(self, rho: HilbertOperator,
                             qudits: tuple[int, ...], dimensions: DimensionsType) -> HilbertOperator:
        """
        Apply this gate to a density matrix in Hilbert space.

        Computes rho_out = U rho U^dagger where U is the gate unitary embedded
        into the full Hilbert space.

        Parameters
        ----------
        rho : HilbertOperator
            The input density matrix.
        qudits : tuple[int, ...]
            Qudit indices on which the gate acts.
        dimensions : DimensionsType
            Local Hilbert space dimensions for each qudit.

        Returns
        -------
        HilbertOperator
            The transformed density matrix.
        """
        dimension = dimensions[qudits[0]]
        unitary = embed_unitary(self.local_unitary(dimension), qudits, dimensions)

        return unitary @ rho @ unitary.conjugate().transpose()

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        """
        Compute the unitary matrix representation of this gate.

        Parameters
        ----------
        dimension : int
            The local Hilbert space dimension (e.g., 2 for qubits, 3 for qutrits).

        Returns
        -------
        HilbertOperator
            The d^n x d^n unitary matrix, where n is the number of qudits the gate acts on.
        """
        # Default implementation for generic gates - subclasses should override
        raise NotImplementedError(
            f"unitary() not implemented for {self.__class__.__name__}. "
            "Override in subclass or use a specific gate class."
        )

    def to_local_hilbert_space(self, dimension: int | None = None) -> HilbertOperator:
        return self.local_unitary(dimension)

    def inverse(self) -> Self:
        """
        Return the inverse of this gate.

        For singleton gates (from GATES), returns the pre-linked inverse.
        Otherwise computes the inverse symplectic matrix.

        Returns
        -------
        Gate
            The inverse gate such that gate @ gate.inverse() = Identity.
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

    def transvection(self, transvection_vector: IntArrayLike, transvection_weight: int = 1) -> Gate:
        """
        Return a new gate that is the transvection of this gate by the given vector.

        Parameters
        ----------
        transvection_vector : IntArrayLike
            A 2n-dimensional vector where n is the number of qudits.
        transvection_weight : int
            Multiplier for the transvection. Default is 1.

        Returns
        -------
        Gate
            A new generic Gate with the transvected symplectic matrix.
        """
        if not isinstance(transvection_weight, int) and not isinstance(transvection_weight, np.int64):
            raise TypeError("Transvection weight must be an integer.")

        transvection_vector = np.asarray(transvection_vector, dtype=int)
        T = transvection_matrix(transvection_vector, multiplier=transvection_weight)
        new_name = self._name if self._name.startswith("T-") else "T-" + self._name

        return _GenericGate(new_name, self._symplectic @ T, self._phase_vector.copy())

    def full_symplectic(self, qudits: tuple[int, ...] | int, n_qudits: int,
                        dimension: int | None = None) -> TableauType:
        """
        Get the full 2n x 2n symplectic matrix for a gate acting on specific qudits.

        Parameters
        ----------
        qudits : tuple[int, ...] | int
            The qudit index(es) the gate acts on
        n_qudits : int
            Total number of qudits in the system
        dimension : int
            The prime dimension for modular arithmetic

        Returns
        -------
        TableauLike
            The full 2n x 2n symplectic matrix mod dimension
        """
        if isinstance(qudits, int):
            qudits = (qudits,)
        symplectic_matrix, _ = embed_symplectic(self.symplectic, self.phase_vector(dimension), qudits, n_qudits)
        if dimension is None:
            return symplectic_matrix

        return symplectic_matrix % dimension


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

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        from sympleq.core.circuits.utils import H_mat
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        U = H_mat(dimension)
        if self._is_inverse:
            return U.conj().T
        return U


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

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        from sympleq.core.circuits.utils import S_mat
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        U = S_mat(dimension)
        if self._is_inverse:
            return U.conj().T
        return U


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

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        """
        CX acts as |j,k⟩ -> |j, (j+k) mod d⟩ (or |j, (k-j) mod d⟩ for inverse).

        Parameters
        ----------
        dimension : int, optional
            The local Hilbert space dimension. Defaults to `DEFAULT_QUDIT_DIMENSION`.

        Returns
        -------
        HilbertOperator
            The d^2 x d^2 unitary matrix for the CX gate.
        """

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
        return HilbertOperator(U)


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

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        """
        SWAP acts as |j,k⟩ -> |k,j⟩.

        Parameters
        ----------
        dimension : int, optional
            The local Hilbert space dimension. Defaults to `DEFAULT_QUDIT_DIMENSION`.

        Returns
        -------
        HilbertOperator
            The d^2 x d^2 unitary matrix for the SWAP gate.
        """
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
        return HilbertOperator(U)

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

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        """
        CZ adds phase ω^{jk} to |j,k⟩ where ω = exp(2πi/d).

        Parameters
        ----------
        dimension : int, optional
            The local Hilbert space dimension. Defaults to `DEFAULT_QUDIT_DIMENSION`.

        Returns
        -------
        HilbertOperator
            The d^2 x d^2 unitary matrix for the CZ gate.
        """
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
        return HilbertOperator(U)

    def inverse(self) -> _CZ:
        # CZ is self-inverse
        return self


class _ZZMax(Gate):
    """ZZ-Phase gate (native for Quantinuum:
    https://docs.quantinuum.com/systems/trainings/helios/getting_started/parameterized_angle_2_qubit_gates.html).
    The angle is set to pi/4 thus the gate is Clifford."""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse

        if is_inverse:
            symplectic = np.array([
                [1, 0, -1, -1],   # image of X0:  X0 -> X0 Z0^{-1} Z1^{-1}
                [0, 1, -1, -1],    # image of X1:  X1 -> Z0 X1^{-1} Z1^{-1}
                [0, 0, 1, 0],    # image of Z0:  Z0 -> Z0
                [0, 0, 0, 1]     # image of Z1:  Z1 -> Z1
            ], dtype=int).T
            exceptional = {2: np.array([-1, -1, 0, 0], dtype=int)}
            name = "ZZP_inv"
        else:
            symplectic = np.array([
                [1, 0, 1, 1],    # image of X0:  X0 -> X0 Z0 Z1
                [0, 1, 1, 1],    # image of X1:  X1 -> Z0 X1 Z1
                [0, 0, 1, 0],    # image of Z0:  Z0 -> Z0
                [0, 0, 0, 1]    # image of Z1:  Z1 -> Z1
            ], dtype=int).T
            exceptional = {2: np.array([1, 1, 0, 0], dtype=int)}
            name = "ZZP"

        super().__init__(name, symplectic, exceptional_phase_vectors=exceptional)

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        """
        Applies phase exp(±iπ(j+k)²/d) to |j,k⟩ (sign flipped for the inverse).

        For qubits this reproduces ZZPhase(±π/4) = exp(∓iπ/4 Z⊗Z) up to a global phase.

        Parameters
        ----------
        dimension : int, optional
            The local Hilbert space dimension. Defaults to `DEFAULT_QUDIT_DIMENSION`.

        Returns
        -------
        HilbertOperator
            The d^2 x d^2 unitary matrix for the ZZMax gate.
        """
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        D = d * d
        sign = -1 if self._is_inverse else 1
        U = np.zeros((D, D), dtype=complex)
        for j in range(d):
            for k in range(d):
                idx = j * d + k
                U[idx, idx] = np.exp(sign * 1j * np.pi * (j + k) ** 2 / d)
        return HilbertOperator(U)


class _V(Gate):
    """V = √X gate: X -> X, Z -> -Y = -XZ. Has special phase vector for qubits.

    V is the X-axis analog of S: V = exp(-iπ/4 X). Together with S and any
    entangling Clifford, V generates the Clifford group, which
    makes ``{S, V, ZZMax}`` a useful generating set on Quantinuum H2 since
    each element maps 1:1 to a single H2 native gate (Rz(0.5), PhasedX(0.5, 0),
    ZZMax respectively).
    """

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse

        # Symplectic image is the same for V and V^{-1}; the phase vector differs.
        # X -> X (image (1, 0))
        # Z -> XZ (image (1, 1)); the sign of the XZ image is what distinguishes V from V^{-1}.
        symplectic = np.array([
            [1, 0],   # image of X:  X -> X
            [1, 1],   # image of Z:  Z -> XZ
        ], dtype=int).T

        if is_inverse:
            # V^{-1} Z V = +σ_y = +i · sympleq_Y, so phase[1] = +1 (=ω).
            exceptional = {2: np.array([0, 1], dtype=int)}
            name = "V_inv"
        else:
            # V Z V^{-1} = -σ_y = -i · sympleq_Y, so phase[1] = -1 (=ω^{-1}).
            exceptional = {2: np.array([0, -1], dtype=int)}
            name = "V"

        super().__init__(name, symplectic, exceptional_phase_vectors=exceptional)

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        if dimension != 2:
            # FIXME: An easy way to define it for qudit is to apply a Hadamard to the S gate:
            # H@S@H_inv
            raise NotImplementedError(
                "V (= √X) is only implemented for qubits (dimension=2)."
            )
        # √X = exp(-iπ/4 X) = (1/√2)(I - iX) = (1/√2)[[1, -i], [-i, 1]]
        sign = 1j if self._is_inverse else -1j
        U = np.array([[1, sign], [sign, 1]], dtype=complex) / np.sqrt(2)
        return HilbertOperator(U)


class _Id(Gate):
    """Identity gate: Id|j⟩ = |j⟩."""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse
        # NOTE: here we define the gate to be single-qudit, but overriding
        # the act method makes it work for any number of qudits.
        symplectic = np.eye(2, dtype=int)
        phase_vector = np.array([0, 0], dtype=int)

        super().__init__("Id", symplectic, phase_vector)

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        # X|j⟩ = |j+1 mod d⟩, X^{-1}|j⟩ = |j-1 mod d⟩
        U = np.eye(d, dtype=complex)
        return HilbertOperator(U)

    def inverse(self) -> _Id:
        # Id is self-inverse
        return self

    def act(self, pauli: PauliObject, qudits: int | tuple[int, ...]) -> PauliObject:
        return pauli


class _X(Gate):
    """Generalized X gate (shift operator): X|j⟩ = |j+1 mod d⟩."""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse
        symplectic = np.eye(2, dtype=int)

        if is_inverse:
            # X^{-1} = X^{d-1} has tableau [-1, 0]
            self._tableau = np.array([-1, 0], dtype=int)
            name = "X_inv"
        else:
            self._tableau = np.array([1, 0], dtype=int)
            name = "X"

        super().__init__(name, symplectic)

    def phase_vector(self, dimension: int | None = None) -> PhasesType:
        # h = 2 * Ω @ tableau, where Ω = [[0, 1], [-1, 0]]
        # Ω @ [x, 0] = [0, -x], so h = [0, -2x]
        x = self._tableau[0]
        if dimension is not None:
            return np.array([0, -2 * x], dtype=int) % (2 * dimension)
        return np.array([0, -2 * x], dtype=int)

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        # X|j⟩ = |j+1 mod d⟩, X^{-1}|j⟩ = |j-1 mod d⟩
        U = np.zeros((d, d), dtype=complex)
        for j in range(d):
            if self._is_inverse:
                U[(j - 1) % d, j] = 1.0
            else:
                U[(j + 1) % d, j] = 1.0
        return HilbertOperator(U)


class _Y(Gate):
    """Generalized Y gate: Y = X * Z."""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse
        symplectic = np.eye(2, dtype=int)

        if is_inverse:
            self._tableau = np.array([-1, -1], dtype=int)
            name = "Y_inv"
        else:
            self._tableau = np.array([1, 1], dtype=int)
            name = "Y"

        super().__init__(name, symplectic)

    def phase_vector(self, dimension: int | None = None) -> PhasesType:
        # h = 2 * Ω @ [x, z] = 2 * [z, -x]
        x, z = self._tableau
        if dimension is not None:
            return np.array([2 * z, -2 * x], dtype=int) % (2 * dimension)
        return np.array([2 * z, -2 * x], dtype=int)

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        omega = np.exp(2j * np.pi / d)
        # Y = X * Z: Y|j⟩ = ω^j |j+1 mod d⟩
        U = np.zeros((d, d), dtype=complex)
        for j in range(d):
            if self._is_inverse:
                # Y^{-1}|j⟩ = ω^{-(j-1)} |j-1 mod d⟩
                U[(j - 1) % d, j] = omega ** (-(j - 1) % d)
            else:
                U[(j + 1) % d, j] = omega ** j
        return HilbertOperator(np.around(U, 10))


class _Z(Gate):
    """Generalized Z gate (clock operator): Z|j⟩ = ω^j |j⟩."""

    def __init__(self, is_inverse: bool = False):
        self._is_inverse = is_inverse
        symplectic = np.eye(2, dtype=int)

        if is_inverse:
            self._tableau = np.array([0, -1], dtype=int)
            name = "Z_inv"
        else:
            self._tableau = np.array([0, 1], dtype=int)
            name = "Z"

        super().__init__(name, symplectic)

    def phase_vector(self, dimension: int | None = None) -> PhasesType:
        # h = 2 * Ω @ [0, z] = 2 * [z, 0]
        z = self._tableau[1]
        if dimension is not None:
            return np.array([2 * z, 0], dtype=int) % (2 * dimension)
        return np.array([2 * z, 0], dtype=int)

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        if dimension is None:
            dimension = DEFAULT_QUDIT_DIMENSION
        d = dimension
        omega = np.exp(2j * np.pi / d)
        # Z|j⟩ = ω^j |j⟩, Z^{-1}|j⟩ = ω^{-j} |j⟩
        if self._is_inverse:
            diag = [omega ** (-j % d) for j in range(d)]
        else:
            diag = [omega ** j for j in range(d)]
        return HilbertOperator(np.around(np.diag(diag), 10))


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

        self._ZZMax = _ZZMax(is_inverse=False)
        self._ZZMax_inv = _ZZMax(is_inverse=True)
        self._ZZMax._inverse = self._ZZMax_inv
        self._ZZMax_inv._inverse = self._ZZMax

        self._V = _V(is_inverse=False)
        self._V_inv = _V(is_inverse=True)
        self._V._inverse = self._V_inv
        self._V_inv._inverse = self._V

        self._SWAP = _SWAP()
        # SWAP is self-inverse, already handled in the class

        self._CZ = _CZ()
        # CZ is self-inverse, already handled in the class

        self._Id = _Id()

        # Pauli gates
        self._X = _X(is_inverse=False)
        self._X_inv = _X(is_inverse=True)
        self._X._inverse = self._X_inv
        self._X_inv._inverse = self._X

        self._Y = _Y(is_inverse=False)
        self._Y_inv = _Y(is_inverse=True)
        self._Y._inverse = self._Y_inv
        self._Y_inv._inverse = self._Y

        self._Z = _Z(is_inverse=False)
        self._Z_inv = _Z(is_inverse=True)
        self._Z._inverse = self._Z_inv
        self._Z_inv._inverse = self._Z

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

    # Quantinuum-ZZ-Phase with angle pi/4
    @property
    def ZZMax(self) -> _ZZMax:
        return self._ZZMax

    @property
    def ZZMax_inv(self) -> _ZZMax:
        return self._ZZMax_inv

    # V = √X (Quantinuum H2 native: PhasedX(0.5, 0))
    @property
    def V(self) -> _V:
        return self._V

    @property
    def V_inv(self) -> _V:
        return self._V_inv

    # SWAP
    @property
    def SWAP(self) -> _SWAP:
        return self._SWAP

    # CZ
    @property
    def CZ(self) -> _CZ:
        return self._CZ

    # Identity
    @property
    def Id(self) -> _Id:
        return self._Id

    # Pauli X
    @property
    def X(self) -> _X:
        return self._X

    @property
    def X_inv(self) -> _X:
        return self._X_inv

    # Pauli Y
    @property
    def Y(self) -> _Y:
        return self._Y

    @property
    def Y_inv(self) -> _Y:
        return self._Y_inv

    # Pauli Z
    @property
    def Z(self) -> _Z:
        return self._Z

    @property
    def Z_inv(self) -> _Z:
        return self._Z_inv


# Global singleton instance
GATES = _Gates()

# All built-in gates (forward and inverse variants).
DEFAULT_GATES_SET: list[Gate] = [
    GATES.Id,
    GATES.H, GATES.H_inv,
    GATES.S, GATES.S_inv,
    GATES.V, GATES.V_inv,
    GATES.X, GATES.X_inv,
    GATES.Y, GATES.Y_inv,
    GATES.Z, GATES.Z_inv,
    GATES.CX, GATES.CX_inv,
    GATES.SWAP,
    GATES.CZ,
    GATES.ZZMax, GATES.ZZMax_inv,
]


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
    def dimensions(self) -> DimensionsType:
        return self._dimensions

    def inverse(self) -> PauliGate:
        # Pauli gates are self-inverse (up to phase)
        return self

    def local_unitary(self, dimension: int | None = None) -> HilbertOperator:
        """
        Compute the unitary for this PauliGate.

        For PauliGate, dimension is optional since it's determined by the stored PauliString.

        Parameters
        ----------
        dimension : int, optional
            Unused; retained for interface compatibility with `Gate.local_unitary`.

        Returns
        -------
        HilbertOperator
            The unitary matrix for this PauliGate.
        """
        from sympleq.core.circuits.utils import pauli_unitary_from_tableau
        # Use the dimension from the PauliString
        d = int(self._dimensions[0])  # assumes uniform dimensions
        x = self.pauli_string.x_exp
        z = self.pauli_string.z_exp
        return pauli_unitary_from_tableau(d, x, z, convention="bare")

    @overload
    def act(self, pauli: PauliSum, qudits: int | tuple[int, ...]) -> PauliSum:
        ...

    @overload
    def act(self, pauli: PauliString, qudits: int | tuple[int, ...]) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliObject, qudits: int | tuple[int, ...]) -> PauliObject:
        ...

    def act(self, pauli, qudits):
        """
        Apply this PauliGate to a Pauli object.

        For PauliGate, qudits defaults to all qudits in order (0, 1, 2, ..., n-1)
        since the gate was constructed for a specific number of qudits.

        Parameters
        ----------
        pauli : Pauli | PauliString | PauliSum
            The Pauli object to transform.
        qudits : int | tuple[int, ...] | None
            The qudit index(es) the gate acts on. If None, defaults to all
            qudits in order.

        Returns
        -------
        PauliObject
            The transformed Pauli object of the same type as the input.
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
