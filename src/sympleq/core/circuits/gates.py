import numpy as np
from sympleq.core.paulis import PauliString, PauliSum, Pauli
from typing import overload
from sympleq.core.circuits.target import find_map_to_target_pauli_sum, get_phase_vector
from sympleq.core.circuits.utils import (transvection_matrix, symplectic_form, tensor, I_mat, H_mat, S_mat, CX_func,
                                         SWAP_func, pauli_unitary_from_tableau)
from sympleq.core.finite_field_solvers import get_linear_dependencies
import scipy.sparse as sp
from .utils import embed_symplectic


class Gate:
    def __init__(self, name: str,
                 qudit_indices: list[int] | np.ndarray,
                 symplectic: np.ndarray,
                 dimensions: int | list[int] | np.ndarray,
                 phase_vector: np.ndarray | list[int]):

        if len(qudit_indices) == 0:
            raise ValueError("Gate must act on at least one qudit_indices.")

        if isinstance(dimensions, int) or isinstance(dimensions, np.signedinteger):
            dimensions = dimensions * np.ones(len(qudit_indices), dtype=int)
        elif len(qudit_indices) != len(dimensions):
            raise ValueError("Dimensions and qudit_indices must have the same length.")

        self.dimensions = np.asarray(dimensions, int)
        self.dimensions.setflags(write=False)
        self.name = name

        if isinstance(qudit_indices, list):
            qudit_indices = np.asarray(qudit_indices, dtype=int)

        self.qudit_indices = qudit_indices
        self.qudit_indices.setflags(write=False)
        self.n_qudits = len(qudit_indices)
        self.symplectic: np.ndarray = symplectic
        self.symplectic.setflags(write=False)
        self.phase_vector = np.asarray(phase_vector, dtype=int)
        self.phase_vector.setflags(write=False)
        self.lcm = np.lcm.reduce(self.dimensions)

        # U = [[0_n, 0_n],
        #      [I_n, 0_n]]
        U = np.zeros((2 * self.n_qudits, 2 * self.n_qudits), dtype=int)
        U[self.n_qudits:, :self.n_qudits] = np.eye(self.n_qudits, dtype=int)
        self.U_symplectic_conjugated = self.symplectic.T @ U @ self.symplectic

        self.V_diag = np.diag(self.U_symplectic_conjugated)
        # This is the part associated with the linear form.
        self.modified_phase_vector = self.phase_vector - self.V_diag  # h - V_diag
        # Remove diagonal part to match definition in Eq.[7] in PHYSICAL REVIEW A 71, 042315 (2005).
        # This is the part associated with the quadratic form.
        self.p_part = 2 * np.triu(self.U_symplectic_conjugated) - np.diag(self.V_diag)

    @classmethod
    def solve_from_target(cls, name: str, input_pauli_sum: PauliSum, target_pauli_sum: PauliSum,
                          dimensions: int | list[int] | np.ndarray):
        """
        Create a gate that maps input_pauli_sum to target_pauli_sum.
        """

        independent_set, dependent_set = get_linear_dependencies(input_pauli_sum.tableau, dimensions)

        if len(dependent_set) != 0:
            raise NotImplementedError("Input PauliSum is not linearly independent. Will be implemented dreckly.")

        symplectic, phase_vector, qudit_indices, dimension = find_map_to_target_pauli_sum(input_pauli_sum,
                                                                                          target_pauli_sum)
        return cls(name, qudit_indices, symplectic.T, dimension, phase_vector)

    @classmethod
    def from_random(cls, n_qudits: int, dimension: int, n_transvection: int = 10, seed: int | None = None):
        if seed is not None:
            np.random.seed(seed)
        seed_vec = np.random.randint(0, 100000, size=n_transvection)
        # if n_qudits < 4 and dimension == 2:
        #     symp_int = np.random.randint(symplectic_group_size(n_qudits))
        #     symplectic = symplectic_gf2(symp_int, n_qudits)
        #     phase_vector = get_phase_vector(symplectic, dimension)
        #     return cls(f"R{symp_int}", list(range(n_qudits)), symplectic.T, dimension, phase_vector)
        # else:
        symplectic = np.eye(2 * n_qudits, dtype=int)

        for i in range(n_transvection):
            np.random.seed(seed_vec[i])
            Tv = transvection_matrix(np.random.randint(0, dimension, size=2 * n_qudits), dimension) % dimension
            symplectic = symplectic @ Tv % dimension

        phase_vector = get_phase_vector(symplectic, dimension)
        return cls(f"R{n_transvection}", list(range(n_qudits)), symplectic, dimension, phase_vector)

    def __repr__(self):
        return f"Gate(name={self.name}, qudit_indices={self.qudit_indices}, " \
            f"dimensions={self.dimensions}, phase_vector={self.phase_vector})"

    @overload
    def act(self, pauli: Pauli) -> Pauli:
        ...

    @overload
    def act(self, pauli: PauliString) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliSum) -> PauliSum:
        ...

    def act(self, pauli):
        """
        Returns the updated tableau and phases acquired by the PauliSum when acted upon by this gate.

        See Eq.[7] in PHYSICAL REVIEW A 71, 042315 (2005)

        """
        if not np.array_equal(self.dimensions, pauli.dimensions[self.qudit_indices]):
            raise ValueError("Gate and Pauli object have different dimensions.")

        T = pauli.tableau

        # Precompute tableau mask. This will be applied to the PauliSum tableau to get
        # the subset of affected columns.
        tableau_mask = np.concatenate([self.qudit_indices, self.qudit_indices + pauli.n_qudits()])

        T_affected = T[:, tableau_mask]
        relevant_dimensions = np.tile(pauli.dimensions[self.qudit_indices], 2)
        updated_tableau = np.mod(T_affected @ self.symplectic.T, relevant_dimensions)
        new_tableau = T.copy()
        new_tableau[:, tableau_mask] = updated_tableau

        # FIXME: should we move this to a separate function?
        linear_terms = T_affected @ self.modified_phase_vector
        quadratic_terms = np.sum(T_affected * (T_affected @ self.p_part), axis=1)

        # FIXME: this is a but of a hack
        dimensional_factor = pauli.lcm // np.lcm.reduce(pauli.dimensions[self.qudit_indices])
        acquired_phases = (linear_terms + quadratic_terms) * dimensional_factor

        new_phases = (pauli.phases + acquired_phases) % (2 * pauli.lcm)

        return pauli.__class__(tableau=new_tableau, dimensions=pauli.dimensions,
                               weights=pauli.weights, phases=new_phases)

    def copy(self) -> 'Gate':
        """
        Returns a copy of the gate.
        """
        return Gate(self.name, self.qudit_indices.copy(), self.symplectic.copy(), self.dimensions,
                    self.phase_vector.copy())

    def inv(self) -> 'Gate':
        n = self.n_qudits
        dims = np.asarray(self.dimensions, dtype=int)
        L = int(np.lcm.reduce(dims))
        modulus = 2 * L

        C = self.symplectic % L

        zero_block = np.zeros((n, n), dtype=int)
        identity_block = np.eye(n, dtype=int)
        Omega = np.block([[zero_block, identity_block], [-identity_block, zero_block]])

        C_inv = (Omega.T @ C.T @ Omega) % L
        C_inv = C_inv.astype(int)

        U = np.zeros((2 * n, 2 * n), dtype=int)
        U[n:, :n] = np.eye(n, dtype=int)

        U_C = (C.T @ U @ C) % modulus
        U_C_inv = (C_inv.T @ U @ C_inv) % modulus

        P_C = (2 * np.triu(U_C) - np.diag(np.diag(U_C))) % modulus
        P_C_inv = (2 * np.triu(U_C_inv) - np.diag(np.diag(U_C_inv))) % modulus

        h = (self.phase_vector % modulus).astype(int)

        term1 = (-h @ C_inv) % modulus
        term2 = (np.diag(U_C.T) % modulus) @ C_inv % modulus
        term3 = np.diag((C_inv.T @ P_C @ C_inv) % modulus) % modulus
        term4 = np.diag(U_C_inv) % modulus
        term5 = np.diag(P_C_inv) % modulus

        h_inv = (term1 + term2 - term3 + term4 - term5) % modulus

        return Gate(self.name + "-inv", self.qudit_indices, C_inv, self.dimensions, h_inv.astype(int))

    def transvection(self, transvection_vector: np.ndarray | list, transvection_weight: int = 1) -> 'Gate':
        """
        Returns a new gate that is the transvection of this gate by the given vector.
        The transvection vector should be a 2n-dimensional vector where n is the number of qudits.
        """
        if not np.all(self.dimensions == self.dimensions[0]):
            raise ValueError("Transvections only implemented for gates with equal dimensions.")
        dimension = self.dimensions[0]
        if transvection_weight >= dimension:
            raise ValueError("Transvection weight must be less than the gate dimension.")
        if not isinstance(transvection_weight, int) and not isinstance(transvection_weight, np.int64):
            raise TypeError("Transvection weight must be an integer.")
        if isinstance(transvection_vector, list):
            transvection_vector = np.array(transvection_vector)

        T = transvection_matrix(transvection_vector, multiplier=transvection_weight, p=dimension)
        name = self.name if self.name.startswith("T-") else "T-" + self.name
        return Gate(name, self.qudit_indices, self.symplectic @ T, self.dimensions, self.phase_vector)

    def unitary(self, dims=None) -> sp.csr_matrix:
        if dims is None:
            dims = self.dimensions
        raise NotImplementedError("Unitary not implemented for generic Gate. Use specific gate subclasses.")

    def full_symplectic(self, n_qudits):
        if n_qudits < max(self.qudit_indices):
            raise ValueError("n_qudits must be greater than or equal to the maximum qudit index.")
        full_symplectic, _ = embed_symplectic(self.symplectic, self.phase_vector, self.qudit_indices, n_qudits)
        return full_symplectic

    def __hash__(self):
        return hash((
            self.name,
            tuple(self.qudit_indices.tobytes()),
            tuple(self.dimensions.tobytes()),
            tuple(self.symplectic.flatten().tobytes()),
            tuple(self.phase_vector.tobytes()),
        ))

    def __eq__(self, other):
        if not isinstance(other, Gate):
            return False
        return np.all(self.qudit_indices == other.qudit_indices) and \
            np.all(self.symplectic == other.symplectic) and np.all(self.dimensions == other.dimensions) and \
            np.all(self.phase_vector == other.phase_vector)


def _scalar_dim(dim):
    return int(np.asarray(dim).reshape(-1)[0])


class SUM(Gate):
    def __init__(self, control, target, dimension):
        symplectic = np.array([
            [1, 1, 0, 0],   # image of X0:  X0 -> X0 X1
            [0, 1, 0, 0],   # image of X1:  X1 -> X1
            [0, 0, 1, 0],   # image of Z0:  Z0 -> Z0
            [0, 0, -1, 1]   # image of Z1:  Z1 -> Z0^-1 Z1
        ], dtype=int).T % dimension

        phase_vector = np.array([0, 0, 0, 0], dtype=int)

        super().__init__("SUM", [control, target], symplectic, dimensions=dimension, phase_vector=phase_vector)

    def unitary(self, dims=None) -> sp.csr_matrix:
        if dims is None:
            dims = self.dimensions
        D = np.prod(dims)
        aa = self.qudit_indices
        a0 = aa[0]
        a1 = aa[1]
        aa2 = np.array([1 for i in range(D)])
        aa3 = np.array([CX_func(i, a0, a1, dims) for i in range(D)])
        aa4 = np.array([i for i in range(D)])
        return sp.csr_matrix((aa2, (aa3, aa4)))

    def copy(self) -> 'Gate':
        d = _scalar_dim(self.dimensions)
        return SUM(int(self.qudit_indices[0]), int(self.qudit_indices[1]), d)


class CZ(Gate):
    def __init__(self, index1, index2, dimension):
        symplectic = np.array([
            [1, 0, 0, 0],  # image of X0:  X0 -> X0
            [0, 1, 0, 0],  # image of X1:  X1 -> X1
            [0, 1, 1, 0],  # image of Z0:  Z0 -> Z0
            [1, 0, 0, 1]   # image of Z1:  Z1 -> Z1
        ], dtype=int).T

        phase_vector = np.array([0, 0, 0, 0], dtype=int)

        super().__init__("CZ", [index1, index2], symplectic, dimensions=dimension, phase_vector=phase_vector)

    def copy(self) -> 'Gate':
        d = _scalar_dim(self.dimensions)
        return CZ(int(self.qudit_indices[0]), int(self.qudit_indices[1]), d)


class SWAP(Gate):
    def __init__(self, index1, index2, dimension):
        symplectic = np.array([
            [0, 1, 0, 0],  # image of X0:  X0 -> X1
            [1, 0, 0, 0],  # image of X1:  X1 -> X0
            [0, 0, 0, 1],  # image of Z0:  Z0 -> Z1
            [0, 0, 1, 0]   # image of Z1:  Z1 -> Z0
        ], dtype=int).T

        phase_vector = np.array([0, 0, 0, 0], dtype=int)

        super().__init__("SWAP", [index1, index2], symplectic, dimensions=dimension, phase_vector=phase_vector)

    def unitary(self, dims=None) -> sp.csr_matrix:
        """
        SWAP on two qudits at positions self.qudit_indices = [a0, a1]
        for an arbitrary mixed-radix register with local dims.
        Builds permutation matrix P with P_{f(i), i} = 1 where f applies the swap.
        """
        if dims is None:
            dims = self.dimensions
        dims = np.asarray(dims, dtype=int)

        a0, a1 = map(int, self.qudit_indices)  # expect 0-based positions
        q = len(dims)
        if not (0 <= a0 < q and 0 <= a1 < q and a0 != a1):
            raise ValueError("Invalid qudit indices to swap.")

        D = int(np.prod(dims))

        # Columns are original indices 0..D-1; rows are mapped indices f(i)
        cols = np.arange(D, dtype=int)
        rows = np.fromiter((SWAP_func(i, a0, a1, dims) for i in cols), count=D, dtype=int)

        data = np.ones(D, dtype=int)
        return sp.csr_matrix((data, (rows, cols)), shape=(D, D))

    def copy(self) -> 'Gate':
        d = _scalar_dim(self.dimensions)
        return SWAP(int(self.qudit_indices[0]), int(self.qudit_indices[1]), d)


class CNOT(Gate):
    def __init__(self, control, target):
        symplectic = np.array([
            [1, 1, 0, 0],   # image of X0:  X0 -> X0 X1
            [0, 1, 0, 0],   # image of X1:  X1 -> X1
            [0, 0, 1, 0],   # image of Z0:  Z0 -> Z0
            [0, 0, 1, 1]   # image of Z1:  Z1 -> Z0^-1 Z1
        ], dtype=int).T

        phase_vector = np.array([0, 0, 0, 0], dtype=int)

        super().__init__("SUM", [control, target], symplectic, dimensions=2, phase_vector=phase_vector)

    def unitary(self, dims=None) -> sp.csr_matrix:
        if dims is None:
            dims = self.dimensions
        D = np.prod(dims)
        aa = self.qudit_indices
        a0 = aa[0]
        a1 = aa[1]
        aa2 = np.array([1 for i in range(D)])
        aa3 = np.array([CX_func(i, a0, a1, dims) for i in range(D)])
        aa4 = np.array([i for i in range(D)])
        return sp.csr_matrix((aa2, (aa3, aa4)))

    def copy(self) -> 'Gate':
        return CNOT(int(self.qudit_indices[0]), int(self.qudit_indices[1]))


class Hadamard(Gate):
    def __init__(self, index: int, dimension: int, inverse: bool = False):
        self.is_inverse = bool(inverse)
        if inverse:
            symplectic = np.array([
                [0, 1],    # image of X:  X -> Z
                [-1, 0]    # image of Z:  Z -> -X
            ], dtype=int) % dimension
        else:
            symplectic = np.array([
                [0, -1],   # image of X:  X -> -Z
                [1, 0]     # image of Z:  Z -> X
            ], dtype=int) % dimension

        phase_vector = np.array([0, 0], dtype=int)

        name = "H" if not inverse else "H_inv"
        super().__init__(name, [index], symplectic, dimensions=dimension, phase_vector=phase_vector)

    def unitary(self, dims=None) -> sp.csr_matrix:
        if dims is None:
            dims = self.dimensions
        return tensor([H_mat(dims[i]) if i in self.qudit_indices else I_mat(dims[i]) for i in range(len(dims))])

    def copy(self) -> 'Gate':
        d = _scalar_dim(self.dimensions)
        return Hadamard(int(self.qudit_indices[0]), d, inverse=self.is_inverse)


class PHASE(Gate):

    def __init__(self, index: int, dimension: int):
        symplectic = np.array([
            [1, 1],  # image of X:  X -> XZ
            [0, 1]   # image of Z:  Z -> Z
        ], dtype=int).T
        if dimension == 2:
            phase_vector = np.array([1, 0], dtype=int)
        else:
            phase_vector = np.array([0, 0], dtype=int)

        super().__init__("S", [index], symplectic, dimensions=dimension, phase_vector=phase_vector)

    def unitary(self, dims=None) -> sp.csr_matrix:
        if dims is None:
            dims = self.dimensions

        unitary = tensor([S_mat(dims[i]) if i in self.qudit_indices else I_mat(dims[i]) for i in range(len(dims))])
        return unitary

    def copy(self) -> 'Gate':
        d = _scalar_dim(self.dimensions)
        return PHASE(int(self.qudit_indices[0]), d)


class PauliGate(Gate):
    def __init__(self, pauli: PauliString):
        self.pauli_string = pauli
        n = pauli.n_qudits()
        lcm = int(pauli.lcm)
        symplectic = np.eye(2 * n, dtype=int)
        phase_vector = (2 * symplectic_form(n, lcm) @ np.concatenate([pauli.x_exp, pauli.z_exp])) % (2 * lcm)
        super().__init__("Pauli", list(range(n)), symplectic, dimensions=pauli.dimensions, phase_vector=phase_vector)

    def copy(self) -> 'Gate':
        return PauliGate(self.pauli_string)

    def unitary(self, dims=None) -> sp.csr_matrix:
        if dims is None:
            dims = self.dimensions
        return pauli_unitary_from_tableau(dims[0], self.pauli_string.x_exp, self.pauli_string.z_exp)
