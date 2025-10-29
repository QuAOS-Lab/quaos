from abc import ABC
import numpy as np
import scipy.sparse as sp
from typing import TypeVar, overload

from sympleq.core.circuits.find_symplectic import map_pauli_sum_to_target_tableau
from sympleq.core.paulis import Pauli, PauliString, PauliSum, PauliObject
from sympleq.core.circuits.target import get_phase_vector
from sympleq.core.circuits.utils import transvection_matrix, symplectic_form, tensor, I_mat, H_mat, S_mat, CX_func

# We define a type using TypeVar to let the type checker know that
# the input and output of the `act` function share the same type.
P = TypeVar("P", bound="PauliObject")


class Gate(ABC):
    def __init__(self, name: str, symplectic: np.ndarray):
        self.name = name
        n_qudits = len(symplectic) // 2
        self._n_qudits = n_qudits
        self.symplectic = symplectic

        # U = [[0_n, 0_n],
        #      [I_n, 0_n]]
        U = np.zeros((2 * self.n_qudits, 2 * self.n_qudits), dtype=int)
        U[self.n_qudits:, :self.n_qudits] = np.eye(self.n_qudits, dtype=int)
        self.U_symplectic_conjugated = self.symplectic.T @ U @ self.symplectic

        self.V_diag = np.diag(self.U_symplectic_conjugated)
        # This is the part associated with the quadratic form.
        # We remove diagonal part to match definition in Eq.[7] in PHYSICAL REVIEW A 71, 042315 (2005).
        self.p_part = 2 * np.triu(self.U_symplectic_conjugated) - np.diag(self.V_diag)

    @classmethod
    def get_phase_vector(cls, dimensions: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    @abstractmethod
    def n_qudits(self) -> int:
        pass

    @classmethod
    def solve_from_target(cls, name: str, input_pauli_sum: PauliSum, target_pauli_sum: PauliSum) -> 'Gate':
        """
        Create a gate that maps input_pauli_sum to target_pauli_sum.
        """

        input_symplectic = input_pauli_sum.tableau()
        target_symplectic = target_pauli_sum.tableau()
        symplectic = map_pauli_sum_to_target_tableau(input_symplectic, target_symplectic)

        return cls(name, symplectic.T)

    @classmethod
    def from_random(cls, n_qudits: int, n_transvection: int = 10, seed: int | None = None):
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
        return cls(f"R{n_transvection}", n_qudits, symplectic, phase_vector)

    def __repr__(self):
        return f"Gate(name={self.name}, n_qudits={self.n_qudits}, phase_vector={self.phase_vector})"

    @overload
    def act(self, pauli: Pauli, qudit_indices: int | list[int]) -> Pauli:
        ...

    @overload
    def act(self, pauli: PauliString, qudit_indices: int | list[int]) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliSum, qudit_indices: int | list[int]) -> PauliSum:
        ...

    def act(self, pauli: P, qudit_indices: int | list[int] | np.ndarray) -> P:
        """
        Returns the updated tableau and phases acquired by the PauliSum when acted upon by this gate.

        See Eq.[7] in PHYSICAL REVIEW A 71, 042315 (2005)

        """
        if isinstance(qudit_indices, int):
            qudit_indices = [qudit_indices]
        qudit_indices = np.asarray(qudit_indices, dtype=int)

        pauli_dimensions = pauli.dimensions()[qudit_indices]

        T = pauli.tableau()

        # Precompute tableau mask. This will be applied to the PauliSum tableau to get
        # the subset of affected columns.
        tableau_mask = np.concatenate([qudit_indices, qudit_indices + pauli.n_qudits()])

        T_affected = T[:, tableau_mask]
        relevant_dimensions = np.tile(pauli_dimensions, 2)
        updated_tableau = np.mod(T_affected @ self.symplectic.T, relevant_dimensions)
        new_tableau = T.copy()
        new_tableau[:, tableau_mask] = updated_tableau

        # FIXME: should we move this to a separate function?
        # This is the part associated with the linear form.
        phase_vector = Gate.get_phase_vector(pauli_dimensions)
        modified_phase_vector = phase_vector - self.V_diag  # h - V_diag
        linear_terms = T_affected @ modified_phase_vector
        quadratic_terms = np.sum(T_affected * (T_affected @ self.p_part), axis=1)

        # FIXME: this is a but of a hack
        dimensional_factor = pauli.lcm() // np.lcm.reduce(pauli.dimensions()[qudit_indices])
        acquired_phases = (linear_terms + quadratic_terms) * dimensional_factor

        new_phases = (pauli.phases() + acquired_phases) % (2 * pauli.lcm())

        return pauli.__class__(tableau=new_tableau, dimensions=pauli.dimensions(),
                               weights=pauli.weights(), phases=new_phases)

    def transvection(self, transvection_vector: np.ndarray | list, transvection_weight: int = 1) -> 'Gate':
        """
        Returns a new gate that is the transvection of this gate by the given vector.
        The transvection vector should be a 2n-dimensional vector where n is the number of qudits.
        """
        if not isinstance(transvection_weight, int) and not isinstance(transvection_weight, np.int64):
            raise TypeError("Transvection weight must be an integer.")

        if isinstance(transvection_vector, list):
            transvection_vector = np.array(transvection_vector)

        T = transvection_matrix(transvection_vector, multiplier=transvection_weight, p=dimension)
        if self.name[0] != "T":
            self.name = "T-" + self.name
        return Gate(self.name, self.symplectic @ T)

    def inv(self) -> 'Gate':
        # TODO: Test for mixed dimensions - not clear that the symplectic form here is correct.
        print("Warning: inverse phase vector not working - PHASES MAY BE INCORRECT.")

        C = self.symplectic.T

        U = np.zeros((2 * self.n_qudits, 2 * self.n_qudits), dtype=int)
        U[self.n_qudits:, :self.n_qudits] = np.eye(self.n_qudits, dtype=int)
        Omega = symplectic_form(int(C.shape[0] / 2), p=self.lcm)

        C_inv = -(Omega.T @ C.T @ Omega) % self.lcm
        U_c = C_inv.T @ U @ C_inv % self.lcm

        p1 = - C_inv.T @ self.phase_vector
        p2 = - np.diag(C.T @ (2 * np.triu(U_c) - np.diag(np.diag(U_c))) @ C)
        p3 = C.T @ np.diag(U_c)

        phase_vector = (p1 + p2 + p3) % (2 * self.lcm)
        return Gate(self.name + "-inv", self.qudit_indices, C_inv.T, phase_vector)

    def inverse(self) -> 'Gate':
        """
        Returns the inverse of this gate.

        The inverse of a gate G is another gate G' such that the composition G'G is the identity.

        The inverse of a gate is computed using the formulae presented in PHYSICAL REVIEW A 71, 042315 (2005)

        :return: A new Gate object, the inverse of this gate.
        """
        n = self.n_qudits
        Id_n = np.eye(n)
        Zero_n = np.zeros((n, n))
        U = np.block([[Zero_n, Zero_n], [Id_n, Zero_n]])
        Omega = (U - U.T)

        C = self.symplectic.T
        C_inv = Omega.T @ C.T @ Omega

        return Gate(self.name + '_inv', C_inv.T)

    def unitary(self, qudit_indices: list[int], dimensions: list[int] | np.ndarray) -> sp.csr_matrix:
        raise NotImplementedError("Unitary not implemented for generic Gate. Use specific gate subclasses.")


class SUM(Gate):
    def __init__(self):
        symplectic = np.array([
            [1, 1, 0, 0],   # image of X0:  X0 -> X0 X1
            [0, 1, 0, 0],   # image of X1:  X1 -> X1
            [0, 0, 1, 0],   # image of Z0:  Z0 -> Z0
            [0, 0, -1, 1]   # image of Z1:  Z1 -> Z0^-1 Z1
        ], dtype=int).T

        super().__init__("SUM", symplectic)

    @classmethod
    def get_phase_vector(cls, dimensions: np.ndarray) -> np.ndarray:
        return np.array([0, 0, 0, 0], dtype=int)

    def unitary(self, qudit_indices: list[int], dimensions: list[int] | np.ndarray) -> sp.csr_matrix:
        D = np.prod(dimensions)
        aa = qudit_indices
        a0 = aa[0]
        a1 = aa[1]
        aa2 = np.array([1 for i in range(D)])
        aa3 = np.array([CX_func(i, a0, a1, dimensions) for i in range(D)])
        aa4 = np.array([i for i in range(D)])
        return sp.csr_matrix((aa2, (aa3, aa4)))


class SWAP(Gate):
    def __init__(self):
        symplectic = np.array([
            [0, 1, 0, 0],  # image of X0:  X0 -> X1
            [1, 0, 0, 0],  # image of X1:  X1 -> X0
            [0, 0, 0, 1],  # image of Z0:  Z0 -> Z1
            [0, 0, 1, 0]   # image of Z1:  Z1 -> Z0
        ], dtype=int).T

        phase_vector = np.array([0, 0, 0, 0], dtype=int)

        super().__init__("SWAP", 2, symplectic, phase_vector=phase_vector)

    def unitary(self, qudit_indices: list[int], dimensions: list[int] | np.ndarray) -> sp.csr_matrix:
        # SWAP on two qudits of equal dimension: |i, j> -> |j, i>.
        # Basis ordering |i>⊗|j> with linear index idx(i, j) = i * d + j.
        aa = qudit_indices
        q = len(dimensions)
        D = np.prod(dimensions)
        a0 = q - 1 - aa[0]
        a1 = q - 1 - aa[1]
        aa2 = np.array([1 for i in range(D)])
        aa3 = np.array([i for i in range(D)])
        aa4 = np.array([SWAP._swap_linear_index(i, a0, a1, dimensions) for i in range(D)])
        return sp.csr_matrix((aa2, (aa3, aa4)))

    @staticmethod
    def _swap_linear_index(i, a0, a1, dims):
        # Convert linear index i to multi-index in row-major order
        multi_idx = []
        rem = i
        for d in reversed(dims):
            multi_idx.append(rem % d)
            rem //= d
        multi_idx = multi_idx[::-1]

        # Swap the two qudits
        multi_idx[a0], multi_idx[a1] = multi_idx[a1], multi_idx[a0]

        # Convert back to linear index (row-major)
        idx = 0
        for j, dim in enumerate(dims):
            idx = idx * dim + multi_idx[j]
        return idx


class CNOT(Gate):
    def __init__(self):
        symplectic = np.array([
            [1, 1, 0, 0],   # image of X0:  X0 -> X0 X1
            [0, 1, 0, 0],   # image of X1:  X1 -> X1
            [0, 0, 1, 0],   # image of Z0:  Z0 -> Z0
            [0, 0, 1, 1]   # image of Z1:  Z1 -> Z0^-1 Z1
        ], dtype=int).T

        phase_vector = np.array([0, 0, 0, 0], dtype=int)

        super().__init__("SUM", 2, symplectic, phase_vector=phase_vector)

    def unitary(self, qudit_indices: list[int], dimensions: list[int] | np.ndarray) -> sp.csr_matrix:
        D = np.prod(dimensions)
        aa = qudit_indices
        a0 = aa[0]
        a1 = aa[1]
        aa2 = np.array([1 for i in range(D)])
        aa3 = np.array([CX_func(i, a0, a1, dimensions) for i in range(D)])
        aa4 = np.array([i for i in range(D)])
        return sp.csr_matrix((aa2, (aa3, aa4)))


class Hadamard(Gate):
    def __init__(self, inverse: bool = False):
        if inverse:
            symplectic = np.array([
                [0, 1],    # image of X:  X -> Z
                [-1, 0]    # image of Z:  Z -> -X
            ], dtype=int)
        else:
            symplectic = np.array([
                [0, -1],   # image of X:  X -> -Z
                [1, 0]     # image of Z:  Z -> X
            ], dtype=int)

        phase_vector = np.array([0, 0], dtype=int)

        name = "H" if not inverse else "H_inv"
        super().__init__(name, 1, symplectic, phase_vector=phase_vector)

    def unitary(self, qudit_indices: list[int], dimensions: list[int] | np.ndarray) -> sp.csr_matrix:
        return tensor(
            [H_mat(dimensions[i]) if i in qudit_indices else I_mat(dimensions[i]) for i in range(len(dimensions))]
        )


class PHASE(Gate):
    def __init__(self):
        symplectic = np.array([
            [1, 1],  # image of X:  X -> XZ
            [0, 1]   # image of Z:  Z -> Z
        ], dtype=int).T

        super().__init__("S", 1, symplectic)

    @classmethod
    def get_phase_vector(cls, dimensions: np.ndarray) -> np.ndarray:
        if dimensions == 2:
            phase_vector = np.array([1, 0], dtype=int)
        else:
            phase_vector = np.array([0, 0], dtype=int)
        return super().get_phase_vector()

    def unitary(self, qudit_indices: list[int], dimensions: list[int] | np.ndarray) -> sp.csr_matrix:
        unitary = tensor([S_mat(dimensions[i]) if i in qudit_indices else I_mat(dimensions[i])
                         for i in range(len(dimensions))])
        return unitary


class _Gates():
    def __init__(self):
        self._sum = SUM()
        self._swap = SWAP()
        self._cnot = CNOT()
        self._hadamard = Hadamard()
        self._phase = PHASE()

    @property
    def sum(self):
        return self._sum

    @property
    def CX(self):
        return self._sum

    @property
    def swap(self):
        return self._swap

    @property
    def cnot(self):
        return self._cnot

    @property
    def hadamard(self):
        return self._hadamard

    @property
    def H(self):
        return self._hadamard

    @property
    def phase(self):
        return self._phase

    @property
    def S(self):
        return self._phase


GATES = _Gates()
