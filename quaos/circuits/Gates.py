import numpy as np
from quaos.paulis.pauli_string import PauliString
from quaos.paulis.pauli_sum import PauliSum
from sympy import Matrix
from typing import List, Optional


class Gate:
    """
    Mapping can be written as set of rules,
    
    e.g. for CNOT
                x1z0*x0z0 -> x1z0*x1z0
                x0z0*x1z0 -> x0z0*x1z0  # doesn't need specifying
                x0z1*x0z0 -> x0z1*x0z0  # doesn't need specifying
                x0z0*x0z1 -> x0z-1*x0z1 # 

    inputs are:

    mapping = ['x1z0*x0z0 -> x1z0*x1z0', 'x0z0*x0z1 -> -x0z1*x0z1']  # (control*target -> control*target)

    """
    def __init__(self, name: str, qudit_indices: list[int], images: list[str], dimension: list[int]):
        self.dimension = dimension
        self.name = name
        self.qudit_indices = qudit_indices
        self.images = images
        self.n_qudits = len(qudit_indices)
        symplectic, phase_function = compute_symplectic_and_phase_function(images, dimension)
        self.symplectic = symplectic
        self.phase_function

    def act_on_pauli_string(self, P: PauliString):
        # Extract the relevant symplectics from the PauliString corresponding to self.qudit_indices
        

        # Act on the symplectic of the PauliString with the symplectic of the gate


        # Replace the symplectics in the positions in self.qudit_indices


        return 
    

    def act_on_pauli_sum(self, P: PauliSum):
        return
    


    def act(self, P: PauliString | PauliSum):
        return
    


    
def solve_underdetermined_symplectic_mapping(A: np.ndarray, B: np.ndarray, d: int) -> Optional[np.ndarray]:
    try:
        n = A.shape[1] // 2
        dim = 2 * n
        Omega = symplectic_form(n, d)

        A_mat = Matrix(A.tolist()).applyfunc(lambda x: x % d)
        B_mat = Matrix(B.tolist()).applyfunc(lambda x: x % d)

        basis_A = A.tolist()
        while len(basis_A) < dim:
            for i in range(dim):
                candidate = np.eye(dim, dtype=int)[i]
                if all(symplectic_inner(candidate, np.array(v), Omega, d) == 0 for v in basis_A):
                    basis_A.append(candidate.tolist())
                    break
            else:
                raise ValueError("Could not extend A to full symplectic basis.")

        extra = len(basis_A) - len(B)
        basis_B = B.tolist() + np.eye(dim, dtype=int)[:extra].tolist()

        A_full = Matrix(basis_A).applyfunc(lambda x: x % d)
        B_full = Matrix(basis_B).applyfunc(lambda x: x % d)
        S = A_full.solve_least_squares(B_full, method='LDL')
        S = np.array(S.applyfunc(lambda x: x % d)).astype(int)

        if (S.T @ Omega @ S % d == Omega).all():  # Check if S is symplectic
            return S
        else:
            return None

    except Exception:
        print("Warning: Could not solve underdetermined symplectic mapping.")
        return None


def compute_symplectic_from_images(images: List[np.ndarray], d: int) -> Optional[np.ndarray]:
    n = len(images) // 2
    dim = 2 * n
    if len(images) != dim:
        raise ValueError("Need exactly 2n image vectors for n-qudit Clifford.")

    standard_basis = np.eye(dim, dtype=int)
    A = np.stack(standard_basis)
    B = np.stack([v % d for v in images])
    return solve_underdetermined_symplectic_mapping(A, B, d)


def compute_symplectic_and_phase_function(images: List[np.ndarray], d: int):
    """
    Computes the symplectic matrix S and a phase function for a Clifford gate defined by Pauli generator images.

    Parameters:
    - images: List of length-2n vectors specifying the images of [X0, ..., Xn-1, Z0, ..., Zn-1] under the gate
    - d: Dimension of the qudit system

    Returns:
    - S: symplectic matrix ∈ Sp(2n, Z_d)
    - phase_fn: function accepting symplectic vector v and returning the induced phase (in ω^k form)
    """
    n = len(images) // 2
    dim = 2 * n
    Omega = symplectic_form(n, d)

    S = compute_symplectic_from_images(images, d)
    if S is None:
        return None, None

    # Define phase function: maps input Pauli vector to its acquired phase (as exponent of ω)
    def phase_fn(pauli_vec: np.ndarray) -> int:
        """
        Computes the phase exponent k such that:
            U P(v) U† = ω^k P(S v)
        where ω = exp(2πi / d), and v is the symplectic representation of the Pauli operator.

        Parameters:
        - pauli_vec: A symplectic vector v ∈ Z_d^{2n}

        Returns:
        - Integer k ∈ Z_d representing the phase exponent ω^k
        """
        v = pauli_vec % d
        Sv = S @ v % d

        # Symplectic inner product (v, Sv) / 2 mod d gives phase
        # This works for standard Clifford gates; phase is linear mod d
        k = symplectic_inner(v, Sv, Omega, d)
        if d % 2 == 0:
            k = (k * pow(2, -1, d)) % d  # divide by 2 mod d
        else:
            k = (k * pow(2, -1, d)) % d
        return k

    return S, phase_fn


class SUM(Gate):
    def __init__(self, control, target, dimension):
        images = [np.array([1, 1, 0, 0]),  # image of X0:  X0 -> X0 X1
                  np.array([0, 1, 0, 0]),  # image of X1:  X1 -> X1
                  np.array([0, 0, 1, 0]),  # image of Z0:  Z0 -> Z0
                  np.array([0, 0, -1, 1])  # image of Z1:  Z1 -> Z0^-1 Z1
                  ]
        super().__init__("SUM", [control, target], images, dimension=dimension)
    

def symplectic_form(n: int, d: int) -> np.ndarray:
    I = np.eye(n, dtype=int)
    O = np.zeros((n, n), dtype=int)
    return np.block([[O, I], [-I % d, O]])


def symplectic_inner(a: np.ndarray, b: np.ndarray, Omega: np.ndarray, d: int) -> int:
    return int(np.dot(a, Omega @ b) % d)
