import numpy as np
from quaos.paulis import PauliString, PauliSum, Pauli


class Gate:
    """
    Mapping can be written as set of rules,
    
    inputs are:

    """
    def __init__(self, name: str, qudit_indices: list[int], images: list[np.ndarray], dimension: int,
                 phase_matrix: np.ndarray | None = None):
        self.dimension = dimension
        self.name = name
        self.qudit_indices = qudit_indices
        self.images = images
        self.n_qudits = len(qudit_indices)
        self.symplectic = np.stack([v % dimension for v in images]).T
        if phase_matrix is None:
            self.phase_matrix = make_phase_matrix_from_symplectic(self.symplectic, self.dimension)
        else:
            self.phase_matrix = phase_matrix
        # self.phase_function = phase_function

    def _act_on_pauli_string(self, P: PauliString) -> tuple[PauliString, int]:
        local_symplectic = np.concatenate([P.x_exp[self.qudit_indices], P.z_exp[self.qudit_indices]])
        acquired_phase = ((local_symplectic.T @ self.phase_matrix @ local_symplectic) // 2) % self.dimension

        local_symplectic = (self.symplectic @ local_symplectic) % self.dimension
        P = P._replace_symplectic(local_symplectic, self.qudit_indices)
        return P, acquired_phase
    
    def _act_on_pauli_sum(self, P: PauliSum):
        pauli_strings, phases = zip(*[self._act_on_pauli_string(p) for p in P.pauli_strings])
        P = PauliSum(pauli_strings, P.weights, np.asarray(phases), P.dimensions, False)
        return P

    def act(self, P: Pauli | PauliString | PauliSum):
        if isinstance(P, Pauli):
            P = P._to_pauli_string()
        if isinstance(P, PauliString):
            return self._act_on_pauli_string(P)[0]
        elif isinstance(P, PauliSum):
            return self._act_on_pauli_sum(P)
        else:
            raise TypeError(f"Unsupported type {type(P)} for Gate.act. Expected Pauli, PauliString or PauliSum.")


def make_phase_matrix_from_symplectic(symplectic: np.ndarray, dimension: int) -> np.ndarray:
    """
    Note that this is not unique! For example putting the SUM gate symplectic in here will not return the phase matrix
    of the SUM gate.

    This gives a consistent choice of the phase matrix only, as it is not uniquely determined by the symplectic.
    """
    n = symplectic.shape[0] // 2
    # Define standard symplectic form Omega
    id = np.eye(n, dtype=int)
    zero_mat = np.zeros((n, n), dtype=int)
    Omega = np.block([[zero_mat, id], [-id, zero_mat]]) % dimension
    
    # Compute modular inverse of 2 mod d (for odd prime d)
    inv2 = pow(2, -1, dimension)
    
    # Define N = Omega / 2 mod d
    N = (inv2 * Omega) % dimension
    
    # Calculate M = S.T @ N @ S - N mod d
    M = (symplectic.T @ N @ symplectic - N) % dimension
    
    # Ensure skew-symmetry: M = (M - M.T) mod d
    M = (M - M.T) % dimension
    
    return M


def symplectic_inner(u, v, d):
    n = len(u) // 2
    x1, z1 = u[:n], u[n:]
    x2, z2 = v[:n], v[n:]
    return (np.dot(x1, z2) - np.dot(x2, z1)) % d
    

class SUM(Gate):
    def __init__(self, control, target, dimension):
        images = [np.array([1, 1, 0, 0]),  # image of X0:  X0 -> X0 X1
                  np.array([0, 1, 0, 0]),  # image of X1:  X1 -> X1
                  np.array([0, 0, 1, 0]),  # image of Z0:  Z0 -> Z0
                  np.array([0, 0, -1, 1])  # image of Z1:  Z1 -> Z0^-1 Z1
                  ]
        phase_matrix = np.array([[0, 0, 0, 0],
                                [0, 0, 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]], dtype=int)
        super().__init__("SUM", [control, target], images, dimension=dimension, phase_matrix=phase_matrix)
    

class SWAP(Gate):
    def __init__(self, index1, index2, dimension):
        images = [np.array([0, 1, 0, 0]),  # image of X0:  X0 -> X0 X1
                  np.array([1, 0, 0, 0]),  # image of X1:  X1 -> X1
                  np.array([0, 0, 0, 1]),  # image of Z0:  Z0 -> Z0
                  np.array([0, 0, 1, 0])  # image of Z1:  Z1 -> Z0^-1 Z1
                  ]
        super().__init__("SWAP", [index1, index2], images, dimension=dimension)
   

class CNOT(Gate):
    def __init__(self, control: int, target: int):
        images = [np.array([1, 1, 0, 0]),  # image of X0:  X0 -> X0 X1
                  np.array([1, 0, 0, 0]),  # image of X1:  X1 -> X1
                  np.array([0, 0, 1, 0]),  # image of Z0:  Z0 -> Z0
                  np.array([0, 0, -1, 1])  # image of Z1:  Z1 -> Z0^-1 Z1
                  ]
        phase_matrix = np.array([[0, 0, 0, 0],
                                [0, 0, 0, 0],
                                [0, 0, 1, 0],
                                [0, 0, 0, 1]], dtype=int)
        super().__init__("SUM", [control, target], images, dimension=2, phase_matrix=phase_matrix)


class Hadamard(Gate):
    def __init__(self, index: int, dimension: int, inverse: bool = False):
        if inverse:
            images = [np.array([0, 1]),  # image of X:  X -> Z
                      np.array([-1, 0]),  # image of Z:  Z -> -X
                      ]
        else:
            images = [np.array([0, -1]),  # image of X:  X -> -Z
                      np.array([1, 0]),  # image of Z:  Z -> X
                      ]
        phase_matrix = np.array([[0, 0.5 if dimension == 2 else pow(2, -1, dimension)],
                                [0.5 if dimension == 2 else pow(2, -1, dimension), 0]]) % dimension
        name = "H" if not inverse else "Hdag"
        super().__init__(name, [index], images, dimension=dimension, phase_matrix=phase_matrix)


class PHASE(Gate):
    
    def __init__(self, index: int, dimension: int):
        images = [np.array([1, 1]),  # image of X:  X -> XZ
                  np.array([0, 1]),  # image of Z:  Z -> Z
            ]
        phase_matrix = np.array([[0, pow(2, -1, dimension)],
                                [pow(2, -1, dimension), 0]]) % dimension
        super().__init__("S", [index], images, dimension=dimension, phase_matrix=phase_matrix)
