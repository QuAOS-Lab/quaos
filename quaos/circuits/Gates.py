import numpy as np
from quaos.paulis import PauliString, PauliSum, Pauli
from typing import Callable


class Gate:

    def __init__(self, name: str,
                 qudit_indices: list[int],
                 images: list[np.ndarray],
                 dimension: int,
                 phase_function: Callable):
        
        self.dimension = dimension
        self.name = name
        self.qudit_indices = qudit_indices
        self.images = images
        self.n_qudits = len(qudit_indices)
        self.symplectic = np.stack([v % dimension for v in images]).T
        self.phase_function = phase_function

    def _act_on_pauli_string(self, P: PauliString) -> tuple[PauliString, int]:
        local_symplectic = np.concatenate([P.x_exp[self.qudit_indices], P.z_exp[self.qudit_indices]])
        acquired_phase = self.phase_function(P)

        local_symplectic = (self.symplectic @ local_symplectic) % self.dimension
        P = P._replace_symplectic(local_symplectic, self.qudit_indices)
        return P, acquired_phase
    
    def _act_on_pauli_sum(self, P: PauliSum):
        pauli_strings, phases = zip(*[self._act_on_pauli_string(p) for p in P.pauli_strings])
        P = PauliSum(pauli_strings, P.weights, np.asarray(phases), P.dimensions, False)
        return P

    def act(self, P: Pauli | PauliString | PauliSum):
        if isinstance(P, Pauli):
            P = PauliString(P, dimensions=[P.dimension])
        if isinstance(P, PauliString):
            return self._act_on_pauli_string(P)[0]
        elif isinstance(P, PauliSum):
            return self._act_on_pauli_sum(P)
        else:
            raise TypeError(f"Unsupported type {type(P)} for Gate.act. Expected Pauli, PauliString or PauliSum.")


class SUM(Gate):
    def __init__(self, control, target, dimension):
        images = [np.array([1, 1, 0, 0]),  # image of X0:  X0 -> X0 X1
                  np.array([0, 1, 0, 0]),  # image of X1:  X1 -> X1
                  np.array([0, 0, 1, 0]),  # image of Z0:  Z0 -> Z0
                  np.array([0, 0, -1, 1])  # image of Z1:  Z1 -> Z0^-1 Z1
                  ]
        
        @staticmethod
        def phase_function(P: PauliString) -> int:
            """
            Returns the phase acquired by the PauliString P when acted upon by this SUM gate.
            """
            return P.x_exp[control] * P.z_exp[target] % P.lcm
        
        super().__init__("SUM", [control, target], images, dimension=dimension, phase_function=phase_function)
    

class SWAP(Gate):
    def __init__(self, index1, index2, dimension):
        images = [np.array([0, 1, 0, 0]),  # image of X0:  X0 -> X0 X1
                  np.array([1, 0, 0, 0]),  # image of X1:  X1 -> X1
                  np.array([0, 0, 0, 1]),  # image of Z0:  Z0 -> Z0
                  np.array([0, 0, 1, 0])  # image of Z1:  Z1 -> Z0^-1 Z1
                  ]
        
        super().__init__("SWAP", [index1, index2], images, dimension=dimension, phase_function=lambda P: 0)
   

class CNOT(Gate):
    def __init__(self, control, target):
        images = [np.array([1, 1, 0, 0]),  # image of X0:  X0 -> X0 X1
                  np.array([0, 1, 0, 0]),  # image of X1:  X1 -> X1
                  np.array([0, 0, 1, 0]),  # image of Z0:  Z0 -> Z0
                  np.array([0, 0, -1, 1])  # image of Z1:  Z1 -> Z0^-1 Z1
                  ]
        
        @staticmethod
        def phase_function(P: PauliString) -> int:
            """
            Returns the phase acquired by the PauliString P when acted upon by this SUM gate.
            """
            return P.x_exp[control] * P.z_exp[target] % 2
        
        super().__init__("SUM", [control, target], images, dimension=2, phase_function=phase_function)


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
            
        @staticmethod
        def phase_function(P: PauliString) -> int:
            """
            Returns the phase acquired by the PauliString P when acted upon by this Hadamard gate.
            """
            return P.x_exp[index] * P.z_exp[index] % P.lcm

        name = "H" if not inverse else "Hdag"
        super().__init__(name, [index], images, dimension=dimension, phase_function=phase_function)


class PHASE(Gate):
    
    def __init__(self, index: int, dimension: int):
        images = [np.array([1, 1]),  # image of X:  X -> XZ
                  np.array([0, 1]),  # image of Z:  Z -> Z
                  ]

        @staticmethod
        def phase_function(P: PauliString) -> int:
            """
            Returns the phase acquired by the PauliString P when acted upon by this Hadamard gate.
            """
            r = P.x_exp[index]
            return r * (r - 1) // 2

        super().__init__("S", [index], images, dimension=dimension, phase_function=phase_function)
