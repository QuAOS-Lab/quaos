

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
import numpy as np


class RMB:
    def __init__(self, n_gates: int | None = None, n_qudits: int | None = None) -> None:
        if n_gates is None:
            n_gates = np.random.randint(2, 10)

        if n_qudits is None:
            n_qudits = np.random.randint(2, 10)

        circuit = Circuit.from_random(
            n_gates, dimensions=[DEFAULT_QUDIT_DIMENSION for _ in range(n_qudits)])

        self.circuit = circuit + circuit.inv()