from __future__ import annotations
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import PHASE, Gate
import numpy as np
from numpy.random import Generator, default_rng
from sympleq.core.paulis.pauli_sum import PauliSum


class RMB:
    def __init__(self, dimensions: list[int] | np.ndarray, n_gates: int, rng: Generator | None) -> None:
        circuit = Circuit.from_random(n_gates, dimensions)

        if rng is None:
            rng = default_rng()
        self.rng = rng

        self.circuit = circuit + circuit.inv()
        self.error_profile: dict[Gate, float] = {}
        for gate in self.circuit.gates:
            if gate not in self.error_profile:
                self.error_profile[gate] = self.rng.random() * 0.05

    @classmethod
    def from_random(cls, dimensions: int | list[int] | np.ndarray, rng: Generator | None = None) -> RMB:
        """
        Create a random RMB object.

        Parameters
        ----------
        dimensions : int | list[int] | np.ndarray
            The dimensions of the qudits. The size of dimensions determines the number of qudits.
        rng : numpy.random.Generator | None = None
            The random number generator. Passing a value can be used to obtain deterministic randomness.

        Returns
        -------
        RMB
            A RMB object.
        """

        if isinstance(dimensions, int):
            dimensions = [dimensions]

        n_gates = np.random.randint(10, 20)

        return cls(dimensions, n_gates, rng)

    @property
    def gates(self) -> list[Gate]:
        return self.circuit.gates

    def act(self, pauli: PauliSum) -> PauliSum:
        for gate in self.circuit.gates:
            error = self.error_profile[gate]
            if self.rng.random() <= error:
                # FIXME: should pass own rng to get deterministic results
                random_gate = Gate.from_random(gate.n_qudits, gate.dimensions[0])
                pauli = random_gate.act(pauli)
            else:
                pauli = gate.act(pauli)

        return pauli


if __name__ == "__main__":
    n_gates = 20
    dimensions = [2] * 10
    rmb = RMB.from_random(dimensions)
    ps = PauliSum.from_random(10, dimensions)

    res = rmb.act(ps)

    print(ps)
    print(res)
