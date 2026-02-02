from __future__ import annotations
from typing import Generator
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import PHASE, SUM, SWAP, Gate, Hadamard
from sympleq.applications.randomized_benchmarking.noise_model import DephasingNoise, NoiseModel, Noiseless
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.core.paulis.pauli_sum import PauliSum


class RMB:  # FIXME: after merging #106 make this s ubclass of circuit and move pretty print to circuit class.
    def __init__(self,
                 circuit: Circuit,
                 random_initial_state: bool,
                 with_random_elimination: float,
                 with_random_insertion: float,
                 noise_model: NoiseModel,
                 rng: RNGGenerator
                 ) -> None:

        self.rng = rng

        self._circuit = circuit + circuit.inv()
        self.noise_model = noise_model
        self._initial_state = RMB.initial_state(circuit.dimensions, random_initial_state, self.rng)

        # Eliminate and insert identity gates from and to the circuit to make it asymmetric.
        # This step is performed without applying errors.
        if with_random_elimination > 0.0:
            gate_to_eliminate_indices = []
            pauli = self._initial_state.copy()
            for idx, gate in enumerate(circuit.gates):
                intermediate = gate.act(pauli)
                if pauli == intermediate:
                    gate_to_eliminate_indices.append(idx)

                pauli = intermediate

            for idx in gate_to_eliminate_indices:
                if self.rng.random() <= with_random_elimination:
                    if self.rng.choice(a=[False, True]):
                        # Remove gate from mirrored circuit (right part)
                        # idx can have max value len(circuit.gates) - 1 == len(self.circuit.gates) / 2 - 1
                        self._circuit.remove_gate(len(self._circuit) - 1 - idx)
                    else:
                        # Remove gate from base circuit (left part)
                        self._circuit.remove_gate(idx)

    @classmethod
    def initial_state(cls,
                      dimensions: list[int] | np.ndarray,
                      random_phases: bool = False,
                      rng: RNGGenerator | None = None) -> PauliSum:
        if rng is None:
            rng = default_rng()

        pauli_strings = []
        for p_idx in range(n_qudits):
            pauli_string = ""
            for q_idx in range(n_qudits):
                if p_idx == q_idx:
                    pauli_string += "x0z1"
                else:
                    pauli_string += "x0z0"
            pauli_strings.append(pauli_string)

        ps = PauliSum.from_string(pauli_strings, dimensions)
        if random_phases:
            ps.set_phases(rng.choice(a=[0, 2], size=ps.n_paulis()))

        return ps

    @classmethod
    def from_circuit(cls,
                     circuit: Circuit,
                     random_initial_state: bool = True,
                     noise_model: NoiseModel = Noiseless(),
                     rng: RNGGenerator | None = None
                     ) -> RMB:
        """
        Create a random RMB object.

        Parameters
        ----------
        circuit: Circuit
            The base circuit to construct the RMB. The circuit will be mirrored, so in a sense this input
            is half the final circuit.
        random_initial_state: bool = True
            Whether the initial state should be set randomly.
        noise_model: NoiseModel
            The noise model to use to get the Kraus operators...
        rng : numpy.random.Generator | None = None
            The random number generator. Passing a value can be used to obtain deterministic randomness.

        Returns
        -------
        RMB
            A RMB object.
        """

        if rng is None:
            rng = default_rng()

        return cls(circuit, random_initial_state, False, False, noise_model, rng)

    @classmethod
    def from_random(cls,
                    dimensions: int | list[int] | np.ndarray,
                    gate_density: float = 1.0,
                    random_initial_state: bool = True,
                    with_random_elimination: float = 0.0,
                    with_random_insertion: float = 0.0,
                    noise_model: NoiseModel = Noiseless(),
                    rng: RNGGenerator | None = None
                    ) -> RMB:
        """
        Create a random RMB object.

        Parameters
        ----------
        dimensions : int | list[int] | np.ndarray
            The dimensions of the qudits. The size of dimensions determines the number of qudits.
        gate_density: float = 1.0
            The average number of gates for each qudit in the fir half of the circuit
            (the second half is the mirror of the first).
        random_initial_state: bool = True
            Whether the initial state should be set randomly.
        with_random_elimination: float = 0.0
            Whether gates acting as identity should be randomly eliminated.
            This is used to break the mirror symmetry of the circuit without
            affecting the output state (in absence of errors).
        with_random_insertion: float = 0.0
            Whether gates acting as identity should be randomly inserted.
                This is used to break the mirror symmetry of the circuit without
                affecting the output state (in absence of errors).
        noise_model: NoiseModel
            The noise model to use to get the Kraus operators...
        rng : numpy.random.Generator | None = None
            The random number generator. Passing a value can be used to obtain deterministic randomness.

        Returns
        -------
        RMB
            A RMB object.
        """

        if isinstance(dimensions, int):
            dimensions = [dimensions]

        if rng is None:
            rng = default_rng()

        n_qudits = len(dimensions)
        n_gates = int(gate_density * n_qudits)
        circuit = Circuit.from_random(n_gates, dimensions, rng=rng)

        return cls(circuit, random_initial_state,
                   with_random_elimination, with_random_insertion, noise_model, rng)

    @property
    def dimensions(self) -> np.ndarray:
        return self._circuit.dimensions

    @property
    def gates(self) -> list[Gate]:
        return self._circuit.gates

    def n_gates(self) -> int:
        return len(self._circuit.gates)

    def n_qudits(self) -> int:
        return len(self._circuit.dimensions)

    def rho_average(self, pauli_sum: PauliSum, n_runs: int = 1000) -> np.ndarray:
        def _get_output_probabilities(pauli_sum: PauliSum, n_runs: int = 1000) -> dict[PauliSum, int]:
            output_probabilities: dict[PauliSum, int] = {}
            for _ in range(n_runs):
                output = self.average_act(pauli_sum)
                if output not in output_probabilities:
                    output_probabilities[output] = 1
                else:
                    output_probabilities[output] += 1

            return output_probabilities

        output = _get_output_probabilities(pauli_sum, n_runs)
        rho: np.ndarray | None = None
        for output_pauli, count in output.items():
            if rho is None:
                rho = pauli_to_rho(output_pauli) * (count / n_runs)
            else:
                rho += pauli_to_rho(output_pauli) * (count / n_runs)

        assert rho is not None
        return rho

    def rho_exact(self, pauli_sum: PauliSum) -> np.ndarray:
        rho = pauli_to_rho(pauli_sum)

        for gate in self.gates:
            rho = self._apply_gate_to_rho_with_error(gate, rho)

        return np.around(rho, 10)

    def average_act(self, pauli_sum: PauliSum, n_runs: int = 1) -> PauliSum:
        """
        Calculate final output state, including random errors.
        Parameters
        ----------
        n_runs: int = 1
            The number times the circuit should be applied. The return value will be the average
            over all these runs.

        Returns
        -------
        PauliSum
            The resulting PauliSum after applying the noisy circuit.
        """
        if n_runs < 1:
            raise ValueError(f"Number of runs must be greater equal to 1 (got {n_runs}).")

        output = self.act(pauli_sum)
        for _ in range(n_runs - 1):
            output += self.act(pauli_sum)

        output = output / n_runs

        return output

    def _apply_gate_to_pauli_with_error(self, gate: Gate, pauli: PauliSum) -> PauliSum:
        # Get the 'correct' pauli
        pauli = gate.act(pauli)

        # Apply noise model
        pauli = self.noise_model.act(pauli, gate.qudit_indices)

        return pauli

    def _apply_gate_to_rho_with_error(self, gate: Gate, rho: np.ndarray) -> np.ndarray:
        # FIXME: after merging #106, clean this up
        if gate.name.endswith("-inv"):
            match gate.name:
                case "H-inv":
                    gatehack = Hadamard(gate.qudit_indices[0], gate.dimensions[0])
                case "S-inv":
                    gatehack = PHASE(gate.qudit_indices[0], gate.dimensions[0])
                case "SWAP-inv":
                    gatehack = SWAP(gate.qudit_indices[0], gate.qudit_indices[1], gate.dimensions[0])
                case "SUM-inv":
                    gatehack = SUM(gate.qudit_indices[0], gate.qudit_indices[1], gate.dimensions[0])
                case _:
                    print(gate)
                    raise ValueError("Unknown gate")

            unitary = gatehack.unitary(self.dimensions).transpose().conjugate()
        else:
            unitary = gate.unitary(self.dimensions)

        rho = unitary @ rho @ unitary.transpose().conjugate()

        # Apply noise using Kraus form: ρ_out = Σ_i K_i ρ K†_i
        output_rho: np.ndarray | None = None
        for K in self.noise_model.kraus_operators(self.dimensions, gate.qudit_indices):
            K_matrix = K.to_hilbert_space()
            term = K_matrix @ rho @ K_matrix.conj().T
            if output_rho is None:
                output_rho = term
            else:
                output_rho += term

        assert output_rho is not None
        return output_rho

    def act(self, pauli_sum: PauliSum) -> PauliSum:
        for gate in self.gates:
            pauli_sum = self._apply_gate_to_pauli_with_error(gate, pauli_sum)

        return pauli_sum

    def act_iter(self, pauli_sum: PauliSum) -> Generator[PauliSum, None, None]:
        for gate in self.gates:
            pauli_sum = self._apply_gate_to_pauli_with_error(gate, pauli_sum)
            yield pauli_sum

    def __str__(self) -> str:
        """
        Returns a more readable string representation of the RMB.

        Returns
        -------
        str
            A string representation of the RMB.
        """

        p_string = f"""
Initial state:
  {self._initial_state}

Circuit:
  {self._circuit}
"""
        return p_string

    def gates_layout(self, with_qudit_indices: bool = False, wrap: bool = True) -> str:
        """
        Returns a visual circuit diagram of the RMB.

        Renders the circuit as ASCII art with gates displayed as boxes
        connected by wires.

        Parameters
        ----------
        with_qudit_indices : bool, default False
            If True, display qudit indices on the left of each wire.
        wrap : bool, default True
            If True, wrap the output to fit the terminal width by splitting
            at gate boundaries.

        Returns
        -------
        str
            A string representation of the circuit diagram.
        """

        def green(s):
            return f"\033[92m{s}\033[0m"

        def red(s):
            return f"\033[91m{s}\033[0m"

        n_qudits = self.n_qudits()
        with_input = self._initial_state
        with_output = self.act(self._initial_state)
        wires = [green("=") if with_input.phases[l_idx] == with_output.phases[l_idx]
                 else red("=") for l_idx in range(n_qudits)]
        return self._circuit.gates_layout(
            with_qudit_indices=with_qudit_indices,
            with_input=with_input,
            with_output=with_output,
            wires=wires,
            wrap=wrap)


def pauli_to_rho(pauli: PauliSum) -> np.ndarray:
    _, states = pauli.ordered_eigenspectrum()
    ground_state = states[0]
    d = ground_state.size
    return np.kron(ground_state.conj(), ground_state).reshape(d, d)


if __name__ == "__main__":
    n_qudits = 4
    gate_density = 4.5
    dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits
    dimension = dimensions[0]
    circuit = Circuit(dimensions, [PHASE(0, dimension)])
    rmb = RMB.from_random(dimensions, gate_density,
                          noise_model=DephasingNoise(0.25),
                          with_random_elimination=False,
                          rng=default_rng())

    print(rmb.gates_layout(with_qudit_indices=True))
