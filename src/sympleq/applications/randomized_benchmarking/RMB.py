from __future__ import annotations
from typing import Generator
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import PHASE, SUM, SWAP, Gate, Hadamard
from .noise_model import DephasingNoise, DepolarizingNoise, NoiseModel, Noiseless
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.core.paulis.pauli_sum import PauliSum


class RMB:  # FIXME: after merging #106 make this s ubclass of circuit and move pretty print to circuit class.
    def __init__(self,
                 circuit: Circuit,
                 random_initial_state: bool,
                 with_random_elimination: float,
                 with_random_insertion: float,
                 noise_models: list[NoiseModel],
                 rng: RNGGenerator
                 ) -> None:

        self.rng = rng

        self.circuit = circuit + circuit.inv()
        self.noise_models = noise_models

        dimensions = circuit.dimensions
        n_qudits = len(dimensions)

        pauli_strings = []
        for p_idx in range(n_qudits):
            pauli_string = ""
            for q_idx in range(n_qudits):
                if p_idx == q_idx:
                    pauli_string += "x0z1"
                else:
                    pauli_string += "x0z0"
            pauli_strings.append(pauli_string)

        self.initial_state = PauliSum.from_string(pauli_strings, dimensions)

        if random_initial_state:
            random_phases = self.rng.choice(a=[0, 2], size=self.initial_state.n_paulis())
            self.initial_state.set_phases(random_phases)

        # Eliminate and insert identity gates from and to the circuit to make it asymmetric.
        # This step is performed without applying errors.
        if with_random_elimination > 0.0:
            pauli = self.initial_state
            gate_to_eliminate_indices = []
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
                        self.circuit.remove_gate(len(self.circuit) - 1 - idx)
                    else:
                        # Remove gate from base circuit (left part)
                        self.circuit.remove_gate(idx)

        # Initialize noise models probabilities.
        # In principle, we don't know if there are gates with n_qudits larger than 2,
        # good enough for now.
        self.noise_models_kraus_probabilities = {}

        for gate_n_qudits in (1, 2, 3):
            probs = np.concatenate([
                model.kraus_probabilities(gate_n_qudits)
                for model in self.noise_models
            ])
            # Normalize probabilities when combining multiple noise models
            probs /= probs.sum()
            # Given probabilities [p0, p1, p2, p3], cumsum gives [p0, p0+p1, p0+p1+p2, 1.0].
            # This allows O(log n) sampling via searchsorted with a uniform random number
            # in _apply_gate_to_pauli_with_error.
            self.noise_models_kraus_probabilities[gate_n_qudits] = np.cumsum(probs)

    @classmethod
    def from_random(cls,
                    dimensions: int | list[int] | np.ndarray,
                    gate_density: float = 1.0,
                    random_initial_state: bool = True,
                    with_random_elimination: float = 0.0,
                    with_random_insertion: float = 0.0,
                    noise_models: NoiseModel | list[NoiseModel] = Noiseless(),
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

        if not isinstance(noise_models, list):
            noise_models = [noise_models]

        n_qudits = len(dimensions)
        n_gates = int(gate_density * n_qudits)
        circuit = Circuit.from_random(n_gates, dimensions, rng=rng)

        return cls(circuit, random_initial_state,
                   with_random_elimination, with_random_insertion, noise_models, rng)

    @property
    def dimensions(self) -> np.ndarray:
        return self.circuit.dimensions

    @property
    def gates(self) -> list[Gate]:
        return self.circuit.gates

    @property
    def n_gates(self) -> int:
        return len(self.circuit.gates)

    @property
    def n_qudits(self) -> int:
        return self.initial_state.n_qudits()

    def average_error(self, n_runs: int = 10) -> float:
        error_runs = 0
        for _ in range(n_runs):
            if self.initial_state != self.get_output():
                error_runs += 1
        return error_runs / n_runs

    def get_output_probabilities(self, n_runs: int = 1000) -> dict[PauliSum, int]:
        output_probabilities: dict[PauliSum, int] = {}
        for _ in range(n_runs):
            output = self.get_output()
            if output not in output_probabilities:
                output_probabilities[output] = 1
            else:
                output_probabilities[output] += 1

        return output_probabilities

    def rho_average(self, n_runs: int = 1000) -> np.ndarray:
        output = self.get_output_probabilities(n_runs)
        rho: np.ndarray | None = None
        for output_pauli, count in output.items():
            if rho is None:
                rho = pauli_to_rho(output_pauli) * (count / n_runs)
            else:
                rho += pauli_to_rho(output_pauli) * (count / n_runs)

        assert rho is not None
        return rho

    def rho_exact(self) -> np.ndarray:
        rho = pauli_to_rho(self.initial_state)

        for gate in self.gates:
            rho = self._apply_gate_to_rho_with_error(gate, rho)

        return np.around(rho, 10)

    def get_output(self, n_runs: int = 1) -> PauliSum:
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

        output = self.act(self.initial_state)
        for _ in range(n_runs - 1):
            output += self.act(self.initial_state)

        output = output / n_runs

        return output

    def _apply_gate_to_pauli_with_error(self, gate: Gate, pauli: PauliSum) -> PauliSum:
        correct_pauli = gate.act(pauli)

        # Get probabilities to select one possible quantum trajectory
        probs = self.noise_models_kraus_probabilities[gate.n_qudits]
        idx = np.searchsorted(probs, self.rng.random())

        # locate the corresponding model and apply Kraus operator
        i = int(idx)
        for model in self.noise_models:
            nk = model.n_kraus_operators() ** gate.n_qudits
            if i < nk:
                # FIXME: way faster just to act with the kraus operator on the PauliSum directly
                # without even initializing the Kraus PauliSum.
                # pauli = model.apply_kraus_operator(correct_pauli, gate.qudit_indices, i)
                k = model.kraus_operators(self.dimensions, gate.qudit_indices)[i]
                prob = model.kraus_probabilities(gate.n_qudits)[i]
                pauli = (k * correct_pauli * k.H()) / prob

                break
            i -= nk

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

        correct_rho = unitary @ rho @ unitary.transpose().conjugate()

        # Apply noise using Kraus form: ρ_out = Σ_i K_i ρ K†_i
        output_rho: np.ndarray | None = None
        for model in self.noise_models:
            for K in model.kraus_operators(self.dimensions, gate.qudit_indices):
                K_matrix = K.to_hilbert_space()
                term = K_matrix @ correct_rho @ K_matrix.conj().T
                if output_rho is None:
                    output_rho = term
                else:
                    output_rho += term

        assert output_rho is not None

        # Normalize when combining multiple noise models
        output_rho /= len(self.noise_models)

        return output_rho

    def act(self, pauli: PauliSum) -> PauliSum:
        for gate in self.gates:
            pauli = self._apply_gate_to_pauli_with_error(gate, pauli)

        return pauli

    def act_iter(self, pauli: PauliSum) -> Generator[PauliSum, None, None]:
        for gate in self.gates:
            pauli = self._apply_gate_to_pauli_with_error(gate, pauli)
            yield pauli

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
  {self.initial_state}

Circuit:
  {self.circuit}
"""
        return p_string

    def fancy_str(self) -> str:
        """
        Returns a fancy string representation of the RMB.

        Returns
        -------
        str
            A string representation of the RMB.
        """

        def green(s):
            return f"\033[92m{s}\033[0m"

        def red(s):
            return f"\033[91m{s}\033[0m"

        def gate_name(gate: Gate) -> str:
            return gate.name.replace("-inv", "*")[:gate_name_len].center(gate_name_len)

        output_state = self.get_output()

        lines: list[str] = ["" for _ in range(3 * self.n_qudits)]
        gate_num: list[int] = [0 for _ in range(self.n_qudits)]
        wires = [green("=") if self.initial_state.phases[l_idx] == output_state.phases[l_idx]
                 else red("=") for l_idx in range(self.n_qudits)]

        gate_name_len = 5
        gate_len = gate_name_len + 4

        # Put initial state phase on the left
        for l_idx in range(self.n_qudits):
            lines[3 * l_idx + 0] = " " * 4
            lines[3 * l_idx + 1] = f"{self.initial_state.phases[l_idx]} " + wires[l_idx] * 2
            lines[3 * l_idx + 2] = " " * 4

        for gate in self.gates:
            if gate.n_qudits == 1:
                l_idx = gate.qudit_indices[0]

                gate_num[l_idx] += 1

                lines[3 * l_idx + 0] += " ┌" + "─" * gate_name_len + "┐ "
                lines[3 * l_idx + 1] += wires[l_idx] + "│" + f"{gate_name(gate)}" + "│" + wires[l_idx]
                lines[3 * l_idx + 2] += " └" + "─" * gate_name_len + "┘ "
            # 2-qudit gate
            else:
                # Get max line length of affected qudits
                max_num_gate_affected_qudits = max(
                    [gate_num[idx] for idx in gate.qudit_indices])

                for l_idx in gate.qudit_indices:
                    while gate_num[l_idx] < max_num_gate_affected_qudits:
                        gate_num[l_idx] += 1
                        lines[3 * l_idx + 0] += " " * gate_len
                        lines[3 * l_idx + 1] += wires[l_idx] * gate_len
                        lines[3 * l_idx + 2] += " " * gate_len

                    gate_num[l_idx] += 1

                    is_top_qudit = l_idx == min(gate.qudit_indices)
                    is_btm_qudit = l_idx == max(gate.qudit_indices)

                    if is_top_qudit:
                        lines[3 * l_idx + 0] += " ┌" + "─" * gate_name_len + "┐ "
                    else:
                        lines[3 * l_idx + 0] += " ┌" + "─" * (gate_name_len // 2) + \
                            "┴" + "─" * (gate_name_len // 2) + "┐ "

                    lines[3 * l_idx + 1] += wires[l_idx] + "│" + f"{gate_name(gate)}" + "│" + wires[l_idx]
                    if is_btm_qudit:
                        lines[3 * l_idx + 2] += " └" + "─" * gate_name_len + "┘ "
                    else:
                        lines[3 * l_idx + 2] += " └" + "─" * (gate_name_len // 2) + \
                            "┬" + "─" * (gate_name_len // 2) + "┘ "

        max_num_gate = max(gate_num)
        for l_idx in range(self.n_qudits):
            while gate_num[l_idx] < max_num_gate:
                gate_num[l_idx] += 1
                lines[3 * l_idx + 0] += " " * gate_len
                lines[3 * l_idx + 1] += wires[l_idx] * gate_len
                lines[3 * l_idx + 2] += " " * gate_len

        # Put final state phase on the right
        for l_idx in range(self.n_qudits):
            lines[3 * l_idx + 0] += " " * 4
            lines[3 * l_idx + 1] += wires[l_idx] * 2 + " " + f"{output_state.phases[l_idx]}"
            lines[3 * l_idx + 2] += " " * 4

        return "\n".join(lines)


def pauli_to_rho(pauli: PauliSum) -> np.ndarray:
    _, states = pauli.ordered_eigenspectrum()
    ground_state = states[0]
    d = ground_state.size
    return np.kron(ground_state.conj(), ground_state).reshape(d, d)


if __name__ == "__main__":
    n_qudits = 3
    gate_density = 2.5
    dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits
    dimension = dimensions[0]
    circuit = Circuit(dimensions, [PHASE(0, dimension)])
    rmb = RMB.from_random(dimensions, gate_density,
                          noise_models=[DephasingNoise(0.25), DepolarizingNoise(0.15)],
                          with_random_elimination=False,
                          rng=default_rng())

    print(rmb.fancy_str())
