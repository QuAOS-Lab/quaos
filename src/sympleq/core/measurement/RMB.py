from __future__ import annotations
import numpy as np
from numpy.random import Generator, default_rng

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import Gate
from sympleq.core.measurement.noise_model import DephasingNoise, DepolarizingNoise, NoiseModel, Noiseless
from sympleq.core.paulis.pauli_sum import PauliSum
from sympleq.core.paulis.utils import ground_state_TMP


class RMB:
    def __init__(self,
                 dimensions: list[int] | np.ndarray,
                 n_gates: int,
                 random_initial_state: bool,
                 with_random_elimination: float,
                 with_random_insertion: float,
                 noise_models: list[NoiseModel],
                 rng: Generator
                 ) -> None:

        self.dimensions = dimensions
        circuit = Circuit.from_random(n_gates, dimensions, rng=rng)
        self.rng = rng

        self.circuit = circuit + circuit.inv()
        self.noise_models = noise_models

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
            random_phases = self.rng.choice(a=[0, 1], size=self.initial_state.n_paulis())
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
                    if self.rng.choice(a=[0, 1]) == 0:
                        # Remove gate from mirrored circuit (right part)
                        # idx can have max value len(circuit.gates) - 1 == len(self.circuit.gates) / 2 - 1
                        self.circuit.remove_gate(len(self.circuit) - 1 - idx)
                    else:
                        # Remove gate from base circuit (left part)
                        self.circuit.remove_gate(idx)

    def get_output(self):
        # Calculate final output state, including random errors.
        return self.act(self.initial_state)

    @classmethod
    def from_random(cls,
                    dimensions: int | list[int] | np.ndarray,
                    n_gates: int | None = None,
                    random_initial_state: bool = True,
                    with_random_elimination: float = 0.0,
                    with_random_insertion: float = 0.0,
                    noise_models: NoiseModel | list[NoiseModel] = Noiseless(),
                    rng: Generator | None = None
                    ) -> RMB:
        """
        Create a random RMB object.

        Parameters
        ----------
        dimensions : int | list[int] | np.ndarray
            The dimensions of the qudits. The size of dimensions determines the number of qudits.
        n_gates: int | None = None
            The number of gates in the first half of the circuit. Since the second half is the mirror of the first,
            this corresponds to half the number of gates of the full circuit.
        random_initial_state: bool = True
            Whether the initial state should be set randomly.
        with_random_elimination: float = 0.0,
        with_random_insertion: float = 0.0,
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

        if n_gates is None:
            n_gates = rng.integers(10, 20)

        if not isinstance(noise_models, list):
            noise_models = [noise_models]

        return cls(dimensions, n_gates, random_initial_state,
                   with_random_elimination, with_random_insertion, noise_models, rng)

    @property
    def gates(self) -> list[Gate]:
        return self.circuit.gates

    @property
    def n_qudits(self) -> int:
        return self.initial_state.n_qudits()

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

    def act(self, pauli: PauliSum) -> PauliSum:
        for gate in self.circuit.gates:
            kraus_operators = [k for model in self.noise_models for k in model.kraus_operators(
                self.dimensions, gate.qudit_indices)]
            correct_pauli = gate.act(pauli)

            # Note: we do not need cross-terms, as we will pick only terms from the diagonal.
            options = [k[1] * correct_pauli * k[1].H() for k in kraus_operators]
            probability_distribution = [np.abs(k[0])**2 / len(self.noise_models) for k in kraus_operators]

            # Pick one possible output with weight given by the Kraus operator weight.
            pauli = self.rng.choice(np.asarray(options), p=probability_distribution)

        return pauli


def luca_check():
    rmb = RMB.from_random(dimensions, n_gates,
                          noise_models=DephasingNoise(0.0))
    output = rmb.get_output()
    output.phase_to_weight()
    print(output)
    _, gs = ground_state_TMP(output)
    rho = np.kron(gs.conj(), gs)

    N = 10
    for _ in range(N - 1):
        rmb = RMB.from_random(dimensions, n_gates,
                              noise_models=DephasingNoise(0.0))
        output = rmb.get_output()
        _, gs = ground_state_TMP(output)
        n_rho = np.kron(gs.conj(), gs)
        rho += n_rho

    rho = rho / N
    print(np.around(rho, decimals=4))

    rmb_exact = RMB.from_random(dimensions, n_gates,
                                noise_models=Noiseless(), rng=default_rng(0))

    ouput_exact = rmb_exact.get_output()
    _, gs = ground_state_TMP(ouput_exact)
    rho_exact = np.kron(gs.conj(), gs)

    print()
    print(np.around(rho_exact, decimals=4))

    print(rmb_exact.fancy_str())

    # Run N samples with fixed RMB (initialize rng).
    # transform output to hilbert space
    # average over all samples
    # start with initial hilbert space
    # apply channel action
    # compare


if __name__ == "__main__":
    n_gates = 8
    dimensions = [2] * 4

    rmb = RMB.from_random(dimensions, n_gates,
                          noise_models=[DephasingNoise(0.05), DepolarizingNoise(0.05)],
                          rng=default_rng(10))
    print(rmb.fancy_str())
