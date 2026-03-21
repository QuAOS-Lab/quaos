from __future__ import annotations
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import Gate
from sympleq.core.paulis._typing import HilbertOperator
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.core.paulis.pauli_sum import PauliSum
from sympleq.core.noise.noise_model import \
    CompositeNoise, DephasingNoise, DepolarizingNoise, NoiseModel, Noiseless
from sympleq.core.bayesian_estimation import BayesianEstimator


class RMB:
    def __init__(self,
                 circuit: Circuit,
                 scrambler: Circuit,
                 with_random_elimination: float,
                 with_random_insertion: float,
                 rng: RNGGenerator,
                 noise_model: NoiseModel | None = None,
                 ) -> None:

        self.rng = rng

        self._scrambler = scrambler
        self._scrambler_inv = scrambler.inverse()

        self._circuit = circuit
        self._circuit_inv = circuit.inverse()

        self._circuit.set_noise(noise_model)
        self._circuit_inv.set_noise(noise_model)

        self._initial_state = RMB.initial_state(circuit.dimensions)

        # Eliminate and insert identity gates from and to the circuit to make it asymmetric.
        # This step is performed without applying errors.
        if with_random_elimination > 0.0:
            gate_to_eliminate_indices = []
            pauli = self._initial_state.copy()
            for idx, (gate, idxs) in enumerate(zip(circuit.gates, circuit.qudit_indices)):
                intermediate = gate.act(pauli, idxs)
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
    def initial_state(cls, dimensions: list[int] | np.ndarray) -> PauliSum:
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

        ps = PauliSum.from_string(pauli_strings, dimensions)
        return ps

    @classmethod
    def from_circuit(cls,
                     circuit: Circuit,
                     rng: RNGGenerator | None = None
                     ) -> RMB:
        """
        Create a random RMB object.

        Parameters
        ----------
        circuit: Circuit
            The base circuit to construct the RMB. The circuit will be mirrored, so in a sense this input
            is half the final circuit.
        rng : numpy.random.Generator | None = None
            The random number generator. Passing a value can be used to obtain deterministic randomness.

        Returns
        -------
        RMB
            A RMB object.
        """

        if rng is None:
            rng = default_rng()

        return cls(circuit, Circuit.empty(circuit.dimensions), False, False, rng)

    @classmethod
    def from_random(cls,
                    dimensions: int | list[int] | np.ndarray,
                    gate_density: float = 1.0,
                    with_random_elimination: float = 0.0,
                    with_random_insertion: float = 0.0,
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
        scrambler = Circuit.from_random(10, dimensions, two_qudit_gate_ratio=0.0, rng=rng)
        circuit = Circuit.from_random(n_gates, dimensions, rng=rng)

        return cls(circuit, Circuit.empty(circuit.dimensions), with_random_elimination, with_random_insertion, rng)\
            .with_scrambler(scrambler)

    def with_noise(self, noise_model: NoiseModel) -> RMB:
        """
        Attach a noise model and return self for chaining.

        Parameters
        ----------
        noise_model : NoiseModel
            The noise model to apply after each gate.

        Returns
        -------
        RMB
            This RMB instance (for method chaining).
        """
        self._circuit.set_noise(noise_model)
        self._circuit_inv.set_noise(noise_model)
        return self

    def with_scrambler(self, scrambler: Circuit) -> RMB:
        """
        Attach a scrambler to randomize the initial state and return self for chaining.

        Parameters
        ----------
        scrambler : Circuit
            The scrambler to prepend to the circuit.

        Returns
        -------
        RMB
            This RMB instance (for method chaining).
        """
        if not np.array_equal(scrambler.dimensions, self._circuit.dimensions):
            raise ValueError("Scrambler and circuit must have the same dimensions.")

        if not all([g.n_qudits == 1 for g in scrambler.gates]):
            raise ValueError("Scrambler gates should be all 1-qudit gates.")

        self._scrambler = scrambler
        self._scrambler_inv = scrambler.inverse()
        return self

    def with_random_scrambler(self, n_gates: int = 10) -> RMB:
        """
        Attach a random scrambler to randomize the initial state and return self for chaining.

        Parameters
        ----------
        n_gates : int | None
            The depth of the random scrambler.

        Returns
        -------
        RMB
            This RMB instance (for method chaining).
        """
        self._scrambler = Circuit.from_random(n_gates, self.dimensions, two_qudit_gate_ratio=0.0, rng=self.rng)
        self._scrambler_inv = self._scrambler.inverse()
        return self

    @property
    def dimensions(self) -> np.ndarray:
        return self._circuit.dimensions

    @property
    def lcm(self) -> int:
        return self._circuit.lcm

    @property
    def gates(self) -> list[Gate]:
        return self._circuit.gates

    @property
    def qudit_indices(self) -> list[tuple[int, ...]]:
        return self._circuit.qudit_indices

    def n_gates(self) -> int:
        return len(self._circuit.gates)

    def n_qudits(self) -> int:
        return len(self._circuit.dimensions)

    def run(self) -> PauliSum:
        ps = self._initial_state
        return self.circuit().act(ps)

    def run_in_hilbert_space(self) -> HilbertOperator:
        rho = self._initial_state.stabilizer_to_hilbert_space()
        return self.circuit().act_in_hilbert_space(rho)

    def circuit(self) -> Circuit:
        return self._scrambler + self._circuit + self._circuit_inv + self._scrambler_inv

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
        with_output = self.run()
        wires = [green("=") if with_input.phases[l_idx] == with_output.phases[l_idx]
                 else red("=") for l_idx in range(n_qudits)]

        return self.circuit().gates_layout(
            with_qudit_indices=with_qudit_indices,
            with_input=with_input,
            with_output=with_output,
            wires=wires,
            wrap=wrap)


def _fidelity(rho: np.matrix, sigma: np.matrix) -> float:
    from scipy.linalg import sqrtm

    sq_rho = sqrtm(rho)
    tmp_matrix = sq_rho @ sigma @ sq_rho
    return np.real(np.trace(sqrtm(tmp_matrix)) ** 2)


def fidelity(estimator: BayesianEstimator, target: HilbertOperator) -> tuple[float, float]:
    sigma = np.sum([estimator.probability(res) * res.stabilizer_to_hilbert_space() for res in results])
    rho = target.todense()
    fidelity = _fidelity(rho, sigma)

    rng = default_rng()

    def _sample_fidelity() -> float:
        sampled_probabilities = {}
        cum_probability = 0
        for res in estimator.results():
            p = estimator.probability(res)
            std = np.sqrt(estimator.variance(res))
            s = max(0.0, rng.normal(p, std))
            cum_probability += s
            sampled_probabilities[res] = s

        sampled_sigma = np.sum([sampled_probabilities[res] / cum_probability *
                                res.stabilizer_to_hilbert_space() for res in results])
        return _fidelity(rho, sampled_sigma)

    N = 100
    sampling = np.empty(N, dtype=float)
    for i in range(N):
        sampling[i] = _sample_fidelity()

    avg = np.mean(sampling)
    std = float(np.std(sampling))

    assert abs(avg - fidelity) <= std, "Inconsistend sampling"
    return fidelity, std


if __name__ == "__main__":
    n_qudits = 4
    gate_density = 6
    dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits
    threshold = 0.05 * 1e-3
    rng = default_rng(11)

    noise_model = CompositeNoise.from_noise_models([DephasingNoise(0.05), DepolarizingNoise(0.01)], rng=rng)
    # noise_model = DepolarizingNoise(1.0, rng)
    # noise_model = DephasingNoise(0.5, rng)
    # noise_model = Noiseless()
    rmb = RMB.from_random(dimensions, gate_density,
                          with_random_elimination=True,
                          rng=rng
                          ).with_noise(noise_model).with_random_scrambler()

    print(rmb.gates_layout(with_qudit_indices=True))

    output_rho = rmb.run_in_hilbert_space()
    output_rho.eliminate_zeros()

    def _callable() -> PauliSum:
        return rmb.run()

    import time
    now = time.time()
    estimator = BayesianEstimator(threshold, min_runs=1000)

    n_printed = 0

    for _ in estimator.run_iter(_callable):
        if n_printed > 0:
            print(f"\033[{n_printed}A", end="")

        n_printed = 1
        print(f"Threshold={threshold} - {time.time() - now:.2f}s")

        results: list[PauliSum] = estimator.results()
        for idx, res in enumerate(results):
            p = estimator.probability(res)
            std = np.sqrt(estimator.variance(res))
            print(f"\033[K{res.phases}: p={p:.5f} ± {std:.5f}")
        n_printed += len(results)

        n_runs = estimator.num_runs()
        print(f"n_runs={n_runs}")
        n_printed += 1

        if n_runs % 1000 == 1:
            average_rho = sum([estimator.probability(res) * res.stabilizer_to_hilbert_space() for res in results])
            m = np.max(np.abs(average_rho - output_rho))
            print(f"Max error={m:.5f}")
            n_printed += 1

    print()
    fid, err = fidelity(estimator, output_rho)
    print(f"Ballistic accuracy={fid:.5f} ± {err:.5f}")

    rho = rmb._initial_state.stabilizer_to_hilbert_space()
    fid, err = fidelity(estimator, rho)
    print(f"Fidelity={fid:.5f} ± {err:.5f}")
