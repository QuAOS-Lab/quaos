from __future__ import annotations
from numpy.random import Generator as RNGGenerator
from pytket.circuit import Circuit as PytketCircuit
import warnings

from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.backends.utils import data_from_pytket_circuit_results
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.core.noise.noise_model import GenericNoise
from sympleq.integrations.quantinuum.utils import (
    NATIVE_GATES_SET,
    fetch_recent_execute_jobs,
    to_pytket_circuit,
    pytket_bare_simulation_cost,
    BASE_SIMULATION_COST
)


MAX_COST_PER_RUN: float = 50.0


class QuantinuumBackend(RMBBackend):
    """
    Quantinuum hardware/emulator backend.

    Builds a random circuit from the config, submits it to the
    Quantinuum device identified by ``device_name`` for ``n_shots``
    shots, and returns whether the most-common measurement outcome is
    the all-zeros bitstring.

    Parameters
    ----------
    device_name : str
        Quantinuum device identifier (e.g. ``"H2-1LE"``).
    n_shots : int
        Number of shots to request from the backend.
    """

    type = "quantinuum"

    def __init__(self, device_name: str = "H2-2E", n_shots: int = 1,
                 project_name: str = "Benchmark", batch_size: int = 11, max_cost_per_run: float = 10.0) -> None:
        self.device_name = device_name
        self.n_shots = n_shots
        self.project_name = project_name
        self.batch_size = batch_size
        if max_cost_per_run > MAX_COST_PER_RUN:
            raise ValueError(
                "Max cost per run should not exceed MAX_COST_PER_RUN ({max_cost_per_run} <= {MAX_COST_PER_RUN}).")
        self.max_cost_per_run = max_cost_per_run

    def fidelity_estimation(self, config: RMBConfig, rng: RNGGenerator) -> list[bool]:
        from sympleq.integrations.quantinuum.workflow import run_circuits_on_device
        from sympleq.integrations.quantinuum.stitching import circuit_stitching, destitch_results, \
            estimate_qasm_program_size, MAX_QASM_PROGRAM_SIZE
        from sympleq.integrations.quantinuum.utils import pytket_simulation_cost

        def generate_compatible_circuit() -> PytketCircuit:
            while True:
                circuit = to_pytket_circuit(config.random_circuit(rng=rng))
                if circuit.n_gates == 0:
                    continue
                if circuit.n_1qb_gates() != config.n_1qb_gates:
                    warnings.warn(
                        Warning(f"Generated circuit \
has mismatching number of 1qb gates ({circuit.n_1qb_gates()} vs {config.n_1qb_gates})."))
                    continue
                if circuit.n_2qb_gates() != config.n_2qb_gates:
                    Warning(f"Generated circuit \
has mismatching number of 1qb gates ({circuit.n_2qb_gates()} vs {config.n_2qb_gates}).")
                    continue
                return circuit

        # We stitch the list circuits into a (smaller) list of stitched circuits.
        stitched_circuits: list[PytketCircuit] = []
        # We keep track of how many underlying circuits were stitched together in the corresponing
        # circuit in stitched_circuits
        stitch_sizes: list[int] = []

        while len(stitched_circuits) < self.batch_size:
            circuits_to_stich = []
            program_sizes = []
            simulation_costs = []
            while True:
                # Approximate extra cost as sum of current stitched circuit and new circuit
                append_circuit = generate_compatible_circuit()
                append_program_size = estimate_qasm_program_size(append_circuit)
                append_simulation_cost = pytket_bare_simulation_cost(append_circuit)
                num_resets = len(circuits_to_stich) * config.n_qubits
                running_cost = BASE_SIMULATION_COST + sum(simulation_costs) + num_resets / 5000
                if sum(program_sizes) + append_program_size > MAX_QASM_PROGRAM_SIZE:
                    break

                if running_cost + append_simulation_cost > self.max_cost_per_run:
                    break

                circuits_to_stich.append(append_circuit)
                program_sizes.append(append_program_size)
                simulation_costs.append(append_simulation_cost)

                print(f"Running cost: {running_cost}")

            stitched_circuit = circuit_stitching(circuits_to_stich)
            stitched_circuits.append(stitched_circuit)
            stitch_sizes.append(len(circuits_to_stich))
            # sub_circuit_configs.append(RMBConfig.from_pytket_circuit(append_circuit))
            append_config = RMBConfig.from_pytket_circuit(append_circuit)
            assert append_config.n_1qb_gates == config.n_1qb_gates
            assert append_config.n_2qb_gates == config.n_2qb_gates
            assert append_config.n_qubits == config.n_qubits

        print(f"Generated {len(stitched_circuits)} stitched circuits with sizes {stitch_sizes}.")

        stitched_costs = [pytket_simulation_cost(c) for c in stitched_circuits]
        native_costs = [pytket_simulation_cost(c) + 5 * (size - 1) for c, size in zip(stitched_circuits, stitch_sizes)]
        print(f"Stitched costs: {stitched_costs}")
        print(f"Native costs:   {native_costs}")
        print(f"You saved {sum(stitched_costs) / sum(native_costs) * 100.0}%.")

        results = run_circuits_on_device(
            stitched_circuits, self.n_shots, self.device_name, self.project_name, verbose=True)

        fidelities = []
        c_idx = 0
        for (res, circuit) in zip(results, stitched_circuits):
            unstitched_results = destitch_results(res, circuit.c_registers)
            for u_res in unstitched_results:
                distribution = u_res.get_empirical_distribution()
                counts = distribution.as_counter()
                if not counts:
                    continue

                total = distribution.total
                for state, count in counts.most_common():
                    print(f"{state}: {count / total:.4f} ({count}/{total})")

                outcome, count = counts.most_common()[0]
                assert count <= self.n_shots
                # Compare to initial state
                res = all(bit == 0 for bit in outcome)
                c_idx += 1
                fidelities.append(res)

        return fidelities

    @classmethod
    def default_config(cls) -> RMBConfig:
        return RMBConfig.default()\
            .with_n_1qb_gates(20 * 6)\
            .with_n_2qb_gates(5 * 6)\
            .with_random_elimination(0.1)\
            .with_n_qubits(6)\
            .with_gates_set(tuple(NATIVE_GATES_SET))

    def default_estimator(self) -> BayesianEstimator:
        """Return a fresh :class:`BayesianEstimator` with ``threshold=0.1`` and ``min_runs=self.batch_size``."""
        return BayesianEstimator(threshold=10**(-1), min_runs=self.batch_size)

    @classmethod
    def default_sympleq_backend(cls) -> SympleqBackend:
        noise_model = GenericNoise.from_paulis([0.000025, 0.000025, 0.000025])
        two_qubit_noise_model = GenericNoise.from_paulis([0.00079, 0.00079, 0.00079])
        return SympleqBackend(
            noise_model=noise_model,
            two_qubit_noise_model=two_qubit_noise_model,
        )

    def populate_from_recent_jobs(self, n: int) -> RMBData:
        """Build :type:`RMBData` from the last ``n`` execute jobs in ``self.project_name``.

        Each fetched circuit's ``n_qubits``, non-measurement gate count, and
        two-qudit gate ratio are used to synthesise an :class:`RMBConfig`;
        remaining fields fall back to defaults. The outcome recorded is
        ``True`` iff the most-common measured bitstring is all-zeros.

        Parameters
        ----------
        n : int
            Maximum number of execute jobs to fetch.

        Returns
        -------
        RMBData
            Mapping of inferred configs to their estimators.
        """
        jobs = fetch_recent_execute_jobs(self.project_name, n, device_name=self.device_name)
        pairs = []
        for circuit, result in jobs:
            counts = result.get_empirical_distribution().as_counter()
            if not counts:
                continue
            outcome, _ = counts.most_common()[0]
            pairs.append(all(bit == 0 for bit in outcome))

        return data_from_pytket_circuit_results(pairs, self.default_estimator)

    def to_dict(self) -> dict:
        return {
            "type": self.type,
            "device_name": self.device_name,
            "n_shots": self.n_shots,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> QuantinuumBackend:
        return cls(device_name=payload["device_name"], n_shots=payload["n_shots"])


# noise_model = GenericNoise.from_paulis([0.000075, 0.000075, 0.000075], rng=rng)
# two_qubit_noise_model = GenericNoise.from_paulis([0.00039, 0.00039, 0.00039], rng=rng)
