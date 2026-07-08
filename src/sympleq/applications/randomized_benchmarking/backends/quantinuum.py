from __future__ import annotations
from numpy.random import Generator as RNGGenerator
from pytket.circuit import Circuit as PytketCircuit
import warnings

from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementOutcomes,
    MeasurementRequest,
    RMBBackend,
    ShotRNG,
)
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


MAX_COST_PER_RUN: float = 25.0


class QuantinuumBackend(RMBBackend):
    """
    Quantinuum hardware/emulator backend.

    Builds the requested random circuits, stitches them into as few
    submissions as the QASM program size limit and ``max_cost_per_run``
    allow, runs them on the Quantinuum device identified by
    ``device_name`` for ``n_shots`` shots each, and reports per circuit
    whether the most-common measurement outcome is the all-zeros
    bitstring.

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

    def fidelity_estimation(self, requests: list[MeasurementRequest], rng: RNGGenerator,
                            shot_rng: ShotRNG | None = None) -> MeasurementOutcomes:
        from sympleq.integrations.quantinuum.workflow import run_circuits_on_device
        from sympleq.integrations.quantinuum.stitching import circuit_stitching, destitch_results, \
            estimate_qasm_program_size, MAX_QASM_PROGRAM_SIZE

        # One independently drawn circuit per requested shot, in request order.
        jobs: list[tuple[RMBConfig, PytketCircuit]] = []
        shot_counts: dict[RMBConfig, int] = {}
        for request in requests:
            for _ in range(max(0, request.shots)):
                index = shot_counts.get(request.config, 0)
                shot_counts[request.config] = index + 1
                circuit_rng = rng if shot_rng is None else shot_rng(request.config, index)
                jobs.append((request.config,
                             self.compatible_circuit(request.config, circuit_rng)))
        if not jobs:
            return MeasurementOutcomes()

        # Pack the circuits greedily into stitched submissions, respecting the
        # QASM program size limit and the per-submission cost cap. Each circuit
        # stitched after another pays the reset cost of its qubits.
        submissions: list[list[tuple[RMBConfig, PytketCircuit]]] = []
        cost = 0.0
        current: list[tuple[RMBConfig, PytketCircuit]] = []
        current_size = 0
        current_cost = float(BASE_SIMULATION_COST)
        for config, circuit in jobs:
            program_size = estimate_qasm_program_size(circuit)
            append_cost = pytket_bare_simulation_cost(circuit)
            if current:
                append_cost += config.n_qubits / 5000
                if (current_size + program_size > MAX_QASM_PROGRAM_SIZE
                        or current_cost + append_cost > self.max_cost_per_run):
                    submissions.append(current)
                    cost += current_cost
                    current = []
                    current_size = 0
                    current_cost = float(BASE_SIMULATION_COST)
                    append_cost = pytket_bare_simulation_cost(circuit)
            current.append((config, circuit))
            current_size += program_size
            current_cost += append_cost
        submissions.append(current)
        cost += current_cost

        stitched_circuits = [circuit_stitching([circuit for _, circuit in submission])
                             for submission in submissions]
        results = run_circuits_on_device(
            stitched_circuits, self.n_shots, self.device_name, self.project_name, verbose=True)

        outcomes: dict[RMBConfig, list[bool]] = {}
        for submission, result, circuit in zip(submissions, results, stitched_circuits):
            # pytket lists registers lexicographically (creg_10 before creg_2);
            # destitching must read them in stitch order to keep each result
            # paired with the config whose circuit wrote it.
            registers = sorted(circuit.c_registers,
                               key=lambda register: int(register.name.removeprefix("creg_")))
            sorted_submission = sorted(submission, key=lambda item: item[1].n_qubits, reverse=True)
            unstitched_results = destitch_results(result, registers)
            print("\n[backend destitched order]", flush=True)
            for i, ((config, _), sub_result) in enumerate(zip(
                sorted_submission,
                unstitched_results)):
                counts = sub_result.get_empirical_distribution().as_counter()
                print(
                    f"destitched_index={i:02d} "
                    f"config_n_qubits={config.n_qubits} "
                    f"config_n_gates={config.n_gates} "
                    f"register={registers[i].name} "
                    f"counts={counts}",
                    flush=True,
                )
            for (config, _), sub_result in zip(
                sorted_submission,
                unstitched_results):
                counts = sub_result.get_empirical_distribution().as_counter()
                if not counts:
                    continue
                # The shot's outcome is whether the most-common measured
                # bitstring is the all-zeros initial state.
                top_outcome, _ = counts.most_common()[0]
                outcomes.setdefault(config, []).append(all(bit == 0 for bit in top_outcome))

        return MeasurementOutcomes(outcomes=outcomes, cost=cost,
                                   n_submissions=len(submissions))

    @staticmethod
    def compatible_circuit(config: RMBConfig, rng: RNGGenerator) -> PytketCircuit:
        """Random pytket circuit of ``config`` whose gate counts match the config."""
        while True:
            circuit = to_pytket_circuit(config.random_circuit(rng=rng))
            if circuit.n_gates == 0:
                continue
            if circuit.n_1qb_gates() != config.n_1qb_gates:
                warnings.warn(f"Generated circuit has mismatching number of 1qb gates "
                              f"({circuit.n_1qb_gates()} vs {config.n_1qb_gates}).")
                continue
            if circuit.n_2qb_gates() != config.n_2qb_gates:
                warnings.warn(f"Generated circuit has mismatching number of 2qb gates "
                              f"({circuit.n_2qb_gates()} vs {config.n_2qb_gates}).")
                continue
            return circuit

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
    def default_sympleq_backend(cls, rng: RNGGenerator | None = None) -> SympleqBackend:
        noise_model = GenericNoise.from_paulis([0.000025, 0.000025, 0.000025], rng)
        two_qubit_noise_model = GenericNoise.from_paulis([0.00079, 0.00079, 0.00079], rng)
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
