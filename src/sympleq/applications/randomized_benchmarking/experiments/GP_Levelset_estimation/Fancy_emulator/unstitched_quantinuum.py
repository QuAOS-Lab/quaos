from __future__ import annotations

import qnexus as qnx
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementOutcomes,
    MeasurementRequest,
    ShotRNG,
)
from sympleq.applications.randomized_benchmarking.backends.quantinuum import (
    QuantinuumBackend,
)
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.integrations.quantinuum.utils import (
    BASE_SIMULATION_COST,
    pytket_bare_simulation_cost,
)


def fancy_noisy_emulator_config(device_name: str) -> qnx.QuantinuumConfig:
    return qnx.QuantinuumConfig(
        device_name=device_name,
        simulator="state-vector",
        noisy_simulation=True,
        no_opt=True,
        allow_implicit_swaps=False,
        leakage_detection=False,
        max_batch_cost=0.0,
        attempt_batching=False,
    )


fancy_ideal_emulator_config = fancy_noisy_emulator_config


class FancyUnstitchedQuantinuumBackend(QuantinuumBackend):
    """Quantinuum emulator backend that submits each circuit as its own program."""

    log_prefix = "[fancy-unstitched]"
    program_name_prefix = "fancy-unstitched"
    simulator = "state-vector"
    noisy_simulation = True

    def to_dict(self) -> dict:
        payload = super().to_dict()
        payload.update(
            {
                "project_name": self.project_name,
                "max_cost_per_run": self.max_cost_per_run,
                "stitching": "unstitched",
                "simulator": self.simulator,
                "noisy_simulation": self.noisy_simulation,
                "fancy_emulator": True,
            }
        )
        return payload

    def fidelity_estimation(
        self,
        requests: list[MeasurementRequest],
        rng: RNGGenerator,
        shot_rng: ShotRNG | None = None,
    ) -> MeasurementOutcomes:
        from sympleq.integrations.quantinuum.workflow import (
            build_and_compile_circuits,
            run_compiled_circuits,
            setup,
        )

        jobs = []
        shot_counts: dict[RMBConfig, int] = {}
        for request_index, request in enumerate(requests):
            for _ in range(max(0, request.shots)):
                shot_index = shot_counts.get(request.config, 0)
                shot_counts[request.config] = shot_index + 1
                circuit_rng = (
                    rng
                    if shot_rng is None
                    else shot_rng(request.config, shot_index)
                )
                circuit = self.compatible_circuit(request.config, circuit_rng)
                jobs.append((request.config, circuit))
                print(
                    f"{self.log_prefix} "
                    f"request={request_index:02d} "
                    f"shot_index={shot_index} "
                    f"qubits={circuit.n_qubits} "
                    f"bits={circuit.n_bits} "
                    f"gates={circuit.n_gates} "
                    f"n_1q={circuit.n_1qb_gates()} "
                    f"n_2q={circuit.n_2qb_gates()}",
                    flush=True,
                )

        if not jobs:
            return MeasurementOutcomes()

        submissions = []
        current = []
        current_cost = float(BASE_SIMULATION_COST)
        for config, circuit in jobs:
            append_cost = pytket_bare_simulation_cost(circuit)
            if current and current_cost + append_cost > self.max_cost_per_run:
                submissions.append(current)
                current = []
                current_cost = float(BASE_SIMULATION_COST)
            current.append((config, circuit))
            current_cost += append_cost
        if current:
            submissions.append(current)

        outcomes: dict[RMBConfig, list[bool]] = {}
        total_cost = 0.0
        backend_config = fancy_noisy_emulator_config(self.device_name)
        setup(self.project_name)
        for submission_index, submission in enumerate(submissions):
            circuits = [circuit for _, circuit in submission]
            submission_cost = float(BASE_SIMULATION_COST) + sum(
                pytket_bare_simulation_cost(circuit)
                for circuit in circuits
            )
            total_cost += submission_cost
            print(
                f"{self.log_prefix} "
                f"submitting batch={submission_index:02d} "
                f"programs={len(circuits)} "
                f"estimated_cost={submission_cost:.4f}",
                flush=True,
            )

            print(
                f"{self.log_prefix} "
                f"backend_config device={self.device_name} "
                f"simulator={self.simulator} "
                f"noisy_simulation={self.noisy_simulation}",
                flush=True,
            )
            ref_circuits = build_and_compile_circuits(
                circuits,
                backend_config=backend_config,
                name=f"{self.program_name_prefix}-{submission_index:02d}",
            )
            results = run_compiled_circuits(
                ref_circuits,
                self.n_shots,
                backend_config,
            )
            if len(results) != len(submission):
                raise RuntimeError(
                    "Unstitched Quantinuum result count mismatch: "
                    f"got {len(results)} results for {len(submission)} circuits."
                )

            for result_index, ((config, _), result) in enumerate(
                zip(submission, results)
            ):
                counts = result.get_empirical_distribution().as_counter()
                if not counts:
                    continue
                top_outcome, _ = counts.most_common()[0]
                success = all(bit == 0 for bit in top_outcome)
                outcomes.setdefault(config, []).append(success)
                print(
                    f"{self.log_prefix} "
                    f"result={result_index:02d} "
                    f"q={config.n_qubits} "
                    f"n_1q={config.n_1qb_gates} "
                    f"n_2q={config.n_2qb_gates} "
                    f"success={success}",
                    flush=True,
                )

        return MeasurementOutcomes(
            outcomes=outcomes,
            cost=total_cost,
            n_submissions=len(submissions),
        )


def _fancy_unstitched_backend(settings, device_name: str, project_name: str):
    return FancyUnstitchedQuantinuumBackend(
        device_name=device_name,
        project_name=project_name,
        batch_size=1,
        max_cost_per_run=settings.max_cost_per_run,
    )


def fancy_unstitched_h21e(settings, rng):
    return _fancy_unstitched_backend(
        settings,
        device_name="H2-1E",
        project_name="Fancy_Emulator_Unstitched_H21E",
    )


def fancy_unstitched_h22e(settings, rng):
    return _fancy_unstitched_backend(
        settings,
        device_name="H2-2E",
        project_name="Fancy_Emulator_Unstitched_H22E",
    )


def fancy_unstitched_emulator_backend_factory(settings, rng):
    """Backward-compatible H2-1E factory."""

    return fancy_unstitched_h21e(settings, rng)
