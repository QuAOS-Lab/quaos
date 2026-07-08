"""Settings handles for running the cost-aware surface method from common_run_models."""

from __future__ import annotations

from dataclasses import dataclass

from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementOutcomes,
    MeasurementRequest,
    RMBBackend,
    ShotRNG,
)
from sympleq.applications.randomized_benchmarking.backends.quantinuum import (
    QuantinuumBackend,
)
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    BASE_SIMULATION_COST,
    CrossingSettings,
    default_backend_factory,
    pytket_bare_simulation_cost,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import (
    CostAwareSurfaceSettings,
    EXPERIMENTS_DIR,
)


def quantinuum_h2_backend_factory(
    settings: CrossingSettings,
    rng: RNGGenerator,
) -> RMBBackend:
    return QuantinuumBackend(
        device_name=H2_DEVICE_NAME,
        project_name=PROJECT_NAME,
        batch_size=1,
        max_cost_per_run=settings.max_cost_per_run,
    )


class GateLimitedQuantinuumBackend(QuantinuumBackend):
    """Quantinuum backend variant that caps actual stitched circuit gates."""

    def __init__(
        self,
        *args,
        max_stitched_gates: int | None = None,
        max_emulator_batch_cost: float | None = None,
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)
        self.max_stitched_gates = max_stitched_gates
        self.max_emulator_batch_cost = max_emulator_batch_cost

    def fidelity_estimation(
        self,
        requests: list[MeasurementRequest],
        rng: RNGGenerator,
        shot_rng: ShotRNG | None = None,
    ) -> MeasurementOutcomes:
        from sympleq.integrations.quantinuum.stitching import (
            MAX_QASM_PROGRAM_SIZE,
            circuit_stitching,
            destitch_results,
            estimate_qasm_program_size,
        )
        from sympleq.integrations.quantinuum.workflow import run_circuits_on_device

        jobs = []
        shot_counts: dict[RMBConfig, int] = {}
        for request in requests:
            for _ in range(max(0, request.shots)):
                index = shot_counts.get(request.config, 0)
                shot_counts[request.config] = index + 1
                circuit_rng = rng if shot_rng is None else shot_rng(request.config, index)
                jobs.append((request.config, self.compatible_circuit(request.config, circuit_rng)))
        if not jobs:
            return MeasurementOutcomes()

        def stitched_gate_count(submission) -> int:
            return int(circuit_stitching([circuit for _, circuit in submission]).n_gates)

        logical_cost = float(BASE_SIMULATION_COST)
        for i, (config, circuit) in enumerate(jobs):
            single = [(config, circuit)]
            logical_cost += pytket_bare_simulation_cost(circuit)
            if i:
                logical_cost += config.n_qubits / 5000

            if (
                self.max_stitched_gates is not None
                and stitched_gate_count(single) > self.max_stitched_gates
            ):
                raise ValueError(
                    "Single emulator circuit exceeds gate_budget "
                    f"({stitched_gate_count(single)} > {self.max_stitched_gates})."
                )

        def submission_cost(submission) -> float:
            cost = float(BASE_SIMULATION_COST)
            for i, (config, circuit) in enumerate(submission):
                cost += pytket_bare_simulation_cost(circuit)
                if i:
                    cost += config.n_qubits / 5000
            return cost

        def submission_size(submission) -> int:
            return int(sum(estimate_qasm_program_size(circuit) for _, circuit in submission))

        def fits(submission) -> bool:
            if not submission:
                return True
            if submission_size(submission) > MAX_QASM_PROGRAM_SIZE:
                return False
            if (
                self.max_emulator_batch_cost is not None
                and submission_cost(submission) > self.max_emulator_batch_cost
            ):
                return False
            return (
                self.max_stitched_gates is None
                or stitched_gate_count(submission) <= self.max_stitched_gates
            )

        def balanced_submissions() -> list[list[tuple[RMBConfig, object]]]:
            if fits(jobs):
                return [jobs]
            if self.max_stitched_gates is None:
                chunks = []
                current = []
                for job in jobs:
                    candidate = current + [job]
                    if current and not fits(candidate):
                        chunks.append(current)
                        current = []
                    current.append(job)
                if current:
                    chunks.append(current)
                return chunks

            total_gates = stitched_gate_count(jobs)
            parts = max(1, -(-total_gates // int(self.max_stitched_gates)))
            while parts <= len(jobs):
                target = total_gates / parts
                chunks = []
                index = 0
                ok = True
                for remaining_parts in range(parts, 0, -1):
                    current = []
                    while index < len(jobs):
                        candidate = current + [jobs[index]]
                        if not fits(candidate):
                            if not current:
                                ok = False
                            break

                        must_leave_one_per_remaining = (
                            len(jobs) - (index + 1) < remaining_parts - 1
                        )
                        if current and not must_leave_one_per_remaining:
                            current_gates = stitched_gate_count(current)
                            candidate_gates = stitched_gate_count(candidate)
                            if abs(candidate_gates - target) > abs(current_gates - target):
                                break

                        current = candidate
                        index += 1
                        if len(jobs) - index == remaining_parts - 1:
                            break

                    if not current:
                        ok = False
                        break
                    chunks.append(current)

                if ok and index == len(jobs):
                    return chunks
                parts += 1

            return [[job] for job in jobs]

        submissions = balanced_submissions()

        stitched_circuits = [
            circuit_stitching([circuit for _, circuit in submission])
            for submission in submissions
        ]
        results = run_circuits_on_device(
            stitched_circuits,
            self.n_shots,
            self.device_name,
            self.project_name,
            verbose=True,
        )

        outcomes: dict[RMBConfig, list[bool]] = {}
        for submission, result, circuit in zip(submissions, results, stitched_circuits):
            registers = sorted(
                circuit.c_registers,
                key=lambda register: int(register.name.removeprefix("creg_")),
            )
            for (config, _), sub_result in zip(submission, destitch_results(result, registers)):
                counts = sub_result.get_empirical_distribution().as_counter()
                if counts:
                    top_outcome, _ = counts.most_common()[0]
                    outcomes.setdefault(config, []).append(all(bit == 0 for bit in top_outcome))

        return MeasurementOutcomes(
            outcomes=outcomes,
            cost=logical_cost,
            n_submissions=1,
        )


def cost_aware_emulator_backend_factory(
    settings: CrossingSettings,
    rng: RNGGenerator,
) -> RMBBackend:
    emulator_max_batch_cost = getattr(settings, "emulator_max_batch_cost", None)
    return GateLimitedQuantinuumBackend(
        device_name="H2-Emulator",
        project_name=PROJECT_NAME,
        batch_size=1,
        max_cost_per_run=(
            settings.max_cost_per_run
            if emulator_max_batch_cost is None
            else emulator_max_batch_cost
        ),
        max_stitched_gates=getattr(settings, "single_batch_gate_budget", None),
        max_emulator_batch_cost=emulator_max_batch_cost,
    )


@dataclass(frozen=True)
class CostAwareSettings(CostAwareSurfaceSettings):
    save_real_checkpoints: bool = True
    save_gp_prediction_grid: bool = False
    single_batch_gate_budget: int | None = None
    emulator_max_batch_cost: float | None = None


# =============================================================================
# CONTROL PANEL
# =============================================================================

BACKEND = 'emulator'  # "emulator"  # "sympleq", "emulator", or "H2"
PROJECT_NAME = "fidelity-benchmark"
H2_DEVICE_NAME = "H2-1"

Q_VALUES = tuple(range(5, 21, 1))
N_QUBITS = round(sum(Q_VALUES) / len(Q_VALUES))
N_GATES_BOUNDS = (10, 3000)
RATIO_BOUNDS = (0.1, 0.9)

HQC_BUDGET = 250.0
MAX_COST_PER_RUN = 35.0

# Only emulator
GATE_BUDGET = 20000
SINGLE_BATCH_GATE_BUDGET = 4200
EMULATOR_MAX_BATCH_COST = 12.0


RNG_SEEDS = [2025]

ACQUISITION_Q_RESOLUTION = 15 # 30
ACQUISITION_RATIO_POINTS = 15
MAX_QUBIT_WINDOW = 5
MIN_DISTINCT_Q_COVERAGE = 10

GRID_RESOLUTION = (11, 11, 7, 7, 1, 5, 1)
BOUNDARY_FIT_RESOLUTION = (17, 17, 9, 9, 1, 7, 1)
CONTINUOUS_REFIT = True

PLOT = False
SURFACE_PLOT_SHOW = True
LIVE_SURFACE_PLOT = True
LIVE_SURFACE_PLOT_SHOW = True
LIVE_VOLUME_PLOT = True
LIVE_VOLUME_PLOT_SHOW = True
LIVE_PLOT_PAUSE = 0.5

USE_SCRAMBLER = True
VERBOSE = True
SAVE_REAL_CHECKPOINTS = True

GP_GRID_SURFACE_PATH = None #(
#    EXPERIMENTS_DIR.parent / "rmb_data" / "FLE_20260702_063751_gp_grid_3d.npz"
#)
GP_GRID_SURFACE_LABEL = None # "Rick/Shreya Grid"


def backend_factory_for_name(name: str):
    if name == "sympleq":
        return default_backend_factory
    if name == "emulator":
        return cost_aware_emulator_backend_factory
    if name == "H2":
        return quantinuum_h2_backend_factory
    raise ValueError(f"Unknown cost-aware backend {name!r}")


def control_panel_settings_kwargs() -> dict:
    return dict(
        q_values=Q_VALUES,
        n_qubits=N_QUBITS,
        n_gates_bounds=N_GATES_BOUNDS,
        ratio_bounds=RATIO_BOUNDS,
        hqc_budget=HQC_BUDGET,
        max_cost_per_run=MAX_COST_PER_RUN,
        gate_budget=GATE_BUDGET,
        single_batch_gate_budget=SINGLE_BATCH_GATE_BUDGET,
        emulator_max_batch_cost=EMULATOR_MAX_BATCH_COST,
        backend_model=BACKEND,
        backend_factory=backend_factory_for_name(BACKEND),
        acquisition_q_resolution=ACQUISITION_Q_RESOLUTION,
        acquisition_ratio_points=ACQUISITION_RATIO_POINTS,
        max_qubit_window=MAX_QUBIT_WINDOW,
        min_distinct_q_coverage=MIN_DISTINCT_Q_COVERAGE,
        grid_resolution=GRID_RESOLUTION,
        boundary_fit_resolution=BOUNDARY_FIT_RESOLUTION,
        continuous_refit=CONTINUOUS_REFIT,
        plot=PLOT,
        surface_plot_show=SURFACE_PLOT_SHOW,
        surface_plot_n_gates_bounds=N_GATES_BOUNDS,
        live_surface_plot=LIVE_SURFACE_PLOT,
        live_surface_plot_show=LIVE_SURFACE_PLOT_SHOW,
        live_surface_plot_pause=LIVE_PLOT_PAUSE,
        live_volume_plot=LIVE_VOLUME_PLOT,
        live_volume_plot_show=LIVE_VOLUME_PLOT_SHOW,
        live_volume_plot_pause=LIVE_PLOT_PAUSE,
        gp_grid_surface_path=GP_GRID_SURFACE_PATH,
        gp_grid_surface_label=GP_GRID_SURFACE_LABEL,
        verbose=VERBOSE,
        use_scrambler=USE_SCRAMBLER,
        save_real_checkpoints=SAVE_REAL_CHECKPOINTS,
        save_gp_prediction_grid=False,
    )
