"""Settings handles for running the cost-aware surface method from common_run_models."""

from __future__ import annotations

from dataclasses import dataclass
from functools import partial
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

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
    BASE_1Q_PAULI_ERROR,
    BASE_2Q_PAULI_ERROR,
    CostAwareSurfaceSettings,
    EXPERIMENTS_DIR,
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


def cost_aware_quantinuum_backend_factory(
    settings: CrossingSettings,
    rng: RNGGenerator,
    *,
    device_name: str,
    gate_limited: bool,
) -> RMBBackend:
    if not gate_limited:
        return QuantinuumBackend(
            device_name=device_name,
            project_name=PROJECT_NAME,
            batch_size=1,
            max_cost_per_run=settings.max_cost_per_run,
        )

    emulator_max_batch_cost = getattr(settings, "emulator_max_batch_cost", None)
    return GateLimitedQuantinuumBackend(
        device_name=device_name,
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

BACKEND = "H2-1"  # "sympleq"  # , "emulator", "H2-1", "H2-2", "H2-1E", or "H2-2E"
PROJECT_NAME = "fidelity-benchmark"

Q_VALUES = tuple(range(27, 57, 1))
N_QUBITS = round(sum(Q_VALUES) / len(Q_VALUES))
N_GATES_BOUNDS = (10, 4000)
RATIO_BOUNDS = (0.1, 0.9)

HQC_BUDGET = 500.0
MAX_COST_PER_RUN = 35.0

# Only gate-limited emulator targets: "emulator", "H2-1E", and "H2-2E".
GATE_BUDGET = 20000
SINGLE_BATCH_GATE_BUDGET = 4200
EMULATOR_MAX_BATCH_COST = 12.0


RNG_SEEDS = [20261]

ACQUISITION_Q_RESOLUTION = 30
ACQUISITION_RATIO_POINTS = 15
MAX_QUBIT_WINDOW = 5
MIN_DISTINCT_Q_COVERAGE = 30
# Tuple order:
#   (L1, L2, m1, m2, nu1, nu2, V)
# where L1/L2 are the base one-/two-qubit rates, m1/m2 are their linear
# Q-slopes, nu1/nu2 are their quadratic Q-curvatures, and V is visibility.
# GRID_RESOLUTION = (11, 11, 7, 7, 1, 1, 1)  
# BOUNDARY_FIT_RESOLUTION = (13, 13, 7, 7, 1, 1, 1)
GRID_RESOLUTION = (15, 15, 9, 9, 1, 9, 1)  
BOUNDARY_FIT_RESOLUTION = (17, 17, 9, 9, 1, 9, 1)
CONTINUOUS_REFIT = True

# Prior centre for the one- and two-qubit Pauli error guesses. The surface
# method converts these to the initial per-gate rate guesses used by the grid.
INITIAL_ONE_Q_PAULI_ERROR = BASE_1Q_PAULI_ERROR
INITIAL_TWO_Q_PAULI_ERROR = BASE_2Q_PAULI_ERROR
INITIAL_ERROR_RELATIVE_UNCERTAINTY = 0.30
INITIAL_ONE_Q_ERROR_RELATIVE_UNCERTAINTY = None  # overrides INITIAL_ERROR_RELATIVE_UNCERTAINTY if necessary
INITIAL_TWO_Q_ERROR_RELATIVE_UNCERTAINTY = None

# Visibility/amplitude nuisance parameter. Leave VISIBILITY_BOUNDS as None for
# the historical physical behaviour: if the V axis is freed, V is log-grid
# sampled and clipped to 0 < V <= 1. For the diagnostic transient-absorption
# test, set e.g. VISIBILITY_BOUNDS = (0.95, 1.05) and make the final entry of
# GRID_RESOLUTION / BOUNDARY_FIT_RESOLUTION greater than 1.
INITIAL_VISIBILITY = 1.0
VISIBILITY_LOG_STD = 0.20
VISIBILITY_BOUNDS = None  # (0.95, 1.05)  # e.g. (0.95, 1.05)

# Optional randomisation of the initial prior centre. When enabled, each
# rng_seed gets a deterministic additive Gaussian perturbation around the base
# values above:
#
#     e_i = e_i_base + Normal(0, relative_width_i * e_i_base).
#
# The realised values are clipped positive, saved in checkpoints, and used by
# replay.
RANDOMIZE_INITIAL_RATE_GUESSES = True
INITIAL_ONE_Q_RATE_RANDOM_RELATIVE_STD = 0.3
INITIAL_TWO_Q_RATE_RANDOM_RELATIVE_STD = 0.3
INITIAL_RATE_RANDOM_SEED_OFFSET = 1729

PLOT = False
SURFACE_PLOT_SHOW = True
SURFACE_PLOT_LINDBLAD = True
LIVE_SURFACE_PLOT = True
LIVE_SURFACE_PLOT_SHOW = True
LIVE_VOLUME_PLOT = True
LIVE_VOLUME_PLOT_SHOW = True
LIVE_PLOT_PAUSE = 0.5
USE_VOXEL_VOLUME = True
VOXEL_VOLUME_N_GATES_GRID = 80
VOXEL_VOLUME_N_RATIO_GRID = 80
VOXEL_VOLUME_N_QUBITS_GRID = 16

USE_SCRAMBLER = True
VERBOSE = True
SAVE_REAL_CHECKPOINTS = True

GP_GRID_SURFACE_PATH = (
    Path("Personal")
    / "CostAware"
    / "reference_grids"
    / "measurement_015_globalsur_20260710_105239_496166_gp_grid_3d.npz"
)
GP_GRID_SURFACE_LABEL = "measurement 015 GlobalSUR grid"

_GATE_LIMITED_QUANTINUUM_BACKENDS = {
    "emulator": "H2-Emulator",
    "H2-1E": "H2-1E",
    "H2-2E": "H2-2E",
}
_QUANTINUUM_HARDWARE_BACKENDS = {"H2-1", "H2-2"}


def backend_factory_for_name(name: str):
    if name == "sympleq":
        return default_backend_factory
    if name in _GATE_LIMITED_QUANTINUUM_BACKENDS:
        factory = partial(
            cost_aware_quantinuum_backend_factory,
            device_name=_GATE_LIMITED_QUANTINUUM_BACKENDS[name],
            gate_limited=True,
        )
        factory.__name__ = f"cost_aware_{name.replace('-', '_')}_backend_factory"
        return factory
    if name in _QUANTINUUM_HARDWARE_BACKENDS:
        factory = partial(
            cost_aware_quantinuum_backend_factory,
            device_name=name,
            gate_limited=False,
        )
        factory.__name__ = f"cost_aware_{name.replace('-', '_')}_backend_factory"
        return factory
    supported = ("sympleq", "emulator", "H2-1", "H2-2", "H2-1E", "H2-2E")
    raise ValueError(
        f"Unknown cost-aware backend {name!r}; choose one of {', '.join(supported)}."
    )


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
        initial_one_q_pauli_error=INITIAL_ONE_Q_PAULI_ERROR,
        initial_two_q_pauli_error=INITIAL_TWO_Q_PAULI_ERROR,
        initial_error_relative_uncertainty=INITIAL_ERROR_RELATIVE_UNCERTAINTY,
        initial_one_q_error_relative_uncertainty=INITIAL_ONE_Q_ERROR_RELATIVE_UNCERTAINTY,
        initial_two_q_error_relative_uncertainty=INITIAL_TWO_Q_ERROR_RELATIVE_UNCERTAINTY,
        initial_visibility=INITIAL_VISIBILITY,
        visibility_log_std=VISIBILITY_LOG_STD,
        visibility_bounds=VISIBILITY_BOUNDS,
        plot=PLOT,
        surface_plot_show=SURFACE_PLOT_SHOW,
        surface_plot_analytic=SURFACE_PLOT_LINDBLAD,
        surface_plot_n_gates_bounds=N_GATES_BOUNDS,
        live_surface_plot=LIVE_SURFACE_PLOT,
        live_surface_plot_show=LIVE_SURFACE_PLOT_SHOW,
        live_surface_plot_pause=LIVE_PLOT_PAUSE,
        live_volume_plot=LIVE_VOLUME_PLOT,
        live_volume_plot_show=LIVE_VOLUME_PLOT_SHOW,
        live_volume_plot_pause=LIVE_PLOT_PAUSE,
        use_voxel_volume=USE_VOXEL_VOLUME,
        voxel_volume_n_gates_grid=VOXEL_VOLUME_N_GATES_GRID,
        voxel_volume_n_ratio_grid=VOXEL_VOLUME_N_RATIO_GRID,
        voxel_volume_n_qubits_grid=VOXEL_VOLUME_N_QUBITS_GRID,
        gp_grid_surface_path=GP_GRID_SURFACE_PATH,
        gp_grid_surface_label=GP_GRID_SURFACE_LABEL,
        verbose=VERBOSE,
        use_scrambler=USE_SCRAMBLER,
        save_real_checkpoints=SAVE_REAL_CHECKPOINTS,
        save_gp_prediction_grid=False,
    )


def with_randomized_initial_rate_guesses(kwargs: dict, rng_seed: int | None) -> dict:
    """Apply deterministic per-seed randomisation to initial rate guesses.

    The randomisation is additive Gaussian in the Pauli-error guesses:

        e_i = e_i_base + Normal(0, relative_width_i * e_i_base).

    ``rng_seed`` plus ``INITIAL_RATE_RANDOM_SEED_OFFSET`` defines the draw, so
    rerunning the same seed gives the same prior centre. The realised positive
    values, relative deltas, and equivalent log multipliers are carried in the
    settings and written to checkpoints.
    """
    updated = dict(kwargs)
    base_one = float(updated["initial_one_q_pauli_error"])
    base_two = float(updated["initial_two_q_pauli_error"])
    random_seed = (
        None
        if rng_seed is None
        else int(rng_seed) + int(INITIAL_RATE_RANDOM_SEED_OFFSET)
    )
    one_log = 0.0
    two_log = 0.0
    one_delta = 0.0
    two_delta = 0.0
    if RANDOMIZE_INITIAL_RATE_GUESSES:
        rng = default_rng(random_seed)
        one_delta = float(rng.normal(0.0, float(INITIAL_ONE_Q_RATE_RANDOM_RELATIVE_STD)))
        two_delta = float(rng.normal(0.0, float(INITIAL_TWO_Q_RATE_RANDOM_RELATIVE_STD)))
        one_error = max(base_one * (1.0 + one_delta), 1e-12)
        two_error = max(base_two * (1.0 + two_delta), 1e-12)
        one_log = float(np.log(one_error / base_one))
        two_log = float(np.log(two_error / base_two))
        updated["initial_one_q_pauli_error"] = one_error
        updated["initial_two_q_pauli_error"] = two_error

    updated["initial_rate_randomization_enabled"] = bool(
        RANDOMIZE_INITIAL_RATE_GUESSES
    )
    updated["initial_rate_random_seed"] = random_seed
    updated["initial_one_q_pauli_error_base"] = base_one
    updated["initial_two_q_pauli_error_base"] = base_two
    updated["initial_one_q_random_log_multiplier"] = one_log
    updated["initial_two_q_random_log_multiplier"] = two_log
    updated["initial_one_q_random_relative_std"] = float(
        INITIAL_ONE_Q_RATE_RANDOM_RELATIVE_STD
    )
    updated["initial_two_q_random_relative_std"] = float(
        INITIAL_TWO_Q_RATE_RANDOM_RELATIVE_STD
    )
    updated["initial_one_q_random_relative_delta"] = one_delta
    updated["initial_two_q_random_relative_delta"] = two_delta
    return updated
