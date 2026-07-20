from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData

from viarregio2 import (
    config_from_parameters,
    fidelity_mean,
    fidelity_variance,
    make_backend,
    spend_measurements,
    total_measurements,
)
from viarregio3 import (
    finite_difference_gradient,
    project_to_contour,
)
from viarregio4 import (
    MonotoneBoundaryExperimentConfig,
    expected_gate_counts,
    fit_monotone_fidelity_surface,
    hqc_cost,
    plot_monotone_fidelity_surface_with_confidence,
    print_experiment_summary,
    print_fit_reports,
)


@dataclass(frozen=True)
class ContourFirstExperimentConfig(MonotoneBoundaryExperimentConfig):
    """
    HQC-aware contour-first boundary estimation.

    The method first finds one p=0.5 anchor on a cheap low-ratio ray, then
    traces the contour by stepping along the current contour direction.  The
    trace phase measures projected contour points directly and only applies
    short local corrections when the projected point is off the boundary.
    """
    ray_ratio_count: int = 5
    ray_probe_shots: int = 2
    ray_bisection_steps: int = 5
    ray_bisection_shots: int = 2
    crossing_decision_confirm_shots: int = 4
    crossing_decision_probability_width: float = 0.20
    crossing_ambiguous_as_failure: bool = True
    initial_crossing_interior_bracket: bool = True
    initial_crossing_center_fraction: float = 0.35
    initial_crossing_half_width_fraction: float = 0.18
    initial_crossing_expand_factor: float = 1.6
    crossing_confirm_candidates: int = 0
    crossing_confirm_shots: int = 4
    initial_anchor_refine: bool = True
    initial_anchor_refine_shots: int = 6
    initial_anchor_refine_steps: int = 3
    initial_anchor_search_fraction: float = 0.06
    initial_anchor_min_runs: int = 8
    batched_initial_anchor_search: bool = True
    initial_anchor_depth_grid_count: int = 17
    initial_anchor_grid_shots: int = 2
    trace_reserve_budget_fraction: float = 0.35
    trace_reserve_min_hqc: float = 25.0
    trace_ratio_step_fraction: float = 0.06
    trace_depth_search_fraction: float = 0.06
    trace_step_shrink_attempts: int = 4
    trace_local_crossing: bool = True
    trace_local_stencil_points: int = 7
    trace_local_stencil_shots: int = 2
    seeded_initial_trace: bool = False
    seeded_trace_ratio_count: int = 4
    seeded_trace_max_batches: int = 2
    seeded_trace_ratio_span_fraction: float = 0.18
    seeded_trace_depth_power: float = 0.55
    high_ratio_projected_trace: bool = False
    high_ratio_projected_threshold: float = 0.55
    high_ratio_low_depth_fraction: float = 0.25
    trace_enforce_monotone_depth: bool = True
    trace_correction_steps: int = 2
    trace_shots: int = 2
    trace_accept_probability_width: float = 0.12
    trace_reject_probability_width: float = 0.25
    trace_directions: tuple[int, ...] = (1,)
    model_projection_after_fit: bool = False
    refine_after_trace: bool = True
    trace_anchor_min_shots: int = 8
    refinement_shots: int = 2
    refinement_boundary_width: float = 0.15
    diagnostics_path: str | Path | None = "viarregio7_diagnostics.json"
    print_diagnostics: bool = False
    save_path: str | Path | None = "viarregio7_boundary.json"
    batching_enabled: bool = True
    batch_max_configs: int = 42
    max_cost_per_batch: float | None = 50.0
    batch_max_hqc_cost: float | None = None
    batch_reset_weight: float = 1.0
    batch_candidate_multiplier: int = 4
    batch_refinement_shots: int = 4
    batch_refinement_fit_passes: int = 3
    batch_acquisition_passes: int = 2
    batch_acquisition_ratio_count: int = 6
    batch_target_fill_fraction: float = 0.88
    batch_fill_repeats: bool = True
    batch_fill_max_shots_per_config: int | None = 8
    batch_discovery_fill_max_shots_per_config: int | None = 4
    batch_fill_all_stages: bool = False
    high_ratio_acquisition_fraction: float = 0.35
    high_ratio_candidate_fraction: float = 0.0
    contour_bracket_probe_shots: int = 2
    contour_bracket_depth_fractions: tuple[float, ...] = (0.04, 0.08, 0.12)
    contour_bracket_max_relative_depth: float = 0.20
    acquisition_require_bracket_straddle: bool = True
    acquisition_follow_trace: bool = True
    acquisition_trace_backtrack_fraction: float = 0.08
    acquisition_trace_extension_fraction: float = 0.16
    acquisition_trace_predict_points: int = 0
    batch_post_trace_reserve_fraction: float = 0.25


@dataclass
class BudgetState:
    remaining_measurements: int
    remaining_hqc: float
    circuit_executions: int = 0
    max_execution_repeats: int = 0
    batched_jobs: int = 0
    max_batched_job_size: int = 0
    stitched_hqc_spent: float = 0.0
    native_hqc_estimate: float = 0.0
    batch_config_ids: dict[str, list[int]] = field(default_factory=dict)

    def can_spend(self, reserve_hqc: float = 0.0) -> bool:
        return self.remaining_measurements > 0 and self.remaining_hqc > reserve_hqc


@dataclass(frozen=True)
class MeasurementRequest:
    config: RMBConfig
    requested_shots: int


@dataclass(frozen=True)
class MeasurementSpend:
    config: RMBConfig
    requested_shots: int
    spent_shots: int


def script_default_settings(**overrides) -> ContourFirstExperimentConfig:
    """
    Default runnable v7 profile used by the standalone script and benchmark defaults.
    """
    params = dict(
        measurement_budget=10000000,
        hqc_budget=500.0,
        n_qubits_values=(50,),
        depth_bounds=(4, 50),
        ratio_bounds=(0.08, 0.8),
        random_elimination=0.1,
        scrambling_probability=0.0,
        min_adaptive_shots_per_config=1,
        max_adaptive_shots_per_config=5,
        max_shots_per_config=10,
        candidate_grid_size=(80, 80),
        boundary_width=0.08,
        shot_boundary_width=0.2,
        surface_smoothing=0.1,
        min_fit_points=8,
        monotone_l2=1e-3,
        contour_ready_probability_width=0.10,
        contour_min_anchors=4,
        contour_step_fraction=0.08,
        contour_projection_fraction=0.12,
        contour_gradient_fraction=0.01,
        hqc_cost_informed_acquisition=True,
        hqc_cost_power=1.0,
        save_path="viarregio7_boundary.json",
        verbose=True,
    )
    params.update(overrides)
    return ContourFirstExperimentConfig(**params)


def max_cost_per_batch(settings: ContourFirstExperimentConfig) -> float | None:
    if settings.max_cost_per_batch is not None:
        return settings.max_cost_per_batch
    return settings.batch_max_hqc_cost


def batch_config_key(config: RMBConfig) -> str:
    return "|".join(
        [
            str(config.n_qubits),
            str(config.n_1qb_gates),
            str(config.n_2qb_gates),
            f"{float(config.scrambling_probability):.12g}",
            f"{float(config.random_elimination):.12g}",
        ]
    )


def target_batch_cost(settings: ContourFirstExperimentConfig, spendable_hqc: float) -> float:
    cap = max_cost_per_batch(settings)
    if cap is None:
        return spendable_hqc
    return min(spendable_hqc, max(0.0, settings.batch_target_fill_fraction) * cap)


DISCOVERY_BATCH_STAGES = {
    "single",
    "probe_batch",
    "initial",
    "initial_grid",
    "initial_refine",
    "seeded_trace",
    "trace",
    "trace_local",
    "trace_projected",
}

POST_TRACE_BATCH_STAGES = {
    "trace_anchor_backfill",
    "refinement",
    "boundary_acquisition",
}


def stage_allows_batch_fill(settings: ContourFirstExperimentConfig, stage: str) -> bool:
    if stage in POST_TRACE_BATCH_STAGES:
        return True
    if settings.batch_fill_all_stages:
        return True
    return stage not in DISCOVERY_BATCH_STAGES


def stage_batch_fill_target_runs(
    settings: ContourFirstExperimentConfig,
    stage: str,
    target_runs: int | None = None,
) -> int:
    if target_runs is not None:
        return target_runs
    if stage in DISCOVERY_BATCH_STAGES:
        cap = settings.batch_discovery_fill_max_shots_per_config
        if cap is None:
            return settings.max_shots_per_config
        return min(settings.max_shots_per_config, max(1, cap))
    cap = settings.batch_fill_max_shots_per_config
    if cap is None:
        return settings.max_shots_per_config
    return max(settings.max_shots_per_config, cap)


def trace_reserve_hqc(settings: ContourFirstExperimentConfig) -> float:
    if settings.hqc_budget is None:
        return 0.0
    return min(
        float(settings.hqc_budget),
        max(
            settings.trace_reserve_min_hqc,
            settings.trace_reserve_budget_fraction * float(settings.hqc_budget),
        ),
    )


def post_trace_reserve_hqc(settings: ContourFirstExperimentConfig) -> float:
    if not settings.refine_after_trace or settings.hqc_budget is None:
        return 0.0
    if settings.batch_post_trace_reserve_fraction <= 0.0:
        return 0.0
    return min(
        float(settings.hqc_budget),
        max(
            settings.hqc_base_cost,
            settings.batch_post_trace_reserve_fraction * float(settings.hqc_budget),
        ),
    )


def template_config(settings: ContourFirstExperimentConfig, n_qubits: int) -> RMBConfig:
    return (
        RMBConfig.default()
        .with_n_qubits(n_qubits)
        .with_random_elimination(settings.random_elimination)
        .with_scrambling_probability(settings.scrambling_probability)
    )


def measured_probability(data: RMBData, config: RMBConfig) -> float | None:
    if config not in data or data[config].num_runs() == 0:
        return None
    return fidelity_mean(data[config])


def measured_std(data: RMBData, config: RMBConfig) -> float:
    if config not in data or data[config].num_runs() == 0:
        return float("inf")
    return float(np.sqrt(fidelity_variance(data[config])))


def add_diagnostic(
    diagnostics: list[dict] | None,
    event: str,
    **values,
) -> None:
    if diagnostics is None:
        return
    clean_values = {}
    for key, value in values.items():
        if isinstance(value, np.generic):
            value = value.item()
        clean_values[key] = value
    diagnostics.append({"event": event, **clean_values})


def weighted_hqc_size(config: RMBConfig, settings: ContourFirstExperimentConfig) -> float:
    """
    Quantinuum-style weighted operation count for one RMB result circuit.
    """
    n_one, n_two, n_meas = expected_gate_counts(config)
    return (
        settings.hqc_one_qubit_weight * n_one
        + settings.hqc_two_qubit_weight * n_two
        + settings.hqc_measurement_weight * n_meas
    )


def batched_hqc_cost(
    requests: list[MeasurementRequest],
    settings: ContourFirstExperimentConfig,
) -> float:
    """
    Estimate the Quantinuum HQC cost if these local requests were stitched.

    Execution in v7 remains local through SympleQ. This model charges the
    Quantinuum flat cost once for a stitched job, adds all subcircuit weighted
    operation costs, and includes reset overhead between subcircuits.
    """
    active_requests = [request for request in requests if request.requested_shots > 0]
    if not active_requests:
        return 0.0
    if not settings.batching_enabled or len(active_requests) == 1:
        return sum(
            hqc_cost(request.config, request.requested_shots, settings)
            for request in active_requests
        )

    n_qubits_values = {request.config.n_qubits for request in active_requests}
    if len(n_qubits_values) != 1:
        return sum(
            hqc_cost(request.config, request.requested_shots, settings)
            for request in active_requests
        )

    total_shots = sum(request.requested_shots for request in active_requests)
    # Each requested shot is one RMB result circuit, so this includes one
    # measurement layer per circuit while charging the flat base only once.
    weighted_size = sum(
        weighted_hqc_size(request.config, settings) * request.requested_shots
        for request in active_requests
    )
    n_qubits = next(iter(n_qubits_values))
    reset_size = settings.batch_reset_weight * n_qubits * max(0, total_shots - 1)
    return settings.hqc_base_cost + (weighted_size + reset_size) / settings.hqc_scale


def native_quantinuum_hqc_cost(
    requests: list[MeasurementRequest],
    settings: ContourFirstExperimentConfig,
) -> float:
    """
    Estimate cost for running every measured circuit as a separate job.

    This is the unstitched baseline for the batching saving: each requested
    shot is one circuit and would pay the flat HQC base cost on its own.
    """
    return sum(
        hqc_cost(request.config, 1, settings) * request.requested_shots
        for request in requests
        if request.requested_shots > 0
    )


def request_marginal_hqc_cost(
    request: MeasurementRequest,
    settings: ContourFirstExperimentConfig,
) -> float:
    """
    Approximate the marginal stitched cost of one request inside a non-empty batch.
    """
    if request.requested_shots <= 0:
        return 0.0
    weighted_size = weighted_hqc_size(request.config, settings) * request.requested_shots
    reset_size = (
        settings.batch_reset_weight
        * request.config.n_qubits
        * request.requested_shots
    )
    return max(1e-12, (weighted_size + reset_size) / settings.hqc_scale)


def split_measurement_batches(
    requests: list[MeasurementRequest],
    settings: ContourFirstExperimentConfig,
    remaining_hqc: float,
    reserve_hqc: float,
) -> list[list[MeasurementRequest]]:
    batches: list[list[MeasurementRequest]] = []
    spendable_hqc = max(0.0, remaining_hqc - reserve_hqc)
    if spendable_hqc <= 0.0:
        return batches

    current: list[MeasurementRequest] = []
    current_n_qubits: int | None = None

    def can_add(batch: list[MeasurementRequest], request: MeasurementRequest) -> bool:
        candidate = batch + [request]
        if len(candidate) > max(1, settings.batch_max_configs):
            return False
        if current_n_qubits is not None and request.config.n_qubits != current_n_qubits:
            return False
        cost = batched_hqc_cost(candidate, settings)
        max_batch_cost = max_cost_per_batch(settings)
        if max_batch_cost is not None and cost > max_batch_cost:
            return False
        return cost <= spendable_hqc

    for request in requests:
        if request.requested_shots <= 0:
            continue
        if not current:
            max_batch_cost = max_cost_per_batch(settings)
            if batched_hqc_cost([request], settings) <= spendable_hqc and (
                max_batch_cost is None
                or batched_hqc_cost([request], settings) <= max_batch_cost
            ):
                current = [request]
                current_n_qubits = request.config.n_qubits
            continue
        if can_add(current, request):
            current.append(request)
            continue
        batches.append(current)
        spendable_hqc -= batched_hqc_cost(current, settings)
        if spendable_hqc <= 0.0:
            current = []
            current_n_qubits = None
            break
        max_batch_cost = max_cost_per_batch(settings)
        if batched_hqc_cost([request], settings) <= spendable_hqc and (
            max_batch_cost is None
            or batched_hqc_cost([request], settings) <= max_batch_cost
        ):
            current = [request]
            current_n_qubits = request.config.n_qubits
        else:
            current = []
            current_n_qubits = None

    if current:
        batches.append(current)
    return batches


def execute_local_sympleq_batch(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    requests: list[MeasurementRequest],
) -> list[MeasurementSpend]:
    """
    Execute a planned batch locally using the existing SympleQ backend.

    The batch is only a planning/accounting unit: no Quantinuum submission is
    made here. Each requested circuit is sampled locally and recorded under
    its RMBConfig.
    """
    spends: list[MeasurementSpend] = []
    for request in requests:
        spent = spend_measurements(
            backend=backend,
            rng=rng,
            data=data,
            config=request.config,
            n_measurements=request.requested_shots,
        )
        spends.append(MeasurementSpend(request.config, request.requested_shots, spent))
    return spends


def spend_configs_batch(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    requests: list[MeasurementRequest],
    budget: BudgetState,
    settings: ContourFirstExperimentConfig,
    reserve_hqc: float = 0.0,
    diagnostics: list[dict] | None = None,
    stage: str = "batch",
) -> list[MeasurementSpend]:
    if not requests or not budget.can_spend(reserve_hqc=reserve_hqc):
        return []

    clipped_requests: list[MeasurementRequest] = []
    remaining_measurements = budget.remaining_measurements
    for request in requests:
        if remaining_measurements <= 0:
            break
        shots = max(0, min(request.requested_shots, remaining_measurements))
        if shots > 0:
            clipped_requests.append(MeasurementRequest(request.config, shots))
            remaining_measurements -= shots

    if not clipped_requests:
        return []

    if settings.hqc_budget is None:
        batches = [clipped_requests]
    else:
        batches = split_measurement_batches(
            clipped_requests,
            settings=settings,
            remaining_hqc=budget.remaining_hqc,
            reserve_hqc=reserve_hqc,
        )
    spends: list[MeasurementSpend] = []

    for batch in batches:
        if not batch or not budget.can_spend(reserve_hqc=reserve_hqc):
            break
        if stage_allows_batch_fill(settings, stage):
            batch = fill_requests_toward_batch_cost(
                batch,
                data=data,
                settings=settings,
                spendable_hqc=max(0.0, budget.remaining_hqc - reserve_hqc),
                target_runs=stage_batch_fill_target_runs(settings, stage),
            )
            if not batch:
                break
        batch_cost_estimate = batched_hqc_cost(batch, settings)
        spendable_hqc = max(0.0, budget.remaining_hqc - reserve_hqc)
        batch_cap = max_cost_per_batch(settings)
        if batch_cost_estimate > spendable_hqc + 1e-9:
            break
        if batch_cap is not None and batch_cost_estimate > batch_cap + 1e-9:
            break
        before_runs = {
            request.config: data[request.config].num_runs() if request.config in data else 0
            for request in batch
        }
        batch_spends = execute_local_sympleq_batch(
            backend=backend,
            rng=rng,
            data=data,
            requests=batch,
        )
        spent_in_batch = sum(spend.spent_shots for spend in batch_spends)
        actual_requests = [
            MeasurementRequest(spend.config, spend.spent_shots)
            for spend in batch_spends
        ]
        spends.extend(batch_spends)

        if spent_in_batch <= 0:
            continue
        actual_cost = batched_hqc_cost(actual_requests, settings)
        native_cost = native_quantinuum_hqc_cost(actual_requests, settings)
        if actual_cost > spendable_hqc + 1e-9:
            raise RuntimeError(
                "Internal batching cost error: executed batch exceeded remaining HQC "
                f"budget ({actual_cost:.6f} > {spendable_hqc:.6f})."
            )
        budget.remaining_measurements -= spent_in_batch
        budget.circuit_executions += 1
        budget.batched_jobs += 1
        batch_id = budget.batched_jobs
        budget.stitched_hqc_spent += actual_cost
        budget.native_hqc_estimate += native_cost
        budget.max_execution_repeats = max(
            budget.max_execution_repeats,
            max((spend.spent_shots for spend in spends), default=0),
        )
        budget.max_batched_job_size = max(
            budget.max_batched_job_size,
            sum(request.requested_shots for request in batch),
        )
        if settings.hqc_budget is not None:
            budget.remaining_hqc = max(0.0, budget.remaining_hqc - actual_cost)
        else:
            budget.remaining_hqc = float(budget.remaining_measurements)
        for spend in batch_spends:
            if spend.spent_shots <= 0:
                continue
            key = batch_config_key(spend.config)
            ids = budget.batch_config_ids.setdefault(key, [])
            if batch_id not in ids:
                ids.append(batch_id)

        add_diagnostic(
            diagnostics,
            "batched_measurement",
            batch_id=batch_id,
            stage=stage,
            configs=len(batch),
            requested_shots=sum(request.requested_shots for request in batch),
            spent_shots=spent_in_batch,
            cost_model="quantinuum_stitched_estimate",
            execution_backend="local_sympleq",
            estimated_hqc=round(float(batch_cost_estimate), 6),
            actual_hqc=round(float(actual_cost), 6),
            native_hqc=round(float(native_cost), 6),
            batch_hqc_cap=None if batch_cap is None else round(float(batch_cap), 6),
            batch_cap_fraction=(
                None
                if batch_cap is None
                else round(float(actual_cost / batch_cap), 6)
            ),
            native_saving_hqc=round(
                float(native_cost - actual_cost),
                6,
            ),
            max_existing_runs=max(before_runs.values(), default=0),
        )

    return spends


def spend_config(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    requested_shots: int,
    budget: BudgetState,
    settings: ContourFirstExperimentConfig,
    reserve_hqc: float = 0.0,
    diagnostics: list[dict] | None = None,
    stage: str = "single",
) -> int:
    spends = spend_configs_batch(
        backend=backend,
        rng=rng,
        data=data,
        requests=[MeasurementRequest(config, requested_shots)],
        budget=budget,
        settings=settings,
        reserve_hqc=reserve_hqc,
        diagnostics=diagnostics,
        stage=stage,
    )
    return sum(spend.spent_shots for spend in spends)


def probe_depth(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    template: RMBConfig,
    depth: float,
    ratio: float,
    shots: int,
    budget: BudgetState,
    settings: ContourFirstExperimentConfig,
    reserve_hqc: float = 0.0,
    diagnostics: list[dict] | None = None,
    stage: str = "probe",
) -> tuple[RMBConfig, float | None]:
    config = config_from_parameters(template=template, depth=depth, ratio=ratio)
    spend_config(
        backend=backend,
        rng=rng,
        data=data,
        config=config,
        requested_shots=shots,
        budget=budget,
        settings=settings,
        reserve_hqc=reserve_hqc,
        diagnostics=diagnostics,
        stage=stage,
    )
    return config, measured_probability(data, config)


def probe_depths_batch(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    template: RMBConfig,
    probes: list[tuple[float, float, int]],
    budget: BudgetState,
    settings: ContourFirstExperimentConfig,
    reserve_hqc: float = 0.0,
    diagnostics: list[dict] | None = None,
    stage: str = "probe_batch",
) -> list[tuple[RMBConfig, float | None]]:
    configs = [
        config_from_parameters(template=template, depth=depth, ratio=ratio)
        for depth, ratio, _ in probes
    ]
    spend_configs_batch(
        backend=backend,
        rng=rng,
        data=data,
        requests=[
            MeasurementRequest(config, shots)
            for config, (_, _, shots) in zip(configs, probes)
        ],
        budget=budget,
        settings=settings,
        reserve_hqc=reserve_hqc,
        diagnostics=diagnostics,
        stage=stage,
    )
    return [(config, measured_probability(data, config)) for config in configs]


def confirm_crossing_decision_if_needed(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    probability: float,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    local: bool,
    reserve_hqc: float = 0.0,
    diagnostics: list[dict] | None = None,
    stage: str = "confirm_crossing",
) -> float:
    if local or settings.crossing_decision_confirm_shots <= 0:
        return probability
    if abs(probability - 0.5) > settings.crossing_decision_probability_width:
        return probability
    if config not in data:
        return probability

    current_runs = data[config].num_runs()
    extra_shots = settings.crossing_decision_confirm_shots - current_runs
    if extra_shots <= 0:
        return probability

    spend_config(
        backend=backend,
        rng=rng,
        data=data,
        config=config,
        requested_shots=extra_shots,
        budget=budget,
        settings=settings,
        reserve_hqc=reserve_hqc,
        diagnostics=diagnostics,
        stage=stage,
    )
    confirmed_probability = measured_probability(data, config)
    return probability if confirmed_probability is None else confirmed_probability


def best_fixed_ratio_bracket(
    measured_points: list[tuple[RMBConfig, float]],
) -> tuple[RMBConfig, float, RMBConfig, float] | None:
    points = sorted(
        measured_points,
        key=lambda item: float(item[0].depth),
    )
    best: tuple[RMBConfig, float, RMBConfig, float] | None = None
    best_width = np.inf
    for (low_config, low_p), (high_config, high_p) in zip(points, points[1:]):
        if low_p >= 0.5 and high_p <= 0.5:
            width = float(high_config.depth) - float(low_config.depth)
            if width < best_width:
                best_width = width
                best = (low_config, low_p, high_config, high_p)
    return best


def find_depth_crossing_at_ratio(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    template: RMBConfig,
    ratio: float,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    center_depth: float | None = None,
    local: bool = False,
    reserve_hqc: float = 0.0,
    diagnostics: list[dict] | None = None,
    stage: str = "crossing",
) -> RMBConfig | None:
    """
    Find a monotone crossing in depth at fixed ratio.
    """
    d_min, d_max = settings.depth_bounds
    if center_depth is None or not local:
        if settings.initial_crossing_interior_bracket and not local:
            depth_span = d_max - d_min
            center = d_min + settings.initial_crossing_center_fraction * depth_span
            half_width = settings.initial_crossing_half_width_fraction * depth_span
            low_depth = max(float(d_min), center - half_width)
            high_depth = min(float(d_max), center + half_width)
        else:
            low_depth = float(d_min)
            high_depth = float(d_max)
    else:
        span = settings.trace_depth_search_fraction * (d_max - d_min)
        low_depth = max(float(d_min), center_depth - span)
        high_depth = min(float(d_max), center_depth + span)

    probe_shots = settings.ray_probe_shots if not local else settings.trace_local_stencil_shots
    if local and settings.trace_local_stencil_points > 2:
        stencil_depths = np.linspace(
            low_depth,
            high_depth,
            max(2, settings.trace_local_stencil_points),
        )
        stencil_results = probe_depths_batch(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            probes=[
                (float(depth), ratio, probe_shots)
                for depth in stencil_depths
            ],
            budget=budget,
            settings=settings,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage=f"{stage}_stencil",
        )
        measured_stencil = [
            (config, p)
            for config, p in stencil_results
            if p is not None
        ]
        bracket = best_fixed_ratio_bracket(measured_stencil)
        if bracket is not None:
            low_config, low_p, high_config, high_p = bracket
            low_depth = float(low_config.depth)
            high_depth = float(high_config.depth)
        elif measured_stencil:
            low_config, low_p = min(
                measured_stencil,
                key=lambda item: float(item[0].depth),
            )
            high_config, high_p = max(
                measured_stencil,
                key=lambda item: float(item[0].depth),
            )
            low_depth = float(low_config.depth)
            high_depth = float(high_config.depth)
        else:
            low_config = high_config = None
            low_p = high_p = None
        crossing_candidates: list[tuple[RMBConfig, float]] = measured_stencil
    else:
        (low_config, low_p), (high_config, high_p) = probe_depths_batch(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            probes=[
                (low_depth, ratio, probe_shots),
                (high_depth, ratio, probe_shots),
            ],
            budget=budget,
            settings=settings,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage=f"{stage}_bracket",
        )
        crossing_candidates = [
            (config, p)
            for config, p in [(low_config, low_p), (high_config, high_p)]
            if p is not None
        ]

    if low_config is None or high_config is None or low_p is None or high_p is None:
        add_diagnostic(
            diagnostics,
            "crossing_failed",
            stage=stage,
            ratio=round(float(ratio), 4),
            local=local,
            reason="missing_probe",
        )
        return None

    if (
        not (low_p >= 0.5 and high_p <= 0.5)
        and (local or settings.initial_crossing_interior_bracket)
    ):
        expand = (
            settings.trace_depth_search_fraction * (d_max - d_min)
            if local
            else settings.initial_crossing_half_width_fraction
            * settings.initial_crossing_expand_factor
            * (d_max - d_min)
        )
        for _ in range(3):
            if not budget.can_spend(reserve_hqc=reserve_hqc):
                return None
            expansion_probes: list[tuple[str, float, float, int]] = []
            if low_p is not None and low_p < 0.5:
                low_depth = max(float(d_min), low_depth - expand)
                expansion_probes.append(("low", low_depth, ratio, settings.trace_shots))
            if high_p is not None and high_p > 0.5:
                high_depth = min(float(d_max), high_depth + expand)
                expansion_probes.append(("high", high_depth, ratio, settings.trace_shots))
            if expansion_probes:
                expansion_results = probe_depths_batch(
                    backend=backend,
                    rng=rng,
                    data=data,
                    template=template,
                    probes=[
                        (depth, probe_ratio, shots)
                        for _, depth, probe_ratio, shots in expansion_probes
                    ],
                    budget=budget,
                    settings=settings,
                    reserve_hqc=reserve_hqc,
                    diagnostics=diagnostics,
                    stage=f"{stage}_expand",
                )
                for (side, _, _, _), (expanded_config, expanded_p) in zip(
                    expansion_probes,
                    expansion_results,
                ):
                    if side == "low":
                        low_config, low_p = expanded_config, expanded_p
                    else:
                        high_config, high_p = expanded_config, expanded_p
                    if expanded_p is not None:
                        crossing_candidates.append((expanded_config, expanded_p))
            if low_p is not None and high_p is not None and low_p >= 0.5 and high_p <= 0.5:
                break
            if low_depth <= d_min and high_depth >= d_max:
                break
            expand *= 1.5

    if low_p is None or high_p is None or not (low_p >= 0.5 and high_p <= 0.5):
        add_diagnostic(
            diagnostics,
            "crossing_failed",
            stage=stage,
            ratio=round(float(ratio), 4),
            local=local,
            low_depth=round(float(low_depth), 4),
            low_p=None if low_p is None else round(float(low_p), 4),
            high_depth=round(float(high_depth), 4),
            high_p=None if high_p is None else round(float(high_p), 4),
            reason="no_bracket",
        )
        return None

    add_diagnostic(
        diagnostics,
        "crossing_bracketed",
        stage=stage,
        ratio=round(float(ratio), 4),
        local=local,
        low_depth=round(float(low_depth), 4),
        low_p=round(float(low_p), 4),
        high_depth=round(float(high_depth), 4),
        high_p=round(float(high_p), 4),
    )

    best_config = low_config if abs(low_p - 0.5) <= abs(high_p - 0.5) else high_config
    best_error = min(abs(low_p - 0.5), abs(high_p - 0.5))
    n_steps = settings.trace_correction_steps if local else settings.ray_bisection_steps
    shots = settings.trace_shots if local else settings.ray_bisection_shots

    for _ in range(n_steps):
        if not budget.can_spend(reserve_hqc=reserve_hqc):
            break
        mid_depth = 0.5 * (low_depth + high_depth)
        mid_config, mid_p = probe_depth(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            depth=mid_depth,
            ratio=ratio,
            shots=shots,
            budget=budget,
            settings=settings,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage=f"{stage}_bisect",
        )
        if mid_p is None:
            break
        mid_p = confirm_crossing_decision_if_needed(
            backend=backend,
            rng=rng,
            data=data,
            config=mid_config,
            probability=mid_p,
            settings=settings,
            budget=budget,
            local=local,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage=f"{stage}_confirm",
        )
        crossing_candidates.append((mid_config, mid_p))
        error = abs(mid_p - 0.5)
        if error < best_error:
            best_error = error
            best_config = mid_config
        ambiguous = abs(mid_p - 0.5) <= settings.crossing_decision_probability_width
        if mid_p >= 0.5 and not (
            not local and settings.crossing_ambiguous_as_failure and ambiguous
        ):
            low_depth = float(mid_config.depth)
            low_config = mid_config
            low_p = mid_p
        else:
            high_depth = float(mid_config.depth)
            high_config = mid_config
            high_p = mid_p

    if (
        not local
        and settings.crossing_confirm_candidates > 0
        and settings.crossing_confirm_shots > 0
    ):
        unique_candidates = {}
        for config, probability in crossing_candidates:
            unique_candidates[config] = probability
        ranked_candidates = sorted(
            unique_candidates.items(),
            key=lambda item: abs(item[1] - 0.5),
        )[: settings.crossing_confirm_candidates]

        for candidate_config, probability in ranked_candidates:
            current_runs = data[candidate_config].num_runs() if candidate_config in data else 0
            extra_shots = settings.crossing_confirm_shots - current_runs
            if extra_shots > 0:
                spend_config(
                    backend=backend,
                    rng=rng,
                    data=data,
                    config=candidate_config,
                    requested_shots=extra_shots,
                    budget=budget,
                    settings=settings,
                    reserve_hqc=reserve_hqc,
                )
            confirmed_probability = measured_probability(data, candidate_config)
            if confirmed_probability is None:
                confirmed_probability = probability
            error = abs(confirmed_probability - 0.5)
            if error < best_error:
                best_error = error
                best_config = candidate_config

    at_depth_bound = (
        best_config.depth <= d_min + 1
        or best_config.depth >= d_max - 1
    )
    if at_depth_bound and best_error > settings.trace_accept_probability_width:
        add_diagnostic(
            diagnostics,
            "crossing_rejected",
            stage=stage,
            ratio=round(float(ratio), 4),
            local=local,
            depth=int(best_config.depth),
            p=round(float(measured_probability(data, best_config)), 4),
            std=round(measured_std(data, best_config), 4),
            runs=data[best_config].num_runs(),
            reason="depth_bound",
        )
        return None
    if best_error <= settings.trace_reject_probability_width:
        add_diagnostic(
            diagnostics,
            "crossing_accepted",
            stage=stage,
            ratio=round(float(ratio), 4),
            local=local,
            depth=int(best_config.depth),
            p=round(float(measured_probability(data, best_config)), 4),
            std=round(measured_std(data, best_config), 4),
            runs=data[best_config].num_runs(),
        )
        return best_config
    add_diagnostic(
        diagnostics,
        "crossing_rejected",
        stage=stage,
        ratio=round(float(ratio), 4),
        local=local,
        depth=int(best_config.depth),
        p=round(float(measured_probability(data, best_config)), 4),
        std=round(measured_std(data, best_config), 4),
        runs=data[best_config].num_runs(),
        reason="outside_reject_width",
    )
    return None


def find_initial_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    n_qubits: int,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> RMBConfig | None:
    template = template_config(settings, n_qubits)
    reserve_hqc = trace_reserve_hqc(settings)
    ratios = np.linspace(
        settings.ratio_bounds[0],
        settings.ratio_bounds[1],
        max(1, settings.ray_ratio_count),
    )
    if settings.batched_initial_anchor_search and budget.can_spend(reserve_hqc=reserve_hqc):
        d_min, d_max = settings.depth_bounds
        ratio = float(ratios[0])
        if settings.initial_crossing_interior_bracket:
            depth_span = d_max - d_min
            center = d_min + settings.initial_crossing_center_fraction * depth_span
            half_width = settings.initial_crossing_half_width_fraction * depth_span
            low_depth = max(float(d_min), center - half_width)
            high_depth = min(float(d_max), center + half_width)
        else:
            low_depth = float(d_min)
            high_depth = float(d_max)

        depth_grid = np.linspace(
            low_depth,
            high_depth,
            max(2, settings.initial_anchor_depth_grid_count),
        )
        grid_configs = [
            config_from_parameters(
                template=template,
                depth=float(depth),
                ratio=ratio,
            )
            for depth in depth_grid
        ]
        requests = [
            MeasurementRequest(config, settings.initial_anchor_grid_shots)
            for config in grid_configs
        ]

        spend_configs_batch(
            backend=backend,
            rng=rng,
            data=data,
            requests=requests,
            budget=budget,
            settings=settings,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage="initial_grid",
        )

        measured_grid: list[tuple[RMBConfig, float]] = []
        for config in grid_configs:
            p = measured_probability(data, config)
            if p is None:
                continue
            measured_grid.append((config, p))
            add_diagnostic(
                diagnostics,
                "initial_grid_probe",
                ratio=round(float(config.min_two_qubit_gate_ratio), 4),
                depth=int(config.depth),
                p=round(float(p), 4),
            )

        anchor_config: RMBConfig | None = None
        anchor_error = np.inf
        measured_grid.sort(key=lambda item: float(item[0].depth))
        for (low_config, low_p), (high_config, high_p) in zip(measured_grid, measured_grid[1:]):
            if low_p >= 0.5 and high_p <= 0.5:
                low_error = abs(low_p - 0.5)
                high_error = abs(high_p - 0.5)
                candidate = low_config if low_error <= high_error else high_config
                error = min(low_error, high_error)
                if error < anchor_error:
                    anchor_config = candidate
                    anchor_error = error
        if anchor_config is not None:
            add_diagnostic(
                diagnostics,
                "initial_grid_anchor",
                depth=int(anchor_config.depth),
                ratio=round(float(anchor_config.min_two_qubit_gate_ratio), 4),
                p=round(float(measured_probability(data, anchor_config)), 4),
                error=round(float(anchor_error), 4),
                selection="fixed_ratio_depth_grid",
            )
            return anchor_config

    for ratio in ratios:
        if not budget.can_spend(reserve_hqc=reserve_hqc):
            return None
        anchor = find_depth_crossing_at_ratio(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            ratio=float(ratio),
            settings=settings,
            budget=budget,
            local=False,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage="initial",
        )
        if anchor is not None:
            return anchor
    return None


def confirm_initial_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchor: RMBConfig,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> RMBConfig:
    reserve_hqc = trace_reserve_hqc(settings)
    if not settings.initial_anchor_refine or not budget.can_spend(reserve_hqc=reserve_hqc):
        return anchor

    template = template_config(settings, anchor.n_qubits)
    confirm_settings = replace(
        settings,
        trace_shots=settings.initial_anchor_refine_shots,
        trace_correction_steps=settings.initial_anchor_refine_steps,
        trace_depth_search_fraction=settings.initial_anchor_search_fraction,
    )
    refined = find_depth_crossing_at_ratio(
        backend=backend,
        rng=rng,
        data=data,
        template=template,
        ratio=float(anchor.min_two_qubit_gate_ratio),
        settings=confirm_settings,
        budget=budget,
        center_depth=float(anchor.depth),
        local=True,
        reserve_hqc=reserve_hqc,
        diagnostics=diagnostics,
        stage="initial_refine",
    )
    confirmed_anchor = refined if refined is not None else anchor

    if confirmed_anchor in data:
        current_runs = data[confirmed_anchor].num_runs()
        extra_shots = settings.initial_anchor_min_runs - current_runs
        if extra_shots > 0:
            spend_config(
                backend=backend,
                rng=rng,
                data=data,
                config=confirmed_anchor,
                requested_shots=extra_shots,
                budget=budget,
                settings=settings,
                reserve_hqc=reserve_hqc,
            )

    add_diagnostic(
        diagnostics,
        "initial_anchor_confirmed",
        depth=int(confirmed_anchor.depth),
        ratio=round(float(confirmed_anchor.min_two_qubit_gate_ratio), 4),
        p=None if measured_probability(data, confirmed_anchor) is None else round(
            float(measured_probability(data, confirmed_anchor)), 4),
        std=round(measured_std(data, confirmed_anchor), 4),
        runs=data[confirmed_anchor].num_runs() if confirmed_anchor in data else 0,
    )
    return confirmed_anchor


def projected_trace_target(
    data: RMBData,
    anchor: RMBConfig,
    direction: int,
    settings: ContourFirstExperimentConfig,
    step_fraction: float | None = None,
) -> tuple[float, float]:
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    if step_fraction is None:
        step_fraction = settings.trace_ratio_step_fraction
    ratio_step = direction * step_fraction * ratio_span
    next_ratio = float(np.clip(
        anchor.min_two_qubit_gate_ratio + ratio_step,
        settings.ratio_bounds[0],
        settings.ratio_bounds[1],
    ))
    next_depth = float(anchor.depth)

    if settings.model_projection_after_fit:
        try:
            surface = fit_monotone_fidelity_surface(data, settings)
        except (RuntimeError, ValueError):
            return next_depth, next_ratio
        point = np.array([float(anchor.depth), float(anchor.min_two_qubit_gate_ratio)])
        gradient = finite_difference_gradient(surface, point, settings)
        tangent = np.array([gradient[1], -gradient[0]], dtype=float)
        if np.linalg.norm(tangent) > 1e-12:
            if np.sign(tangent[1]) != np.sign(direction):
                tangent *= -1.0
            depth_span = settings.depth_bounds[1] - settings.depth_bounds[0]
            scaled_step = np.array([depth_span, ratio_span], dtype=float)
            tangent_scaled = tangent * scaled_step
            tangent_scaled /= np.linalg.norm(tangent_scaled)
            trial = point + tangent_scaled * step_fraction * scaled_step
            trial[1] = next_ratio
            trial[0] = np.clip(trial[0], settings.depth_bounds[0], settings.depth_bounds[1])
            projected = project_to_contour(surface, trial, gradient, settings)
            next_depth = float(projected[0])

    if settings.trace_enforce_monotone_depth:
        if direction > 0:
            next_depth = min(next_depth, float(anchor.depth))
        elif direction < 0:
            next_depth = max(next_depth, float(anchor.depth))

    return next_depth, next_ratio


def trace_projected_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    current: RMBConfig,
    depth_guess: float,
    ratio: float,
    direction: int,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    reserve_hqc: float = 0.0,
    diagnostics: list[dict] | None = None,
) -> RMBConfig | None:
    """
    Measure the projected contour point.  Optionally solve a local fixed-ratio
    depth crossing, but keep the cheaper projected trace as the default because
    local bracketing can spend too much budget at low HQC.
    """
    template = template_config(settings, current.n_qubits)
    d_min, d_max = settings.depth_bounds
    lower_depth = float(d_min)
    upper_depth = float(d_max)
    if settings.trace_enforce_monotone_depth:
        if direction > 0:
            upper_depth = min(upper_depth, float(current.depth))
        elif direction < 0:
            lower_depth = max(lower_depth, float(current.depth))

    if lower_depth > upper_depth:
        add_diagnostic(
            diagnostics,
            "trace_candidate_rejected",
            reason="empty_depth_bounds",
            current_depth=int(current.depth),
            current_ratio=round(float(current.min_two_qubit_gate_ratio), 4),
            proposed_ratio=round(float(ratio), 4),
            depth_guess=round(float(depth_guess), 4),
        )
        return None

    depth = float(np.clip(depth_guess, lower_depth, upper_depth))
    projected_config = config_from_parameters(template=template, depth=depth, ratio=ratio)
    use_projected_trace = (
        settings.high_ratio_projected_trace
        and settings.trace_local_crossing
        and direction > 0
        and high_ratio_region(projected_config, settings)
    )
    if settings.trace_local_crossing and not use_projected_trace:
        return find_depth_crossing_at_ratio(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            ratio=float(ratio),
            budget=budget,
            settings=settings,
            center_depth=depth,
            local=True,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage="trace_local",
        )
    if use_projected_trace:
        add_diagnostic(
            diagnostics,
            "trace_high_ratio_projected",
            depth=int(projected_config.depth),
            ratio=round(float(projected_config.min_two_qubit_gate_ratio), 4),
            reason="high_ratio_low_depth",
        )

    best_config, best_p = probe_depth(
        backend=backend,
        rng=rng,
        data=data,
        template=template,
        depth=depth,
        ratio=ratio,
        shots=settings.trace_shots,
        budget=budget,
        settings=settings,
        reserve_hqc=reserve_hqc,
        diagnostics=diagnostics,
        stage="trace_projected_probe",
    )
    if best_p is None:
        add_diagnostic(
            diagnostics,
            "trace_candidate_rejected",
            reason="no_measurement",
            proposed_depth=round(float(depth), 4),
            proposed_ratio=round(float(ratio), 4),
        )
        return None

    best_error = abs(best_p - 0.5)
    if best_error <= settings.trace_accept_probability_width:
        add_diagnostic(
            diagnostics,
            "trace_candidate_accepted",
            reason="projected_accept",
            depth=int(best_config.depth),
            ratio=round(float(best_config.min_two_qubit_gate_ratio), 4),
            p=round(float(best_p), 4),
            std=round(measured_std(data, best_config), 4),
            runs=data[best_config].num_runs(),
        )
        return best_config

    correction_step = settings.trace_depth_search_fraction * (d_max - d_min)
    for _ in range(settings.trace_correction_steps):
        if not budget.can_spend(reserve_hqc=reserve_hqc):
            break
        if best_p >= 0.5:
            depth = min(upper_depth, float(best_config.depth) + correction_step)
        else:
            depth = max(lower_depth, float(best_config.depth) - correction_step)
        if abs(depth - best_config.depth) < 1e-9:
            break

        candidate, p = probe_depth(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            depth=depth,
            ratio=ratio,
            shots=settings.trace_shots,
            budget=budget,
            settings=settings,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage="trace_projected_correction",
        )
        if p is None:
            break
        error = abs(p - 0.5)
        if error < best_error:
            best_config = candidate
            best_p = p
            best_error = error
        if error <= settings.trace_accept_probability_width:
            add_diagnostic(
                diagnostics,
                "trace_candidate_accepted",
                reason="correction_accept",
                depth=int(candidate.depth),
                ratio=round(float(candidate.min_two_qubit_gate_ratio), 4),
                p=round(float(p), 4),
                std=round(measured_std(data, candidate), 4),
                runs=data[candidate].num_runs(),
            )
            return candidate
        correction_step *= 0.5

    add_diagnostic(
        diagnostics,
        "trace_candidate_accepted",
        reason="best_available",
        depth=int(best_config.depth),
        ratio=round(float(best_config.min_two_qubit_gate_ratio), 4),
        p=round(float(best_p), 4),
        std=round(measured_std(data, best_config), 4),
        runs=data[best_config].num_runs(),
        error=round(float(best_error), 4),
    )
    return best_config


def trace_from_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchor: RMBConfig,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> list[RMBConfig]:
    anchors = [anchor]
    reserve_hqc = post_trace_reserve_hqc(settings)
    for direction in settings.trace_directions:
        current = anchor
        while budget.can_spend(reserve_hqc=reserve_hqc):
            next_anchor = None
            n_attempts = max(1, settings.trace_step_shrink_attempts if settings.trace_local_crossing else 1)
            for attempt in range(n_attempts):
                step_fraction = settings.trace_ratio_step_fraction * (0.5 ** attempt)
                depth_guess, ratio = projected_trace_target(
                    data,
                    current,
                    direction,
                    settings,
                    step_fraction=step_fraction,
                )
                if abs(ratio - current.min_two_qubit_gate_ratio) < 1e-9:
                    add_diagnostic(
                        diagnostics,
                        "trace_stopped",
                        reason="ratio_bounds",
                        direction=direction,
                        current_depth=int(current.depth),
                        current_ratio=round(float(current.min_two_qubit_gate_ratio), 4),
                    )
                    break
                add_diagnostic(
                    diagnostics,
                    "trace_proposed",
                    direction=direction,
                    attempt=attempt,
                    step_fraction=round(float(step_fraction), 4),
                    current_depth=int(current.depth),
                    current_ratio=round(float(current.min_two_qubit_gate_ratio), 4),
                    depth_guess=round(float(depth_guess), 4),
                    proposed_ratio=round(float(ratio), 4),
                )
                next_anchor = trace_projected_anchor(
                    backend=backend,
                    rng=rng,
                    data=data,
                    current=current,
                    depth_guess=depth_guess,
                    ratio=float(ratio),
                    direction=direction,
                    settings=settings,
                    budget=budget,
                    reserve_hqc=reserve_hqc,
                    diagnostics=diagnostics,
                )
                if next_anchor is not None:
                    break
            if next_anchor is None:
                add_diagnostic(
                    diagnostics,
                    "trace_stopped",
                    reason="candidate_failed",
                    direction=direction,
                    current_depth=int(current.depth),
                    current_ratio=round(float(current.min_two_qubit_gate_ratio), 4),
                )
                break
            add_diagnostic(
                diagnostics,
                "trace_accepted",
                direction=direction,
                depth=int(next_anchor.depth),
                ratio=round(float(next_anchor.min_two_qubit_gate_ratio), 4),
                p=None if measured_probability(data, next_anchor) is None else round(
                    float(measured_probability(data, next_anchor)), 4),
                std=round(measured_std(data, next_anchor), 4),
                runs=data[next_anchor].num_runs() if next_anchor in data else 0,
            )
            anchors.append(next_anchor)
            current = next_anchor
    return anchors


def seeded_batched_trace_from_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchor: RMBConfig,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> list[RMBConfig]:
    """
    Seed the contour with a small number of same-ratio depth stacks.

    A rough inverse-depth prior lets us exploit batching early: instead of
    adaptively spending one batch per ratio while discovering the contour, we
    measure several likely contour ratios in one or two stitched batches.
    """
    if (
        not settings.seeded_initial_trace
        or settings.seeded_trace_ratio_count <= 0
        or not budget.can_spend(reserve_hqc=post_trace_reserve_hqc(settings))
    ):
        return []

    d_min, d_max = settings.depth_bounds
    ratio_min, ratio_max = settings.ratio_bounds
    current_ratio = float(anchor.min_two_qubit_gate_ratio)
    ratio_span = ratio_max - ratio_min
    upper_ratio = min(
        float(ratio_max),
        current_ratio + max(0.0, settings.seeded_trace_ratio_span_fraction) * ratio_span,
    )
    if upper_ratio <= current_ratio + 1e-9:
        return []

    target_ratios = np.linspace(
        current_ratio,
        upper_ratio,
        settings.seeded_trace_ratio_count + 1,
    )[1:]
    template = template_config(settings, anchor.n_qubits)
    depth_power = max(0.0, settings.seeded_trace_depth_power)
    request_seen: set[RMBConfig] = set()
    requests: list[MeasurementRequest] = []
    predicted_configs: list[RMBConfig] = []
    last_depth = float(anchor.depth)

    for ratio in target_ratios:
        predicted_depth = float(anchor.depth) * (
            max(current_ratio, 1e-6) / max(float(ratio), 1e-6)
        ) ** depth_power
        if settings.trace_enforce_monotone_depth:
            predicted_depth = min(predicted_depth, last_depth)
        predicted_depth = float(np.clip(predicted_depth, d_min, d_max))
        predicted_config = config_from_parameters(
            template=template,
            depth=predicted_depth,
            ratio=float(ratio),
        )
        predicted_configs.append(predicted_config)
        requests.extend(
            contour_bracket_requests(
                config=predicted_config,
                data=data,
                settings=settings,
                template=template,
                seen=request_seen,
            )
        )
        last_depth = predicted_depth

    if not requests:
        return []

    max_batches = max(1, settings.seeded_trace_max_batches)
    chunk_size = max(1, int(np.ceil(len(requests) / max_batches)))
    reserve_hqc = post_trace_reserve_hqc(settings)
    before = total_measurements(data)
    for start in range(0, len(requests), chunk_size):
        if not budget.can_spend(reserve_hqc=reserve_hqc):
            break
        spend_configs_batch(
            backend=backend,
            rng=rng,
            data=data,
            requests=requests[start:start + chunk_size],
            budget=budget,
            settings=settings,
            reserve_hqc=reserve_hqc,
            diagnostics=diagnostics,
            stage="seeded_trace",
        )
    if total_measurements(data) == before:
        return []

    anchors: list[RMBConfig] = []
    used_ratios: set[float] = set()
    for predicted_config in predicted_configs:
        ratio = float(predicted_config.min_two_qubit_gate_ratio)
        ratio_key = round(ratio, 8)
        candidates = [
            config
            for config in data
            if config.n_qubits == predicted_config.n_qubits
            and round(float(config.min_two_qubit_gate_ratio), 8) == ratio_key
            and data[config].num_runs() > 0
        ]
        if not candidates or ratio_key in used_ratios:
            continue
        probabilities = [fidelity_mean(data[config]) for config in candidates]
        if not (any(p >= 0.5 for p in probabilities) and any(p <= 0.5 for p in probabilities)):
            continue
        best_config = min(
            candidates,
            key=lambda config: abs(fidelity_mean(data[config]) - 0.5),
        )
        best_p = fidelity_mean(data[best_config])
        if abs(best_p - 0.5) > settings.trace_reject_probability_width:
            continue
        if anchors and settings.trace_enforce_monotone_depth:
            if float(best_config.depth) > float(anchors[-1].depth):
                continue
        anchors.append(best_config)
        used_ratios.add(ratio_key)
        add_diagnostic(
            diagnostics,
            "seeded_trace_anchor",
            depth=int(best_config.depth),
            ratio=round(float(best_config.min_two_qubit_gate_ratio), 4),
            p=round(float(best_p), 4),
            std=round(measured_std(data, best_config), 4),
            runs=data[best_config].num_runs(),
        )

    return anchors


def refinement_score(
    config: RMBConfig,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
    surface=None,
) -> float:
    if config not in data:
        return 0.0
    estimator = data[config]
    if estimator.num_runs() >= settings.max_shots_per_config:
        return 0.0

    variance = fidelity_variance(estimator)
    variance_score = min(1.0, variance / (1.0 / 12.0))
    if surface is None:
        p_model = fidelity_mean(estimator)
    else:
        p_model = float(surface.probability(np.array([[
            float(config.depth),
            float(config.min_two_qubit_gate_ratio),
        ]]))[0])

    boundary_score = np.exp(-((abs(p_model - 0.5) / settings.refinement_boundary_width) ** 2))
    return float(variance_score * boundary_score)


def refinement_request_efficiency(
    config: RMBConfig,
    score: float,
    shots: int,
    settings: ContourFirstExperimentConfig,
) -> float:
    request = MeasurementRequest(config, shots)
    marginal_cost = request_marginal_hqc_cost(request, settings)
    return score / (marginal_cost ** max(0.0, settings.hqc_cost_power))


def batched_topup_shots(
    config: RMBConfig,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
    target_runs: int | None = None,
) -> int:
    if config not in data:
        current_runs = 0
    else:
        current_runs = data[config].num_runs()
    if target_runs is None:
        target_runs = settings.max_shots_per_config
    missing = max(0, min(target_runs, settings.max_shots_per_config) - current_runs)
    if settings.batching_enabled:
        increment = max(settings.refinement_shots, settings.batch_refinement_shots)
    else:
        increment = settings.refinement_shots
    return min(missing, max(1, increment))


def contour_bracket_depths(
    *,
    depth: float,
    settings: ContourFirstExperimentConfig,
) -> list[float]:
    d_min, d_max = settings.depth_bounds
    depth_span = max(1.0, float(d_max - d_min))
    candidate_depths = [float(depth)]
    for fraction in settings.contour_bracket_depth_fractions:
        span_offset = abs(float(fraction)) * depth_span
        relative_offset = (
            max(0.0, settings.contour_bracket_max_relative_depth)
            * max(1.0, float(depth))
        )
        offset = max(1.0, min(span_offset, max(1.0, relative_offset)))
        candidate_depths.extend([
            float(depth) - offset,
            float(depth) + offset,
        ])
    return [
        candidate_depth
        for candidate_depth in candidate_depths
        if float(d_min) <= candidate_depth <= float(d_max)
    ]


def contour_bracket_requests(
    *,
    config: RMBConfig,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
    template: RMBConfig | None = None,
    seen: set[RMBConfig] | None = None,
) -> list[MeasurementRequest]:
    """
    Request a same-ratio depth stack around a contour candidate.

    The monotone fit is much better constrained when each ratio has measured
    depths on both sides of the 0.5 crossing, not just an isolated near-contour
    point.
    """
    if seen is None:
        seen = set()
    if template is None:
        template = template_config(settings, config.n_qubits)

    ratio = float(config.min_two_qubit_gate_ratio)
    candidate_depths = contour_bracket_depths(
        depth=float(config.depth),
        settings=settings,
    )

    requests: list[MeasurementRequest] = []
    for depth in candidate_depths:
        bracket_config = config_from_parameters(
            template=template,
            depth=depth,
            ratio=ratio,
        )
        if bracket_config in seen:
            continue
        seen.add(bracket_config)
        current_runs = data[bracket_config].num_runs() if bracket_config in data else 0
        missing = max(0, settings.max_shots_per_config - current_runs)
        shots = min(max(1, settings.contour_bracket_probe_shots), missing)
        if shots > 0:
            requests.append(MeasurementRequest(bracket_config, shots))
    return requests


def fill_requests_toward_batch_cost(
    requests: list[MeasurementRequest],
    *,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
    spendable_hqc: float,
    target_runs: int | None = None,
) -> list[MeasurementRequest]:
    """
    Increase repeats on selected configs until the stitched batch is near full.
    """
    if not settings.batch_fill_repeats or not requests:
        return requests
    max_target_cost = target_batch_cost(settings, spendable_hqc)
    if max_target_cost <= 0.0:
        return requests

    filled = [
        MeasurementRequest(request.config, max(0, request.requested_shots))
        for request in requests
        if request.requested_shots > 0
    ]
    if not filled:
        return []

    while True:
        current_cost = batched_hqc_cost(filled, settings)
        if current_cost >= max_target_cost:
            break
        best_index: int | None = None
        best_score = -np.inf
        for idx, request in enumerate(filled):
            config = request.config
            current_runs = data[config].num_runs() if config in data else 0
            default_target = (
                settings.max_shots_per_config
                if settings.batch_fill_max_shots_per_config is None
                else max(settings.max_shots_per_config, settings.batch_fill_max_shots_per_config)
            )
            max_runs = default_target if target_runs is None else max(1, target_runs)
            if current_runs + request.requested_shots >= max_runs:
                continue
            trial = list(filled)
            trial[idx] = MeasurementRequest(config, request.requested_shots + 1)
            trial_cost = batched_hqc_cost(trial, settings)
            if trial_cost > spendable_hqc + 1e-9 or trial_cost > max_target_cost + 1e-9:
                continue
            marginal_cost = max(1e-12, trial_cost - current_cost)
            variance = fidelity_variance(data[config]) if config in data else 1.0 / 12.0
            score = variance / (marginal_cost ** max(0.0, settings.hqc_cost_power))
            if score > best_score:
                best_score = score
                best_index = idx
        if best_index is None:
            break
        request = filled[best_index]
        filled[best_index] = MeasurementRequest(request.config, request.requested_shots + 1)

    return filled


def scaled_depth_ratio_point(
    *,
    depth: float,
    ratio: float,
    settings: ContourFirstExperimentConfig,
) -> np.ndarray:
    depth_span = max(1e-12, settings.depth_bounds[1] - settings.depth_bounds[0])
    ratio_span = max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0])
    return np.array([
        (depth - settings.depth_bounds[0]) / depth_span,
        (ratio - settings.ratio_bounds[0]) / ratio_span,
    ])


def high_ratio_region(config: RMBConfig, settings: ContourFirstExperimentConfig) -> bool:
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    depth_span = settings.depth_bounds[1] - settings.depth_bounds[0]
    ratio_threshold = settings.ratio_bounds[0] + settings.high_ratio_projected_threshold * ratio_span
    depth_threshold = settings.depth_bounds[0] + settings.high_ratio_low_depth_fraction * depth_span
    return (
        float(config.min_two_qubit_gate_ratio) >= ratio_threshold
        and float(config.depth) <= depth_threshold
    )


def high_ratio_acquisition_bonus(config: RMBConfig, settings: ContourFirstExperimentConfig) -> float:
    if settings.high_ratio_acquisition_fraction <= 0.0:
        return 1.0
    ratio_span = max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0])
    depth_span = max(1e-12, settings.depth_bounds[1] - settings.depth_bounds[0])
    ratio_position = (float(config.min_two_qubit_gate_ratio) - settings.ratio_bounds[0]) / ratio_span
    low_depth_position = (settings.depth_bounds[1] - float(config.depth)) / depth_span
    shape_score = np.clip(ratio_position * low_depth_position, 0.0, 1.0)
    return float(1.0 + settings.high_ratio_acquisition_fraction * shape_score)


def acquisition_bracket_straddles_surface(
    *,
    config: RMBConfig,
    surface,
    settings: ContourFirstExperimentConfig,
) -> bool:
    if not settings.acquisition_require_bracket_straddle:
        return True
    ratio = float(config.min_two_qubit_gate_ratio)
    depths = contour_bracket_depths(
        depth=float(config.depth),
        settings=settings,
    )
    if not depths:
        return False
    points = np.array([[depth, ratio] for depth in depths], dtype=float)
    probabilities = surface.probability(points)
    return bool(
        float(np.min(probabilities)) <= 0.5 <= float(np.max(probabilities))
    )


def acquisition_trace_window(
    anchors: list[RMBConfig] | None,
    *,
    pass_index: int,
    settings: ContourFirstExperimentConfig,
) -> tuple[float, float] | None:
    if not settings.acquisition_follow_trace or not anchors:
        return None
    ratios = [
        float(anchor.min_two_qubit_gate_ratio)
        for anchor in anchors
        if settings.ratio_bounds[0] <= float(anchor.min_two_qubit_gate_ratio) <= settings.ratio_bounds[1]
    ]
    if not ratios:
        return None

    ratio_span = max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0])
    last_ratio = max(ratios)
    backtrack = max(0.0, settings.acquisition_trace_backtrack_fraction) * ratio_span
    extension = max(0.0, settings.acquisition_trace_extension_fraction) * ratio_span
    lower = max(settings.ratio_bounds[0], last_ratio - backtrack)
    upper = min(settings.ratio_bounds[1], last_ratio + (pass_index + 1) * extension)
    if upper < lower:
        return None
    return lower, upper


def in_acquisition_trace_window(
    config: RMBConfig,
    window: tuple[float, float] | None,
) -> bool:
    if window is None:
        return True
    ratio = float(config.min_two_qubit_gate_ratio)
    return window[0] <= ratio <= window[1]


def trace_follow_acquisition_candidates(
    anchors: list[RMBConfig] | None,
    *,
    surface,
    pass_index: int,
    n_qubits_values: list[int],
    settings: ContourFirstExperimentConfig,
    candidate_count: int,
) -> list[RMBConfig]:
    if (
        not settings.acquisition_follow_trace
        or not anchors
        or candidate_count <= 0
        or settings.acquisition_trace_predict_points <= 0
    ):
        return []

    sorted_anchors = sorted(
        {
            (
                anchor.n_qubits,
                float(anchor.min_two_qubit_gate_ratio),
                float(anchor.depth),
            )
            for anchor in anchors
        },
        key=lambda item: item[1],
    )
    if not sorted_anchors:
        return []

    ratios = np.asarray([item[1] for item in sorted_anchors], dtype=float)
    depths = np.asarray([item[2] for item in sorted_anchors], dtype=float)
    ratio_span = max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0])
    last_ratio = float(ratios[-1])
    if last_ratio >= settings.ratio_bounds[1]:
        return []

    fit_points = min(
        max(2, settings.acquisition_trace_predict_points),
        len(ratios),
    )
    recent_ratios = ratios[-fit_points:]
    recent_depths = depths[-fit_points:]
    if len(recent_ratios) >= 2 and float(np.ptp(recent_ratios)) > 1e-12:
        slope, intercept = np.polyfit(recent_ratios, recent_depths, deg=1)
    else:
        slope = -(
            settings.depth_bounds[1] - settings.depth_bounds[0]
        ) / ratio_span
        intercept = float(recent_depths[-1]) - slope * last_ratio
    slope = min(0.0, float(slope))

    step = max(
        settings.trace_ratio_step_fraction * ratio_span,
        ratio_span / max(2, settings.candidate_grid_size[1] - 1),
    )
    start_ratio = last_ratio + step
    window = acquisition_trace_window(
        anchors,
        pass_index=pass_index,
        settings=settings,
    )
    upper_ratio = settings.ratio_bounds[1] if window is None else window[1]
    if start_ratio > upper_ratio + 1e-12:
        return []

    target_ratios = np.linspace(
        start_ratio,
        upper_ratio,
        num=max(1, candidate_count),
    )
    configs: list[RMBConfig] = []
    seen: set[RMBConfig] = set()
    for n_qubits in n_qubits_values:
        template = template_config(settings, n_qubits)
        for ratio in target_ratios:
            depth = float(intercept + slope * ratio)
            if len(recent_depths) > 0:
                depth = min(depth, float(recent_depths[-1]))
            depth = float(np.clip(depth, settings.depth_bounds[0], settings.depth_bounds[1]))
            config = config_from_parameters(template=template, depth=depth, ratio=float(ratio))
            if config in seen:
                continue
            if not acquisition_bracket_straddles_surface(
                config=config,
                surface=surface,
                settings=settings,
            ):
                continue
            seen.add(config)
            configs.append(config)
    return configs


def boundary_acquisition_score(
    *,
    config: RMBConfig,
    probability: float,
    existing_points: np.ndarray,
    selected_points: list[np.ndarray],
    data: RMBData,
    settings: ContourFirstExperimentConfig,
) -> float:
    n_runs = data[config].num_runs() if config in data else 0
    if n_runs >= settings.max_shots_per_config:
        return 0.0

    boundary_distance = abs(probability - 0.5)
    if boundary_distance > settings.max_boundary_probability_distance:
        return 0.0
    boundary_relevance = np.exp(-((boundary_distance / settings.boundary_width) ** 2))
    boundary_relevance = boundary_relevance ** settings.boundary_focus_power

    point = scaled_depth_ratio_point(
        depth=float(config.depth),
        ratio=float(config.min_two_qubit_gate_ratio),
        settings=settings,
    )
    if len(existing_points) > 0:
        sparsity = min(1.0, float(np.min(np.linalg.norm(existing_points - point, axis=1))) / settings.sparsity_radius)
    else:
        sparsity = 1.0
    if selected_points:
        diversity = min(1.0, float(np.min(np.linalg.norm(np.asarray(selected_points) - point, axis=1))
                                   ) / settings.batch_diversity_radius)
    else:
        diversity = 1.0

    undersampled = 1.0 - n_runs / max(1, settings.max_shots_per_config)
    sparsity_multiplier = 1.0 + settings.exploration_weight * sparsity
    diversity_multiplier = settings.diversity_floor + (1.0 - settings.diversity_floor) * diversity
    request = MeasurementRequest(config, batched_topup_shots(config, data, settings))
    cost = request_marginal_hqc_cost(request, settings)
    return float(
        boundary_relevance
        * sparsity_multiplier
        * diversity_multiplier
        * undersampled
        * high_ratio_acquisition_bonus(config, settings)
        / (cost ** max(0.0, settings.hqc_cost_power))
    )


def refine_uncertain_boundary_points(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> int:
    """
    Spend leftover budget on uncertain measured points near the fitted contour.
    """
    refinements = 0
    fit_passes = 0
    while budget.can_spend():
        if fit_passes >= max(1, settings.batch_refinement_fit_passes):
            break
        fit_passes += 1
        measured_configs = [
            config
            for config, estimator in data.items()
            if 0 < estimator.num_runs() < settings.max_shots_per_config
        ]
        if not measured_configs:
            break

        try:
            surface = fit_monotone_fidelity_surface(data, settings)
        except (RuntimeError, ValueError):
            surface = None

        scored_configs = [
            (config, refinement_score(config, data, settings, surface=surface))
            for config in measured_configs
        ]
        scored_configs = [
            (config, score)
            for config, score in scored_configs
            if score > 0.0
        ]
        if not scored_configs:
            break

        ranked_configs = sorted(
            scored_configs,
            key=lambda item: refinement_request_efficiency(
                item[0],
                item[1],
                batched_topup_shots(item[0], data, settings),
                settings,
            ),
            reverse=True,
        )
        best_score = max(score for _, score in ranked_configs)
        if best_score <= 0.0:
            break

        before = total_measurements(data)
        candidate_count = max(
            1,
            settings.batch_max_configs * max(1, settings.batch_candidate_multiplier),
        )
        selected_configs = [config for config, _ in ranked_configs[:candidate_count]]
        requests = [
            MeasurementRequest(config, batched_topup_shots(config, data, settings))
            for config in selected_configs
        ]
        requests = fill_requests_toward_batch_cost(
            requests,
            data=data,
            settings=settings,
            spendable_hqc=max(0.0, budget.remaining_hqc),
        )
        spends = spend_configs_batch(
            backend=backend,
            rng=rng,
            data=data,
            requests=requests,
            budget=budget,
            settings=settings,
            diagnostics=diagnostics,
            stage="refinement",
        )
        if not spends or total_measurements(data) == before:
            break
        for spend in spends:
            if spend.spent_shots <= 0:
                continue
            add_diagnostic(
                diagnostics,
                "refinement",
                depth=int(spend.config.depth),
                ratio=round(float(spend.config.min_two_qubit_gate_ratio), 4),
                p=round(float(measured_probability(data, spend.config)), 4),
                std=round(measured_std(data, spend.config), 4),
                runs=data[spend.config].num_runs(),
                score=round(
                    float(refinement_score(spend.config, data, settings, surface=surface)),
                    6,
                ),
                spent=spend.spent_shots,
            )
            refinements += 1

    return refinements


def confirm_traced_anchors(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchors: list[RMBConfig],
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> int:
    """
    Spend leftover budget by repeating already traced contour anchors.

    This reduces uncertainty along the discovered curve without introducing new
    fixed-depth or fixed-ratio exploration bands.
    """
    unique_anchors = []
    seen = set()
    for anchor in anchors:
        if anchor in seen:
            continue
        seen.add(anchor)
        unique_anchors.append(anchor)

    confirmations = 0
    target_runs = min(settings.trace_anchor_min_shots, settings.max_shots_per_config)
    while budget.can_spend():
        candidates = [
            anchor
            for anchor in unique_anchors
            if anchor in data and data[anchor].num_runs() < target_runs
        ]
        if not candidates:
            break

        ranked_candidates = sorted(
            candidates,
            key=lambda config: (
                fidelity_variance(data[config])
                / (
                    request_marginal_hqc_cost(
                        MeasurementRequest(
                            config,
                            batched_topup_shots(
                                config,
                                data,
                                settings,
                                target_runs=target_runs,
                            ),
                        ),
                        settings,
                    )
                    ** max(0.0, settings.hqc_cost_power)
                ),
                -data[config].num_runs(),
            ),
            reverse=True,
        )
        candidate_count = max(
            1,
            settings.batch_max_configs * max(1, settings.batch_candidate_multiplier),
        )
        selected_configs = ranked_candidates[:candidate_count]
        requests = [
            MeasurementRequest(
                config,
                batched_topup_shots(
                    config,
                    data,
                    settings,
                    target_runs=target_runs,
                ),
            )
            for config in selected_configs
        ]
        requests = fill_requests_toward_batch_cost(
            requests,
            data=data,
            settings=settings,
            spendable_hqc=max(0.0, budget.remaining_hqc),
            target_runs=target_runs,
        )
        spends = spend_configs_batch(
            backend=backend,
            rng=rng,
            data=data,
            requests=requests,
            budget=budget,
            settings=settings,
            diagnostics=diagnostics,
            stage="trace_anchor_backfill",
        )
        if not spends:
            break
        for spend in spends:
            if spend.spent_shots <= 0:
                continue
            add_diagnostic(
                diagnostics,
                "trace_anchor_backfill",
                depth=int(spend.config.depth),
                ratio=round(float(spend.config.min_two_qubit_gate_ratio), 4),
                p=round(float(measured_probability(data, spend.config)), 4),
                std=round(measured_std(data, spend.config), 4),
                runs=data[spend.config].num_runs(),
                spent=spend.spent_shots,
            )
            confirmations += 1

    return confirmations


def acquire_boundary_candidates(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
    anchors: list[RMBConfig] | None = None,
    diagnostics: list[dict] | None = None,
) -> int:
    """
    Spend remaining stitched budget on new sparse candidates near the fitted contour.
    """
    acquired = 0
    for pass_index in range(max(0, settings.batch_acquisition_passes)):
        if not budget.can_spend():
            break
        try:
            surface = fit_monotone_fidelity_surface(data, settings)
        except (RuntimeError, ValueError):
            break

        depth_grid, ratio_grid, probabilities = surface.probability_grid(settings)
        measured_points = [
            scaled_depth_ratio_point(
                depth=float(config.depth),
                ratio=float(config.min_two_qubit_gate_ratio),
                settings=settings,
            )
            for config in data
        ]
        existing_points = (
            np.asarray(measured_points, dtype=float)
            if measured_points
            else np.empty((0, 2), dtype=float)
        )

        trace_window = acquisition_trace_window(
            anchors,
            pass_index=pass_index,
            settings=settings,
        )
        scored: list[tuple[float, RMBConfig, np.ndarray]] = []
        fallback_scored: list[tuple[float, RMBConfig, np.ndarray]] = []
        selected_points: list[np.ndarray] = []
        seen: set[RMBConfig] = set()
        n_qubits_values = sorted({config.n_qubits for config in data})
        if not n_qubits_values:
            n_qubits_values = list(settings.n_qubits_values)

        for n_qubits in n_qubits_values:
            template = template_config(settings, n_qubits)
            for depth, ratio, probability in zip(
                depth_grid.ravel(),
                ratio_grid.ravel(),
                probabilities.ravel(),
            ):
                config = config_from_parameters(
                    template=template,
                    depth=float(depth),
                    ratio=float(ratio),
                )
                if config in seen:
                    continue
                seen.add(config)
                score = boundary_acquisition_score(
                    config=config,
                    probability=float(probability),
                    existing_points=existing_points,
                    selected_points=selected_points,
                    data=data,
                    settings=settings,
                )
                if score <= 0.0:
                    continue
                point = scaled_depth_ratio_point(
                    depth=float(config.depth),
                    ratio=float(config.min_two_qubit_gate_ratio),
                    settings=settings,
                )
                if in_acquisition_trace_window(config, trace_window):
                    scored.append((score, config, point))
                else:
                    fallback_scored.append((score, config, point))

        if not scored and not fallback_scored:
            break

        scored.sort(key=lambda item: item[0], reverse=True)
        fallback_scored.sort(key=lambda item: item[0], reverse=True)
        candidate_count = max(1, settings.batch_acquisition_ratio_count)
        selected_configs: list[RMBConfig] = []
        selected_seen: set[RMBConfig] = set()
        high_ratio_threshold = (
            settings.ratio_bounds[0]
            + 0.55 * (settings.ratio_bounds[1] - settings.ratio_bounds[0])
        )
        high_ratio_slots = min(
            candidate_count,
            max(
                0,
                int(round(candidate_count * max(0.0, settings.high_ratio_candidate_fraction))),
            ),
        )

        def try_select_candidate(config: RMBConfig, point: np.ndarray) -> bool:
            if config in selected_seen:
                return False
            if not acquisition_bracket_straddles_surface(
                config=config,
                surface=surface,
                settings=settings,
            ):
                return False
            diversity_score = boundary_acquisition_score(
                config=config,
                probability=float(surface.probability(np.array([[
                    float(config.depth),
                    float(config.min_two_qubit_gate_ratio),
                ]]))[0]),
                existing_points=existing_points,
                selected_points=selected_points,
                data=data,
                settings=settings,
            )
            if diversity_score <= 0.0:
                return False
            selected_configs.append(config)
            selected_seen.add(config)
            selected_points.append(point)
            return True

        predicted_configs = trace_follow_acquisition_candidates(
            anchors,
            surface=surface,
            pass_index=pass_index,
            n_qubits_values=list(n_qubits_values),
            settings=settings,
            candidate_count=candidate_count,
        )
        for config in predicted_configs:
            point = scaled_depth_ratio_point(
                depth=float(config.depth),
                ratio=float(config.min_two_qubit_gate_ratio),
                settings=settings,
            )
            try_select_candidate(config, point)
            if len(selected_configs) >= candidate_count:
                break

        for score, config, point in scored:
            if len(selected_configs) >= high_ratio_slots:
                break
            if float(config.min_two_qubit_gate_ratio) < high_ratio_threshold:
                continue
            try_select_candidate(config, point)

        for score, config, point in scored:
            if config in selected_seen:
                continue
            try_select_candidate(config, point)
            if len(selected_configs) >= candidate_count:
                break

        if not selected_configs and fallback_scored:
            add_diagnostic(
                diagnostics,
                "boundary_acquisition_trace_window_fallback",
                pass_index=pass_index,
                window=None if trace_window is None else [round(trace_window[0], 4), round(trace_window[1], 4)],
                candidates=len(fallback_scored),
            )
            for score, config, point in fallback_scored:
                try_select_candidate(config, point)
                if len(selected_configs) >= candidate_count:
                    break

        if not selected_configs:
            break

        before = total_measurements(data)
        request_seen: set[RMBConfig] = set()
        requests: list[MeasurementRequest] = []
        templates_by_qubits = {
            n_qubits: template_config(settings, n_qubits)
            for n_qubits in sorted({config.n_qubits for config in selected_configs})
        }
        for config in selected_configs:
            requests.extend(
                contour_bracket_requests(
                    config=config,
                    data=data,
                    settings=settings,
                    template=templates_by_qubits[config.n_qubits],
                    seen=request_seen,
                )
            )
        requests = fill_requests_toward_batch_cost(
            requests,
            data=data,
            settings=settings,
            spendable_hqc=max(0.0, budget.remaining_hqc),
        )
        spends = spend_configs_batch(
            backend=backend,
            rng=rng,
            data=data,
            requests=requests,
            budget=budget,
            settings=settings,
            diagnostics=diagnostics,
            stage="boundary_acquisition",
        )
        if not spends or total_measurements(data) == before:
            break
        for spend in spends:
            if spend.spent_shots <= 0:
                continue
            add_diagnostic(
                diagnostics,
                "boundary_acquisition",
                pass_index=pass_index,
                trace_window=None if trace_window is None else [round(trace_window[0], 4), round(trace_window[1], 4)],
                depth=int(spend.config.depth),
                ratio=round(float(spend.config.min_two_qubit_gate_ratio), 4),
                p=round(float(measured_probability(data, spend.config)), 4),
                std=round(measured_std(data, spend.config), 4),
                runs=data[spend.config].num_runs(),
                spent=spend.spent_shots,
            )
            acquired += 1

    return acquired


def plot_v7_boundary_with_batch_labels(
    rmb: RMB,
    settings: ContourFirstExperimentConfig,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
    show: bool = True,
    batch_label_limit: int = 300,
) -> list:
    axes = plot_monotone_fidelity_surface_with_confidence(
        rmb._data,
        settings,
        n_bootstrap=n_bootstrap,
        seed=seed,
        show=False,
    )
    batch_ids_by_config = getattr(rmb, "_batch_config_ids", {})
    summary = getattr(rmb, "_batch_cost_summary", {})
    batch_count = int(summary.get("batched_jobs", 0)) if summary else 0
    if batch_ids_by_config and 0 < batch_count <= batch_label_limit:
        groups = sorted(grouped_by_n_qubits_local(rmb._data).items())
        for ax, (_, group) in zip(axes, groups):
            for config, estimator in group.items():
                if estimator.num_runs() <= 0:
                    continue
                ids = batch_ids_by_config.get(batch_config_key(config))
                if not ids:
                    continue
                ax.text(
                    float(config.depth),
                    float(config.min_two_qubit_gate_ratio),
                    str(min(int(batch_id) for batch_id in ids)),
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white",
                    weight="bold",
                    zorder=6,
                    clip_on=True,
                )

    if show:
        import matplotlib.pyplot as plt
        plt.show()
    return axes


def grouped_by_n_qubits_local(data: RMBData) -> dict[int, RMBData]:
    groups: dict[int, RMBData] = {}
    for config, estimator in data.items():
        groups.setdefault(config.n_qubits, {})[config] = estimator
    return groups


def estimate_boundary(settings: ContourFirstExperimentConfig) -> RMB:
    rng = default_rng(settings.rng_seed)
    backend = make_backend()
    rmb = RMB.default(rng).with_backend(backend)
    data: RMBData = rmb._data
    budget = BudgetState(
        remaining_measurements=settings.measurement_budget,
        remaining_hqc=(
            float(settings.hqc_budget)
            if settings.hqc_budget is not None
            else float(settings.measurement_budget)
        ),
    )
    stop_reason = "budget not exhausted"
    batch_count = 0
    traced_anchors: list[RMBConfig] = []
    diagnostics: list[dict] = []

    for n_qubits in settings.n_qubits_values:
        if not budget.can_spend():
            break
        anchor = find_initial_anchor(
            backend=backend,
            rng=rng,
            data=data,
            n_qubits=n_qubits,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        if anchor is None:
            stop_reason = "no initial contour crossing found"
            continue
        anchor = confirm_initial_anchor(
            backend=backend,
            rng=rng,
            data=data,
            anchor=anchor,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        if settings.verbose:
            print(
                "\nInitial contour anchor: "
                f"n_qubits={anchor.n_qubits}, depth={anchor.depth}, "
                f"ratio={anchor.min_two_qubit_gate_ratio:.2f}, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

        seeded_anchors = seeded_batched_trace_from_anchor(
            backend=backend,
            rng=rng,
            data=data,
            anchor=anchor,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        if seeded_anchors:
            traced_anchors.extend(seeded_anchors)
            batch_count += len(seeded_anchors)
            if settings.verbose:
                print(
                    f"\nSeeded trace complete for n_qubits={n_qubits}: "
                    f"{len(seeded_anchors)} anchors, measurements {total_measurements(data)} / "
                    f"{settings.measurement_budget}."
                )
                print_fit_reports(data, settings)

        trace_start = seeded_anchors[-1] if seeded_anchors else anchor
        anchors = trace_from_anchor(
            backend=backend,
            rng=rng,
            data=data,
            anchor=trace_start,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        if seeded_anchors and anchors and anchors[0] == seeded_anchors[-1]:
            anchors = anchors[1:]
        traced_anchors.extend(anchors)
        batch_count += max(0, len(anchors) - 1)
        if settings.verbose:
            print(
                f"\nContour trace complete for n_qubits={n_qubits}: "
                f"{len(anchors)} anchors, measurements {total_measurements(data)} / "
                f"{settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

    refinements = 0
    if settings.refine_after_trace and traced_anchors and budget.can_spend():
        acquisitions = acquire_boundary_candidates(
            backend=backend,
            rng=rng,
            data=data,
            settings=settings,
            budget=budget,
            anchors=traced_anchors,
            diagnostics=diagnostics,
        )
        refinements += acquisitions
        batch_count += acquisitions
        if acquisitions > 0 and settings.verbose:
            print(
                f"\nBoundary acquisition complete: {acquisitions} extra circuit executions, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

        confirmations = confirm_traced_anchors(
            backend=backend,
            rng=rng,
            data=data,
            anchors=traced_anchors,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        refinements += confirmations
        batch_count += confirmations
        if confirmations > 0 and settings.verbose:
            print(
                f"\nTrace anchor confirmation complete: {confirmations} extra circuit executions, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

        extra_refinements = refine_uncertain_boundary_points(
            backend=backend,
            rng=rng,
            data=data,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        refinements += extra_refinements
        batch_count += extra_refinements
        if extra_refinements > 0 and settings.verbose:
            print(
                f"\nRefinement complete: {extra_refinements} extra circuit executions, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

    if budget.remaining_measurements <= 0:
        stop_reason = "MEASUREMENT BUDGET HIT BEFORE HQC BUDGET"
    elif settings.hqc_budget is not None and budget.remaining_hqc <= 0.0:
        stop_reason = "HQC budget hit"
    elif (
        settings.hqc_budget is not None
        and budget.remaining_hqc <= settings.hqc_base_cost
    ):
        stop_reason = "HQC budget effectively hit: remaining HQC is below the base circuit cost"
    elif refinements > 0:
        stop_reason = "trace reached bounds, then refinement stopped with remaining budget"
    elif stop_reason == "budget not exhausted":
        stop_reason = "contour trace reached requested parameter bounds before exhausting budget"

    batch_saving = budget.native_hqc_estimate - budget.stitched_hqc_spent
    batch_saving_fraction = (
        batch_saving / budget.native_hqc_estimate
        if budget.native_hqc_estimate > 0.0
        else 0.0
    )
    rmb._batch_cost_summary = {
        "stitched_hqc_spent": budget.stitched_hqc_spent,
        "native_hqc_estimate": budget.native_hqc_estimate,
        "estimated_batching_saving_hqc": batch_saving,
        "estimated_batching_saving_fraction": batch_saving_fraction,
        "batched_jobs": budget.batched_jobs,
        "max_batched_job_size": budget.max_batched_job_size,
        "max_cost_per_batch": max_cost_per_batch(settings),
    }
    rmb._batch_config_ids = budget.batch_config_ids

    if settings.verbose:
        print_experiment_summary(
            data=data,
            settings=settings,
            stop_reason=stop_reason,
            batch_count=batch_count,
            circuit_executions=budget.circuit_executions,
            max_execution_repeats=budget.max_execution_repeats,
            remaining_measurements=budget.remaining_measurements,
            remaining_hqc=budget.remaining_hqc,
        )
        if settings.batching_enabled:
            print(f"  batched jobs: {budget.batched_jobs}")
            print(f"  max stitched subcircuits in one job: {budget.max_batched_job_size}")
            print(f"  stitched HQC spent: {budget.stitched_hqc_spent:.3f}")
            print(f"  native separate-job HQC estimate: {budget.native_hqc_estimate:.3f}")
            print(
                f"  estimated batching saving: {batch_saving:.3f} "
                f"({batch_saving_fraction:.1%})"
            )

    if settings.save_path is not None:
        rmb.save(settings.save_path)
        if settings.verbose:
            print(f"\nSaved boundary data to {settings.save_path}")

    if settings.diagnostics_path is not None:
        payload = {
            "settings": {
                "measurement_budget": settings.measurement_budget,
                "hqc_budget": settings.hqc_budget,
                "n_qubits_values": settings.n_qubits_values,
                "depth_bounds": settings.depth_bounds,
                "ratio_bounds": settings.ratio_bounds,
                "ray_ratio_count": settings.ray_ratio_count,
                "batched_initial_anchor_search": settings.batched_initial_anchor_search,
                "initial_anchor_depth_grid_count": settings.initial_anchor_depth_grid_count,
                "initial_anchor_grid_shots": settings.initial_anchor_grid_shots,
                "trace_ratio_step_fraction": settings.trace_ratio_step_fraction,
                "trace_depth_search_fraction": settings.trace_depth_search_fraction,
                "seeded_initial_trace": settings.seeded_initial_trace,
                "seeded_trace_ratio_count": settings.seeded_trace_ratio_count,
                "seeded_trace_max_batches": settings.seeded_trace_max_batches,
                "seeded_trace_ratio_span_fraction": settings.seeded_trace_ratio_span_fraction,
                "seeded_trace_depth_power": settings.seeded_trace_depth_power,
                "high_ratio_projected_trace": settings.high_ratio_projected_trace,
                "high_ratio_projected_threshold": settings.high_ratio_projected_threshold,
                "high_ratio_low_depth_fraction": settings.high_ratio_low_depth_fraction,
                "trace_accept_probability_width": settings.trace_accept_probability_width,
                "trace_reject_probability_width": settings.trace_reject_probability_width,
                "trace_local_stencil_points": settings.trace_local_stencil_points,
                "trace_local_stencil_shots": settings.trace_local_stencil_shots,
                "trace_anchor_min_shots": settings.trace_anchor_min_shots,
                "batching_enabled": settings.batching_enabled,
                "batch_max_configs": settings.batch_max_configs,
                "batch_max_hqc_cost": settings.batch_max_hqc_cost,
                "max_cost_per_batch": max_cost_per_batch(settings),
                "batch_reset_weight": settings.batch_reset_weight,
                "batch_candidate_multiplier": settings.batch_candidate_multiplier,
                "batch_refinement_shots": settings.batch_refinement_shots,
                "batch_refinement_fit_passes": settings.batch_refinement_fit_passes,
                "batch_acquisition_passes": settings.batch_acquisition_passes,
                "batch_acquisition_ratio_count": settings.batch_acquisition_ratio_count,
                "batch_target_fill_fraction": settings.batch_target_fill_fraction,
                "batch_fill_repeats": settings.batch_fill_repeats,
                "batch_fill_max_shots_per_config": settings.batch_fill_max_shots_per_config,
                "batch_discovery_fill_max_shots_per_config": (
                    settings.batch_discovery_fill_max_shots_per_config
                ),
                "batch_fill_all_stages": settings.batch_fill_all_stages,
                "high_ratio_acquisition_fraction": settings.high_ratio_acquisition_fraction,
                "high_ratio_candidate_fraction": settings.high_ratio_candidate_fraction,
                "contour_bracket_probe_shots": settings.contour_bracket_probe_shots,
                "contour_bracket_depth_fractions": settings.contour_bracket_depth_fractions,
                "contour_bracket_max_relative_depth": settings.contour_bracket_max_relative_depth,
                "acquisition_require_bracket_straddle": settings.acquisition_require_bracket_straddle,
                "acquisition_follow_trace": settings.acquisition_follow_trace,
                "acquisition_trace_backtrack_fraction": settings.acquisition_trace_backtrack_fraction,
                "acquisition_trace_extension_fraction": settings.acquisition_trace_extension_fraction,
                "acquisition_trace_predict_points": settings.acquisition_trace_predict_points,
                "batch_post_trace_reserve_fraction": settings.batch_post_trace_reserve_fraction,
                "cost_model": "quantinuum_stitched_estimate",
                "execution_backend": "local_sympleq",
            },
            "stop_reason": stop_reason,
            "total_measurements": total_measurements(data),
            "n_configs": len(data),
            "batched_jobs": budget.batched_jobs,
            "max_batched_job_size": budget.max_batched_job_size,
            "batch_config_ids": budget.batch_config_ids,
            **rmb._batch_cost_summary,
            "events": diagnostics,
        }
        diagnostics_path = Path(settings.diagnostics_path)
        diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
        diagnostics_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        if settings.verbose or settings.print_diagnostics:
            print(f"\nSaved diagnostics to {diagnostics_path}")
            if settings.print_diagnostics:
                for event in diagnostics:
                    print(event)

    return rmb


if __name__ == "__main__":
    settings = script_default_settings()

    rmb = estimate_boundary(settings)
    print_fit_reports(rmb._data, settings)
    plot_v7_boundary_with_batch_labels(
        rmb,
        settings,
        n_bootstrap=100,
        seed=settings.rng_seed,
    )
