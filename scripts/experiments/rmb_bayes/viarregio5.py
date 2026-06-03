from __future__ import annotations

import json
from dataclasses import dataclass, replace
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
    affordable_shot_count,
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
    trace_reserve_budget_fraction: float = 0.35
    trace_reserve_min_hqc: float = 25.0
    trace_ratio_step_fraction: float = 0.06
    trace_depth_search_fraction: float = 0.12
    trace_step_shrink_attempts: int = 4
    trace_local_crossing: bool = False
    trace_enforce_monotone_depth: bool = True
    trace_correction_steps: int = 4
    trace_shots: int = 2
    trace_accept_probability_width: float = 0.12
    trace_reject_probability_width: float = 0.25
    trace_directions: tuple[int, ...] = (1,)
    model_projection_after_fit: bool = True
    refine_after_trace: bool = True
    trace_anchor_min_shots: int = 8
    refinement_shots: int = 2
    refinement_boundary_width: float = 0.15
    diagnostics_path: str | Path | None = "viarregio5_diagnostics.json"
    print_diagnostics: bool = False
    save_path: str | Path | None = "viarregio5_boundary.json"


@dataclass
class BudgetState:
    remaining_measurements: int
    remaining_hqc: float
    circuit_executions: int = 0
    max_execution_repeats: int = 0

    def can_spend(self, reserve_hqc: float = 0.0) -> bool:
        return self.remaining_measurements > 0 and self.remaining_hqc > reserve_hqc


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
) -> int:
    if requested_shots <= 0 or not budget.can_spend(reserve_hqc=reserve_hqc):
        return 0
    shots = min(requested_shots, budget.remaining_measurements)
    spendable_hqc = max(0.0, budget.remaining_hqc - reserve_hqc)
    shots = affordable_shot_count(
        config=config,
        requested_shots=shots,
        remaining_hqc=spendable_hqc,
        settings=settings,
    )
    if shots <= 0:
        return 0

    spent = spend_measurements(
        backend=backend,
        rng=rng,
        data=data,
        config=config,
        n_measurements=shots,
    )
    budget.remaining_measurements -= spent
    if spent > 0:
        budget.circuit_executions += 1
        budget.max_execution_repeats = max(budget.max_execution_repeats, spent)
    if settings.hqc_budget is not None:
        budget.remaining_hqc -= hqc_cost(config, spent, settings)
    else:
        budget.remaining_hqc = float(budget.remaining_measurements)
    return spent


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
    )
    return config, measured_probability(data, config)


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
    )
    confirmed_probability = measured_probability(data, config)
    return probability if confirmed_probability is None else confirmed_probability


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

    low_config, low_p = probe_depth(
        backend=backend,
        rng=rng,
        data=data,
        template=template,
        depth=low_depth,
        ratio=ratio,
        shots=settings.ray_probe_shots if not local else settings.trace_shots,
        budget=budget,
        settings=settings,
        reserve_hqc=reserve_hqc,
    )
    high_config, high_p = probe_depth(
        backend=backend,
        rng=rng,
        data=data,
        template=template,
        depth=high_depth,
        ratio=ratio,
        shots=settings.ray_probe_shots if not local else settings.trace_shots,
        budget=budget,
        settings=settings,
        reserve_hqc=reserve_hqc,
    )
    if low_p is None or high_p is None:
        add_diagnostic(
            diagnostics,
            "crossing_failed",
            stage=stage,
            ratio=round(float(ratio), 4),
            local=local,
            reason="missing_probe",
        )
        return None
    crossing_candidates: list[tuple[RMBConfig, float]] = [(low_config, low_p), (high_config, high_p)]

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
            if low_p is not None and low_p < 0.5:
                low_depth = max(float(d_min), low_depth - expand)
                low_config, low_p = probe_depth(
                    backend=backend,
                    rng=rng,
                    data=data,
                    template=template,
                    depth=low_depth,
                    ratio=ratio,
                    shots=settings.trace_shots,
                    budget=budget,
                    settings=settings,
                    reserve_hqc=reserve_hqc,
                )
                if low_p is not None:
                    crossing_candidates.append((low_config, low_p))
            if high_p is not None and high_p > 0.5:
                high_depth = min(float(d_max), high_depth + expand)
                high_config, high_p = probe_depth(
                    backend=backend,
                    rng=rng,
                    data=data,
                    template=template,
                    depth=high_depth,
                    ratio=ratio,
                    shots=settings.trace_shots,
                    budget=budget,
                    settings=settings,
                    reserve_hqc=reserve_hqc,
                )
                if high_p is not None:
                    crossing_candidates.append((high_config, high_p))
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
        p=None if measured_probability(data, confirmed_anchor) is None else round(float(measured_probability(data, confirmed_anchor)), 4),
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
    if settings.trace_local_crossing:
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
            diagnostics=diagnostics,
            stage="trace_local",
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
        if not budget.can_spend():
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
    for direction in settings.trace_directions:
        current = anchor
        while budget.can_spend():
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
                p=None if measured_probability(data, next_anchor) is None else round(float(measured_probability(data, next_anchor)), 4),
                std=round(measured_std(data, next_anchor), 4),
                runs=data[next_anchor].num_runs() if next_anchor in data else 0,
            )
            anchors.append(next_anchor)
            current = next_anchor
    return anchors


def refinement_score(
    config: RMBConfig,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
) -> float:
    if config not in data:
        return 0.0
    estimator = data[config]
    if estimator.num_runs() >= settings.max_shots_per_config:
        return 0.0

    variance = fidelity_variance(estimator)
    variance_score = min(1.0, variance / (1.0 / 12.0))
    try:
        surface = fit_monotone_fidelity_surface(data, settings)
        p_model = float(surface.probability(np.array([[
            float(config.depth),
            float(config.min_two_qubit_gate_ratio),
        ]]))[0])
    except (RuntimeError, ValueError):
        p_model = fidelity_mean(estimator)

    boundary_score = np.exp(-((abs(p_model - 0.5) / settings.refinement_boundary_width) ** 2))
    return float(variance_score * boundary_score)


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
    while budget.can_spend():
        measured_configs = [
            config
            for config, estimator in data.items()
            if 0 < estimator.num_runs() < settings.max_shots_per_config
        ]
        if not measured_configs:
            break

        best_config = max(
            measured_configs,
            key=lambda config: refinement_score(config, data, settings),
        )
        best_score = refinement_score(best_config, data, settings)
        if best_score <= 0.0:
            break

        before = total_measurements(data)
        spent = spend_config(
            backend=backend,
            rng=rng,
            data=data,
            config=best_config,
            requested_shots=settings.refinement_shots,
            budget=budget,
            settings=settings,
        )
        if spent <= 0 or total_measurements(data) == before:
            break
        add_diagnostic(
            diagnostics,
            "refinement",
            depth=int(best_config.depth),
            ratio=round(float(best_config.min_two_qubit_gate_ratio), 4),
            p=round(float(measured_probability(data, best_config)), 4),
            std=round(measured_std(data, best_config), 4),
            runs=data[best_config].num_runs(),
            score=round(float(best_score), 6),
            spent=spent,
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

        selected = max(
            candidates,
            key=lambda config: (
                fidelity_variance(data[config]),
                -data[config].num_runs(),
            ),
        )
        requested_shots = min(
            settings.refinement_shots,
            target_runs - data[selected].num_runs(),
        )
        spent = spend_config(
            backend=backend,
            rng=rng,
            data=data,
            config=selected,
            requested_shots=requested_shots,
            budget=budget,
            settings=settings,
        )
        if spent <= 0:
            break
        add_diagnostic(
            diagnostics,
            "trace_anchor_backfill",
            depth=int(selected.depth),
            ratio=round(float(selected.min_two_qubit_gate_ratio), 4),
            p=round(float(measured_probability(data, selected)), 4),
            std=round(measured_std(data, selected), 4),
            runs=data[selected].num_runs(),
            spent=spent,
        )
        confirmations += 1

    return confirmations


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

        anchors = trace_from_anchor(
            backend=backend,
            rng=rng,
            data=data,
            anchor=anchor,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
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
                "trace_ratio_step_fraction": settings.trace_ratio_step_fraction,
                "trace_depth_search_fraction": settings.trace_depth_search_fraction,
                "trace_accept_probability_width": settings.trace_accept_probability_width,
                "trace_reject_probability_width": settings.trace_reject_probability_width,
                "trace_anchor_min_shots": settings.trace_anchor_min_shots,
            },
            "stop_reason": stop_reason,
            "total_measurements": total_measurements(data),
            "n_configs": len(data),
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
    settings = ContourFirstExperimentConfig(
        measurement_budget=10000000,
        hqc_budget=250.0,
        n_qubits_values=(10,),
        depth_bounds=(4, 300),
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
        ray_ratio_count=5,
        ray_probe_shots=2,
        ray_bisection_steps=5,
        ray_bisection_shots=2,
        trace_ratio_step_fraction=0.06,
        trace_depth_search_fraction=0.06,
        trace_correction_steps=2,
        trace_shots=2,
        trace_accept_probability_width=0.12,
        trace_directions=(1,),
        refine_after_trace=True,
        refinement_shots=2,
        refinement_boundary_width=0.15,
        hqc_cost_informed_acquisition=True,
        hqc_cost_power=1.0,
        save_path="viarregio5_boundary.json",
        verbose=True,
    )

    rmb = estimate_boundary(settings)
    print_fit_reports(rmb._data, settings)
    plot_monotone_fidelity_surface_with_confidence(
        rmb._data,
        settings,
        n_bootstrap=100,
        seed=settings.rng_seed,
    )
