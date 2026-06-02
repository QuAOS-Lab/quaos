from __future__ import annotations

from dataclasses import dataclass
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
    trace_ratio_step_fraction: float = 0.06
    trace_depth_search_fraction: float = 0.12
    trace_correction_steps: int = 4
    trace_shots: int = 2
    trace_decision_min_shots: int = 4
    crossing_confirm_candidates: int = 3
    trace_accept_probability_width: float = 0.12
    trace_directions: tuple[int, ...] = (1,)
    model_projection_after_fit: bool = True
    refine_after_trace: bool = True
    trace_anchor_min_shots: int = 6
    refinement_shots: int = 2
    refinement_boundary_width: float = 0.15
    refinement_empirical_boundary_width: float = 0.20
    save_path: str | Path | None = "viarregio5_boundary.json"


@dataclass
class BudgetState:
    remaining_measurements: int
    remaining_hqc: float
    circuit_executions: int = 0
    max_execution_repeats: int = 0

    def can_spend(self) -> bool:
        return self.remaining_measurements > 0 and self.remaining_hqc > 0.0


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


def spend_config(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    requested_shots: int,
    budget: BudgetState,
    settings: ContourFirstExperimentConfig,
) -> int:
    if requested_shots <= 0 or not budget.can_spend():
        return 0
    shots = min(requested_shots, budget.remaining_measurements)
    shots = affordable_shot_count(
        config=config,
        requested_shots=shots,
        remaining_hqc=budget.remaining_hqc,
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
    )
    return config, measured_probability(data, config)


def confirm_trace_candidate_if_needed(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    probability: float,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
) -> float:
    """
    Re-measure an apparently near-boundary trace candidate before accepting it.
    """
    if abs(probability - 0.5) > settings.trace_accept_probability_width:
        return probability
    if config not in data:
        return probability

    target = min(settings.trace_decision_min_shots, settings.max_shots_per_config)
    current_runs = data[config].num_runs()
    if current_runs >= target:
        return probability

    additional_shots = target - current_runs
    spend_config(
        backend=backend,
        rng=rng,
        data=data,
        config=config,
        requested_shots=additional_shots,
        budget=budget,
        settings=settings,
    )
    updated_probability = measured_probability(data, config)
    return probability if updated_probability is None else updated_probability


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
) -> RMBConfig | None:
    """
    Find a monotone crossing in depth at fixed ratio.
    """
    d_min, d_max = settings.depth_bounds
    if center_depth is None or not local:
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
    )
    if low_p is None or high_p is None:
        return None
    candidates: list[tuple[RMBConfig, float]] = [(low_config, low_p), (high_config, high_p)]

    if local and not (low_p >= 0.5 and high_p <= 0.5):
        expand = settings.trace_depth_search_fraction * (d_max - d_min)
        for _ in range(3):
            if not budget.can_spend():
                return None
            if low_p < 0.5:
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
                )
                if low_p is not None:
                    candidates.append((low_config, low_p))
            if high_p > 0.5:
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
                )
                if high_p is not None:
                    candidates.append((high_config, high_p))
            if low_p is not None and high_p is not None and low_p >= 0.5 and high_p <= 0.5:
                break
            if low_depth <= d_min and high_depth >= d_max:
                break
            expand *= 1.5

    if not (low_p >= 0.5 and high_p <= 0.5):
        return None

    best_config = low_config if abs(low_p - 0.5) <= abs(high_p - 0.5) else high_config
    best_error = min(abs(low_p - 0.5), abs(high_p - 0.5))
    n_steps = settings.trace_correction_steps if local else settings.ray_bisection_steps
    shots = settings.trace_shots if local else settings.ray_bisection_shots

    for _ in range(n_steps):
        if not budget.can_spend():
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
        )
        if mid_p is None:
            break
        candidates.append((mid_config, mid_p))
        error = abs(mid_p - 0.5)
        if error < best_error:
            best_error = error
            best_config = mid_config
        if mid_p >= 0.5:
            low_depth = float(mid_config.depth)
            low_config = mid_config
            low_p = mid_p
        else:
            high_depth = float(mid_config.depth)
            high_config = mid_config
            high_p = mid_p

    confirmed_config = best_config
    confirmed_error = best_error
    ranked_candidates = sorted(
        candidates,
        key=lambda item: abs(item[1] - 0.5),
    )[: max(1, settings.crossing_confirm_candidates)]
    for candidate_config, probability in ranked_candidates:
        confirmed_probability = confirm_trace_candidate_if_needed(
            backend=backend,
            rng=rng,
            data=data,
            config=candidate_config,
            probability=probability,
            settings=settings,
            budget=budget,
        )
        error = abs(confirmed_probability - 0.5)
        if error < confirmed_error:
            confirmed_config = candidate_config
            confirmed_error = error

    best_config = confirmed_config
    return best_config


def find_initial_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    n_qubits: int,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
) -> RMBConfig | None:
    template = template_config(settings, n_qubits)
    ratios = np.linspace(
        settings.ratio_bounds[0],
        settings.ratio_bounds[1],
        max(1, settings.ray_ratio_count),
    )
    for ratio in ratios:
        if not budget.can_spend():
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
        )
        if anchor is not None:
            return anchor
    return None


def projected_trace_target(
    data: RMBData,
    anchor: RMBConfig,
    direction: int,
    settings: ContourFirstExperimentConfig,
) -> tuple[float, float]:
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    ratio_step = direction * settings.trace_ratio_step_fraction * ratio_span
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
            trial = point + tangent_scaled * settings.trace_ratio_step_fraction * scaled_step
            trial[1] = next_ratio
            trial[0] = np.clip(trial[0], settings.depth_bounds[0], settings.depth_bounds[1])
            projected = project_to_contour(surface, trial, gradient, settings)
            next_depth = float(projected[0])

    return next_depth, next_ratio


def trace_projected_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    current: RMBConfig,
    depth_guess: float,
    ratio: float,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
) -> RMBConfig | None:
    """
    Measure the projected contour point and make short depth corrections only if needed.
    """
    template = template_config(settings, current.n_qubits)
    depth = float(np.clip(depth_guess, settings.depth_bounds[0], settings.depth_bounds[1]))
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
        return None
    best_p = confirm_trace_candidate_if_needed(
        backend=backend,
        rng=rng,
        data=data,
        config=best_config,
        probability=best_p,
        settings=settings,
        budget=budget,
    )

    best_error = abs(best_p - 0.5)
    if best_error <= settings.trace_accept_probability_width:
        return best_config

    d_min, d_max = settings.depth_bounds
    correction_step = settings.trace_depth_search_fraction * (d_max - d_min)
    for _ in range(settings.trace_correction_steps):
        if not budget.can_spend():
            break
        if best_p >= 0.5:
            depth = min(float(d_max), float(best_config.depth) + correction_step)
        else:
            depth = max(float(d_min), float(best_config.depth) - correction_step)
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
        p = confirm_trace_candidate_if_needed(
            backend=backend,
            rng=rng,
            data=data,
            config=candidate,
            probability=p,
            settings=settings,
            budget=budget,
        )
        error = abs(p - 0.5)
        if error < best_error:
            best_config = candidate
            best_p = p
            best_error = error
        if error <= settings.trace_accept_probability_width:
            return candidate
        correction_step *= 0.5

    return best_config


def trace_from_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchor: RMBConfig,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
) -> list[RMBConfig]:
    anchors = [anchor]
    for direction in settings.trace_directions:
        current = anchor
        while budget.can_spend():
            depth_guess, ratio = projected_trace_target(data, current, direction, settings)
            if abs(ratio - current.min_two_qubit_gate_ratio) < 1e-9:
                break
            next_anchor = trace_projected_anchor(
                backend=backend,
                rng=rng,
                data=data,
                current=current,
                depth_guess=depth_guess,
                ratio=float(ratio),
                settings=settings,
                budget=budget,
            )
            if next_anchor is None:
                break
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
    mean = fidelity_mean(estimator)
    empirical_boundary_score = np.exp(
        -((abs(mean - 0.5) / settings.refinement_empirical_boundary_width) ** 2)
    )
    try:
        surface = fit_monotone_fidelity_surface(data, settings)
        p_model = float(surface.probability(np.array([[
            float(config.depth),
            float(config.min_two_qubit_gate_ratio),
        ]]))[0])
    except (RuntimeError, ValueError):
        p_model = mean

    model_boundary_score = np.exp(-((abs(p_model - 0.5) / settings.refinement_boundary_width) ** 2))
    boundary_score = max(model_boundary_score, empirical_boundary_score)
    return float(variance_score * boundary_score)


def refine_uncertain_boundary_points(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: ContourFirstExperimentConfig,
    budget: BudgetState,
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

        ranked_configs = sorted(
            (
                (refinement_score(config, data, settings), config)
                for config in measured_configs
            ),
            key=lambda item: item[0],
            reverse=True,
        )

        selected_config = None
        selected_shots = 0
        for score, config in ranked_configs:
            if score <= 0.0:
                break
            remaining_for_config = settings.max_shots_per_config - data[config].num_runs()
            requested_shots = min(settings.refinement_shots, remaining_for_config)
            affordable_shots = affordable_shot_count(
                config=config,
                requested_shots=requested_shots,
                remaining_hqc=budget.remaining_hqc,
                settings=settings,
            )
            if affordable_shots > 0:
                selected_config = config
                selected_shots = affordable_shots
                break

        if selected_config is None:
            break

        before = total_measurements(data)
        spent = spend_config(
            backend=backend,
            rng=rng,
            data=data,
            config=selected_config,
            requested_shots=selected_shots,
            budget=budget,
            settings=settings,
        )
        if spent <= 0 or total_measurements(data) == before:
            break
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
) -> int:
    """
    Revisit traced contour anchors until each has a minimum repeat count.
    """
    unique_anchors = []
    seen = set()
    for anchor in anchors:
        if anchor in seen:
            continue
        seen.add(anchor)
        unique_anchors.append(anchor)

    confirmations = 0
    while budget.can_spend():
        candidates = [
            anchor
            for anchor in unique_anchors
            if anchor in data
            and data[anchor].num_runs() < min(
                settings.trace_anchor_min_shots,
                settings.max_shots_per_config,
            )
        ]
        if not candidates:
            break

        candidates.sort(
            key=lambda config: (
                data[config].num_runs(),
                -fidelity_variance(data[config]),
            )
        )

        selected_config = None
        selected_shots = 0
        for config in candidates:
            target = min(settings.trace_anchor_min_shots, settings.max_shots_per_config)
            remaining_to_target = target - data[config].num_runs()
            requested_shots = min(settings.refinement_shots, remaining_to_target)
            affordable_shots = affordable_shot_count(
                config=config,
                requested_shots=requested_shots,
                remaining_hqc=budget.remaining_hqc,
                settings=settings,
            )
            if affordable_shots > 0:
                selected_config = config
                selected_shots = affordable_shots
                break

        if selected_config is None:
            break

        spent = spend_config(
            backend=backend,
            rng=rng,
            data=data,
            config=selected_config,
            requested_shots=selected_shots,
            budget=budget,
            settings=settings,
        )
        if spent <= 0:
            break
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
        )
        if anchor is None:
            stop_reason = "no initial contour crossing found"
            continue
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
    if settings.refine_after_trace and budget.can_spend():
        confirmations = confirm_traced_anchors(
            backend=backend,
            rng=rng,
            data=data,
            anchors=traced_anchors,
            settings=settings,
            budget=budget,
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

    return rmb


if __name__ == "__main__":
    settings = ContourFirstExperimentConfig(
        measurement_budget=1000,
        hqc_budget=200.0,
        n_qubits_values=(20,),
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
        trace_decision_min_shots=4,
        crossing_confirm_candidates=3,
        trace_accept_probability_width=0.12,
        trace_directions=(1,),
        refine_after_trace=True,
        trace_anchor_min_shots=6,
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
