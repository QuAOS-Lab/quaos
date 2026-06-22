from __future__ import annotations

import json
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
    template_config,
    total_measurements,
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
class SimpleContourExperimentConfig(MonotoneBoundaryExperimentConfig):
    """
    Minimal contour-first estimator.

    The algorithm has three stages:
      1. Search in depth at a fixed low two-qubit ratio until a confident
         p=0.5 crossing is found.
      2. Trace the contour by stepping in ratio and solving only local depth
         crossings around the predicted next point.
      3. Spend leftover budget repeating uncertain near-boundary points.
    """

    search_ratio: float | None = None
    search_depth_center_fraction: float = 0.35
    search_depth_half_width_fraction: float = 0.22
    search_expand_factor: float = 1.6
    search_max_expansions: int = 3
    crossing_bisection_steps: int = 6
    crossing_probe_shots: int = 2
    crossing_bisection_shots: int = 2
    crossing_min_shots: int = 10
    crossing_target_std: float = 0.14
    crossing_probability_width: float = 0.18

    trace_ratio_step_fraction: float = 0.05
    trace_depth_window_fraction: float = 0.10
    trace_max_expansions: int = 1
    trace_bisection_steps: int = 4
    trace_probe_shots: int = 2
    trace_bisection_shots: int = 2
    trace_min_shots: int = 4
    trace_accept_width: float = 0.22
    trace_directions: tuple[int, ...] = (1,)

    backfill_min_shots: int = 8
    backfill_shots: int = 2
    backfill_boundary_width: float = 0.18
    diagnostics_path: str | Path | None = "viarregio6_diagnostics.json"
    print_diagnostics: bool = False
    save_path: str | Path | None = "viarregio6_boundary.json"


@dataclass
class BudgetState:
    remaining_measurements: int
    remaining_hqc: float
    circuit_executions: int = 0
    max_execution_repeats: int = 0

    def can_spend(self) -> bool:
        return self.remaining_measurements > 0 and self.remaining_hqc > 0.0


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
    settings: SimpleContourExperimentConfig,
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


def measure_depth(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    template: RMBConfig,
    depth: float,
    ratio: float,
    shots: int,
    budget: BudgetState,
    settings: SimpleContourExperimentConfig,
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


def top_up_config(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    target_shots: int,
    budget: BudgetState,
    settings: SimpleContourExperimentConfig,
) -> None:
    if config not in data:
        return
    missing = min(target_shots, settings.max_shots_per_config) - data[config].num_runs()
    if missing > 0:
        spend_config(
            backend=backend,
            rng=rng,
            data=data,
            config=config,
            requested_shots=missing,
            budget=budget,
            settings=settings,
        )


def interior_depth_bracket(
    settings: SimpleContourExperimentConfig,
) -> tuple[float, float]:
    d_min, d_max = settings.depth_bounds
    span = d_max - d_min
    center = d_min + settings.search_depth_center_fraction * span
    half_width = settings.search_depth_half_width_fraction * span
    return (
        max(float(d_min), center - half_width),
        min(float(d_max), center + half_width),
    )


def bracket_crossing(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    template: RMBConfig,
    ratio: float,
    low_depth: float,
    high_depth: float,
    probe_shots: int,
    max_expansions: int,
    expand_fraction: float,
    budget: BudgetState,
    settings: SimpleContourExperimentConfig,
) -> tuple[RMBConfig, float, float, RMBConfig, float, float] | None:
    d_min, d_max = settings.depth_bounds
    low_config, low_p = measure_depth(
        backend=backend,
        rng=rng,
        data=data,
        template=template,
        depth=low_depth,
        ratio=ratio,
        shots=probe_shots,
        budget=budget,
        settings=settings,
    )
    high_config, high_p = measure_depth(
        backend=backend,
        rng=rng,
        data=data,
        template=template,
        depth=high_depth,
        ratio=ratio,
        shots=probe_shots,
        budget=budget,
        settings=settings,
    )
    if low_p is None or high_p is None:
        return None

    expand = expand_fraction * (d_max - d_min)
    for _ in range(max_expansions + 1):
        if low_p >= 0.5 and high_p <= 0.5:
            return low_config, low_depth, low_p, high_config, high_depth, high_p

        if not budget.can_spend():
            return None
        if low_p < 0.5 and low_depth > d_min:
            low_depth = max(float(d_min), low_depth - expand)
            low_config, low_p = measure_depth(
                backend=backend,
                rng=rng,
                data=data,
                template=template,
                depth=low_depth,
                ratio=ratio,
                shots=probe_shots,
                budget=budget,
                settings=settings,
            )
        if high_p > 0.5 and high_depth < d_max:
            high_depth = min(float(d_max), high_depth + expand)
            high_config, high_p = measure_depth(
                backend=backend,
                rng=rng,
                data=data,
                template=template,
                depth=high_depth,
                ratio=ratio,
                shots=probe_shots,
                budget=budget,
                settings=settings,
            )
        if low_p is None or high_p is None:
            return None
        if low_depth <= d_min and high_depth >= d_max:
            break
        expand *= settings.search_expand_factor

    return None


def find_depth_crossing(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    template: RMBConfig,
    ratio: float,
    settings: SimpleContourExperimentConfig,
    budget: BudgetState,
    center_depth: float | None = None,
    local: bool = False,
    diagnostics: list[dict] | None = None,
    stage: str = "crossing",
) -> RMBConfig | None:
    if center_depth is None:
        low_depth, high_depth = interior_depth_bracket(settings)
        max_expansions = settings.search_max_expansions
        expand_fraction = settings.search_depth_half_width_fraction
        probe_shots = settings.crossing_probe_shots
        bisection_steps = settings.crossing_bisection_steps
        bisection_shots = settings.crossing_bisection_shots
        target_shots = settings.crossing_min_shots
        target_std = settings.crossing_target_std
    else:
        d_min, d_max = settings.depth_bounds
        half_width = settings.trace_depth_window_fraction * (d_max - d_min)
        low_depth = max(float(d_min), center_depth - half_width)
        high_depth = min(float(d_max), center_depth + half_width)
        max_expansions = settings.trace_max_expansions
        expand_fraction = settings.trace_depth_window_fraction
        probe_shots = settings.trace_probe_shots
        bisection_steps = settings.trace_bisection_steps
        bisection_shots = settings.trace_bisection_shots
        target_shots = settings.trace_min_shots
        target_std = settings.crossing_target_std

    bracket = bracket_crossing(
        backend=backend,
        rng=rng,
        data=data,
        template=template,
        ratio=ratio,
        low_depth=low_depth,
        high_depth=high_depth,
        probe_shots=probe_shots,
        max_expansions=max_expansions,
        expand_fraction=expand_fraction,
        budget=budget,
        settings=settings,
    )
    if bracket is None:
        add_diagnostic(
            diagnostics,
            "crossing_failed",
            stage=stage,
            ratio=round(float(ratio), 4),
            center_depth=None if center_depth is None else round(float(center_depth), 4),
            local=local,
            reason="no_bracket",
        )
        return None

    low_config, low_depth, low_p, high_config, high_depth, high_p = bracket
    add_diagnostic(
        diagnostics,
        "crossing_bracketed",
        stage=stage,
        ratio=round(float(ratio), 4),
        low_depth=round(float(low_depth), 4),
        low_p=round(float(low_p), 4),
        high_depth=round(float(high_depth), 4),
        high_p=round(float(high_p), 4),
        local=local,
    )
    best_config = low_config if abs(low_p - 0.5) <= abs(high_p - 0.5) else high_config
    best_error = min(abs(low_p - 0.5), abs(high_p - 0.5))

    for _ in range(bisection_steps):
        if not budget.can_spend():
            break
        mid_depth = 0.5 * (low_depth + high_depth)
        mid_config, mid_p = measure_depth(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            depth=mid_depth,
            ratio=ratio,
            shots=bisection_shots,
            budget=budget,
            settings=settings,
        )
        if mid_p is None:
            break

        error = abs(mid_p - 0.5)
        if error < best_error:
            best_config = mid_config
            best_error = error

        if mid_p >= 0.5:
            low_depth = float(mid_config.depth)
            low_config = mid_config
            low_p = mid_p
        else:
            high_depth = float(mid_config.depth)
            high_config = mid_config
            high_p = mid_p

    top_up_config(
        backend=backend,
        rng=rng,
        data=data,
        config=best_config,
        target_shots=target_shots,
        budget=budget,
        settings=settings,
    )

    p_best = measured_probability(data, best_config)
    if p_best is None:
        add_diagnostic(
            diagnostics,
            "crossing_failed",
            stage=stage,
            ratio=round(float(ratio), 4),
            local=local,
            reason="no_measurement",
        )
        return None
    if local and abs(p_best - 0.5) > settings.trace_accept_width:
        add_diagnostic(
            diagnostics,
            "crossing_rejected",
            stage=stage,
            ratio=round(float(ratio), 4),
            depth=int(best_config.depth),
            p=round(float(p_best), 4),
            std=round(measured_std(data, best_config), 4),
            local=local,
            reason="outside_trace_accept_width",
        )
        return None
    if not local and (
        abs(p_best - 0.5) > settings.crossing_probability_width
        and measured_std(data, best_config) > target_std
    ):
        add_diagnostic(
            diagnostics,
            "crossing_rejected",
            stage=stage,
            ratio=round(float(ratio), 4),
            depth=int(best_config.depth),
            p=round(float(p_best), 4),
            std=round(measured_std(data, best_config), 4),
            local=local,
            reason="low_confidence_initial_crossing",
        )
        return None
    add_diagnostic(
        diagnostics,
        "crossing_accepted",
        stage=stage,
        ratio=round(float(ratio), 4),
        depth=int(best_config.depth),
        p=round(float(p_best), 4),
        std=round(measured_std(data, best_config), 4),
        runs=data[best_config].num_runs(),
        local=local,
    )
    return best_config


def contour_depth_prediction(
    data: RMBData,
    current: RMBConfig,
    next_ratio: float,
    settings: SimpleContourExperimentConfig,
) -> float:
    try:
        surface = fit_monotone_fidelity_surface(data, settings)
    except (RuntimeError, ValueError):
        return float(current.depth)

    d_min, d_max = settings.depth_bounds
    depths = np.linspace(d_min, d_max, settings.candidate_grid_size[0])
    points = np.column_stack([
        depths,
        np.full_like(depths, next_ratio, dtype=float),
    ])
    probabilities = surface.probability(points)
    if np.min(probabilities) > 0.5 or np.max(probabilities) < 0.5:
        return float(current.depth)
    return float(depths[int(np.argmin(np.abs(probabilities - 0.5)))])


def trace_contour(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchor: RMBConfig,
    settings: SimpleContourExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> list[RMBConfig]:
    anchors = [anchor]
    template = template_config(settings, anchor.n_qubits)
    ratio_min, ratio_max = settings.ratio_bounds
    ratio_step = settings.trace_ratio_step_fraction * (ratio_max - ratio_min)

    for direction in settings.trace_directions:
        current = anchor
        while budget.can_spend():
            next_ratio = float(current.min_two_qubit_gate_ratio + direction * ratio_step)
            if next_ratio < ratio_min or next_ratio > ratio_max:
                add_diagnostic(
                    diagnostics,
                    "trace_stopped",
                    reason="ratio_bounds",
                    direction=direction,
                    current_depth=int(current.depth),
                    current_ratio=round(float(current.min_two_qubit_gate_ratio), 4),
                )
                break

            depth_guess = contour_depth_prediction(data, current, next_ratio, settings)
            add_diagnostic(
                diagnostics,
                "trace_proposed",
                direction=direction,
                current_depth=int(current.depth),
                current_ratio=round(float(current.min_two_qubit_gate_ratio), 4),
                next_ratio=round(float(next_ratio), 4),
                depth_guess=round(float(depth_guess), 4),
            )
            next_anchor = find_depth_crossing(
                backend=backend,
                rng=rng,
                data=data,
                template=template,
                ratio=next_ratio,
                settings=settings,
                budget=budget,
                center_depth=depth_guess,
                local=True,
                diagnostics=diagnostics,
                stage="trace",
            )
            if next_anchor is None:
                add_diagnostic(
                    diagnostics,
                    "trace_stopped",
                    reason="local_crossing_failed",
                    direction=direction,
                    current_depth=int(current.depth),
                    current_ratio=round(float(current.min_two_qubit_gate_ratio), 4),
                    proposed_ratio=round(float(next_ratio), 4),
                    depth_guess=round(float(depth_guess), 4),
                )
                break

            add_diagnostic(
                diagnostics,
                "trace_accepted",
                direction=direction,
                depth=int(next_anchor.depth),
                ratio=round(float(next_anchor.min_two_qubit_gate_ratio), 4),
                p=round(float(measured_probability(data, next_anchor)), 4),
                std=round(measured_std(data, next_anchor), 4),
                runs=data[next_anchor].num_runs(),
            )
            anchors.append(next_anchor)
            current = next_anchor

    return anchors


def backfill_uncertain_points(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchors: list[RMBConfig],
    settings: SimpleContourExperimentConfig,
    budget: BudgetState,
    diagnostics: list[dict] | None = None,
) -> int:
    refinements = 0
    anchor_set = set(anchors)

    while budget.can_spend():
        candidates = [
            config
            for config, estimator in data.items()
            if estimator.num_runs() < settings.max_shots_per_config
        ]
        if not candidates:
            break

        def score(config: RMBConfig) -> float:
            variance = fidelity_variance(data[config])
            if config in anchor_set:
                boundary_score = 1.0
            else:
                mean = fidelity_mean(data[config])
                boundary_score = np.exp(
                    -((abs(mean - 0.5) / settings.backfill_boundary_width) ** 2)
                )
            return float(boundary_score * variance)

        selected = max(candidates, key=score)
        if score(selected) <= 0.0:
            break

        target = (
            settings.backfill_min_shots
            if selected in anchor_set
            else min(settings.backfill_min_shots, settings.max_shots_per_config)
        )
        missing = min(target, settings.max_shots_per_config) - data[selected].num_runs()
        if missing <= 0 and selected in anchor_set:
            anchor_set.remove(selected)
            continue

        spent = spend_config(
            backend=backend,
            rng=rng,
            data=data,
            config=selected,
            requested_shots=min(settings.backfill_shots, max(1, missing)),
            budget=budget,
            settings=settings,
        )
        if spent <= 0:
            break
        add_diagnostic(
            diagnostics,
            "backfill",
            depth=int(selected.depth),
            ratio=round(float(selected.min_two_qubit_gate_ratio), 4),
            p=round(float(measured_probability(data, selected)), 4),
            std=round(measured_std(data, selected), 4),
            runs=data[selected].num_runs(),
            spent=spent,
        )
        refinements += 1

    return refinements


def estimate_boundary(settings: SimpleContourExperimentConfig) -> RMB:
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
    phase_count = 0
    all_anchors: list[RMBConfig] = []
    diagnostics: list[dict] = []

    for n_qubits in settings.n_qubits_values:
        if not budget.can_spend():
            break

        template = template_config(settings, n_qubits)
        ratio = settings.search_ratio
        if ratio is None:
            ratio = settings.ratio_bounds[0]

        anchor = find_depth_crossing(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            ratio=float(ratio),
            settings=settings,
            budget=budget,
            local=False,
            diagnostics=diagnostics,
            stage="initial",
        )
        if anchor is None:
            stop_reason = "no confident initial crossing found"
            continue

        if settings.verbose:
            print(
                "\nStage 1 complete: "
                f"anchor n={anchor.n_qubits}, depth={anchor.depth}, "
                f"ratio={anchor.min_two_qubit_gate_ratio:.2f}, "
                f"p={measured_probability(data, anchor):.3f}, "
                f"std={measured_std(data, anchor):.3f}."
            )
            print_fit_reports(data, settings)

        anchors = trace_contour(
            backend=backend,
            rng=rng,
            data=data,
            anchor=anchor,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        all_anchors.extend(anchors)
        phase_count += max(0, len(anchors) - 1)
        if settings.verbose:
            print(
                f"\nStage 2 complete: traced {len(anchors)} anchors, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

    refinements = 0
    if all_anchors and budget.can_spend():
        refinements = backfill_uncertain_points(
            backend=backend,
            rng=rng,
            data=data,
            anchors=all_anchors,
            settings=settings,
            budget=budget,
            diagnostics=diagnostics,
        )
        phase_count += refinements
        if settings.verbose:
            print(
                f"\nStage 3 complete: {refinements} backfill executions, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

    if budget.remaining_measurements <= 0:
        stop_reason = "MEASUREMENT BUDGET HIT BEFORE HQC BUDGET"
    elif settings.hqc_budget is not None and budget.remaining_hqc <= settings.hqc_base_cost:
        stop_reason = "HQC budget effectively hit"
    elif stop_reason == "budget not exhausted":
        stop_reason = "contour traced/backfilled until no useful spend remained"

    if settings.verbose:
        print_experiment_summary(
            data=data,
            settings=settings,
            stop_reason=stop_reason,
            batch_count=phase_count,
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
        diagnostics_payload = {
            "settings": {
                "measurement_budget": settings.measurement_budget,
                "hqc_budget": settings.hqc_budget,
                "n_qubits_values": settings.n_qubits_values,
                "depth_bounds": settings.depth_bounds,
                "ratio_bounds": settings.ratio_bounds,
                "search_ratio": settings.search_ratio,
                "trace_ratio_step_fraction": settings.trace_ratio_step_fraction,
                "trace_depth_window_fraction": settings.trace_depth_window_fraction,
                "crossing_min_shots": settings.crossing_min_shots,
                "trace_min_shots": settings.trace_min_shots,
                "backfill_min_shots": settings.backfill_min_shots,
            },
            "stop_reason": stop_reason,
            "total_measurements": total_measurements(data),
            "n_configs": len(data),
            "events": diagnostics,
        }
        diagnostics_path = Path(settings.diagnostics_path)
        diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
        diagnostics_path.write_text(json.dumps(diagnostics_payload, indent=2), encoding="utf-8")
        if settings.verbose or settings.print_diagnostics:
            print(f"\nSaved diagnostics to {diagnostics_path}")
            if settings.print_diagnostics:
                for event in diagnostics:
                    print(event)

    return rmb


if __name__ == "__main__":
    settings = SimpleContourExperimentConfig(
        measurement_budget=10_000_000,
        hqc_budget=250.0,
        n_qubits_values=(10,),
        depth_bounds=(4, 300),
        ratio_bounds=(0.08, 0.8),
        random_elimination=0.1,
        scrambling_probability=0.0,
        max_shots_per_config=12,
        candidate_grid_size=(100, 100),
        min_fit_points=6,
        monotone_l2=1e-3,
        search_ratio=None,
        crossing_min_shots=10,
        trace_min_shots=4,
        trace_ratio_step_fraction=0.05,
        trace_depth_window_fraction=0.10,
        backfill_min_shots=8,
        hqc_cost_informed_acquisition=True,
        save_path="viarregio6_boundary.json",
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
