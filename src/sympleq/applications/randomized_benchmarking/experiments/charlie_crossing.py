"""
Contour-first tracing of the fidelity = 0.5 line in (total gates, two-qubit
gate ratio) space at fixed n_qubits.

Stripped-down port of ``viarregio7.py`` (charlie_rmb branch). The method
first finds one p = 0.5 anchor on a gate-count grid at the lowest ratio,
then traces the contour by stepping the ratio up and solving a local
fixed-ratio crossing around the previous anchor. The remaining budget is
spent on a monotone-surface-guided acquisition, anchor confirmation, and
refinement of uncertain near-contour points.

Measurements are planned as stitched batches and priced with the same HQC
model as ``level_crossing.py`` (one base submission cost per batch plus the
pytket bare cost of every circuit), and the run stops when the HQC budget is
exhausted. Circuits run on the backend built by ``settings.backend_factory``:
by default the SympleQ emulation of Quantinuum hardware, so nothing is
submitted to Quantinuum.

Code paths that the viarregio7 default profile disabled (seeded traces,
model-projected trace targets, high-ratio projected traces, crossing
candidate confirmation, diagnostics, the measurement-count budget) were
removed in this port; the live behaviour is unchanged.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.integrations.quantinuum.utils import BASE_SIMULATION_COST
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    MeasurementRequest,
    MonotoneFidelitySurface,
    batch_hqc_cost,
    config_points,
    marginal_hqc_cost,
    print_crossing,
    print_experiment_summary,
    print_fit_reports,
    print_progress,
    save_crossings,
    single_circuit_bare_hqc,
    spend_request_batch,
    start_run,
    try_fit_monotone_fidelity_surface,
)


@dataclass(frozen=True)
class CharlieCrossingSettings(CrossingSettings):
    """
    Settings for a contour-first fidelity = 0.5 crossing experiment.

    The problem definition lives in :class:`~.common.CrossingSettings` so the
    contour-first and level-crossing experiments can be compared on the same
    footing; the fields here are the method's own knobs with the viarregio7
    default profile values.
    """
    save_path: str | Path | None = "charlie_crossing.json"

    # Per-config sampling caps.
    max_shots_per_config: int = 10

    # Initial anchor search.
    initial_anchor_grid_count: int = 17
    initial_anchor_grid_shots: int = 2
    initial_crossing_center_fraction: float = 0.8
    initial_crossing_half_width_fraction: float = 0.18
    initial_crossing_expand_factor: float = 1.4
    ray_ratio_count: int = 5
    ray_probe_shots: int = 2
    ray_bisection_steps: int = 5
    ray_bisection_shots: int = 2
    crossing_decision_confirm_shots: int = 4
    crossing_decision_probability_width: float = 0.20
    initial_anchor_refine_shots: int = 6
    initial_anchor_refine_steps: int = 3
    initial_anchor_search_fraction: float = 0.06
    initial_anchor_min_runs: int = 8

    # Contour trace.
    trace_reserve_budget_fraction: float = 0.35
    trace_reserve_min_hqc: float = 25.0
    trace_ratio_step_fraction: float = 0.03
    trace_gates_search_fraction: float = 0.03
    trace_step_shrink_attempts: int = 4
    trace_local_stencil_points: int = 7
    trace_local_stencil_shots: int = 2
    trace_correction_steps: int = 3
    trace_shots: int = 2
    trace_accept_probability_width: float = 0.1
    trace_reject_probability_width: float = 0.25

    # Stitched batch planning.
    batch_max_configs: int = 42
    batch_candidate_multiplier: int = 4
    batch_target_fill_fraction: float = 0.88
    batch_fill_max_shots_per_config: int = 8

    # Post-trace refinement and acquisition.
    refine_after_trace: bool = True
    batch_post_trace_reserve_fraction: float = 0.25
    trace_anchor_min_shots: int = 8
    refinement_shots: int = 2
    batch_refinement_shots: int = 4
    batch_refinement_fit_passes: int = 3
    refinement_boundary_width: float = 0.15
    batch_acquisition_passes: int = 2
    batch_acquisition_ratio_count: int = 6
    boundary_width: float = 0.08
    boundary_focus_power: float = 2.0
    max_boundary_probability_distance: float = 0.25
    exploration_weight: float = 0.15
    diversity_floor: float = 0.5
    sparsity_radius: float = 0.15
    batch_diversity_radius: float = 0.12
    high_ratio_acquisition_fraction: float = 0.35
    contour_bracket_probe_shots: int = 2
    contour_bracket_gates_fractions: tuple[float, ...] = (0.04, 0.08, 0.12)
    contour_bracket_max_relative_gates: float = 0.20
    acquisition_trace_backtrack_fraction: float = 0.08
    acquisition_trace_extension_fraction: float = 0.16


def measured_probability(data: RMBData, config: RMBConfig) -> float | None:
    if config not in data or data[config].num_runs() == 0:
        return None
    return data[config].posterior_mean()


def surface_probability(surface: MonotoneFidelitySurface, config: RMBConfig) -> float:
    """Fitted fidelity of the surface at one config's coordinates."""
    return float(surface.probability(config_points([config]))[0])


def trace_reserve_hqc(settings: CharlieCrossingSettings) -> float:
    """HQC held back during the anchor search so the trace can run."""
    return min(
        settings.hqc_budget,
        max(settings.trace_reserve_min_hqc,
            settings.trace_reserve_budget_fraction * settings.hqc_budget),
    )


def post_trace_reserve_hqc(settings: CharlieCrossingSettings) -> float:
    """HQC held back during the trace so refinement can run."""
    if not settings.refine_after_trace:
        return 0.0
    return min(
        settings.hqc_budget,
        max(float(BASE_SIMULATION_COST),
            settings.batch_post_trace_reserve_fraction * settings.hqc_budget),
    )


def split_measurement_batches(
    requests: list[MeasurementRequest],
    settings: CharlieCrossingSettings,
    remaining_hqc: float,
    reserve_hqc: float,
) -> list[list[MeasurementRequest]]:
    """Greedily pack requests into stitched batches that fit cap and budget."""
    batches: list[list[MeasurementRequest]] = []
    spendable_hqc = max(0.0, remaining_hqc - reserve_hqc)
    if spendable_hqc <= 0.0:
        return batches

    current: list[MeasurementRequest] = []
    for request in requests:
        if request.shots <= 0:
            continue
        if current:
            candidate_cost = batch_hqc_cost(current + [request])
            if (len(current) + 1 <= settings.batch_max_configs
                    and candidate_cost <= spendable_hqc
                    and candidate_cost <= settings.max_cost_per_run):
                current.append(request)
                continue
            batches.append(current)
            spendable_hqc -= batch_hqc_cost(current)
            current = []
            if spendable_hqc <= 0.0:
                break
        cost = batch_hqc_cost([request])
        if cost <= spendable_hqc and cost <= settings.max_cost_per_run:
            current = [request]

    if current:
        batches.append(current)
    return batches


def fill_requests_toward_batch_cost(
    requests: list[MeasurementRequest],
    *,
    data: RMBData,
    settings: CharlieCrossingSettings,
    spendable_hqc: float,
    target_runs: int | None = None,
) -> list[MeasurementRequest]:
    """
    Increase repeats on selected configs until the stitched batch is near full.

    Extra shots go to the config with the best posterior-variance per marginal
    HQC, so a batch that is paying the base submission cost anyway is used to
    shrink the most uncertainty.
    """
    max_target_cost = min(spendable_hqc,
                          settings.batch_target_fill_fraction * settings.max_cost_per_run)
    if max_target_cost <= 0.0:
        return requests

    filled = [request for request in requests if request.shots > 0]
    if not filled:
        return []
    if target_runs is None:
        max_runs = max(settings.max_shots_per_config, settings.batch_fill_max_shots_per_config)
    else:
        max_runs = max(1, target_runs)

    # The candidates' runs, per-shot cost, and score stay fixed while filling;
    # only the chosen request's shot count and the running cost change.
    current_runs = [data[request.config].num_runs() if request.config in data else 0
                    for request in filled]
    shot_costs = [single_circuit_bare_hqc(request.config) for request in filled]
    scores = [
        (data[request.config].posterior_variance() if request.config in data else 1.0 / 12.0)
        / (marginal_hqc_cost(request.config, 1))
        for request in filled
    ]

    current_cost = batch_hqc_cost(filled)
    while current_cost < max_target_cost:
        best_index: int | None = None
        best_score = -np.inf
        for idx, request in enumerate(filled):
            if current_runs[idx] + request.shots >= max_runs:
                continue
            trial_cost = current_cost + shot_costs[idx]
            if trial_cost > spendable_hqc + 1e-9 or trial_cost > max_target_cost + 1e-9:
                continue
            if scores[idx] > best_score:
                best_score = scores[idx]
                best_index = idx
        if best_index is None:
            break
        request = filled[best_index]
        filled[best_index] = MeasurementRequest(request.config, request.shots + 1)
        current_cost += shot_costs[best_index]

    return filled


def spend_configs_batch(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    requests: list[MeasurementRequest],
    budget: Budget,
    settings: CharlieCrossingSettings,
    reserve_hqc: float = 0.0,
    fill: bool = True,
) -> int:
    """
    Execute measurement requests as stitched batches and charge the budget.

    Returns the number of circuits actually run. ``fill`` tops batches up
    with extra repeats toward the target fill fraction of the batch cost cap.
    """
    requests = [request for request in requests if request.shots > 0]
    if not requests or budget.remaining_hqc <= reserve_hqc:
        return 0

    batches = split_measurement_batches(requests, settings, budget.remaining_hqc, reserve_hqc)
    spent_total = 0

    for batch in batches:
        if not batch or budget.remaining_hqc <= reserve_hqc:
            break
        spendable_hqc = max(0.0, budget.remaining_hqc - reserve_hqc)
        if fill:
            batch = fill_requests_toward_batch_cost(
                batch, data=data, settings=settings, spendable_hqc=spendable_hqc)
            if not batch:
                break
        cost = batch_hqc_cost(batch)
        if cost > spendable_hqc + 1e-9 or cost > settings.max_cost_per_run + 1e-9:
            break
        spend_request_batch(backend, rng, data, batch, seed=settings.rng_seed)
        batch_shots = sum(request.shots for request in batch)
        spent_total += batch_shots
        budget.spend_batch(cost, batch_shots)

    return spent_total


def spend_config(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    shots: int,
    budget: Budget,
    settings: CharlieCrossingSettings,
    reserve_hqc: float = 0.0,
    fill: bool = True,
) -> int:
    return spend_configs_batch(
        backend=backend, rng=rng, data=data,
        requests=[MeasurementRequest(config, shots)],
        budget=budget, settings=settings, reserve_hqc=reserve_hqc, fill=fill)


def probe_points_batch(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: CharlieCrossingSettings,
    probes: list[tuple[float, float, int]],
    budget: Budget,
    reserve_hqc: float = 0.0,
) -> list[tuple[RMBConfig, float | None]]:
    """Probe (n_gates, ratio, shots) points in one stitched batch plan."""
    configs = [settings.make_config(n_gates, ratio) for n_gates, ratio, _ in probes]
    spend_configs_batch(
        backend=backend, rng=rng, data=data,
        requests=[MeasurementRequest(config, shots)
                  for config, (_, _, shots) in zip(configs, probes)],
        budget=budget, settings=settings, reserve_hqc=reserve_hqc)
    return [(config, measured_probability(data, config)) for config in configs]


def confirm_crossing_decision(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    probability: float,
    settings: CharlieCrossingSettings,
    budget: Budget,
    local: bool,
    reserve_hqc: float = 0.0,
) -> float:
    """Top up shots on a near-0.5 bisection point before branching on it."""
    if local or settings.crossing_decision_confirm_shots <= 0:
        return probability
    if abs(probability - 0.5) > settings.crossing_decision_probability_width:
        return probability
    if config not in data:
        return probability

    extra_shots = settings.crossing_decision_confirm_shots - data[config].num_runs()
    if extra_shots <= 0:
        return probability
    spend_config(backend=backend, rng=rng, data=data, config=config, shots=extra_shots,
                 budget=budget, settings=settings, reserve_hqc=reserve_hqc)
    confirmed = measured_probability(data, config)
    return probability if confirmed is None else confirmed


def best_fixed_ratio_bracket(
    measured_points: list[tuple[RMBConfig, float]],
) -> tuple[RMBConfig, float, RMBConfig, float] | None:
    """Tightest adjacent pair with p >= 0.5 below and p <= 0.5 above."""
    points = sorted(measured_points, key=lambda item: float(item[0].n_gates))
    best: tuple[RMBConfig, float, RMBConfig, float] | None = None
    best_width = np.inf
    for (low_config, low_p), (high_config, high_p) in zip(points, points[1:]):
        if low_p >= 0.5 and high_p <= 0.5:
            width = float(high_config.n_gates) - float(low_config.n_gates)
            if width < best_width:
                best_width = width
                best = (low_config, low_p, high_config, high_p)
    return best


def initial_gates_window(settings: CharlieCrossingSettings) -> tuple[float, float]:
    """Interior gate-count bracket where the initial crossing is searched."""
    g_min, g_max = settings.n_gates_bounds
    gates_span = g_max - g_min
    center = g_min + settings.initial_crossing_center_fraction * gates_span
    half_width = settings.initial_crossing_half_width_fraction * gates_span
    return max(float(g_min), center - half_width), min(float(g_max), center + half_width)


def find_crossing_at_ratio(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    ratio: float,
    settings: CharlieCrossingSettings,
    budget: Budget,
    center_gates: float | None = None,
    reserve_hqc: float = 0.0,
) -> RMBConfig | None:
    """
    Find a monotone crossing in total gates at fixed ratio.

    Without ``center_gates`` the search brackets an interior window of the
    gate bounds and bisects; with it, a local search probes a short stencil
    around ``center_gates`` and refines. Both expand the bracket when the
    endpoints sit on the same side.
    """
    g_min, g_max = settings.n_gates_bounds
    local = center_gates is not None
    if center_gates is None:
        low_gates, high_gates = initial_gates_window(settings)
    else:
        span = settings.trace_gates_search_fraction * (g_max - g_min)
        low_gates = max(float(g_min), center_gates - span)
        high_gates = min(float(g_max), center_gates + span)

    probe_shots = settings.ray_probe_shots if not local else settings.trace_local_stencil_shots
    if local and settings.trace_local_stencil_points > 2:
        stencil_gates = np.linspace(low_gates, high_gates,
                                    max(2, settings.trace_local_stencil_points))
        stencil_results = probe_points_batch(
            backend=backend, rng=rng, data=data, settings=settings,
            probes=[(float(gates), ratio, probe_shots) for gates in stencil_gates],
            budget=budget, reserve_hqc=reserve_hqc)
        measured_stencil = [(config, p) for config, p in stencil_results if p is not None]
        bracket = best_fixed_ratio_bracket(measured_stencil)
        if bracket is not None:
            low_config, low_p, high_config, high_p = bracket
            low_gates = float(low_config.n_gates)
            high_gates = float(high_config.n_gates)
        elif measured_stencil:
            low_config, low_p = min(measured_stencil, key=lambda item: float(item[0].n_gates))
            high_config, high_p = max(measured_stencil, key=lambda item: float(item[0].n_gates))
            low_gates = float(low_config.n_gates)
            high_gates = float(high_config.n_gates)
        else:
            low_config = high_config = None
            low_p = high_p = None
    else:
        (low_config, low_p), (high_config, high_p) = probe_points_batch(
            backend=backend, rng=rng, data=data, settings=settings,
            probes=[(low_gates, ratio, probe_shots), (high_gates, ratio, probe_shots)],
            budget=budget, reserve_hqc=reserve_hqc)

    if low_config is None or high_config is None or low_p is None or high_p is None:
        return None

    if not (low_p >= 0.5 and high_p <= 0.5):
        expand = (
            settings.trace_gates_search_fraction * (g_max - g_min)
            if local
            else settings.initial_crossing_half_width_fraction
            * settings.initial_crossing_expand_factor
            * (g_max - g_min)
        )
        for _ in range(3):
            if budget.remaining_hqc <= reserve_hqc:
                return None
            expansion_probes: list[tuple[str, float]] = []
            if low_p is not None and low_p < 0.5:
                low_gates = max(float(g_min), low_gates - expand)
                expansion_probes.append(("low", low_gates))
            if high_p is not None and high_p > 0.5:
                high_gates = min(float(g_max), high_gates + expand)
                expansion_probes.append(("high", high_gates))
            if expansion_probes:
                expansion_results = probe_points_batch(
                    backend=backend, rng=rng, data=data, settings=settings,
                    probes=[(gates, ratio, settings.trace_shots)
                            for _, gates in expansion_probes],
                    budget=budget, reserve_hqc=reserve_hqc)
                for (side, _), (expanded_config, expanded_p) in zip(
                        expansion_probes, expansion_results):
                    if side == "low":
                        low_config, low_p = expanded_config, expanded_p
                    else:
                        high_config, high_p = expanded_config, expanded_p
            if low_p is not None and high_p is not None and low_p >= 0.5 and high_p <= 0.5:
                break
            if low_gates <= g_min and high_gates >= g_max:
                break
            expand *= 1.5

    if low_p is None or high_p is None or not (low_p >= 0.5 and high_p <= 0.5):
        return None

    best_config = low_config if abs(low_p - 0.5) <= abs(high_p - 0.5) else high_config
    best_error = min(abs(low_p - 0.5), abs(high_p - 0.5))
    n_steps = settings.trace_correction_steps if local else settings.ray_bisection_steps
    shots = settings.trace_shots if local else settings.ray_bisection_shots

    for _ in range(n_steps):
        if budget.remaining_hqc <= reserve_hqc:
            break
        mid_gates = 0.5 * (low_gates + high_gates)
        mid_config, mid_p = probe_points_batch(
            backend=backend, rng=rng, data=data, settings=settings,
            probes=[(mid_gates, ratio, shots)],
            budget=budget, reserve_hqc=reserve_hqc)[0]
        if mid_p is None:
            break
        mid_p = confirm_crossing_decision(
            backend=backend, rng=rng, data=data, config=mid_config, probability=mid_p,
            settings=settings, budget=budget, local=local, reserve_hqc=reserve_hqc)
        error = abs(mid_p - 0.5)
        if error < best_error:
            best_error = error
            best_config = mid_config
        # Non-local searches treat a statistically ambiguous midpoint as the
        # high side, so the bracket only branches on confident decisions.
        ambiguous = error <= settings.crossing_decision_probability_width
        if mid_p >= 0.5 and (local or not ambiguous):
            low_gates = float(mid_config.n_gates)
            low_p = mid_p
        else:
            high_gates = float(mid_config.n_gates)
            high_p = mid_p

    at_gates_bound = (best_config.n_gates <= g_min + 1 or best_config.n_gates >= g_max - 1)
    if at_gates_bound and best_error > settings.trace_accept_probability_width:
        return None
    if best_error <= settings.trace_reject_probability_width:
        return best_config
    return None


def find_initial_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: CharlieCrossingSettings,
    budget: Budget,
) -> RMBConfig | None:
    """
    Find one p = 0.5 anchor at the cheapest (lowest) ratio.

    A gate-count grid over the interior bracket is measured in stitched
    batches first; if no adjacent pair straddles 0.5, fall back to bracketed
    ray searches over increasing ratios.
    """
    reserve_hqc = trace_reserve_hqc(settings)
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1],
                         max(1, settings.ray_ratio_count))

    if budget.remaining_hqc > reserve_hqc:
        ratio = float(ratios[0])
        low_gates, high_gates = initial_gates_window(settings)
        gates_grid = np.linspace(low_gates, high_gates,
                                 max(2, settings.initial_anchor_grid_count))
        grid_configs = [settings.make_config(float(gates), ratio) for gates in gates_grid]
        spend_configs_batch(
            backend=backend, rng=rng, data=data,
            requests=[MeasurementRequest(config, settings.initial_anchor_grid_shots)
                      for config in grid_configs],
            budget=budget, settings=settings, reserve_hqc=reserve_hqc, fill=False)

        measured_grid = [(config, p) for config in grid_configs
                         if (p := measured_probability(data, config)) is not None]
        measured_grid.sort(key=lambda item: float(item[0].n_gates))

        anchor_config: RMBConfig | None = None
        anchor_error = np.inf
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
            return anchor_config

    for ratio in ratios:
        if budget.remaining_hqc <= reserve_hqc:
            return None
        anchor = find_crossing_at_ratio(
            backend=backend, rng=rng, data=data, ratio=float(ratio),
            settings=settings, budget=budget, reserve_hqc=reserve_hqc)
        if anchor is not None:
            return anchor
    return None


def confirm_initial_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchor: RMBConfig,
    settings: CharlieCrossingSettings,
    budget: Budget,
) -> RMBConfig:
    """Refine the anchor with a short local search and top up its shots."""
    reserve_hqc = trace_reserve_hqc(settings)
    if budget.remaining_hqc <= reserve_hqc:
        return anchor

    refine_settings = replace(
        settings,
        trace_shots=settings.initial_anchor_refine_shots,
        trace_correction_steps=settings.initial_anchor_refine_steps,
        trace_gates_search_fraction=settings.initial_anchor_search_fraction,
    )
    refined = find_crossing_at_ratio(
        backend=backend, rng=rng, data=data, ratio=float(anchor.ratio_2_qb_gates),
        settings=refine_settings, budget=budget,
        center_gates=float(anchor.n_gates), reserve_hqc=reserve_hqc)
    confirmed_anchor = refined if refined is not None else anchor

    if confirmed_anchor in data:
        extra_shots = settings.initial_anchor_min_runs - data[confirmed_anchor].num_runs()
        if extra_shots > 0:
            spend_config(backend=backend, rng=rng, data=data, config=confirmed_anchor,
                         shots=extra_shots, budget=budget, settings=settings,
                         reserve_hqc=reserve_hqc, fill=False)
    return confirmed_anchor


def trace_from_anchor(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchor: RMBConfig,
    settings: CharlieCrossingSettings,
    budget: Budget,
) -> list[RMBConfig]:
    """
    Trace the contour up in ratio from the anchor.

    Each step proposes the previous anchor's gate count at a higher ratio and
    solves a local fixed-ratio crossing there; failed steps retry with the
    ratio step halved. The crossing gate count is kept non-increasing in the
    ratio, as the fidelity is monotone in both coordinates.
    """
    anchors = [anchor]
    reserve_hqc = post_trace_reserve_hqc(settings)
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    current = anchor

    while budget.remaining_hqc > reserve_hqc:
        next_anchor = None
        for attempt in range(max(1, settings.trace_step_shrink_attempts)):
            step = settings.trace_ratio_step_fraction * (0.5 ** attempt) * ratio_span
            ratio = float(np.clip(current.ratio_2_qb_gates + step,
                                  settings.ratio_bounds[0], settings.ratio_bounds[1]))
            if abs(ratio - current.ratio_2_qb_gates) < 1e-9:
                break
            next_anchor = find_crossing_at_ratio(
                backend=backend, rng=rng, data=data, ratio=ratio,
                settings=settings, budget=budget,
                center_gates=float(current.n_gates), reserve_hqc=reserve_hqc)
            if next_anchor is not None:
                break
        if next_anchor is None:
            break
        if next_anchor.ratio_2_qb_gates <= current.ratio_2_qb_gates:
            # The gate-count rounding could not realize a larger ratio at this
            # step size, so the trace cannot advance any further.
            break
        anchors.append(next_anchor)
        current = next_anchor
        print_crossing(settings, budget, "Trace anchor", next_anchor)
    return anchors


def refinement_score(
    config: RMBConfig,
    data: RMBData,
    settings: CharlieCrossingSettings,
    surface: MonotoneFidelitySurface | None = None,
) -> float:
    """Variance-weighted closeness of a measured config to the fitted contour."""
    if config not in data:
        return 0.0
    estimator = data[config]
    if estimator.num_runs() >= settings.max_shots_per_config:
        return 0.0

    variance_score = min(1.0, estimator.posterior_variance() / (1.0 / 12.0))
    if surface is None:
        p_model = estimator.posterior_mean()
    else:
        p_model = surface_probability(surface, config)

    boundary_score = np.exp(-((abs(p_model - 0.5) / settings.refinement_boundary_width) ** 2))
    return float(variance_score * boundary_score)


def batched_topup_shots(
    config: RMBConfig,
    data: RMBData,
    settings: CharlieCrossingSettings,
    target_runs: int | None = None,
) -> int:
    """Shots needed to push a config toward its target runs, in one increment."""
    current_runs = data[config].num_runs() if config in data else 0
    if target_runs is None:
        target_runs = settings.max_shots_per_config
    missing = max(0, min(target_runs, settings.max_shots_per_config) - current_runs)
    increment = max(settings.refinement_shots, settings.batch_refinement_shots)
    return min(missing, max(1, increment))


def contour_bracket_gates(gates: float, settings: CharlieCrossingSettings) -> list[float]:
    """Gate counts of a fixed-ratio stack bracketing a contour candidate."""
    g_min, g_max = settings.n_gates_bounds
    gates_span = max(1.0, float(g_max - g_min))
    candidate_gates = [float(gates)]
    for fraction in settings.contour_bracket_gates_fractions:
        span_offset = abs(float(fraction)) * gates_span
        relative_offset = settings.contour_bracket_max_relative_gates * max(1.0, float(gates))
        offset = max(1.0, min(span_offset, relative_offset))
        candidate_gates.extend([float(gates) - offset, float(gates) + offset])
    return [gates for gates in candidate_gates if float(g_min) <= gates <= float(g_max)]


def contour_bracket_requests(
    config: RMBConfig,
    data: RMBData,
    settings: CharlieCrossingSettings,
    seen: set[RMBConfig] | None = None,
) -> list[MeasurementRequest]:
    """
    Request a same-ratio gate-count stack around a contour candidate.

    The monotone fit is much better constrained when each ratio has measured
    points on both sides of the 0.5 crossing, not just an isolated
    near-contour point.
    """
    if seen is None:
        seen = set()
    ratio = float(config.ratio_2_qb_gates)
    requests: list[MeasurementRequest] = []
    for gates in contour_bracket_gates(float(config.n_gates), settings):
        bracket_config = settings.make_config(gates, ratio)
        if bracket_config in seen:
            continue
        seen.add(bracket_config)
        current_runs = data[bracket_config].num_runs() if bracket_config in data else 0
        missing = max(0, settings.max_shots_per_config - current_runs)
        shots = min(max(1, settings.contour_bracket_probe_shots), missing)
        if shots > 0:
            requests.append(MeasurementRequest(bracket_config, shots))
    return requests


def scaled_point(gates: float, ratio: float, settings: CharlieCrossingSettings) -> np.ndarray:
    gates_span = max(1e-12, settings.n_gates_bounds[1] - settings.n_gates_bounds[0])
    ratio_span = max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0])
    return np.array([
        (gates - settings.n_gates_bounds[0]) / gates_span,
        (ratio - settings.ratio_bounds[0]) / ratio_span,
    ])


def high_ratio_acquisition_bonus(config: RMBConfig, settings: CharlieCrossingSettings) -> float:
    """Mild score bonus for the cheap high-ratio, low-gate-count corner."""
    if settings.high_ratio_acquisition_fraction <= 0.0:
        return 1.0
    ratio_span = max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0])
    gates_span = max(1e-12, settings.n_gates_bounds[1] - settings.n_gates_bounds[0])
    ratio_position = (float(config.ratio_2_qb_gates) - settings.ratio_bounds[0]) / ratio_span
    low_gates_position = (settings.n_gates_bounds[1] - float(config.n_gates)) / gates_span
    shape_score = np.clip(ratio_position * low_gates_position, 0.0, 1.0)
    return float(1.0 + settings.high_ratio_acquisition_fraction * shape_score)


def acquisition_bracket_straddles_surface(
    config: RMBConfig,
    surface: MonotoneFidelitySurface,
    settings: CharlieCrossingSettings,
) -> bool:
    """Whether the candidate's bracket stack straddles the fitted 0.5 contour."""
    gates = contour_bracket_gates(float(config.n_gates), settings)
    if not gates:
        return False
    points = np.array([[g, float(config.ratio_2_qb_gates)] for g in gates], dtype=float)
    probabilities = surface.probability(points)
    return bool(float(np.min(probabilities)) <= 0.5 <= float(np.max(probabilities)))


def acquisition_trace_window(
    anchors: list[RMBConfig] | None,
    *,
    pass_index: int,
    settings: CharlieCrossingSettings,
) -> tuple[float, float] | None:
    """Ratio window around the end of the trace where acquisition focuses."""
    if not anchors:
        return None
    ratios = [
        float(anchor.ratio_2_qb_gates)
        for anchor in anchors
        if settings.ratio_bounds[0] <= float(anchor.ratio_2_qb_gates) <= settings.ratio_bounds[1]
    ]
    if not ratios:
        return None

    ratio_span = max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0])
    last_ratio = max(ratios)
    backtrack = settings.acquisition_trace_backtrack_fraction * ratio_span
    extension = settings.acquisition_trace_extension_fraction * ratio_span
    lower = max(settings.ratio_bounds[0], last_ratio - backtrack)
    upper = min(settings.ratio_bounds[1], last_ratio + (pass_index + 1) * extension)
    if upper < lower:
        return None
    return lower, upper


def boundary_acquisition_score(
    *,
    config: RMBConfig,
    probability: float,
    point: np.ndarray,
    existing_points: np.ndarray,
    selected_points: list[np.ndarray],
    data: RMBData,
    settings: CharlieCrossingSettings,
) -> float:
    """Score a candidate at its scaled ``point`` by contour closeness, sparsity,
    diversity, and cost."""
    n_runs = data[config].num_runs() if config in data else 0
    if n_runs >= settings.max_shots_per_config:
        return 0.0

    boundary_distance = abs(probability - 0.5)
    if boundary_distance > settings.max_boundary_probability_distance:
        return 0.0
    boundary_relevance = np.exp(-((boundary_distance / settings.boundary_width) ** 2))
    boundary_relevance = boundary_relevance ** settings.boundary_focus_power

    if len(existing_points) > 0:
        sparsity = min(1.0, float(np.min(np.linalg.norm(existing_points - point, axis=1)))
                       / settings.sparsity_radius)
    else:
        sparsity = 1.0
    if selected_points:
        diversity = min(1.0, float(np.min(np.linalg.norm(np.asarray(selected_points) - point,
                                                         axis=1)))
                        / settings.batch_diversity_radius)
    else:
        diversity = 1.0

    undersampled = 1.0 - n_runs / max(1, settings.max_shots_per_config)
    sparsity_multiplier = 1.0 + settings.exploration_weight * sparsity
    diversity_multiplier = settings.diversity_floor + (1.0 - settings.diversity_floor) * diversity
    cost = marginal_hqc_cost(config, batched_topup_shots(config, data, settings))
    return float(
        boundary_relevance
        * sparsity_multiplier
        * diversity_multiplier
        * undersampled
        * high_ratio_acquisition_bonus(config, settings)
        / (cost)
    )


def acquire_boundary_candidates(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: CharlieCrossingSettings,
    budget: Budget,
    anchors: list[RMBConfig] | None = None,
) -> int:
    """Spend remaining stitched budget on new sparse candidates near the contour."""
    acquired = 0
    for pass_index in range(max(0, settings.batch_acquisition_passes)):
        if budget.remaining_hqc <= 0.0:
            break
        # Rebinding the narrowed fit result keeps the closure below typed
        # against a non-optional surface.
        maybe_surface = try_fit_monotone_fidelity_surface(data, settings)
        if maybe_surface is None:
            break
        surface = maybe_surface

        gates_grid, ratio_grid, probabilities = surface.probability_grid(
            settings.candidate_grid_size)
        measured_points = [
            scaled_point(float(config.n_gates), float(config.ratio_2_qb_gates), settings)
            for config in data
        ]
        existing_points = (np.asarray(measured_points, dtype=float)
                           if measured_points else np.empty((0, 2), dtype=float))

        trace_window = acquisition_trace_window(anchors, pass_index=pass_index,
                                                settings=settings)
        scored: list[tuple[float, RMBConfig, np.ndarray]] = []
        fallback_scored: list[tuple[float, RMBConfig, np.ndarray]] = []
        selected_points: list[np.ndarray] = []
        seen: set[RMBConfig] = set()

        # Points far from the contour would score 0 anyway; masking them first
        # skips building configs for most of the grid.
        near_boundary = (np.abs(probabilities.ravel() - 0.5)
                         <= settings.max_boundary_probability_distance)
        for gates, ratio, probability in zip(gates_grid.ravel()[near_boundary],
                                             ratio_grid.ravel()[near_boundary],
                                             probabilities.ravel()[near_boundary]):
            config = settings.make_config(float(gates), float(ratio))
            if config in seen:
                continue
            seen.add(config)
            point = scaled_point(float(config.n_gates), float(config.ratio_2_qb_gates),
                                 settings)
            score = boundary_acquisition_score(
                config=config, probability=float(probability), point=point,
                existing_points=existing_points, selected_points=selected_points,
                data=data, settings=settings)
            if score <= 0.0:
                continue
            in_window = (trace_window is None or trace_window[0] <= float(config.ratio_2_qb_gates) <= trace_window[1])
            if in_window:
                scored.append((score, config, point))
            else:
                fallback_scored.append((score, config, point))

        if not scored and not fallback_scored:
            break

        scored.sort(key=lambda item: item[0], reverse=True)
        fallback_scored.sort(key=lambda item: item[0], reverse=True)
        candidate_count = max(1, settings.batch_acquisition_ratio_count)
        selected_configs: list[RMBConfig] = []

        def try_select_candidate(config: RMBConfig, point: np.ndarray) -> None:
            if config in selected_configs:
                return
            if not acquisition_bracket_straddles_surface(config, surface, settings):
                return
            diversity_score = boundary_acquisition_score(
                config=config,
                probability=surface_probability(surface, config),
                point=point,
                existing_points=existing_points, selected_points=selected_points,
                data=data, settings=settings)
            if diversity_score <= 0.0:
                return
            selected_configs.append(config)
            selected_points.append(point)

        for _, config, point in scored:
            try_select_candidate(config, point)
            if len(selected_configs) >= candidate_count:
                break
        if not selected_configs:
            for _, config, point in fallback_scored:
                try_select_candidate(config, point)
                if len(selected_configs) >= candidate_count:
                    break
        if not selected_configs:
            break

        request_seen: set[RMBConfig] = set()
        requests: list[MeasurementRequest] = []
        for config in selected_configs:
            requests.extend(contour_bracket_requests(config, data, settings,
                                                     seen=request_seen))
        requests = fill_requests_toward_batch_cost(
            requests, data=data, settings=settings,
            spendable_hqc=max(0.0, budget.remaining_hqc))
        spent = spend_configs_batch(
            backend=backend, rng=rng, data=data, requests=requests,
            budget=budget, settings=settings)
        if spent == 0:
            break
        acquired += spent

    return acquired


def confirm_traced_anchors(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    anchors: list[RMBConfig],
    settings: CharlieCrossingSettings,
    budget: Budget,
) -> int:
    """
    Spend leftover budget by repeating already traced contour anchors.

    This reduces uncertainty along the discovered curve without introducing
    new exploration bands.
    """
    unique_anchors = list(dict.fromkeys(anchors))
    confirmations = 0
    target_runs = min(settings.trace_anchor_min_shots, settings.max_shots_per_config)

    while budget.remaining_hqc > 0.0:
        candidates = [anchor for anchor in unique_anchors
                      if anchor in data and data[anchor].num_runs() < target_runs]
        if not candidates:
            break

        ranked_candidates = sorted(
            candidates,
            key=lambda config: (
                data[config].posterior_variance()
                / (marginal_hqc_cost(config,
                                     batched_topup_shots(config, data, settings,
                                                         target_runs=target_runs))),
                -data[config].num_runs(),
            ),
            reverse=True,
        )
        candidate_count = max(1, settings.batch_max_configs * settings.batch_candidate_multiplier)
        requests = [
            MeasurementRequest(config, batched_topup_shots(config, data, settings,
                                                           target_runs=target_runs))
            for config in ranked_candidates[:candidate_count]
        ]
        requests = fill_requests_toward_batch_cost(
            requests, data=data, settings=settings,
            spendable_hqc=max(0.0, budget.remaining_hqc), target_runs=target_runs)
        spent = spend_configs_batch(
            backend=backend, rng=rng, data=data, requests=requests,
            budget=budget, settings=settings)
        if spent == 0:
            break
        confirmations += spent

    return confirmations


def refine_uncertain_boundary_points(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: CharlieCrossingSettings,
    budget: Budget,
) -> int:
    """Spend leftover budget on uncertain measured points near the fitted contour."""
    refinements = 0
    fit_passes = 0
    while budget.remaining_hqc > 0.0 and fit_passes < max(1, settings.batch_refinement_fit_passes):
        fit_passes += 1
        measured_configs = [
            config for config, estimator in data.items()
            if 0 < estimator.num_runs() < settings.max_shots_per_config
        ]
        if not measured_configs:
            break

        surface = try_fit_monotone_fidelity_surface(data, settings)
        scored_configs = [
            (config, refinement_score(config, data, settings, surface=surface))
            for config in measured_configs
        ]
        scored_configs = [(config, score) for config, score in scored_configs if score > 0.0]
        if not scored_configs:
            break

        ranked_configs = sorted(
            scored_configs,
            key=lambda item: item[1] / (marginal_hqc_cost(
                item[0], batched_topup_shots(item[0], data, settings))),
            reverse=True,
        )
        candidate_count = max(1, settings.batch_max_configs * settings.batch_candidate_multiplier)
        requests = [
            MeasurementRequest(config, batched_topup_shots(config, data, settings))
            for config, _ in ranked_configs[:candidate_count]
        ]
        requests = fill_requests_toward_batch_cost(
            requests, data=data, settings=settings,
            spendable_hqc=max(0.0, budget.remaining_hqc))
        spent = spend_configs_batch(
            backend=backend, rng=rng, data=data, requests=requests,
            budget=budget, settings=settings)
        if spent == 0:
            break
        refinements += spent

    return refinements


def run(settings: CharlieCrossingSettings) -> tuple[RMB, list[RMBConfig]]:
    """
    Run the contour-first experiment.

    Returns
    -------
    tuple[RMB, list[RMBConfig]]
        The RMB holding all recorded data and the traced contour anchors in
        increasing ratio order.
    """
    rng, rmb, budget = start_run(settings)
    data = rmb._data
    stop_reason = None
    anchors: list[RMBConfig] = []

    anchor = find_initial_anchor(backend=rmb.backend, rng=rng, data=data,
                                 settings=settings, budget=budget)
    if anchor is None:
        stop_reason = "no initial contour crossing found"
    else:
        anchor = confirm_initial_anchor(backend=rmb.backend, rng=rng, data=data,
                                        anchor=anchor, settings=settings, budget=budget)
        print_crossing(settings, budget, "\nInitial contour anchor", anchor)
        print_fit_reports(data, settings)

        anchors = trace_from_anchor(backend=rmb.backend, rng=rng, data=data,
                                    anchor=anchor, settings=settings, budget=budget)
        print_progress(settings, budget, f"\nContour trace complete: {len(anchors)} anchors")
        print_fit_reports(data, settings)

    refinements = 0
    if settings.refine_after_trace and anchors and budget.remaining_hqc > 0.0:
        acquired = acquire_boundary_candidates(
            backend=rmb.backend, rng=rng, data=data, settings=settings,
            budget=budget, anchors=anchors)
        refinements += acquired
        if acquired > 0:
            print_progress(settings, budget,
                           f"\nBoundary acquisition complete: {acquired} extra circuits")

        confirmed = confirm_traced_anchors(
            backend=rmb.backend, rng=rng, data=data, anchors=anchors,
            settings=settings, budget=budget)
        refinements += confirmed
        if confirmed > 0:
            print_progress(settings, budget,
                           f"\nTrace anchor confirmation complete: {confirmed} extra circuits")

        extra = refine_uncertain_boundary_points(
            backend=rmb.backend, rng=rng, data=data, settings=settings, budget=budget)
        refinements += extra
        if extra > 0:
            print_progress(settings, budget, f"\nRefinement complete: {extra} extra circuits")

    if budget.remaining_hqc <= 0.0:
        stop_reason = "HQC budget exhausted"
    elif budget.remaining_hqc <= BASE_SIMULATION_COST:
        stop_reason = "HQC budget effectively exhausted: remaining HQC below the base cost"
    elif refinements > 0:
        stop_reason = "trace reached bounds, then refinement stopped with remaining budget"
    elif stop_reason is None:
        stop_reason = "contour trace reached the requested bounds before exhausting the budget"

    print_experiment_summary(data, settings, budget, stop_reason=stop_reason)
    print_fit_reports(data, settings)

    base_path = None
    if settings.save_path is not None:
        base_path = save_crossings(rmb, settings, budget, anchors)

    if settings.plot:
        from sympleq.applications.randomized_benchmarking.experiments.plots import plot_crossing_results
        plot_crossing_results(data, settings, anchors, base_path=base_path)

    return rmb, anchors


if __name__ == "__main__":
    rmb, anchors = run(CharlieCrossingSettings())
