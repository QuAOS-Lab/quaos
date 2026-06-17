"""
Deterministic monotone tracing of the fidelity = 0.5 boundary.

This is the cleaned-up standalone version of the current ``mc_tracing.py``
behaviour. It uses the local SympleQ backend, finds an initial contour point
with a cheap low-to-high gate search, and then traces neighbouring ratios with
deterministic brackets, bisection, a monotone surface hint, and an explicit
local monotonicity constraint.
"""
from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    new_estimator,
    posterior_above,
    print_crossing,
    print_experiment_summary,
    print_progress,
    save_crossings,
    single_circuit_bare_hqc,
    spend_measurements,
    start_run,
    stitch_batch_size,
    stitched_batch_hqc,
    try_fit_monotone_fidelity_surface,
)
from sympleq.core.noise.noise_model import GenericNoise


BASE_1Q_PAULI_ERROR = 0.000025
BASE_2Q_PAULI_ERROR = 0.00079


def sympleq_backend_factory(
    settings: CrossingSettings,
    rng: RNGGenerator,
) -> RMBBackend:
    """Local SympleQ backend with the same default noise as the other demos."""
    noise_model = GenericNoise.from_paulis([BASE_1Q_PAULI_ERROR] * 3, rng)
    two_qubit_noise_model = GenericNoise.from_paulis([BASE_2Q_PAULI_ERROR] * 3, rng)
    return SympleqBackend(
        noise_model=noise_model,
        two_qubit_noise_model=two_qubit_noise_model,
    )


@dataclass(frozen=True)
class MonotoneTracingSettings(CrossingSettings):
    """
    Settings for deterministic local contour tracing.

    The contour is traced in fixed-ratio slices. Each slice searches from low
    gate counts upward until it brackets the fidelity = 0.5 crossing, then
    bisects that bracket. Previous crossings and the fitted monotone surface
    guide the next bracket, while a local monotonicity constraint rejects
    crossings that clearly move in the wrong direction.
    """

    ratio_step: float = 0.15
    start_ratio: float | None = None
    trace_direction: Literal["up", "down", "both"] = "down"
    save_path: str | Path | None = "monotone_tracing.json"
    backend_factory: Callable[[CrossingSettings, RNGGenerator], RMBBackend] = (
        sympleq_backend_factory
    )

    adaptive_ratio_step: bool = False
    min_ratio_step: float = 0.04
    max_ratio_step: float = 0.2
    target_hqc_per_ratio: float = 35.0

    n_gates_resolution: float = 0.08
    max_shots_per_config: int = 7
    decision_confidence: float = 0.75
    max_hqc_per_ratio: float | None = 45.0
    accept_uncertain_crossing: bool = True
    min_crossing_confidence: float = 0.6
    extra_midpoint_shots: int = 4

    initial_bracket_fraction: float = 0.1
    high_ratio_gate_bias_power: float = 1.0
    low_to_high_growth_factor: float = 1.7
    bracket_half_width_fraction: float = 0.24
    trace_gate_growth: float = 1.7
    trace_gate_shrink: float = 0.6
    use_surface_bracket_hint: bool = True
    surface_min_configs: int = 8

    enforce_ratio_monotonicity: bool = True
    monotonic_min_slack_gates: int = 2
    monotonic_slack_fraction: float | None = 0.5


def implied_above(
    data: RMBData,
    config: RMBConfig,
    settings: MonotoneTracingSettings,
) -> float | None:
    """Infer a side of 0.5 from already confident monotone comparisons."""
    for other, estimator in data.items():
        if other.n_qubits != config.n_qubits:
            continue
        easier = (
            config.n_1qb_gates <= other.n_1qb_gates
            and config.n_2qb_gates <= other.n_2qb_gates
        )
        harder = (
            config.n_1qb_gates >= other.n_1qb_gates
            and config.n_2qb_gates >= other.n_2qb_gates
        )
        if not easier and not harder:
            continue

        above = posterior_above(estimator)
        if easier and above >= settings.decision_confidence:
            return 1.0
        if harder and 1.0 - above >= settings.decision_confidence:
            return 0.0
    return None


def measured_config_count(data: RMBData) -> int:
    return sum(1 for estimator in data.values() if estimator.num_runs() > 0)


def spendable_ratio_budget(global_budget: Budget, ratio_budget: Budget) -> float:
    return max(0.0, min(global_budget.remaining_hqc, ratio_budget.remaining_hqc))


def confidently_above(above: float, settings: MonotoneTracingSettings) -> bool:
    return above >= settings.decision_confidence


def confidently_below(above: float, settings: MonotoneTracingSettings) -> bool:
    return 1.0 - above >= settings.decision_confidence


def narrow_bracket(lo: int, hi: int, settings: MonotoneTracingSettings) -> bool:
    return hi - lo <= max(2, settings.n_gates_resolution * hi)


def accept_uncertain_crossing(
    above: float,
    settings: MonotoneTracingSettings,
) -> bool:
    return (
        settings.accept_uncertain_crossing
        and max(above, 1.0 - above) >= settings.min_crossing_confidence
    )


def resolved_start_ratio(settings: MonotoneTracingSettings) -> float:
    low, high = settings.ratio_bounds
    if settings.start_ratio is None:
        if settings.trace_direction == "down":
            return high
        if settings.trace_direction == "both":
            return 0.5 * (low + high)
        return low
    if not low <= settings.start_ratio <= high:
        raise ValueError(
            f"start_ratio={settings.start_ratio} outside ratio_bounds={settings.ratio_bounds}."
        )
    return settings.start_ratio


def ratio_sweep(start: float, stop: float, step: float) -> list[float]:
    """Inclusive ratio sweep from start to stop using the sign implied by stop."""
    if step <= 0:
        raise ValueError(f"ratio_step must be positive, got {step}.")
    direction = 1.0 if stop >= start else -1.0
    ratios: list[float] = []
    current = start
    while direction * (current - stop) <= 1e-9:
        ratios.append(float(current))
        current += direction * step
    if ratios and abs(ratios[-1] - stop) > 1e-9:
        ratios.append(float(stop))
    return ratios


def adaptive_ratio_sweep(
    start: float,
    stop: float,
    settings: MonotoneTracingSettings,
    budget: Budget,
) -> Iterator[float]:
    """Yield ratios, widening steps only when explicitly requested."""
    if settings.ratio_step <= 0:
        raise ValueError(f"ratio_step must be positive, got {settings.ratio_step}.")
    if settings.min_ratio_step <= 0 or settings.max_ratio_step <= 0:
        raise ValueError("adaptive ratio step bounds must be positive.")
    if settings.max_ratio_step < settings.min_ratio_step:
        raise ValueError("max_ratio_step must be at least min_ratio_step.")

    direction = 1.0 if stop >= start else -1.0
    current = start
    while direction * (current - stop) <= 1e-9:
        yield float(current)
        if abs(current - stop) <= 1e-9:
            break

        if settings.adaptive_ratio_step:
            remaining_span = abs(stop - current)
            expected_per_ratio = max(1.0, settings.target_hqc_per_ratio)
            affordable_remaining = max(1, int(budget.remaining_hqc // expected_per_ratio))
            step = remaining_span / affordable_remaining
            step = float(np.clip(step, settings.min_ratio_step, settings.max_ratio_step))
        else:
            step = settings.ratio_step

        next_ratio = current + direction * step
        if direction * (next_ratio - stop) > 0.0:
            next_ratio = stop
        if abs(next_ratio - current) <= 1e-12:
            break
        current = next_ratio


def ratio_sweeps(settings: MonotoneTracingSettings) -> list[list[float]]:
    low, high = settings.ratio_bounds
    start = resolved_start_ratio(settings)
    if settings.trace_direction == "up":
        return [ratio_sweep(start, high, settings.ratio_step)]
    if settings.trace_direction == "down":
        return [ratio_sweep(start, low, settings.ratio_step)]
    if settings.trace_direction == "both":
        upward = ratio_sweep(start, high, settings.ratio_step)
        downward_start = start - settings.ratio_step
        downward = (
            []
            if downward_start < low - 1e-9
            else ratio_sweep(downward_start, low, settings.ratio_step)
        )
        return [upward, downward]
    raise ValueError(f"Unsupported trace_direction={settings.trace_direction!r}.")


def ratio_sweep_iterators(
    settings: MonotoneTracingSettings,
    budget: Budget,
) -> Iterator[Iterator[float] | list[float]]:
    if not settings.adaptive_ratio_step:
        yield from ratio_sweeps(settings)
        return

    low, high = settings.ratio_bounds
    start = resolved_start_ratio(settings)
    if settings.trace_direction == "up":
        yield adaptive_ratio_sweep(start, high, settings, budget)
    elif settings.trace_direction == "down":
        yield adaptive_ratio_sweep(start, low, settings, budget)
    elif settings.trace_direction == "both":
        yield adaptive_ratio_sweep(start, high, settings, budget)
        downward_start = start - settings.ratio_step
        if downward_start >= low - 1e-9:
            yield adaptive_ratio_sweep(downward_start, low, settings, budget)
    else:
        raise ValueError(f"Unsupported trace_direction={settings.trace_direction!r}.")


def ratio_budget_for(
    settings: MonotoneTracingSettings,
    budget: Budget,
    crossings: list[RMBConfig],
) -> Budget:
    """Cap later ratios without restricting the first anchor search."""
    if crossings and settings.max_hqc_per_ratio is not None:
        return Budget(remaining_hqc=min(settings.max_hqc_per_ratio, budget.remaining_hqc))
    return Budget(remaining_hqc=budget.remaining_hqc)


def probe_fidelity(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    budget: Budget,
    ratio_budget: Budget,
    settings: MonotoneTracingSettings,
) -> tuple[float | None, float]:
    """
    Probe one config, bounded by both the global budget and this ratio's budget.
    """
    estimator = data.get(config, new_estimator())
    above = posterior_above(estimator)
    if max(above, 1.0 - above) < settings.decision_confidence:
        implied = implied_above(data, config, settings)
        if implied is not None:
            return estimator.posterior_mean(), implied

    bare_hqc = single_circuit_bare_hqc(config)
    while estimator.num_runs() < settings.max_shots_per_config:
        above = posterior_above(estimator)
        if max(above, 1.0 - above) >= settings.decision_confidence:
            break

        n_circuits = min(
            stitch_batch_size(
                bare_hqc,
                settings.max_cost_per_run,
                settings.max_shots_per_config,
            ),
            settings.max_shots_per_config - estimator.num_runs(),
        )
        while n_circuits > 0:
            cost = stitched_batch_hqc(n_circuits * bare_hqc)
            if cost <= spendable_ratio_budget(budget, ratio_budget):
                break
            n_circuits -= 1
        if n_circuits <= 0:
            break

        data[config] = estimator
        spend_measurements(backend, rng, data, config, n_circuits, seed=settings.rng_seed)
        cost = stitched_batch_hqc(n_circuits * bare_hqc)
        budget.spend_batch(cost, n_circuits)
        ratio_budget.spend_batch(cost, n_circuits)

    above = posterior_above(estimator)
    decided = max(above, 1.0 - above) >= settings.decision_confidence
    if not decided and estimator.num_runs() < settings.max_shots_per_config:
        return None, above
    return estimator.posterior_mean(), above


def inverse_line_fit(crossings: list[RMBConfig]) -> tuple[float, float] | None:
    """Fit 1/n*(ratio) = intercept + slope * ratio from previous crossings."""
    if len(crossings) < 2:
        return None
    ratios = [crossing.ratio_2_qb_gates for crossing in crossings]
    inverse_sizes = [1.0 / crossing.n_gates for crossing in crossings]
    slope, intercept = np.polyfit(ratios, inverse_sizes, 1)
    return float(slope), float(intercept)


def line_fit_prediction(crossings: list[RMBConfig], ratio: float) -> int | None:
    fit = inverse_line_fit(crossings)
    if fit is None:
        return None
    slope, intercept = fit
    inverse = intercept + slope * ratio
    if inverse <= 0:
        return None
    return round(1.0 / inverse)


def surface_prediction(
    data: RMBData,
    settings: MonotoneTracingSettings,
    ratio: float,
) -> int | None:
    """Predict crossing gate count from the shared monotone surface, if usable."""
    if not settings.use_surface_bracket_hint:
        return None
    if measured_config_count(data) < settings.surface_min_configs:
        return None
    surface = try_fit_monotone_fidelity_surface(data, settings)
    if surface is None:
        return None

    gates_axis = np.linspace(
        settings.n_gates_bounds[0],
        settings.n_gates_bounds[1],
        settings.candidate_grid_size[0],
    )
    points = np.column_stack([gates_axis, np.full_like(gates_axis, ratio)])
    delta = surface.probability(points) - 0.5
    crossing = np.where(delta[:-1] * delta[1:] <= 0)[0]
    if len(crossing) == 0:
        return None
    i = int(crossing[0])
    denom = abs(delta[i]) + abs(delta[i + 1])
    if denom <= 0:
        return round(float(gates_axis[i]))
    t = abs(delta[i]) / denom
    return round(float((1.0 - t) * gates_axis[i] + t * gates_axis[i + 1]))


def next_bracket(
    data: RMBData,
    crossings: list[RMBConfig],
    settings: MonotoneTracingSettings,
    ratio: float,
) -> tuple[int, int]:
    """Bracket the next fixed-ratio search from line and monotone hints."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    last = crossings[-1] if crossings else None
    predictions = [
        prediction
        for prediction in (
            line_fit_prediction(crossings, ratio),
            surface_prediction(data, settings, ratio),
        )
        if prediction is not None
    ]
    if predictions:
        center = int(np.median(predictions))
        width = max(
            int(settings.bracket_half_width_fraction * max(center, 1)),
            max(abs(center - prediction) for prediction in predictions),
            2,
        )
        lo = center - width
        hi = center + width
        if last is not None:
            if ratio < last.ratio_2_qb_gates:
                lo = min(lo, round(last.n_gates * settings.trace_gate_shrink))
                hi = max(hi, round(last.n_gates * settings.trace_gate_growth))
            elif ratio > last.ratio_2_qb_gates:
                lo = min(lo, round(last.n_gates / settings.trace_gate_growth))
                hi = max(hi, round(last.n_gates / settings.trace_gate_shrink))
        return max(n_gates_min, lo), min(n_gates_max, hi)

    if last is not None:
        if ratio < last.ratio_2_qb_gates:
            hi = round(last.n_gates * settings.trace_gate_growth)
        elif ratio > last.ratio_2_qb_gates:
            hi = round(last.n_gates / settings.trace_gate_shrink)
        else:
            hi = last.n_gates
        return n_gates_min, min(n_gates_max, max(n_gates_min + 2, hi))
    return n_gates_min, n_gates_max


def initial_sweep_bracket(
    settings: MonotoneTracingSettings,
    sweep: list[float],
) -> tuple[int, int]:
    """Initial bracket for a sweep, biased low until the trace has an anchor."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    if not sweep:
        return n_gates_min, n_gates_max
    hi = round(settings.initial_bracket_fraction * n_gates_max)
    return n_gates_min, min(n_gates_max, max(n_gates_min + 2, hi))


def ratio_biased_gate_cap(settings: MonotoneTracingSettings, ratio: float) -> int:
    """
    Maximum total-gate endpoint to try at this ratio.

    High-ratio circuits are expensive because most gates are two-qubit gates, so
    avoid jumping to the global high-gate bound there. At low ratio this returns
    the full upper bound; at high ratio it approaches the cheap initial bracket.
    """
    n_gates_min, n_gates_max = settings.n_gates_bounds
    low_ratio, high_ratio = settings.ratio_bounds
    span = max(high_ratio - low_ratio, 1e-12)
    normalized_ratio = np.clip((ratio - low_ratio) / span, 0.0, 1.0)
    expensive_corner_weight = normalized_ratio ** settings.high_ratio_gate_bias_power
    low_gate_cap = max(
        n_gates_min + 2,
        round(settings.initial_bracket_fraction * n_gates_max),
    )
    cap = (
        (1.0 - expensive_corner_weight) * n_gates_max
        + expensive_corner_weight * low_gate_cap
    )
    return int(min(n_gates_max, max(low_gate_cap, round(cap))))


def fallback_bracket(
    settings: MonotoneTracingSettings,
    ratio: float,
) -> tuple[int, int]:
    """Fallback bracket that avoids high-gate/high-ratio probes."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    return n_gates_min, min(n_gates_max, ratio_biased_gate_cap(settings, ratio))


def monotonic_slack(settings: MonotoneTracingSettings, anchor: RMBConfig) -> int:
    fraction = (
        settings.n_gates_resolution
        if settings.monotonic_slack_fraction is None
        else settings.monotonic_slack_fraction
    )
    return max(settings.monotonic_min_slack_gates, round(fraction * anchor.n_gates))


def monotonic_gate_bounds(
    settings: MonotoneTracingSettings,
    crossings: list[RMBConfig],
    ratio: float,
) -> tuple[int, int] | None:
    """
    Gate bounds implied by adjacent accepted crossings.

    The expected contour is monotone decreasing as two-qubit ratio increases.
    Equivalently, moving down in ratio should move to the same or larger gate
    count. Only the nearest accepted crossing on each side is used, with a
    slack band, so one distant noisy anchor cannot overconstrain the trace.
    """
    if not settings.enforce_ratio_monotonicity or not crossings:
        return None

    n_gates_min, n_gates_max = settings.n_gates_bounds
    lower = n_gates_min
    upper = n_gates_max
    nearest_higher_ratio = min(
        (anchor for anchor in crossings if anchor.ratio_2_qb_gates > ratio),
        key=lambda anchor: anchor.ratio_2_qb_gates - ratio,
        default=None,
    )
    nearest_lower_ratio = min(
        (anchor for anchor in crossings if anchor.ratio_2_qb_gates < ratio),
        key=lambda anchor: ratio - anchor.ratio_2_qb_gates,
        default=None,
    )

    if nearest_higher_ratio is not None:
        slack = monotonic_slack(settings, nearest_higher_ratio)
        lower = max(
            lower,
            nearest_higher_ratio.n_gates - slack - settings.monotonic_min_slack_gates,
        )
    if nearest_lower_ratio is not None:
        slack = monotonic_slack(settings, nearest_lower_ratio)
        upper = min(
            upper,
            nearest_lower_ratio.n_gates + slack + settings.monotonic_min_slack_gates,
        )

    if lower > upper:
        return None
    return lower, upper


def apply_monotonic_bounds(
    settings: MonotoneTracingSettings,
    crossings: list[RMBConfig],
    ratio: float,
    lo: int,
    hi: int,
) -> tuple[int, int]:
    bounds = monotonic_gate_bounds(settings, crossings, ratio)
    if bounds is None:
        return lo, hi
    lower, upper = bounds
    return max(lo, lower), min(hi, upper)


def monotonicity_violation_message(
    settings: MonotoneTracingSettings,
    crossings: list[RMBConfig],
    crossing: RMBConfig,
) -> str | None:
    bounds = monotonic_gate_bounds(settings, crossings, crossing.ratio_2_qb_gates)
    if bounds is None:
        return None
    lower, upper = bounds
    if lower <= crossing.n_gates <= upper:
        return None
    return (
        f"rejected non-monotone crossing at ratio={crossing.ratio_2_qb_gates:.3f}: "
        f"n_gates={crossing.n_gates} outside allowed [{lower}, {upper}]"
    )


def find_crossing(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    budget: Budget,
    ratio_budget: Budget,
    settings: MonotoneTracingSettings,
    ratio: float,
    lo: int,
    hi: int,
) -> RMBConfig | None:
    """Search low-to-high, then bisect in total gates at fixed ratio."""
    config_lo = settings.make_config(lo, ratio)
    p_lo, above_lo = probe_fidelity(
        backend=backend,
        rng=rng,
        data=data,
        config=config_lo,
        budget=budget,
        ratio_budget=ratio_budget,
        settings=settings,
    )
    if p_lo is None:
        return None
    if not confidently_above(above_lo, settings):
        if (
            not confidently_below(above_lo, settings)
            and accept_uncertain_crossing(above_lo, settings)
        ):
            return config_lo
        return None

    previous = lo
    probe_hi = max(
        lo + 2,
        2 * round(max(lo + 2, lo * settings.low_to_high_growth_factor) / 2),
    )
    while probe_hi < hi:
        config_hi = settings.make_config(probe_hi, ratio)
        p_hi, above_hi = probe_fidelity(
            backend=backend,
            rng=rng,
            data=data,
            config=config_hi,
            budget=budget,
            ratio_budget=ratio_budget,
            settings=settings,
        )
        if p_hi is None:
            return None
        print_progress(
            settings,
            budget,
            f"  bracket ratio={ratio:.3f}: n_gates={probe_hi} p={p_hi:.2f}",
        )
        if not confidently_above(above_hi, settings):
            if not confidently_below(above_hi, settings):
                if accept_uncertain_crossing(above_hi, settings):
                    return config_hi
                previous = probe_hi
                probe_hi = max(
                    probe_hi + 2,
                    2 * round((probe_hi * settings.low_to_high_growth_factor) / 2),
                )
                continue
            lo = previous
            hi = probe_hi
            break
        previous = probe_hi
        probe_hi = max(
            probe_hi + 2,
            2 * round((probe_hi * settings.low_to_high_growth_factor) / 2),
        )
    else:
        lo = previous

    config_hi = settings.make_config(hi, ratio)
    p_hi, above_hi = probe_fidelity(
        backend=backend,
        rng=rng,
        data=data,
        config=config_hi,
        budget=budget,
        ratio_budget=ratio_budget,
        settings=settings,
    )
    if p_hi is None:
        return None
    if not confidently_below(above_hi, settings):
        if not confidently_above(above_hi, settings) and (
            narrow_bracket(lo, hi, settings)
            or accept_uncertain_crossing(above_hi, settings)
        ):
            return config_hi
        return None

    while not narrow_bracket(lo, hi, settings):
        mid = 2 * round((lo + hi) / 4)
        if mid in (lo, hi):
            break
        p_mid, above_mid = probe_fidelity(
            backend=backend,
            rng=rng,
            data=data,
            config=settings.make_config(mid, ratio),
            budget=budget,
            ratio_budget=ratio_budget,
            settings=settings,
        )
        if p_mid is None:
            break
        print_progress(
            settings,
            budget,
            f"  bisect ratio={ratio:.3f}: n_gates={mid} p={p_mid:.2f}",
        )
        if confidently_above(above_mid, settings):
            lo = mid
            continue
        if confidently_below(above_mid, settings):
            hi = mid
            continue

        if settings.extra_midpoint_shots > 0:
            extra_settings = replace(
                settings,
                max_shots_per_config=(
                    settings.max_shots_per_config + settings.extra_midpoint_shots
                ),
            )
            p_mid, above_mid = probe_fidelity(
                backend=backend,
                rng=rng,
                data=data,
                config=settings.make_config(mid, ratio),
                budget=budget,
                ratio_budget=ratio_budget,
                settings=extra_settings,
            )
            if p_mid is not None:
                print_progress(
                    settings,
                    budget,
                    f"  refine midpoint ratio={ratio:.3f}: n_gates={mid} p={p_mid:.2f}",
                )

        if confidently_above(above_mid, settings):
            lo = mid
        elif confidently_below(above_mid, settings):
            hi = mid
        elif accept_uncertain_crossing(above_mid, settings):
            return settings.make_config(mid, ratio)
        else:
            break

    if not narrow_bracket(lo, hi, settings):
        return None
    return settings.make_config((lo + hi) // 2, ratio)


def search_one_ratio(
    rmb: RMB,
    rng: RNGGenerator,
    data: RMBData,
    crossings: list[RMBConfig],
    budget: Budget,
    settings: MonotoneTracingSettings,
    ratio: float,
) -> RMBConfig | None:
    """Find one fixed-ratio crossing with deterministic bracketing."""
    ratio_budget = ratio_budget_for(settings, budget, crossings)
    if crossings:
        lo, hi = next_bracket(data, crossings, settings, ratio)
        lo, hi = apply_monotonic_bounds(settings, crossings, ratio, lo, hi)
        if hi <= lo:
            bounds = monotonic_gate_bounds(settings, crossings, ratio)
            lo, hi = bounds if bounds is not None else settings.n_gates_bounds
        if hi <= lo:
            print_progress(
                settings,
                budget,
                f"No monotone search interval at ratio={ratio:.3f} "
                f"(allowed n_gates {lo} to {hi})",
            )
            return None
    else:
        lo, hi = initial_sweep_bracket(settings, [ratio])

    crossing = find_crossing(
        backend=rmb.backend,
        rng=rng,
        data=data,
        budget=budget,
        ratio_budget=ratio_budget,
        settings=settings,
        ratio=ratio,
        lo=lo,
        hi=hi,
    )

    fallback_lo, fallback_hi = fallback_bracket(settings, ratio)
    fallback_lo, fallback_hi = apply_monotonic_bounds(
        settings,
        crossings,
        ratio,
        fallback_lo,
        fallback_hi,
    )
    if crossing is None and fallback_hi > fallback_lo and (lo, hi) != (
        fallback_lo,
        fallback_hi,
    ):
        ratio_budget.remaining_hqc = min(ratio_budget.remaining_hqc, budget.remaining_hqc)
        crossing = find_crossing(
            backend=rmb.backend,
            rng=rng,
            data=data,
            budget=budget,
            ratio_budget=ratio_budget,
            settings=settings,
            ratio=ratio,
            lo=fallback_lo,
            hi=fallback_hi,
        )

    if crossing is not None:
        message = monotonicity_violation_message(settings, crossings, crossing)
        if message is not None:
            print_progress(settings, budget, message)
            return None
    return crossing


def run_with_budget(
    settings: MonotoneTracingSettings,
) -> tuple[RMB, list[RMBConfig], Budget]:
    """Run the deterministic local boundary tracer."""
    rng, rmb, budget = start_run(settings)
    if not isinstance(rmb.backend, SympleqBackend):
        raise TypeError(
            "MonotoneTracingSettings must use SympleqBackend; refusing non-local backend."
        )

    data = rmb._data
    crossings: list[RMBConfig] = []

    for sweep in ratio_sweep_iterators(settings, budget):
        for ratio in sweep:
            if budget.remaining_hqc <= 0:
                break
            crossing = search_one_ratio(
                rmb,
                rng,
                data,
                crossings,
                budget,
                settings,
                ratio,
            )
            if crossing is None:
                reason = (
                    "budget exhausted"
                    if budget.remaining_hqc <= 0
                    else "no local boundary found"
                )
                print_progress(
                    settings,
                    budget,
                    f"No crossing found at ratio={ratio:.3f} ({reason})",
                )
                continue
            crossings.append(crossing)
            print_crossing(settings, budget, "Crossing", crossing)

    print_progress(settings, budget, f"\nTraced {len(crossings)} crossings")
    if settings.verbose:
        for config in crossings:
            print(f"  ratio={config.ratio_2_qb_gates:.3f} n_gates={config.n_gates:>6}")
    print_experiment_summary(data, settings, budget)

    base_path = None
    if settings.save_path is not None:
        base_path = save_crossings(rmb, settings, budget, crossings)

    if settings.plot:
        from sympleq.applications.randomized_benchmarking.experiments.plots import (
            plot_crossing_results,
        )

        plot_crossing_results(data, settings, crossings, base_path=base_path)

    return rmb, crossings, budget


def run(settings: MonotoneTracingSettings) -> tuple[RMB, list[RMBConfig]]:
    """Run the deterministic local boundary tracer and return data plus crossings."""
    rmb, crossings, _ = run_with_budget(settings)
    return rmb, crossings


if __name__ == "__main__":
    rmb, crossings = run(MonotoneTracingSettings(rng_seed=None))
