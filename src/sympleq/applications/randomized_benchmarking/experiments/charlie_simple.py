"""
Simple hybrid crossing experiment.

This keeps the transparent ratio-by-ratio bisection structure of
``level_crossing.py`` and borrows only small, low-complexity ideas from
``charlie_crossing.py``:

* reserve budget by limiting HQC spent at later ratios, after the first
  crossing has bootstrapped the trace;
* use a fitted monotone surface, when available, as a bracket hint;
* spend a few local refinement probes around each found crossing.

It deliberately avoids full boundary acquisition, diversity scoring,
and multi-pass refinement machinery.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    print_crossing,
    print_experiment_summary,
    print_progress,
    save_crossings,
    single_circuit_bare_hqc,
    start_run,
    stitch_batch_size,
    stitched_batch_hqc,
    try_fit_monotone_fidelity_surface,
)
from sympleq.applications.randomized_benchmarking.experiments.monotone_tracing import (
    new_estimator,
    posterior_above,
    spend_measurements,
)


@dataclass(frozen=True)
class CharlieSimpleSettings(CrossingSettings):
    """
    Settings for a simple hybrid of level crossing and Charlie-style tracing.

    The algorithm is still a fixed-ratio bisection sweep. The extra knobs only
    constrain budget use and improve the next bracket when enough data exists.
    """
    ratio_step: float = 0.1
    start_ratio: float | None = None
    trace_direction: Literal["up", "down", "both"] = "down"
    n_gates_resolution: float = 0.1
    max_shots_per_config: int = 5
    decision_confidence: float = 0.8
    save_path: str | Path | None = "charlie_simple.json"
    one_q_noise_scale: float = 1.0
    two_q_noise_scale: float = 1.0

    # Optional protection against one difficult ratio consuming the whole
    # experiment after the trace has started. Disabled by default because the
    # plain level-crossing baseline sometimes legitimately spends more than a
    # small cap to resolve an early ratio.
    max_hqc_per_ratio: float | None = None

    # Small deterministic refinement stack around accepted crossings.
    refine_crossings: bool = False
    refinement_fraction: float = 0.06
    refinement_shots: int = 2

    # Use the shared monotone surface as a bracket hint after enough data exists.
    use_surface_bracket_hint: bool = True
    surface_min_configs: int = 8
    bracket_half_width_fraction: float = 0.25

    # High-gate endpoints are costly. Before any crossing has been found, start
    # every sweep with a cheap low-gate bracket and fall back to the full search
    # range only when that bracket fails.
    initial_bracket_fraction: float = 0.1
    high_ratio_gate_bias_power: float = 1.0

    # Use low-to-high endpoint discovery. This avoids direct jumps to very
    # large, confidently low-fidelity high endpoints, which otherwise create
    # horizontal bands of wasted low-fidelity measurements.
    low_to_high_growth_factor: float = 1.9

    # Directional tracing: when the ratio decreases, the boundary is expected
    # to move to more gates; when it increases, to fewer gates.
    trace_gate_growth: float = 1.5
    trace_gate_shrink: float = 0.75

    # Do not turn broad p ~= 0.5 vertical stacks into accepted crossings too
    # early. Uncertain points only become crossings when the bracket is narrow,
    # unless this explicit escape hatch is enabled.
    accept_uncertain_crossing: bool = True
    min_crossing_confidence: float = 0.65
    extra_midpoint_shots: int = 4


def implied_above(data: RMBData, config: RMBConfig,
                  settings: CharlieSimpleSettings) -> float | None:
    """Side of 0.5 implied by already confident monotone comparisons."""
    for other, estimator in data.items():
        if other.n_qubits != config.n_qubits:
            continue
        easier = (config.n_1qb_gates <= other.n_1qb_gates
                  and config.n_2qb_gates <= other.n_2qb_gates)
        harder = (config.n_1qb_gates >= other.n_1qb_gates
                  and config.n_2qb_gates >= other.n_2qb_gates)
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


def confidently_above(above: float, settings: CharlieSimpleSettings) -> bool:
    return above >= settings.decision_confidence


def confidently_below(above: float, settings: CharlieSimpleSettings) -> bool:
    return 1.0 - above >= settings.decision_confidence


def narrow_bracket(lo: int, hi: int, settings: CharlieSimpleSettings) -> bool:
    return hi - lo <= max(2, settings.n_gates_resolution * hi)


def accept_uncertain_crossing(above: float, settings: CharlieSimpleSettings) -> bool:
    return (
        settings.accept_uncertain_crossing
        and max(above, 1.0 - above) >= settings.min_crossing_confidence
    )


def resolved_start_ratio(settings: CharlieSimpleSettings) -> float:
    """Start ratio for the trace, defaulting by trace direction."""
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


def ratio_sweeps(settings: CharlieSimpleSettings) -> list[list[float]]:
    """One or two ratio sweeps requested by the trace settings."""
    low, high = settings.ratio_bounds
    start = resolved_start_ratio(settings)
    if settings.trace_direction == "up":
        return [ratio_sweep(start, high, settings.ratio_step)]
    if settings.trace_direction == "down":
        return [ratio_sweep(start, low, settings.ratio_step)]
    if settings.trace_direction == "both":
        upward = ratio_sweep(start, high, settings.ratio_step)
        downward_start = start - settings.ratio_step
        downward = ([] if downward_start < low - 1e-9
                    else ratio_sweep(downward_start, low, settings.ratio_step))
        return [upward, downward]
    raise ValueError(f"Unsupported trace_direction={settings.trace_direction!r}.")


def probe_fidelity(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    budget: Budget,
    ratio_budget: Budget,
    settings: CharlieSimpleSettings,
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
            stitch_batch_size(bare_hqc, settings.max_cost_per_run,
                              settings.max_shots_per_config),
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
    ratios = [c.ratio_2_qb_gates for c in crossings]
    inverse_sizes = [1.0 / c.n_gates for c in crossings]
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


def surface_prediction(data: RMBData, settings: CharlieSimpleSettings,
                       ratio: float) -> int | None:
    """Predict crossing gate count from the shared monotone surface, if usable."""
    if not settings.use_surface_bracket_hint:
        return None
    if measured_config_count(data) < settings.surface_min_configs:
        return None
    surface = try_fit_monotone_fidelity_surface(data, settings)
    if surface is None:
        return None

    gates_axis = np.linspace(settings.n_gates_bounds[0], settings.n_gates_bounds[1],
                             settings.candidate_grid_size[0])
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
    settings: CharlieSimpleSettings,
    ratio: float,
) -> tuple[int, int]:
    """Bracket the next fixed-ratio search from line and monotone hints."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    last = crossings[-1] if crossings else None
    predictions = [
        prediction for prediction in (
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
    data: RMBData,
    crossings: list[RMBConfig],
    settings: CharlieSimpleSettings,
    sweep: list[float],
) -> tuple[int, int]:
    """Initial bracket for a sweep, biased low until the trace has an anchor."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    if not sweep:
        return n_gates_min, n_gates_max
    if not crossings:
        hi = round(settings.initial_bracket_fraction * n_gates_max)
        return n_gates_min, min(n_gates_max, max(n_gates_min + 2, hi))
    return next_bracket(data, crossings, settings, sweep[0])


def ratio_biased_gate_cap(settings: CharlieSimpleSettings, ratio: float) -> int:
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
    low_gate_cap = max(n_gates_min + 2, round(settings.initial_bracket_fraction * n_gates_max))
    cap = (1.0 - expensive_corner_weight) * n_gates_max + expensive_corner_weight * low_gate_cap
    return int(min(n_gates_max, max(low_gate_cap, round(cap))))


def fallback_bracket(settings: CharlieSimpleSettings, ratio: float) -> tuple[int, int]:
    """Fallback bracket that avoids high-gate/high-ratio probes."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    return n_gates_min, min(n_gates_max, ratio_biased_gate_cap(settings, ratio))


def find_crossing(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    budget: Budget,
    ratio_budget: Budget,
    settings: CharlieSimpleSettings,
    ratio: float,
    lo: int,
    hi: int,
) -> RMBConfig | None:
    """Bisect in total gates for a crossing at fixed ratio."""
    config_lo = settings.make_config(lo, ratio)
    p_lo, above_lo = probe_fidelity(
        backend=backend, rng=rng, data=data, config=config_lo,
        budget=budget, ratio_budget=ratio_budget, settings=settings)
    if p_lo is None:
        return None
    if not confidently_above(above_lo, settings):
        if not confidently_below(above_lo, settings) and accept_uncertain_crossing(above_lo, settings):
            return config_lo
        return None

    previous = lo
    probe_hi = max(lo + 2, 2 * round(max(lo + 2, lo * settings.low_to_high_growth_factor) / 2))
    while probe_hi < hi:
        config_hi = settings.make_config(probe_hi, ratio)
        p_hi, above_hi = probe_fidelity(
            backend=backend, rng=rng, data=data, config=config_hi,
            budget=budget, ratio_budget=ratio_budget, settings=settings)
        if p_hi is None:
            return None
        print_progress(settings, budget,
                       f"  bracket ratio={ratio:.3f}: n_gates={probe_hi} p={p_hi:.2f}")
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
        probe_hi = max(probe_hi + 2, 2 * round((probe_hi * settings.low_to_high_growth_factor) / 2))
    else:
        lo = previous

    config_hi = settings.make_config(hi, ratio)
    p_hi, above_hi = probe_fidelity(
        backend=backend, rng=rng, data=data, config=config_hi,
        budget=budget, ratio_budget=ratio_budget, settings=settings)
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
            backend=backend, rng=rng, data=data, config=settings.make_config(mid, ratio),
            budget=budget, ratio_budget=ratio_budget, settings=settings)
        if p_mid is None:
            break
        print_progress(settings, budget, f"  bisect ratio={ratio:.3f}: n_gates={mid} p={p_mid:.2f}")
        if confidently_above(above_mid, settings):
            lo = mid
        elif confidently_below(above_mid, settings):
            hi = mid
        else:
            if settings.extra_midpoint_shots > 0:
                extra_settings = replace(
                    settings,
                    max_shots_per_config=settings.max_shots_per_config
                    + settings.extra_midpoint_shots,
                )
                p_mid, above_mid = probe_fidelity(
                    backend=backend, rng=rng, data=data,
                    config=settings.make_config(mid, ratio),
                    budget=budget, ratio_budget=ratio_budget,
                    settings=extra_settings)
                if p_mid is not None:
                    print_progress(
                        settings, budget,
                        f"  refine midpoint ratio={ratio:.3f}: n_gates={mid} p={p_mid:.2f}")
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


def refine_crossing(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    budget: Budget,
    ratio_budget: Budget,
    settings: CharlieSimpleSettings,
    crossing: RMBConfig,
) -> None:
    """Take a small local stack around an accepted crossing."""
    if not settings.refine_crossings or settings.refinement_shots <= 0:
        return
    ratio = crossing.ratio_2_qb_gates
    span = max(2, round(settings.refinement_fraction * crossing.n_gates))
    gates = [
        crossing.n_gates - span,
        crossing.n_gates,
        crossing.n_gates + span,
    ]
    local_settings = replace(settings, max_shots_per_config=settings.refinement_shots)
    for n_gates in gates:
        if budget.remaining_hqc <= 0 or ratio_budget.remaining_hqc <= 0:
            break
        if not settings.n_gates_bounds[0] <= n_gates <= settings.n_gates_bounds[1]:
            continue
        probe_fidelity(
            backend=backend, rng=rng, data=data,
            config=settings.make_config(n_gates, ratio),
            budget=budget, ratio_budget=ratio_budget, settings=local_settings)


def run_with_budget(settings: CharlieSimpleSettings) -> tuple[RMB, list[RMBConfig], Budget]:
    """
    Trace the fidelity = 0.5 line with a simple level/Charlie hybrid.
    """
    rng, rmb, budget = start_run(settings)
    data = rmb._data
    crossings: list[RMBConfig] = []

    n_gates_min, n_gates_max = settings.n_gates_bounds
    for sweep in ratio_sweeps(settings):
        if not sweep:
            continue
        lo, hi = initial_sweep_bracket(data, crossings, settings, sweep)
        if hi <= lo:
            lo, hi = n_gates_min, n_gates_max

        for ratio in sweep:
            if budget.remaining_hqc <= 0:
                break
            if crossings and settings.max_hqc_per_ratio is not None:
                ratio_hqc = min(settings.max_hqc_per_ratio, budget.remaining_hqc)
            else:
                ratio_hqc = budget.remaining_hqc
            ratio_budget = Budget(remaining_hqc=ratio_hqc)
            crossing = find_crossing(
                backend=rmb.backend, rng=rng, data=data, budget=budget,
                ratio_budget=ratio_budget, settings=settings, ratio=ratio, lo=lo, hi=hi)
            fallback_lo, fallback_hi = fallback_bracket(settings, ratio)
            if crossing is None and (lo, hi) != (fallback_lo, fallback_hi):
                lo, hi = fallback_lo, fallback_hi
                ratio_budget.remaining_hqc = min(ratio_budget.remaining_hqc,
                                                 budget.remaining_hqc)
                crossing = find_crossing(
                    backend=rmb.backend, rng=rng, data=data, budget=budget,
                    ratio_budget=ratio_budget, settings=settings, ratio=ratio, lo=lo, hi=hi)

            if crossing is not None:
                crossings.append(crossing)
                print_crossing(settings, budget, "Crossing", crossing)
                refine_crossing(
                    backend=rmb.backend, rng=rng, data=data, budget=budget,
                    ratio_budget=ratio_budget, settings=settings, crossing=crossing)
            else:
                reason = "budget exhausted" if budget.remaining_hqc <= 0 else "ratio budget or endpoints unresolved"
                print_progress(settings, budget, f"No crossing found at ratio={ratio:.3f} ({reason})")

            next_ratio = ratio + (settings.ratio_step if sweep[-1] >= sweep[0]
                                  else -settings.ratio_step)
            lo, hi = next_bracket(data, crossings, settings, next_ratio)
            if hi <= lo:
                lo, hi = n_gates_min, n_gates_max

    print_progress(settings, budget, f"\nTraced {len(crossings)} crossings")
    if settings.verbose:
        for config in crossings:
            print(f"  ratio={config.ratio_2_qb_gates:.3f} n_gates={config.n_gates:>6}")
    print_experiment_summary(data, settings, budget)

    base_path = None
    if settings.save_path is not None:
        base_path = save_crossings(rmb, settings, budget, crossings)

    if settings.plot:
        from sympleq.applications.randomized_benchmarking.experiments.plots import plot_crossing_results
        plot_crossing_results(data, settings, crossings, base_path=base_path)

    return rmb, crossings, budget


def run(settings: CharlieSimpleSettings) -> tuple[RMB, list[RMBConfig]]:
    """Trace the fidelity = 0.5 line and return the run data and crossings."""
    rmb, crossings, _ = run_with_budget(settings)
    return rmb, crossings


if __name__ == "__main__":
    rmb, crossings = run(CharlieSimpleSettings(rng_seed=None))
