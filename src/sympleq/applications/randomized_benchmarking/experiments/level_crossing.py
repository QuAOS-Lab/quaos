"""
Trace the fidelity = 0.5 line in (total gates, two-qubit gate ratio) space at
fixed n_qubits.

The experiment finds a first crossing by starting at the high two-qubit ratio
and growing upward from a minimal gate count until the contour is bracketed.
It then tracks the line top-down in ratio, reusing previous crossings to seed
small low-to-high search windows instead of probing the global max gate count.
Configs are built from the transformed parameters on demand. Every probe is
metered in Quantinuum credits (HQC) via ``pytket_bare_simulation_cost`` with
the base submission cost paid once per stitched batch, and the whole run stops
when the HQC budget is exhausted.

Circuits run on the backend built by ``settings.backend_factory``: by default
the SympleQ emulation of Quantinuum hardware, so nothing is submitted to
Quantinuum.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.base import MeasurementRequest
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    print_crossing,
    print_experiment_summary,
    print_progress,
    save_crossings,
    single_circuit_bare_hqc,
    spend_request_batch,
    start_run,
    stitch_batch_size,
    stitched_batch_hqc,
)


@dataclass(frozen=True)
class LevelCrossingSettings(CrossingSettings):
    """
    Settings for a fidelity = 0.5 level-crossing experiment.

    The problem definition (qubits, bounds, budget, ...) lives in
    :class:`~.common.CrossingSettings`; by default crossings are traced from
    ``ratio_bounds[1]`` down to ``ratio_bounds[0]``.

    Parameters
    ----------
    ratio_step : float
        Increment of the ratio between consecutive traced crossings.
    n_gates_resolution : float
        Bisection stops once the bracket is narrower than this fraction of
        its upper edge.
    max_shots_per_config : int
        Hard cap on circuits spent on a single config probe.
    decision_confidence : float
        Posterior mass required on one side of 0.5 to end a probe early.
        Values at or below 0.875 let one or two unanimous shots decide a
        side, which makes bracket checks unreliable; keep it above that.
    """
    ratio_step: float = 0.1
    start_ratio: float | None = None
    trace_direction: Literal["up", "down"] = "down"
    n_gates_resolution: float = 0.05
    max_shots_per_config: int = 24
    decision_confidence: float = 0.975
    save_path: str | Path | None = "level_crossing.json"

    # Avoid immediately probing the expensive high-gate endpoint. Each ratio
    # starts from a low gate count and grows upward until the crossing is
    # bracketed or this capped search window is exhausted.
    initial_bracket_fraction: float = 0.1
    low_to_high_growth_factor: float = 1.9
    trace_gate_growth: float = 1.5
    trace_gate_shrink: float = 0.75


def implied_above(data: RMBData, config: RMBConfig,
                  settings: LevelCrossingSettings) -> float | None:
    """
    Side of 0.5 implied for ``config`` by monotonicity of already-decided data.

    Fidelity decreases when gates of either kind are added, so a config with
    at least as many one- and two-qubit gates as a confidently-below config
    is below, and one with at most as many of each as a confidently-above
    config is above. Such configs need no measurements at all.

    Returns
    -------
    float | None
        1.0 (above), 0.0 (below), or ``None`` when existing data does not
        determine the side.
    """
    for other, estimator in data.items():
        if other.n_qubits != config.n_qubits:
            continue
        # The cheap gate-count comparisons gate the Beta-posterior evaluation.
        easier = (config.n_1qb_gates <= other.n_1qb_gates
                  and config.n_2qb_gates <= other.n_2qb_gates)
        harder = (config.n_1qb_gates >= other.n_1qb_gates
                  and config.n_2qb_gates >= other.n_2qb_gates)
        if not easier and not harder:
            continue
        above = estimator.posterior_above()
        if easier and above >= settings.decision_confidence:
            return 1.0
        if harder and 1.0 - above >= settings.decision_confidence:
            return 0.0
    return None


def probe_fidelity(*, backend, rng: RNGGenerator, data: RMBData, config: RMBConfig,
                   budget: Budget, settings: LevelCrossingSettings) -> tuple[float | None, float]:
    """
    Estimate the fidelity of ``config`` with as few circuits as possible.

    Records circuit outcomes into ``data``, in stitched batches priced like
    the Quantinuum backend, until the Beta posterior places
    ``decision_confidence`` mass on one side of 0.5, the shot cap is reached,
    or the budget runs out. Configs whose side is already implied by
    monotonicity from existing data are not measured at all.

    Returns
    -------
    tuple[float | None, float]
        Posterior mean fidelity and posterior probability that the fidelity
        is above 0.5. The mean is ``None`` when the budget stopped the probe
        before it could either decide a side or reach the shot cap.
    """
    estimator = data.get(config, BayesianEstimator.default())
    above = estimator.posterior_above()
    if max(above, 1.0 - above) < settings.decision_confidence:
        implied = implied_above(data, config, settings)
        if implied is not None:
            return estimator.posterior_mean(), implied

    bare_hqc = single_circuit_bare_hqc(config)

    while estimator.num_runs() < settings.max_shots_per_config:
        above = estimator.posterior_above()
        if max(above, 1.0 - above) >= settings.decision_confidence:
            break
        n_circuits = min(
            stitch_batch_size(bare_hqc, settings.max_cost_per_run, settings.max_shots_per_config),
            settings.max_shots_per_config - estimator.num_runs())
        while n_circuits > 0 and not budget.can_afford(stitched_batch_hqc(n_circuits * bare_hqc)):
            n_circuits -= 1
        if n_circuits <= 0:
            break
        # Only configs with recorded outcomes enter the data set; implied or
        # unaffordable probes must not leave empty estimators behind.
        data[config] = estimator
        spend_request_batch(backend, rng, data,
                            [MeasurementRequest(config, n_circuits)], seed=settings.rng_seed)
        budget.spend_batch(stitched_batch_hqc(n_circuits * bare_hqc), n_circuits)

    above = estimator.posterior_above()
    decided = max(above, 1.0 - above) >= settings.decision_confidence
    if not decided and estimator.num_runs() < settings.max_shots_per_config:
        return None, above
    return estimator.posterior_mean(), above


def inverse_line_fit(crossings: list[RMBConfig]) -> tuple[float, float] | None:
    """
    Fit the level line as 1/n*(r) = intercept + slope * r.

    At the crossing the accumulated noise is roughly constant, so the inverse
    crossing size is approximately linear in the two-qubit gate ratio. Only
    used to predict the bisection bracket of the next traced crossing; the
    shared plots overlay the monotone-surface contour instead.

    Returns
    -------
    tuple[float, float] | None
        ``(slope, intercept)`` of the inverse-size fit, or ``None`` with
        fewer than two crossings.
    """
    if len(crossings) < 2:
        return None
    ratios = [c.ratio_2_qb_gates for c in crossings]
    inverse_sizes = [1.0 / c.n_gates for c in crossings]
    slope, intercept = np.polyfit(ratios, inverse_sizes, 1)
    return float(slope), float(intercept)


def predict_crossing(crossings: list[RMBConfig], ratio: float) -> int | None:
    """
    Predict the crossing total gates at ``ratio`` from previous crossings.

    Returns
    -------
    int | None
        Predicted total gates, or ``None`` with fewer than two crossings or
        when the fit does not cross at this ratio.
    """
    fit = inverse_line_fit(crossings)
    if fit is None:
        return None
    slope, intercept = fit
    inverse = intercept + slope * ratio
    if inverse <= 0:
        return None
    return round(1.0 / inverse)


def even_gate_count(value: float) -> int:
    return max(2, 2 * round(value / 2))


def resolved_start_ratio(settings: LevelCrossingSettings) -> float:
    low, high = settings.ratio_bounds
    if settings.start_ratio is None:
        return high if settings.trace_direction == "down" else low
    if not low <= settings.start_ratio <= high:
        raise ValueError(
            f"start_ratio={settings.start_ratio} outside ratio_bounds={settings.ratio_bounds}."
        )
    return settings.start_ratio


def ratio_sweep(settings: LevelCrossingSettings) -> list[float]:
    """Inclusive ratio sweep, descending by default from the high-ratio end."""
    if settings.ratio_step <= 0:
        raise ValueError(f"ratio_step must be positive, got {settings.ratio_step}.")
    low, high = settings.ratio_bounds
    start = resolved_start_ratio(settings)
    stop = low if settings.trace_direction == "down" else high
    direction = -1.0 if settings.trace_direction == "down" else 1.0

    ratios: list[float] = []
    current = start
    while direction * (current - stop) <= 1e-9:
        ratios.append(float(current))
        current += direction * settings.ratio_step
    if ratios and abs(ratios[-1] - stop) > 1e-9:
        ratios.append(float(stop))
    return ratios


def initial_search_cap(settings: LevelCrossingSettings) -> int:
    """Low initial high endpoint; deliberately below the global max bound."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    cap = even_gate_count(settings.initial_bracket_fraction * n_gates_max)
    return min(n_gates_max - 2, max(n_gates_min + 2, cap))


def next_bracket(
    crossings: list[RMBConfig],
    settings: LevelCrossingSettings,
    ratio: float,
) -> tuple[int, int]:
    """
    Search window for the next ratio, anchored by the traced contour so far.

    Moving downward in two-qubit ratio should require more total gates, so the
    next high endpoint grows from the nearest previous crossing instead of
    jumping to ``n_gates_bounds[1]``.
    """
    n_gates_min, n_gates_max = settings.n_gates_bounds
    cap = n_gates_max - 2
    if not crossings:
        return n_gates_min, initial_search_cap(settings)

    predicted = predict_crossing(crossings, ratio)
    if predicted is not None:
        lo = even_gate_count(settings.trace_gate_shrink * predicted)
        hi = even_gate_count(settings.trace_gate_growth * predicted)
    else:
        anchor = crossings[-1]
        if ratio < anchor.ratio_2_qb_gates:
            lo = even_gate_count(settings.trace_gate_shrink * anchor.n_gates)
            hi = even_gate_count(settings.trace_gate_growth * anchor.n_gates)
        else:
            lo = n_gates_min
            hi = even_gate_count(anchor.n_gates)

    lo = max(n_gates_min, min(cap - 2, lo))
    hi = min(cap, max(lo + 2, hi))
    return lo, hi


def find_crossing(*, backend, rng: RNGGenerator, data: RMBData, budget: Budget,
                  settings: LevelCrossingSettings,
                  ratio: float, lo: int, hi: int) -> RMBConfig | None:
    """
    Search low-to-high, then bisect in total gates at fixed ratio.

    Fidelity decreases monotonically with gate count, so the bracket needs
    fidelity confidently above 0.5 at ``lo`` and confidently below at a grown
    high endpoint. The global maximum gate bound is not probed directly.
    The bisection only branches on confident side decisions; any point whose
    shot cap cannot tell it from 0.5, endpoints included, is statistically on
    the line and is returned as the crossing directly.

    Returns
    -------
    RMBConfig | None
        The crossing config, or ``None`` when the line lies confidently
        outside the bracket or the budget runs out.
    """
    config_lo = settings.make_config(lo, ratio)
    p_lo, above_lo = probe_fidelity(backend=backend, rng=rng, data=data,
                                    config=config_lo, budget=budget, settings=settings)
    if p_lo is None:
        return None
    if above_lo < settings.decision_confidence:
        if 1.0 - above_lo < settings.decision_confidence:
            return config_lo
        return None

    previous = lo
    probe_hi = even_gate_count(
        max(lo + 2, lo * settings.low_to_high_growth_factor)
    )
    while probe_hi < hi:
        config_hi = settings.make_config(probe_hi, ratio)
        p_hi, above_hi = probe_fidelity(backend=backend, rng=rng, data=data,
                                        config=config_hi, budget=budget,
                                        settings=settings)
        if p_hi is None:
            return None
        print_progress(
            settings,
            budget,
            f"  bracket ratio={ratio:.3f}: n_gates={probe_hi} p={p_hi:.2f}",
        )
        if above_hi < settings.decision_confidence:
            if 1.0 - above_hi < settings.decision_confidence:
                return config_hi
            lo = previous
            hi = probe_hi
            break
        previous = probe_hi
        probe_hi = even_gate_count(
            max(probe_hi + 2, probe_hi * settings.low_to_high_growth_factor)
        )
    else:
        lo = previous

    config_hi = settings.make_config(hi, ratio)
    p_hi, above_hi = probe_fidelity(backend=backend, rng=rng, data=data,
                                    config=config_hi, budget=budget, settings=settings)
    if p_hi is None:
        return None
    if 1.0 - above_hi < settings.decision_confidence:
        if above_hi < settings.decision_confidence:
            return config_hi
        return None

    while hi - lo > max(2, settings.n_gates_resolution * hi):
        mid = 2 * round((lo + hi) / 4)
        if mid in (lo, hi):
            break
        p_mid, above_mid = probe_fidelity(backend=backend, rng=rng, data=data,
                                          config=settings.make_config(mid, ratio),
                                          budget=budget, settings=settings)
        if p_mid is None:
            break
        print_progress(settings, budget, f"  bisect ratio={ratio:.3f}: n_gates={mid} p={p_mid:.2f}")
        if above_mid >= settings.decision_confidence:
            lo = mid
        elif 1.0 - above_mid >= settings.decision_confidence:
            hi = mid
        else:
            return settings.make_config(mid, ratio)

    if hi - lo > max(2, settings.n_gates_resolution * hi):
        # The budget died before the bracket converged; the midpoint would be
        # a guess, not a measurement.
        return None
    return settings.make_config((lo + hi) // 2, ratio)


def run_with_budget(settings: LevelCrossingSettings) -> tuple[RMB, list[RMBConfig], Budget]:
    """
    Find the fidelity = 0.5 line and trace it through the ratio sweep.

    Returns
    -------
    tuple[RMB, list[RMBConfig], Budget]
        The RMB holding all recorded data, the crossing configs in traced
        ratio order, and the final budget state.
    """
    rng, rmb, budget = start_run(settings)
    data = rmb._data
    crossings: list[RMBConfig] = []

    lo, hi = next_bracket(crossings, settings, resolved_start_ratio(settings))
    for ratio in ratio_sweep(settings):
        if budget.remaining_hqc <= 0:
            break
        crossing = find_crossing(backend=rmb.backend, rng=rng, data=data, budget=budget,
                                 settings=settings, ratio=ratio, lo=lo, hi=hi)
        fallback_lo, fallback_hi = next_bracket([], settings, ratio)
        if crossing is None and (lo, hi) != (fallback_lo, fallback_hi):
            # The predicted bracket missed the line; retry from the low-gate
            # initial window, still avoiding the global max endpoint.
            lo, hi = fallback_lo, fallback_hi
            crossing = find_crossing(backend=rmb.backend, rng=rng, data=data, budget=budget,
                                     settings=settings, ratio=ratio, lo=lo, hi=hi)

        if crossing is not None:
            crossings.append(crossing)
            print_crossing(settings, budget, "Crossing", crossing)
        else:
            # A crossing in total gates exists at every ratio, so a miss means
            # either the budget ran dry or unlucky endpoint reads; the next
            # ratio is still worth trying with the remaining budget.
            reason = "budget exhausted" if budget.remaining_hqc <= 0 else "endpoints unresolved"
            print_progress(settings, budget, f"No crossing found at ratio={ratio:.3f} ({reason})")

        next_ratio = (
            ratio - settings.ratio_step
            if settings.trace_direction == "down"
            else ratio + settings.ratio_step
        )
        lo, hi = next_bracket(crossings, settings, next_ratio)

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


def run(settings: LevelCrossingSettings) -> tuple[RMB, list[RMBConfig]]:
    """Trace the fidelity = 0.5 line and return the run data and crossings."""
    rmb, crossings, _ = run_with_budget(settings)
    return rmb, crossings


if __name__ == "__main__":
    rmb, crossings = run(LevelCrossingSettings())
