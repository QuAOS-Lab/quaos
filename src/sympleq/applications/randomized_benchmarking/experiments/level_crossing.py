"""
Trace the fidelity = 0.5 line in (total gates, two-qubit gate ratio) space at
fixed n_qubits.

The experiment finds a first crossing by bisecting in the total gate count at
a small fixed two-qubit ratio, then tracks the line by stepping the ratio up
and re-bisecting inside a bracket predicted from the previous crossings.
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

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.RMB import RMB
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
)


@dataclass(frozen=True)
class LevelCrossingSettings(CrossingSettings):
    """
    Settings for a fidelity = 0.5 level-crossing experiment.

    The problem definition (qubits, bounds, budget, ...) lives in
    :class:`~.common.CrossingSettings`; crossings are traced from
    ``ratio_bounds[0]`` up to ``ratio_bounds[1]``.

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
    ratio_step: float = 0.05
    n_gates_resolution: float = 0.1
    max_shots_per_config: int = 8
    decision_confidence: float = 0.75
    save_path: str | Path | None = "level_crossing.json"


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
        above = posterior_above(estimator)
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
            stitch_batch_size(bare_hqc, settings.max_cost_per_run, settings.max_shots_per_config),
            settings.max_shots_per_config - estimator.num_runs())
        while n_circuits > 0 and not budget.can_afford(stitched_batch_hqc(n_circuits * bare_hqc)):
            n_circuits -= 1
        if n_circuits <= 0:
            break
        # Only configs with recorded outcomes enter the data set; implied or
        # unaffordable probes must not leave empty estimators behind.
        data[config] = estimator
        spend_measurements(backend, rng, data, config, n_circuits, seed=settings.rng_seed)
        budget.spend_batch(stitched_batch_hqc(n_circuits * bare_hqc), n_circuits)

    above = posterior_above(estimator)
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


def find_crossing(*, backend, rng: RNGGenerator, data: RMBData, budget: Budget,
                  settings: LevelCrossingSettings,
                  ratio: float, lo: int, hi: int) -> RMBConfig | None:
    """
    Bisect in total gates for the fidelity = 0.5 crossing at fixed ratio.

    Fidelity decreases monotonically with gate count, so the bracket needs
    fidelity confidently above 0.5 at ``lo`` and confidently below at ``hi``.
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


def run(settings: LevelCrossingSettings) -> tuple[RMB, list[RMBConfig]]:
    """
    Find the fidelity = 0.5 line and trace it up in the two-qubit gate ratio.

    Returns
    -------
    tuple[RMB, list[RMBConfig]]
        The RMB holding all recorded data and the crossing configs in
        increasing ratio order.
    """
    rng, rmb, budget = start_run(settings)
    data = rmb._data
    crossings: list[RMBConfig] = []

    n_gates_min, n_gates_max = settings.n_gates_bounds
    lo, hi = n_gates_min, n_gates_max
    step_index = 0
    ratio = settings.ratio_bounds[0]

    while ratio <= settings.ratio_bounds[1] + 1e-9 and budget.remaining_hqc > 0:
        crossing = find_crossing(backend=rmb.backend, rng=rng, data=data, budget=budget,
                                 settings=settings, ratio=ratio, lo=lo, hi=hi)
        if crossing is None and (lo, hi) != (n_gates_min, n_gates_max):
            # The predicted bracket missed the line; retry with full bounds.
            lo, hi = n_gates_min, n_gates_max
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

        step_index += 1
        ratio = settings.ratio_bounds[0] + step_index * settings.ratio_step

        # Bracket the next crossing around the inverse-linear fit prediction;
        # before the fit is possible, search below the previous crossing.
        predicted = predict_crossing(crossings, ratio)
        if predicted is not None:
            lo = max(n_gates_min, predicted // 2)
            hi = min(n_gates_max, 2 * predicted)
        elif crossings:
            lo, hi = n_gates_min, crossings[-1].n_gates
        else:
            lo, hi = n_gates_min, n_gates_max
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

    return rmb, crossings


if __name__ == "__main__":
    rmb, crossings = run(LevelCrossingSettings())
