"""
High-confidence ground-truth baseline for the fidelity = 0.5 line in (total
gates, two-qubit gate ratio) space at fixed n_qubits.

This method ignores HQC cost and trades it for confidence. It fits the
physical randomized-benchmarking decay surface
``p(N1, N2) = B + A exp(-(gamma_1 N1 + gamma_2 N2))`` (see
:class:`~.common.PhysicalDecaySurface`) and refines it iteratively:

1. seed - bisect the full gate range at a few ratios to find a first line and
   fit the decay surface;
2. cover - sample tight brackets straddling the predicted crossing at many
   ratios, relocating the line precisely where the fit is most sensitive;
3. anchor - sample a coarse low-to-high gate rake at a few ratios so the off-
   line decay shape pins the amplitude, floor, and per-gate rates;
4. refit and repeat until the four physical parameters converge.

The traced line is then read off the converged surface in closed form, so it
is smooth rather than quantized onto the bisection lattice. The bisection
primitive is reused from :mod:`level_crossing`; only the orchestration and the
parametric fit differ.
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
    MeasurementRequest,
    PhysicalDecaySurface,
    batch_hqc_cost,
    print_crossing,
    print_experiment_summary,
    print_progress,
    save_crossings,
    spend_request_batch,
    start_run,
    try_fit_physical_decay,
)
from sympleq.applications.randomized_benchmarking.experiments.level_crossing import (
    LevelCrossingSettings,
    find_crossing,
)


@dataclass(frozen=True)
class BaselineCrossingSettings(LevelCrossingSettings):
    """
    Settings for the iterative physical-decay baseline.

    Inherits the level-crossing bisection knobs (it reuses
    :func:`~.level_crossing.find_crossing` to locate crossings) but spends
    without a budget cap and at high confidence, so the achievable fidelity
    tolerance is set only by ``max_shots_per_config``. The scatter plots are
    drawn at full resolution because the run samples densely along the line.

    Parameters
    ----------
    seed_ratio_count : int
        Ratios at which the first full-range bisection seeds the fit.
    cover_ratio_count : int
        Ratios covered each iteration with a tight bracket on the line.
    cover_bracket_fraction : float
        Bracket half-width as a fraction of the predicted crossing gate count.
    anchor_ratio_count : int
        Ratios that get an off-line decay rake each iteration.
    anchor_gate_points : int
        Gate counts in each rake, spread across ``n_gates_bounds``.
    anchor_shots : int
        Shots measured on each rake point, fixed so the fidelity is pinned even
        far from 0.5 (the side-deciding probe would stop too early there).
    max_iterations : int
        Hard cap on cover/anchor/refit rounds.
    line_tol : float
        Relative move of the predicted crossing line, compared across the ratio
        grid, below which the fit is declared converged. Measuring the line
        (counts in the hundreds) rather than the raw parameters keeps the test
        immune to the relative jitter of a near-zero ``gamma_1`` or ``floor``.
    line_ratio_count : int
        Ratios at which the converged surface's 0.5 line is reported.
    fit_restarts : int
        Gamma-scale restarts of the physical-decay fit.
    batch_max_configs : int
        Cap on configs packed into one stitched anchor submission.
    """
    hqc_budget: float = float("inf")
    decision_confidence: float = 0.999
    max_shots_per_config: int = 200
    n_gates_resolution: float = 0.02
    scatter_merge_bins: tuple[int, int] | None = None
    save_path: str | Path | None = "baseline_crossing.json"

    seed_ratio_count: int = 6
    cover_ratio_count: int = 41
    cover_bracket_fraction: float = 0.3
    anchor_ratio_count: int = 4
    anchor_gate_points: int = 6
    anchor_shots: int = 48
    max_iterations: int = 5
    line_tol: float = 0.02
    line_ratio_count: int = 60
    fit_restarts: int = 6
    batch_max_configs: int = 60


def _ratio_grid(settings: BaselineCrossingSettings, count: int) -> np.ndarray:
    """``count`` evenly spaced ratios across the search bounds."""
    return np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], max(1, count))


def spend_fixed_shots(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    configs: list[RMBConfig],
    shots: int,
    budget: Budget,
    settings: BaselineCrossingSettings,
) -> int:
    """
    Top each config up to ``shots`` recorded outcomes in stitched batches.

    Unlike the side-deciding probe, this measures a fixed number of shots
    regardless of which side of 0.5 a config is on, so off-line anchor points
    get tight fidelities the decay fit can use.
    """
    requests: list[MeasurementRequest] = []
    seen: set[RMBConfig] = set()
    for config in configs:
        if config in seen:
            continue
        seen.add(config)
        current = data[config].num_runs() if config in data else 0
        needed = shots - current
        if needed > 0:
            requests.append(MeasurementRequest(config, needed))

    spent = 0
    for start in range(0, len(requests), max(1, settings.batch_max_configs)):
        batch = requests[start:start + max(1, settings.batch_max_configs)]
        spend_request_batch(backend, rng, data, batch, seed=settings.rng_seed)
        batch_shots = sum(request.shots for request in batch)
        budget.spend_batch(batch_hqc_cost(batch), batch_shots)
        spent += batch_shots
    return spent


def cover_line(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    surface: PhysicalDecaySurface,
    settings: BaselineCrossingSettings,
    budget: Budget,
) -> int:
    """Relocate the crossing at many ratios inside a tight predicted bracket."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    located = 0
    for ratio in _ratio_grid(settings, settings.cover_ratio_count):
        predicted = surface.crossing_n_gates(float(ratio))
        if predicted is None:
            continue
        half_width = settings.cover_bracket_fraction * predicted
        lo = max(float(n_gates_min), predicted - half_width)
        hi = min(float(n_gates_max), predicted + half_width)
        if hi - lo < 2.0:
            continue
        crossing = find_crossing(backend=backend, rng=rng, data=data, budget=budget,
                                 settings=settings, ratio=float(ratio), lo=lo, hi=hi)
        if crossing is not None:
            located += 1
    return located


def anchor_decay(
    *,
    backend,
    rng: RNGGenerator,
    data: RMBData,
    settings: BaselineCrossingSettings,
    budget: Budget,
) -> int:
    """Measure off-line gate rakes that pin the decay amplitude, floor, and rates."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    rake = np.linspace(n_gates_min, n_gates_max, max(2, settings.anchor_gate_points))
    configs = [settings.make_config(float(gates), float(ratio))
               for ratio in _ratio_grid(settings, settings.anchor_ratio_count)
               for gates in rake]
    return spend_fixed_shots(backend=backend, rng=rng, data=data, configs=configs,
                             shots=settings.anchor_shots, budget=budget, settings=settings)


def fit(data: RMBData, settings: BaselineCrossingSettings) -> PhysicalDecaySurface | None:
    """Fit the physical decay surface with the run's restart count."""
    return try_fit_physical_decay(data, settings)


def predicted_line(surface: PhysicalDecaySurface,
                   settings: BaselineCrossingSettings) -> np.ndarray:
    """Predicted crossing gate count over the ratio grid; NaN where none lies in bounds."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    ratios = _ratio_grid(settings, settings.cover_ratio_count)
    line = np.full(len(ratios), np.nan)
    for index, ratio in enumerate(ratios):
        n_gates = surface.crossing_n_gates(float(ratio))
        if n_gates is not None and n_gates_min <= n_gates <= n_gates_max:
            line[index] = n_gates
    return line


def line_relative_change(previous: np.ndarray | None, current: np.ndarray) -> float:
    """
    Largest relative move of the predicted line between iterations.

    Compared only at ratios where both iterations place a crossing in bounds,
    with a one-gate floor in the denominator so the metric stays well-scaled
    (the crossing counts are in the hundreds) and is immune to the relative
    jitter of near-zero fit parameters. ``inf`` on the first fit.
    """
    if previous is None:
        return float("inf")
    mask = np.isfinite(previous) & np.isfinite(current)
    if not np.any(mask):
        return float("inf")
    return float(np.max(np.abs(current[mask] - previous[mask])
                        / (np.abs(previous[mask]) + 1.0)))


def params_line(surface: PhysicalDecaySurface) -> str:
    """One-line ``(gamma_1, gamma_2, A, B)`` summary for the convergence trace."""
    return (f"g1={surface.gamma_1:.3e}  g2={surface.gamma_2:.3e}  "
            f"A={surface.amplitude:.4f}  B={surface.floor:.4f}")


def line_from_surface(surface: PhysicalDecaySurface,
                      settings: BaselineCrossingSettings) -> list[RMBConfig]:
    """Configs on the converged surface's fidelity = 0.5 line, in ratio order."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    configs: list[RMBConfig] = []
    seen: set[RMBConfig] = set()
    for ratio in _ratio_grid(settings, settings.line_ratio_count):
        n_gates = surface.crossing_n_gates(float(ratio))
        if n_gates is None or not (n_gates_min <= n_gates <= n_gates_max):
            continue
        config = settings.make_config(float(n_gates), float(ratio))
        if config not in seen:
            seen.add(config)
            configs.append(config)
    return sorted(configs, key=lambda c: (c.ratio_2_qb_gates, c.n_gates))


def load_fit(path: str | Path | None = None,
             settings: BaselineCrossingSettings | None = None) -> PhysicalDecaySurface | None:
    """
    Refit and return the physical decay surface from a saved baseline run.

    The run saves the measured outcomes, not the fit itself, so the converged
    surface is recovered by refitting the loaded data (the fit is fast and
    deterministic). ``path`` defaults to the settings' ``save_path``; print the
    returned surface's :meth:`~..common.PhysicalDecaySurface.report` for the
    line equation, or read its coefficients with ``crossing_line()``.
    """
    settings = settings or BaselineCrossingSettings()
    if path is None:
        path = settings.save_path
    return try_fit_physical_decay(RMB.load(path)._data, settings)


def run(settings: BaselineCrossingSettings) -> tuple[RMB, list[RMBConfig]]:
    """
    Trace the fidelity = 0.5 line as a converged physical-decay fit.

    Returns
    -------
    tuple[RMB, list[RMBConfig]]
        The RMB holding all recorded data and the configs on the converged
        surface's 0.5 line, in increasing ratio order.
    """
    rng, rmb, budget = start_run(settings)
    data = rmb._data
    n_gates_min, n_gates_max = settings.n_gates_bounds

    # Seed: full-range bisection at a few ratios for a first line and fit.
    for ratio in _ratio_grid(settings, settings.seed_ratio_count):
        find_crossing(backend=rmb.backend, rng=rng, data=data, budget=budget,
                      settings=settings, ratio=float(ratio),
                      lo=n_gates_min, hi=n_gates_max)
    surface = fit(data, settings)
    if settings.verbose:
        print(f"\n{'fit':>11}   "
              f"{'g1':<11} {'g2':<11} {'A':<8} {'B':<8} {'Δline':<8}")
        if surface is not None:
            print(f"{'seed':>11}   {params_line(surface)}")

    best_surface = surface
    best_change = float("inf")
    previous_line = predicted_line(surface, settings) if surface is not None else None
    stop_reason = "fit failed on the seed data"
    for iteration in range(settings.max_iterations):
        if surface is None:
            break
        cover_line(backend=rmb.backend, rng=rng, data=data, surface=surface,
                   settings=settings, budget=budget)
        anchor_decay(backend=rmb.backend, rng=rng, data=data,
                     settings=settings, budget=budget)
        surface = fit(data, settings)
        if surface is None:
            stop_reason = "fit failed during refinement"
            break

        current_line = predicted_line(surface, settings)
        change = line_relative_change(previous_line, current_line)
        if settings.verbose:
            print(f"{f'iter {iteration + 1}':>11}   {params_line(surface)}  {change:.4f}")
        if change <= best_change:
            best_change = change
            best_surface = surface
        previous_line = current_line
        if change < settings.line_tol:
            stop_reason = f"line converged after {iteration + 1} iterations"
            break
    else:
        stop_reason = (f"reached the {settings.max_iterations}-iteration cap; "
                       "kept the most line-stable fit")

    surface = best_surface
    if settings.verbose and surface is not None:
        print(f"\n{surface.report()}")

    crossings = line_from_surface(surface, settings) if surface is not None else []
    if settings.verbose:
        print_progress(settings, budget, f"\nTraced {len(crossings)} line points")
        for config in crossings:
            print(f"  ratio={config.ratio_2_qb_gates:.3f} n_gates={config.n_gates:>6}")
    print_experiment_summary(data, settings, budget, stop_reason=stop_reason)

    base_path = None
    if settings.save_path is not None:
        base_path = save_crossings(rmb, settings, budget, crossings)

    print(try_fit_physical_decay(rmb._data, settings).report())

    if settings.plot:
        import matplotlib.pyplot as plt
        from sympleq.applications.randomized_benchmarking.experiments.plots import (
            plot_level_line,
        )
        bins = settings.scatter_merge_bins or (None, None)
        png_path = None if base_path is None else base_path.parent / f"{base_path.stem}.png"
        plot_level_line(data, crossings,
                        n_1qb_gates_bin=bins[0], n_2qb_gates_bin=bins[1],
                        png_path=png_path, show=False)
        plt.show()

    return rmb, crossings


if __name__ == "__main__":
    settings = BaselineCrossingSettings(n_qubits=20)
    rmb, crossings = run(settings)
    # print(load_fit(settings=settings).report())
