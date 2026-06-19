"""
Deterministic batched-window monotone tracing of the fidelity = 0.5 boundary.

This variant uses the local SympleQ backend and probes small stitched batches
inside predicted contour windows. A crossing is accepted only when a batch
contains a confident above-to-below bracket, and the bracket midpoint is used
as the estimate.
"""
from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Literal

_SRC_ROOT = Path(__file__).resolve().parents[4]
sys.path = [path for path in sys.path if path != str(_SRC_ROOT)]
sys.path.insert(0, str(_SRC_ROOT))

import numpy as np
from numpy.random import Generator as RNGGenerator
from scipy.special import betainc

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementRequest,
    RMBBackend,
)
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    batch_hqc_cost,
    Budget,
    CrossingSettings,
    print_crossing,
    print_experiment_summary,
    print_progress,
    save_crossings,
    measurement_rng,
    start_run,
    try_fit_monotone_fidelity_surface,
)
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.core.noise.noise_model import GenericNoise


BASE_1Q_PAULI_ERROR = 0.000025
BASE_2Q_PAULI_ERROR = 0.00079


def new_estimator() -> BayesianEstimator:
    """Estimator constructor local to this experiment."""
    return BayesianEstimator()


def posterior_above(estimator: BayesianEstimator) -> float:
    """Posterior probability that the Boolean success probability is above 0.5."""
    alpha, beta = estimator.posterior_alpha_beta()
    return float(1.0 - betainc(alpha, beta, 0.5))


def spend_request_batch_local(
    backend: RMBBackend,
    rng: RNGGenerator,
    data: RMBData,
    requests: list[MeasurementRequest],
    budget: Budget,
    ratio_budget: Budget,
    *,
    seed: int | None,
) -> None:
    """Record a stitched batch using this module's estimator constructor."""
    requests = [request for request in requests if request.shots > 0]
    if not requests:
        return
    cost = batch_hqc_cost(requests)
    if cost > spendable_ratio_budget(budget, ratio_budget):
        return

    offsets: dict[RMBConfig, int] = {}
    for request in requests:
        estimator = data.setdefault(request.config, new_estimator())
        offsets.setdefault(request.config, estimator.num_runs())

    shot_rng = None
    if seed is not None:
        def shot_rng(config: RMBConfig, index: int) -> RNGGenerator:
            return measurement_rng(seed, config, offsets[config] + index)

    outcomes = backend.fidelity_estimation(requests, rng, shot_rng=shot_rng).outcomes
    for config, results in outcomes.items():
        for outcome in results:
            data[config].record(bool(outcome))
    n_circuits = sum(request.shots for request in requests)
    budget.spend_batch(cost, n_circuits)
    ratio_budget.spend_batch(cost, n_circuits)


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
class BatchedMonotoneTracingSettings(CrossingSettings):
    """
    Settings for deterministic local contour tracing.

    The contour is traced in fixed-ratio slices. Each slice submits a small
    stitched batch inside a predicted gate-count window. If the batch contains
    a confident above-to-below transition, the midpoint of that local bracket
    is used as the crossing estimate.
    """

    ratio_step: float = 0.1
    start_ratio: float | None = None
    trace_direction: Literal["up", "down", "both"] = "down"
    save_path: str | Path | None = "batched_monotone_tracing.json"
    backend_factory: Callable[[CrossingSettings, RNGGenerator], RMBBackend] = (
        sympleq_backend_factory
    )

    adaptive_ratio_step: bool = True
    min_ratio_step: float = 0.1
    max_ratio_step: float = 0.2
    target_hqc_per_ratio: float = 35.0

    n_gates_resolution: float = 0.1
    decision_confidence: float = 0.70
    max_hqc_per_ratio: float | None = 35.0

    initial_bracket_fraction: float = 0.35
    low_to_high_growth_factor: float = 1.7
    bracket_half_width_fraction: float = 0.24
    trace_gate_growth: float = 2.4
    trace_gate_shrink: float = 0.65
    bracket_hint: Literal["line", "surface", "both"] = "surface"
    use_surface_bracket_hint: bool = True
    surface_min_configs: int = 12
    initial_batch_points: int = 5
    initial_batch_shots: int = 3
    initial_prediction_window_fraction: float = 0.35
    initial_prediction_min_width: int = 80
    trace_batch_points: int = 5
    trace_batch_shots: int = 3
    trace_batch_window_fraction: float = 0.14
    trace_batch_min_width: int = 40
    downward_prediction_quantile: float = 0.75
    upward_prediction_quantile: float = 0.25

    enforce_ratio_monotonicity: bool = True
    monotonic_min_slack_gates: int = 2
    monotonic_slack_fraction: float | None = 0.65
    max_trace_gate_growth_fraction: float | None = 0.9
    max_trace_gate_shrink_fraction: float | None = 0.7
    trace_gate_jump_min_slack: int = 120
    trace_gate_jump_step_scale_cap: float = 1.5
    trace_gate_jump_min_anchors: int = 2
    downward_growth_floor_fraction: float = 0.45
    downward_growth_floor_min_gates: int = 24


def implied_above(
    data: RMBData,
    config: RMBConfig,
    settings: BatchedMonotoneTracingSettings,
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


def confidently_above(above: float, settings: BatchedMonotoneTracingSettings) -> bool:
    return above >= settings.decision_confidence


def confidently_below(above: float, settings: BatchedMonotoneTracingSettings) -> bool:
    return 1.0 - above >= settings.decision_confidence


def resolved_start_ratio(settings: BatchedMonotoneTracingSettings) -> float:
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
    settings: BatchedMonotoneTracingSettings,
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


def ratio_sweeps(settings: BatchedMonotoneTracingSettings) -> list[list[float]]:
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
    settings: BatchedMonotoneTracingSettings,
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
    settings: BatchedMonotoneTracingSettings,
    budget: Budget,
    crossings: list[RMBConfig],
) -> Budget:
    """Cap later ratios without restricting the first anchor search."""
    if crossings and settings.max_hqc_per_ratio is not None:
        return Budget(remaining_hqc=min(settings.max_hqc_per_ratio, budget.remaining_hqc))
    return Budget(remaining_hqc=budget.remaining_hqc)


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
    settings: BatchedMonotoneTracingSettings,
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
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
) -> tuple[int, int]:
    """Bracket the next fixed-ratio search from line and monotone hints."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    last = crossings[-1] if crossings else None
    predictions = []
    if settings.bracket_hint in ("line", "both"):
        prediction = line_fit_prediction(crossings, ratio)
        if prediction is not None:
            predictions.append(prediction)
    if settings.bracket_hint in ("surface", "both"):
        prediction = surface_prediction(data, settings, ratio)
        if prediction is not None:
            predictions.append(prediction)
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
    settings: BatchedMonotoneTracingSettings,
    sweep: list[float],
) -> tuple[int, int]:
    """Initial bracket for a sweep, biased low until the trace has an anchor."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    if not sweep:
        return n_gates_min, n_gates_max
    hi = round(settings.initial_bracket_fraction * n_gates_max)
    return n_gates_min, min(n_gates_max, max(n_gates_min + 2, hi))


def analytic_prediction(
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
) -> int | None:
    """Predict a crossing from the analytic Lindblad contour."""
    from sympleq.applications.randomized_benchmarking.experiments.plots import (
        analytic_gate_counts,
        analytic_noise_scales,
    )

    one_q_noise_scale, two_q_noise_scale = analytic_noise_scales(settings)
    prediction = analytic_gate_counts(
        np.array([ratio]),
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )[0]
    if not np.isfinite(prediction) or prediction <= 0:
        return None
    return round(float(prediction))


def initial_prediction_bracket(
    data: RMBData,
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
) -> tuple[int, int, str]:
    """
    Choose the initial-anchor batch window from the best current estimate.

    Before any accepted crossing exists, the surface fit can still become
    useful after failed high-ratio probes because the accumulated configs carry
    shape information. If it is not ready, use the analytic contour as a cheap
    prior. Only fall back to the broad low-to-high sweep when neither estimate
    is usable.
    """
    n_gates_min, n_gates_max = settings.n_gates_bounds
    prediction = surface_prediction(data, settings, ratio)
    source = "surface"
    if prediction is None:
        prediction = analytic_prediction(settings, ratio)
        source = "analytic"
    if prediction is None:
        lo, hi = initial_sweep_bracket(settings, [ratio])
        return lo, hi, "broad"

    center = int(np.clip(prediction, n_gates_min, n_gates_max))
    half_width = max(
        2,
        round(settings.initial_prediction_window_fraction * center),
        settings.initial_prediction_min_width // 2,
    )
    lo = max(n_gates_min, 2 * round((center - half_width) / 2))
    hi = min(n_gates_max, 2 * round((center + half_width) / 2))
    if hi <= lo:
        lo, hi = initial_sweep_bracket(settings, [ratio])
        return lo, hi, "broad"
    return lo, hi, source


def centered_prediction_window(
    *,
    center: int,
    lo: int,
    hi: int,
    fraction: float,
    min_width: int,
) -> tuple[int, int] | None:
    """Return an even-gate window centered on a prediction and clamped to bounds."""
    if hi <= lo:
        return None
    center = int(np.clip(center, lo, hi))
    half_width = max(2, round(fraction * max(center, 1)), min_width // 2)
    window_lo = max(lo, 2 * round((center - half_width) / 2))
    window_hi = min(hi, 2 * round((center + half_width) / 2))
    if window_hi <= window_lo:
        return None
    return window_lo, window_hi


def trace_batch_window(
    data: RMBData,
    crossings: list[RMBConfig],
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
    lo: int,
    hi: int,
) -> tuple[int, int, str]:
    """Focused batch window inside the wider safe search interval."""
    predictions: list[tuple[str, int]] = []
    surface = surface_prediction(data, settings, ratio)
    if surface is not None:
        predictions.append(("surface", surface))
    line = line_fit_prediction(crossings, ratio)
    if line is not None:
        predictions.append(("line", line))
    analytic = analytic_prediction(settings, ratio)
    if analytic is not None:
        predictions.append(("analytic", analytic))
    if not predictions:
        return lo, hi, "trace"

    values = np.array([prediction for _, prediction in predictions], dtype=float)
    last = crossings[-1] if crossings else None
    if last is not None and ratio < last.ratio_2_qb_gates:
        quantile = settings.downward_prediction_quantile
    elif last is not None and ratio > last.ratio_2_qb_gates:
        quantile = settings.upward_prediction_quantile
    else:
        quantile = 0.5
    prediction = round(float(np.quantile(values, quantile)))
    source = "+".join(source for source, _ in predictions)

    window = centered_prediction_window(
        center=prediction,
        lo=lo,
        hi=hi,
        fraction=settings.trace_batch_window_fraction,
        min_width=settings.trace_batch_min_width,
    )
    if window is None:
        return lo, hi, "trace"
    batch_lo, batch_hi = window
    return batch_lo, batch_hi, source


def initial_batch_gate_counts(
    settings: BatchedMonotoneTracingSettings,
    lo: int,
    hi: int,
    *,
    geometric: bool,
    n_points: int | None = None,
) -> list[int]:
    """Gate-count stack for initial anchor discovery."""
    n_points = max(1, settings.initial_batch_points if n_points is None else n_points)
    if not geometric:
        gates = [2 * round(value / 2) for value in np.linspace(lo, hi, n_points)]
        return sorted({int(np.clip(gate, lo, hi)) for gate in gates})

    gates = [lo]
    current = lo
    while len(gates) < n_points and current < hi:
        current = max(
            current + 2,
            2 * round((current * settings.low_to_high_growth_factor) / 2),
        )
        gates.append(min(current, hi))
    return sorted(set(gates))


def submit_initial_batch(
    *,
    backend: RMBBackend,
    rng: RNGGenerator,
    data: RMBData,
    budget: Budget,
    ratio_budget: Budget,
    settings: BatchedMonotoneTracingSettings,
    configs: list[RMBConfig],
    shots: int,
) -> None:
    """Probe candidate configs in stitched batches."""
    requests = [MeasurementRequest(config, shots)
                for config in configs
                if shots > 0]
    batch: list[MeasurementRequest] = []
    for request in requests:
        candidate = [*batch, request]
        if batch and batch_hqc_cost(candidate) > settings.max_cost_per_run:
            spend_request_batch_local(
                backend, rng, data, batch, budget, ratio_budget,
                seed=settings.rng_seed)
            batch = [request]
        else:
            batch = candidate

        if batch_hqc_cost(batch) > spendable_ratio_budget(budget, ratio_budget):
            return

    if batch:
        spend_request_batch_local(
            backend, rng, data, batch, budget, ratio_budget,
            seed=settings.rng_seed)


def batch_window_bracket(
    *,
    backend: RMBBackend,
    rng: RNGGenerator,
    data: RMBData,
    budget: Budget,
    ratio_budget: Budget,
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
    lo: int,
    hi: int,
    source: str,
    n_points: int,
    shots: int,
) -> tuple[int, int] | None:
    """
    Batch-probe a vertical slice and return a confident bracket.

    This amortizes the stitched submission base cost across several configs and
    requires a confident above-to-below transition within the probed window.
    """
    gates = initial_batch_gate_counts(
        settings,
        lo,
        hi,
        geometric=source == "broad",
        n_points=n_points,
    )
    configs = [settings.make_config(gate, ratio) for gate in gates]
    print_progress(
        settings,
        budget,
        f"  batch {source} window ratio={ratio:.3f}: "
        f"n_gates {lo} to {hi}, probes={gates}",
    )
    submit_initial_batch(
        backend=backend,
        rng=rng,
        data=data,
        budget=budget,
        ratio_budget=ratio_budget,
        settings=settings,
        configs=configs,
        shots=shots,
    )

    previous_above: RMBConfig | None = None
    summaries = []
    for config in configs:
        estimator = data.get(config)
        if estimator is None or estimator.num_runs() <= 0:
            continue
        mean = estimator.posterior_mean()
        above = posterior_above(estimator)
        summaries.append(f"n={config.n_gates} p={mean:.2f}")
        if confidently_above(above, settings):
            previous_above = config
            continue
        if previous_above is not None and confidently_below(above, settings):
            print_progress(
                settings,
                budget,
                f"  batch {source} bracket: "
                + ", ".join(summaries)
                + f" -> [{previous_above.n_gates}, {config.n_gates}]",
            )
            return previous_above.n_gates, config.n_gates

    print_progress(
        settings,
        budget,
        f"  batch {source} no bracket: " + ", ".join(summaries),
    )
    return None


def initial_batch_bracket(
    *,
    backend: RMBBackend,
    rng: RNGGenerator,
    data: RMBData,
    budget: Budget,
    ratio_budget: Budget,
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
    lo: int,
    hi: int,
    source: str,
) -> tuple[int, int] | None:
    """Batch-probe the pre-anchor slice and require a visible bracket."""
    return batch_window_bracket(
        backend=backend,
        rng=rng,
        data=data,
        budget=budget,
        ratio_budget=ratio_budget,
        settings=settings,
        ratio=ratio,
        lo=lo,
        hi=hi,
        source=source,
        n_points=settings.initial_batch_points,
        shots=settings.initial_batch_shots,
    )


def midpoint_crossing_from_bracket(
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
    lo: int,
    hi: int,
) -> RMBConfig:
    """Use the midpoint of a batched bracket as the crossing estimate."""
    n_gates = 2 * round(((lo + hi) / 2) / 2)
    n_gates = int(np.clip(n_gates, lo, hi))
    return settings.make_config(n_gates, ratio)


def monotonic_slack(settings: BatchedMonotoneTracingSettings, anchor: RMBConfig) -> int:
    fraction = (
        settings.n_gates_resolution
        if settings.monotonic_slack_fraction is None
        else settings.monotonic_slack_fraction
    )
    return max(settings.monotonic_min_slack_gates, round(fraction * anchor.n_gates))


def trace_gate_jump_slack(
    settings: BatchedMonotoneTracingSettings,
    anchor: RMBConfig,
    ratio: float,
    *,
    fraction: float | None,
) -> int | None:
    """Allowed local gate-count movement from a neighbouring accepted crossing."""
    if fraction is None:
        return None
    step_scale = min(
        settings.trace_gate_jump_step_scale_cap,
        max(
            1.0,
            abs(ratio - anchor.ratio_2_qb_gates) / max(settings.ratio_step, 1e-12),
        ),
    )
    return max(
        settings.trace_gate_jump_min_slack,
        round(fraction * step_scale * anchor.n_gates),
    )


def downward_growth_floor(
    settings: BatchedMonotoneTracingSettings,
    anchor: RMBConfig,
    ratio: float,
) -> int:
    """Minimum expected gate-count increase when tracing to a lower ratio."""
    if ratio >= anchor.ratio_2_qb_gates:
        return 0
    current_prediction = analytic_prediction(settings, ratio)
    anchor_prediction = analytic_prediction(settings, anchor.ratio_2_qb_gates)
    if current_prediction is None or anchor_prediction is None:
        return settings.downward_growth_floor_min_gates
    expected_increase = current_prediction - anchor_prediction
    if expected_increase <= 0:
        return settings.downward_growth_floor_min_gates
    return max(
        settings.downward_growth_floor_min_gates,
        round(settings.downward_growth_floor_fraction * expected_increase),
    )


def monotonic_gate_bounds(
    settings: BatchedMonotoneTracingSettings,
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

    use_jump_guard = len(crossings) >= settings.trace_gate_jump_min_anchors

    if nearest_higher_ratio is not None:
        slack = monotonic_slack(settings, nearest_higher_ratio)
        lower = max(
            lower,
            nearest_higher_ratio.n_gates - slack - settings.monotonic_min_slack_gates,
        )
        growth_floor = downward_growth_floor(settings, nearest_higher_ratio, ratio)
        if growth_floor > 0:
            lower = max(lower, nearest_higher_ratio.n_gates + growth_floor)
        if use_jump_guard:
            growth_slack = trace_gate_jump_slack(
                settings,
                nearest_higher_ratio,
                ratio,
                fraction=settings.max_trace_gate_growth_fraction,
            )
            if growth_slack is not None:
                upper = min(upper, nearest_higher_ratio.n_gates + growth_slack)
    if nearest_lower_ratio is not None:
        slack = monotonic_slack(settings, nearest_lower_ratio)
        upper = min(
            upper,
            nearest_lower_ratio.n_gates + slack + settings.monotonic_min_slack_gates,
        )
        if use_jump_guard:
            shrink_slack = trace_gate_jump_slack(
                settings,
                nearest_lower_ratio,
                ratio,
                fraction=settings.max_trace_gate_shrink_fraction,
            )
            if shrink_slack is not None:
                lower = max(lower, nearest_lower_ratio.n_gates - shrink_slack)

    if lower > upper:
        return None
    return lower, upper


def apply_monotonic_bounds(
    settings: BatchedMonotoneTracingSettings,
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
    settings: BatchedMonotoneTracingSettings,
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


def search_one_ratio(
    rmb: RMB,
    rng: RNGGenerator,
    data: RMBData,
    crossings: list[RMBConfig],
    budget: Budget,
    settings: BatchedMonotoneTracingSettings,
    ratio: float,
) -> RMBConfig | None:
    """Find one fixed-ratio crossing using only batched local brackets."""
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
        lo, hi, source = trace_batch_window(
            data,
            crossings,
            settings,
            ratio,
            lo,
            hi,
        )
        bracket = batch_window_bracket(
            backend=rmb.backend,
            rng=rng,
            data=data,
            budget=budget,
            ratio_budget=ratio_budget,
            settings=settings,
            ratio=ratio,
            lo=lo,
            hi=hi,
            source=source,
            n_points=settings.trace_batch_points,
            shots=settings.trace_batch_shots,
        )
    else:
        lo, hi, source = initial_prediction_bracket(data, settings, ratio)
        bracket = initial_batch_bracket(
            backend=rmb.backend,
            rng=rng,
            data=data,
            budget=budget,
            ratio_budget=ratio_budget,
            settings=settings,
            ratio=ratio,
            lo=lo,
            hi=hi,
            source=source,
        )

    if bracket is None:
        return None

    lo, hi = bracket
    crossing = midpoint_crossing_from_bracket(settings, ratio, lo, hi)
    message = monotonicity_violation_message(settings, crossings, crossing)
    if message is not None:
        print_progress(settings, budget, message)
        return None
    return crossing


def run_with_budget(
    settings: BatchedMonotoneTracingSettings,
) -> tuple[RMB, list[RMBConfig], Budget]:
    """Run the deterministic local boundary tracer."""
    rng, rmb, budget = start_run(settings)
    if not isinstance(rmb.backend, SympleqBackend):
        raise TypeError(
            "BatchedMonotoneTracingSettings must use SympleqBackend; refusing non-local backend."
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

    if settings.verbose:
        print_run_score(rmb, crossings, budget, settings)

    if settings.plot:
        from sympleq.applications.randomized_benchmarking.experiments.plots import (
            plot_crossing_results,
        )

        plot_crossing_results(data, settings, crossings, base_path=base_path)

    return rmb, crossings, budget


def run(settings: BatchedMonotoneTracingSettings) -> tuple[RMB, list[RMBConfig]]:
    """Run the deterministic local boundary tracer and return data plus crossings."""
    rmb, crossings, _ = run_with_budget(settings)
    return rmb, crossings


def fitted_gate_counts(data: RMBData, ratios: np.ndarray) -> np.ndarray | None:
    """Parametric inverse-boundary gate counts on ``ratios``."""
    from sympleq.applications.randomized_benchmarking.experiments.plots import (
        parametric_boundary_fit,
    )

    fit = parametric_boundary_fit(data)
    if fit is None:
        return None
    q, slope, _ = fit
    gates = 1.0 / (q + slope * ratios)
    gates[~np.isfinite(gates)] = np.nan
    gates[gates <= 0.0] = np.nan
    return gates


def fit_score(data: RMBData, settings: BatchedMonotoneTracingSettings) -> float:
    """Score the fitted boundary against the analytic Lindblad boundary."""
    from sympleq.applications.randomized_benchmarking.experiments.plots import (
        analytic_gate_counts,
        analytic_noise_scales,
    )

    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 300)
    fitted = fitted_gate_counts(data, ratios)
    if fitted is None:
        return 0.0
    one_q_noise_scale, two_q_noise_scale = analytic_noise_scales(settings)
    analytic = analytic_gate_counts(
        ratios,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    mask = (
        np.isfinite(fitted)
        & np.isfinite(analytic)
        & (fitted > 0.0)
        & (analytic > 0.0)
    )
    if not np.any(mask):
        return 0.0
    log_error = np.log(fitted[mask] / analytic[mask])
    return float(np.exp(-np.sqrt(np.mean(log_error**2))))


def nan_quantiles(
    samples: np.ndarray,
    quantiles: tuple[float, ...],
) -> tuple[np.ndarray, ...]:
    """Column-wise nan-safe quantiles without all-NaN warnings."""
    output = [np.full(samples.shape[1], np.nan, dtype=float) for _ in quantiles]
    for column_index in range(samples.shape[1]):
        column = samples[:, column_index]
        column = column[np.isfinite(column)]
        if len(column) == 0:
            continue
        for output_array, quantile in zip(output, quantiles):
            output_array[column_index] = float(np.quantile(column, quantile))
    return tuple(output)


def bootstrap_coverage(
    data: RMBData,
    settings: BatchedMonotoneTracingSettings,
    *,
    n_bootstrap: int = 100,
) -> tuple[float, float]:
    """Fraction of analytic Lindblad points inside bootstrap 50% and 90% bands."""
    from sympleq.applications.randomized_benchmarking.experiments.plots import (
        analytic_gate_counts,
        analytic_noise_scales,
        parametric_boundary_bootstrap,
        parametric_boundary_gate_samples,
    )

    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 200)
    fits = parametric_boundary_bootstrap(
        data,
        n_bootstrap=n_bootstrap,
        seed=None if settings.rng_seed is None else settings.rng_seed + 500_000,
    )
    samples = parametric_boundary_gate_samples(fits, ratios)
    if len(samples) == 0:
        return 0.0, 0.0

    one_q_noise_scale, two_q_noise_scale = analytic_noise_scales(settings)
    analytic = analytic_gate_counts(
        ratios,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    q05, q25, q75, q95 = nan_quantiles(samples, (0.05, 0.25, 0.75, 0.95))
    valid50 = np.isfinite(q25) & np.isfinite(q75) & np.isfinite(analytic)
    valid90 = np.isfinite(q05) & np.isfinite(q95) & np.isfinite(analytic)
    coverage50 = (
        np.mean((q25[valid50] <= analytic[valid50])
                & (analytic[valid50] <= q75[valid50]))
        if np.any(valid50)
        else 0.0
    )
    coverage90 = (
        np.mean((q05[valid90] <= analytic[valid90])
                & (analytic[valid90] <= q95[valid90]))
        if np.any(valid90)
        else 0.0
    )
    return float(coverage50), float(coverage90)


def print_run_score(
    rmb: RMB,
    crossings: list[RMBConfig],
    budget: Budget,
    settings: BatchedMonotoneTracingSettings,
) -> None:
    """Print the direct-run benchmark metrics for one monotone tracing run."""
    score = fit_score(rmb._data, settings)
    coverage50, coverage90 = bootstrap_coverage(rmb._data, settings)
    print("\nRun score")
    print(f"  score: {score:.3f}")
    print(f"  analytic line inside bootstrap 50% band: {100.0 * coverage50:.1f}%")
    print(f"  analytic line inside bootstrap 90% band: {100.0 * coverage90:.1f}%")
    print(f"  spent: {budget.spent_hqc:.1f} HQC")
    print(f"  crossings: {len(crossings)}")


if __name__ == "__main__":
    settings = BatchedMonotoneTracingSettings(rng_seed=None)
    rmb, crossings, budget = run_with_budget(settings)
