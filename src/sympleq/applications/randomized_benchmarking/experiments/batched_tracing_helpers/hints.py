"""Prediction and window hint helpers for batched monotone tracing."""
from __future__ import annotations

from typing import Any

import numpy as np

from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    try_fit_monotone_fidelity_surface,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    analytic_gate_counts,
    analytic_noise_scales,
)


def measured_config_count(data: RMBData) -> int:
    return sum(1 for estimator in data.values() if estimator.num_runs() > 0)


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
    settings: Any,
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


def analytic_prediction(settings: Any, ratio: float) -> int | None:
    """Predict a crossing from the analytic Lindblad contour."""
    one_q_noise_scale, two_q_noise_scale = analytic_noise_scales(settings)
    prediction = analytic_gate_counts(
        np.array([ratio]),
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )[0]
    if not np.isfinite(prediction) or prediction <= 0:
        return None
    return round(float(prediction))


def next_bracket(
    data: RMBData,
    crossings: list[RMBConfig],
    settings: Any,
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


def initial_sweep_bracket(settings: Any, sweep: list[float]) -> tuple[int, int]:
    """Initial bracket for a sweep, biased low until the trace has an anchor."""
    n_gates_min, n_gates_max = settings.n_gates_bounds
    if not sweep:
        return n_gates_min, n_gates_max
    hi = round(settings.initial_bracket_fraction * n_gates_max)
    return n_gates_min, min(n_gates_max, max(n_gates_min + 2, hi))


def initial_prediction_bracket(
    data: RMBData,
    settings: Any,
    ratio: float,
) -> tuple[int, int, str]:
    """
    Choose the initial-anchor batch window from the best current estimate.

    If the surface is not ready, use the analytic contour as a cheap prior.
    Only fall back to the broad low-to-high sweep when neither estimate works.
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
    settings: Any,
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

