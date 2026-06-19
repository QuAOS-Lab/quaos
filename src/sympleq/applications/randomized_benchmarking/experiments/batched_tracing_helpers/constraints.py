"""Monotonicity constraints for batched monotone tracing."""
from __future__ import annotations

from typing import Any

from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.batched_tracing_helpers.hints import (
    analytic_prediction,
)


def monotonic_slack(settings: Any, anchor: RMBConfig) -> int:
    fraction = (
        settings.n_gates_resolution
        if settings.monotonic_slack_fraction is None
        else settings.monotonic_slack_fraction
    )
    return max(settings.monotonic_min_slack_gates, round(fraction * anchor.n_gates))


def trace_gate_jump_slack(
    settings: Any,
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


def downward_growth_floor(settings: Any, anchor: RMBConfig, ratio: float) -> int:
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
    settings: Any,
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
    settings: Any,
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
    settings: Any,
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

