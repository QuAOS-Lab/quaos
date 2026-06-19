"""Batch window probing for batched monotone tracing."""
from __future__ import annotations

from typing import Any

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementRequest,
    RMBBackend,
)
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.batched_tracing_helpers.constraints import (
    apply_monotonic_bounds,
    monotonic_gate_bounds,
    monotonicity_violation_message,
)
from sympleq.applications.randomized_benchmarking.experiments.batched_tracing_helpers.hints import (
    initial_prediction_bracket,
    next_bracket,
    trace_batch_window,
)
from sympleq.applications.randomized_benchmarking.experiments.common import (
    batch_hqc_cost,
    Budget,
    print_progress,
    spend_request_batch,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    posterior_above,
)


def spendable_ratio_budget(global_budget: Budget, ratio_budget: Budget) -> float:
    return max(0.0, min(global_budget.remaining_hqc, ratio_budget.remaining_hqc))


def confidently_above(above: float, settings: Any) -> bool:
    return above >= settings.decision_confidence


def confidently_below(above: float, settings: Any) -> bool:
    return 1.0 - above >= settings.decision_confidence


def spend_request_batch_with_budget(
    backend: RMBBackend,
    rng: RNGGenerator,
    data: RMBData,
    requests: list[MeasurementRequest],
    budget: Budget,
    ratio_budget: Budget,
    *,
    seed: int | None,
) -> None:
    """Record a stitched batch and charge both global and ratio budgets."""
    requests = [request for request in requests if request.shots > 0]
    if not requests:
        return
    cost = batch_hqc_cost(requests)
    if cost > spendable_ratio_budget(budget, ratio_budget):
        return

    spend_request_batch(backend, rng, data, requests, seed=seed)
    n_circuits = sum(request.shots for request in requests)
    budget.spend_batch(cost, n_circuits)
    ratio_budget.spend_batch(cost, n_circuits)


def initial_batch_gate_counts(
    settings: Any,
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
    settings: Any,
    configs: list[RMBConfig],
    shots: int,
) -> None:
    """Probe candidate configs in stitched batches."""
    requests = [MeasurementRequest(config, shots) for config in configs if shots > 0]
    batch: list[MeasurementRequest] = []
    for request in requests:
        candidate = [*batch, request]
        if batch and batch_hqc_cost(candidate) > settings.max_cost_per_run:
            spend_request_batch_with_budget(
                backend,
                rng,
                data,
                batch,
                budget,
                ratio_budget,
                seed=settings.rng_seed,
            )
            batch = [request]
        else:
            batch = candidate

        if batch_hqc_cost(batch) > spendable_ratio_budget(budget, ratio_budget):
            return

    if batch:
        spend_request_batch_with_budget(
            backend,
            rng,
            data,
            batch,
            budget,
            ratio_budget,
            seed=settings.rng_seed,
        )


def batch_window_bracket(
    *,
    backend: RMBBackend,
    rng: RNGGenerator,
    data: RMBData,
    budget: Budget,
    ratio_budget: Budget,
    settings: Any,
    ratio: float,
    lo: int,
    hi: int,
    source: str,
    n_points: int,
    shots: int,
    probe_gates: list[int] | None = None,
) -> tuple[int, int] | None:
    """
    Batch-probe a vertical slice and return a confident bracket.

    This amortizes the stitched submission base cost across several configs and
    requires a confident above-to-below transition within the probed window.
    """
    gates = (
        sorted({
            int(np.clip(2 * round(gate / 2), lo, hi))
            for gate in probe_gates
        })
        if probe_gates is not None
        else initial_batch_gate_counts(
            settings,
            lo,
            hi,
            geometric=source == "broad",
            n_points=n_points,
        )
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
    settings: Any,
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
    settings: Any,
    ratio: float,
    lo: int,
    hi: int,
) -> RMBConfig:
    """Use the midpoint of a batched bracket as the crossing estimate."""
    n_gates = 2 * round(((lo + hi) / 2) / 2)
    n_gates = int(np.clip(n_gates, lo, hi))
    return settings.make_config(n_gates, ratio)


def ratio_budget_for(settings: Any, budget: Budget, crossings: list[RMBConfig]) -> Budget:
    """Cap later ratios without restricting the first anchor search."""
    if crossings and settings.max_hqc_per_ratio is not None:
        return Budget(remaining_hqc=min(settings.max_hqc_per_ratio, budget.remaining_hqc))
    return Budget(remaining_hqc=budget.remaining_hqc)


def search_one_ratio(
    rmb: RMB,
    rng: RNGGenerator,
    data: RMBData,
    crossings: list[RMBConfig],
    budget: Budget,
    settings: Any,
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
        lo, hi, source = trace_batch_window(data, crossings, settings, ratio, lo, hi)
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
