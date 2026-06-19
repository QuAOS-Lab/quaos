"""
Posterior-window batched monotone tracing of the fidelity = 0.5 boundary.

This variant probes one batch per ratio. After each batch it updates a particle
posterior over an inverse-form decay boundary from all measured Bernoulli
outcomes, then uses posterior uncertainty to choose the next lower-ratio batch
window.
"""
from __future__ import annotations

from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from statistics import NormalDist
from typing import Literal

import numpy as np
from numpy.random import Generator as RNGGenerator
from scipy.special import expit

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.batched_tracing_helpers.hints import (
    initial_prediction_bracket,
    measured_config_count,
)
from sympleq.applications.randomized_benchmarking.experiments.batched_tracing_helpers.windows import (
    batch_window_bracket,
    midpoint_crossing_from_bracket,
)
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    measured_items,
    print_crossing,
    print_experiment_summary,
    print_progress,
    save_crossings,
    start_run,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    REFERENCE_OFFSET,
    REFERENCE_SLOPE,
    analytic_gate_counts,
    parametric_bootstrap_analytic_coverage,
    parametric_boundary_fit_score,
)
from sympleq.core.noise.noise_model import GenericNoise


BASE_1Q_PAULI_ERROR = 0.000025
BASE_2Q_PAULI_ERROR = 0.00079


def sympleq_backend_factory(
    settings: CrossingSettings,
    rng: RNGGenerator,
) -> RMBBackend:
    """Local SympleQ backend with the default benchmarking noise model."""
    noise_model = GenericNoise.from_paulis([BASE_1Q_PAULI_ERROR] * 3, rng)
    two_qubit_noise_model = GenericNoise.from_paulis([BASE_2Q_PAULI_ERROR] * 3, rng)
    return SympleqBackend(
        noise_model=noise_model,
        two_qubit_noise_model=two_qubit_noise_model,
    )


@dataclass(frozen=True)
class BootstrapBatchedMonotoneTracingSettings(CrossingSettings):
    """
    Settings for posterior-window batched contour tracing.

    The first anchor window is centered on an analytic prior computed from
    ``initial_one_q_pauli_error`` and ``initial_two_q_pauli_error``. Its width
    is the analytic contour interval implied by
    ``initial_error_relative_uncertainty`` at ``initial_window_confidence``.
    Subsequent fixed-ratio windows are chosen from posterior samples of the
    inverse-form decay boundary fitted to all measured data so far. A confident
    bracket is reported when one is seen, but the next window is driven by the
    updated decay-rate posterior rather than by accepted crossings.
    """

    save_path: str | Path | None = "bootstrap_batched_monotone_tracing.json"
    backend_factory: Callable[[CrossingSettings, RNGGenerator], RMBBackend] = (
        sympleq_backend_factory
    )

    ratio_step: float = 0.1
    start_ratio: float | None = None
    trace_direction: Literal["up", "down", "both"] = "down"
    adaptive_ratio_step: bool = True
    min_ratio_step: float = 0.1
    max_ratio_step: float = 0.2
    target_hqc_per_ratio: float = 35.0

    decision_confidence: float = 0.70
    initial_bracket_fraction: float = 0.35
    low_to_high_growth_factor: float = 1.7
    use_surface_bracket_hint: bool = True
    surface_min_configs: int = 12
    initial_prediction_window_fraction: float = 0.35
    initial_prediction_min_width: int = 80

    initial_one_q_pauli_error: float = BASE_1Q_PAULI_ERROR
    initial_two_q_pauli_error: float = BASE_2Q_PAULI_ERROR
    initial_error_relative_uncertainty: float = 0.30
    initial_one_q_error_relative_uncertainty: float | None = None
    initial_two_q_error_relative_uncertainty: float | None = None
    initial_window_confidence: float = 0.50
    initial_window_min_width: int = 100
    initial_batch_points: int = 6
    initial_batch_shots: int = 4

    trace_batch_points: int = 6
    trace_batch_shots: int = 4
    max_hqc_per_ratio: float | None = 35.0

    particle_count: int = 10000
    particle_sharpness: float = 1.0
    particle_sharpness_log_std: float = 0.75
    particle_seed_offset: int = 800_000
    bootstrap_window_quantiles: tuple[float, float] = (0.10, 0.90)
    posterior_probe_quantiles: tuple[float, ...] = (0.10, 0.30, 0.50, 0.70, 0.90)
    bootstrap_window_min_valid_fraction: float = 0.35
    bootstrap_window_padding_fraction: float = 0.15
    bootstrap_window_min_width: int = 40
    bootstrap_window_max_width_fraction: float | None = 0.45
    bootstrap_min_configs: int = 3

    def boundary_fit(self, data) -> tuple[float, float, float] | None:
        """Prior-aware boundary fit used by plots and scores."""
        return decay_boundary_fit(data, self)

    def boundary_bootstrap(
        self,
        data,
        *,
        n_bootstrap: int = 100,
        seed: int | None = None,
    ) -> list[tuple[float, float, float]]:
        """Prior-aware boundary bootstrap used by plots and scores."""
        return decay_boundary_bootstrap(
            data,
            self,
            n_bootstrap=n_bootstrap,
            seed=seed,
        )


def bootstrap_window_seed(
    settings: BootstrapBatchedMonotoneTracingSettings,
    ratio: float,
    step_index: int,
) -> int | None:
    """Stable per-ratio posterior seed when the run itself is seeded."""
    if settings.rng_seed is None:
        return None
    ratio_key = int(round(1_000_000 * ratio))
    return settings.rng_seed + settings.particle_seed_offset + ratio_key + 997 * step_index


def resolved_start_ratio(settings: BootstrapBatchedMonotoneTracingSettings) -> float:
    """Starting ratio for the configured trace direction."""
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
    """Inclusive ratio sweep from start to stop."""
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
    settings: BootstrapBatchedMonotoneTracingSettings,
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


def ratio_sweeps(settings: BootstrapBatchedMonotoneTracingSettings) -> list[list[float]]:
    """Non-adaptive ratio sweeps for the configured trace direction."""
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
    settings: BootstrapBatchedMonotoneTracingSettings,
    budget: Budget,
) -> Iterator[Iterator[float] | list[float]]:
    """Ratio iterators for posterior-window tracing."""
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


def initial_error_relative_uncertainties(
    settings: BootstrapBatchedMonotoneTracingSettings,
) -> tuple[float, float]:
    """Independent relative uncertainty for the initial one- and two-qubit errors."""
    one_q = (
        settings.initial_error_relative_uncertainty
        if settings.initial_one_q_error_relative_uncertainty is None
        else settings.initial_one_q_error_relative_uncertainty
    )
    two_q = (
        settings.initial_error_relative_uncertainty
        if settings.initial_two_q_error_relative_uncertainty is None
        else settings.initial_two_q_error_relative_uncertainty
    )
    if one_q < 0.0 or two_q < 0.0:
        raise ValueError("initial error relative uncertainties must be non-negative.")
    return one_q, two_q


def initial_decay_prior(settings: BootstrapBatchedMonotoneTracingSettings) -> tuple[float, float]:
    """Prior center for ``n(r) = 1 / (q + slope * r)`` from initial errors."""
    one_q_scale = settings.initial_one_q_pauli_error / BASE_1Q_PAULI_ERROR
    two_q_scale = settings.initial_two_q_pauli_error / BASE_2Q_PAULI_ERROR
    return (
        REFERENCE_OFFSET * one_q_scale / np.log(2.0),
        REFERENCE_SLOPE * two_q_scale / np.log(2.0),
    )


def measured_count_arrays(data) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Measured configs as ratio, gate count, success count, failure count arrays."""
    measured = measured_items(data)
    ratios = np.asarray([config.ratio_2_qb_gates for config, _ in measured], dtype=float)
    n_gates = np.asarray([config.n_gates for config, _ in measured], dtype=float)
    successes = []
    failures = []
    for _, estimator in measured:
        counts = estimator.counts()
        successes.append(float(counts.get(True, 0)))
        failures.append(float(counts.get(False, 0)))
    return (
        ratios,
        n_gates,
        np.asarray(successes, dtype=float),
        np.asarray(failures, dtype=float),
    )


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: float | np.ndarray,
) -> np.ndarray:
    """Weighted quantile with normalized non-negative weights."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    quantiles_array = np.asarray(quantiles, dtype=float)
    finite = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(finite):
        return np.full_like(quantiles_array, np.nan, dtype=float)

    values = values[finite]
    weights = weights[finite]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cumulative /= cumulative[-1]
    result = np.interp(quantiles_array, cumulative, values)
    return result


def prior_decay_particles(
    settings: BootstrapBatchedMonotoneTracingSettings,
    *,
    seed: int | None,
    n_particles: int | None = None,
) -> np.ndarray:
    """Sample prior particles for q, slope, and logistic sharpness."""
    count = max(1, settings.particle_count if n_particles is None else n_particles)
    rng = np.random.default_rng(seed)
    q0, slope0 = initial_decay_prior(settings)
    one_q_uncertainty, two_q_uncertainty = initial_error_relative_uncertainties(settings)
    log_center = np.log([
        max(q0, 1e-12),
        max(slope0, 1e-12),
        max(settings.particle_sharpness, 1e-6),
    ])
    log_sigma = np.array(
        [
            max(0.05, np.log1p(one_q_uncertainty)),
            max(0.05, np.log1p(two_q_uncertainty)),
            max(0.05, settings.particle_sharpness_log_std),
        ],
        dtype=float,
    )
    particles = rng.normal(log_center, log_sigma, size=(count, 3))
    particles = np.exp(particles)
    particles[:, 0] = np.clip(particles[:, 0], 1e-8, 1.0)
    particles[:, 1] = np.clip(particles[:, 1], 1e-8, 1.0)
    particles[:, 2] = np.clip(particles[:, 2], 0.05, 100.0)
    return particles


def particle_log_likelihood(
    particles: np.ndarray,
    ratios: np.ndarray,
    n_gates: np.ndarray,
    successes: np.ndarray,
    failures: np.ndarray,
) -> np.ndarray:
    """Bernoulli log likelihood for all particles against measured outcomes."""
    if len(ratios) == 0:
        return np.zeros(len(particles), dtype=float)

    q = particles[:, 0][:, None]
    slope = particles[:, 1][:, None]
    sharpness = particles[:, 2][:, None]
    boundary = 1.0 / (q + slope * ratios[None, :])
    probabilities = expit(sharpness * (boundary - n_gates[None, :]) / boundary)
    eps = 1e-12
    return np.sum(
        successes[None, :] * np.log(probabilities + eps)
        + failures[None, :] * np.log(1.0 - probabilities + eps),
        axis=1,
    )


def decay_particle_posterior(
    data,
    settings: BootstrapBatchedMonotoneTracingSettings,
    *,
    seed: int | None,
    n_particles: int | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Particle posterior over inverse-boundary decay parameters."""
    particles = prior_decay_particles(settings, seed=seed, n_particles=n_particles)
    ratios, n_gates, successes, failures = measured_count_arrays(data)
    log_weights = particle_log_likelihood(
        particles,
        ratios,
        n_gates,
        successes,
        failures,
    )
    log_weights -= np.max(log_weights)
    weights = np.exp(log_weights)
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0.0:
        weights = np.full(len(particles), 1.0 / len(particles), dtype=float)
    else:
        weights /= total
    ess = 1.0 / float(np.sum(weights**2))
    return particles, weights, ess


def decay_boundary_bootstrap(
    data,
    settings: BootstrapBatchedMonotoneTracingSettings,
    *,
    n_bootstrap: int,
    seed: int | None,
) -> list[tuple[float, float, float]]:
    """Posterior draws of inverse-boundary parameters for plots and scores."""
    if measured_config_count(data) < settings.bootstrap_min_configs:
        return []
    rng = np.random.default_rng(seed)
    particles, weights, _ = decay_particle_posterior(
        data,
        settings,
        seed=seed,
        n_particles=max(settings.particle_count, n_bootstrap),
    )
    indices = rng.choice(len(particles), size=n_bootstrap, replace=True, p=weights)
    return [tuple(float(x) for x in particles[index]) for index in indices]


def decay_boundary_fit(
    data,
    settings: BootstrapBatchedMonotoneTracingSettings,
) -> tuple[float, float, float] | None:
    """Posterior-median inverse-boundary fit from measured Bernoulli counts."""
    if measured_config_count(data) < settings.bootstrap_min_configs:
        return None
    particles, weights, _ = decay_particle_posterior(
        data,
        settings,
        seed=bootstrap_window_seed(settings, 0.0, measured_config_count(data)),
    )
    medians = np.asarray(
        [weighted_quantile(particles[:, index], weights, 0.5) for index in range(3)],
        dtype=float,
    )
    if np.any(~np.isfinite(medians)):
        return None
    return tuple(float(value) for value in medians)


def decay_gate_samples(
    fits: list[tuple[float, float, float]],
    ratio: float,
) -> np.ndarray:
    """Boundary gate-count samples at one ratio."""
    values = []
    for q, slope, _ in fits:
        gate_count = 1.0 / (q + slope * ratio)
        if np.isfinite(gate_count) and gate_count > 0.0:
            values.append(float(gate_count))
    return np.asarray(values, dtype=float)


def posterior_probe_gates(
    samples: np.ndarray,
    weights: np.ndarray,
    settings: BootstrapBatchedMonotoneTracingSettings,
    lo: int,
    hi: int,
    n_points: int,
) -> list[int]:
    """Probe posterior quantiles, filling gaps if duplicate rounding occurs."""
    quantiles = settings.posterior_probe_quantiles
    if len(quantiles) != n_points:
        quantiles = tuple(np.linspace(
            settings.bootstrap_window_quantiles[0],
            settings.bootstrap_window_quantiles[1],
            n_points,
        ))
    gates = [
        2 * round(float(weighted_quantile(samples, weights, quantile)) / 2)
        for quantile in quantiles
    ]
    gates = [int(np.clip(gate, lo, hi)) for gate in gates]
    unique = sorted(set(gates))
    if len(unique) >= n_points:
        return unique[:n_points]

    fill = [2 * round(value / 2) for value in np.linspace(lo, hi, n_points)]
    unique = sorted(set([*unique, *(int(np.clip(gate, lo, hi)) for gate in fill)]))
    return unique[:n_points]


def decay_model_window(
    data,
    settings: BootstrapBatchedMonotoneTracingSettings,
    ratio: float,
    step_index: int,
) -> tuple[int, int, str, list[int]] | None:
    """
    Choose a batch window from the updated inverse-decay posterior.

    Particles are drawn from the configured analytic-error prior and weighted
    by the Bernoulli likelihood of every measurement so far. Evaluating those
    weighted particles at the next ratio gives the gate-count distribution that
    defines the next batch window.
    """
    if measured_config_count(data) < settings.bootstrap_min_configs:
        return None
    q_low, q_high = settings.bootstrap_window_quantiles
    if not 0.0 <= q_low < q_high <= 1.0:
        raise ValueError(
            "bootstrap_window_quantiles must satisfy 0 <= low < high <= 1."
        )

    particles, weights, ess = decay_particle_posterior(
        data,
        settings,
        seed=bootstrap_window_seed(settings, ratio, step_index),
    )
    samples = 1.0 / (particles[:, 0] + particles[:, 1] * ratio)
    finite = np.isfinite(samples) & (samples > 0.0) & np.isfinite(weights) & (weights > 0.0)
    valid_weight = float(np.sum(weights[finite]))
    valid_fraction = valid_weight / max(float(np.sum(weights)), 1e-12)
    if valid_fraction < settings.bootstrap_window_min_valid_fraction:
        return None

    finite_samples = samples[finite]
    finite_weights = weights[finite] / valid_weight
    lower = float(weighted_quantile(finite_samples, finite_weights, q_low))
    upper = float(weighted_quantile(finite_samples, finite_weights, q_high))
    median = float(weighted_quantile(finite_samples, finite_weights, 0.5))
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        return None

    width = upper - lower
    pad = max(2.0, settings.bootstrap_window_padding_fraction * width)
    lower -= pad
    upper += pad

    min_width = max(2, settings.bootstrap_window_min_width)
    if upper - lower < min_width:
        center = median
        lower = center - 0.5 * min_width
        upper = center + 0.5 * min_width

    if settings.bootstrap_window_max_width_fraction is not None:
        max_width = max(
            min_width,
            settings.bootstrap_window_max_width_fraction * max(median, 1.0),
        )
        if upper - lower > max_width:
            center = median
            lower = center - 0.5 * max_width
            upper = center + 0.5 * max_width

    n_gates_min, n_gates_max = settings.n_gates_bounds
    window_lo = max(n_gates_min, 2 * round(lower / 2))
    window_hi = min(n_gates_max, 2 * round(upper / 2))
    if window_hi <= window_lo:
        return None
    probes = posterior_probe_gates(
        finite_samples,
        finite_weights,
        settings,
        int(window_lo),
        int(window_hi),
        settings.trace_batch_points,
    )

    source = (
        f"particle posterior q{q_low:.2f}-q{q_high:.2f} "
        f"ess={ess:.0f}/{len(particles)} valid={valid_fraction:.2f}"
    )
    return int(window_lo), int(window_hi), source, probes


def initial_analytic_gate_count(
    settings: BootstrapBatchedMonotoneTracingSettings,
    ratio: float,
    *,
    one_q_pauli_error: float,
    two_q_pauli_error: float,
) -> float | None:
    """Analytic total-gate crossing for explicit one- and two-qubit errors."""
    if one_q_pauli_error <= 0.0 or two_q_pauli_error <= 0.0:
        return None
    prediction = analytic_gate_counts(
        np.array([ratio]),
        one_q_noise_scale=one_q_pauli_error / BASE_1Q_PAULI_ERROR,
        two_q_noise_scale=two_q_pauli_error / BASE_2Q_PAULI_ERROR,
    )[0]
    if not np.isfinite(prediction) or prediction <= 0.0:
        return None
    return float(prediction)


def initial_analytic_prior_window(
    settings: BootstrapBatchedMonotoneTracingSettings,
    ratio: float,
) -> tuple[int, int] | None:
    """
    First-anchor gate window from uncertainty in the analytic prior errors.

    The error inputs are treated as central estimates with independent relative
    normal uncertainty. Since the analytic crossing is monotone decreasing in
    both error rates, the low-gate side comes from the high-error contour and
    the high-gate side from the low-error contour.
    """
    if not 0.0 < settings.initial_window_confidence < 1.0:
        raise ValueError("initial_window_confidence must be between 0 and 1.")
    one_q_uncertainty, two_q_uncertainty = initial_error_relative_uncertainties(settings)

    center = initial_analytic_gate_count(
        settings,
        ratio,
        one_q_pauli_error=settings.initial_one_q_pauli_error,
        two_q_pauli_error=settings.initial_two_q_pauli_error,
    )
    if center is None:
        return None

    z_score = NormalDist().inv_cdf(0.5 + 0.5 * settings.initial_window_confidence)
    one_q_radius = z_score * one_q_uncertainty
    two_q_radius = z_score * two_q_uncertainty
    lower_one_q_factor = max(1e-9, 1.0 - one_q_radius)
    upper_one_q_factor = 1.0 + one_q_radius
    lower_two_q_factor = max(1e-9, 1.0 - two_q_radius)
    upper_two_q_factor = 1.0 + two_q_radius

    high_error_gate_count = initial_analytic_gate_count(
        settings,
        ratio,
        one_q_pauli_error=settings.initial_one_q_pauli_error * upper_one_q_factor,
        two_q_pauli_error=settings.initial_two_q_pauli_error * upper_two_q_factor,
    )
    low_error_gate_count = initial_analytic_gate_count(
        settings,
        ratio,
        one_q_pauli_error=settings.initial_one_q_pauli_error * lower_one_q_factor,
        two_q_pauli_error=settings.initial_two_q_pauli_error * lower_two_q_factor,
    )
    if high_error_gate_count is None or low_error_gate_count is None:
        return None

    lo_float = min(high_error_gate_count, low_error_gate_count, center)
    hi_float = max(high_error_gate_count, low_error_gate_count, center)
    if hi_float - lo_float < settings.initial_window_min_width:
        half_width = 0.5 * settings.initial_window_min_width
        lo_float = center - half_width
        hi_float = center + half_width

    n_gates_min, n_gates_max = settings.n_gates_bounds
    lo = max(n_gates_min, 2 * round(lo_float / 2))
    hi = min(n_gates_max, 2 * round(hi_float / 2))
    if hi <= lo:
        return None
    return int(lo), int(hi)


def initial_anchor_window(
    data,
    settings: BootstrapBatchedMonotoneTracingSettings,
    ratio: float,
) -> tuple[int, int, str]:
    """Choose the first-anchor window from the configured analytic prior."""
    window = initial_analytic_prior_window(settings, ratio)
    if window is None:
        return initial_prediction_bracket(data, settings, ratio)
    lo, hi = window
    return lo, hi, (
        f"initial analytic prior {100.0 * settings.initial_window_confidence:.0f}%"
    )


def search_one_ratio(
    rmb: RMB,
    rng: RNGGenerator,
    data,
    budget: Budget,
    settings: BootstrapBatchedMonotoneTracingSettings,
    ratio: float,
) -> RMBConfig | None:
    """Probe one ratio using the current decay-rate posterior window."""
    ratio_budget = Budget(
        remaining_hqc=(
            budget.remaining_hqc
            if measured_config_count(data) == 0 or settings.max_hqc_per_ratio is None
            else min(settings.max_hqc_per_ratio, budget.remaining_hqc)
        )
    )
    if measured_config_count(data) == 0:
        lo, hi, source = initial_anchor_window(data, settings, ratio)
        n_points = settings.initial_batch_points
        shots = settings.initial_batch_shots
        probe_gates = None
    else:
        window = decay_model_window(data, settings, ratio, measured_config_count(data))
        if window is None:
            lo, hi, source = initial_anchor_window(data, settings, ratio)
            probe_gates = None
        else:
            lo, hi, source, probe_gates = window
        n_points = settings.trace_batch_points
        shots = settings.trace_batch_shots

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
        n_points=n_points,
        shots=shots,
        probe_gates=probe_gates,
    )
    if bracket is None:
        return None

    lo, hi = bracket
    return midpoint_crossing_from_bracket(settings, ratio, lo, hi)


def run_with_budget(
    settings: BootstrapBatchedMonotoneTracingSettings,
) -> tuple[RMB, list[RMBConfig], Budget]:
    """Run the posterior-window batched boundary tracer."""
    rng, rmb, budget = start_run(settings)
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


def run(
    settings: BootstrapBatchedMonotoneTracingSettings,
) -> tuple[RMB, list[RMBConfig]]:
    """Run the posterior-window batched boundary tracer."""
    rmb, crossings, _ = run_with_budget(settings)
    return rmb, crossings


def print_run_score(
    rmb: RMB,
    crossings: list[RMBConfig],
    budget: Budget,
    settings: BootstrapBatchedMonotoneTracingSettings,
) -> None:
    """Print the direct-run benchmark metrics for one posterior-window run."""
    score = parametric_boundary_fit_score(rmb._data, settings)
    coverage50, coverage90 = parametric_bootstrap_analytic_coverage(
        rmb._data,
        settings,
        seed=None if settings.rng_seed is None else settings.rng_seed + 500_000,
    )
    print("\nRun score")
    print(f"  score: {score:.3f}")
    print(f"  analytic line inside bootstrap 50% band: {100.0 * coverage50:.1f}%")
    print(f"  analytic line inside bootstrap 90% band: {100.0 * coverage90:.1f}%")
    print(f"  spent: {budget.spent_hqc:.1f} HQC")
    print(f"  crossings: {len(crossings)}")


if __name__ == "__main__":
    settings = BootstrapBatchedMonotoneTracingSettings(rng_seed=None)
    rmb, crossings, budget = run_with_budget(settings)
