from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from scipy.optimize import differential_evolution, minimize
from scipy.special import expit

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData

from viarregio2 import (
    BoundaryExperimentConfig,
    adaptive_shot_count,
    budgeted_initial_design,
    config_from_parameters,
    fidelity_mean,
    fidelity_variance,
    grouped_by_n_qubits,
    make_backend,
    nearest_scaled_distance,
    random_candidate_configs,
    scaled_config_point,
    scaled_data_points,
    spend_measurements,
    total_measurements,
)
from viarregio3 import (
    contour_follow_candidates,
    contour_points_from_surface,
    near_contour_anchors,
)


@dataclass(frozen=True)
class MonotoneBoundaryExperimentConfig(BoundaryExperimentConfig):
    """
    Hybrid boundary-estimation settings with a monotone fidelity model.

    The fitted surface is constrained so that expected fidelity cannot increase
    as either depth or two-qubit gate ratio increases.
    """
    monotone_l2: float = 1e-3
    contour_follow_fraction: float = 0.50
    contour_ready_probability_width: float = 0.10
    contour_min_anchors: int = 4
    contour_step_fraction: float = 0.08
    contour_projection_fraction: float = 0.12
    contour_gradient_fraction: float = 0.01
    contour_candidate_multiplier: int = 5
    adaptive_batch_size: bool = True
    min_batch_size: int = 1
    late_contour_follow_fraction: float = 0.85
    late_exploration_weight: float = 0.05
    adaptive_batch_late_budget_fraction: float = 0.50
    adaptive_batch_anchor_multiplier: int = 2
    hqc_budget: float | None = None
    hqc_base_cost: float = 5.0
    hqc_one_qubit_weight: float = 1.0
    hqc_two_qubit_weight: float = 10.0
    hqc_measurement_weight: float = 5.0
    hqc_scale: float = 5000.0
    hqc_cost_informed_acquisition: bool = True
    hqc_cost_power: float = 1.0
    save_path: str | Path | None = "viarregio4_boundary.json"


@dataclass
class MonotoneFidelitySurface:
    alpha: float
    coefficients: np.ndarray
    depth_bounds: tuple[int, int]
    ratio_bounds: tuple[float, float]
    n_points: int

    def _scale(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        lower = np.array([self.depth_bounds[0], self.ratio_bounds[0]], dtype=float)
        upper = np.array([self.depth_bounds[1], self.ratio_bounds[1]], dtype=float)
        return np.clip((points - lower) / (upper - lower), 0.0, 1.0)

    def features(self, points: np.ndarray) -> np.ndarray:
        scaled = self._scale(points)
        depth = scaled[:, 0]
        ratio = scaled[:, 1]
        return np.column_stack([
            depth,
            ratio,
            depth**2,
            ratio**2,
            depth * ratio,
            depth**3,
            ratio**3,
            depth**2 * ratio,
            depth * ratio**2,
        ])

    def probability(self, points: np.ndarray) -> np.ndarray:
        features = self.features(points)
        return expit(self.alpha - features @ self.coefficients)

    def probability_grid(
        self,
        settings: BoundaryExperimentConfig,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n_depths, n_ratios = settings.candidate_grid_size
        depths = np.linspace(settings.depth_bounds[0], settings.depth_bounds[1], n_depths)
        ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], n_ratios)
        depth_grid, ratio_grid = np.meshgrid(depths, ratios)
        points = np.column_stack([depth_grid.ravel(), ratio_grid.ravel()])
        probabilities = self.probability(points).reshape(depth_grid.shape)
        return depth_grid, ratio_grid, probabilities

    def report(self, settings: BoundaryExperimentConfig) -> str:
        _, _, probabilities = self.probability_grid(settings)
        min_probability = float(np.min(probabilities))
        max_probability = float(np.max(probabilities))
        has_contour = min_probability <= 0.5 <= max_probability
        return (
            "Fitted monotone fidelity surface:\n"
            f"  fit points: {self.n_points}\n"
            f"  predicted fidelity range on grid: "
            f"{min_probability:.4f} to {max_probability:.4f}\n"
            f"  contains fidelity=0.5 contour: {has_contour}"
        )


def fit_monotone_surface_from_arrays(
    points: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    settings: MonotoneBoundaryExperimentConfig,
    n_points: int,
) -> MonotoneFidelitySurface:
    points = np.asarray(points, dtype=float)
    y = np.asarray(y, dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights = weights / np.mean(weights)

    probe_surface = MonotoneFidelitySurface(
        alpha=0.0,
        coefficients=np.zeros(9, dtype=float),
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=n_points,
    )
    features = probe_surface.features(points)

    def loss_and_grad(params: np.ndarray) -> tuple[float, np.ndarray]:
        alpha = params[0]
        coefficients = params[1:]
        z = alpha - features @ coefficients
        p = expit(z)

        eps = 1e-12
        loss = -np.sum(weights * (y * np.log(p + eps) + (1.0 - y) * np.log(1.0 - p + eps)))
        loss += 0.5 * settings.monotone_l2 * np.sum(coefficients**2)

        residual = weights * (p - y)
        grad_alpha = np.sum(residual)
        grad_coefficients = -features.T @ residual + settings.monotone_l2 * coefficients
        return loss, np.concatenate([[grad_alpha], grad_coefficients])

    x0 = np.zeros(features.shape[1] + 1, dtype=float)
    bounds = [(None, None)] + [(0.0, None)] * features.shape[1]
    result = minimize(
        fun=lambda params: loss_and_grad(params)[0],
        x0=x0,
        jac=lambda params: loss_and_grad(params)[1],
        bounds=bounds,
        method="L-BFGS-B",
    )
    if not result.success:
        raise RuntimeError(result.message)

    return MonotoneFidelitySurface(
        alpha=float(result.x[0]),
        coefficients=result.x[1:],
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=n_points,
    )


def fit_monotone_surface_from_counts(
    points: np.ndarray,
    successes: np.ndarray,
    failures: np.ndarray,
    settings: MonotoneBoundaryExperimentConfig,
    n_points: int,
) -> MonotoneFidelitySurface:
    """
    Fit the monotone surface with the binomial likelihood for Boolean outcomes.
    """
    points = np.asarray(points, dtype=float)
    successes = np.asarray(successes, dtype=float)
    failures = np.asarray(failures, dtype=float)
    if np.any(successes < 0.0) or np.any(failures < 0.0):
        raise ValueError("Success and failure counts must be non-negative.")

    probe_surface = MonotoneFidelitySurface(
        alpha=0.0,
        coefficients=np.zeros(9, dtype=float),
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=n_points,
    )
    features = probe_surface.features(points)

    def loss_and_grad(params: np.ndarray) -> tuple[float, np.ndarray]:
        alpha = params[0]
        coefficients = params[1:]
        z = alpha - features @ coefficients
        p = expit(z)

        eps = 1e-12
        loss = -np.sum(
            successes * np.log(p + eps)
            + failures * np.log(1.0 - p + eps)
        )
        loss += 0.5 * settings.monotone_l2 * np.sum(coefficients**2)

        residual = p * (successes + failures) - successes
        grad_alpha = np.sum(residual)
        grad_coefficients = -features.T @ residual + settings.monotone_l2 * coefficients
        return loss, np.concatenate([[grad_alpha], grad_coefficients])

    x0 = np.zeros(features.shape[1] + 1, dtype=float)
    bounds = [(None, None)] + [(0.0, None)] * features.shape[1]
    result = minimize(
        fun=lambda params: loss_and_grad(params)[0],
        x0=x0,
        jac=lambda params: loss_and_grad(params)[1],
        bounds=bounds,
        method="L-BFGS-B",
    )
    if not result.success:
        raise RuntimeError(result.message)

    return MonotoneFidelitySurface(
        alpha=float(result.x[0]),
        coefficients=result.x[1:],
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=n_points,
    )


def fit_monotone_fidelity_surface(
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
) -> MonotoneFidelitySurface:
    configs = [config for config, estimator in data.items() if estimator.num_runs() > 0]
    if len(configs) < settings.min_fit_points:
        raise ValueError(
            f"Need at least {settings.min_fit_points} data points to fit a surface, "
            f"got {len(configs)}."
        )

    points = np.array(
        [[float(config.depth), float(config.min_two_qubit_gate_ratio)] for config in configs],
        dtype=float,
    )
    successes = []
    failures = []
    for config in configs:
        counts = data[config].counts()
        successes.append(float(counts.get(True, 0)))
        failures.append(float(counts.get(False, 0)))

    return fit_monotone_surface_from_counts(
        points=points,
        successes=np.asarray(successes, dtype=float),
        failures=np.asarray(failures, dtype=float),
        settings=settings,
        n_points=len(configs),
    )


def fit_monotone_surface_from_values(
    configs: list[RMBConfig],
    fidelities: np.ndarray,
    weights: np.ndarray,
    settings: MonotoneBoundaryExperimentConfig,
) -> MonotoneFidelitySurface:
    points = np.array(
        [[float(config.depth), float(config.min_two_qubit_gate_ratio)] for config in configs],
        dtype=float,
    )
    return fit_monotone_surface_from_arrays(
        points=points,
        y=np.asarray(fidelities, dtype=float),
        weights=np.asarray(weights, dtype=float),
        settings=settings,
        n_points=len(configs),
    )


def bootstrap_contours(
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
) -> list[np.ndarray]:
    """
    Draw posterior monotone surfaces and return their p=0.5 contours.
    """
    rng = default_rng(seed)
    configs = [config for config, estimator in data.items() if estimator.num_runs() > 0]
    if len(configs) < settings.min_fit_points:
        return []

    alpha = []
    beta = []
    weights = []
    for config in configs:
        counts = data[config].counts()
        alpha.append(counts.get(True, 0) + 1.0)
        beta.append(counts.get(False, 0) + 1.0)
        weights.append(max(1, data[config].num_runs()))

    alpha_array = np.asarray(alpha, dtype=float)
    beta_array = np.asarray(beta, dtype=float)
    weights_array = np.asarray(weights, dtype=float)
    contours = []

    for _ in range(n_bootstrap):
        sampled_fidelities = rng.beta(alpha_array, beta_array)
        try:
            surface = fit_monotone_surface_from_values(
                configs,
                sampled_fidelities,
                weights_array,
                settings,
            )
        except (RuntimeError, ValueError):
            continue
        contour = contour_points_from_surface(surface, settings)
        if len(contour) > 0:
            contours.append(contour)

    return contours


def print_fit_reports(data: RMBData, settings: MonotoneBoundaryExperimentConfig) -> None:
    for n_qubits, group in sorted(grouped_by_n_qubits(data).items()):
        try:
            surface = fit_monotone_fidelity_surface(group, settings)
        except (RuntimeError, ValueError):
            print(f"n_qubits={n_qubits}: not enough data for a monotone surface fit yet.")
            continue

        print(f"\nn_qubits={n_qubits}")
        print(surface.report(settings))


def expected_gate_counts(config: RMBConfig) -> tuple[int, int, int]:
    """
    Expected (N_1, N_2, N_m) for the RB circuit represented by a config.

    `Circuit.from_depth` creates `n_qubits * depth` gate slots before the
    inverse. A two-qubit gate occupies two of those slots. The RMB construction
    appends the inverse circuit, so the random-circuit gate counts are doubled.
    Scrambler X gates are included in expectation. Identity scrambler gates are
    not charged as physical one-qubit gates.
    """
    n_slots = config.n_qubits * config.depth
    ratio = 0.5 * (config.min_two_qubit_gate_ratio + config.max_two_qubit_gate_ratio)
    two_qubit_base = int(ratio * n_slots) // 2
    one_qubit_base = n_slots - 2 * two_qubit_base
    scrambler_one_qubit = int(round(2 * config.n_qubits * config.scrambling_probability))
    n_one_qubit = 2 * one_qubit_base + scrambler_one_qubit
    n_two_qubit = 2 * two_qubit_base
    n_measurements = config.n_qubits
    return n_one_qubit, n_two_qubit, n_measurements


def hqc_cost(
    config: RMBConfig,
    shot_count: int,
    settings: MonotoneBoundaryExperimentConfig,
) -> float:
    """
    Hardware quantum credits for one circuit at `config` repeated `shot_count` times.
    """
    if shot_count <= 0:
        return 0.0
    n_one, n_two, n_meas = expected_gate_counts(config)
    weighted_size = (
        settings.hqc_one_qubit_weight * n_one
        + settings.hqc_two_qubit_weight * n_two
        + settings.hqc_measurement_weight * n_meas
    )
    return settings.hqc_base_cost + (weighted_size / settings.hqc_scale) * shot_count


def affordable_shot_count(
    config: RMBConfig,
    requested_shots: int,
    remaining_hqc: float,
    settings: MonotoneBoundaryExperimentConfig,
) -> int:
    """
    Clip a requested shot count to what the remaining HQC budget can pay for.
    """
    if settings.hqc_budget is None:
        return requested_shots
    if requested_shots <= 0 or remaining_hqc <= settings.hqc_base_cost:
        return 0

    n_one, n_two, n_meas = expected_gate_counts(config)
    weighted_size = (
        settings.hqc_one_qubit_weight * n_one
        + settings.hqc_two_qubit_weight * n_two
        + settings.hqc_measurement_weight * n_meas
    )
    unit_hqc = weighted_size / settings.hqc_scale
    if unit_hqc <= 0.0:
        return requested_shots

    max_affordable = int(np.floor((remaining_hqc - settings.hqc_base_cost) / unit_hqc))
    return max(0, min(requested_shots, max_affordable))


def planned_candidate_shots(
    config: RMBConfig,
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
) -> int:
    """
    Estimate the shot count that would be requested for acquisition-cost scoring.
    """
    return max(
        1,
        adaptive_shot_count(
            config=config,
            data=data,
            settings=settings,
            remaining=settings.measurement_budget,
        ),
    )


def hqc_efficiency_multiplier(
    config: RMBConfig,
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
) -> float:
    """
    Penalize acquisition scores by expected HQC cost.

    The baseline is the fixed circuit cost, so a minimum-cost circuit receives
    a multiplier near one and larger/deeper circuits receive smaller values.
    """
    if not settings.hqc_cost_informed_acquisition:
        return 1.0
    shot_count = planned_candidate_shots(config, data, settings)
    cost = hqc_cost(config, shot_count, settings)
    if cost <= 0.0:
        return 1.0
    normalized_cost = max(1.0, cost / max(1e-12, settings.hqc_base_cost))
    return float(normalized_cost ** (-settings.hqc_cost_power))


def sort_configs_by_hqc_efficiency(
    configs: list[RMBConfig],
    settings: MonotoneBoundaryExperimentConfig,
    shot_count: int,
) -> list[RMBConfig]:
    """
    Prefer cheaper configurations when no fitted surface can score information.
    """
    if not settings.hqc_cost_informed_acquisition:
        return configs
    return sorted(
        configs,
        key=lambda config: (
            hqc_cost(config, max(1, shot_count), settings),
            config.n_qubits,
            config.depth,
            config.min_two_qubit_gate_ratio,
        ),
    )


def hqc_aware_acquisition_score(
    *,
    theta: np.ndarray,
    surface: MonotoneFidelitySurface,
    existing_points: np.ndarray,
    selected_points: list[np.ndarray],
    data: RMBData,
    template: RMBConfig,
    settings: MonotoneBoundaryExperimentConfig,
) -> float:
    """
    Boundary acquisition score normalized by expected HQC cost.
    """
    depth = float(theta[0])
    ratio = float(theta[1])
    if not settings.depth_bounds[0] <= depth <= settings.depth_bounds[1]:
        return 0.0
    if not settings.ratio_bounds[0] <= ratio <= settings.ratio_bounds[1]:
        return 0.0

    probability = float(surface.probability(np.array([[depth, ratio]]))[0])
    boundary_distance = abs(probability - 0.5)
    if boundary_distance > settings.max_boundary_probability_distance:
        return 0.0
    boundary_relevance = np.exp(-((boundary_distance / settings.boundary_width) ** 2))
    boundary_relevance = boundary_relevance ** settings.boundary_focus_power

    scaled_point = scaled_config_point(depth=depth, ratio=ratio, settings=settings)
    data_distance = nearest_scaled_distance(scaled_point, existing_points)
    sparsity = min(1.0, data_distance / settings.sparsity_radius)

    if selected_points:
        selected_distance = nearest_scaled_distance(scaled_point, np.asarray(selected_points))
        diversity = min(1.0, selected_distance / settings.batch_diversity_radius)
    else:
        diversity = 1.0

    config = config_from_parameters(template=template, depth=depth, ratio=ratio)
    n_runs = data.get(config).num_runs() if config in data else 0
    if n_runs >= settings.max_shots_per_config:
        return 0.0
    undersampled = 1.0 - n_runs / max(1, settings.max_shots_per_config)

    sparsity_multiplier = 1.0 + settings.exploration_weight * sparsity
    diversity_multiplier = settings.diversity_floor + (1.0 - settings.diversity_floor) * diversity
    information_score = boundary_relevance * sparsity_multiplier * diversity_multiplier * undersampled
    return float(information_score * hqc_efficiency_multiplier(config, data, settings))


def optimize_hqc_aware_candidate_config(
    *,
    surface: MonotoneFidelitySurface,
    existing_points: np.ndarray,
    selected_points: list[np.ndarray],
    data: RMBData,
    template: RMBConfig,
    settings: MonotoneBoundaryExperimentConfig,
    rng: RNGGenerator,
) -> tuple[float, RMBConfig, np.ndarray]:
    bounds = [
        settings.depth_bounds,
        settings.ratio_bounds,
    ]

    def objective(theta: np.ndarray) -> float:
        score = hqc_aware_acquisition_score(
            theta=theta,
            surface=surface,
            existing_points=existing_points,
            selected_points=selected_points,
            data=data,
            template=template,
            settings=settings,
        )
        return -score

    best_score = -np.inf
    best_theta: np.ndarray | None = None
    for _ in range(settings.optimizer_attempts_per_config):
        seed = int(rng.integers(0, np.iinfo(np.int32).max))
        result = differential_evolution(
            objective,
            bounds=bounds,
            seed=seed,
            maxiter=settings.optimizer_maxiter,
            popsize=settings.optimizer_popsize,
            polish=True,
            updating="immediate",
            workers=1,
        )
        score = -float(result.fun)
        if score > best_score:
            best_score = score
            best_theta = np.asarray(result.x, dtype=float)

    if best_theta is None:
        raise RuntimeError("Candidate optimizer did not return a point.")

    config = config_from_parameters(
        template=template,
        depth=float(best_theta[0]),
        ratio=float(best_theta[1]),
    )
    rounded_point = scaled_config_point(
        depth=float(config.depth),
        ratio=float(config.min_two_qubit_gate_ratio),
        settings=settings,
    )
    return best_score, config, rounded_point


def adaptive_runtime_settings(
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
    remaining: float,
) -> MonotoneBoundaryExperimentConfig:
    """
    Derive the current batch behavior from boundary maturity diagnostics.

    Early updates keep the configured exploratory batch.  Once a fitted contour
    exists and enough near-boundary anchors have been measured, later updates
    shrink the batch and spend a larger fraction of it following the contour.
    """
    if not settings.adaptive_batch_size:
        return settings

    base_batch_size = max(1, settings.batch_size)
    min_batch_size = max(1, min(settings.min_batch_size, base_batch_size))
    budget = settings.hqc_budget if settings.hqc_budget is not None else settings.measurement_budget
    spent = budget - remaining
    progress = spent / max(1.0, float(budget))

    contour_found = False
    max_anchor_count = 0
    groups = grouped_by_n_qubits(data)
    for n_qubits in settings.n_qubits_values:
        group = groups.get(n_qubits, {})
        try:
            surface = fit_monotone_fidelity_surface(group, settings)
        except (RuntimeError, ValueError):
            continue

        _, _, probabilities = surface.probability_grid(settings)
        contour_found = contour_found or (
            float(np.min(probabilities)) <= 0.5 <= float(np.max(probabilities))
        )
        anchors = near_contour_anchors(group, surface, settings)
        max_anchor_count = max(max_anchor_count, len(anchors))

    if not contour_found:
        return settings

    ready_anchor_count = settings.contour_min_anchors
    mature_anchor_count = max(
        ready_anchor_count,
        settings.adaptive_batch_anchor_multiplier * ready_anchor_count,
    )
    contour_fraction = settings.contour_follow_fraction
    exploration_weight = settings.exploration_weight
    batch_size = base_batch_size

    if max_anchor_count >= mature_anchor_count and progress >= settings.adaptive_batch_late_budget_fraction:
        batch_size = min_batch_size
        contour_fraction = max(contour_fraction, settings.late_contour_follow_fraction)
        exploration_weight = min(exploration_weight, settings.late_exploration_weight)
    elif max_anchor_count >= ready_anchor_count:
        batch_size = max(min_batch_size, int(np.ceil(0.5 * (base_batch_size + min_batch_size))))
        contour_fraction = max(contour_fraction, 0.65)
        exploration_weight = min(exploration_weight, 0.10)

    return replace(
        settings,
        batch_size=batch_size,
        contour_follow_fraction=contour_fraction,
        exploration_weight=exploration_weight,
    )


def print_experiment_summary(
    *,
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
    stop_reason: str,
    batch_count: int,
    circuit_executions: int,
    max_execution_repeats: int,
    remaining_measurements: int,
    remaining_hqc: float,
) -> None:
    measured_configs = [
        config
        for config, estimator in data.items()
        if estimator.num_runs() > 0
    ]
    repeats = [
        data[config].num_runs()
        for config in measured_configs
    ]
    total_repeats = total_measurements(data)
    max_repeats_same_config = max(repeats, default=0)
    mean_repeats = float(np.mean(repeats)) if repeats else 0.0

    if measured_configs:
        depths = [config.depth for config in measured_configs]
        ratios = [config.min_two_qubit_gate_ratio for config in measured_configs]
        n_qubits = sorted({config.n_qubits for config in measured_configs})
        n_one_total = 0
        n_two_total = 0
        n_meas_total = 0
        for config in measured_configs:
            n_one, n_two, n_meas = expected_gate_counts(config)
            config_repeats = data[config].num_runs()
            n_one_total += n_one * config_repeats
            n_two_total += n_two * config_repeats
            n_meas_total += n_meas * config_repeats
    else:
        depths = []
        ratios = []
        n_qubits = []
        n_one_total = n_two_total = n_meas_total = 0

    print("\nExperiment summary")
    print(f"  stop reason: {stop_reason}")
    print(f"  batch proposals made: {batch_count}")
    print(f"  circuit executions: {circuit_executions}")
    print(f"  distinct configs measured: {len(measured_configs)}")
    print(f"  total repeats / Boolean outcomes: {total_repeats}")
    print(f"  max repeats in one execution: {max_execution_repeats}")
    print(f"  max repeats of the same config: {max_repeats_same_config}")
    print(f"  mean repeats per measured config: {mean_repeats:.2f}")

    if settings.hqc_budget is not None:
        hqc_spent = settings.hqc_budget - remaining_hqc
        print(f"  HQC spent: {hqc_spent:.3f} / {settings.hqc_budget:.3f}")
        print(f"  HQC remaining: {max(0.0, remaining_hqc):.3f}")
    print(
        f"  measurement cap used: {total_repeats} / {settings.measurement_budget} "
        f"(remaining {max(0, remaining_measurements)})"
    )

    if measured_configs:
        print(f"  n_qubits sampled: {n_qubits}")
        print(f"  sampled depth range: {min(depths)} to {max(depths)}")
        print(f"  sampled two-qubit ratio range: {min(ratios):.2f} to {max(ratios):.2f}")
        print(f"  expected one-qubit gate applications over repeats: {n_one_total}")
        print(f"  expected two-qubit gate applications over repeats: {n_two_total}")
        print(f"  expected measurement operations over repeats: {n_meas_total}")

    print(f"  adaptive batch size enabled: {settings.adaptive_batch_size}")
    print(f"  HQC-informed acquisition enabled: {settings.hqc_cost_informed_acquisition}")
    print(f"  HQC cost power: {settings.hqc_cost_power:.2f}")
    print(f"  base batch size: {settings.batch_size}")
    print(f"  min runtime batch size: {settings.min_batch_size}")
    print(f"  initial contour-follow fraction: {settings.contour_follow_fraction:.2f}")
    print(f"  late contour-follow fraction: {settings.late_contour_follow_fraction:.2f}")


def propose_monotone_hybrid_configs(
    *,
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
    rng: RNGGenerator,
) -> list[RMBConfig]:
    groups = grouped_by_n_qubits(data)
    candidates = []
    selected_points = []
    fitted_surface_found = False

    n_contour = int(round(settings.batch_size * settings.contour_follow_fraction))
    n_contour = max(0, min(settings.batch_size, n_contour))
    n_global = settings.batch_size - n_contour

    for n_qubits in settings.n_qubits_values:
        group = groups.get(n_qubits, {})
        try:
            surface = fit_monotone_fidelity_surface(group, settings)
        except (RuntimeError, ValueError):
            continue
        fitted_surface_found = True

        contour_candidates = contour_follow_candidates(
            group=group,
            surface=surface,
            settings=settings,
            rng=rng,
            n_qubits=n_qubits,
            n_candidates=max(n_contour, 1) * settings.contour_candidate_multiplier,
        )
        contour_candidates = [
            (
                score * hqc_efficiency_multiplier(config, data, settings),
                config,
                selected_point,
            )
            for score, config, selected_point in contour_candidates
        ]
        candidates.extend(contour_candidates)

        existing_points = scaled_data_points(group, settings)
        for _ in range(max(n_global, settings.batch_size - len(contour_candidates))):
            try:
                score, config, selected_point = optimize_hqc_aware_candidate_config(
                    surface=surface,
                    existing_points=existing_points,
                    selected_points=selected_points,
                    data=data,
                    template=RMBConfig.default()
                    .with_n_qubits(n_qubits)
                    .with_random_elimination(settings.random_elimination)
                    .with_scrambling_probability(settings.scrambling_probability),
                    settings=settings,
                    rng=rng,
                )
            except RuntimeError:
                continue
            candidates.append((score, config, selected_point))
            selected_points.append(selected_point)

    candidates.sort(key=lambda item: item[0], reverse=True)
    unique_configs = []
    seen = set()
    for _, config, _ in candidates:
        if config in seen:
            continue
        seen.add(config)
        unique_configs.append(config)
        if len(unique_configs) >= settings.batch_size:
            break

    if not fitted_surface_found and len(unique_configs) < settings.batch_size:
        random_configs = random_candidate_configs(
            settings=settings,
            rng=rng,
            n_candidates=max(
                settings.batch_size - len(unique_configs),
                settings.initial_candidate_multiplier * (settings.batch_size - len(unique_configs)),
            ),
        )
        random_configs = sort_configs_by_hqc_efficiency(
            random_configs,
            settings,
            shot_count=settings.max_adaptive_shots_per_config,
        )
        unique_configs.extend(random_configs[:settings.batch_size - len(unique_configs)])

    return unique_configs


def plot_monotone_fidelity_surface_contours(
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
    *,
    show: bool = True,
    axis_padding: float = 0.05,
) -> list:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    groups = sorted(grouped_by_n_qubits(data).items())
    if not groups:
        return []

    _, axes_arr = plt.subplots(1, len(groups), figsize=(5 * len(groups), 4), squeeze=False)
    axes = list(axes_arr[0])
    cmap = LinearSegmentedColormap.from_list(
        "darkred_to_lime",
        ["darkred", "red", "orange", "lime", "green"],
    )

    depth_span = settings.depth_bounds[1] - settings.depth_bounds[0]
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    x_pad = axis_padding * depth_span
    y_pad = axis_padding * ratio_span

    for ax, (n_qubits, group) in zip(axes, groups):
        depths = np.array([config.depth for config in group], dtype=float)
        ratios = np.array([config.min_two_qubit_gate_ratio for config in group], dtype=float)
        fidelities = np.array([fidelity_mean(estimator) for estimator in group.values()], dtype=float)
        variances = np.array([fidelity_variance(estimator) for estimator in group.values()], dtype=float)
        certainty = 1.0 - np.clip(variances / (1.0 / 12.0), 0.0, 1.0)
        marker_sizes = 30.0 + 160.0 * certainty

        scatter = ax.scatter(
            depths,
            ratios,
            c=fidelities,
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            s=marker_sizes,
            edgecolors="black",
            linewidths=0.4,
            zorder=3,
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Fidelity")
        ax.scatter([], [], s=30, facecolors="none", edgecolors="black", label="uncertain point")
        ax.scatter([], [], s=190, facecolors="none", edgecolors="black", label="certain point")

        try:
            surface = fit_monotone_fidelity_surface(group, settings)
        except (RuntimeError, ValueError):
            surface = None

        if surface is not None:
            depth_grid, ratio_grid, probabilities = surface.probability_grid(settings)
            if float(np.min(probabilities)) <= 0.5 <= float(np.max(probabilities)):
                ax.contour(
                    depth_grid,
                    ratio_grid,
                    probabilities,
                    levels=[0.5],
                    colors="black",
                    linewidths=2,
                    zorder=4,
                )
                ax.plot([], [], color="black", linewidth=2, label="monotone E[fidelity] = 0.5")

        ax.set_xlabel("# Gates")
        ax.set_ylabel("Two-qudit gate ratio")
        ax.set_title(f"# Qubits = {n_qubits}")
        ax.set_xlim(settings.depth_bounds[0] - x_pad, settings.depth_bounds[1] + x_pad)
        ax.set_ylim(
            max(0.0, settings.ratio_bounds[0] - y_pad),
            min(1.0, settings.ratio_bounds[1] + y_pad),
        )
        ax.legend(loc="best")

    if show:
        plt.show()

    return axes


def plot_monotone_fidelity_surface_with_confidence(
    data: RMBData,
    settings: MonotoneBoundaryExperimentConfig,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
    show: bool = True,
    mode: str = "heatmap",
) -> list:
    """
    Plot the monotone p=0.5 line with monotone bootstrap uncertainty.
    """
    import matplotlib.pyplot as plt

    axes = plot_monotone_fidelity_surface_contours(data, settings, show=False)
    groups = sorted(grouped_by_n_qubits(data).items())

    for ax, (_, group) in zip(axes, groups):
        contours = bootstrap_contours(
            group,
            settings,
            n_bootstrap=n_bootstrap,
            seed=seed,
        )
        if not contours:
            continue

        if mode == "lines":
            for contour in contours:
                ax.plot(
                    contour[:, 0],
                    contour[:, 1],
                    color="tab:blue",
                    alpha=min(0.2, 8.0 / max(1, n_bootstrap)),
                    linewidth=0.8,
                    zorder=2,
                )
            ax.plot([], [], color="tab:blue", alpha=0.6, linewidth=1.0,
                    label=f"monotone bootstrap contours ({len(contours)})")
            ax.legend(loc="best")
            continue

        depth_bins = np.linspace(
            settings.depth_bounds[0],
            settings.depth_bounds[1],
            settings.candidate_grid_size[0],
        )
        ratio_bins = np.linspace(
            settings.ratio_bounds[0],
            settings.ratio_bounds[1],
            settings.candidate_grid_size[1],
        )
        occupancy = np.zeros((len(ratio_bins) - 1, len(depth_bins) - 1), dtype=float)

        for contour in contours:
            depth_idx = np.searchsorted(depth_bins, contour[:, 0], side="right") - 1
            ratio_idx = np.searchsorted(ratio_bins, contour[:, 1], side="right") - 1
            valid = (
                (0 <= depth_idx)
                & (depth_idx < occupancy.shape[1])
                & (0 <= ratio_idx)
                & (ratio_idx < occupancy.shape[0])
            )
            cells = set(zip(ratio_idx[valid], depth_idx[valid]))
            for row, col in cells:
                occupancy[row, col] += 1.0

        occupancy /= len(contours)
        mesh = ax.pcolormesh(
            depth_bins,
            ratio_bins,
            occupancy,
            cmap="Blues",
            vmin=0.0,
            vmax=max(0.05, float(np.max(occupancy))),
            shading="auto",
            alpha=0.55,
            zorder=1,
        )
        cbar = plt.colorbar(mesh, ax=ax)
        cbar.set_label("Monotone bootstrap contour occupancy")
        ax.plot([], [], color="tab:blue", alpha=0.6, linewidth=6,
                label=f"monotone contour uncertainty ({len(contours)}/{n_bootstrap})")
        ax.legend(loc="best")

    if show:
        plt.show()

    return axes


def estimate_boundary(settings: MonotoneBoundaryExperimentConfig) -> RMB:
    rng = default_rng(settings.rng_seed)
    backend = make_backend()
    rmb = RMB.default(rng).with_backend(backend)
    data: RMBData = rmb._data

    remaining_measurements = settings.measurement_budget
    remaining_hqc = (
        float(settings.hqc_budget)
        if settings.hqc_budget is not None
        else float(settings.measurement_budget)
    )
    circuit_executions = 0
    max_execution_repeats = 0
    stop_reason = "budget not exhausted"
    initial_configs, initial_shots_per_config = budgeted_initial_design(settings, rng)
    initial_configs = sort_configs_by_hqc_efficiency(
        initial_configs,
        settings,
        shot_count=initial_shots_per_config,
    )

    for config in initial_configs:
        if remaining_measurements <= 0 or remaining_hqc <= 0.0:
            break
        requested_shots = min(initial_shots_per_config, remaining_measurements)
        requested_shots = affordable_shot_count(
            config=config,
            requested_shots=requested_shots,
            remaining_hqc=remaining_hqc,
            settings=settings,
        )
        if requested_shots <= 0:
            continue
        spent = spend_measurements(
            backend=backend,
            rng=rng,
            data=data,
            config=config,
            n_measurements=requested_shots,
        )
        remaining_measurements -= spent
        if spent > 0:
            circuit_executions += 1
            max_execution_repeats = max(max_execution_repeats, spent)
        if settings.hqc_budget is not None:
            remaining_hqc -= hqc_cost(config, spent, settings)
        else:
            remaining_hqc = float(remaining_measurements)

    if settings.verbose:
        budget_message = (
            f"{settings.hqc_budget - remaining_hqc:.3f} / {settings.hqc_budget:.3f} HQC"
            if settings.hqc_budget is not None
            else f"{total_measurements(data)} / {settings.measurement_budget} measurements"
        )
        print(
            f"Initial design complete: {budget_message}, "
            f"{total_measurements(data)} measurements across {len(data)} configs "
            f"({len(initial_configs)} initial configs, {initial_shots_per_config} shots each)."
        )
        print_fit_reports(data, settings)

    batch_index = 0
    while remaining_measurements > 0 and remaining_hqc > 0.0:
        batch_index += 1
        remaining_for_progress = (
            remaining_hqc
            if settings.hqc_budget is not None
            else float(remaining_measurements)
        )
        runtime_settings = adaptive_runtime_settings(data, settings, remaining_for_progress)
        batch = propose_monotone_hybrid_configs(data=data, settings=runtime_settings, rng=rng)
        if not batch:
            stop_reason = "no candidate batch proposed"
            break

        spent_this_batch = 0
        hqc_spent_this_batch = 0.0
        for config in batch:
            if remaining_measurements <= 0 or remaining_hqc <= 0.0:
                break
            n_measurements = adaptive_shot_count(
                config=config,
                data=data,
                settings=runtime_settings,
                remaining=remaining_measurements,
            )
            n_measurements = affordable_shot_count(
                config=config,
                requested_shots=n_measurements,
                remaining_hqc=remaining_hqc,
                settings=runtime_settings,
            )
            if n_measurements <= 0:
                continue
            spent = spend_measurements(
                backend=backend,
                rng=rng,
                data=data,
                config=config,
                n_measurements=n_measurements,
            )
            remaining_measurements -= spent
            spent_this_batch += spent
            if spent > 0:
                circuit_executions += 1
                max_execution_repeats = max(max_execution_repeats, spent)
            if settings.hqc_budget is not None:
                spent_hqc = hqc_cost(config, spent, runtime_settings)
                remaining_hqc -= spent_hqc
                hqc_spent_this_batch += spent_hqc
            else:
                remaining_hqc = float(remaining_measurements)

        if settings.verbose:
            if settings.hqc_budget is not None:
                budget_message = (
                    f"spent {hqc_spent_this_batch:.3f} HQC, "
                    f"total {settings.hqc_budget - remaining_hqc:.3f} / "
                    f"{settings.hqc_budget:.3f} HQC, "
                    f"measurements {total_measurements(data)} / {settings.measurement_budget}"
                )
            else:
                budget_message = (
                    f"spent {spent_this_batch}, "
                    f"total {total_measurements(data)} / {settings.measurement_budget}"
                )
            print(
                f"\nMonotone hybrid batch {batch_index}: {budget_message}, "
                f"configs {len(data)}, "
                f"batch_size {runtime_settings.batch_size}, "
                f"contour_fraction {runtime_settings.contour_follow_fraction:.2f}, "
                f"exploration_weight {runtime_settings.exploration_weight:.2f}."
            )
            print_fit_reports(data, settings)

        if spent_this_batch == 0:
            stop_reason = "no affordable measurements in proposed batch"
            break

    if remaining_measurements <= 0:
        stop_reason = "MEASUREMENT BUDGET HIT BEFORE HQC BUDGET"
    elif settings.hqc_budget is not None and remaining_hqc <= 0.0:
        stop_reason = "HQC budget hit"
    elif (
        settings.hqc_budget is not None
        and remaining_hqc <= settings.hqc_base_cost
        and stop_reason == "no affordable measurements in proposed batch"
    ):
        stop_reason = "HQC budget effectively hit: remaining HQC is below the base circuit cost"
    elif settings.hqc_budget is None and remaining_measurements <= 0:
        stop_reason = "MEASUREMENT BUDGET HIT"

    if settings.verbose:
        print_experiment_summary(
            data=data,
            settings=settings,
            stop_reason=stop_reason,
            batch_count=batch_index,
            circuit_executions=circuit_executions,
            max_execution_repeats=max_execution_repeats,
            remaining_measurements=remaining_measurements,
            remaining_hqc=remaining_hqc,
        )

    if settings.save_path is not None:
        rmb.save(settings.save_path)
        if settings.verbose:
            print(f"\nSaved boundary data to {settings.save_path}")

    return rmb


if __name__ == "__main__":
    settings = MonotoneBoundaryExperimentConfig(
        measurement_budget=1000,
        hqc_budget=1000.0,
        n_qubits_values=(50,),
        depth_bounds=(4, 100),
        ratio_bounds=(0.08, 0.6),
        initial_depths=4,
        initial_ratios=4,
        initial_shots_per_config=2,
        initial_budget_fraction=0.25,
        reserve_adaptive_measurements=60,
        initial_edge_margin=0.08,
        initial_candidate_multiplier=8,
        include_bracketing_points=True,
        batch_size=4,
        min_adaptive_shots_per_config=1,
        max_adaptive_shots_per_config=5,
        max_shots_per_config=10,
        candidate_grid_size=(80, 80),
        optimizer_maxiter=25,
        optimizer_popsize=8,
        optimizer_attempts_per_config=2,
        boundary_width=0.08,
        boundary_focus_power=2.0,
        max_boundary_probability_distance=0.25,
        exploration_weight=0.15,
        diversity_floor=0.5,
        shot_boundary_width=0.2,
        sparsity_radius=0.15,
        batch_diversity_radius=0.12,
        surface_smoothing=0.1,
        min_fit_points=8,
        monotone_l2=1e-3,
        contour_follow_fraction=0.50,
        contour_ready_probability_width=0.10,
        contour_min_anchors=4,
        contour_step_fraction=0.08,
        contour_projection_fraction=0.12,
        contour_gradient_fraction=0.01,
        contour_candidate_multiplier=5,
        adaptive_batch_size=True,
        min_batch_size=1,
        late_contour_follow_fraction=0.85,
        late_exploration_weight=0.05,
        adaptive_batch_late_budget_fraction=0.50,
        adaptive_batch_anchor_multiplier=2,
        hqc_cost_informed_acquisition=True,
        hqc_cost_power=1.0,
        save_path="viarregio4_boundary.json",
        verbose=True,
    )

    rmb = estimate_boundary(settings)
    print_fit_reports(rmb._data, settings)
    plot_monotone_fidelity_surface_with_confidence(
        rmb._data,
        settings,
        n_bootstrap=100,
        seed=settings.rng_seed,
    )
