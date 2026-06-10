from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from scipy.interpolate import RBFInterpolator
from scipy.optimize import differential_evolution
from scipy.special import expit, logit

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.core.noise.noise_model import GenericNoise


@dataclass(frozen=True)
class BoundaryExperimentConfig:
    """
    Settings for an adaptive fidelity=0.5 boundary-estimation experiment.

    `measurement_budget` is the total number of boolean fidelity outcomes to
    record across all configurations. The loop never asks estimators to run to
    convergence; it spends this global budget in small batches.
    """
    measurement_budget: int = 1_000
    n_qubits_values: tuple[int, ...] = (3,)
    depth_bounds: tuple[int, int] = (4, 80)
    ratio_bounds: tuple[float, float] = (0.0, 0.8)
    random_elimination: float = 0.1
    scrambling_probability: float = 0.0

    initial_depths: int = 4
    initial_ratios: int = 4
    initial_shots_per_config: int = 2
    initial_budget_fraction: float = 0.25
    reserve_adaptive_measurements: int = 60
    initial_edge_margin: float = 0.08
    initial_candidate_multiplier: int = 8
    include_bracketing_points: bool = True

    batch_size: int = 4
    min_adaptive_shots_per_config: int = 1
    max_adaptive_shots_per_config: int = 5
    max_shots_per_config: int = 10
    candidate_grid_size: tuple[int, int] = (80, 80)
    optimizer_maxiter: int = 25
    optimizer_popsize: int = 8
    optimizer_attempts_per_config: int = 2
    boundary_width: float = 0.08
    boundary_focus_power: float = 2.0
    max_boundary_probability_distance: float = 0.25
    exploration_weight: float = 0.15
    diversity_floor: float = 0.5
    shot_boundary_width: float = 0.2
    sparsity_radius: float = 0.15
    batch_diversity_radius: float = 0.12

    surface_smoothing: float = 0.05
    min_fit_points: int = 8
    rng_seed: int | None = 1234
    save_path: str | Path | None = "viarregio2_boundary.json"
    verbose: bool = True


@dataclass
class FlexibleFidelitySurface:
    """
    Nonlinear fitted surface for E[fidelity | depth, ratio].

    The interpolator is fit to logit-smoothed fidelity estimates. The
    fidelity=0.5 contour is therefore the zero level set of the interpolated
    logit surface, but callers can work directly with probabilities.
    """
    interpolator: RBFInterpolator
    depth_bounds: tuple[int, int]
    ratio_bounds: tuple[float, float]
    n_points: int

    def _scale(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        lower = np.array([self.depth_bounds[0], self.ratio_bounds[0]], dtype=float)
        upper = np.array([self.depth_bounds[1], self.ratio_bounds[1]], dtype=float)
        return (points - lower) / (upper - lower)

    def logit(self, points: np.ndarray) -> np.ndarray:
        return np.asarray(self.interpolator(self._scale(points))).reshape(-1)

    def probability(self, points: np.ndarray) -> np.ndarray:
        return expit(self.logit(points))

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
            "Fitted flexible fidelity surface:\n"
            f"  fit points: {self.n_points}\n"
            f"  predicted fidelity range on grid: "
            f"{min_probability:.4f} to {max_probability:.4f}\n"
            f"  contains fidelity=0.5 contour: {has_contour}"
        )


def make_backend() -> SympleqBackend:
    """
    Use the same SympleQBackend noise model as viarregio1.py.
    """
    noise_model = GenericNoise.from_paulis([0.000075, 0.000075, 0.000075])
    two_qubit_noise_model = GenericNoise.from_paulis([0.00039, 0.00039, 0.00039])
    return require_sympleq_backend(SympleqBackend(
        noise_model=noise_model,
        two_qubit_noise_model=two_qubit_noise_model,
    ))


def require_sympleq_backend(backend: object) -> SympleqBackend:
    """
    Hard guard for experiment scripts: never execute non-SympleQ backends here.
    """
    if type(backend) is not SympleqBackend:
        raise RuntimeError(
            "RMB Bayes experiments are locked to the local SympleqBackend; "
            f"refusing to run backend {type(backend).__module__}.{type(backend).__name__}."
        )
    return backend


def config_from_parameters(
    *,
    template: RMBConfig,
    depth: float,
    ratio: float,
    ratio_digits: int = 2,
) -> RMBConfig:
    ratio = round(float(np.clip(ratio, 0.0, 1.0)), ratio_digits)
    depth = max(1, int(round(depth)))
    n_slots = template.n_qubits * depth
    n_2qb_before_inverse = int(ratio * n_slots) // 2
    n_1qb_before_inverse = n_slots - 2 * n_2qb_before_inverse

    return (
        template
        .with_n_1qb_gates(2 * (n_1qb_before_inverse + template.n_qubits))
        .with_n_2qb_gates(2 * n_2qb_before_inverse)
    )


def config_depth(config: RMBConfig) -> int:
    """
    Return the legacy depth coordinate represented by a gate-count RMBConfig.
    """
    n_1qb_before_inverse = max(0, config.n_1qb_gates // 2 - config.n_qubits)
    n_2qb_before_inverse = max(0, config.n_2qb_gates // 2)
    n_slots = n_1qb_before_inverse + 2 * n_2qb_before_inverse
    return max(1, int(round(n_slots / config.n_qubits)))


def config_two_qubit_gate_ratio(config: RMBConfig) -> float:
    """
    Return the legacy two-qubit slot ratio represented by a gate-count config.
    """
    n_1qb_before_inverse = max(0, config.n_1qb_gates // 2 - config.n_qubits)
    n_2qb_before_inverse = max(0, config.n_2qb_gates // 2)
    n_slots = n_1qb_before_inverse + 2 * n_2qb_before_inverse
    if n_slots <= 0:
        return 0.0
    return float(2 * n_2qb_before_inverse / n_slots)


def _legacy_with_depth(config: RMBConfig, depth: int) -> RMBConfig:
    return config_from_parameters(
        template=config,
        depth=depth,
        ratio=config_two_qubit_gate_ratio(config),
    )


def _legacy_with_two_qubit_gate_ratio(
    config: RMBConfig,
    min_ratio: float,
    max_ratio: float | None = None,
) -> RMBConfig:
    if max_ratio is None:
        ratio = min_ratio
    else:
        ratio = 0.5 * (min_ratio + max_ratio)
    return config_from_parameters(
        template=config,
        depth=config_depth(config),
        ratio=ratio,
    )


if not hasattr(RMBConfig, "depth"):
    RMBConfig.depth = property(config_depth)  # type: ignore[attr-defined]
if not hasattr(RMBConfig, "min_two_qubit_gate_ratio"):
    RMBConfig.min_two_qubit_gate_ratio = property(config_two_qubit_gate_ratio)  # type: ignore[attr-defined]
if not hasattr(RMBConfig, "max_two_qubit_gate_ratio"):
    RMBConfig.max_two_qubit_gate_ratio = property(config_two_qubit_gate_ratio)  # type: ignore[attr-defined]
if not hasattr(RMBConfig, "with_depth"):
    RMBConfig.with_depth = _legacy_with_depth  # type: ignore[attr-defined]
if not hasattr(RMBConfig, "with_two_qubit_gate_ratio"):
    RMBConfig.with_two_qubit_gate_ratio = _legacy_with_two_qubit_gate_ratio  # type: ignore[attr-defined]


def template_config(settings: BoundaryExperimentConfig, n_qubits: int) -> RMBConfig:
    return (
        RMBConfig.default()
        .with_n_qubits(n_qubits)
        .with_random_elimination(settings.random_elimination)
        .with_scrambling_probability(settings.scrambling_probability)
    )


def interior_bounds(
    bounds: tuple[float, float] | tuple[int, int],
    margin_fraction: float,
) -> tuple[float, float]:
    lower, upper = float(bounds[0]), float(bounds[1])
    margin = margin_fraction * (upper - lower)
    if 2.0 * margin >= upper - lower:
        return lower, upper
    return lower + margin, upper - margin


def latin_hypercube_points(
    n_points: int,
    n_dimensions: int,
    rng: RNGGenerator,
) -> np.ndarray:
    points = np.empty((n_points, n_dimensions), dtype=float)
    for dim in range(n_dimensions):
        points[:, dim] = (np.arange(n_points) + rng.random(n_points)) / n_points
        rng.shuffle(points[:, dim])
    return points


def initial_design(settings: BoundaryExperimentConfig, rng: RNGGenerator) -> list[RMBConfig]:
    n_candidates = max(
        settings.initial_depths * settings.initial_ratios,
        settings.min_fit_points * settings.initial_candidate_multiplier,
    )
    depth_bounds = interior_bounds(settings.depth_bounds, settings.initial_edge_margin)
    ratio_bounds = interior_bounds(settings.ratio_bounds, settings.initial_edge_margin)
    configs = []

    for n_qubits in settings.n_qubits_values:
        template = template_config(settings, n_qubits)
        points = latin_hypercube_points(n_candidates, 2, rng)
        for depth_unit, ratio_unit in points:
            depth = depth_bounds[0] + depth_unit * (depth_bounds[1] - depth_bounds[0])
            ratio = ratio_bounds[0] + ratio_unit * (ratio_bounds[1] - ratio_bounds[0])
            configs.append(config_from_parameters(
                template=template,
                depth=depth,
                ratio=ratio,
            ))

    return sorted(
        set(configs),
        key=lambda c: (c.n_qubits, c.depth, c.min_two_qubit_gate_ratio),
    )


def bracketing_design(settings: BoundaryExperimentConfig) -> list[RMBConfig]:
    """
    Include a few deliberate probes that help the fitted surface bracket 0.5.

    These are not meant to estimate the whole line. They reduce the chance that
    a low-budget run only observes one side of the transition and therefore
    cannot produce a p=0.5 contour at all.
    """
    if not settings.include_bracketing_points:
        return []

    d_min, d_max = settings.depth_bounds
    r_min, r_max = settings.ratio_bounds
    d_inner_min, d_inner_max = interior_bounds(settings.depth_bounds, settings.initial_edge_margin)
    r_inner_min, r_inner_max = interior_bounds(settings.ratio_bounds, settings.initial_edge_margin)
    d_mid = 0.5 * (d_min + d_max)
    r_mid = 0.5 * (r_min + r_max)

    parameter_points = [
        (d_mid, r_mid),                    # central transition check
        (d_inner_max, r_mid),              # depth stress, away from corner
        (d_mid, r_inner_max),              # ratio stress, away from corner
        (d_inner_max, r_inner_min),        # depth-dominant hard point
        (d_inner_min, r_inner_max),        # ratio-dominant hard point
    ]

    configs = []
    for n_qubits in settings.n_qubits_values:
        template = template_config(settings, n_qubits)
        for depth, ratio in parameter_points:
            configs.append(config_from_parameters(
                template=template,
                depth=depth,
                ratio=ratio,
            ))

    return sorted(
        set(configs),
        key=lambda c: (c.n_qubits, c.depth, c.min_two_qubit_gate_ratio),
    )


def config_design_point(config: RMBConfig, settings: BoundaryExperimentConfig) -> np.ndarray:
    lower = np.array([settings.depth_bounds[0], settings.ratio_bounds[0]], dtype=float)
    upper = np.array([settings.depth_bounds[1], settings.ratio_bounds[1]], dtype=float)
    point = np.array([float(config.depth), float(config.min_two_qubit_gate_ratio)], dtype=float)
    return (point - lower) / (upper - lower)


def select_space_filling_configs(
    configs: list[RMBConfig],
    n_configs: int,
    settings: BoundaryExperimentConfig,
) -> list[RMBConfig]:
    """
    Deterministically thin candidate configs with farthest-point sampling.
    """
    if n_configs >= len(configs):
        return configs

    points = np.asarray([config_design_point(config, settings) for config in configs])
    center = np.array([0.5, 0.5], dtype=float)
    first_idx = int(np.argmin(np.linalg.norm(points - center, axis=1)))

    selected_indices = [first_idx]
    remaining_indices = set(range(len(configs))) - {first_idx}

    while len(selected_indices) < n_configs and remaining_indices:
        selected_points = points[selected_indices]
        best_idx = max(
            remaining_indices,
            key=lambda idx: float(np.min(np.linalg.norm(selected_points - points[idx], axis=1))),
        )
        selected_indices.append(best_idx)
        remaining_indices.remove(best_idx)

    selected = [configs[idx] for idx in selected_indices]
    selected.sort(key=lambda c: (c.n_qubits, c.depth, c.min_two_qubit_gate_ratio))
    return selected


def budgeted_initial_design(
    settings: BoundaryExperimentConfig,
    rng: RNGGenerator,
) -> tuple[list[RMBConfig], int]:
    """
    Choose an initial design that leaves measurements for adaptive updates.
    """
    configs = initial_design(settings, rng)
    if not configs or settings.measurement_budget <= 0:
        return [], 0

    reserved = max(0, settings.reserve_adaptive_measurements)
    adaptive_respecting_budget = max(1, settings.measurement_budget - reserved)
    fractional_budget = max(1, int(settings.initial_budget_fraction * settings.measurement_budget))
    initial_budget = min(adaptive_respecting_budget, fractional_budget)

    shots_per_config = max(1, min(settings.initial_shots_per_config, initial_budget))
    max_configs = max(1, initial_budget // shots_per_config)

    bracket_configs = bracketing_design(settings)
    if max_configs <= len(bracket_configs):
        configs = select_space_filling_configs(bracket_configs, max_configs, settings)
    elif max_configs < len(configs) + len(bracket_configs):
        bracket_set = set(bracket_configs)
        remaining_candidates = [config for config in configs if config not in bracket_set]
        n_fill = max_configs - len(bracket_configs)
        fill_configs = select_space_filling_configs(remaining_candidates, n_fill, settings)
        configs = sorted(
            set(bracket_configs + fill_configs),
            key=lambda c: (c.n_qubits, c.depth, c.min_two_qubit_gate_ratio),
        )
    else:
        configs = sorted(
            set(bracket_configs + configs),
            key=lambda c: (c.n_qubits, c.depth, c.min_two_qubit_gate_ratio),
        )

    return configs, shots_per_config


def estimator_for(data: RMBData, config: RMBConfig) -> BayesianEstimator:
    return data.setdefault(config, BayesianEstimator(threshold=0.0, min_runs=0))


def spend_measurements(
    *,
    backend: SympleqBackend,
    rng: RNGGenerator,
    data: RMBData,
    config: RMBConfig,
    n_measurements: int,
) -> int:
    """
    Record up to `n_measurements` boolean outcomes for one config.
    """
    backend = require_sympleq_backend(backend)
    estimator = estimator_for(data, config)
    spent = 0

    while spent < n_measurements:
        outcomes = backend.fidelity_estimation(config, rng)
        if not isinstance(outcomes, list):
            outcomes = [outcomes]

        for outcome in outcomes:
            if spent >= n_measurements:
                break
            record_estimator = estimator
            fidelity = outcome
            if isinstance(outcome, tuple) and len(outcome) == 2:
                outcome_config, fidelity = outcome
                record_estimator = estimator_for(data, outcome_config)
            record_estimator.record(bool(fidelity))
            spent += 1

    return spent


def total_measurements(data: RMBData) -> int:
    return sum(estimator.num_runs() for estimator in data.values())


def fidelity_mean(estimator: BayesianEstimator) -> float:
    """
    Beta(1, 1)-smoothed posterior mean for the boolean fidelity outcome.
    """
    counts = estimator.counts()
    true_count = counts.get(True, 0)
    false_count = counts.get(False, 0)
    return (true_count + 1.0) / (true_count + false_count + 2.0)


def fidelity_variance(estimator: BayesianEstimator) -> float:
    """
    Beta(1, 1)-posterior variance for the boolean fidelity probability.
    """
    counts = estimator.counts()
    true_count = counts.get(True, 0)
    false_count = counts.get(False, 0)
    alpha = true_count + 1.0
    beta = false_count + 1.0
    total = alpha + beta
    return alpha * beta / (total * total * (total + 1.0))


def grouped_by_n_qubits(data: RMBData) -> dict[int, RMBData]:
    groups: dict[int, RMBData] = {}
    for config, estimator in data.items():
        if estimator.num_runs() == 0:
            continue
        groups.setdefault(config.n_qubits, {})[config] = estimator
    return groups


def fit_fidelity_surface(
    data: RMBData,
    settings: BoundaryExperimentConfig,
) -> FlexibleFidelitySurface:
    """
    Fit a nonlinear smooth fidelity surface for one fixed n_qubits group.
    """
    points = []
    values = []

    for config, estimator in data.items():
        if estimator.num_runs() == 0:
            continue
        points.append([float(config.depth), float(config.min_two_qubit_gate_ratio)])
        values.append(fidelity_mean(estimator))

    if len(points) < settings.min_fit_points:
        raise ValueError(
            f"Need at least {settings.min_fit_points} data points to fit a surface, "
            f"got {len(points)}."
        )

    points_array = np.asarray(points, dtype=float)
    values_array = np.clip(np.asarray(values, dtype=float), 1e-4, 1.0 - 1e-4)
    lower = np.array([settings.depth_bounds[0], settings.ratio_bounds[0]], dtype=float)
    upper = np.array([settings.depth_bounds[1], settings.ratio_bounds[1]], dtype=float)
    scaled_points = (points_array - lower) / (upper - lower)

    interpolator = RBFInterpolator(
        scaled_points,
        logit(values_array),
        kernel="thin_plate_spline",
        degree=1,
        smoothing=settings.surface_smoothing,
    )
    return FlexibleFidelitySurface(
        interpolator=interpolator,
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=len(points),
    )


def random_candidate_configs(
    *,
    settings: BoundaryExperimentConfig,
    rng: RNGGenerator,
    n_candidates: int,
) -> list[RMBConfig]:
    configs = []
    for _ in range(n_candidates):
        n_qubits = int(rng.choice(settings.n_qubits_values))
        template = template_config(settings, n_qubits)
        depth = rng.uniform(settings.depth_bounds[0], settings.depth_bounds[1])
        ratio = rng.uniform(settings.ratio_bounds[0], settings.ratio_bounds[1])
        configs.append(config_from_parameters(
            template=template,
            depth=depth,
            ratio=ratio,
        ))
    return configs


def scaled_config_point(
    *,
    depth: float,
    ratio: float,
    settings: BoundaryExperimentConfig,
) -> np.ndarray:
    lower = np.array([settings.depth_bounds[0], settings.ratio_bounds[0]], dtype=float)
    upper = np.array([settings.depth_bounds[1], settings.ratio_bounds[1]], dtype=float)
    return (np.array([depth, ratio], dtype=float) - lower) / (upper - lower)


def scaled_data_points(data: RMBData, settings: BoundaryExperimentConfig) -> np.ndarray:
    points = [
        scaled_config_point(
            depth=float(config.depth),
            ratio=float(config.min_two_qubit_gate_ratio),
            settings=settings,
        )
        for config, estimator in data.items()
        if estimator.num_runs() > 0
    ]
    if not points:
        return np.empty((0, 2), dtype=float)
    return np.asarray(points, dtype=float)


def nearest_scaled_distance(point: np.ndarray, points: np.ndarray) -> float:
    if len(points) == 0:
        return 1.0
    distances = np.linalg.norm(points - point, axis=1)
    return float(np.min(distances))


def acquisition_score(
    *,
    theta: np.ndarray,
    surface: FlexibleFidelitySurface,
    existing_points: np.ndarray,
    selected_points: list[np.ndarray],
    data: RMBData,
    template: RMBConfig,
    settings: BoundaryExperimentConfig,
) -> float:
    """
    Boundary-focused acquisition score for continuous candidate parameters.

    The score is dominated by closeness to the fitted p=0.5 contour. Sparsity
    and batch diversity are weak tie-breakers so low-budget runs do not spend
    many shots far from the line.
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
    return float(boundary_relevance * sparsity_multiplier * diversity_multiplier * undersampled)


def optimize_candidate_config(
    *,
    surface: FlexibleFidelitySurface,
    existing_points: np.ndarray,
    selected_points: list[np.ndarray],
    data: RMBData,
    template: RMBConfig,
    settings: BoundaryExperimentConfig,
    rng: RNGGenerator,
) -> tuple[float, RMBConfig, np.ndarray]:
    bounds = [
        settings.depth_bounds,
        settings.ratio_bounds,
    ]

    def objective(theta: np.ndarray) -> float:
        score = acquisition_score(
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


def adaptive_shot_count(
    *,
    config: RMBConfig,
    data: RMBData,
    settings: BoundaryExperimentConfig,
    remaining: int,
) -> int:
    """
    Allocate more shots to configs that are both relevant and uncertain.
    """
    current_runs = data.get(config).num_runs() if config in data else 0
    available_for_config = settings.max_shots_per_config - current_runs
    if available_for_config <= 0 or remaining <= 0:
        return 0

    if config in data:
        estimator = data[config]
        mean = fidelity_mean(estimator)
        variance = fidelity_variance(estimator)
    else:
        mean = 0.5
        variance = 1.0 / 12.0

    boundary_relevance = np.exp(-((mean - 0.5) / settings.shot_boundary_width) ** 2)
    uncertainty = min(1.0, variance / (1.0 / 12.0))
    need = float(boundary_relevance * uncertainty)

    shot_span = settings.max_adaptive_shots_per_config - settings.min_adaptive_shots_per_config
    requested = settings.min_adaptive_shots_per_config + int(np.ceil(shot_span * need))
    return max(0, min(requested, available_for_config, remaining))


def propose_adaptive_configs(
    *,
    data: RMBData,
    settings: BoundaryExperimentConfig,
    rng: RNGGenerator,
) -> list[RMBConfig]:
    """
    Propose configs by optimizing a boundary-focused acquisition rule.

    The optimized score is high where the fitted expected fidelity is close
    to 0.5, existing data are sparse, the rounded config is undersampled, and
    the point is separated from candidates already selected for this batch.
    When no surface fit is available yet, random candidates fill the batch.
    """
    candidates: list[tuple[float, RMBConfig, np.ndarray]] = []
    groups = grouped_by_n_qubits(data)
    selected_points: list[np.ndarray] = []
    fitted_surface_found = False

    for n_qubits in settings.n_qubits_values:
        group = groups.get(n_qubits, {})
        template = template_config(settings, n_qubits)

        try:
            surface = fit_fidelity_surface(group, settings)
        except (RuntimeError, ValueError):
            continue
        fitted_surface_found = True

        existing_points = scaled_data_points(group, settings)
        for _ in range(settings.batch_size):
            try:
                score, config, selected_point = optimize_candidate_config(
                    surface=surface,
                    existing_points=existing_points,
                    selected_points=selected_points,
                    data=data,
                    template=template,
                    settings=settings,
                    rng=rng,
                )
            except RuntimeError:
                continue

            jitter = 1e-9 * float(rng.random())
            candidates.append((score + jitter, config, selected_point))
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
        unique_configs.extend(random_candidate_configs(
            settings=settings,
            rng=rng,
            n_candidates=settings.batch_size - len(unique_configs),
        ))

    return unique_configs


def print_fit_reports(data: RMBData, settings: BoundaryExperimentConfig) -> None:
    for n_qubits, group in sorted(grouped_by_n_qubits(data).items()):
        try:
            surface = fit_fidelity_surface(group, settings)
        except (RuntimeError, ValueError):
            print(f"n_qubits={n_qubits}: not enough data for a surface fit yet.")
            continue

        print(f"\nn_qubits={n_qubits}")
        print(surface.report(settings))


def plot_fidelity_surface_contours(
    data: RMBData,
    settings: BoundaryExperimentConfig,
    *,
    show: bool = True,
    axis_padding: float = 0.05,
) -> list:
    """
    Plot data and the nonlinear fitted E[fidelity]=0.5 contour per n_qubits.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    groups = sorted(grouped_by_n_qubits(data).items())
    if not groups:
        return []

    _, axes_arr = plt.subplots(
        1,
        len(groups),
        figsize=(5 * len(groups), 4),
        squeeze=False,
    )
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
        ratios = np.array(
            [config.min_two_qubit_gate_ratio for config in group],
            dtype=float,
        )
        fidelities = np.array(
            [fidelity_mean(estimator) for estimator in group.values()],
            dtype=float,
        )

        scatter = ax.scatter(
            depths,
            ratios,
            c=fidelities,
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            s=90,
            edgecolors="black",
            linewidths=0.4,
            zorder=3,
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Fidelity")

        try:
            surface = fit_fidelity_surface(group, settings)
        except (RuntimeError, ValueError):
            surface = None

        if surface is not None:
            depth_grid, ratio_grid, probabilities = surface.probability_grid(settings)
            min_probability = float(np.min(probabilities))
            max_probability = float(np.max(probabilities))

            if min_probability <= 0.5 <= max_probability:
                ax.contour(
                    depth_grid,
                    ratio_grid,
                    probabilities,
                    levels=[0.5],
                    colors="black",
                    linewidths=2,
                    zorder=4,
                )
                ax.plot([], [], color="black", linewidth=2, label="E[fidelity] = 0.5")
                ax.legend(loc="best")
            else:
                closeness = np.abs(probabilities - 0.5)
                closest = np.argsort(closeness.ravel())[:settings.candidate_grid_size[0]]
                ax.scatter(
                    depth_grid.ravel()[closest],
                    ratio_grid.ravel()[closest],
                    color="black",
                    s=8,
                    marker="x",
                    label="closest fitted points to 0.5",
                    zorder=4,
                )
                ax.legend(loc="best")

        ax.set_xlabel("# Gates")
        ax.set_ylabel("Two-qudit gate ratio")
        ax.set_title(f"# Qubits = {n_qubits}")
        ax.set_xlim(settings.depth_bounds[0] - x_pad, settings.depth_bounds[1] + x_pad)
        ax.set_ylim(
            max(0.0, settings.ratio_bounds[0] - y_pad),
            min(1.0, settings.ratio_bounds[1] + y_pad),
        )

    if show:
        plt.show()

    return axes


def estimate_boundary(settings: BoundaryExperimentConfig) -> RMB:
    rng = default_rng(settings.rng_seed)
    backend = make_backend()
    rmb = RMB.default(rng).with_backend(backend)
    data: RMBData = rmb._data

    remaining = settings.measurement_budget
    initial_configs, initial_shots_per_config = budgeted_initial_design(settings, rng)

    for config in initial_configs:
        if remaining <= 0:
            break
        n_measurements = min(initial_shots_per_config, remaining)
        spent = spend_measurements(
            backend=backend,
            rng=rng,
            data=data,
            config=config,
            n_measurements=n_measurements,
        )
        remaining -= spent

    if settings.verbose:
        print(
            f"Initial design complete: {total_measurements(data)} / "
            f"{settings.measurement_budget} measurements across {len(data)} configs "
            f"({len(initial_configs)} initial configs, {initial_shots_per_config} shots each)."
        )
        print_fit_reports(data, settings)

    batch_index = 0
    while remaining > 0:
        batch_index += 1
        batch = propose_adaptive_configs(data=data, settings=settings, rng=rng)
        if not batch:
            break

        spent_this_batch = 0
        for config in batch:
            if remaining <= 0:
                break

            n_measurements = adaptive_shot_count(
                config=config,
                data=data,
                settings=settings,
                remaining=remaining,
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
            remaining -= spent
            spent_this_batch += spent

        if settings.verbose:
            print(
                f"\nBatch {batch_index}: spent {spent_this_batch}, "
                f"total {total_measurements(data)} / {settings.measurement_budget}, "
                f"configs {len(data)}."
            )
            print_fit_reports(data, settings)

        if spent_this_batch == 0:
            break

    if settings.save_path is not None:
        rmb.save(settings.save_path)
        if settings.verbose:
            print(f"\nSaved boundary data to {settings.save_path}")

    return rmb


if __name__ == "__main__":
    settings = BoundaryExperimentConfig(
        measurement_budget=100,
        n_qubits_values=(10,),
        depth_bounds=(4, 80),
        ratio_bounds=(0.0, 0.8),
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
        save_path="viarregio2_boundary.json",
        verbose=True,
    )

    rmb = estimate_boundary(settings)
    print_fit_reports(rmb._data, settings)
    plot_fidelity_surface_contours(rmb._data, settings)
