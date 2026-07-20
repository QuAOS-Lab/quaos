from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from scipy.optimize import minimize_scalar

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData

from viarregio2 import (
    BoundaryExperimentConfig,
    FlexibleFidelitySurface,
    adaptive_shot_count,
    budgeted_initial_design,
    config_from_parameters,
    fit_fidelity_surface,
    grouped_by_n_qubits,
    make_backend,
    optimize_candidate_config,
    plot_fidelity_surface_contours,
    print_fit_reports,
    random_candidate_configs,
    scaled_config_point,
    scaled_data_points,
    spend_measurements,
    template_config,
    total_measurements,
)


def contour_points_from_surface(
    surface: FlexibleFidelitySurface,
    settings: BoundaryExperimentConfig,
    level: float = 0.5,
) -> np.ndarray:
    """
    Approximate contour points by linearly interpolating grid-edge crossings.
    """
    depth_grid, ratio_grid, probabilities = surface.probability_grid(settings)
    points = []

    for row in range(probabilities.shape[0]):
        for col in range(probabilities.shape[1] - 1):
            p0 = probabilities[row, col] - level
            p1 = probabilities[row, col + 1] - level
            if p0 == 0.0:
                points.append([depth_grid[row, col], ratio_grid[row, col]])
            if p0 * p1 < 0.0:
                t = abs(p0) / (abs(p0) + abs(p1))
                points.append([
                    (1.0 - t) * depth_grid[row, col] + t * depth_grid[row, col + 1],
                    ratio_grid[row, col],
                ])

    for row in range(probabilities.shape[0] - 1):
        for col in range(probabilities.shape[1]):
            p0 = probabilities[row, col] - level
            p1 = probabilities[row + 1, col] - level
            if p0 == 0.0:
                points.append([depth_grid[row, col], ratio_grid[row, col]])
            if p0 * p1 < 0.0:
                t = abs(p0) / (abs(p0) + abs(p1))
                points.append([
                    depth_grid[row, col],
                    (1.0 - t) * ratio_grid[row, col] + t * ratio_grid[row + 1, col],
                ])

    if not points:
        return np.empty((0, 2), dtype=float)
    return np.unique(np.asarray(points, dtype=float), axis=0)


def fit_surface_from_values(
    configs: list[RMBConfig],
    fidelities: np.ndarray,
    settings: BoundaryExperimentConfig,
) -> FlexibleFidelitySurface:
    from scipy.interpolate import RBFInterpolator
    from scipy.special import logit

    points = np.array(
        [[float(config.depth), float(config.min_two_qubit_gate_ratio)] for config in configs],
        dtype=float,
    )
    values = np.clip(np.asarray(fidelities, dtype=float), 1e-4, 1.0 - 1e-4)
    lower = np.array([settings.depth_bounds[0], settings.ratio_bounds[0]], dtype=float)
    upper = np.array([settings.depth_bounds[1], settings.ratio_bounds[1]], dtype=float)
    scaled_points = (points - lower) / (upper - lower)

    interpolator = RBFInterpolator(
        scaled_points,
        logit(values),
        kernel="thin_plate_spline",
        degree=1,
        smoothing=settings.surface_smoothing,
    )
    return FlexibleFidelitySurface(
        interpolator=interpolator,
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=len(configs),
    )


def bootstrap_contours(
    data: RMBData,
    settings: BoundaryExperimentConfig,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
) -> list[np.ndarray]:
    """
    Draw posterior fidelity surfaces and return their p=0.5 contours.

    For each measured config, sample a fidelity probability from its
    Beta(1 + successes, 1 + failures) posterior, refit the surface, and
    extract a contour.
    """
    rng = default_rng(seed)
    configs = [config for config, estimator in data.items() if estimator.num_runs() > 0]
    if len(configs) < settings.min_fit_points:
        return []

    alpha = []
    beta = []
    for config in configs:
        counts = data[config].counts()
        alpha.append(counts.get(True, 0) + 1.0)
        beta.append(counts.get(False, 0) + 1.0)

    alpha_array = np.asarray(alpha, dtype=float)
    beta_array = np.asarray(beta, dtype=float)
    contours = []

    for _ in range(n_bootstrap):
        sampled_fidelities = rng.beta(alpha_array, beta_array)
        try:
            surface = fit_surface_from_values(configs, sampled_fidelities, settings)
        except (RuntimeError, ValueError):
            continue
        contour = contour_points_from_surface(surface, settings)
        if len(contour) > 0:
            contours.append(contour)

    return contours


def plot_fidelity_surface_with_confidence(
    data: RMBData,
    settings: BoundaryExperimentConfig,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
    show: bool = True,
    mode: str = "heatmap",
) -> list:
    """
    Plot the fitted p=0.5 line with bootstrap uncertainty.

    `mode="heatmap"` shows where bootstrap contours pass most often, which is
    clearer than raw contour lines when uncertainty is large. `mode="lines"`
    overlays the raw bootstrap contour samples.
    """
    import matplotlib.pyplot as plt

    axes = plot_fidelity_surface_contours(data, settings, show=False)
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
                    label=f"bootstrap contours ({len(contours)})")
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
        cbar.set_label("Bootstrap contour occupancy")
        ax.plot([], [], color="tab:blue", alpha=0.6, linewidth=6,
                label=f"contour uncertainty ({len(contours)}/{n_bootstrap})")
        ax.legend(loc="best")

    if show:
        plt.show()

    return axes


@dataclass(frozen=True)
class HybridBoundaryExperimentConfig(BoundaryExperimentConfig):
    """
    Boundary-estimation settings with a hybrid update rule.

    The update rule starts with global boundary acquisition. Once the fitted
    surface has enough near-contour anchors, most proposals follow the current
    continuous p=0.5 contour by stepping along local tangent directions and
    projecting back to the contour.
    """
    contour_follow_fraction: float = 0.50
    contour_ready_probability_width: float = 0.10
    contour_min_anchors: int = 4
    contour_step_fraction: float = 0.08
    contour_projection_fraction: float = 0.12
    contour_gradient_fraction: float = 0.01
    contour_candidate_multiplier: int = 5
    save_path: str | Path | None = "viarregio3_boundary.json"


def finite_difference_gradient(
    surface: FlexibleFidelitySurface,
    point: np.ndarray,
    settings: HybridBoundaryExperimentConfig,
) -> np.ndarray:
    """
    Estimate grad p(depth, ratio) at one point.
    """
    depth_span = settings.depth_bounds[1] - settings.depth_bounds[0]
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    h_depth = settings.contour_gradient_fraction * depth_span
    h_ratio = settings.contour_gradient_fraction * ratio_span

    depth, ratio = float(point[0]), float(point[1])
    p_depth_plus = surface.probability(np.array([[depth + h_depth, ratio]]))[0]
    p_depth_minus = surface.probability(np.array([[depth - h_depth, ratio]]))[0]
    p_ratio_plus = surface.probability(np.array([[depth, ratio + h_ratio]]))[0]
    p_ratio_minus = surface.probability(np.array([[depth, ratio - h_ratio]]))[0]

    return np.array([
        (p_depth_plus - p_depth_minus) / (2.0 * h_depth),
        (p_ratio_plus - p_ratio_minus) / (2.0 * h_ratio),
    ])


def project_to_contour(
    surface: FlexibleFidelitySurface,
    point: np.ndarray,
    normal: np.ndarray,
    settings: HybridBoundaryExperimentConfig,
) -> np.ndarray:
    """
    Move along the local normal to get close to the fitted p=0.5 contour.
    """
    depth_span = settings.depth_bounds[1] - settings.depth_bounds[0]
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    scale = np.array([depth_span, ratio_span], dtype=float)
    max_step = settings.contour_projection_fraction

    normal_scaled = normal * scale
    norm = np.linalg.norm(normal_scaled)
    if norm < 1e-12:
        return point
    direction = normal_scaled / norm

    def objective(t: float) -> float:
        candidate = point + t * direction * scale
        candidate[0] = np.clip(candidate[0], settings.depth_bounds[0], settings.depth_bounds[1])
        candidate[1] = np.clip(candidate[1], settings.ratio_bounds[0], settings.ratio_bounds[1])
        probability = float(surface.probability(candidate.reshape(1, 2))[0])
        return abs(probability - 0.5)

    result = minimize_scalar(objective, bounds=(-max_step, max_step), method="bounded")
    projected = point + float(result.x) * direction * scale
    projected[0] = np.clip(projected[0], settings.depth_bounds[0], settings.depth_bounds[1])
    projected[1] = np.clip(projected[1], settings.ratio_bounds[0], settings.ratio_bounds[1])
    return projected


def near_contour_anchors(
    data: RMBData,
    surface: FlexibleFidelitySurface,
    settings: HybridBoundaryExperimentConfig,
) -> list[RMBConfig]:
    anchors = []
    for config, estimator in data.items():
        if estimator.num_runs() == 0:
            continue
        point = np.array([[float(config.depth), float(config.min_two_qubit_gate_ratio)]])
        probability = float(surface.probability(point)[0])
        if abs(probability - 0.5) <= settings.contour_ready_probability_width:
            anchors.append(config)
    return anchors


def contour_follow_candidates(
    *,
    group: RMBData,
    surface: FlexibleFidelitySurface,
    settings: HybridBoundaryExperimentConfig,
    rng: RNGGenerator,
    n_qubits: int,
    n_candidates: int,
) -> list[tuple[float, RMBConfig, np.ndarray]]:
    """
    Follow the current contour from near-boundary anchors.
    """
    anchors = near_contour_anchors(group, surface, settings)
    if len(anchors) < settings.contour_min_anchors:
        return []

    template = template_config(settings, n_qubits)
    existing_points = scaled_data_points(group, settings)
    selected_points: list[np.ndarray] = []
    candidates: list[tuple[float, RMBConfig, np.ndarray]] = []

    depth_span = settings.depth_bounds[1] - settings.depth_bounds[0]
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    tangent_step = settings.contour_step_fraction * np.array([depth_span, ratio_span])

    anchor_order = list(rng.permutation(len(anchors)))
    for anchor_index in anchor_order:
        anchor = anchors[int(anchor_index)]
        anchor_point = np.array([float(anchor.depth), float(anchor.min_two_qubit_gate_ratio)])
        gradient = finite_difference_gradient(surface, anchor_point, settings)
        if np.linalg.norm(gradient) < 1e-12:
            continue

        tangent = np.array([gradient[1], -gradient[0]], dtype=float)
        tangent_scaled = tangent * np.array([depth_span, ratio_span])
        tangent_norm = np.linalg.norm(tangent_scaled)
        if tangent_norm < 1e-12:
            continue
        tangent_unit = tangent_scaled / tangent_norm

        for direction in (-1.0, 1.0):
            trial = anchor_point + direction * tangent_unit * tangent_step
            trial[0] = np.clip(trial[0], settings.depth_bounds[0], settings.depth_bounds[1])
            trial[1] = np.clip(trial[1], settings.ratio_bounds[0], settings.ratio_bounds[1])
            projected = project_to_contour(surface, trial, gradient, settings)

            config = config_from_parameters(
                template=template,
                depth=float(projected[0]),
                ratio=float(projected[1]),
            )
            if config in group and group[config].num_runs() >= settings.max_shots_per_config:
                continue

            selected_point = scaled_config_point(
                depth=float(config.depth),
                ratio=float(config.min_two_qubit_gate_ratio),
                settings=settings,
            )
            probability = float(surface.probability(
                np.array([[float(config.depth), float(config.min_two_qubit_gate_ratio)]])
            )[0])
            contour_score = np.exp(-((abs(probability - 0.5) / settings.boundary_width) ** 2))
            if len(existing_points) > 0:
                distance = float(np.min(np.linalg.norm(existing_points - selected_point, axis=1)))
            else:
                distance = 1.0
            score = contour_score * (1.0 + settings.exploration_weight * distance)
            candidates.append((score, config, selected_point))
            selected_points.append(selected_point)

            if len(candidates) >= n_candidates:
                return candidates

    return candidates


def propose_hybrid_configs(
    *,
    data: RMBData,
    settings: HybridBoundaryExperimentConfig,
    rng: RNGGenerator,
) -> list[RMBConfig]:
    """
    Propose a batch using contour following plus global boundary acquisition.
    """
    groups = grouped_by_n_qubits(data)
    candidates: list[tuple[float, RMBConfig, np.ndarray]] = []
    selected_points: list[np.ndarray] = []
    fitted_surface_found = False

    n_contour = int(round(settings.batch_size * settings.contour_follow_fraction))
    n_contour = max(0, min(settings.batch_size, n_contour))
    n_global = settings.batch_size - n_contour

    for n_qubits in settings.n_qubits_values:
        group = groups.get(n_qubits, {})
        template = template_config(settings, n_qubits)
        try:
            surface = fit_fidelity_surface(group, settings)
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
        candidates.extend(contour_candidates)

        existing_points = scaled_data_points(group, settings)
        for _ in range(max(n_global, settings.batch_size - len(contour_candidates))):
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
        unique_configs.extend(random_candidate_configs(
            settings=settings,
            rng=rng,
            n_candidates=settings.batch_size - len(unique_configs),
        ))

    return unique_configs


def estimate_boundary(settings: HybridBoundaryExperimentConfig) -> RMB:
    rng = default_rng(settings.rng_seed)
    backend = make_backend()
    rmb = RMB.default(rng).with_backend(backend)
    data: RMBData = rmb._data

    remaining = settings.measurement_budget
    initial_configs, initial_shots_per_config = budgeted_initial_design(settings, rng)

    for config in initial_configs:
        if remaining <= 0:
            break
        spent = spend_measurements(
            backend=backend,
            rng=rng,
            data=data,
            config=config,
            n_measurements=min(initial_shots_per_config, remaining),
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
        batch = propose_hybrid_configs(data=data, settings=settings, rng=rng)
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
                f"\nHybrid batch {batch_index}: spent {spent_this_batch}, "
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
    settings = HybridBoundaryExperimentConfig(
        measurement_budget=150,
        n_qubits_values=(5,),
        depth_bounds=(4, 200),
        ratio_bounds=(0.0, 0.8),
        initial_depths=4,
        initial_ratios=4,
        initial_shots_per_config=2,
        initial_budget_fraction=0.10,
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
        contour_follow_fraction=0.50,
        contour_ready_probability_width=0.10,
        contour_min_anchors=4,
        contour_step_fraction=0.08,
        contour_projection_fraction=0.12,
        contour_gradient_fraction=0.01,
        contour_candidate_multiplier=5,
        save_path="viarregio3_boundary.json",
        verbose=True,
    )

    rmb = estimate_boundary(settings)
    print_fit_reports(rmb._data, settings)
    plot_fidelity_surface_with_confidence(
        rmb._data,
        settings,
        n_bootstrap=100,
        seed=settings.rng_seed,
    )
