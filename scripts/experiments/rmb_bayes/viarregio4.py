from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from scipy.optimize import minimize
from scipy.special import expit

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData

from viarregio2 import (
    BoundaryExperimentConfig,
    adaptive_shot_count,
    budgeted_initial_design,
    fidelity_mean,
    grouped_by_n_qubits,
    make_backend,
    optimize_candidate_config,
    random_candidate_configs,
    scaled_data_points,
    spend_measurements,
    total_measurements,
)
from viarregio3 import (
    contour_follow_candidates,
    contour_points_from_surface,
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
    y = np.array([fidelity_mean(data[config]) for config in configs], dtype=float)
    weights = np.array([max(1, data[config].num_runs()) for config in configs], dtype=float)

    return fit_monotone_surface_from_arrays(
        points=points,
        y=y,
        weights=weights,
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
        candidates.extend(contour_candidates)

        existing_points = scaled_data_points(group, settings)
        for _ in range(max(n_global, settings.batch_size - len(contour_candidates))):
            try:
                score, config, selected_point = optimize_candidate_config(
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
        unique_configs.extend(random_candidate_configs(
            settings=settings,
            rng=rng,
            n_candidates=settings.batch_size - len(unique_configs),
        ))

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
        batch = propose_monotone_hybrid_configs(data=data, settings=settings, rng=rng)
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
                f"\nMonotone hybrid batch {batch_index}: spent {spent_this_batch}, "
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
    settings = MonotoneBoundaryExperimentConfig(
        measurement_budget=500,
        n_qubits_values=(50,),
        depth_bounds=(4, 300),
        ratio_bounds=(0.08, 0.8),
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
