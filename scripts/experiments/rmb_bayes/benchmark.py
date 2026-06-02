from __future__ import annotations

import argparse
from dataclasses import replace
from pathlib import Path
from typing import Callable

import numpy as np

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData

from viarregio2 import (
    BoundaryExperimentConfig,
    FlexibleFidelitySurface,
    config_from_parameters,
    fit_fidelity_surface,
    make_backend,
    spend_measurements,
)
from viarregio2 import estimate_boundary as estimate_boundary_v2
from viarregio3 import HybridBoundaryExperimentConfig
from viarregio3 import estimate_boundary as estimate_boundary_v3


StrategyRunner = Callable[[BoundaryExperimentConfig], RMB]


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

    # Horizontal grid edges.
    for row in range(probabilities.shape[0]):
        for col in range(probabilities.shape[1] - 1):
            p0 = probabilities[row, col] - level
            p1 = probabilities[row, col + 1] - level
            if p0 == 0.0:
                points.append([depth_grid[row, col], ratio_grid[row, col]])
            if p0 * p1 < 0.0:
                t = abs(p0) / (abs(p0) + abs(p1))
                depth = (1.0 - t) * depth_grid[row, col] + t * depth_grid[row, col + 1]
                ratio = ratio_grid[row, col]
                points.append([depth, ratio])

    # Vertical grid edges.
    for row in range(probabilities.shape[0] - 1):
        for col in range(probabilities.shape[1]):
            p0 = probabilities[row, col] - level
            p1 = probabilities[row + 1, col] - level
            if p0 == 0.0:
                points.append([depth_grid[row, col], ratio_grid[row, col]])
            if p0 * p1 < 0.0:
                t = abs(p0) / (abs(p0) + abs(p1))
                depth = depth_grid[row, col]
                ratio = (1.0 - t) * ratio_grid[row, col] + t * ratio_grid[row + 1, col]
                points.append([depth, ratio])

    if not points:
        return np.empty((0, 2), dtype=float)
    return np.unique(np.asarray(points, dtype=float), axis=0)


def scale_points(points: np.ndarray, settings: BoundaryExperimentConfig) -> np.ndarray:
    lower = np.array([settings.depth_bounds[0], settings.ratio_bounds[0]], dtype=float)
    upper = np.array([settings.depth_bounds[1], settings.ratio_bounds[1]], dtype=float)
    return (points - lower) / (upper - lower)


def mean_nearest_distance(source: np.ndarray, target: np.ndarray) -> float:
    if len(source) == 0 or len(target) == 0:
        return float("nan")
    distances = np.linalg.norm(source[:, None, :] - target[None, :, :], axis=2)
    return float(np.mean(np.min(distances, axis=1)))


def chamfer_distance(
    estimated_points: np.ndarray,
    reference_points: np.ndarray,
    settings: BoundaryExperimentConfig,
) -> float:
    if len(estimated_points) == 0 or len(reference_points) == 0:
        return float("nan")
    estimated_scaled = scale_points(estimated_points, settings)
    reference_scaled = scale_points(reference_points, settings)
    return 0.5 * (
        mean_nearest_distance(estimated_scaled, reference_scaled)
        + mean_nearest_distance(reference_scaled, estimated_scaled)
    )


def calibration_error(
    estimated_points: np.ndarray,
    reference_surface: FlexibleFidelitySurface,
) -> float:
    if len(estimated_points) == 0:
        return float("nan")
    reference_probabilities = reference_surface.probability(estimated_points)
    return float(np.mean(np.abs(reference_probabilities - 0.5)))


def build_reference_data(
    settings: BoundaryExperimentConfig,
    *,
    grid_size: tuple[int, int],
    shots_per_config: int,
    seed: int,
) -> RMBData:
    rng = np.random.default_rng(seed)
    backend = make_backend()
    data: RMBData = {}

    depths = np.linspace(settings.depth_bounds[0], settings.depth_bounds[1], grid_size[0])
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], grid_size[1])

    for n_qubits in settings.n_qubits_values:
        template = (
            RMBConfig.default()
            .with_n_qubits(n_qubits)
            .with_random_elimination(settings.random_elimination)
            .with_scrambling_probability(settings.scrambling_probability)
        )
        for depth in depths:
            for ratio in ratios:
                config = config_from_parameters(template=template, depth=depth, ratio=ratio)
                spend_measurements(
                    backend=backend,
                    rng=rng,
                    data=data,
                    config=config,
                    n_measurements=shots_per_config,
                )

    return data


def reference_cache_path(args: argparse.Namespace) -> Path:
    if args.reference_cache is not None:
        return Path(args.reference_cache).expanduser().resolve()
    name = (
        f"benchmark_reference_"
        f"{args.n_qubits}q_"
        f"d{args.depth_min}-{args.depth_max}_"
        f"r{args.ratio_min:g}-{args.ratio_max:g}_"
        f"g{args.reference_grid}_"
        f"s{args.reference_shots}_"
        f"seed{args.seed + 100_000}.json"
    )
    return Path(__file__).resolve().parent / "benchmark_cache" / name


def load_or_build_reference_data(
    reference_settings: BoundaryExperimentConfig,
    args: argparse.Namespace,
) -> RMBData:
    path = reference_cache_path(args)
    if path.exists() and not args.rebuild_reference:
        print(f"Loading reference data from {path}...")
        return RMB.load(path)._data

    print(
        "Building reference surface "
        f"({args.reference_grid}x{args.reference_grid}, "
        f"{args.reference_shots} shots/config)..."
    )
    data = build_reference_data(
        reference_settings,
        grid_size=(args.reference_grid, args.reference_grid),
        shots_per_config=args.reference_shots,
        seed=args.seed + 100_000,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    rmb = RMB.default().with_backend(make_backend())
    rmb._data = data
    rmb.save(path)
    print(f"Saved reference data to {path}.")
    return data


def make_v2_settings(seed: int, budget: int, args: argparse.Namespace) -> BoundaryExperimentConfig:
    return BoundaryExperimentConfig(
        measurement_budget=budget,
        n_qubits_values=(args.n_qubits,),
        depth_bounds=(args.depth_min, args.depth_max),
        ratio_bounds=(args.ratio_min, args.ratio_max),
        initial_depths=4,
        initial_ratios=4,
        initial_shots_per_config=2,
        initial_budget_fraction=0.25,
        reserve_adaptive_measurements=max(0, int(0.6 * budget)),
        initial_edge_margin=0.08,
        initial_candidate_multiplier=8,
        include_bracketing_points=True,
        batch_size=4,
        min_adaptive_shots_per_config=1,
        max_adaptive_shots_per_config=5,
        max_shots_per_config=10,
        candidate_grid_size=(args.contour_grid, args.contour_grid),
        optimizer_maxiter=args.optimizer_maxiter,
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
        surface_smoothing=args.surface_smoothing,
        min_fit_points=8,
        rng_seed=seed,
        save_path=None,
        verbose=False,
    )


def make_v3_settings(seed: int, budget: int, args: argparse.Namespace) -> HybridBoundaryExperimentConfig:
    base = make_v2_settings(seed, budget, args)
    params = base.__dict__.copy()
    params.update({
        "contour_follow_fraction": 0.50,
        "contour_ready_probability_width": 0.10,
        "contour_min_anchors": 4,
        "contour_step_fraction": 0.08,
        "contour_projection_fraction": 0.12,
        "contour_gradient_fraction": 0.01,
        "contour_candidate_multiplier": 5,
        "save_path": None,
    })
    return HybridBoundaryExperimentConfig(**params)


def evaluate_strategy(
    *,
    name: str,
    runner: StrategyRunner,
    settings: BoundaryExperimentConfig,
    reference_surface: FlexibleFidelitySurface,
    reference_contour: np.ndarray,
) -> dict[str, float | int | str]:
    rmb = runner(settings)
    try:
        estimated_surface = fit_fidelity_surface(rmb._data, settings)
        estimated_contour = contour_points_from_surface(estimated_surface, settings)
    except (RuntimeError, ValueError):
        estimated_contour = np.empty((0, 2), dtype=float)

    return {
        "strategy": name,
        "seed": settings.rng_seed if settings.rng_seed is not None else -1,
        "budget": settings.measurement_budget,
        "measurements": sum(estimator.num_runs() for estimator in rmb._data.values()),
        "configs": len(rmb._data),
        "contour_points": len(estimated_contour),
        "calibration_error": calibration_error(estimated_contour, reference_surface),
        "chamfer_distance": chamfer_distance(estimated_contour, reference_contour, settings),
        "failure": int(len(estimated_contour) == 0),
    }


def print_table(rows: list[dict[str, float | int | str]]) -> None:
    headers = [
        "strategy",
        "seed",
        "budget",
        "measurements",
        "configs",
        "contour_points",
        "calibration_error",
        "chamfer_distance",
        "failure",
    ]
    print(",".join(headers))
    for row in rows:
        values = []
        for header in headers:
            value = row[header]
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        print(",".join(values))


def print_summary(rows: list[dict[str, float | int | str]]) -> None:
    print("\nSummary")
    for strategy in sorted({str(row["strategy"]) for row in rows}):
        strategy_rows = [row for row in rows if row["strategy"] == strategy]
        calibration = np.array([float(row["calibration_error"]) for row in strategy_rows])
        chamfer = np.array([float(row["chamfer_distance"]) for row in strategy_rows])
        failures = np.array([int(row["failure"]) for row in strategy_rows])

        print(
            f"{strategy}: "
            f"calibration median={np.nanmedian(calibration):.5f}, "
            f"calibration mean={np.nanmean(calibration):.5f}, "
            f"chamfer median={np.nanmedian(chamfer):.5f}, "
            f"failure rate={np.mean(failures):.2f}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark viarregio2 and viarregio3 against a high-shot reference contour."
    )
    parser.add_argument("--budget", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--n-qubits", type=int, default=10)
    parser.add_argument("--depth-min", type=int, default=4)
    parser.add_argument("--depth-max", type=int, default=80)
    parser.add_argument("--ratio-min", type=float, default=0.0)
    parser.add_argument("--ratio-max", type=float, default=0.8)
    parser.add_argument("--reference-grid", type=int, default=14)
    parser.add_argument("--reference-shots", type=int, default=80)
    parser.add_argument("--contour-grid", type=int, default=80)
    parser.add_argument("--optimizer-maxiter", type=int, default=25)
    parser.add_argument("--surface-smoothing", type=float, default=0.1)
    parser.add_argument(
        "--reference-cache",
        type=str,
        default=None,
        help="Path to save/load the high-shot reference data. Defaults to scripts/personal/benchmark_cache/...",
    )
    parser.add_argument(
        "--rebuild-reference",
        action="store_true",
        help="Rebuild the high-shot reference even if a cache file exists.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reference_settings = make_v2_settings(args.seed, args.budget, args)
    reference_settings = replace(
        reference_settings,
        min_fit_points=min(8, args.reference_grid * args.reference_grid),
        candidate_grid_size=(args.contour_grid, args.contour_grid),
    )

    reference_data = load_or_build_reference_data(reference_settings, args)
    reference_surface = fit_fidelity_surface(reference_data, reference_settings)
    reference_contour = contour_points_from_surface(reference_surface, reference_settings)
    print(f"Reference contour points: {len(reference_contour)}\n")

    rows = []
    for repeat in range(args.repeats):
        seed = args.seed + repeat
        v2_settings = make_v2_settings(seed, args.budget, args)
        v3_settings = make_v3_settings(seed, args.budget, args)

        rows.append(evaluate_strategy(
            name="viarregio2_global",
            runner=estimate_boundary_v2,
            settings=v2_settings,
            reference_surface=reference_surface,
            reference_contour=reference_contour,
        ))
        rows.append(evaluate_strategy(
            name="viarregio3_hybrid",
            runner=estimate_boundary_v3,
            settings=v3_settings,
            reference_surface=reference_surface,
            reference_contour=reference_contour,
        ))

    print_table(rows)
    print_summary(rows)


if __name__ == "__main__":
    main()
