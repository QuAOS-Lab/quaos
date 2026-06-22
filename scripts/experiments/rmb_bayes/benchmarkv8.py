from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from dataclasses import replace
from math import ceil
from pathlib import Path

import numpy as np
from numpy.random import default_rng

from sympleq.applications.randomized_benchmarking.RMB import RMB

from viarregio2 import config_from_parameters, fidelity_mean, make_backend, total_measurements
from viarregio4 import fidelity_colormap, hqc_cost
from viarregio7 import (
    BudgetState,
    MeasurementRequest,
    batch_config_key,
    spend_configs_batch,
    template_config,
)
from viarregio8 import (
    MonotoneFidelityVolume,
    VolumeExperimentConfig,
    depth_domain_for_qubits,
    estimate_boundary as estimate_boundary_v8,
    extract_volume_grid,
    fit_monotone_fidelity_volume,
    script_default_settings,
    volume_metrics,
)


def configure_plot_caches() -> None:
    cache_root = Path(tempfile.gettempdir()) / "sympleq_benchmarkv8_cache"
    mpl_cache = cache_root / "matplotlib"
    xdg_cache = cache_root / "xdg"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))


configure_plot_caches()


def parse_budget_list(value: str) -> list[float]:
    budgets = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not budgets:
        raise ValueError("At least one budget is required.")
    return budgets


def parse_optional_float(value: str) -> float | None:
    if value.strip().lower() in {"none", "null", "off"}:
        return None
    return float(value)


def make_settings(
    *,
    budget: float,
    seed: int,
    args: argparse.Namespace,
) -> VolumeExperimentConfig:
    if args.budget_mode == "hqc":
        measurement_budget = args.measurement_cap
        hqc_budget = float(budget)
    else:
        measurement_budget = int(round(budget))
        hqc_budget = args.hqc_budget

    base = script_default_settings(
        save_path=None,
        volume_save_path=None,
        diagnostics_path=(
            Path(args.diagnostics_dir).expanduser().resolve()
            / f"v8_{args.budget_mode}_{budget:g}_seed{seed}_diagnostics.json"
            if args.diagnostics_dir is not None
            else None
        ),
        print_diagnostics=args.print_diagnostics,
        verbose=args.verbose,
    )
    return replace(
        base,
        measurement_budget=measurement_budget,
        hqc_budget=hqc_budget,
        n_qubits_values=tuple(args.n_qubits_values),
        seed_n_qubits_values=tuple(args.seed_n_qubits_values),
        depth_bounds=(args.depth_min, args.depth_max),
        ratio_bounds=(args.ratio_min, args.ratio_max),
        random_elimination=args.random_elimination,
        scrambling_probability=args.scrambling_probability,
        max_shots_per_config=args.max_shots_per_config,
        min_fit_points=args.min_fit_points,
        min_volume_fit_points=args.min_volume_fit_points,
        monotone_l2=args.monotone_l2,
        volume_grid_size=(args.volume_depth_grid, args.volume_ratio_grid, args.volume_qubit_grid),
        acquisition_ratio_count=args.acquisition_ratio_count,
        acquisition_qubit_count=args.acquisition_qubit_count,
        acquisition_passes=args.acquisition_passes,
        acquisition_candidates_per_pass=args.acquisition_candidates_per_pass,
        initial_seed_ratio_count=args.initial_seed_ratio_count,
        initial_seed_depth_grid_count=args.initial_seed_depth_grid_count,
        initial_low_ratio_depth_fractions=tuple(args.initial_low_ratio_depth_fractions),
        initial_seed_shots=args.initial_seed_shots,
        qubit_depth_cap_enabled=args.qubit_depth_cap,
        qubit_depth_cap_plateau=args.qubit_depth_cap_plateau,
        qubit_depth_cap_power=args.qubit_depth_cap_power,
        qubit_depth_cap_margin_fraction=args.qubit_depth_cap_margin_fraction,
        qubit_depth_cap_acquisition=args.qubit_depth_cap_acquisition,
        contour_bracket_probe_shots=args.contour_bracket_probe_shots,
        contour_bracket_depth_fractions=tuple(args.contour_bracket_depth_fractions),
        contour_bracket_max_relative_depth=args.contour_bracket_max_relative_depth,
        batching_enabled=not args.no_batching,
        batch_max_configs=args.batch_max_configs,
        max_cost_per_batch=args.max_cost_per_batch,
        batch_target_fill_fraction=args.batch_target_fill_fraction,
        batch_fill_repeats=not args.no_batch_fill_repeats,
        batch_fill_max_shots_per_config=args.batch_fill_max_shots_per_config,
        batch_discovery_fill_max_shots_per_config=args.batch_discovery_fill_max_shots_per_config,
        batch_fill_all_stages=args.batch_fill_all_stages,
        volume_sparsity_radius=args.volume_sparsity_radius,
        volume_uncertainty_weight=args.volume_uncertainty_weight,
        volume_sparsity_weight=args.volume_sparsity_weight,
        volume_cost_power=args.volume_cost_power,
        volume_guardrails_enabled=not args.no_volume_guardrails,
        volume_guardrail_min_successes_per_qubit=args.volume_guardrail_min_successes_per_qubit,
        volume_guardrail_min_failures_per_qubit=args.volume_guardrail_min_failures_per_qubit,
        volume_guardrail_shots=args.volume_guardrail_shots,
        volume_guardrail_max_configs_per_pass=args.volume_guardrail_max_configs_per_pass,
        bracket_completion_enabled=not args.no_bracket_completion,
        bracket_completion_max_groups_per_pass=args.bracket_completion_max_groups_per_pass,
        bracket_completion_ratio_decimals=args.bracket_completion_ratio_decimals,
        bracket_completion_probability_width=args.bracket_completion_probability_width,
        bracket_completion_depth_fractions=tuple(args.bracket_completion_depth_fractions),
        hqc_cost_informed_acquisition=(args.budget_mode == "hqc" or hqc_budget is not None),
        hqc_cost_power=1.0,
        rng_seed=seed,
    )


def cache_stem(label: str, budget: float, seed: int, args: argparse.Namespace) -> str:
    q = "-".join(str(value) for value in args.n_qubits_values)
    seed_q = "-".join(str(value) for value in args.seed_n_qubits_values)
    reference_suffix = ""
    if label == "true":
        reference_suffix = (
            f"_rr{args.reference_ratios}"
            f"_rb{args.reference_bisection_steps}"
            f"_rs{args.reference_shots}"
            f"_rc{args.reference_confirm_shots}"
        )
    name = (
        f"v8_{label}_{args.budget_mode}{budget:g}_q{q}_seedq{seed_q}_"
        f"d{args.depth_min}-{args.depth_max}_r{args.ratio_min:g}-{args.ratio_max:g}_"
        f"seed{seed}{reference_suffix}_vg{args.volume_depth_grid}x{args.volume_ratio_grid}x{args.volume_qubit_grid}_"
        f"ap{args.acquisition_passes}ar{args.acquisition_ratio_count}aq{args.acquisition_qubit_count}"
        f"ac{args.acquisition_candidates_per_pass}"
        f"_il{','.join(f'{value:g}' for value in args.initial_low_ratio_depth_fractions)}"
        f"_cb{args.contour_bracket_probe_shots}"
        f"_dc{'y' if args.qubit_depth_cap else 'n'}"
        f"p{args.qubit_depth_cap_plateau}"
        f"w{args.qubit_depth_cap_power:g}"
        f"m{args.qubit_depth_cap_margin_fraction:g}"
        f"a{'y' if args.qubit_depth_cap_acquisition else 'n'}"
        f"_mc{args.max_cost_per_batch}"
        f"_vgd{'n' if args.no_volume_guardrails else 'y'}"
        f"s{args.volume_guardrail_min_successes_per_qubit}"
        f"f{args.volume_guardrail_min_failures_per_qubit}"
        f"h{args.volume_guardrail_shots}"
        f"c{args.volume_guardrail_max_configs_per_pass}"
        f"_bc{'n' if args.no_bracket_completion else 'y'}"
        f"g{args.bracket_completion_max_groups_per_pass}"
        f"r{args.bracket_completion_ratio_decimals}"
        f"p{args.bracket_completion_probability_width:g}"
    )
    if len(name) > 190:
        digest = hashlib.sha1(name.encode("utf-8")).hexdigest()[:12]
        name = f"v8_{label}_{args.budget_mode}{budget:g}_seed{seed}_{digest}"
    return name


def cache_path(label: str, budget: float, seed: int, args: argparse.Namespace) -> Path:
    return Path(args.cache_dir).expanduser().resolve() / f"{cache_stem(label, budget, seed, args)}.json"


def batch_summary_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".batch_summary.json")


def batch_config_ids_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".batch_config_ids.json")


def volume_grid_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".volume_grid.json")


def metrics_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".metrics.json")


def save_sidecars(path: Path, rmb: RMB) -> None:
    summary = getattr(rmb, "_batch_cost_summary", None)
    if summary:
        batch_summary_path(path).write_text(json.dumps(summary, indent=2), encoding="utf-8")
    ids = getattr(rmb, "_batch_config_ids", None)
    if ids:
        batch_config_ids_path(path).write_text(json.dumps(ids, indent=2), encoding="utf-8")
    grid = getattr(rmb, "_volume_grid", None)
    if grid:
        volume_grid_path(path).write_text(json.dumps(grid, indent=2), encoding="utf-8")


def load_sidecars(path: Path, rmb: RMB) -> None:
    if batch_summary_path(path).exists():
        rmb._batch_cost_summary = json.loads(batch_summary_path(path).read_text(encoding="utf-8"))
    if batch_config_ids_path(path).exists():
        rmb._batch_config_ids = json.loads(batch_config_ids_path(path).read_text(encoding="utf-8"))
    if volume_grid_path(path).exists():
        rmb._volume_grid = json.loads(volume_grid_path(path).read_text(encoding="utf-8"))


def spend_reference_request(
    *,
    rmb: RMB,
    config,
    shots: int,
    budget: BudgetState,
    settings: VolumeExperimentConfig,
) -> float | None:
    spends = spend_configs_batch(
        backend=rmb.backend,
        rng=rmb.rng,
        data=rmb._data,
        requests=[MeasurementRequest(config, shots)],
        budget=budget,
        settings=settings,
        stage="initial_grid",
    )
    if not spends or config not in rmb._data:
        return None
    return fidelity_mean(rmb._data[config])


def dense_reference_run(
    *,
    settings: VolumeExperimentConfig,
    seed: int,
    args: argparse.Namespace,
) -> RMB:
    rng = default_rng(seed)
    rmb = RMB.default(rng).with_backend(make_backend())
    budget = BudgetState(
        remaining_measurements=settings.measurement_budget,
        remaining_hqc=float(settings.hqc_budget) if settings.hqc_budget is not None else float(settings.measurement_budget),
    )
    ratios = np.linspace(args.ratio_min, args.ratio_max, args.reference_ratios)
    for n_qubits in settings.n_qubits_values:
        template = template_config(settings, n_qubits)
        for ratio in ratios:
            if not budget.can_spend():
                break
            low_depth, high_depth = depth_domain_for_qubits(n_qubits, settings)
            low_config = config_from_parameters(template=template, depth=low_depth, ratio=float(ratio))
            high_config = config_from_parameters(template=template, depth=high_depth, ratio=float(ratio))
            low_p = spend_reference_request(
                rmb=rmb,
                config=low_config,
                shots=args.reference_edge_shots,
                budget=budget,
                settings=settings,
            )
            high_p = spend_reference_request(
                rmb=rmb,
                config=high_config,
                shots=args.reference_edge_shots,
                budget=budget,
                settings=settings,
            )
            if low_p is None or high_p is None:
                break
            if low_p < 0.5 or high_p > 0.5:
                continue

            best_config = low_config if abs(low_p - 0.5) <= abs(high_p - 0.5) else high_config
            best_error = min(abs(low_p - 0.5), abs(high_p - 0.5))
            for _ in range(args.reference_bisection_steps):
                if not budget.can_spend():
                    break
                mid_depth = 0.5 * (low_depth + high_depth)
                mid_config = config_from_parameters(template=template, depth=mid_depth, ratio=float(ratio))
                mid_p = spend_reference_request(
                    rmb=rmb,
                    config=mid_config,
                    shots=args.reference_shots,
                    budget=budget,
                    settings=settings,
                )
                if mid_p is None:
                    break
                error = abs(mid_p - 0.5)
                if error < best_error:
                    best_config = mid_config
                    best_error = error
                if mid_p >= 0.5:
                    low_depth = float(mid_config.depth)
                else:
                    high_depth = float(mid_config.depth)
            if args.reference_confirm_shots > 0:
                spend_reference_request(
                    rmb=rmb,
                    config=best_config,
                    shots=args.reference_confirm_shots,
                    budget=budget,
                    settings=settings,
                )

    saving = budget.native_hqc_estimate - budget.stitched_hqc_spent
    rmb._batch_cost_summary = {
        "stitched_hqc_spent": budget.stitched_hqc_spent,
        "native_hqc_estimate": budget.native_hqc_estimate,
        "estimated_batching_saving_hqc": saving,
        "estimated_batching_saving_fraction": saving / budget.native_hqc_estimate if budget.native_hqc_estimate > 0 else 0.0,
        "batched_jobs": budget.batched_jobs,
        "max_batched_job_size": budget.max_batched_job_size,
        "max_cost_per_batch": settings.max_cost_per_batch,
    }
    rmb._batch_config_ids = budget.batch_config_ids
    try:
        surface = fit_monotone_fidelity_volume(rmb._data, settings)
    except (RuntimeError, ValueError):
        pass
    else:
        rmb._volume_grid = extract_volume_grid(surface, settings)
    return rmb


def run_or_load(
    *,
    label: str,
    budget: float,
    seed: int,
    args: argparse.Namespace,
) -> tuple[RMB, VolumeExperimentConfig]:
    settings = make_settings(budget=budget, seed=seed, args=args)
    path = cache_path(label, budget, seed, args)
    if path.exists() and not args.rebuild:
        print(f"Loading {label} budget {budget:g} from {path}")
        rmb = RMB.load(path).with_backend(make_backend())
        load_sidecars(path, rmb)
        return rmb, settings

    print(f"Running {label} budget {budget:g}...")
    if label == "true":
        rmb = dense_reference_run(settings=settings, seed=seed, args=args)
    else:
        rmb = estimate_boundary_v8(settings)
    path.parent.mkdir(parents=True, exist_ok=True)
    rmb.save(path)
    save_sidecars(path, rmb)
    print(f"Saved {label} budget {budget:g} to {path}")
    summary = getattr(rmb, "_batch_cost_summary", None)
    if summary:
        print(
            f"{label} budget {budget:g}: stitched HQC {summary['stitched_hqc_spent']:.1f}, "
            f"native {summary['native_hqc_estimate']:.1f}, "
            f"saving {summary['estimated_batching_saving_fraction']:.1%}"
        )
    return rmb, settings


def fitted_volume_and_grid(rmb: RMB, settings: VolumeExperimentConfig) -> tuple[MonotoneFidelityVolume | None, dict | None]:
    grid = getattr(rmb, "_volume_grid", None)
    try:
        surface = fit_monotone_fidelity_volume(rmb._data, settings)
    except (RuntimeError, ValueError):
        return None, grid
    if grid is None:
        grid = extract_volume_grid(surface, settings)
        rmb._volume_grid = grid
    return surface, grid


def hqc_spent(rmb: RMB, settings: VolumeExperimentConfig) -> float:
    summary = getattr(rmb, "_batch_cost_summary", None)
    if summary:
        return float(summary["stitched_hqc_spent"])
    if settings.hqc_budget is None:
        return float("nan")
    return float(sum(
        hqc_cost(config, estimator.num_runs(), settings)
        for config, estimator in rmb._data.items()
        if estimator.num_runs() > 0
    ))


def point_arrays(rmb: RMB) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    configs = [config for config, estimator in rmb._data.items() if estimator.num_runs() > 0]
    depths = np.asarray([config.depth for config in configs], dtype=float)
    ratios = np.asarray([config.min_two_qubit_gate_ratio for config in configs], dtype=float)
    qubits = np.asarray([config.n_qubits for config in configs], dtype=float)
    fidelities = np.asarray([fidelity_mean(rmb._data[config]) for config in configs], dtype=float)
    return depths, ratios, qubits, fidelities


def plot_3d_surface(ax, *, rmb: RMB, settings: VolumeExperimentConfig, title: str) -> None:
    surface, grid = fitted_volume_and_grid(rmb, settings)
    if grid is None:
        ax.set_title(f"{title}\nno 3D fit")
        return
    ratios = np.asarray(grid["ratios"], dtype=float)
    qubits = np.asarray(grid["n_qubits"], dtype=float)
    depth50 = np.asarray(grid["depth50"], dtype=float)
    ratio_grid, qubit_grid = np.meshgrid(ratios, qubits)
    ax.plot_surface(depth50, ratio_grid, qubit_grid, color="#4c78a8", alpha=0.42, linewidth=0, antialiased=True)
    depths, point_ratios, point_qubits, fidelities = point_arrays(rmb)
    ax.scatter(depths, point_ratios, point_qubits, c=fidelities, cmap=fidelity_colormap(), vmin=0, vmax=1, s=12, edgecolors="black", linewidths=0.2)
    ax.set_xlabel("Depth")
    ax.set_ylabel("Two-qubit ratio")
    ax.set_zlabel("n_qubits")
    ax.set_title(f"{title}\nvolume={grid['normalized_volume']:.4f}, HQC={hqc_spent(rmb, settings):.1f}")


def plot_slice_panel(ax, *, rmb: RMB, settings: VolumeExperimentConfig, n_qubits: int, true_grid: dict | None) -> None:
    surface, grid = fitted_volume_and_grid(rmb, settings)
    cmap = fidelity_colormap()
    d_min, d_max = depth_domain_for_qubits(n_qubits, settings)
    depths = np.linspace(d_min, d_max, settings.volume_grid_size[0])
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], settings.volume_grid_size[1])
    depth_grid, ratio_grid = np.meshgrid(depths, ratios)
    if surface is not None:
        points = np.column_stack([
            depth_grid.ravel(),
            ratio_grid.ravel(),
            np.full(depth_grid.size, n_qubits, dtype=float),
        ])
        probabilities = surface.probability(points).reshape(depth_grid.shape)
        ax.contourf(depth_grid, ratio_grid, probabilities, levels=np.linspace(0, 1, 16), cmap=cmap, alpha=0.28, vmin=0, vmax=1)
    if grid is not None:
        qubits = np.asarray(grid["n_qubits"], dtype=float)
        ratios_line = np.asarray(grid["ratios"], dtype=float)
        depth50 = np.asarray(grid["depth50"], dtype=float)
        q_idx = int(np.argmin(np.abs(qubits - n_qubits)))
        ax.plot(depth50[q_idx], ratios_line, color="black", linewidth=2, label="estimate")
    if true_grid is not None:
        qubits = np.asarray(true_grid["n_qubits"], dtype=float)
        ratios_line = np.asarray(true_grid["ratios"], dtype=float)
        depth50 = np.asarray(true_grid["depth50"], dtype=float)
        q_idx = int(np.argmin(np.abs(qubits - n_qubits)))
        ax.plot(depth50[q_idx], ratios_line, color="red", linewidth=1.5, label="reference")
    configs = [
        config
        for config, estimator in rmb._data.items()
        if estimator.num_runs() > 0 and config.n_qubits == n_qubits
    ]
    if configs:
        ax.scatter(
            [config.depth for config in configs],
            [config.min_two_qubit_gate_ratio for config in configs],
            c=[fidelity_mean(rmb._data[config]) for config in configs],
            cmap=cmap,
            vmin=0,
            vmax=1,
            s=28,
            edgecolors="black",
            linewidths=0.3,
        )
    ax.set_title(f"{n_qubits} qubits")
    ax.set_xlabel("Depth")
    ax.set_ylabel("Two-qubit ratio")
    ax.set_xlim(d_min, d_max)
    ax.set_ylim(settings.ratio_bounds)
    ax.legend(loc="best", fontsize=7)


def make_plots(
    *,
    runs: list[tuple[str, RMB, VolumeExperimentConfig, dict]],
    true_run: tuple[str, RMB, VolumeExperimentConfig],
    args: argparse.Namespace,
) -> None:
    import matplotlib.pyplot as plt

    true_label, true_rmb, true_settings = true_run
    true_surface, true_grid = fitted_volume_and_grid(true_rmb, true_settings)
    n_rows = len(runs) + 1
    n_slice_cols = len(args.slice_n_qubits)
    fig = plt.figure(figsize=(5.0 * (n_slice_cols + 1), 4.2 * n_rows), constrained_layout=True)

    plot_rows = [(true_label, true_rmb, true_settings, {})] + runs
    for row, (label, rmb, settings, metrics) in enumerate(plot_rows):
        ax3d = fig.add_subplot(n_rows, n_slice_cols + 1, row * (n_slice_cols + 1) + 1, projection="3d")
        plot_3d_surface(ax3d, rmb=rmb, settings=settings, title=label)
        for col, n_qubits in enumerate(args.slice_n_qubits, start=2):
            ax = fig.add_subplot(n_rows, n_slice_cols + 1, row * (n_slice_cols + 1) + col)
            plot_slice_panel(ax, rmb=rmb, settings=settings, n_qubits=n_qubits, true_grid=true_grid if label != true_label else None)
            if metrics and col == n_slice_cols + 1:
                text = (
                    f"vol={metrics.get('normalized_volume', float('nan')):.4f}\n"
                    f"vol err={metrics.get('absolute_volume_error', float('nan')):.4f}\n"
                    f"depth err={metrics.get('mean_abs_depth50_error', float('nan')):.2f}\n"
                    f"cal={metrics.get('calibration_error', float('nan')):.3f}"
                )
                ax.text(0.02, 0.02, text, transform=ax.transAxes, fontsize=8, va="bottom", bbox={"facecolor": "white", "alpha": 0.75})

    output = Path(args.output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=args.dpi)
    print(f"Saved plot to {output}")
    if not args.no_show:
        plt.show()
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    defaults = script_default_settings()
    parser = argparse.ArgumentParser(description="Benchmark viarregio8 3D HQC-aware fidelity volume estimation.")
    parser.add_argument("--budgets", type=str, default="250,500,1000")
    parser.add_argument("--true-budget", type=float, default=5000.0)
    parser.add_argument("--budget-mode", choices=("hqc", "measurements"), default="hqc")
    parser.add_argument("--hqc-budget", type=float, default=None)
    parser.add_argument("--measurement-cap", type=int, default=1_000_000_000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--true-seed-offset", type=int, default=100_000)
    parser.add_argument("--n-qubits-values", type=int, nargs="+", default=list(defaults.n_qubits_values))
    parser.add_argument("--seed-n-qubits-values", type=int, nargs="+", default=list(defaults.seed_n_qubits_values))
    parser.add_argument("--slice-n-qubits", type=int, nargs="+", default=[10, 30, 50])
    parser.add_argument("--depth-min", type=int, default=defaults.depth_bounds[0])
    parser.add_argument("--depth-max", type=int, default=defaults.depth_bounds[1])
    parser.add_argument("--ratio-min", type=float, default=defaults.ratio_bounds[0])
    parser.add_argument("--ratio-max", type=float, default=defaults.ratio_bounds[1])
    parser.add_argument("--random-elimination", type=float, default=defaults.random_elimination)
    parser.add_argument("--scrambling-probability", type=float, default=defaults.scrambling_probability)
    parser.add_argument("--max-shots-per-config", type=int, default=defaults.max_shots_per_config)
    parser.add_argument("--min-fit-points", type=int, default=defaults.min_fit_points)
    parser.add_argument("--min-volume-fit-points", type=int, default=defaults.min_volume_fit_points)
    parser.add_argument("--monotone-l2", type=float, default=defaults.monotone_l2)
    parser.add_argument("--volume-depth-grid", type=int, default=defaults.volume_grid_size[0])
    parser.add_argument("--volume-ratio-grid", type=int, default=defaults.volume_grid_size[1])
    parser.add_argument("--volume-qubit-grid", type=int, default=defaults.volume_grid_size[2])
    parser.add_argument("--acquisition-ratio-count", type=int, default=defaults.acquisition_ratio_count)
    parser.add_argument("--acquisition-qubit-count", type=int, default=defaults.acquisition_qubit_count)
    parser.add_argument("--acquisition-passes", type=int, default=defaults.acquisition_passes)
    parser.add_argument("--acquisition-candidates-per-pass", type=int, default=defaults.acquisition_candidates_per_pass)
    parser.add_argument("--initial-seed-ratio-count", type=int, default=defaults.initial_seed_ratio_count)
    parser.add_argument("--initial-seed-depth-grid-count", type=int, default=defaults.initial_seed_depth_grid_count)
    parser.add_argument("--initial-low-ratio-depth-fractions", type=float, nargs="+", default=list(defaults.initial_low_ratio_depth_fractions))
    parser.add_argument("--initial-seed-shots", type=int, default=defaults.initial_seed_shots)
    parser.add_argument("--qubit-depth-cap", action="store_true", default=defaults.qubit_depth_cap_enabled)
    parser.add_argument("--no-qubit-depth-cap", dest="qubit_depth_cap", action="store_false")
    parser.add_argument("--qubit-depth-cap-plateau", type=int, default=defaults.qubit_depth_cap_plateau)
    parser.add_argument("--qubit-depth-cap-power", type=float, default=defaults.qubit_depth_cap_power)
    parser.add_argument("--qubit-depth-cap-margin-fraction", type=float, default=defaults.qubit_depth_cap_margin_fraction)
    parser.add_argument("--qubit-depth-cap-acquisition", action="store_true")
    parser.add_argument("--reference-ratios", type=int, default=12)
    parser.add_argument("--reference-edge-shots", type=int, default=3)
    parser.add_argument("--reference-bisection-steps", type=int, default=5)
    parser.add_argument("--reference-shots", type=int, default=4)
    parser.add_argument("--reference-confirm-shots", type=int, default=6)
    parser.add_argument("--contour-bracket-probe-shots", type=int, default=defaults.contour_bracket_probe_shots)
    parser.add_argument("--contour-bracket-depth-fractions", type=float, nargs="+", default=list(defaults.contour_bracket_depth_fractions))
    parser.add_argument("--contour-bracket-max-relative-depth", type=float, default=defaults.contour_bracket_max_relative_depth)
    parser.add_argument("--no-batching", action="store_true")
    parser.add_argument("--batch-max-configs", type=int, default=defaults.batch_max_configs)
    parser.add_argument("--max-cost-per-batch", type=parse_optional_float, default=defaults.max_cost_per_batch)
    parser.add_argument("--batch-target-fill-fraction", type=float, default=defaults.batch_target_fill_fraction)
    parser.add_argument("--no-batch-fill-repeats", action="store_true")
    parser.add_argument("--batch-fill-max-shots-per-config", type=int, default=defaults.batch_fill_max_shots_per_config)
    parser.add_argument("--batch-discovery-fill-max-shots-per-config", type=int, default=defaults.batch_discovery_fill_max_shots_per_config)
    parser.add_argument("--batch-fill-all-stages", action="store_true")
    parser.add_argument("--volume-sparsity-radius", type=float, default=defaults.volume_sparsity_radius)
    parser.add_argument("--volume-uncertainty-weight", type=float, default=defaults.volume_uncertainty_weight)
    parser.add_argument("--volume-sparsity-weight", type=float, default=defaults.volume_sparsity_weight)
    parser.add_argument("--volume-cost-power", type=float, default=defaults.volume_cost_power)
    parser.add_argument("--no-volume-guardrails", action="store_true")
    parser.add_argument("--volume-guardrail-min-successes-per-qubit", type=int, default=defaults.volume_guardrail_min_successes_per_qubit)
    parser.add_argument("--volume-guardrail-min-failures-per-qubit", type=int, default=defaults.volume_guardrail_min_failures_per_qubit)
    parser.add_argument("--volume-guardrail-shots", type=int, default=defaults.volume_guardrail_shots)
    parser.add_argument("--volume-guardrail-max-configs-per-pass", type=int, default=defaults.volume_guardrail_max_configs_per_pass)
    parser.add_argument("--no-bracket-completion", action="store_true")
    parser.add_argument("--bracket-completion-max-groups-per-pass", type=int, default=defaults.bracket_completion_max_groups_per_pass)
    parser.add_argument("--bracket-completion-ratio-decimals", type=int, default=defaults.bracket_completion_ratio_decimals)
    parser.add_argument("--bracket-completion-probability-width", type=float, default=defaults.bracket_completion_probability_width)
    parser.add_argument("--bracket-completion-depth-fractions", type=float, nargs="+", default=list(defaults.bracket_completion_depth_fractions))
    parser.add_argument("--diagnostics-dir", type=str, default=None)
    parser.add_argument("--print-diagnostics", action="store_true")
    parser.add_argument("--cache-dir", type=str, default=str(Path(__file__).resolve().parent / "benchmark_cache"))
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--output", type=str, default=str(Path(__file__).resolve().parent / "benchmarkv8.png"))
    parser.add_argument("--metrics-output", type=str, default=str(Path(__file__).resolve().parent / "benchmarkv8_metrics.json"))
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main(**overrides) -> None:
    args = parse_args()
    for name, value in overrides.items():
        if not hasattr(args, name):
            raise ValueError(f"Unknown benchmarkv8 option: {name}")
        setattr(args, name, value)
    if args.overwrite:
        args.rebuild = True
    budgets = parse_budget_list(args.budgets)

    true_seed = args.seed + args.true_seed_offset
    true_run = (
        f"true {args.budget_mode} {args.true_budget:g}",
        *run_or_load(label="true", budget=args.true_budget, seed=true_seed, args=args),
    )
    true_surface, true_grid = fitted_volume_and_grid(true_run[1], true_run[2])

    runs = []
    metrics_by_label = {}
    for budget in budgets:
        label = f"{args.budget_mode} {budget:g}"
        rmb, settings = run_or_load(label="budget", budget=budget, seed=args.seed, args=args)
        surface, grid = fitted_volume_and_grid(rmb, settings)
        metrics = {}
        if grid is not None:
            metrics = volume_metrics(grid, reference_grid=true_grid, reference_surface=true_surface)
            metrics["hqc_spent"] = hqc_spent(rmb, settings)
            metrics["measurements"] = total_measurements(rmb._data)
            metrics["configs"] = len(rmb._data)
            metrics_by_label[label] = metrics
            path = metrics_path(cache_path("budget", budget, args.seed, args))
            path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
            print(
                f"{label}: volume={metrics.get('normalized_volume', float('nan')):.4f}, "
                f"volume_error={metrics.get('absolute_volume_error', float('nan')):.4f}, "
                f"depth_error={metrics.get('mean_abs_depth50_error', float('nan')):.3f}, "
                f"calibration={metrics.get('calibration_error', float('nan')):.4f}"
            )
        runs.append((label, rmb, settings, metrics))

    Path(args.metrics_output).expanduser().resolve().write_text(json.dumps(metrics_by_label, indent=2), encoding="utf-8")
    make_plots(runs=runs, true_run=true_run, args=args)


if __name__ == "__main__":
    main()
