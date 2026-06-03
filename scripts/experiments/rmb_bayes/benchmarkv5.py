from __future__ import annotations

import argparse
from math import ceil
from pathlib import Path

import numpy as np
from numpy.random import default_rng

from sympleq.applications.randomized_benchmarking.RMB import RMB

from viarregio2 import (
    config_from_parameters,
    fidelity_mean,
    fidelity_variance,
    make_backend,
    spend_measurements,
    template_config,
    total_measurements,
)
from viarregio3 import contour_points_from_surface
from viarregio4 import (
    affordable_shot_count,
    fidelity_colormap,
    fit_monotone_fidelity_surface,
    hqc_cost,
)
from viarregio5 import ContourFirstExperimentConfig
from viarregio5 import estimate_boundary as estimate_boundary_v5


def parse_budget_list(value: str) -> list[float]:
    budgets = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not budgets:
        raise ValueError("At least one budget is required.")
    if any(budget <= 0.0 for budget in budgets):
        raise ValueError("Budgets must be positive.")
    return budgets


def make_settings(
    *,
    budget: float,
    seed: int,
    args: argparse.Namespace,
) -> ContourFirstExperimentConfig:
    if args.budget_mode == "hqc":
        measurement_budget = args.measurement_cap
        hqc_budget = float(budget)
    else:
        measurement_budget = int(round(budget))
        hqc_budget = args.hqc_budget

    return ContourFirstExperimentConfig(
        measurement_budget=measurement_budget,
        hqc_budget=hqc_budget,
        n_qubits_values=(args.n_qubits,),
        depth_bounds=(args.depth_min, args.depth_max),
        ratio_bounds=(args.ratio_min, args.ratio_max),
        random_elimination=args.random_elimination,
        scrambling_probability=args.scrambling_probability,
        min_adaptive_shots_per_config=1,
        max_adaptive_shots_per_config=args.max_adaptive_shots,
        max_shots_per_config=args.max_shots_per_config,
        candidate_grid_size=(args.contour_grid, args.contour_grid),
        boundary_width=0.08,
        shot_boundary_width=0.2,
        surface_smoothing=args.surface_smoothing,
        min_fit_points=8,
        monotone_l2=1e-3,
        contour_ready_probability_width=0.10,
        contour_min_anchors=4,
        contour_step_fraction=0.08,
        contour_projection_fraction=0.12,
        contour_gradient_fraction=0.01,
        ray_ratio_count=args.ray_ratio_count,
        ray_probe_shots=2,
        ray_bisection_steps=5,
        ray_bisection_shots=2,
        crossing_decision_confirm_shots=args.crossing_decision_confirm_shots,
        crossing_decision_probability_width=args.crossing_decision_probability_width,
        crossing_ambiguous_as_failure=not args.no_crossing_ambiguous_as_failure,
        initial_crossing_interior_bracket=not args.no_initial_crossing_interior_bracket,
        initial_crossing_center_fraction=args.initial_crossing_center_fraction,
        initial_crossing_half_width_fraction=args.initial_crossing_half_width_fraction,
        initial_crossing_expand_factor=args.initial_crossing_expand_factor,
        crossing_confirm_candidates=args.crossing_confirm_candidates,
        crossing_confirm_shots=args.crossing_confirm_shots,
        initial_anchor_refine=not args.no_initial_anchor_refine,
        initial_anchor_refine_shots=args.initial_anchor_refine_shots,
        initial_anchor_refine_steps=args.initial_anchor_refine_steps,
        initial_anchor_search_fraction=args.initial_anchor_search_fraction,
        initial_anchor_min_runs=args.initial_anchor_min_runs,
        trace_reserve_budget_fraction=args.trace_reserve_budget_fraction,
        trace_reserve_min_hqc=args.trace_reserve_min_hqc,
        trace_ratio_step_fraction=args.trace_ratio_step_fraction,
        trace_depth_search_fraction=0.06,
        trace_local_crossing=args.trace_local_crossing,
        trace_enforce_monotone_depth=not args.no_trace_monotone_depth,
        trace_correction_steps=args.trace_correction_steps,
        trace_shots=2,
        trace_accept_probability_width=0.12,
        trace_reject_probability_width=args.trace_reject_probability_width,
        trace_directions=(1,),
        model_projection_after_fit=True,
        refine_after_trace=True,
        trace_anchor_min_shots=args.trace_anchor_min_shots,
        refinement_shots=args.refinement_shots,
        refinement_boundary_width=0.15,
        hqc_cost_informed_acquisition=(args.budget_mode == "hqc" or hqc_budget is not None),
        hqc_cost_power=1.0,
        rng_seed=seed,
        save_path=None,
        diagnostics_path=(
            Path(args.diagnostics_dir).expanduser().resolve()
            / f"{args.budget_mode}_{budget:g}_seed{seed}_diagnostics.json"
            if args.diagnostics_dir is not None
            else None
        ),
        print_diagnostics=args.print_diagnostics,
        verbose=args.verbose,
    )


def spend_for_reference(
    *,
    rmb: RMB,
    config,
    requested_shots: int,
    settings: ContourFirstExperimentConfig,
    remaining_measurements: int,
    remaining_hqc: float,
) -> tuple[int, int, float]:
    if requested_shots <= 0 or remaining_measurements <= 0 or remaining_hqc <= 0.0:
        return 0, remaining_measurements, remaining_hqc

    shots = min(requested_shots, remaining_measurements)
    if settings.hqc_budget is not None:
        shots = affordable_shot_count(
            config=config,
            requested_shots=shots,
            remaining_hqc=remaining_hqc,
            settings=settings,
        )
    if shots <= 0:
        return 0, remaining_measurements, remaining_hqc

    spent = spend_measurements(
        backend=rmb.backend,
        rng=rmb.rng,
        data=rmb._data,
        config=config,
        n_measurements=shots,
    )
    remaining_measurements -= spent
    if settings.hqc_budget is not None:
        remaining_hqc -= hqc_cost(config, spent, settings)
    else:
        remaining_hqc = float(remaining_measurements)
    return spent, remaining_measurements, remaining_hqc


def dense_reference_run(
    *,
    settings: ContourFirstExperimentConfig,
    seed: int,
    args: argparse.Namespace,
) -> RMB:
    """
    Build a high-quality reference from fixed-ratio depth bisections.

    This avoids using the same adaptive contour-following stopping rule as the
    candidate method.  The resulting reference is still fit with the same
    monotone surface, but the data are deliberately spread along the boundary.
    """
    rng = default_rng(seed)
    rmb = RMB.default(rng).with_backend(make_backend())
    template = template_config(settings, args.n_qubits)

    remaining_measurements = settings.measurement_budget
    remaining_hqc = (
        float(settings.hqc_budget)
        if settings.hqc_budget is not None
        else float(settings.measurement_budget)
    )

    def measure(depth: float, ratio: float, shots: int) -> tuple[object, float | None]:
        nonlocal remaining_measurements, remaining_hqc
        config = config_from_parameters(template=template, depth=depth, ratio=ratio)
        spent, remaining_measurements, remaining_hqc = spend_for_reference(
            rmb=rmb,
            config=config,
            requested_shots=shots,
            settings=settings,
            remaining_measurements=remaining_measurements,
            remaining_hqc=remaining_hqc,
        )
        if spent <= 0 and config not in rmb._data:
            return config, None
        return config, fidelity_mean(rmb._data[config])

    ratios = np.linspace(args.ratio_min, args.ratio_max, args.reference_ratios)
    for ratio in ratios:
        if remaining_measurements <= 0 or remaining_hqc <= settings.hqc_base_cost:
            break

        low_depth = float(args.depth_min)
        high_depth = float(args.depth_max)
        low_config, low_p = measure(low_depth, float(ratio), args.reference_edge_shots)
        high_config, high_p = measure(high_depth, float(ratio), args.reference_edge_shots)
        if low_p is None or high_p is None:
            break

        if low_p < 0.5 or high_p > 0.5:
            continue

        best_config = low_config if abs(low_p - 0.5) <= abs(high_p - 0.5) else high_config
        best_error = min(abs(low_p - 0.5), abs(high_p - 0.5))

        for _ in range(args.reference_bisection_steps):
            if remaining_measurements <= 0 or remaining_hqc <= settings.hqc_base_cost:
                break
            mid_depth = 0.5 * (low_depth + high_depth)
            mid_config, mid_p = measure(mid_depth, float(ratio), args.reference_shots)
            if mid_p is None:
                break
            error = abs(mid_p - 0.5)
            if error < best_error:
                best_config = mid_config
                best_error = error
            if mid_p >= 0.5:
                low_depth = float(mid_config.depth)
                low_p = mid_p
            else:
                high_depth = float(mid_config.depth)
                high_p = mid_p

        if args.reference_confirm_shots > 0:
            measure(
                float(best_config.depth),
                float(best_config.min_two_qubit_gate_ratio),
                args.reference_confirm_shots,
            )

    return rmb


def cache_path(
    *,
    label: str,
    budget: float,
    seed: int,
    args: argparse.Namespace,
) -> Path:
    cache_dir = Path(args.cache_dir).expanduser().resolve()
    reference_suffix = ""
    if label == "true":
        reference_suffix = (
            f"_{args.true_method}"
            f"_rr{args.reference_ratios}"
            f"_rs{args.reference_shots}"
            f"_rb{args.reference_bisection_steps}"
            f"_rc{args.reference_confirm_shots}"
        )
    name = (
        f"v5_{label}_"
        f"{args.budget_mode}{budget:g}_"
        f"{args.n_qubits}q_"
        f"d{args.depth_min}-{args.depth_max}_"
        f"r{args.ratio_min:g}-{args.ratio_max:g}_"
        f"seed{seed}{reference_suffix}_trace"
        f"{'local' if args.trace_local_crossing else 'projected'}"
        f"{'mono' if not args.no_trace_monotone_depth else 'free'}"
        f"_dc{args.crossing_decision_confirm_shots}"
        f"w{args.crossing_decision_probability_width:g}"
        f"{'af' if not args.no_crossing_ambiguous_as_failure else 'ap'}"
        f"_ib{'y' if not args.no_initial_crossing_interior_bracket else 'n'}"
        f"c{args.initial_crossing_center_fraction:g}"
        f"h{args.initial_crossing_half_width_fraction:g}"
        f"_ia{'r' if not args.no_initial_anchor_refine else 'n'}"
        f"{args.initial_anchor_refine_shots}"
        f"m{args.initial_anchor_min_runs}"
        f"_tr{args.trace_reserve_budget_fraction:g}m{args.trace_reserve_min_hqc:g}"
        f"_tc{args.trace_correction_steps}rw{args.trace_reject_probability_width:g}"
        f"_ta{args.trace_anchor_min_shots}"
        f"_cc{args.crossing_confirm_candidates}s{args.crossing_confirm_shots}7.json"
    )
    return cache_dir / name


def requested_budget_value(settings: ContourFirstExperimentConfig) -> float:
    if settings.hqc_budget is not None:
        return float(settings.hqc_budget)
    return float(settings.measurement_budget)


def spent_budget_value(rmb: RMB, settings: ContourFirstExperimentConfig) -> float:
    if settings.hqc_budget is not None:
        return hqc_spent(rmb, settings)
    return float(total_measurements(rmb._data))


def print_budget_warning(
    *,
    label: str,
    rmb: RMB,
    settings: ContourFirstExperimentConfig,
    args: argparse.Namespace,
) -> None:
    requested = requested_budget_value(settings)
    spent = spent_budget_value(rmb, settings)
    if requested <= 0.0 or not np.isfinite(spent):
        return
    fraction = spent / requested
    if fraction < args.warn_spend_fraction:
        print(
            f"WARNING: {label} spent only {spent:.1f} / {requested:.1f} "
            f"({fraction:.1%}) of the requested "
            f"{'HQC' if settings.hqc_budget is not None else 'measurement'} budget."
        )


def run_or_load(
    *,
    label: str,
    budget: float,
    seed: int,
    args: argparse.Namespace,
) -> tuple[RMB, ContourFirstExperimentConfig]:
    settings = make_settings(budget=budget, seed=seed, args=args)
    path = cache_path(label=label, budget=budget, seed=seed, args=args)
    if path.exists() and not args.rebuild:
        print(f"Loading {label} budget {budget:g} from {path}")
        rmb = RMB.load(path)
        print_budget_warning(label=label, rmb=rmb, settings=settings, args=args)
        return rmb, settings

    print(f"Running {label} budget {budget:g}...")
    if label == "true" and args.true_method == "dense":
        rmb = dense_reference_run(settings=settings, seed=seed, args=args)
    else:
        rmb = estimate_boundary_v5(settings)
    path.parent.mkdir(parents=True, exist_ok=True)
    rmb.save(path)
    print(f"Saved {label} budget {budget:g} to {path}")
    print_budget_warning(label=label, rmb=rmb, settings=settings, args=args)
    return rmb, settings


def extract_contour(
    rmb: RMB,
    settings: ContourFirstExperimentConfig,
) -> tuple[np.ndarray, object | None]:
    try:
        surface = fit_monotone_fidelity_surface(rmb._data, settings)
    except (RuntimeError, ValueError):
        return np.empty((0, 2), dtype=float), None
    return contour_points_from_surface(surface, settings), surface


def hqc_spent(rmb: RMB, settings: ContourFirstExperimentConfig) -> float:
    if settings.hqc_budget is None:
        return float("nan")
    return float(sum(
        hqc_cost(config, estimator.num_runs(), settings)
        for config, estimator in rmb._data.items()
        if estimator.num_runs() > 0
    ))


def sorted_contour(contour: np.ndarray) -> np.ndarray:
    if len(contour) == 0:
        return contour
    order = np.lexsort((contour[:, 0], contour[:, 1]))
    return contour[order]


def point_arrays(rmb: RMB) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    configs = [config for config, estimator in rmb._data.items() if estimator.num_runs() > 0]
    depths = np.array([config.depth for config in configs], dtype=float)
    ratios = np.array([config.min_two_qubit_gate_ratio for config in configs], dtype=float)
    fidelities = np.array([fidelity_mean(rmb._data[config]) for config in configs], dtype=float)
    variances = np.array([fidelity_variance(rmb._data[config]) for config in configs], dtype=float)
    certainty = 1.0 - np.clip(variances / (1.0 / 12.0), 0.0, 1.0)
    sizes = 25.0 + 140.0 * certainty
    return depths, ratios, fidelities, sizes


def plot_surface_panel(
    ax,
    *,
    rmb: RMB,
    settings: ContourFirstExperimentConfig,
    budget_label: str,
    true_contour: np.ndarray,
    show_points: bool = True,
) -> None:
    import matplotlib.pyplot as plt

    contour, surface = extract_contour(rmb, settings)
    cmap = fidelity_colormap()
    scatter = None
    if show_points:
        depths, ratios, fidelities, sizes = point_arrays(rmb)
        scatter = ax.scatter(
            depths,
            ratios,
            c=fidelities,
            s=sizes,
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            edgecolors="black",
            linewidths=0.35,
            zorder=3,
        )

    if surface is not None:
        depth_grid, ratio_grid, probabilities = surface.probability_grid(settings)
        ax.contourf(
            depth_grid,
            ratio_grid,
            probabilities,
            levels=np.linspace(0.0, 1.0, 16),
            cmap=cmap,
            alpha=0.25,
            vmin=0.0,
            vmax=1.0,
            zorder=1,
        )
        if len(contour) > 0:
            c = sorted_contour(contour)
            ax.plot(c[:, 0], c[:, 1], color="black", linewidth=2.0, label="estimate", zorder=4)

    if len(true_contour) > 0:
        t = sorted_contour(true_contour)
        ax.plot(t[:, 0], t[:, 1], color="white", linewidth=3.0, alpha=0.8, zorder=4)
        ax.plot(t[:, 0], t[:, 1], color="red", linewidth=1.5, label="true high-budget", zorder=5)

    used = total_measurements(rmb._data)
    hqc = hqc_spent(rmb, settings)
    if np.isfinite(hqc):
        requested = requested_budget_value(settings)
        hqc_text = f", HQC {hqc:.1f}/{requested:.0f}"
    else:
        hqc_text = ""
    ax.set_title(f"{budget_label}\n{used} outcomes, {len(rmb._data)} configs{hqc_text}")
    ax.set_xlabel("# Gates")
    ax.set_ylabel("Two-qubit gate ratio")
    ax.set_xlim(settings.depth_bounds)
    ax.set_ylim(settings.ratio_bounds)
    ax.legend(loc="best", fontsize=8)
    return scatter


def plot_overlay_panel(
    ax,
    *,
    contours: list[tuple[str, np.ndarray]],
    true_label: str,
    true_contour: np.ndarray,
    settings: ContourFirstExperimentConfig,
) -> None:
    cmap = __import__("matplotlib").colormaps.get_cmap("plasma")
    plotted = 0
    for idx, (label, contour) in enumerate(contours):
        if len(contour) == 0:
            continue
        color = cmap(idx / max(1, len(contours) - 1))
        c = sorted_contour(contour)
        ax.plot(c[:, 0], c[:, 1], linewidth=1.8, color=color, label=label)
        plotted += 1

    if len(true_contour) > 0:
        t = sorted_contour(true_contour)
        ax.plot(t[:, 0], t[:, 1], color="black", linewidth=3.0, label=true_label)

    ax.set_title("Fidelity = 0.5 contour comparison")
    ax.set_xlabel("# Gates")
    ax.set_ylabel("Two-qubit gate ratio")
    ax.set_xlim(settings.depth_bounds)
    ax.set_ylim(settings.ratio_bounds)
    if plotted or len(true_contour) > 0:
        ax.legend(loc="best", fontsize=8)


def make_plots(
    *,
    runs: list[tuple[str, RMB, ContourFirstExperimentConfig]],
    true_run: tuple[str, RMB, ContourFirstExperimentConfig],
    args: argparse.Namespace,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    true_label, true_rmb, true_settings = true_run
    true_contour, _ = extract_contour(true_rmb, true_settings)

    surface_runs = runs + [true_run]
    n_panels = len(surface_runs) + 1
    n_cols = min(args.columns, n_panels)
    n_rows = ceil(n_panels / n_cols)
    fig = plt.figure(figsize=(5.2 * n_cols + 0.8, 4.4 * n_rows), constrained_layout=True)
    grid = GridSpec(
        n_rows,
        n_cols + 1,
        figure=fig,
        width_ratios=[1.0] * n_cols + [0.06],
    )
    flat_axes = [
        fig.add_subplot(grid[row, col])
        for row in range(n_rows)
        for col in range(n_cols)
    ]
    cbar_ax = fig.add_subplot(grid[:, -1])

    scatter = None
    contours = []
    for ax, (label, rmb, settings) in zip(flat_axes, surface_runs):
        panel_scatter = plot_surface_panel(
            ax,
            rmb=rmb,
            settings=settings,
            budget_label=label,
            true_contour=true_contour,
            show_points=(label != true_label or args.show_reference_points),
        )
        if panel_scatter is not None:
            scatter = panel_scatter
        contour, _ = extract_contour(rmb, settings)
        if label != true_label:
            contours.append((label, contour))

    overlay_ax = flat_axes[len(surface_runs)]
    plot_overlay_panel(
        overlay_ax,
        contours=contours,
        true_label=true_label,
        true_contour=true_contour,
        settings=true_settings,
    )

    for ax in flat_axes[n_panels:]:
        ax.axis("off")

    if scatter is not None:
        cbar = fig.colorbar(scatter, cax=cbar_ax)
        cbar.set_label("Posterior mean fidelity")
    else:
        cbar_ax.axis("off")

    if args.output is not None:
        output = Path(args.output).expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists() and args.overwrite:
            output.unlink()
            print(f"Overwriting plot at {output}")
        fig.savefig(output, dpi=args.dpi)
        print(f"Saved plot to {output}")
    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot viarregio5 contour estimates across budgets against a high-budget reference."
    )
    parser.add_argument("--budgets", type=str, default="100,250,500,1000")
    parser.add_argument("--true-budget", type=float, default=5000.0)
    parser.add_argument(
        "--true-method",
        choices=("dense", "adaptive"),
        default="dense",
        help="How to build the high-budget reference contour.",
    )
    parser.add_argument("--reference-ratios", type=int, default=28)
    parser.add_argument("--reference-bisection-steps", type=int, default=7)
    parser.add_argument("--reference-edge-shots", type=int, default=4)
    parser.add_argument("--reference-shots", type=int, default=6)
    parser.add_argument("--reference-confirm-shots", type=int, default=10)
    parser.add_argument("--budget-mode", choices=("hqc", "measurements"), default="hqc")
    parser.add_argument("--hqc-budget", type=float, default=None)
    parser.add_argument("--measurement-cap", type=int, default=1_000_000_000)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--true-seed-offset", type=int, default=100_000)
    parser.add_argument("--n-qubits", type=int, default=10)
    parser.add_argument("--depth-min", type=int, default=4)
    parser.add_argument("--depth-max", type=int, default=300)
    parser.add_argument("--ratio-min", type=float, default=0.08)
    parser.add_argument("--ratio-max", type=float, default=0.8)
    parser.add_argument("--random-elimination", type=float, default=0.1)
    parser.add_argument("--scrambling-probability", type=float, default=0.0)
    parser.add_argument("--contour-grid", type=int, default=100)
    parser.add_argument("--surface-smoothing", type=float, default=0.1)
    parser.add_argument("--max-adaptive-shots", type=int, default=5)
    parser.add_argument("--max-shots-per-config", type=int, default=10)
    parser.add_argument("--ray-ratio-count", type=int, default=5)
    parser.add_argument("--crossing-decision-confirm-shots", type=int, default=4)
    parser.add_argument("--crossing-decision-probability-width", type=float, default=0.20)
    parser.add_argument(
        "--no-crossing-ambiguous-as-failure",
        action="store_true",
        help="Let ambiguous initial bisection points update the pass side if their mean is above 0.5.",
    )
    parser.add_argument(
        "--no-initial-crossing-interior-bracket",
        action="store_true",
        help="Start initial crossing searches at min/max depth instead of an interior bracket.",
    )
    parser.add_argument("--initial-crossing-center-fraction", type=float, default=0.35)
    parser.add_argument("--initial-crossing-half-width-fraction", type=float, default=0.18)
    parser.add_argument("--initial-crossing-expand-factor", type=float, default=1.6)
    parser.add_argument("--crossing-confirm-candidates", type=int, default=0)
    parser.add_argument("--crossing-confirm-shots", type=int, default=4)
    parser.add_argument("--no-initial-anchor-refine", action="store_true")
    parser.add_argument("--initial-anchor-refine-shots", type=int, default=6)
    parser.add_argument("--initial-anchor-refine-steps", type=int, default=3)
    parser.add_argument("--initial-anchor-search-fraction", type=float, default=0.06)
    parser.add_argument("--initial-anchor-min-runs", type=int, default=8)
    parser.add_argument("--trace-reserve-budget-fraction", type=float, default=0.35)
    parser.add_argument("--trace-reserve-min-hqc", type=float, default=25.0)
    parser.add_argument("--trace-ratio-step-fraction", type=float, default=0.06)
    parser.add_argument(
        "--trace-local-crossing",
        action="store_true",
        help="Use local fixed-ratio depth crossings for trace steps. More principled, but costly at low budgets.",
    )
    parser.add_argument(
        "--no-trace-monotone-depth",
        action="store_true",
        help="Allow the traced contour depth to increase when the two-qubit ratio increases.",
    )
    parser.add_argument("--trace-correction-steps", type=int, default=2)
    parser.add_argument("--trace-reject-probability-width", type=float, default=0.25)
    parser.add_argument("--trace-anchor-min-shots", type=int, default=8)
    parser.add_argument("--refinement-shots", type=int, default=2)
    parser.add_argument("--diagnostics-dir", type=str, default=None)
    parser.add_argument("--print-diagnostics", action="store_true")
    parser.add_argument(
        "--cache-dir",
        type=str,
        default=str(Path(__file__).resolve().parent / "benchmark_cache"),
    )
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Rebuild all cached benchmark runs and overwrite the output plot.",
    )
    parser.add_argument("--output", type=str, default="scripts/experiments/rmb_bayes/benchmarkv5.png")
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--columns", type=int, default=3)
    parser.add_argument(
        "--show-reference-points",
        dest="show_reference_points",
        action="store_true",
        help="Show dense-reference probe points as well as the fitted reference contour.",
    )
    parser.add_argument(
        "--hide-reference-points",
        dest="show_reference_points",
        action="store_false",
        help="Hide dense-reference probe points and show only the fitted reference contour.",
    )
    parser.set_defaults(show_reference_points=False)
    parser.add_argument("--no-show", action="store_true", help="Save the figure without opening a window.")
    parser.add_argument("--warn-spend-fraction", type=float, default=0.8)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main(**overrides) -> None:
    args = parse_args()
    for name, value in overrides.items():
        if not hasattr(args, name):
            raise ValueError(f"Unknown benchmarkv5 option: {name}")
        setattr(args, name, value)
    if args.overwrite:
        args.rebuild = True
    budgets = parse_budget_list(args.budgets)

    true_seed = args.seed + args.true_seed_offset
    true_run = (
        f"true {args.true_method} {args.budget_mode} {args.true_budget:g}",
        *run_or_load(label="true", budget=args.true_budget, seed=true_seed, args=args),
    )

    runs = []
    for budget in budgets:
        label = f"{args.budget_mode} {budget:g}"
        rmb, settings = run_or_load(label="budget", budget=budget, seed=args.seed, args=args)
        runs.append((label, rmb, settings))

    make_plots(runs=runs, true_run=true_run, args=args)


if __name__ == "__main__":
    main(overwrite=True)
