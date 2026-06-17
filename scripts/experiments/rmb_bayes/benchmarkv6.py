from __future__ import annotations

import argparse
from pathlib import Path

from sympleq.applications.randomized_benchmarking.RMB import RMB

from benchmarkv5 import (
    dense_reference_run,
    make_plots,
    parse_budget_list,
    print_budget_warning,
)
from viarregio6 import GeometryAnchorExperimentConfig
from viarregio6 import estimate_boundary as estimate_boundary_v6
from viarregio2 import make_backend


def make_settings(
    *,
    budget: float,
    seed: int,
    args: argparse.Namespace,
) -> GeometryAnchorExperimentConfig:
    if args.budget_mode == "hqc":
        measurement_budget = args.measurement_cap
        hqc_budget = float(budget)
    else:
        measurement_budget = int(round(budget))
        hqc_budget = args.hqc_budget

    return GeometryAnchorExperimentConfig(
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
        anchor_ratio_count=args.anchor_ratio_count,
        anchor_max_ratio_fraction=args.anchor_max_ratio_fraction,
        anchor_probe_shots=args.anchor_probe_shots,
        anchor_probe_bisection_steps=args.anchor_probe_bisection_steps,
        first_anchor_full_search=not args.no_first_anchor_full_search,
        anchor_accept_score=args.anchor_accept_score,
        anchor_min_probe_candidates=args.anchor_min_probe_candidates,
        anchor_refine_selected=args.anchor_refine_selected,
        anchor_use_local_prediction=not args.no_anchor_local_prediction,
        anchor_depth_margin_weight=args.anchor_depth_margin_weight,
        anchor_ratio_cost_weight=args.anchor_ratio_cost_weight,
        anchor_boundary_weight=args.anchor_boundary_weight,
        anchor_slope_weight=args.anchor_slope_weight,
        anchor_target_scaled_slope=args.anchor_target_scaled_slope,
        ray_ratio_count=args.ray_ratio_count,
        ray_probe_shots=2,
        ray_bisection_steps=5,
        ray_bisection_shots=2,
        trace_ratio_step_fraction=args.trace_ratio_step_fraction,
        trace_depth_search_fraction=0.06,
        trace_local_crossing=args.trace_local_crossing,
        trace_correction_steps=4,
        trace_shots=2,
        trace_accept_probability_width=0.12,
        trace_directions=tuple(args.trace_directions),
        model_projection_after_fit=True,
        refine_after_trace=True,
        refinement_shots=args.refinement_shots,
        refinement_boundary_width=0.15,
        hqc_cost_informed_acquisition=(args.budget_mode == "hqc" or hqc_budget is not None),
        hqc_cost_power=1.0,
        rng_seed=seed,
        save_path=None,
        verbose=args.verbose,
    )


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
    trace = "local" if args.trace_local_crossing else "projected"
    directions = "".join("m" if direction < 0 else "p" for direction in args.trace_directions)
    name = (
        f"v6_{label}_"
        f"{args.budget_mode}{budget:g}_"
        f"{args.n_qubits}q_"
        f"d{args.depth_min}-{args.depth_max}_"
        f"r{args.ratio_min:g}-{args.ratio_max:g}_"
        f"seed{seed}{reference_suffix}_"
        f"a{args.anchor_ratio_count}_"
        f"am{args.anchor_max_ratio_fraction:g}_"
        f"apb{args.anchor_probe_bisection_steps}_"
        f"{'ffs' if not args.no_first_anchor_full_search else 'fps'}_"
        f"aa{args.anchor_accept_score:g}_"
        f"{'arl' if not args.no_anchor_local_prediction else 'arg'}_"
        f"trace{trace}_{directions}.json"
    )
    return cache_dir / name


def run_or_load(
    *,
    label: str,
    budget: float,
    seed: int,
    args: argparse.Namespace,
) -> tuple[RMB, GeometryAnchorExperimentConfig]:
    settings = make_settings(budget=budget, seed=seed, args=args)
    path = cache_path(label=label, budget=budget, seed=seed, args=args)
    if path.exists() and not args.rebuild:
        print(f"Loading {label} budget {budget:g} from {path}")
        rmb = RMB.load(path).with_backend(make_backend())
        print_budget_warning(label=label, rmb=rmb, settings=settings, args=args)
        return rmb, settings

    print(f"Running {label} budget {budget:g}...")
    if label == "true" and args.true_method == "dense":
        rmb = dense_reference_run(settings=settings, seed=seed, args=args)
    else:
        rmb = estimate_boundary_v6(settings)
    path.parent.mkdir(parents=True, exist_ok=True)
    rmb.save(path)
    print(f"Saved {label} budget {budget:g} to {path}")
    print_budget_warning(label=label, rmb=rmb, settings=settings, args=args)
    return rmb, settings


def parse_trace_directions(value: str) -> tuple[int, ...]:
    directions = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        direction = int(part)
        if direction not in (-1, 1):
            raise argparse.ArgumentTypeError("trace directions must be -1 or 1")
        directions.append(direction)
    if not directions:
        raise argparse.ArgumentTypeError("at least one trace direction is required")
    return tuple(directions)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot viarregio6 geometry-anchor contour estimates across budgets."
    )
    parser.add_argument("--budgets", type=str, default="100,250,500,1000")
    parser.add_argument("--true-budget", type=float, default=5000.0)
    parser.add_argument("--true-method", choices=("dense", "adaptive"), default="dense")
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
    parser.add_argument("--anchor-ratio-count", type=int, default=5)
    parser.add_argument("--anchor-max-ratio-fraction", type=float, default=0.45)
    parser.add_argument("--anchor-probe-shots", type=int, default=1)
    parser.add_argument("--anchor-probe-bisection-steps", type=int, default=1)
    parser.add_argument(
        "--no-first-anchor-full-search",
        action="store_true",
        help="Use the cheap anchor probe even for the first ratio. Usually worse at low budgets.",
    )
    parser.add_argument("--anchor-accept-score", type=float, default=1.15)
    parser.add_argument("--anchor-min-probe-candidates", type=int, default=1)
    parser.add_argument("--anchor-refine-selected", action="store_true")
    parser.add_argument("--no-anchor-local-prediction", action="store_true")
    parser.add_argument("--anchor-depth-margin-weight", type=float, default=1.0)
    parser.add_argument("--anchor-ratio-cost-weight", type=float, default=0.35)
    parser.add_argument("--anchor-boundary-weight", type=float, default=0.6)
    parser.add_argument("--anchor-slope-weight", type=float, default=0.8)
    parser.add_argument("--anchor-target-scaled-slope", type=float, default=0.8)
    parser.add_argument("--ray-ratio-count", type=int, default=5)
    parser.add_argument("--trace-ratio-step-fraction", type=float, default=0.06)
    parser.add_argument("--trace-directions", type=parse_trace_directions, default=(-1, 1))
    parser.add_argument("--trace-local-crossing", action="store_true")
    parser.add_argument("--refinement-shots", type=int, default=2)
    parser.add_argument(
        "--cache-dir",
        type=str,
        default=str(Path(__file__).resolve().parent / "benchmark_cache"),
    )
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--output", type=str, default="scripts/experiments/rmb_bayes/benchmarkv6.png")
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--columns", type=int, default=3)
    parser.add_argument(
        "--hide-reference-points",
        dest="show_reference_points",
        action="store_false",
        help="Hide dense-reference probe points and show only the fitted reference contour.",
    )
    parser.set_defaults(show_reference_points=True)
    parser.add_argument("--no-show", action="store_true")
    parser.add_argument("--warn-spend-fraction", type=float, default=0.8)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def main(**overrides) -> None:
    args = parse_args()
    for name, value in overrides.items():
        if not hasattr(args, name):
            raise ValueError(f"Unknown benchmarkv6 option: {name}")
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
