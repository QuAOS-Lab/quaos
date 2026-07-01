"""
Benchmark contour-finding methods with the S1/S2 log-contour scores.

This is a scoring-focused variant of ``benchmarkbenchmark.py``.  It runs the
same configured approaches and noise realisations, but scores the fitted
fidelity-0.5 line using the conventions in ``scores.py``:

* S1: integrated absolute error in log gate count, normalized by contour area.
* S2: integrated fitted-contour uncertainty in log gate count, normalized by
  contour area.

Lower S1 and S2 are better.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from statistics import mean

_SRC_ROOT = Path(__file__).resolve().parents[4]
sys.path = [path for path in sys.path if path != str(_SRC_ROOT)]
sys.path.insert(0, str(_SRC_ROOT))

import matplotlib.pyplot as plt
import numpy as np

from sympleq.applications.randomized_benchmarking.experiments.benchmarkbenchmark import (
    APPROACHES,
    HQC_BUDGET,
    MAX_COST_PER_RUN,
    MAX_WORKERS,
    N_BOOTSTRAP,
    N_NOISE_VALUES,
    N_REALISATIONS,
    NOISE_SPREAD,
    approach_settings,
    noise_values,
    run_approach,
)
from sympleq.applications.randomized_benchmarking.experiments.common import (
    CrossingSettings,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    boundary_bootstrap_for_settings,
    boundary_fit_for_settings,
    gate_counts_from_boundary_fit,
    parametric_boundary_gate_samples,
    plot_uncertainty_diagnostics,
)
from sympleq.applications.randomized_benchmarking.experiments.scores import (
    integrate_trapezoid,
    true_log_gates,
)
from sympleq.applications.randomized_benchmarking.experiments.score_summary_histogram import (
    plot_score_summaries,
)


SCORE_RATIO_POINTS = int(os.environ.get("BENCHMARKBENCHMARK_SCORE_RATIO_POINTS", "300"))
SAVE_PLOTS = os.environ.get("BENCHMARKBENCHMARK_SCORE_PLOTS", "1") != "0"
SHOW_SUMMARY_PLOTS = os.environ.get("BENCHMARKBENCHMARK_SCORE_SHOW", "1") != "0"
FIG_DIR = Path(__file__).resolve().parent / "figs" / "benchmarkbenchmark_scores"
RESULTS_JSON = FIG_DIR / "results.json"
RUNS_CSV = FIG_DIR / "runs.csv"


def finite_nanstd_columns(samples: np.ndarray) -> np.ndarray:
    """Column-wise standard deviation, ignoring NaNs without all-NaN warnings."""
    output = np.full(samples.shape[1], np.nan, dtype=float)
    for column_index in range(samples.shape[1]):
        column = samples[:, column_index]
        column = column[np.isfinite(column)]
        if len(column) > 1:
            output[column_index] = float(np.std(column, ddof=1))
        elif len(column) == 1:
            output[column_index] = 0.0
    return output


def contour_scores(
    data,
    settings: CrossingSettings,
    *,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
    seed: int,
) -> dict[str, float | int | None]:
    """Score an approach-provided boundary fit with the S1/S2 convention."""
    ratios = np.linspace(
        settings.ratio_bounds[0],
        settings.ratio_bounds[1],
        SCORE_RATIO_POINTS,
    )
    fit = boundary_fit_for_settings(data, settings)
    if fit is None:
        return unavailable_scores()

    gates = gate_counts_from_boundary_fit(fit, ratios)
    log_gates = np.full_like(gates, np.nan, dtype=float)
    valid_gates = np.isfinite(gates) & (gates > 0.0)
    log_gates[valid_gates] = np.log(gates[valid_gates])
    analytic_log_gates = true_log_gates(
        ratios,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    valid = (
        np.isfinite(ratios)
        & np.isfinite(log_gates)
        & np.isfinite(analytic_log_gates)
    )
    if np.count_nonzero(valid) < 2:
        return unavailable_scores()

    scored_ratios = ratios[valid]
    scored_log_gates = log_gates[valid]
    scored_analytic = analytic_log_gates[valid]
    area = float(integrate_trapezoid(scored_log_gates, scored_ratios))
    if not np.isfinite(area) or abs(area) <= 1e-12:
        return unavailable_scores()

    delta = np.abs(scored_log_gates - scored_analytic)
    s1 = float(integrate_trapezoid(delta, scored_ratios) / abs(area))

    fits = boundary_bootstrap_for_settings(
        data,
        settings,
        n_bootstrap=N_BOOTSTRAP,
        seed=seed,
    )
    gate_samples = parametric_boundary_gate_samples(fits, ratios)
    s2: float | None = None
    mean_sigma: float | None = None
    if len(gate_samples) > 0:
        log_samples = np.full_like(gate_samples, np.nan, dtype=float)
        valid_samples = np.isfinite(gate_samples) & (gate_samples > 0.0)
        log_samples[valid_samples] = np.log(gate_samples[valid_samples])
        sigma = finite_nanstd_columns(log_samples)
        valid_sigma = valid & np.isfinite(sigma)
        if np.count_nonzero(valid_sigma) >= 2:
            s2 = float(
                integrate_trapezoid(sigma[valid_sigma], ratios[valid_sigma])
                / abs(area)
            )
            mean_sigma = float(np.mean(sigma[valid_sigma]))

    return {
        "S1": s1,
        "S2": s2,
        "S_total": None if s2 is None else s1 + s2,
        "A_fit": area,
        "mean_delta_log_gates": float(np.mean(delta)),
        "mean_sigma_contour": mean_sigma,
        "n_contour_points": int(np.count_nonzero(valid)),
    }


def unavailable_scores() -> dict[str, float | int | None]:
    return {
        "S1": None,
        "S2": None,
        "S_total": None,
        "A_fit": None,
        "mean_delta_log_gates": None,
        "mean_sigma_contour": None,
        "n_contour_points": 0,
    }


def format_optional(value: float | int | None, precision: int = 6) -> str:
    if value is None:
        return "unavailable"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{precision}g}"


def print_row(row: dict) -> None:
    print(
        f"approach={row['approach']} "
        f"noise={row['noise_index']:02d} realisation={row['realisation']:02d}: "
        f"S1={format_optional(row['S1'])}, "
        f"S2={format_optional(row['S2'])}, "
        f"spent={row['spent_hqc']:.1f} HQC, "
        f"crossings={row['n_crossings']}"
    )


def run_one(
    *,
    approach: str,
    noise_index: int,
    realisation: int,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> dict:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    seed = 10_000 * noise_index + realisation
    settings = approach_settings(
        approach,
        seed=seed,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )

    rmb, crossings, budget = run_approach(approach, settings)
    scores = contour_scores(
        rmb._data,
        settings,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
        seed=seed + 500_000,
    )

    out_path = None
    if SAVE_PLOTS:
        out_path = (
            FIG_DIR
            / f"{approach}_noise_{noise_index:02d}_realisation_{realisation:02d}.png"
        )
        plot_uncertainty_diagnostics(
            rmb._data,
            settings,
            png_path=out_path,
            show=False,
        )
        plt.close("all")

    return {
        "approach": approach,
        "noise_index": noise_index,
        "realisation": realisation,
        "seed": seed,
        **scores,
        "spent_hqc": budget.spent_hqc,
        "n_crossings": len(crossings),
        "one_q_noise_scale": one_q_noise_scale,
        "two_q_noise_scale": two_q_noise_scale,
        "figure": "" if out_path is None else str(out_path),
    }


def benchmark_tasks() -> list[dict]:
    tasks = []
    for noise_index, (one_q_scale, two_q_scale) in enumerate(noise_values()):
        for approach in APPROACHES:
            for realisation in range(N_REALISATIONS):
                tasks.append({
                    "approach": approach,
                    "noise_index": noise_index,
                    "realisation": realisation,
                    "one_q_noise_scale": one_q_scale,
                    "two_q_noise_scale": two_q_scale,
                })
    return tasks


def summarize_optional(values: list[float | None]) -> dict[str, float | None]:
    finite = np.asarray([value for value in values if value is not None], dtype=float)
    if len(finite) == 0:
        return {
            "average": None,
            "variance": None,
            "best": None,
            "worst": None,
        }
    return {
        "average": float(np.mean(finite)),
        "variance": float(np.var(finite)),
        "best": float(np.min(finite)),
        "worst": float(np.max(finite)),
    }


def mean_optional(values: list[float | None]) -> float | None:
    finite = [value for value in values if value is not None]
    return None if not finite else float(mean(finite))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    tasks = benchmark_tasks()
    max_workers = (
        min(os.process_cpu_count() or 1, max(1, len(tasks)))
        if MAX_WORKERS is None
        else min(MAX_WORKERS, max(1, len(tasks)))
    )
    print(
        f"Running {len(APPROACHES)} approaches, {N_NOISE_VALUES} noise values, "
        f"{N_REALISATIONS} realisations per noise value with {max_workers} workers "
        f"(budget {HQC_BUDGET:.1f} HQC, max cost/run {MAX_COST_PER_RUN:.1f} HQC)"
    )
    print("Scoring with S1/S2 log-contour metrics from scores.py; lower is better.")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_one, **task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print_row(row)

    rows.sort(key=lambda row: (row["approach"], row["noise_index"], row["realisation"]))

    by_approach = {}
    for approach in APPROACHES:
        approach_rows = [row for row in rows if row["approach"] == approach]
        spent = [row["spent_hqc"] for row in approach_rows]
        by_approach[approach] = {
            "n_runs": len(approach_rows),
            "S1": summarize_optional([row["S1"] for row in approach_rows]),
            "S2": summarize_optional([row["S2"] for row in approach_rows]),
            "S_total": summarize_optional([row["S_total"] for row in approach_rows]),
            "mean_delta_log_gates": mean_optional(
                [row["mean_delta_log_gates"] for row in approach_rows]
            ),
            "mean_sigma_contour": mean_optional(
                [row["mean_sigma_contour"] for row in approach_rows]
            ),
            "average_spent_hqc": mean(spent) if spent else None,
        }

    summary = {
        "n_noise_values": N_NOISE_VALUES,
        "n_realisations": N_REALISATIONS,
        "approaches": list(APPROACHES),
        "n_runs": len(rows),
        "noise_spread": NOISE_SPREAD,
        "n_bootstrap": N_BOOTSTRAP,
        "score_ratio_points": SCORE_RATIO_POINTS,
        "hqc_budget": HQC_BUDGET,
        "max_cost_per_run": MAX_COST_PER_RUN,
        "save_plots": SAVE_PLOTS,
        "show_summary_plots": SHOW_SUMMARY_PLOTS,
        "by_approach": by_approach,
    }
    RESULTS_JSON.write_text(
        json.dumps({"summary": summary, "runs": rows}, indent=2),
        encoding="utf-8",
    )
    fieldnames = [
        "approach",
        "noise_index",
        "realisation",
        "seed",
        "S1",
        "S2",
        "S_total",
        "A_fit",
        "mean_delta_log_gates",
        "mean_sigma_contour",
        "n_contour_points",
        "spent_hqc",
        "n_crossings",
        "one_q_noise_scale",
        "two_q_noise_scale",
        "figure",
    ]
    write_csv(RUNS_CSV, rows, fieldnames)
    score_plot_paths = plot_score_summaries(
        RESULTS_JSON,
        output_dir=FIG_DIR,
        show=SHOW_SUMMARY_PLOTS,
    )

    print("\nSummary")
    for approach, approach_summary in summary["by_approach"].items():
        print(f"  {approach}")
        for score_name in ("S1", "S2", "S_total"):
            score_summary = approach_summary[score_name]
            print(
                f"    {score_name}: "
                f"avg={format_optional(score_summary['average'])}, "
                f"var={format_optional(score_summary['variance'])}, "
                f"best={format_optional(score_summary['best'])}, "
                f"worst={format_optional(score_summary['worst'])}"
            )
        print(
            "    mean delta log gates: "
            f"{format_optional(approach_summary['mean_delta_log_gates'])}"
        )
        print(
            "    mean sigma contour: "
            f"{format_optional(approach_summary['mean_sigma_contour'])}"
        )
        spent = approach_summary["average_spent_hqc"]
        print(f"    average spent: {format_optional(spent)} HQC")
    print("  score plots:")
    for name, path in score_plot_paths.items():
        print(f"    {name}: {path}")
    print(f"  wrote results to {FIG_DIR}")


if __name__ == "__main__":
    main()
