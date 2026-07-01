"""
Benchmark monotone tracing over a sweep of parameter regimes.

This is the monotone-only analogue of ``benchmarkbenchmark.py``. It evaluates
each parameter regime on the same deterministic noise values and realisations,
prints a summary as soon as that regime completes, and reports the best five
regimes at the end. The regimes are a refined second sweep around the best
regions found in the first 50-regime sweep.
"""
from __future__ import annotations

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
from numpy.random import default_rng

from sympleq.applications.randomized_benchmarking.experiments.benchmarkbenchmark import (
    N_REALISATIONS,
    NOISE_SPREAD,
    bootstrap_coverage,
    fit_score,
    noise_values,
    scaled_noise_backend_factory,
    summarize,
    write_csv,
)
from sympleq.applications.randomized_benchmarking.experiments.monotone_tracing import (
    MonotoneTracingSettings,
    run_with_budget,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    plot_uncertainty_diagnostics,
)


N_PARAMETER_SETS = int(os.environ.get("MONOTONE_SWEEP_PARAMETER_SETS", "50"))
MAX_WORKERS = int(os.environ.get("MONOTONE_SWEEP_WORKERS", "25")) or None
SAVE_FIGURES = bool(int(os.environ.get("MONOTONE_SWEEP_SAVE_FIGURES", "0")))

SWEEP_NAME = os.environ.get("MONOTONE_SWEEP_NAME", "monotone_sweep_refined")
FIG_DIR = Path(__file__).resolve().parent / "figs" / SWEEP_NAME
RESULTS_JSON = FIG_DIR / "results.json"
RUNS_CSV = FIG_DIR / "runs.csv"
SUMMARY_CSV = FIG_DIR / "summary.csv"

PARAMETER_FIELDS = [
    "initial_bracket_fraction",
    "high_ratio_gate_bias_power",
    "low_to_high_growth_factor",
    "bracket_half_width_fraction",
    "trace_gate_growth",
    "trace_gate_shrink",
    "bracket_hint",
    "n_gates_resolution",
    "max_hqc_per_ratio",
    "max_shots_per_config",
    "decision_confidence",
    "min_crossing_confidence",
    "extra_midpoint_shots",
    "min_ratio_step",
    "max_ratio_step",
    "initial_refine_hqc",
    "initial_refine_half_width_fraction",
    "initial_refine_extra_shots",
    "initial_refine_decision_confidence",
    "initial_refine_min_crossing_confidence",
    "adaptive_ratio_step",
    "target_hqc_per_ratio",
    "surface_min_configs",
    "monotonic_min_slack_gates",
    "monotonic_slack_fraction",
]


def current_parameters() -> dict:
    settings = MonotoneTracingSettings()
    return {field: getattr(settings, field) for field in PARAMETER_FIELDS}


def clipped_float(value: float, low: float, high: float) -> float:
    return float(np.clip(value, low, high))


def refined_anchor_parameters() -> list[dict]:
    """
    Anchor regimes from the first sweep, normalized to include new fields.

    These cover the current baseline, the best average-score regime, the best
    composite/coverage regime, and the strongest low-tail regimes.
    """
    base = current_parameters()
    anchors = [dict(base)]
    for updates in [
        # First-sweep #48: best average score with lower spend.
        {
            "initial_bracket_fraction": 0.39062521897393,
            "high_ratio_gate_bias_power": 2.0309491699785296,
            "low_to_high_growth_factor": 1.861600925440569,
            "bracket_half_width_fraction": 0.23851132485188098,
            "trace_gate_growth": 2.1050784610211704,
            "trace_gate_shrink": 1.0083045454001958,
            "max_hqc_per_ratio": 46.30457302306524,
            "max_shots_per_config": 9,
            "decision_confidence": 0.8594322413636332,
            "min_crossing_confidence": 0.5655538990407013,
            "initial_refine_hqc": 53.2748480466087,
            "initial_refine_half_width_fraction": 0.3000547541090491,
            "initial_refine_extra_shots": 3,
            "initial_refine_decision_confidence": 0.8789248223691047,
            "initial_refine_min_crossing_confidence": 0.9632439672236701,
            "adaptive_ratio_step": False,
            "target_hqc_per_ratio": 49.99488056551102,
        },
        # First-sweep #32: best composite score and best bootstrap coverage.
        {
            "initial_bracket_fraction": 0.35120211893211256,
            "high_ratio_gate_bias_power": 1.5895314157491474,
            "low_to_high_growth_factor": 1.951831382179539,
            "bracket_half_width_fraction": 0.22965933928225818,
            "trace_gate_growth": 2.598813231512525,
            "trace_gate_shrink": 0.6792152913745838,
            "max_hqc_per_ratio": 34.76250596293123,
            "max_shots_per_config": 7,
            "decision_confidence": 0.6964607539235972,
            "min_crossing_confidence": 0.5620043155511613,
            "initial_refine_hqc": 48.72949747823284,
            "initial_refine_half_width_fraction": 0.3499729753896317,
            "initial_refine_extra_shots": 4,
            "initial_refine_decision_confidence": 0.8970156961050784,
            "initial_refine_min_crossing_confidence": 0.8636327348878575,
            "adaptive_ratio_step": True,
            "target_hqc_per_ratio": 43.755902491515215,
        },
        # First-sweep #23: strong low-tail score with good average.
        {
            "initial_bracket_fraction": 0.23766849169546267,
            "high_ratio_gate_bias_power": 1.3503225212469359,
            "low_to_high_growth_factor": 1.7917462299966245,
            "bracket_half_width_fraction": 0.3056568485559962,
            "trace_gate_growth": 2.1801830852960054,
            "trace_gate_shrink": 0.8826983155055188,
            "max_hqc_per_ratio": 62.06351737691256,
            "max_shots_per_config": 9,
            "decision_confidence": 0.6900413596787146,
            "min_crossing_confidence": 0.5617793297538064,
            "initial_refine_hqc": 32.36695054179802,
            "initial_refine_half_width_fraction": 0.2474593547659775,
            "initial_refine_extra_shots": 7,
            "initial_refine_decision_confidence": 0.9340074448343959,
            "initial_refine_min_crossing_confidence": 0.9174594964007631,
            "adaptive_ratio_step": True,
            "target_hqc_per_ratio": 48.06805664805335,
        },
        # First-sweep #44: best low-tail score, kept to test if coverage can be rescued.
        {
            "initial_bracket_fraction": 0.236663520325022,
            "high_ratio_gate_bias_power": 2.3978458499027493,
            "low_to_high_growth_factor": 2.0190009254627213,
            "bracket_half_width_fraction": 0.35848228317559144,
            "trace_gate_growth": 2.3958927775852,
            "trace_gate_shrink": 0.69020651678864,
            "max_hqc_per_ratio": 63.08578703597753,
            "max_shots_per_config": 10,
            "decision_confidence": 0.8450842608721691,
            "min_crossing_confidence": 0.6374890433933627,
            "initial_refine_hqc": 28.751413566074255,
            "initial_refine_half_width_fraction": 0.1991116709688927,
            "initial_refine_extra_shots": 3,
            "initial_refine_decision_confidence": 0.893080791038469,
            "initial_refine_min_crossing_confidence": 0.9493181179042972,
            "adaptive_ratio_step": True,
            "target_hqc_per_ratio": 47.588026770020534,
        },
    ]:
        regime = dict(base)
        regime.update(updates)
        anchors.append(regime)
    return anchors


def perturb_regime(anchor: dict, rng) -> dict:
    regime = dict(anchor)
    regime.update({
        "initial_bracket_fraction": clipped_float(
            anchor["initial_bracket_fraction"] * rng.uniform(0.85, 1.18),
            0.20,
            0.45,
        ),
        "high_ratio_gate_bias_power": clipped_float(
            anchor["high_ratio_gate_bias_power"] * rng.uniform(0.8, 1.25),
            0.8,
            2.8,
        ),
        "low_to_high_growth_factor": clipped_float(
            anchor["low_to_high_growth_factor"] * rng.uniform(0.9, 1.15),
            1.5,
            2.25,
        ),
        "bracket_half_width_fraction": clipped_float(
            anchor["bracket_half_width_fraction"] * rng.uniform(0.8, 1.25),
            0.18,
            0.42,
        ),
        "trace_gate_growth": clipped_float(
            anchor["trace_gate_growth"] * rng.uniform(0.9, 1.18),
            1.8,
            2.8,
        ),
        "trace_gate_shrink": clipped_float(
            anchor["trace_gate_shrink"] * rng.uniform(0.85, 1.12),
            0.55,
            1.08,
        ),
        "max_hqc_per_ratio": clipped_float(
            anchor["max_hqc_per_ratio"] * rng.uniform(0.85, 1.25),
            32.0,
            75.0,
        ),
        "decision_confidence": clipped_float(
            anchor["decision_confidence"] + rng.uniform(-0.05, 0.05),
            0.68,
            0.88,
        ),
        "min_crossing_confidence": clipped_float(
            anchor["min_crossing_confidence"] + rng.uniform(-0.03, 0.06),
            0.54,
            0.66,
        ),
        "initial_refine_hqc": clipped_float(
            anchor["initial_refine_hqc"] * rng.uniform(0.85, 1.25),
            25.0,
            70.0,
        ),
        "initial_refine_half_width_fraction": clipped_float(
            anchor["initial_refine_half_width_fraction"] * rng.uniform(0.85, 1.25),
            0.18,
            0.42,
        ),
        "initial_refine_decision_confidence": clipped_float(
            anchor["initial_refine_decision_confidence"] + rng.uniform(-0.04, 0.04),
            0.84,
            0.96,
        ),
        "initial_refine_min_crossing_confidence": clipped_float(
            anchor["initial_refine_min_crossing_confidence"] + rng.uniform(-0.08, 0.03),
            0.82,
            0.98,
        ),
        "target_hqc_per_ratio": clipped_float(
            anchor["target_hqc_per_ratio"] * rng.uniform(0.85, 1.2),
            30.0,
            60.0,
        ),
        "n_gates_resolution": float(rng.choice([0.06, 0.08, 0.10])),
        "min_ratio_step": float(rng.choice([0.03, 0.04, 0.05])),
        "max_ratio_step": float(rng.choice([0.14, 0.18, 0.20, 0.24])),
        "surface_min_configs": int(rng.choice([6, 8, 10, 12])),
        "monotonic_min_slack_gates": int(rng.choice([2, 4, 6])),
        "monotonic_slack_fraction": float(rng.choice([0.35, 0.5, 0.65, 0.8])),
    })
    regime["max_shots_per_config"] = int(rng.choice([7, 8, 9, 10]))
    regime["extra_midpoint_shots"] = int(rng.choice([3, 4, 5, 6]))
    regime["initial_refine_extra_shots"] = int(rng.choice([3, 4, 5, 6, 7]))
    regime["adaptive_ratio_step"] = bool(rng.choice([True, True, True, False]))
    regime["bracket_hint"] = str(rng.choice(["line", "both", "surface"]))
    return regime


def parameter_regimes(n: int = N_PARAMETER_SETS) -> list[dict]:
    """
    Deterministic refined regimes around first-sweep winners.

    The first regimes are fixed anchors from the first sweep. The remaining
    regimes are local perturbations, including parameters not explored in the
    first sweep: ``bracket_hint``, ratio-step bounds, gate resolution,
    ``surface_min_configs``, and monotonicity slack.
    """
    anchors = refined_anchor_parameters()
    regimes = anchors[:min(n, len(anchors))]
    rng = default_rng(20260619)

    while len(regimes) < n:
        anchor = anchors[len(regimes) % len(anchors)]
        regimes.append(perturb_regime(anchor, rng))

    return regimes


def run_one(
    *,
    parameter_index: int,
    parameters: dict,
    noise_index: int,
    realisation: int,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> dict:
    seed = 10_000 * noise_index + realisation
    settings = MonotoneTracingSettings(
        rng_seed=seed,
        plot=False,
        verbose=False,
        save_path=None,
        backend_factory=scaled_noise_backend_factory(
            one_q_noise_scale,
            two_q_noise_scale,
        ),
        **parameters,
    )
    object.__setattr__(settings, "one_q_noise_scale", one_q_noise_scale)
    object.__setattr__(settings, "two_q_noise_scale", two_q_noise_scale)

    rmb, crossings, budget = run_with_budget(settings)
    score = fit_score(
        rmb._data,
        settings,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    coverage50, coverage90 = bootstrap_coverage(
        rmb._data,
        settings,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
        seed=seed + 500_000,
    )

    out_path = ""
    if SAVE_FIGURES:
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        png_path = (
            FIG_DIR
            / f"params_{parameter_index:02d}_noise_{noise_index:02d}_realisation_{realisation:02d}.png"
        )
        plot_uncertainty_diagnostics(
            rmb._data,
            settings,
            png_path=png_path,
            show=False,
        )
        out_path = str(png_path)
    plt.close("all")

    return {
        "parameter_index": parameter_index,
        "noise_index": noise_index,
        "realisation": realisation,
        "seed": seed,
        "score": score,
        "bootstrap_50_coverage": coverage50,
        "bootstrap_90_coverage": coverage90,
        "spent_hqc": budget.spent_hqc,
        "n_crossings": len(crossings),
        "one_q_noise_scale": one_q_noise_scale,
        "two_q_noise_scale": two_q_noise_scale,
        "figure": out_path,
        **parameters,
    }


def summarize_rows(parameter_index: int, parameters: dict, rows: list[dict]) -> dict:
    scores = [row["score"] for row in rows]
    coverage50 = [row["bootstrap_50_coverage"] for row in rows]
    coverage90 = [row["bootstrap_90_coverage"] for row in rows]
    spent = [row["spent_hqc"] for row in rows]
    crossings = [row["n_crossings"] for row in rows]
    return {
        "parameter_index": parameter_index,
        "n_runs": len(rows),
        "average_score": summarize(scores)["average"],
        "score_variance": summarize(scores)["variance"],
        "lowest_score": summarize(scores)["lowest"],
        "bootstrap_50_coverage_percent": 100.0 * mean(coverage50),
        "bootstrap_90_coverage_percent": 100.0 * mean(coverage90),
        "average_spent_hqc": mean(spent),
        "average_crossings": mean(crossings),
        **parameters,
    }


def print_summary(summary: dict) -> None:
    print(
        f"params={summary['parameter_index']:02d}: "
        f"avg_score={summary['average_score']:.3f}, "
        f"lowest={summary['lowest_score']:.3f}, "
        f"var={summary['score_variance']:.5f}, "
        f"coverage50={summary['bootstrap_50_coverage_percent']:.1f}%, "
        f"coverage90={summary['bootstrap_90_coverage_percent']:.1f}%, "
        f"spent={summary['average_spent_hqc']:.1f}, "
        f"crossings={summary['average_crossings']:.1f}"
    )


def parameter_tasks(parameter_index: int, parameters: dict) -> list[dict]:
    tasks = []
    for noise_index, (one_q_scale, two_q_scale) in enumerate(noise_values()):
        for realisation in range(N_REALISATIONS):
            tasks.append({
                "parameter_index": parameter_index,
                "parameters": parameters,
                "noise_index": noise_index,
                "realisation": realisation,
                "one_q_noise_scale": one_q_scale,
                "two_q_noise_scale": two_q_scale,
            })
    return tasks


def write_outputs(all_rows: list[dict], summaries: list[dict]) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    all_rows.sort(key=lambda row: (
        row["parameter_index"],
        row["noise_index"],
        row["realisation"],
    ))
    summaries.sort(key=lambda row: row["parameter_index"])
    RESULTS_JSON.write_text(
        json.dumps({
            "summary": {
                "n_parameter_sets": len(summaries),
                "n_noise_values": len(noise_values()),
                "n_realisations": N_REALISATIONS,
                "noise_spread": NOISE_SPREAD,
                "save_figures": SAVE_FIGURES,
            },
            "parameter_summaries": summaries,
            "runs": all_rows,
        }, indent=2),
        encoding="utf-8",
    )
    row_fields = [
        "parameter_index",
        "noise_index",
        "realisation",
        "seed",
        "score",
        "bootstrap_50_coverage",
        "bootstrap_90_coverage",
        "spent_hqc",
        "n_crossings",
        "one_q_noise_scale",
        "two_q_noise_scale",
        "figure",
        *PARAMETER_FIELDS,
    ]
    summary_fields = [
        "parameter_index",
        "n_runs",
        "average_score",
        "score_variance",
        "lowest_score",
        "bootstrap_50_coverage_percent",
        "bootstrap_90_coverage_percent",
        "average_spent_hqc",
        "average_crossings",
        *PARAMETER_FIELDS,
    ]
    write_csv(RUNS_CSV, all_rows, row_fields)
    write_csv(SUMMARY_CSV, summaries, summary_fields)


def main() -> None:
    regimes = parameter_regimes()
    all_rows: list[dict] = []
    summaries: list[dict] = []
    jobs_per_regime = len(noise_values()) * N_REALISATIONS
    max_workers = (
        min(os.process_cpu_count() or 1, jobs_per_regime)
        if MAX_WORKERS is None
        else min(MAX_WORKERS, jobs_per_regime)
    )

    print(
        f"Running {len(regimes)} monotone regimes, {len(noise_values())} noise values, "
        f"{N_REALISATIONS} realisations with {max_workers} workers per regime"
    )

    for parameter_index, parameters in enumerate(regimes):
        tasks = parameter_tasks(parameter_index, parameters)
        rows = []
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(run_one, **task) for task in tasks]
            for future in as_completed(futures):
                rows.append(future.result())
        rows.sort(key=lambda row: (row["noise_index"], row["realisation"]))
        all_rows.extend(rows)
        summary = summarize_rows(parameter_index, parameters, rows)
        summaries.append(summary)
        print_summary(summary)
        write_outputs(all_rows, summaries)

    best = sorted(
        summaries,
        key=lambda row: (
            row["average_score"],
            row["lowest_score"],
            row["bootstrap_90_coverage_percent"],
        ),
        reverse=True,
    )[:5]

    print("\nBest 5")
    for rank, summary in enumerate(best, start=1):
        print(
            f"{rank}. params={summary['parameter_index']:02d} "
            f"avg_score={summary['average_score']:.3f} "
            f"lowest={summary['lowest_score']:.3f} "
            f"coverage90={summary['bootstrap_90_coverage_percent']:.1f}% "
            f"spent={summary['average_spent_hqc']:.1f}"
        )
        print(
            "   "
            f"initial_bracket_fraction={summary['initial_bracket_fraction']:.3f}, "
            f"trace_gate_growth={summary['trace_gate_growth']:.3f}, "
            f"trace_gate_shrink={summary['trace_gate_shrink']:.3f}, "
            f"bracket_hint={summary['bracket_hint']}, "
            f"n_gates_resolution={summary['n_gates_resolution']:.3f}, "
            f"max_hqc_per_ratio={summary['max_hqc_per_ratio']:.1f}, "
            f"initial_refine_hqc={summary['initial_refine_hqc']:.1f}, "
            f"surface_min_configs={summary['surface_min_configs']}, "
            f"monotonic_slack_fraction={summary['monotonic_slack_fraction']:.2f}"
        )

    print(f"\nWrote sweep results to {FIG_DIR}")


if __name__ == "__main__":
    main()
