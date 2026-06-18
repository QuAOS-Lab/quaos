"""
Evaluate the current ``monotone_tracing.py`` defaults over noisy simulations.

The script runs 20 measurement realisations for each of 20 fixed noise values,
where each one- and two-qubit noise scale is sampled within +/-20% of the
default SympleQ noise.  For each run it saves only the two-column uncertainty
diagnostic figure used by ``monotone_tracing.py``.

Each run is scored by comparing the fitted inverse parametric boundary against
the analytic Lindblad boundary for that run's noise scales.  Bootstrap coverage
reports how much of the analytic boundary lies inside the parametric fit's
50% and 90% bootstrap bands.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.experiments.common import (
    CrossingSettings,
)
from sympleq.applications.randomized_benchmarking.experiments.monotone_tracing import (
    BASE_1Q_PAULI_ERROR,
    BASE_2Q_PAULI_ERROR,
    MonotoneTracingSettings,
    run_with_budget,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    analytic_gate_counts,
    parametric_boundary_bootstrap,
    parametric_boundary_fit,
    parametric_boundary_gate_samples,
    plot_uncertainty_diagnostics,
)
from sympleq.core.noise.noise_model import GenericNoise


N_NOISE_VALUES = 5
N_REALISATIONS = 10
NOISE_SPREAD = 0.20
N_BOOTSTRAP = 100
FIG_DIR = Path(__file__).resolve().parent / "figs" / "monotone_tracing"
RESULTS_JSON = FIG_DIR / "results.json"
RUNS_CSV = FIG_DIR / "runs.csv"


def scaled_noise_backend_factory(
    one_q_noise_scale: float,
    two_q_noise_scale: float,
):
    """Build a local SympleQ backend with fixed scaled Lindblad noise."""
    def factory(settings: CrossingSettings, rng: RNGGenerator) -> RMBBackend:
        noise_model = GenericNoise.from_paulis(
            [BASE_1Q_PAULI_ERROR * one_q_noise_scale] * 3,
            rng,
        )
        two_qubit_noise_model = GenericNoise.from_paulis(
            [BASE_2Q_PAULI_ERROR * two_q_noise_scale] * 3,
            rng,
        )
        return SympleqBackend(
            noise_model=noise_model,
            two_qubit_noise_model=two_qubit_noise_model,
        )

    return factory


def noise_values() -> list[tuple[float, float]]:
    """Twenty deterministic noise-scale pairs within +/-NOISE_SPREAD."""
    rng = default_rng(20260618)
    return [
        tuple(float(x) for x in rng.uniform(
            1.0 - NOISE_SPREAD,
            1.0 + NOISE_SPREAD,
            size=2,
        ))
        for _ in range(N_NOISE_VALUES)
    ]


def fitted_gate_counts(data, ratios: np.ndarray) -> np.ndarray | None:
    """Parametric inverse-boundary gate counts on ``ratios``."""
    fit = parametric_boundary_fit(data)
    if fit is None:
        return None
    q, slope, _ = fit
    gates = 1.0 / (q + slope * ratios)
    gates[~np.isfinite(gates)] = np.nan
    gates[gates <= 0.0] = np.nan
    return gates


def fit_score(
    data,
    settings: MonotoneTracingSettings,
    *,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> float:
    """
    Score the fitted boundary against the analytic Lindblad boundary.

    The score is ``exp(-RMSE(log(fit / analytic)))`` over the configured ratio
    range, so 1 is an identical curve and lower values mean larger average
    multiplicative error.
    """
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 300)
    fitted = fitted_gate_counts(data, ratios)
    if fitted is None:
        return 0.0
    analytic = analytic_gate_counts(
        ratios,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    mask = (
        np.isfinite(fitted)
        & np.isfinite(analytic)
        & (fitted > 0.0)
        & (analytic > 0.0)
    )
    if not np.any(mask):
        return 0.0
    log_error = np.log(fitted[mask] / analytic[mask])
    return float(np.exp(-np.sqrt(np.mean(log_error**2))))


def bootstrap_coverage(
    data,
    settings: MonotoneTracingSettings,
    *,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
    seed: int,
) -> tuple[float, float]:
    """
    Fraction of ratio-grid points where the analytic line is inside the
    parametric bootstrap 50% and 90% bands.
    """
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 200)
    fits = parametric_boundary_bootstrap(data, n_bootstrap=N_BOOTSTRAP, seed=seed)
    samples = parametric_boundary_gate_samples(fits, ratios)
    if len(samples) == 0:
        return 0.0, 0.0

    analytic = analytic_gate_counts(
        ratios,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    q05, q25, q75, q95 = nan_quantiles(samples, (0.05, 0.25, 0.75, 0.95))
    valid50 = np.isfinite(q25) & np.isfinite(q75) & np.isfinite(analytic)
    valid90 = np.isfinite(q05) & np.isfinite(q95) & np.isfinite(analytic)
    if not np.any(valid50) and not np.any(valid90):
        return 0.0, 0.0

    coverage50 = (
        np.mean((q25[valid50] <= analytic[valid50])
                & (analytic[valid50] <= q75[valid50]))
        if np.any(valid50)
        else 0.0
    )
    coverage90 = (
        np.mean((q05[valid90] <= analytic[valid90])
                & (analytic[valid90] <= q95[valid90]))
        if np.any(valid90)
        else 0.0
    )
    return float(coverage50), float(coverage90)


def nan_quantiles(
    samples: np.ndarray,
    quantiles: tuple[float, ...],
) -> tuple[np.ndarray, ...]:
    """Column-wise nan-safe quantiles without all-NaN warnings."""
    output = [np.full(samples.shape[1], np.nan, dtype=float) for _ in quantiles]
    for column_index in range(samples.shape[1]):
        column = samples[:, column_index]
        column = column[np.isfinite(column)]
        if len(column) == 0:
            continue
        for output_array, quantile in zip(output, quantiles):
            output_array[column_index] = float(np.quantile(column, quantile))
    return tuple(output)


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "average": float(np.mean(array)) if len(array) else 0.0,
        "variance": float(np.var(array)) if len(array) else 0.0,
        "lowest": float(np.min(array)) if len(array) else 0.0,
    }


def run_one(
    *,
    noise_index: int,
    realisation: int,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> dict:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
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
    )
    # Plot helpers read these optional attributes when drawing analytic lines.
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

    out_path = FIG_DIR / f"noise_{noise_index:02d}_realisation_{realisation:02d}.png"
    plot_uncertainty_diagnostics(
        rmb._data,
        settings,
        png_path=out_path,
        show=False,
    )
    plt.close("all")

    return {
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
        "figure": str(out_path),
    }


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for noise_index, (one_q_scale, two_q_scale) in enumerate(noise_values()):
        for realisation in range(N_REALISATIONS):
            row = run_one(
                noise_index=noise_index,
                realisation=realisation,
                one_q_noise_scale=one_q_scale,
                two_q_noise_scale=two_q_scale,
            )
            rows.append(row)
            print(
                f"noise={noise_index:02d} realisation={realisation:02d}: "
                f"score={row['score']:.3f}, "
                f"coverage50={100.0 * row['bootstrap_50_coverage']:.1f}%, "
                f"coverage90={100.0 * row['bootstrap_90_coverage']:.1f}%, "
                f"spent={row['spent_hqc']:.1f} HQC, "
                f"crossings={row['n_crossings']}"
            )

    scores = [row["score"] for row in rows]
    coverage50 = [row["bootstrap_50_coverage"] for row in rows]
    coverage90 = [row["bootstrap_90_coverage"] for row in rows]
    spent = [row["spent_hqc"] for row in rows]
    summary = {
        "n_noise_values": N_NOISE_VALUES,
        "n_realisations": N_REALISATIONS,
        "n_runs": len(rows),
        "noise_spread": NOISE_SPREAD,
        "n_bootstrap": N_BOOTSTRAP,
        "score": summarize(scores),
        "bootstrap_50_coverage_percent": 100.0 * mean(coverage50),
        "bootstrap_90_coverage_percent": 100.0 * mean(coverage90),
        "average_spent_hqc": mean(spent),
    }
    RESULTS_JSON.write_text(
        json.dumps({"summary": summary, "runs": rows}, indent=2),
        encoding="utf-8",
    )
    write_csv(
        RUNS_CSV,
        rows,
        [
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
        ],
    )

    print("\nSummary")
    print(f"  average score: {summary['score']['average']:.3f}")
    print(f"  score variance: {summary['score']['variance']:.5f}")
    print(f"  lowest score: {summary['score']['lowest']:.3f}")
    print(
        "  analytic line inside bootstrap 50% band: "
        f"{summary['bootstrap_50_coverage_percent']:.1f}%"
    )
    print(
        "  analytic line inside bootstrap 90% band: "
        f"{summary['bootstrap_90_coverage_percent']:.1f}%"
    )
    print(f"  average spent: {summary['average_spent_hqc']:.1f} HQC")
    print(f"  wrote figures and results to {FIG_DIR}")


if __name__ == "__main__":
    main()
