"""
Benchmark ``cost_aware_surface_design.py`` on full 3-D Q-dependent surfaces.

This is the 3-D analogue of ``benchmarkbenchmark_scores.py``.  Each noise value
draws base one-/two-qubit Lindblad rates as in the original benchmark, then
adds slow linear Q-dependence across the chosen register-size range.  A
lightweight synthetic backend samples Bernoulli RMB outcomes from that truth,
so the benchmark can test whether the cost-aware design recovers the full
fidelity-0.5 surface

    n_*(r,Q) = ln(2) / (lambda_1(Q) + delta_2(Q) r).

Lower S1 and S2 are better:

* S1: integrated absolute error in log gate count over (r,Q), normalized by
  the fitted log-boundary surface scale.
* S2: integrated posterior standard deviation of log gate count over (r,Q),
  using the same normalization.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Callable

_SRC_ROOT = Path(__file__).resolve().parents[4]
sys.path = [path for path in sys.path if path != str(_SRC_ROOT)]
sys.path.insert(0, str(_SRC_ROOT))

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementOutcomes,
    MeasurementRequest,
    RMBBackend,
    ShotRNG,
)
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_design import (
    BASE_1Q_PAULI_ERROR,
    BASE_2Q_PAULI_ERROR,
    CostAwareSurfaceSettings,
    _I_L1,
    _I_L2,
    _I_M1,
    _I_M2,
    _LN2,
    _q_reference,
    _score_grid,
    _stateless_posterior,
    _surface_log_uncertainty,
    _wquantile,
    run_with_budget,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    REFERENCE_OFFSET,
    REFERENCE_SLOPE,
)
from sympleq.applications.randomized_benchmarking.experiments.scores import (
    integrate_trapezoid,
)
from sympleq.core.bayesian_estimation import BayesianEstimator


N_NOISE_VALUES = int(os.environ.get("COST_AWARE_3D_N_NOISE_VALUES", "15"))
N_REALISATIONS = int(os.environ.get("COST_AWARE_3D_N_REALISATIONS", "15"))
NOISE_SEED = int(os.environ.get("COST_AWARE_3D_NOISE_SEED", "20260618"))
NOISE_SPREAD = float(os.environ.get("COST_AWARE_3D_NOISE_SPREAD", "0.5"))
# Fractional rate change from Qref to either end of the Q range.
Q_SLOPE_SPREAD = float(os.environ.get("COST_AWARE_3D_Q_SLOPE_SPREAD", "0.12"))
Q_VALUES = tuple(range(5, 51, 5))
SCORE_RATIO_POINTS = int(os.environ.get("COST_AWARE_3D_SCORE_RATIO_POINTS", "60"))
HQC_BUDGET = float(os.environ.get("COST_AWARE_3D_HQC_BUDGET", "500.0"))
MAX_COST_PER_RUN = float(os.environ.get("COST_AWARE_3D_MAX_COST_PER_RUN", "45.0"))
MAX_WORKERS = int(os.environ.get("COST_AWARE_3D_WORKERS", "12")) or None

FIG_DIR = Path(__file__).resolve().parent / "figs" / "cost_aware_surface_3d"
RESULTS_JSON = FIG_DIR / "results.json"
RUNS_CSV = FIG_DIR / "runs.csv"


@dataclass(frozen=True)
class LinearQTruth:
    """Linear-Q truth for the RB boundary rates."""

    one_q_noise_scale: float
    two_q_noise_scale: float
    one_q_rate_q_slope: float
    two_q_rate_q_slope: float
    q_values: tuple[int, ...] = Q_VALUES
    visibility: float = 1.0

    @property
    def q_ref(self) -> float:
        return 0.5 * (min(self.q_values) + max(self.q_values))

    @property
    def q_half_span(self) -> float:
        return max(max(self.q_values) - self.q_ref, self.q_ref - min(self.q_values), 1.0)

    @property
    def lambda_1_ref(self) -> float:
        return REFERENCE_OFFSET * self.one_q_noise_scale

    @property
    def delta_2_ref(self) -> float:
        return REFERENCE_SLOPE * self.two_q_noise_scale

    def q_coordinate(self, q: float) -> float:
        return float((q - self.q_ref) / self.q_half_span)

    def rates(self, q: float) -> tuple[float, float]:
        x = self.q_coordinate(q)
        lambda_1 = self.lambda_1_ref * max(0.05, 1.0 + self.one_q_rate_q_slope * x)
        delta_2 = self.delta_2_ref * max(0.05, 1.0 + self.two_q_rate_q_slope * x)
        return float(lambda_1), float(delta_2)

    def denominator(self, ratio: float, q: float) -> float:
        lambda_1, delta_2 = self.rates(q)
        return (lambda_1 + delta_2 * ratio) / _LN2

    def boundary(self, ratio: float, q: float) -> float:
        d = self.denominator(ratio, q)
        return float("inf") if d <= 0.0 else 1.0 / d

    def probability(self, n_gates: float, ratio: float, q: float) -> float:
        b = 2.0 ** (-float(q))
        p = self.visibility * (1.0 - b) * 2.0 ** (-n_gates * self.denominator(ratio, q)) + b
        return float(np.clip(p, 1e-12, 1.0 - 1e-12))

    def as_row(self) -> dict[str, float]:
        return {
            "one_q_noise_scale": self.one_q_noise_scale,
            "two_q_noise_scale": self.two_q_noise_scale,
            "one_q_rate_q_slope": self.one_q_rate_q_slope,
            "two_q_rate_q_slope": self.two_q_rate_q_slope,
        }


class LinearQTruthBackend(RMBBackend):
    """Synthetic RMB backend for controlled Q-dependent rate benchmarks."""

    type = "linear_q_truth"

    def __init__(self, truth: LinearQTruth) -> None:
        self.truth = truth

    def fidelity_estimation(
        self,
        requests: list[MeasurementRequest],
        rng: RNGGenerator,
        shot_rng: ShotRNG | None = None,
    ) -> MeasurementOutcomes:
        outcomes: dict[RMBConfig, list[bool]] = {}
        shot_counts: dict[RMBConfig, int] = {}
        for request in requests:
            for _ in range(max(0, request.shots)):
                index = shot_counts.get(request.config, 0)
                shot_counts[request.config] = index + 1
                circuit_rng = rng if shot_rng is None else shot_rng(request.config, index)
                p = self.truth.probability(
                    request.config.n_gates,
                    request.config.ratio_2_qb_gates,
                    request.config.n_qubits,
                )
                outcomes.setdefault(request.config, []).append(bool(circuit_rng.random() < p))
        return MeasurementOutcomes(outcomes=outcomes)

    @classmethod
    def default_config(cls) -> RMBConfig:
        return RMBConfig.default()

    def default_estimator(self) -> BayesianEstimator:
        return BayesianEstimator.default()

    def to_dict(self) -> dict:
        return {"type": self.type, **self.truth.as_row()}

    @classmethod
    def from_dict(cls, payload: dict) -> "LinearQTruthBackend":
        truth = LinearQTruth(
            one_q_noise_scale=float(payload["one_q_noise_scale"]),
            two_q_noise_scale=float(payload["two_q_noise_scale"]),
            one_q_rate_q_slope=float(payload["one_q_rate_q_slope"]),
            two_q_rate_q_slope=float(payload["two_q_rate_q_slope"]),
        )
        return cls(truth)


def truth_values() -> list[LinearQTruth]:
    """Deterministic base-rate disorder plus slow linear Q-rate disorder."""
    rng = default_rng(NOISE_SEED)
    truths = []
    for _ in range(N_NOISE_VALUES):
        one_q_scale, two_q_scale = rng.uniform(
            1.0 - NOISE_SPREAD,
            1.0 + NOISE_SPREAD,
            size=2,
        )
        one_q_slope, two_q_slope = rng.uniform(
            -Q_SLOPE_SPREAD,
            Q_SLOPE_SPREAD,
            size=2,
        )
        truths.append(
            LinearQTruth(
                one_q_noise_scale=float(one_q_scale),
                two_q_noise_scale=float(two_q_scale),
                one_q_rate_q_slope=float(one_q_slope),
                two_q_rate_q_slope=float(two_q_slope),
            )
        )
    return truths


def backend_factory_for_truth(
    truth: LinearQTruth,
) -> Callable[[CostAwareSurfaceSettings, RNGGenerator], RMBBackend]:
    def factory(settings: CostAwareSurfaceSettings, rng: RNGGenerator) -> RMBBackend:
        return LinearQTruthBackend(truth)

    return factory


def make_settings(
    *,
    seed: int,
    truth: LinearQTruth,
) -> CostAwareSurfaceSettings:
    return CostAwareSurfaceSettings(
        q_values=Q_VALUES,
        n_qubits=int(round(truth.q_ref)),
        rng_seed=seed,
        plot=False,
        verbose=False,
        save_path=None,
        surface_plot_path=None,
        hqc_budget=HQC_BUDGET,
        max_cost_per_run=MAX_COST_PER_RUN,
        backend_factory=backend_factory_for_truth(truth),
        truth_boundary=truth.boundary,
        initial_one_q_pauli_error=BASE_1Q_PAULI_ERROR * truth.one_q_noise_scale,
        initial_two_q_pauli_error=BASE_2Q_PAULI_ERROR * truth.two_q_noise_scale,
    )


def _surface_average(values: np.ndarray, ratios: np.ndarray, qubits: np.ndarray) -> float:
    values = np.nan_to_num(values)
    r_span = max(float(max(ratios) - min(ratios)), 1e-12)
    if len(qubits) == 1:
        return float(integrate_trapezoid(values[0], ratios) / r_span)
    q_span = max(float(max(qubits) - min(qubits)), 1e-12)
    return float(
        integrate_trapezoid(integrate_trapezoid(values, ratios, axis=1), qubits)
        / (r_span * q_span)
    )


def _surface_integral(values: np.ndarray, ratios: np.ndarray, qubits: np.ndarray) -> float:
    values = np.nan_to_num(values)
    if len(qubits) == 1:
        return float(integrate_trapezoid(values[0], ratios))
    return float(integrate_trapezoid(integrate_trapezoid(values, ratios, axis=1), qubits))


def surface_scores(
    data,
    settings: CostAwareSurfaceSettings,
    truth: LinearQTruth,
) -> dict[str, float | int | None]:
    """Full 3-D S1/S2 scores against the generated Q-dependent truth."""
    posterior = _stateless_posterior(data, settings)
    if posterior[0] is None:
        return unavailable_scores()
    params, weights = posterior
    mean_var, _ = _surface_log_uncertainty(params, weights, _score_grid(settings), settings)
    log_rms_uncertainty = float(np.sqrt(mean_var))
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], SCORE_RATIO_POINTS)
    qubits = np.asarray(settings.q_values, dtype=float)
    rr, qq = np.meshgrid(ratios, qubits)
    q_ref = _q_reference(settings)
    dq = qq - q_ref
    med = np.asarray([float(_wquantile(params[:, j], weights, 0.5)) for j in range(5)])
    lam1 = med[_I_L1] + med[_I_M1] * dq
    lam2 = med[_I_L2] + med[_I_M2] * dq
    d_fit = (lam1 * (1.0 - rr) + lam2 * rr) / _LN2
    fit_log = np.where(d_fit > 0.0, -np.log(d_fit), np.nan)
    truth_gates = np.vectorize(truth.boundary)(rr, qq)
    truth_log = np.where(truth_gates > 0.0, np.log(truth_gates), np.nan)
    mean_d_gates = np.full_like(fit_log, np.nan, dtype=float)

    sigma = np.full_like(fit_log, np.nan, dtype=float)
    for q_index, q in enumerate(qubits):
        for r_index, r in enumerate(ratios):
            dq_i = q - q_ref
            lam1_i = params[:, _I_L1] + params[:, _I_M1] * dq_i
            lam2_i = params[:, _I_L2] + params[:, _I_M2] * dq_i
            d = (lam1_i * (1.0 - r) + lam2_i * r) / _LN2
            ok = (d > 0.0) & (weights > 0.0)
            if not np.any(ok):
                continue
            w = weights[ok] / np.sum(weights[ok])
            mean_d = float(np.sum(w * d[ok]))
            if mean_d > 0.0:
                mean_d_gates[q_index, r_index] = 1.0 / mean_d
            log_n = -np.log(d[ok])
            mean_log = float(np.sum(w * log_n))
            var_log = float(np.sum(w * log_n * log_n) - mean_log * mean_log)
            sigma[q_index, r_index] = float(np.sqrt(max(var_log, 0.0)))

    valid = np.isfinite(fit_log) & np.isfinite(truth_log) & np.isfinite(sigma)
    if np.count_nonzero(valid) < 2:
        return unavailable_scores()
    fit_log = np.where(valid, fit_log, np.nan)
    delta = np.where(valid, np.abs(fit_log - truth_log), np.nan)
    sigma = np.where(valid, sigma, np.nan)
    normalizer = max(abs(_surface_average(fit_log, ratios, qubits)), 1e-12)
    s1 = float(_surface_average(delta, ratios, qubits) / normalizer)
    s2 = float(_surface_average(sigma, ratios, qubits) / normalizer)
    valid_volume = np.isfinite(mean_d_gates) & np.isfinite(truth_gates)
    volume_fit = _surface_integral(np.where(valid_volume, mean_d_gates, np.nan), ratios, qubits)
    volume_truth = _surface_integral(np.where(valid_volume, truth_gates, np.nan), ratios, qubits)
    volume_ratio = (
        float(volume_fit / volume_truth)
        if np.isfinite(volume_fit) and np.isfinite(volume_truth) and abs(volume_truth) > 1e-12
        else None
    )
    return {
        "S1": s1,
        "S2": s2,
        "S_total": s1 + s2,
        "log_rms_uncertainty": log_rms_uncertainty,
        "volume_fit": float(volume_fit),
        "volume_truth": float(volume_truth),
        "volume_ratio": volume_ratio,
        "mean_delta_log_gates": float(np.nanmean(delta)),
        "mean_sigma_log_gates": float(np.nanmean(sigma)),
        "n_surface_points": int(np.count_nonzero(valid)),
    }


def unavailable_scores() -> dict[str, float | int | None]:
    return {
        "S1": None,
        "S2": None,
        "S_total": None,
        "log_rms_uncertainty": None,
        "volume_fit": None,
        "volume_truth": None,
        "volume_ratio": None,
        "mean_delta_log_gates": None,
        "mean_sigma_log_gates": None,
        "n_surface_points": 0,
    }


def format_optional(value: float | int | None, precision: int = 6) -> str:
    if value is None:
        return "unavailable"
    if isinstance(value, int):
        return str(value)
    return f"{value:.{precision}g}"


def print_row(row: dict) -> None:
    print(
        f"noise={row['noise_index']:02d} realisation={row['realisation']:02d}: "
        f"S1={format_optional(row['S1'])}, "
        f"S2={format_optional(row['S2'])}, "
        f"log-rms={format_optional(row['log_rms_uncertainty'])}, "
        f"volume ratio={format_optional(row['volume_ratio'])}, "
        f"spent={row['spent_hqc']:.1f} HQC, "
        f"configs={row['n_training_configs']}"
    )


def run_one(
    *,
    noise_index: int,
    realisation: int,
    truth: LinearQTruth,
) -> dict:
    seed = 10_000 * noise_index + realisation
    settings = make_settings(seed=seed, truth=truth)
    rmb, submitted, budget = run_with_budget(settings)
    scores = surface_scores(rmb._data, settings, truth)
    return {
        "approach": "cost_aware_surface_design",
        "noise_index": noise_index,
        "realisation": realisation,
        "seed": seed,
        **truth.as_row(),
        **scores,
        "spent_hqc": budget.spent_hqc,
        "n_training_configs": len(rmb._data),
        "n_submitted_configs": len(submitted),
        "q_values": list(Q_VALUES),
    }


def benchmark_tasks() -> list[dict]:
    tasks = []
    for noise_index, truth in enumerate(truth_values()):
        for realisation in range(N_REALISATIONS):
            tasks.append({
                "noise_index": noise_index,
                "realisation": realisation,
                "truth": truth,
            })
    return tasks


def summarize_optional(values: list[float | None]) -> dict[str, float | None]:
    finite = np.asarray([value for value in values if value is not None], dtype=float)
    if len(finite) == 0:
        return {"average": None, "variance": None, "best": None, "worst": None}
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
        f"Running cost_aware_surface_design 3D benchmark: "
        f"{N_NOISE_VALUES} noise values, {N_REALISATIONS} realisations each, "
        f"{max_workers} workers, budget {HQC_BUDGET:.1f} HQC"
    )
    print(
        f"Base-rate spread +/-{100 * NOISE_SPREAD:.1f}%, "
        f"linear Q-rate slope spread +/-{100 * Q_SLOPE_SPREAD:.1f}%"
    )

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_one, **task) for task in tasks]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print_row(row)

    rows.sort(key=lambda row: (row["noise_index"], row["realisation"]))
    summary = {
        "approach": "cost_aware_surface_design",
        "n_noise_values": N_NOISE_VALUES,
        "n_realisations": N_REALISATIONS,
        "n_runs": len(rows),
        "noise_seed": NOISE_SEED,
        "noise_spread": NOISE_SPREAD,
        "q_slope_spread": Q_SLOPE_SPREAD,
        "q_values": list(Q_VALUES),
        "score_ratio_points": SCORE_RATIO_POINTS,
        "hqc_budget": HQC_BUDGET,
        "max_cost_per_run": MAX_COST_PER_RUN,
        "S1": summarize_optional([row["S1"] for row in rows]),
        "S2": summarize_optional([row["S2"] for row in rows]),
        "S_total": summarize_optional([row["S_total"] for row in rows]),
        "log_rms_uncertainty": summarize_optional(
            [row["log_rms_uncertainty"] for row in rows]
        ),
        "volume_ratio": summarize_optional([row["volume_ratio"] for row in rows]),
        "mean_delta_log_gates": mean_optional([row["mean_delta_log_gates"] for row in rows]),
        "mean_sigma_log_gates": mean_optional([row["mean_sigma_log_gates"] for row in rows]),
        "average_spent_hqc": mean([row["spent_hqc"] for row in rows]) if rows else None,
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
        "one_q_noise_scale",
        "two_q_noise_scale",
        "one_q_rate_q_slope",
        "two_q_rate_q_slope",
        "S1",
        "S2",
        "S_total",
        "log_rms_uncertainty",
        "volume_fit",
        "volume_truth",
        "volume_ratio",
        "mean_delta_log_gates",
        "mean_sigma_log_gates",
        "n_surface_points",
        "spent_hqc",
        "n_training_configs",
        "n_submitted_configs",
        "q_values",
    ]
    write_csv(RUNS_CSV, rows, fieldnames)

    print("\nSummary")
    for score_name in ("S1", "S2", "S_total", "log_rms_uncertainty", "volume_ratio"):
        score_summary = summary[score_name]
        print(
            f"  {score_name}: "
            f"avg={format_optional(score_summary['average'])}, "
            f"var={format_optional(score_summary['variance'])}, "
            f"best={format_optional(score_summary['best'])}, "
            f"worst={format_optional(score_summary['worst'])}"
        )
    print(f"  mean delta log gates: {format_optional(summary['mean_delta_log_gates'])}")
    print(f"  mean sigma log gates: {format_optional(summary['mean_sigma_log_gates'])}")
    print(f"  average spent: {format_optional(summary['average_spent_hqc'])} HQC")
    print(f"  wrote results to {FIG_DIR}")


if __name__ == "__main__":
    main()
