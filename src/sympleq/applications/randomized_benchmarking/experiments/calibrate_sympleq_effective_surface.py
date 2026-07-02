"""
Calibrate the effective SympleQ RB decay surface as a function of register size.

This script directly measures dense depth sweeps on the local backend and
extracts effective rates

    log2 F(n, r, Q) = a0(r,Q) - D(r,Q) n,
    ln(2) D(r,Q) = (1-r) lambda_1(Q) + r lambda_2(Q),

where F is the asymptote-renormalised raw survival probability.  The resulting
smooth lambda_i(Q) curves are saved as a reference surface that
``cost_aware_surface_design.py`` can use for plots and scoring via
``reference_surface_path=...``.
"""
from __future__ import annotations

import csv
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from numpy.random import default_rng

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.base import MeasurementRequest
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    measured_items,
    spend_request_batch,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_design import (
    BASE_1Q_PAULI_ERROR,
    BASE_2Q_PAULI_ERROR,
    CostAwareSurfaceSettings,
    _LN2,
    _analytic_lindblad_gates,
    _asymptote,
    make_config_q,
)


EXPERIMENTS_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = EXPERIMENTS_DIR / "figs" / "sympleq_effective_surface_calibration"


def _seed_label(settings: "SympleQEffectiveSurfaceCalibrationSettings") -> str:
    return f"seed{settings.rng_seed}" if settings.rng_seed is not None else "seednone"


def _output_paths(
    output_dir: Path,
    settings: "SympleQEffectiveSurfaceCalibrationSettings",
) -> dict[str, Path]:
    seed = _seed_label(settings)
    return {
        "reference": output_dir / f"effective_surface_reference_{seed}.json",
        "measurements": output_dir / f"measurement_rows_{seed}.csv",
        "decay_slopes": output_dir / f"decay_slope_rows_{seed}.csv",
        "lambdas": output_dir / f"lambda_rows_{seed}.csv",
        "slice_fits_plot": output_dir / f"decay_slice_fits_{seed}.png",
        "decay_slopes_plot": output_dir / f"effective_decay_slopes_{seed}.png",
        "lambdas_plot": output_dir / f"effective_lambda_q_{seed}.png",
    }


def _migrate_legacy_outputs(
    output_dir: Path,
    paths: dict[str, Path],
    *,
    progress: bool,
) -> None:
    legacy = {
        "reference": output_dir / "effective_surface_reference.json",
        "measurements": output_dir / "measurement_rows.csv",
        "decay_slopes": output_dir / "decay_slope_rows.csv",
        "lambdas": output_dir / "lambda_rows.csv",
        "slice_fits_plot": output_dir / "decay_slice_fits.png",
        "decay_slopes_plot": output_dir / "effective_decay_slopes.png",
        "lambdas_plot": output_dir / "effective_lambda_q.png",
    }
    for key, old_path in legacy.items():
        new_path = paths[key]
        if old_path.exists() and not new_path.exists():
            old_path.rename(new_path)
            if progress:
                print(f"renamed legacy output: {old_path.name} -> {new_path.name}")


@dataclass(frozen=True)
class SympleQEffectiveSurfaceCalibrationSettings:
    q_values: tuple[int, ...] = tuple(range(5, 51, 5))
    ratios: tuple[float, ...] = (0.03, 0.05, 0.07, 0.10, 0.20, 0.30, 0.45, 0.60, 0.75, 0.85, 0.87, 0.90)
    depth_multiples: tuple[float, ...] = (
        0.35, 0.50, 0.65, 0.80, 0.95, 1.10, 1.30, 1.55, 1.85, 2.20,
    )
    shots_per_point: int = 500
    rng_seed: int | None = 12345
    n_gates_bounds: tuple[int, int] = (100, 15000)
    ratio_bounds: tuple[float, float] = (0.05, 0.90)
    random_elimination: float = 0.1
    use_scrambler: bool = True
    fit_f_range: tuple[float, float] = (0.15, 0.85)
    # Optional transient guard: shallow-depth points can sit off the asymptotic
    # RB line.  These settings exclude them from the straight-line decay fit
    # while keeping them visible as grey points in the diagnostic plot.
    min_fit_depth_multiple: float | None = None
    drop_shallowest_fit_points: int = 1
    min_points_per_slice: int = 4
    realized_ratio_tolerance: float | None = None
    lambda_poly_degree: int = 2
    output_dir: str | Path = DEFAULT_OUTPUT_DIR
    show_plots: bool = False
    progress: bool = True
    reuse_existing_measurements: bool = True
    # Parallelises the expensive measurement phase over independent (Q, r)
    # depth sweeps.  Keep 1 for deterministic serial debugging.
    n_workers: int = 36


def _round_even(value: float) -> int:
    return max(2, int(2 * round(value / 2)))


def _build_design_settings(
    settings: SympleQEffectiveSurfaceCalibrationSettings,
) -> CostAwareSurfaceSettings:
    q_values = tuple(sorted({int(q) for q in settings.q_values}))
    return CostAwareSurfaceSettings(
        q_values=q_values,
        n_qubits=int(round(float(np.median(q_values)))),
        n_gates_bounds=settings.n_gates_bounds,
        ratio_bounds=settings.ratio_bounds,
        rng_seed=settings.rng_seed,
        random_elimination=settings.random_elimination,
        use_scrambler=settings.use_scrambler,
        plot=False,
        save_path=None,
        verbose=False,
    )


def _local_rmb(design: CostAwareSurfaceSettings) -> tuple[np.random.Generator, RMB]:
    rng = default_rng(design.rng_seed)
    rmb = RMB.default(rng).with_backend(design.backend_factory(design, rng))
    return rng, rmb


def _measure_requests(
    requests: list[MeasurementRequest],
    design: CostAwareSurfaceSettings,
) -> RMB:
    rng, rmb = _local_rmb(design)
    data = rmb._data
    backend = rmb.backend
    spend_request_batch(backend, rng, data, requests, seed=design.rng_seed)
    return rmb


def _probe_requests_for_slice(
    settings: SympleQEffectiveSurfaceCalibrationSettings,
    design: CostAwareSurfaceSettings,
    q: int,
    ratio: float,
) -> tuple[list[MeasurementRequest], dict[RMBConfig, dict]]:
    requests: list[MeasurementRequest] = []
    metadata: dict[RMBConfig, dict] = {}
    n_star = float(
        _analytic_lindblad_gates(
            np.array([ratio]),
            np.array([float(min(settings.q_values))]),
            design,
        ).ravel()[0]
    )
    for multiple in settings.depth_multiples:
        n_gates = _round_even(np.clip(multiple * n_star, *settings.n_gates_bounds))
        config = make_config_q(design, n_gates, float(ratio), int(q))
        if config in metadata:
            continue
        metadata[config] = {
            "target_ratio": float(ratio),
            "target_q": int(q),
            "target_n_gates": int(n_gates),
            "depth_multiple": float(multiple),
            "realized_ratio": float(config.ratio_2_qb_gates),
            "realized_n_gates": int(config.n_gates),
        }
        requests.append(MeasurementRequest(config, int(settings.shots_per_point)))
    return requests, metadata


def _slice_seed(
    base_seed: int | None,
    q: int,
    ratio: float,
) -> int | None:
    if base_seed is None:
        return None
    return int(base_seed + 1_000_003 * int(q) + round(1_000_000 * float(ratio)))


def _measure_slice_worker(payload: tuple) -> tuple[list[dict], dict]:
    settings, q, ratio = payload
    slice_settings = replace(settings, rng_seed=_slice_seed(settings.rng_seed, q, ratio))
    design = _build_design_settings(slice_settings)
    requests, metadata = _probe_requests_for_slice(slice_settings, design, int(q), float(ratio))
    rmb = _measure_requests(requests, design)
    rows = _measurement_rows(rmb, metadata, design, slice_settings)
    return rows, {
        "q": int(q),
        "ratio": float(ratio),
        "configs": len(rows),
        "shots": sum(int(row["shots"]) for row in rows),
    }


def _measure_slices(
    settings: SympleQEffectiveSurfaceCalibrationSettings,
    slices: list[tuple[int, float]] | None = None,
) -> list[dict]:
    if slices is None:
        slices = [(int(q), float(r)) for q in settings.q_values for r in settings.ratios]
    tasks = [(settings, int(q), float(r)) for q, r in slices]
    rows: list[dict] = []
    total = len(tasks)

    def record_result(slice_rows: list[dict], stats: dict, done: int) -> None:
        rows.extend(slice_rows)
        if settings.progress:
            print(
                f"  slice {done:>3d}/{total:<3d} "
                f"Q={stats['q']:<3d} r={stats['ratio']:.2f}: "
                f"{stats['configs']} configs, {stats['shots']} shots",
                flush=True,
            )

    if int(settings.n_workers) <= 1:
        for index, task in enumerate(tasks, start=1):
            slice_rows, stats = _measure_slice_worker(task)
            record_result(slice_rows, stats, index)
    else:
        with ProcessPoolExecutor(max_workers=max(1, int(settings.n_workers))) as pool:
            futures = [pool.submit(_measure_slice_worker, task) for task in tasks]
            for done, future in enumerate(as_completed(futures), start=1):
                slice_rows, stats = future.result()
                record_result(slice_rows, stats, done)
    return rows


def _fit_weighted_line(x: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0.0)
    x, y, sigma = x[keep], y[keep], sigma[keep]
    if len(x) < 2:
        raise ValueError("Need at least two points for a weighted line fit.")
    X = np.column_stack([np.ones_like(x), x])
    w = 1.0 / np.maximum(sigma, 1e-12) ** 2
    xtw = X.T * w
    cov0 = np.linalg.pinv(xtw @ X)
    beta = cov0 @ (xtw @ y)
    residuals = y - X @ beta
    dof = max(len(x) - 2, 1)
    chi2 = float(np.sum(w * residuals * residuals))
    scale = max(chi2 / dof, 1.0)
    cov = cov0 * scale
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    ss_res = float(np.sum(residuals * residuals))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    return {
        "intercept": float(beta[0]),
        "slope": float(beta[1]),
        "intercept_stderr": float(np.sqrt(max(cov[0, 0], 0.0))),
        "slope_stderr": float(np.sqrt(max(cov[1, 1], 0.0))),
        "r2": float(r2),
        "chi2": chi2,
        "dof": int(dof),
    }


def _measurement_rows(
    rmb,
    metadata: dict[RMBConfig, dict],
    design: CostAwareSurfaceSettings,
    settings: SympleQEffectiveSurfaceCalibrationSettings,
) -> list[dict]:
    rows: list[dict] = []
    for config, estimator in measured_items(rmb._data):
        meta = metadata.get(config)
        if meta is None:
            continue
        counts = estimator.counts()
        successes = int(counts.get(True, 0))
        failures = int(counts.get(False, 0))
        shots = successes + failures
        if shots <= 0:
            continue
        b = _asymptote(design, float(config.n_qubits))
        # Jeffreys-style finite-sample smoothing keeps log errors finite without
        # moving high-shot estimates noticeably.
        p_hat = (successes + 0.5) / (shots + 1.0)
        f_hat = (p_hat - b) / max(1.0 - b, 1e-12)
        f_hat = float(np.clip(f_hat, 1e-12, 1.0 - 1e-12))
        p_var = max(p_hat * (1.0 - p_hat) / max(shots, 1), 1e-12)
        y_sigma = np.sqrt(p_var) / max(np.log(2.0) * abs(p_hat - b), 1e-12)
        row = {
            **meta,
            "n_qubits": int(config.n_qubits),
            "n_1qb_gates": int(config.n_1qb_gates),
            "n_2qb_gates": int(config.n_2qb_gates),
            "n_gates": int(config.n_gates),
            "ratio": float(config.ratio_2_qb_gates),
            "successes": successes,
            "failures": failures,
            "shots": shots,
            "p_hat": float(p_hat),
            "asymptote": float(b),
            "f_hat": f_hat,
            "log2_f_hat": float(np.log2(f_hat)),
            "log2_f_sigma": float(y_sigma),
            "in_fit_range": bool(settings.fit_f_range[0] <= f_hat <= settings.fit_f_range[1]),
        }
        rows.append(row)
    return rows


def _select_fit_points(
    points: list[dict],
    ratio: float,
    settings: SympleQEffectiveSurfaceCalibrationSettings,
) -> list[dict]:
    fit_points = [row for row in points if row["in_fit_range"]]
    if settings.realized_ratio_tolerance is not None:
        fit_points = [
            row for row in fit_points
            if abs(float(row["realized_ratio"]) - ratio)
            <= settings.realized_ratio_tolerance
        ]
    if settings.min_fit_depth_multiple is not None:
        fit_points = [
            row for row in fit_points
            if float(row["depth_multiple"]) >= settings.min_fit_depth_multiple
        ]
    drop = max(0, int(settings.drop_shallowest_fit_points))
    if drop > 0 and len(fit_points) > drop:
        fit_points = sorted(
            fit_points,
            key=lambda row: (float(row["depth_multiple"]), float(row["n_gates"])),
        )[drop:]
    return fit_points


def _fit_decay_slices(
    measurement_rows: list[dict],
    settings: SympleQEffectiveSurfaceCalibrationSettings,
    design: CostAwareSurfaceSettings,
) -> list[dict]:
    slope_rows: list[dict] = []
    total_slices = len(settings.q_values) * len(settings.ratios)
    completed_slices = 0
    for q in settings.q_values:
        for ratio in settings.ratios:
            completed_slices += 1
            points = [
                row for row in measurement_rows
                if int(row["target_q"]) == int(q)
                and abs(float(row["target_ratio"]) - float(ratio)) < 1e-12
            ]
            fit_points = _select_fit_points(points, float(ratio), settings)
            if len(fit_points) < settings.min_points_per_slice:
                slope_rows.append({
                    "target_ratio": float(ratio),
                    "q": int(q),
                    "ok": False,
                    "reason": f"too few fit points ({len(fit_points)})",
                    "n_points": len(fit_points),
                    "n_measured_points": len(points),
                })
                if settings.progress:
                    print(
                        f"  fit slice {completed_slices:>3d}/{total_slices:<3d} "
                        f"Q={int(q):<3d} r={ratio:.2f}: skipped "
                        f"({len(fit_points)} fit points, {len(points)} measured)",
                        flush=True,
                    )
                continue
            n = np.asarray([row["n_gates"] for row in fit_points], dtype=float)
            y = np.asarray([row["log2_f_hat"] for row in fit_points], dtype=float)
            sigma = np.asarray([row["log2_f_sigma"] for row in fit_points], dtype=float)
            fit = _fit_weighted_line(n, y, sigma)
            d_eff = -fit["slope"]
            d_stderr = max(fit["slope_stderr"], 1e-12)
            r_real = float(np.mean([row["ratio"] for row in fit_points]))
            n_star_analytic = float(
                _analytic_lindblad_gates(
                    np.array([r_real]), np.array([float(q)]), design
                ).ravel()[0]
            )
            d_analytic = 1.0 / n_star_analytic if n_star_analytic > 0.0 else float("nan")
            slope_rows.append({
                "target_ratio": float(ratio),
                "realized_ratio_mean": r_real,
                "realized_ratio_min": float(min(row["ratio"] for row in points)),
                "realized_ratio_max": float(max(row["ratio"] for row in points)),
                "q": int(q),
                "ok": bool(d_eff > 0.0 and np.isfinite(d_eff)),
                "reason": "",
                "n_points": len(fit_points),
                "n_measured_points": len(points),
                "min_fit_depth_multiple": settings.min_fit_depth_multiple,
                "drop_shallowest_fit_points": int(settings.drop_shallowest_fit_points),
                "D": float(d_eff),
                "D_stderr": float(d_stderr),
                "n_star": float(1.0 / d_eff) if d_eff > 0.0 else float("nan"),
                "intercept": fit["intercept"],
                "intercept_stderr": fit["intercept_stderr"],
                "r2": fit["r2"],
                "chi2": fit["chi2"],
                "dof": fit["dof"],
                "D_analytic": float(d_analytic),
                "D_over_analytic": float(d_eff / d_analytic)
                if d_analytic > 0.0 else float("nan"),
            })
            if settings.progress:
                print(
                    f"  fit slice {completed_slices:>3d}/{total_slices:<3d} "
                    f"Q={int(q):<3d} r={ratio:.2f}: "
                    f"D={d_eff:.4e}, n*={1.0 / d_eff:.1f}, "
                    f"pts={len(fit_points)}/{len(points)}, r_real={r_real:.3f}",
                    flush=True,
                )
    return slope_rows


def _fit_lambda_rows(slope_rows: list[dict]) -> list[dict]:
    lambda_rows: list[dict] = []
    q_values = sorted({int(row["q"]) for row in slope_rows})
    for q in q_values:
        rows = [row for row in slope_rows if row.get("ok") and int(row["q"]) == q]
        if len(rows) < 2:
            lambda_rows.append({
                "q": q,
                "ok": False,
                "reason": f"too few decay slopes ({len(rows)})",
            })
            continue
        r = np.asarray([row["realized_ratio_mean"] for row in rows], dtype=float)
        y = _LN2 * np.asarray([row["D"] for row in rows], dtype=float)
        sigma = _LN2 * np.asarray([row["D_stderr"] for row in rows], dtype=float)
        X = np.column_stack([1.0 - r, r])
        w = 1.0 / np.maximum(sigma, 1e-12) ** 2
        xtw = X.T * w
        cov0 = np.linalg.pinv(xtw @ X)
        beta = cov0 @ (xtw @ y)
        residuals = y - X @ beta
        dof = max(len(rows) - 2, 1)
        chi2 = float(np.sum(w * residuals * residuals))
        scale = max(chi2 / dof, 1.0)
        cov = cov0 * scale
        lambda_rows.append({
            "q": q,
            "ok": bool(beta[0] > 0.0 and beta[1] > 0.0),
            "reason": "",
            "lambda1": float(beta[0]),
            "lambda2": float(beta[1]),
            "lambda1_stderr": float(np.sqrt(max(cov[0, 0], 0.0))),
            "lambda2_stderr": float(np.sqrt(max(cov[1, 1], 0.0))),
            "n_slopes": len(rows),
            "chi2": chi2,
            "dof": int(dof),
        })
    return lambda_rows


def _fit_lambda_polynomials(
    lambda_rows: list[dict],
    settings: SympleQEffectiveSurfaceCalibrationSettings,
) -> dict:
    ok_rows = [row for row in lambda_rows if row.get("ok")]
    if len(ok_rows) < 2:
        reasons = "; ".join(
            f"Q={row.get('q')}: {row.get('reason', 'invalid fit')}"
            for row in lambda_rows
            if not row.get("ok")
        )
        raise RuntimeError(f"Not enough Q slices to fit lambda(Q). {reasons}")
    q = np.asarray([row["q"] for row in ok_rows], dtype=float)
    q_ref = float(0.5 * (min(settings.q_values) + max(settings.q_values)))
    x = q - q_ref
    degree = min(int(settings.lambda_poly_degree), len(ok_rows) - 1)
    l1 = np.asarray([row["lambda1"] for row in ok_rows], dtype=float)
    l2 = np.asarray([row["lambda2"] for row in ok_rows], dtype=float)
    s1 = np.asarray([max(row["lambda1_stderr"], 1e-12) for row in ok_rows], dtype=float)
    s2 = np.asarray([max(row["lambda2_stderr"], 1e-12) for row in ok_rows], dtype=float)
    c1 = np.polyfit(x, l1, degree, w=1.0 / s1)
    c2 = np.polyfit(x, l2, degree, w=1.0 / s2)
    return {
        "q_reference": q_ref,
        "lambda_poly_degree": int(degree),
        "lambda1_poly_coefficients": [float(v) for v in c1],
        "lambda2_poly_coefficients": [float(v) for v in c2],
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


_INT_ROW_FIELDS = {
    "target_q",
    "target_n_gates",
    "realized_n_gates",
    "n_qubits",
    "n_1qb_gates",
    "n_2qb_gates",
    "n_gates",
    "successes",
    "failures",
    "shots",
}
_FLOAT_ROW_FIELDS = {
    "target_ratio",
    "depth_multiple",
    "realized_ratio",
    "ratio",
    "p_hat",
    "asymptote",
    "f_hat",
    "log2_f_hat",
    "log2_f_sigma",
}
_BOOL_ROW_FIELDS = {"in_fit_range"}


def _read_measurement_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for raw in reader:
            row: dict = {}
            for key, value in raw.items():
                if value is None or value == "":
                    row[key] = value
                elif key in _INT_ROW_FIELDS:
                    row[key] = int(float(value))
                elif key in _FLOAT_ROW_FIELDS:
                    row[key] = float(value)
                elif key in _BOOL_ROW_FIELDS:
                    row[key] = value.strip().lower() in {"1", "true", "yes"}
                else:
                    row[key] = value
            rows.append(row)
    return rows


def _same_slice(row: dict, q: int, ratio: float) -> bool:
    return (
        int(row.get("target_q", -1)) == int(q)
        and abs(float(row.get("target_ratio", float("nan"))) - float(ratio)) < 1e-12
    )


def _row_config_key(row: dict) -> tuple[int, int, int]:
    return (
        int(row["n_qubits"]),
        int(row["n_1qb_gates"]),
        int(row["n_2qb_gates"]),
    )


def _expected_slice_keys(
    settings: SympleQEffectiveSurfaceCalibrationSettings,
    design: CostAwareSurfaceSettings,
    q: int,
    ratio: float,
) -> set[tuple[int, int, int]]:
    requests, _ = _probe_requests_for_slice(settings, design, int(q), float(ratio))
    return {
        (
            int(request.config.n_qubits),
            int(request.config.n_1qb_gates),
            int(request.config.n_2qb_gates),
        )
        for request in requests
    }


def _slice_cache_complete(
    rows: list[dict],
    settings: SympleQEffectiveSurfaceCalibrationSettings,
    design: CostAwareSurfaceSettings,
    q: int,
    ratio: float,
) -> bool:
    expected = _expected_slice_keys(settings, design, q, ratio)
    if not expected:
        return True
    existing: dict[tuple[int, int, int], dict] = {}
    for row in rows:
        if _same_slice(row, q, ratio):
            existing[_row_config_key(row)] = row
    for key in expected:
        row = existing.get(key)
        if row is None:
            return False
        if int(row.get("shots", 0)) < int(settings.shots_per_point):
            return False
    return True


def _missing_cached_slices(
    rows: list[dict],
    settings: SympleQEffectiveSurfaceCalibrationSettings,
    design: CostAwareSurfaceSettings,
) -> list[tuple[int, float]]:
    missing: list[tuple[int, float]] = []
    for q in settings.q_values:
        for ratio in settings.ratios:
            if not _slice_cache_complete(rows, settings, design, int(q), float(ratio)):
                missing.append((int(q), float(ratio)))
    return missing


def _drop_slices(
    rows: list[dict],
    slices: list[tuple[int, float]],
) -> list[dict]:
    if not slices:
        return rows
    return [
        row for row in rows
        if not any(_same_slice(row, q, ratio) for q, ratio in slices)
    ]


def _make_plots(
    paths: dict[str, Path],
    measurement_rows: list[dict],
    slope_rows: list[dict],
    lambda_rows: list[dict],
    polynomial: dict,
    settings: SympleQEffectiveSurfaceCalibrationSettings,
) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    import matplotlib
    if not settings.show_plots:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    q_values = tuple(sorted({int(q) for q in settings.q_values}))
    ratios = tuple(float(r) for r in settings.ratios)

    fig, axes = plt.subplots(
        len(q_values),
        len(ratios),
        figsize=(max(3.0 * len(ratios), 7.0), max(2.4 * len(q_values), 4.0)),
        sharex=False,
        sharey=False,
        squeeze=False,
    )
    slope_by_slice = {
        (int(row["q"]), float(row["target_ratio"])): row
        for row in slope_rows
    }
    for q_index, q in enumerate(q_values):
        for r_index, ratio in enumerate(ratios):
            ax = axes[q_index, r_index]
            points = [
                row for row in measurement_rows
                if int(row["target_q"]) == q
                and abs(float(row["target_ratio"]) - ratio) < 1e-12
            ]
            fit_points = _select_fit_points(points, ratio, settings)
            excluded = [row for row in points if row not in fit_points]
            if excluded:
                ax.scatter(
                    [row["n_gates"] for row in excluded],
                    [row["log2_f_hat"] for row in excluded],
                    color="0.7",
                    s=13,
                    label="excluded" if q_index == 0 and r_index == 0 else None,
                )
            if fit_points:
                ax.errorbar(
                    [row["n_gates"] for row in fit_points],
                    [row["log2_f_hat"] for row in fit_points],
                    yerr=[row["log2_f_sigma"] for row in fit_points],
                    marker="o",
                    linestyle="",
                    markersize=3.0,
                    capsize=1.5,
                    color="tab:blue",
                    label="fit points" if q_index == 0 and r_index == 0 else None,
                )
            slope_row = slope_by_slice.get((q, ratio))
            if slope_row is not None and slope_row.get("ok") and points:
                xs = np.asarray([row["n_gates"] for row in points], dtype=float)
                x_line = np.linspace(float(xs.min()), float(xs.max()), 100)
                y_line = float(slope_row["intercept"]) - float(slope_row["D"]) * x_line
                ax.plot(
                    x_line,
                    y_line,
                    color="tab:red",
                    linewidth=1.1,
                    label="fit" if q_index == 0 and r_index == 0 else None,
                )
                ax.text(
                    0.03,
                    0.05,
                    f"a0={float(slope_row['intercept']):+.2f}\n"
                    f"D={float(slope_row['D']):.2e}\n"
                    f"R2={float(slope_row['r2']):.2f}",
                    transform=ax.transAxes,
                    fontsize=7,
                    va="bottom",
                )
            ax.set_title(f"Q={q}, r={ratio:.2f}", fontsize=8)
            ax.tick_params(labelsize=7)
            if q_index == len(q_values) - 1:
                ax.set_xlabel("n", fontsize=8)
            if r_index == 0:
                ax.set_ylabel(r"$\log_2 F$", fontsize=8)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=3)
    fig.suptitle(r"Per-slice fits: $\log_2 F = a_0 - D n$")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(paths["slice_fits_plot"], dpi=150)

    ok_slopes = [row for row in slope_rows if row.get("ok")]
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    for ratio in settings.ratios:
        rows = [row for row in ok_slopes if abs(row["target_ratio"] - ratio) < 1e-12]
        if not rows:
            continue
        rows.sort(key=lambda row: row["q"])
        ax.errorbar(
            [row["q"] for row in rows],
            [row["D"] for row in rows],
            yerr=[row["D_stderr"] for row in rows],
            marker="o",
            capsize=2,
            label=f"r={ratio:.2f}",
        )
    ax.set_xlabel("Q", fontsize=13)
    ax.set_ylabel(r"$D(r,Q)=-d\log_2 F/dn$", fontsize=13)
    ax.set_title("Measured SympleQ effective decay slopes", fontsize=15)
    ax.tick_params(labelsize=11)
    ax.legend(ncols=2, fontsize=10)
    fig.tight_layout()
    fig.savefig(paths["decay_slopes_plot"], dpi=150)

    ok_lambda = [row for row in lambda_rows if row.get("ok")]
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.0), sharex=True)
    q = np.asarray([row["q"] for row in ok_lambda], dtype=float)
    x = q - float(polynomial["q_reference"])
    q_dense = np.linspace(min(settings.q_values), max(settings.q_values), 300)
    x_dense = q_dense - float(polynomial["q_reference"])
    for ax, name, color, ylabel in (
        (
            axes[0],
            "lambda1",
            "tab:blue",
            r"$\lambda_1(Q)=\ln 2\,D(0,Q)$",
        ),
        (
            axes[1],
            "lambda2",
            "tab:orange",
            r"$\lambda_2(Q)=\ln 2\,D(1,Q)$",
        ),
    ):
        y = np.asarray([row[name] for row in ok_lambda], dtype=float)
        yerr = np.asarray([row[f"{name}_stderr"] for row in ok_lambda], dtype=float)
        coeff = np.asarray(polynomial[f"{name}_poly_coefficients"], dtype=float)
        ax.errorbar(q, y, yerr=yerr, marker="o", linestyle="", capsize=2,
                    color=color, label="measured")
        ax.plot(q_dense, np.polyval(coeff, x_dense), color=color,
                label="polynomial")
        ax.set_ylabel(ylabel, fontsize=12)
        ax.tick_params(labelsize=11)
        ax.legend(fontsize=11)
    axes[1].set_xlabel("Q", fontsize=13)
    fig.suptitle(
        r"Extracted rates from: $\ln 2\,D(r,Q)=(1-r)\lambda_1(Q)+r\lambda_2(Q)$",
        fontsize=16,
    )
    fig.tight_layout()
    fig.savefig(paths["lambdas_plot"], dpi=150)
    if settings.show_plots:
        plt.show()
    plt.close("all")


def run_calibration(
    settings: SympleQEffectiveSurfaceCalibrationSettings
    = SympleQEffectiveSurfaceCalibrationSettings(),
) -> dict:
    output_dir = Path(settings.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = _output_paths(output_dir, settings)
    _migrate_legacy_outputs(output_dir, paths, progress=settings.progress)
    measurement_path = paths["measurements"]
    design = _build_design_settings(settings)
    n_configs = sum(
        len(_probe_requests_for_slice(settings, design, int(q), float(ratio))[0])
        for q in settings.q_values
        for ratio in settings.ratios
    )
    print(
        f"calibrating local SympleQ: {n_configs} configs, "
        f"{settings.shots_per_point} shots each, Q={settings.q_values}, "
        f"ratios={settings.ratios}, workers={max(1, int(settings.n_workers))}"
    )
    if settings.reuse_existing_measurements and measurement_path.exists():
        measurement_rows = _read_measurement_rows(measurement_path)
        print(f"loaded existing measurements: {measurement_path} ({len(measurement_rows)} rows)")
        missing_slices = _missing_cached_slices(measurement_rows, settings, design)
        if missing_slices:
            print(f"measuring {len(missing_slices)} missing/incomplete slices")
            new_rows = _measure_slices(settings, missing_slices)
            measurement_rows = _drop_slices(measurement_rows, missing_slices) + new_rows
        else:
            print("all requested slices found in cached measurements")
    else:
        measurement_rows = _measure_slices(settings)
    if settings.progress:
        print(f"fitting {len(settings.q_values) * len(settings.ratios)} decay slices")
    slope_rows = _fit_decay_slices(measurement_rows, settings, design)
    if settings.progress:
        print("fitting lambda_1(Q), lambda_2(Q)")
    lambda_rows = _fit_lambda_rows(slope_rows)
    if settings.progress:
        print("fitting smooth lambda_i(Q) reference curves")
    polynomial = _fit_lambda_polynomials(lambda_rows, settings)

    payload = {
        "version": 1,
        "label": "calibrated SympleQ",
        "model": "effective_lambda_polynomial_q",
        "polynomial_variable": "q_minus_reference",
        "backend": "sympleq",
        "base_one_q_pauli_error": BASE_1Q_PAULI_ERROR,
        "base_two_q_pauli_error": BASE_2Q_PAULI_ERROR,
        "settings": {
            key: (str(value) if isinstance(value, Path) else value)
            for key, value in asdict(settings).items()
        },
        **polynomial,
        "lambda_rows": lambda_rows,
        "slope_rows": slope_rows,
    }

    reference_path = paths["reference"]
    with open(reference_path, "w") as f:
        json.dump(payload, f, indent=2)
    _write_csv(measurement_path, measurement_rows)
    _write_csv(paths["decay_slopes"], slope_rows)
    _write_csv(paths["lambdas"], lambda_rows)
    _make_plots(paths, measurement_rows, slope_rows, lambda_rows, polynomial, settings)

    ok_slope_count = sum(1 for row in slope_rows if row.get("ok"))
    ok_lambda_count = sum(1 for row in lambda_rows if row.get("ok"))
    print("\nCalibration summary")
    print(f"  measured configs: {len(measurement_rows)}")
    print(f"  usable decay slices: {ok_slope_count} / {len(slope_rows)}")
    print(f"  usable Q lambda fits: {ok_lambda_count} / {len(lambda_rows)}")
    print(f"  reference surface: {reference_path}")
    print(f"  use with run_surface(reference_surface_path={str(reference_path)!r})")
    return payload


if __name__ == "__main__":
    show = os.environ.get("SYMPLEQ_EFFECTIVE_CALIBRATION_SHOW", "0") == "1"
    run_calibration(
        SympleQEffectiveSurfaceCalibrationSettings(
            output_dir=DEFAULT_OUTPUT_DIR,
            show_plots=True,
        )
    )
