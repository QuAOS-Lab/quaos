"""Gaussian process level-set estimator for fidelity target contour.

This script estimates the contour where fidelity(depth, ratio) = TARGET over a
smooth ratio range using a Gaussian process model. Fidelity uncertainties from
shot-based Bayesian estimation are baked into the GP observation noise.

The estimator uses:
- a GP surrogate over (ratio, depth)
- heteroscedastic observation noise from fidelity_std
- adaptive sampling near the estimated target contour
- posterior uncertainty propagation to depth error bars
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from datetime import datetime
import json
import os
from pathlib import Path
import warnings

os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES
from sympleq.core.noise.noise_model import GenericNoise


# -----------------------------
# Main knobs
# -----------------------------

TARGET = 0.5

STUDY_SEEDS = (5,10,15,20,30,40,50,60,70,)

N_QUBITS = 50

DEPTH_MIN = 10
DEPTH_MAX = 150

RATIO_MIN = 0.2
RATIO_MAX = 0.8
N_RATIO_POINTS = 16
RATIO_EVAL_POINTS = tuple(np.linspace(RATIO_MIN, RATIO_MAX, N_RATIO_POINTS))

INITIAL_RATIO_POINTS = 10
INITIAL_DEPTH_POINTS = 10
N_ADAPTIVE_TRIALS = 50

BISECTION_STEPS = 0  # not used in GP-level-set search
RUNS_PER_CONFIG = 10

EVALUATE_ENDPOINTS = False
TARGET_OVERLAP_SIGMA = 1.0

OUTPUT_DIR = Path(__file__).parent
RESULTS_DIR = OUTPUT_DIR / "monotone_bisection_results"

SINGLE_QUBIT_ERROR = 0.000025
TWO_QUBIT_GATE_ERROR = 0.00079
SINGLE_QUBIT_PAULI_ERROR = SINGLE_QUBIT_ERROR / 3
TWO_QUBIT_PER_QUBIT_ERROR = 1 - np.sqrt(1 - TWO_QUBIT_GATE_ERROR)
TWO_QUBIT_PAULI_ERROR = TWO_QUBIT_PER_QUBIT_ERROR / 3

TWO_QUBIT_GATE_REFERENCE = 100
PLOT_CMAPS = {
    "levelset_points": "viridis",
    "levelsets_colored": "coolwarm",
    "evaluations": "plasma",
}


@dataclass
class EvalRecord:
    seed: int
    ratio: float
    depth: int
    fidelity: float
    fidelity_std: float
    target: float
    signed_error: float
    absolute_error: float
    squared_error: float
    runs: int
    success_counts: dict
    requested_two_qubit_gates: int
    actual_two_qubit_gates: int
    two_qubit_gate_shortfall: int
    two_qubit_gate_warning: bool


@dataclass
class CrossingRecord:
    seed: int
    ratio: float
    crossing_depth: float
    crossing_depth_std: float
    crossing_method: str
    gp_pred_mean: float
    gp_pred_std: float
    gp_depth_slope: float
    closest_depth: int
    closest_fidelity: float
    closest_fidelity_std: float


def write_json(path: Path, data: dict | list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(data, f, indent=2)


def requested_two_qubit_gate_depth(
    ratios: np.ndarray,
    gate_count: int = TWO_QUBIT_GATE_REFERENCE,
    n_qubits: int = N_QUBITS,
) -> np.ndarray:
    return 2.0 * gate_count / (n_qubits * ratios)


def add_two_qubit_gate_reference_line(
    ax,
    gate_count: int = TWO_QUBIT_GATE_REFERENCE,
    n_qubits: int = N_QUBITS,
) -> None:
    ratios = np.linspace(max(RATIO_MIN, 1e-9), RATIO_MAX, 400)
    depths = requested_two_qubit_gate_depth(ratios, gate_count=gate_count, n_qubits=n_qubits)
    y_min = DEPTH_MIN - 4
    y_max = DEPTH_MAX + 4
    visible = (depths >= y_min) & (depths <= y_max)
    if not np.any(visible):
        return

    ax.plot(
        ratios[visible],
        depths[visible],
        color="#2458a6",
        linestyle="--",
        linewidth=1.7,
        alpha=0.95,
        label=f"{gate_count} requested two-qubit gates",
        zorder=4,
    )


def json_safe_counts(counts: dict) -> dict:
    return {str(k): int(v) for k, v in counts.items()}


def make_backend(rng: np.random.Generator) -> SympleqBackend:
    noise_model = GenericNoise.from_paulis([SINGLE_QUBIT_PAULI_ERROR] * 3, rng=rng)
    two_qubit_noise_model = GenericNoise.from_paulis([TWO_QUBIT_PAULI_ERROR] * 3, rng=rng)
    return SympleqBackend(noise_model, two_qubit_noise_model)


def make_config(depth: int, two_qubit_ratio: float) -> RMBConfig:
    return (
        RMBConfig.default()
        .with_depth(depth)
        .with_n_qubits(N_QUBITS)
        .with_two_qubit_gate_ratio(two_qubit_ratio)
        .with_scrambling_probability(0.85)
        .with_random_elimination(0.3)
        .with_gates_set((GATES.H, GATES.S, GATES.CX))
    )


def two_qubit_gate_stats(config: RMBConfig, two_qubit_ratio: float, seed: int) -> tuple[int, int, bool]:
    requested = int(two_qubit_ratio * config.n_qubits * config.depth) // 2

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        circuit = Circuit.from_depth(
            config.depth,
            config.dimensions,
            gates_set=config.gates_set,
            two_qudit_gate_ratio=two_qubit_ratio,
            rng=np.random.default_rng(seed),
        )

    actual = sum(gate.n_qudits == 2 for gate in circuit.gates)
    had_warning = any(
        "Could not satisfy the required number of 2-qudit gates" in str(w.message)
        for w in caught
    )
    return requested, actual, had_warning


def estimate_fidelity(
    config: RMBConfig,
    backend: SympleqBackend,
    rng: np.random.Generator,
) -> tuple[float, float, int, dict]:
    estimator = BayesianEstimator(
        threshold=2e-3,
        min_runs=RUNS_PER_CONFIG,
        max_runs=RUNS_PER_CONFIG,
    )

    for _ in estimator.run_iter(lambda c=config: backend.fidelity_estimation(c, rng)):
        pass

    fidelity = float(estimator.probability(True))
    fidelity_std = float(np.sqrt(estimator.variance(True)))
    runs = int(estimator.num_runs())
    counts = json_safe_counts(estimator.counts())
    return fidelity, fidelity_std, runs, counts


def evaluate_depth(
    seed: int,
    ratio: float,
    depth: int,
    backend: SympleqBackend,
    rng: np.random.Generator,
) -> EvalRecord:
    config = make_config(depth=depth, two_qubit_ratio=ratio)

    requested, actual, had_warning = two_qubit_gate_stats(
        config=config,
        two_qubit_ratio=ratio,
        seed=seed * 1_000_000 + int(round(10_000 * ratio)) + depth,
    )

    fidelity, fidelity_std, runs, counts = estimate_fidelity(
        config=config,
        backend=backend,
        rng=rng,
    )

    signed_error = fidelity - TARGET

    return EvalRecord(
        seed=int(seed),
        ratio=float(ratio),
        depth=int(depth),
        fidelity=float(fidelity),
        fidelity_std=float(fidelity_std),
        target=float(TARGET),
        signed_error=float(signed_error),
        absolute_error=float(abs(signed_error)),
        squared_error=float(signed_error**2),
        runs=runs,
        success_counts=counts,
        requested_two_qubit_gates=int(requested),
        actual_two_qubit_gates=int(actual),
        two_qubit_gate_shortfall=int(requested - actual),
        two_qubit_gate_warning=bool(had_warning),
    )


def build_gp(X: np.ndarray, y: np.ndarray, y_std: np.ndarray) -> GaussianProcessRegressor:
    kernel = ConstantKernel(1.0, constant_value_bounds=(1e-3, 1e3)) * Matern(
        length_scale=[0.2, 20.0], length_scale_bounds=[(1e-2, 100.0), (1.0, 100.0)], nu=2.5
    ) + WhiteKernel(noise_level=1e-6, noise_level_bounds=(1e-8, 1e0))
    gp = GaussianProcessRegressor(kernel=kernel, alpha=np.maximum(y_std**2, 1e-8), normalize_y=True, n_restarts_optimizer=5)
    gp.fit(X, y)
    return gp


def predict_gp(gp: GaussianProcessRegressor, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mu, sigma = gp.predict(X, return_std=True)
    return mu, sigma


def make_initial_design() -> np.ndarray:
    ratios = np.linspace(RATIO_MIN, RATIO_MAX, INITIAL_RATIO_POINTS)
    depths = np.linspace(DEPTH_MIN, DEPTH_MAX, INITIAL_DEPTH_POINTS)
    grid = np.array([[r, d] for r in ratios for d in depths], dtype=float)
    return grid


def adaptive_acquisition(mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    proximity = np.exp(-0.5 * ((mu - TARGET) / 0.08)**2)
    score = proximity * (sigma + 0.05)
    return score


def find_crossing_depth_for_ratio(
    gp: GaussianProcessRegressor,
    ratio: float,
    depth_grid: np.ndarray,
    mean_curve: np.ndarray,
) -> tuple[float | None, float | None, float | None, float | None, float | None]:
    target = TARGET
    if mean_curve[0] < target and mean_curve[-1] < target:
        return None, None, None, None, None
    if mean_curve[0] > target and mean_curve[-1] > target:
        return None, None, None, None, None

    sign = np.sign(mean_curve - target)
    crossings = np.where(sign[:-1] * sign[1:] <= 0)[0]
    if len(crossings) == 0:
        closest_idx = np.argmin(np.abs(mean_curve - target))
        return float(depth_grid[closest_idx]), float(mean_curve[closest_idx]), None, None, None

    idx = crossings[0]
    d0, d1 = depth_grid[idx], depth_grid[idx + 1]
    f0, f1 = mean_curve[idx], mean_curve[idx + 1]
    if abs(f1 - f0) < 1e-12:
        depth_cross = 0.5 * (d0 + d1)
    else:
        depth_cross = d0 + (target - f0) * (d1 - d0) / (f1 - f0)

    delta = 1.0
    depth_lo = max(DEPTH_MIN, depth_cross - delta)
    depth_hi = min(DEPTH_MAX, depth_cross + delta)
    eval_pts = np.array([[ratio, depth_lo], [ratio, depth_hi]], dtype=float)
    mu_vals, _ = gp.predict(eval_pts, return_std=True)
    slope = (mu_vals[1] - mu_vals[0]) / (depth_hi - depth_lo)

    pred_mu, pred_sigma = gp.predict(np.array([[ratio, depth_cross]]), return_std=True)
    grad = abs(slope)
    depth_std = float(pred_sigma[0] / np.maximum(grad, 1e-3)) if grad > 0 else float(pred_sigma[0] * 10.0)
    depth_std = float(max(depth_std, 0.5))

    return float(depth_cross), float(pred_mu[0]), float(pred_sigma[0]), float(depth_std), float(slope)


def run_gp_levelset(seed: int, run_dir: Path) -> tuple[list[EvalRecord], list[CrossingRecord]]:
    seed_dir = run_dir / f"seed_{seed}"
    seed_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    backend = make_backend(rng)
    print(f"[seed {seed}] Backend initialized")

    train_X = []
    train_y = []
    train_y_std = []
    records: list[EvalRecord] = []

    design = make_initial_design()
    if EVALUATE_ENDPOINTS:
        ratios = np.linspace(RATIO_MIN, RATIO_MAX, INITIAL_RATIO_POINTS)
        depths = np.array([DEPTH_MIN, DEPTH_MAX], dtype=float)
        endpoints = np.array([[r, d] for r in ratios for d in depths], dtype=float)
        design = np.vstack((design, endpoints))

    design = np.unique(design, axis=0)
    print(f"[seed {seed}] Initial design prepared: {len(design)} unique points")

    for ratio, depth in design:
        record = evaluate_depth(seed=seed, ratio=float(ratio), depth=int(depth), backend=backend, rng=rng)
        records.append(record)
        train_X.append([ratio, depth])
        train_y.append(record.fidelity)
        train_y_std.append(record.fidelity_std)

    print(f"[seed {seed}] Initial design grid evaluated: {len(records)} points")

    for trial in range(N_ADAPTIVE_TRIALS):
        X = np.array(train_X, dtype=float)
        y = np.array(train_y, dtype=float)
        y_std = np.array(train_y_std, dtype=float)
        gp = build_gp(X, y, y_std)

        candidate_ratios = np.linspace(RATIO_MIN, RATIO_MAX, 41)
        candidate_depths = np.linspace(DEPTH_MIN, DEPTH_MAX, 41)
        candidates = np.array([[r, d] for r in candidate_ratios for d in candidate_depths], dtype=float)

        mu, sigma = predict_gp(gp, candidates)
        score = adaptive_acquisition(mu, sigma)

        sampled = np.array(train_X, dtype=float)
        if sampled.size > 0:
            diff = np.abs(candidates[:, np.newaxis, :] - sampled[np.newaxis, :, :])
            all_close = np.all(diff < 1e-6, axis=2)
            already_seen = np.any(all_close, axis=1)
            if np.all(already_seen):
                print(f"[seed {seed}] All candidates exhausted at trial {trial+1}")
                break
            score = np.where(already_seen, -np.inf, score)

        best_idx = int(np.argmax(score))
        best_ratio, best_depth = float(candidates[best_idx, 0]), float(candidates[best_idx, 1])

        record = evaluate_depth(seed=seed, ratio=best_ratio, depth=int(best_depth), backend=backend, rng=rng)
        records.append(record)
        train_X.append([best_ratio, best_depth])
        train_y.append(record.fidelity)
        train_y_std.append(record.fidelity_std)

        print(
            f"[seed {seed} | trial {trial+1} | ratio={best_ratio:.4f} | depth={best_depth:.1f}] "
            f"fidelity={record.fidelity:.4f} ± {record.fidelity_std:.4f}"
        )

    print(f"[seed {seed}] Adaptive trials complete: {len(records)} total points evaluated")

    X = np.array(train_X, dtype=float)
    y = np.array(train_y, dtype=float)
    y_std = np.array(train_y_std, dtype=float)
    gp = build_gp(X, y, y_std)
    print(f"[seed {seed}] Final GP fit complete using {len(train_X)} training points")

    crossings: list[CrossingRecord] = []
    ratio_grid = np.array(RATIO_EVAL_POINTS, dtype=float)
    depth_grid = np.linspace(DEPTH_MIN, DEPTH_MAX, 201, dtype=float)

    for ratio in ratio_grid:
        query = np.array([[ratio, d] for d in depth_grid], dtype=float)
        mean_curve, _ = predict_gp(gp, query)
        crossing_depth, pred_mu, pred_sigma, depth_std, slope = find_crossing_depth_for_ratio(gp, ratio, depth_grid, mean_curve)

        closest_idx = np.argmin(np.abs(mean_curve - TARGET))
        closest_depth = int(depth_grid[closest_idx])
        closest_fidelity = float(mean_curve[closest_idx])
        closest_fidelity_std = float(_[closest_idx])

        method = "gp_mean_root" if crossing_depth is not None else "gp_closest_observed"
        if crossing_depth is None:
            crossing_depth = float(closest_depth)
            depth_std = float(max(0.5 * (DEPTH_MAX - DEPTH_MIN), 1.0))
            pred_mu = float(closest_fidelity)
            pred_sigma = float(closest_fidelity_std)
            slope = 0.0

        crossings.append(CrossingRecord(
            seed=int(seed),
            ratio=float(ratio),
            crossing_depth=float(crossing_depth),
            crossing_depth_std=float(depth_std),
            crossing_method=method,
            gp_pred_mean=float(pred_mu),
            gp_pred_std=float(pred_sigma),
            gp_depth_slope=float(slope),
            closest_depth=closest_depth,
            closest_fidelity=closest_fidelity,
            closest_fidelity_std=closest_fidelity_std,
        ))

    print(f"[seed {seed}] Crossing extraction complete: {len(crossings)} ratios")

    write_json(seed_dir / "records.json", [asdict(r) for r in records])
    print(f"[seed {seed}] Saved evaluation records: {seed_dir / 'records.json'}")
    write_json(seed_dir / "crossings.json", [asdict(c) for c in crossings])
    print(f"[seed {seed}] Saved crossing records: {seed_dir / 'crossings.json'}")
    write_json(seed_dir / "gp_training_data.json", {
        "X": train_X,
        "y": train_y,
        "y_std": train_y_std,
    })
    print(f"[seed {seed}] Saved GP training data: {seed_dir / 'gp_training_data.json'}")

    return records, crossings


def plot_gp_levelset(crossings_by_seed: dict[int, list[CrossingRecord]], records_by_seed: dict[int, list[EvalRecord]], run_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    total_rb_runs = sum(record.runs for records in records_by_seed.values() for record in records)

    seed_colors = {seed: plt.get_cmap("tab10")(idx % 10) for idx, seed in enumerate(sorted(crossings_by_seed))}

    for seed, crossings in sorted(crossings_by_seed.items()):
        crossings_sorted = sorted(crossings, key=lambda c: c.ratio)
        ratios = np.array([c.ratio for c in crossings_sorted], dtype=float)
        depths = np.array([c.crossing_depth for c in crossings_sorted], dtype=float)
        depth_stds = np.array([c.crossing_depth_std for c in crossings_sorted], dtype=float)

        ax.plot(ratios, depths, color=seed_colors[seed], alpha=0.85, linewidth=2.0, label=f"seed {seed} contour")
        ax.fill_between(
            ratios,
            depths - depth_stds,
            depths + depth_stds,
            color=seed_colors[seed],
            alpha=0.18,
            linewidth=0,
        )

    # Overlay measured evaluation points
    all_ratio = []
    all_depth = []
    all_fidelity = []
    all_std = []
    for records in records_by_seed.values():
        for record in records:
            all_ratio.append(record.ratio)
            all_depth.append(record.depth)
            all_fidelity.append(record.fidelity)
            all_std.append(record.fidelity_std)

    if all_ratio:
        all_ratio = np.array(all_ratio)
        all_depth = np.array(all_depth)
        all_fidelity = np.array(all_fidelity)
        all_std = np.array(all_std)
        size = 24 + 240 * np.clip(all_std, 0.0, 0.25)
        scatter = ax.scatter(
            all_ratio,
            all_depth,
            c=all_fidelity,
            s=size,
            edgecolors="black",
            linewidths=0.3,
            alpha=0.6,
            cmap=PLOT_CMAPS["levelset_points"],
            vmin=0.0,
            vmax=1.0,
            label="measured evaluations",
            zorder=2,
        )
        colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
        colorbar.set_label("observed fidelity")

    ax.set_xlabel("two_qubit_ratio")
    ax.set_ylabel("depth at fidelity=0.5")
    ax.set_title("GP level-set estimate with measured evaluation points")
    ax.text(0.02, 0.98, f"total RB runs: {total_rb_runs}", transform=ax.transAxes, va="top", ha="left", fontsize=8, bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 3})
    ax.set_xlim(RATIO_MIN - 0.02, RATIO_MAX + 0.02)
    ax.set_ylim(DEPTH_MIN - 4, DEPTH_MAX + 4)
    add_two_qubit_gate_reference_line(ax)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", framealpha=0.92, fontsize=8)

    out_dir = run_dir / "images"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "monotone_bisection_gp_levelset_shadow.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Saved GP contour plot: {out_path}")
    return out_path


def plot_gp_levelsets_colored(crossings_by_seed: dict[int, list[CrossingRecord]], records_by_seed: dict[int, list[EvalRecord]], run_dir: Path) -> Path:
    """Plot fidelity level sets with measured evaluations colored by fidelity value.
    
    Each contour of constant fidelity is shown as a band of the same color,
    making it easy to see level-set structure across (ratio, depth) space.
    The target contour (fidelity=TARGET) is overlaid as a red line.
    """
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    total_rb_runs = sum(record.runs for records in records_by_seed.values() for record in records)
    
    all_depth = []
    all_ratio = []
    all_fidelity = []
    all_std = []
    for records in records_by_seed.values():
        for record in records:
            all_depth.append(record.depth)
            all_ratio.append(record.ratio)
            all_fidelity.append(record.fidelity)
            all_std.append(record.fidelity_std)

    if not all_ratio:
        plt.close(fig)
        return None

    all_depth = np.array(all_depth)
    all_ratio = np.array(all_ratio)
    all_fidelity = np.array(all_fidelity)
    all_std = np.array(all_std)

    size = 28 + 260 * np.clip(all_std, 0.0, 0.25)
    target_band = (all_fidelity >= TARGET - 0.1) & (all_fidelity <= TARGET + 0.1)
    edgecolors = np.where(target_band, "black", "none")
    linewidths = np.where(target_band, 1.25, 0.0)

    scatter = ax.scatter(
        all_ratio,
        all_depth,
        c=all_fidelity,
        s=size,
        edgecolors=edgecolors,
        linewidths=linewidths,
        alpha=0.75,
        cmap=PLOT_CMAPS["levelsets_colored"],
        vmin=0.3,
        vmax=0.9,
    )
    
    seed_colors = {seed: plt.get_cmap("tab10")(idx % 10) for idx, seed in enumerate(sorted(crossings_by_seed))}

    for seed, crossings in sorted(crossings_by_seed.items()):
        crossings_sorted = sorted(crossings, key=lambda c: c.ratio)
        ratios = np.array([c.ratio for c in crossings_sorted], dtype=float)
        depths = np.array([c.crossing_depth for c in crossings_sorted], dtype=float)
        depth_stds = np.array([c.crossing_depth_std for c in crossings_sorted], dtype=float)
        
        ax.plot(ratios, depths, color=seed_colors[seed], alpha=0.9, linewidth=2.5, label=f"seed {seed} target contour (fidelity={TARGET})", zorder=3)
        ax.fill_between(ratios, depths - depth_stds, depths + depth_stds, color=seed_colors[seed], alpha=0.14, linewidth=0, label=f"seed {seed} depth uncertainty")
    
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("observed fidelity")
    
    ax.set_xlabel("two_qubit_ratio")
    ax.set_ylabel("depth")
    ax.set_title("Fidelity level sets: points colored by observed fidelity value")
    ax.text(0.02, 0.98, f"total RB runs: {total_rb_runs}", transform=ax.transAxes, va="top", ha="left", fontsize=8, bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 3})
    ax.set_xlim(RATIO_MIN - 0.02, RATIO_MAX + 0.02)
    ax.set_ylim(DEPTH_MIN - 4, DEPTH_MAX + 4)
    add_two_qubit_gate_reference_line(ax)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", framealpha=0.92, fontsize=8)

    out_dir = run_dir / "images"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "monotone_bisection_gp_levelsets_colored.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Saved colored level-set plot: {out_path}")
    return out_path


def plot_gp_evaluations(records_by_seed: dict[int, list[EvalRecord]], run_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    total_rb_runs = sum(record.runs for records in records_by_seed.values() for record in records)
    all_depth = []
    all_ratio = []
    all_fidelity = []
    all_std = []

    for records in records_by_seed.values():
        for record in records:
            all_depth.append(record.depth)
            all_ratio.append(record.ratio)
            all_fidelity.append(record.fidelity)
            all_std.append(record.fidelity_std)

    all_depth = np.array(all_depth)
    all_ratio = np.array(all_ratio)
    all_fidelity = np.array(all_fidelity)
    all_std = np.array(all_std)

    size = 24 + 240 * np.clip(all_std, 0.0, 0.25)

    scatter = ax.scatter(
        all_ratio,
        all_depth,
        c=all_fidelity,
        s=size,
        edgecolors="black",
        linewidths=0.3,
        alpha=0.85,
        cmap=PLOT_CMAPS["evaluations"],
        vmin=0.0,
        vmax=1.0,
    )
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("observed fidelity")

    ax.set_xlabel("two_qubit_ratio")
    ax.set_ylabel("depth")
    ax.set_title("GP-level-set measurement locations; marker size ~ fidelity uncertainty")
    ax.text(0.02, 0.98, f"total RB runs: {total_rb_runs}", transform=ax.transAxes, va="top", ha="left", fontsize=8, bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none", "pad": 3})
    ax.set_xlim(RATIO_MIN - 0.02, RATIO_MAX + 0.02)
    ax.set_ylim(DEPTH_MIN - 4, DEPTH_MAX + 4)
    add_two_qubit_gate_reference_line(ax)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", framealpha=0.92, fontsize=8)

    out_dir = run_dir / "images"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "monotone_bisection_gp_levelset_evaluations.png"
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(f"Saved evaluation plot: {out_path}")
    return out_path


def run_all() -> Path:
    run_dir = RESULTS_DIR / (datetime.now().strftime("%Y%m%d_%H%M%S") + "_gp_levelset")
    run_dir.mkdir(parents=True, exist_ok=True)

    metadata = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "method": "gp_levelset_estimator",
        "target": TARGET,
        "study_seeds": list(STUDY_SEEDS),
        "n_qubits": N_QUBITS,
        "depth_range": [DEPTH_MIN, DEPTH_MAX],
        "ratio_range": [RATIO_MIN, RATIO_MAX],
        "n_ratio_points": N_RATIO_POINTS,
        "initial_ratio_points": INITIAL_RATIO_POINTS,
        "initial_depth_points": INITIAL_DEPTH_POINTS,
        "n_adaptive_trials": N_ADAPTIVE_TRIALS,
        "runs_per_config": RUNS_PER_CONFIG,
        "evaluate_endpoints": EVALUATE_ENDPOINTS,
        "target_overlap_sigma": TARGET_OVERLAP_SIGMA,
        "approx_total_rb_runs": (
            len(STUDY_SEEDS) * (len(np.unique(np.array(make_initial_design(), dtype=float), axis=0)) + N_ADAPTIVE_TRIALS) * RUNS_PER_CONFIG
        ),
        "single_qubit_error": SINGLE_QUBIT_ERROR,
        "two_qubit_gate_error": TWO_QUBIT_GATE_ERROR,
        "assumption": "fidelity is smooth in ratio and depth and measurement uncertainty is heteroscedastic",
        "uncertainty_note": (
            "Fidelity uncertainty is baked into the GP observation noise via fidelity_std. "
            "Depth error bars are estimated from the posterior fidelity std and local slope." 
        ),
    }
    write_json(run_dir / "run_config.json", metadata)
    print(f"Saved run config: {run_dir / 'run_config.json'}")

    records_by_seed: dict[int, list[EvalRecord]] = {}
    crossings_by_seed: dict[int, list[CrossingRecord]] = {}

    for seed in STUDY_SEEDS:
        records, crossings = run_gp_levelset(seed, run_dir)
        records_by_seed[seed] = records
        crossings_by_seed[seed] = crossings

        write_json(run_dir / "all_records_by_seed.json", {str(s): [asdict(r) for r in recs] for s, recs in records_by_seed.items()})
        print(f"[seed {seed}] Updated aggregate records")
        write_json(run_dir / "crossings_by_seed.json", {str(s): [asdict(c) for c in crs] for s, crs in crossings_by_seed.items()})
        print(f"[seed {seed}] Updated aggregate crossings")

    plot_gp_levelset(crossings_by_seed, records_by_seed, run_dir)
    plot_gp_levelsets_colored(crossings_by_seed, records_by_seed, run_dir)
    plot_gp_evaluations(records_by_seed, run_dir)

    print("\nSaved")
    print("-----")
    print(f"run folder: {run_dir}")
    return run_dir


def main() -> None:
    run_all()


if __name__ == "__main__":
    main()
