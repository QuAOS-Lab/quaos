"""Hybrid contour-first + GP-uncertainty boundary estimator.

Goal
----
Estimate the fidelity = 0.5 boundary in (two_qubit_ratio, depth), with an
uncertainty band on the crossing depth.

Design
------
This is deliberately a hybrid of the two approaches we discussed:

1. Use the contour-first method from viarregio5 to cheaply get onto the
   fid = 0.5 boundary and trace it.
2. Replace the final "dumb" refinement with a GP-style, uncertainty-aware
   refinement. Candidate points are scored by:

       probability of being close to fid = 0.5
       × uncertainty / cost

   where closeness is normalized by total uncertainty. Therefore, for example,
   fid = 0.60 ± 0.20 is treated as more plausibly boundary-relevant than
   fid = 0.60 ± 0.01.
3. Fit a GP to all measured data using heteroscedastic observation noise from
   the fidelity estimator variance.
4. Extract the final fid = TARGET contour and convert posterior fidelity
   uncertainty into depth uncertainty via local slope.

Assumptions
-----------
This file is meant to sit beside your existing viarregio files. It imports the
contour-first search utilities from viarregio5 and the monotone-surface/cost
utilities from viarregio4.

If your contour-first file has a different module name, change the viarregio5
import block below.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import warnings
import json
import os
from typing import Iterable

SCRIPT_DIR = Path(__file__).resolve().parent

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from numpy.random import default_rng

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.isotonic import IsotonicRegression

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData

from viarregio2 import (
    config_from_parameters,
    fidelity_mean,
    fidelity_variance,
    make_backend,
    total_measurements,
)
from viarregio4 import (
    affordable_shot_count,
    hqc_cost,
    plot_monotone_fidelity_surface_with_confidence,
    print_experiment_summary,
    print_fit_reports,
)

# This is the contour-first file. Rename this import if your file has another name.
import viarregio5 as _viarregio5
from viarregio5 import (
    BudgetState,
    ContourFirstExperimentConfig,
    confirm_initial_anchor,
    find_initial_anchor,
    find_depth_crossing_at_ratio,
    spend_config,
    template_config,
    trace_from_anchor,
)


TARGET = 0.5


@dataclass(frozen=True)
class HybridContourGPConfig(ContourFirstExperimentConfig):
    """Settings for contour-first search plus GP uncertainty refinement."""

    # GP refinement budget strategy. These shots are NOT extra by design; they
    # consume the same remaining measurement/HQC budget after the contour trace.
    gp_refine_after_trace: bool = True
    # Reserve this fraction of the HQC budget for the GP uncertainty refinement stage.
    # With the default 0.25, contour search can use at most 75% of HQC.
    gp_reserved_hqc_fraction: float = 0.25
    gp_refinement_shots: int = 2
    gp_max_refinement_executions: int | None = None

    # Candidate grid around the traced contour anchors.
    gp_candidate_grid_size: tuple[int, int] = (80, 80)  # depth grid, ratio grid
    gp_local_depth_radius_fraction: float = 0.10

    # Ask GP refinement to focus near measured contour-like points.
    # A point is contour-like if measured fidelity is within TARGET +/- band.
    gp_focus_on_measured_target_band: bool = True
    gp_focus_target_band_width: float = 0.20

    # Candidate band around those measured near-target points.
    gp_focus_ratio_radius: float = 0.06
    gp_focus_depth_radius_fraction: float = 0.12

    # Before GP refinement, force a few fixed-ratio crossing searches across
    # the full ratio range. This prevents the GP from extrapolating a boundary
    # from one low-ratio cluster only.
    force_ratio_coverage_before_gp: bool = True
    forced_ratio_coverage_count: int = 6
    forced_ratio_min_separation: float = 0.04

    # Do not spend contour-search budget on an exact config already measured.
    contour_skip_existing_configs: bool = False

    # viarregio2.config_from_parameters rounds ratios to two decimal places.
    # Keep contour-trace moves at least this large so a projected target like
    # 0.204 does not collapse back to the already-measured 0.20 anchor.
    contour_trace_min_ratio_step: float = 0.01

    # Cap contour-stage measurements per ratio bin. This prevents the first
    # successful ratio from consuming the whole contour-search budget.
    # None means no cap. This cap does not apply to GP uncertainty refinement.
    max_contour_executions_per_ratio: int | None = None
    contour_ratio_bin_width: float = 0.01

    # Plot only the contour segments supported by measured ratios nearby.
    plot_only_supported_ratios: bool = True
    plot_support_ratio_radius: float = 0.08
    gp_include_existing_configs: bool = True

    # Score = boundary_probability * uncertainty / cost.
    # boundary_probability uses z = |mu - target| / total_sigma.
    gp_target_width_floor: float = 0.03
    gp_measurement_variance_floor: float = 1e-4
    gp_cost_power: float = 1.0
    gp_min_fit_points: int = 8

    # GP model details.
    gp_depth_length_scale: float = 30.0
    gp_ratio_length_scale: float = 0.15
    gp_n_restarts_optimizer: int = 4
    # If False, do not let sklearn optimize the GP kernel hyperparameters.
    # This bypasses repeated ConvergenceWarning messages by using the supplied
    # gp_ratio_length_scale / gp_depth_length_scale directly.
    gp_optimize_kernel: bool = False
    gp_suppress_convergence_warnings: bool = True

    # Output.
    gp_ratio_eval_points: int = 41
    gp_depth_eval_points: int = 301

    # Final contour cleanup. The physical boundary depth should not increase
    # as the two-qubit ratio increases. Apply a weighted isotonic projection
    # to the extracted GP crossing depths before saving/plotting.
    enforce_monotone_crossing_depth: bool = False

    hybrid_save_path: str | Path | None = "viarregio6_hybrid_boundary.json"
    hybrid_output_dir: str | Path = "viarregio6_hybrid_outputs"


@dataclass
class HybridCrossingRecord:
    ratio: float
    crossing_depth: float
    crossing_depth_std: float
    gp_pred_mean: float
    gp_pred_std: float
    gp_depth_slope: float
    method: str


@dataclass
class HybridRefinementRecord:
    depth: int
    ratio: float
    predicted_fidelity: float
    model_std: float
    total_std: float
    z_to_target: float
    score: float
    shots: int
    cost: float


@dataclass
class MeasurementTrialRecord:
    trial_index: int
    phase: str
    n_qubits: int
    depth: int
    ratio: float
    requested_shots: int
    spent_shots: int
    total_runs_at_config: int
    fidelity: float | None
    fidelity_std: float | None
    abs_error_from_target: float | None
    remaining_measurements: int
    remaining_hqc: float


_ORIGINAL_VIARREGIO5_SPEND_CONFIG = _viarregio5.spend_config


def ensure_trial_log(budget: BudgetState) -> list[MeasurementTrialRecord]:
    if not hasattr(budget, "trial_log"):
        setattr(budget, "trial_log", [])
    return getattr(budget, "trial_log")


def ensure_contour_ratio_counts(budget: BudgetState) -> dict[float, int]:
    if not hasattr(budget, "contour_ratio_counts"):
        setattr(budget, "contour_ratio_counts", {})
    return getattr(budget, "contour_ratio_counts")


def contour_ratio_key(ratio: float, settings: HybridContourGPConfig) -> float:
    width = max(float(settings.contour_ratio_bin_width), 1e-12)
    return round(float(ratio) / width) * width


def is_contour_phase(phase: str) -> bool:
    return phase.startswith("contour") or phase in {
        "initial_anchor",
        "anchor_bisection",
        "local_crossing",
        "local_bisection",
        "crossing_confirmation",
        "trace_correction",
    }


def logged_spend_config(
    *,
    backend,
    rng,
    data: RMBData,
    config: RMBConfig,
    requested_shots: int,
    budget: BudgetState,
    settings: HybridContourGPConfig,
    phase: str = "contour_search",
) -> int:
    """Wrapper around viarregio5.spend_config that logs every measured config."""
    spent = _ORIGINAL_VIARREGIO5_SPEND_CONFIG(
        backend=backend,
        rng=rng,
        data=data,
        config=config,
        requested_shots=requested_shots,
        budget=budget,
        settings=settings,
    )

    if spent <= 0:
        return spent

    if (
        settings.max_contour_executions_per_ratio is not None
        and is_contour_phase(phase)
    ):
        counts = ensure_contour_ratio_counts(budget)
        r_key = contour_ratio_key(float(config.min_two_qubit_gate_ratio), settings)
        counts[r_key] = counts.get(r_key, 0) + 1

    trials = ensure_trial_log(budget)
    p = None
    p_std = None
    abs_err = None
    total_runs = 0
    if config in data and data[config].num_runs() > 0:
        total_runs = int(data[config].num_runs())
        p = float(fidelity_mean(data[config]))
        p_std = float(np.sqrt(max(fidelity_variance(data[config]), 0.0)))
        abs_err = float(abs(p - TARGET))

    record = MeasurementTrialRecord(
        trial_index=len(trials) + 1,
        phase=phase,
        n_qubits=int(config.n_qubits),
        depth=int(config.depth),
        ratio=float(config.min_two_qubit_gate_ratio),
        requested_shots=int(requested_shots),
        spent_shots=int(spent),
        total_runs_at_config=int(total_runs),
        fidelity=p,
        fidelity_std=p_std,
        abs_error_from_target=abs_err,
        remaining_measurements=int(budget.remaining_measurements),
        remaining_hqc=float(budget.remaining_hqc),
    )
    trials.append(record)

    if p is None:
        print(
            f"[trial {record.trial_index:03d} | {phase}] "
            f"n={record.n_qubits}, ratio={record.ratio:.4f}, depth={record.depth}, "
            f"shots={spent}/{requested_shots}, fidelity=None, "
            f"HQC_left={budget.remaining_hqc:.3f}"
        )
    else:
        print(
            f"[trial {record.trial_index:03d} | {phase}] "
            f"n={record.n_qubits}, ratio={record.ratio:.4f}, depth={record.depth}, "
            f"shots={spent}/{requested_shots}, runs_here={total_runs}, "
            f"fidelity={p:.4f} ± {p_std:.4f}, |err|={abs_err:.4f}, "
            f"HQC_left={budget.remaining_hqc:.3f}"
        )

    return spent


# Patch viarregio5 so imported contour-first routines also get logged.
_viarregio5.spend_config = logged_spend_config
spend_config = logged_spend_config


def print_trial_log(trials: list[MeasurementTrialRecord]) -> None:
    print()
    print("Trial / fidelity estimate log")
    print("-----------------------------")
    if not trials:
        print("No trials recorded.")
        return
    for r in trials:
        if r.fidelity is None:
            fid_text = "fidelity=None"
        else:
            fid_text = f"fidelity={r.fidelity:.4f} ± {r.fidelity_std:.4f}, |err|={r.abs_error_from_target:.4f}"
        print(
            f"[{r.trial_index:03d} | {r.phase}] "
            f"n={r.n_qubits}, ratio={r.ratio:.4f}, depth={r.depth}, "
            f"shots={r.spent_shots}/{r.requested_shots}, runs_here={r.total_runs_at_config}, "
            f"{fid_text}, HQC_left={r.remaining_hqc:.3f}"
        )


def save_trial_log(trials: list[MeasurementTrialRecord], out_dir: Path) -> None:
    write_json(out_dir / "hybrid_trial_log.json", [asdict(t) for t in trials])


def json_safe(obj):
    """Make run summaries robust to Path and NumPy objects."""
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return str(obj)


def write_json(path: str | Path, data: dict | list) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(data, f, indent=2, default=json_safe)


def config_depth_ratio(config: RMBConfig) -> tuple[float, float]:
    """Return coordinates as (depth, ratio)."""
    return float(config.depth), float(config.min_two_qubit_gate_ratio)


def measured_fidelity_std(data: RMBData, config: RMBConfig) -> float:
    if config not in data or data[config].num_runs() == 0:
        return np.inf
    return float(np.sqrt(max(fidelity_variance(data[config]), 0.0)))


def training_arrays_from_data(data: RMBData) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build GP training arrays from all measured configs.

    X columns are [ratio, depth], matching the older GP file.
    """
    X: list[list[float]] = []
    y: list[float] = []
    y_std: list[float] = []

    for config, estimator in data.items():
        if estimator.num_runs() <= 0:
            continue
        depth, ratio = config_depth_ratio(config)
        X.append([ratio, depth])
        y.append(float(fidelity_mean(estimator)))
        y_std.append(float(np.sqrt(max(fidelity_variance(estimator), 0.0))))

    if not X:
        return np.empty((0, 2)), np.empty((0,)), np.empty((0,))
    return np.asarray(X, dtype=float), np.asarray(y, dtype=float), np.asarray(y_std, dtype=float)


def measured_target_band_points(
    data: RMBData,
    settings: HybridContourGPConfig,
) -> list[tuple[float, float]]:
    """Return measured (ratio, depth) points with fidelity within TARGET +/- band."""
    points: list[tuple[float, float]] = []
    band = float(settings.gp_focus_target_band_width)

    for config, estimator in data.items():
        if estimator.num_runs() <= 0:
            continue

        p = float(fidelity_mean(estimator))
        if abs(p - TARGET) <= band:
            points.append(
                (
                    float(config.min_two_qubit_gate_ratio),
                    float(config.depth),
                )
            )

    return points


def expected_two_qubit_gates_from_depth_ratio(
    *,
    depth: float,
    ratio: float,
    n_qubits: int,
) -> float:
    """Expected/requested two-qubit gate count for a depth-ratio point.

    This matches the convention used elsewhere in the RMB scripts:
        requested = int(ratio * n_qubits * depth) // 2
    For plotting a smooth contour, keep the continuous value.
    """
    return float(ratio * n_qubits * depth / 2.0)


def representative_n_qubits(data: RMBData, settings: HybridContourGPConfig) -> int:
    for config, estimator in data.items():
        if estimator.num_runs() > 0:
            return int(config.n_qubits)
    return int(settings.n_qubits_values[0])


def build_gp_from_data(data: RMBData, settings: HybridContourGPConfig) -> GaussianProcessRegressor:
    X, y, y_std = training_arrays_from_data(data)
    if len(y) < settings.gp_min_fit_points:
        raise RuntimeError(
            f"Need at least {settings.gp_min_fit_points} measured points for GP fit; got {len(y)}."
        )

    kernel = ConstantKernel(1.0, constant_value_bounds=(1e-3, 1e3)) * Matern(
        length_scale=[settings.gp_ratio_length_scale, settings.gp_depth_length_scale],
        length_scale_bounds=[(1e-3, 10.0), (1.0, 1e3)],
        nu=2.5,
    ) + WhiteKernel(noise_level=1e-6, noise_level_bounds=(1e-8, 1e0))

    alpha = np.maximum(y_std**2, settings.gp_measurement_variance_floor)
    optimizer = "fmin_l_bfgs_b" if settings.gp_optimize_kernel else None
    n_restarts = settings.gp_n_restarts_optimizer if settings.gp_optimize_kernel else 0
    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=alpha,
        normalize_y=True,
        optimizer=optimizer,
        n_restarts_optimizer=n_restarts,
    )

    if settings.gp_suppress_convergence_warnings:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", ConvergenceWarning)
            gp.fit(X, y)
    else:
        gp.fit(X, y)
    return gp


def candidate_configs_near_anchors(
    *,
    anchors: Iterable[RMBConfig],
    data: RMBData,
    settings: HybridContourGPConfig,
) -> list[RMBConfig]:
    """Generate candidate configs near traced anchors or measured near-target points.

    If gp_focus_on_measured_target_band=True, the candidate band is built around
    measured points whose fidelity is already close to TARGET, e.g. TARGET +/- 0.2.
    """
    anchors = list(anchors)

    if not anchors and not data:
        return []

    n_depth, n_ratio = settings.gp_candidate_grid_size
    d_min, d_max = settings.depth_bounds
    r_min, r_max = settings.ratio_bounds
    depth_span = float(d_max - d_min)

    if anchors:
        n_qubits = int(anchors[0].n_qubits)
    else:
        n_qubits = int(settings.n_qubits_values[0])

    template = template_config(settings, n_qubits)

    anchor_points = sorted(
        [(float(a.min_two_qubit_gate_ratio), float(a.depth)) for a in anchors],
        key=lambda x: x[0],
    )

    if settings.gp_focus_on_measured_target_band:
        band_points = measured_target_band_points(data, settings)
        if band_points:
            anchor_points = sorted(band_points, key=lambda x: x[0])

    if not anchor_points:
        return []

    anchor_ratios = np.asarray([p[0] for p in anchor_points], dtype=float)
    anchor_depths = np.asarray([p[1] for p in anchor_points], dtype=float)

    ratio_low = max(
        float(r_min),
        float(anchor_ratios.min()) - settings.gp_focus_ratio_radius,
    )
    ratio_high = min(
        float(r_max),
        float(anchor_ratios.max()) + settings.gp_focus_ratio_radius,
    )

    if ratio_high <= ratio_low:
        return []

    ratio_grid = np.linspace(ratio_low, ratio_high, n_ratio)

    if settings.gp_focus_on_measured_target_band:
        local_radius = settings.gp_focus_depth_radius_fraction * depth_span
    else:
        local_radius = settings.gp_local_depth_radius_fraction * depth_span

    configs: dict[tuple[int, float], RMBConfig] = {}

    for ratio in ratio_grid:
        center_depth = float(np.interp(ratio, anchor_ratios, anchor_depths))
        depth_low = max(float(d_min), center_depth - local_radius)
        depth_high = min(float(d_max), center_depth + local_radius)

        for depth in np.linspace(depth_low, depth_high, n_depth):
            config = config_from_parameters(
                template=template,
                depth=float(depth),
                ratio=float(ratio),
            )
            key = (
                int(config.depth),
                round(float(config.min_two_qubit_gate_ratio), 8),
            )
            configs[key] = config

    return list(configs.values())


def estimate_prospective_measurement_std(mu: float, shots: int, settings: HybridContourGPConfig) -> float:
    """Approximate future Bernoulli/Bayesian measurement uncertainty.

    This is only for acquisition scoring before the candidate is measured.
    After measurement, the real estimator variance is used in the GP alpha.
    """
    shots = max(1, int(shots))
    bern_var = max(mu * (1.0 - mu), settings.gp_measurement_variance_floor)
    return float(np.sqrt(bern_var / shots))


def score_candidate(
    *,
    gp: GaussianProcessRegressor,
    data: RMBData,
    config: RMBConfig,
    requested_shots: int,
    budget: BudgetState,
    settings: HybridContourGPConfig,
) -> HybridRefinementRecord | None:
    if requested_shots <= 0 or not budget.can_spend():
        return None

    shots = min(requested_shots, budget.remaining_measurements)
    shots = affordable_shot_count(
        config=config,
        requested_shots=shots,
        remaining_hqc=budget.remaining_hqc,
        settings=settings,
    )
    if shots <= 0:
        return None

    current_runs = data[config].num_runs() if config in data else 0
    if current_runs >= settings.max_shots_per_config:
        return None

    ratio = float(config.min_two_qubit_gate_ratio)
    depth = float(config.depth)
    mu, model_std = gp.predict(np.array([[ratio, depth]], dtype=float), return_std=True)
    mu = float(mu[0])
    model_std = float(model_std[0])

    if config in data and data[config].num_runs() > 0:
        measurement_std = measured_fidelity_std(data, config)
    else:
        measurement_std = estimate_prospective_measurement_std(mu, shots, settings)

    target_width = max(settings.gp_target_width_floor, 1e-8)
    total_std = float(np.sqrt(model_std**2 + measurement_std**2 + target_width**2))
    z = abs(mu - TARGET) / total_std

    # This is the key idea Shreya wanted:
    # closeness to 0.5 is judged relative to uncertainty.
    boundary_probability = float(np.exp(-0.5 * z**2))

    # Prefer points that are both plausibly on the boundary and still uncertain.
    information_value = boundary_probability * total_std

    cost = hqc_cost(config, shots, settings) if settings.hqc_budget is not None else float(shots)
    cost = max(float(cost), 1e-12)
    score = information_value / (cost**settings.gp_cost_power)

    return HybridRefinementRecord(
        depth=int(config.depth),
        ratio=float(config.min_two_qubit_gate_ratio),
        predicted_fidelity=mu,
        model_std=model_std,
        total_std=total_std,
        z_to_target=float(z),
        score=float(score),
        shots=int(shots),
        cost=float(cost),
    )


def force_ratio_coverage_anchors(
    *,
    backend,
    rng,
    data: RMBData,
    anchors: list[RMBConfig],
    settings: HybridContourGPConfig,
    budget: BudgetState,
) -> list[RMBConfig]:
    """Force fixed-ratio crossing searches across the full ratio range.

    This prevents a failure mode where all measured data sit near the first
    low-ratio anchor and the GP extrapolates the rest of the fid=0.5 line.
    """
    if not settings.force_ratio_coverage_before_gp or not budget.can_spend():
        return []

    n_qubits = int(anchors[0].n_qubits) if anchors else int(settings.n_qubits_values[0])
    template = template_config(settings, n_qubits)
    r_min, r_max = settings.ratio_bounds
    target_ratios = np.linspace(float(r_min), float(r_max), max(2, settings.forced_ratio_coverage_count))

    existing = np.array([float(a.min_two_qubit_gate_ratio) for a in anchors], dtype=float) if anchors else np.array([])
    added: list[RMBConfig] = []

    for ratio in target_ratios:
        if not budget.can_spend():
            break
        if existing.size and np.min(np.abs(existing - ratio)) < settings.forced_ratio_min_separation:
            continue

        if settings.verbose:
            print(f"\nForced ratio-coverage crossing search at ratio={ratio:.4f}")

        anchor = find_depth_crossing_at_ratio(
            backend=backend,
            rng=rng,
            data=data,
            template=template,
            ratio=float(ratio),
            settings=settings,
            budget=budget,
            local=False,
        )
        if anchor is not None:
            added.append(anchor)
            existing = np.append(existing, float(anchor.min_two_qubit_gate_ratio))
            if settings.verbose:
                p = fidelity_mean(data[anchor]) if anchor in data else float("nan")
                print(
                    "  added coverage anchor: "
                    f"ratio={anchor.min_two_qubit_gate_ratio:.4f}, "
                    f"depth={anchor.depth}, fidelity={p:.4f}"
                )
        elif settings.verbose:
            print(f"  no crossing found at ratio={ratio:.4f}")

    return added


def nudged_trace_ratio(
    *,
    target_ratio: float,
    current_ratio: float,
    direction: int,
    settings: HybridContourGPConfig,
) -> float | None:
    """Return a ratio that really advances after config rounding."""
    if direction == 0:
        return None

    r_min, r_max = settings.ratio_bounds
    min_step = max(float(settings.contour_trace_min_ratio_step), 0.0)
    ratio = float(np.clip(target_ratio, r_min, r_max))

    rounded_current = round(float(current_ratio), 2)
    rounded_target = round(ratio, 2)
    if direction > 0 and rounded_target <= rounded_current:
        ratio = float(current_ratio) + min_step
    elif direction < 0 and rounded_target >= rounded_current:
        ratio = float(current_ratio) - min_step

    ratio = float(np.clip(ratio, r_min, r_max))
    if round(ratio, 2) == rounded_current:
        return None
    return ratio


def trace_from_anchor(
    *,
    backend,
    rng,
    data: RMBData,
    anchor: RMBConfig,
    settings: HybridContourGPConfig,
    budget: BudgetState,
) -> list[RMBConfig]:
    """Trace the contour, forcing projected targets to move across ratio bins."""
    anchors = [anchor]
    for direction in settings.trace_directions:
        current = anchor
        while budget.can_spend():
            next_anchor = None
            n_attempts = max(1, settings.trace_step_shrink_attempts if settings.trace_local_crossing else 1)
            for attempt in range(n_attempts):
                step_fraction = settings.trace_ratio_step_fraction * (0.5 ** attempt)
                depth_guess, ratio = _viarregio5.projected_trace_target(
                    data,
                    current,
                    direction,
                    settings,
                    step_fraction=step_fraction,
                )
                ratio = nudged_trace_ratio(
                    target_ratio=float(ratio),
                    current_ratio=float(current.min_two_qubit_gate_ratio),
                    direction=int(direction),
                    settings=settings,
                )
                if ratio is None:
                    break

                if settings.verbose and abs(float(ratio) - float(current.min_two_qubit_gate_ratio)) >= 1e-9:
                    print(
                        "  contour trace moving: "
                        f"ratio {current.min_two_qubit_gate_ratio:.4f} -> {ratio:.4f}"
                    )

                next_anchor = _viarregio5.trace_projected_anchor(
                    backend=backend,
                    rng=rng,
                    data=data,
                    current=current,
                    depth_guess=depth_guess,
                    ratio=float(ratio),
                    direction=direction,
                    settings=settings,
                    budget=budget,
                )
                if next_anchor is not None:
                    break
            if next_anchor is None:
                break
            anchors.append(next_anchor)
            current = next_anchor
    return anchors


def gp_uncertainty_refinement(
    *,
    backend,
    rng,
    data: RMBData,
    anchors: list[RMBConfig],
    settings: HybridContourGPConfig,
    budget: BudgetState,
) -> list[HybridRefinementRecord]:
    """Spend remaining budget using File-2-style uncertainty-aware scoring."""
    if not settings.gp_refine_after_trace or not budget.can_spend():
        return []

    history: list[HybridRefinementRecord] = []
    executions = 0

    while budget.can_spend():
        if settings.gp_max_refinement_executions is not None and executions >= settings.gp_max_refinement_executions:
            break

        try:
            gp = build_gp_from_data(data, settings)
        except RuntimeError:
            break

        candidates = candidate_configs_near_anchors(
            anchors=anchors,
            data=data,
            settings=settings,
        )
        if settings.gp_include_existing_configs:
            candidates.extend([c for c, est in data.items() if est.num_runs() > 0])

        best_config: RMBConfig | None = None
        best_record: HybridRefinementRecord | None = None

        seen: set[tuple[int, float]] = set()
        for candidate in candidates:
            key = (int(candidate.depth), float(candidate.min_two_qubit_gate_ratio))
            if key in seen:
                continue
            seen.add(key)

            record = score_candidate(
                gp=gp,
                data=data,
                config=candidate,
                requested_shots=settings.gp_refinement_shots,
                budget=budget,
                settings=settings,
            )
            if record is None:
                continue
            if best_record is None or record.score > best_record.score:
                best_record = record
                best_config = candidate

        if best_config is None or best_record is None or best_record.score <= 0.0:
            break

        spent = spend_config(
            backend=backend,
            rng=rng,
            data=data,
            config=best_config,
            requested_shots=best_record.shots,
            budget=budget,
            settings=settings,
            phase="gp_uncertainty_refinement",
        )
        if spent <= 0:
            break

        # Store the pre-measurement score for auditability.
        history.append(best_record)
        executions += 1

    return history


def extract_gp_crossings(
    *,
    gp: GaussianProcessRegressor,
    settings: HybridContourGPConfig,
) -> list[HybridCrossingRecord]:
    d_min, d_max = settings.depth_bounds
    r_min, r_max = settings.ratio_bounds
    depth_grid = np.linspace(float(d_min), float(d_max), settings.gp_depth_eval_points)
    ratio_grid = np.linspace(float(r_min), float(r_max), settings.gp_ratio_eval_points)

    records: list[HybridCrossingRecord] = []
    for ratio in ratio_grid:
        X = np.array([[ratio, d] for d in depth_grid], dtype=float)
        mean_curve, std_curve = gp.predict(X, return_std=True)

        centered = mean_curve - TARGET
        sign = np.sign(centered)
        crossing_indices = np.where(sign[:-1] * sign[1:] <= 0)[0]

        if len(crossing_indices) == 0:
            closest_idx = int(np.argmin(np.abs(centered)))
            crossing_depth = float(depth_grid[closest_idx])
            pred_mu = float(mean_curve[closest_idx])
            pred_std = float(std_curve[closest_idx])
            slope = 0.0
            depth_std = float(max(0.5 * (d_max - d_min), 1.0))
            method = "closest_no_bracket"
        else:
            idx = int(crossing_indices[0])
            d0, d1 = float(depth_grid[idx]), float(depth_grid[idx + 1])
            f0, f1 = float(mean_curve[idx]), float(mean_curve[idx + 1])
            if abs(f1 - f0) < 1e-12:
                crossing_depth = 0.5 * (d0 + d1)
            else:
                crossing_depth = d0 + (TARGET - f0) * (d1 - d0) / (f1 - f0)

            pred_mu, pred_std = gp.predict(
                np.array([[ratio, crossing_depth]], dtype=float),
                return_std=True,
            )
            pred_mu = float(pred_mu[0])
            pred_std = float(pred_std[0])

            delta = max(1.0, 0.005 * (d_max - d_min))
            dlo = max(float(d_min), crossing_depth - delta)
            dhi = min(float(d_max), crossing_depth + delta)
            if dhi <= dlo:
                slope = 0.0
            else:
                vals, _ = gp.predict(np.array([[ratio, dlo], [ratio, dhi]], dtype=float), return_std=True)
                slope = float((vals[1] - vals[0]) / (dhi - dlo))

            grad = max(abs(slope), 1e-3)
            depth_std = float(max(pred_std / grad, 0.5))
            method = "gp_mean_root"

        records.append(
            HybridCrossingRecord(
                ratio=float(ratio),
                crossing_depth=float(crossing_depth),
                crossing_depth_std=float(depth_std),
                gp_pred_mean=float(pred_mu),
                gp_pred_std=float(pred_std),
                gp_depth_slope=float(slope),
                method=method,
            )
        )

    return records



def enforce_monotone_crossing_depths(
    crossings: list[HybridCrossingRecord],
    settings: HybridContourGPConfig,
) -> list[HybridCrossingRecord]:
    """Project extracted crossing depths onto a non-increasing curve in ratio.

    The GP is unconstrained, so extracting one root independently at each ratio
    can produce a jagged/non-monotone boundary. For RB-like decay, increasing
    the two-qubit ratio should not require a larger depth to hit the same
    fidelity target. This applies a weighted isotonic regression to the final
    depth-vs-ratio curve. It does not change the raw measured data.
    """
    if not settings.enforce_monotone_crossing_depth or len(crossings) < 2:
        return crossings

    ordered = sorted(crossings, key=lambda c: c.ratio)
    ratios = np.asarray([c.ratio for c in ordered], dtype=float)
    depths = np.asarray([c.crossing_depth for c in ordered], dtype=float)
    stds = np.asarray([max(c.crossing_depth_std, 1e-6) for c in ordered], dtype=float)
    weights = 1.0 / (stds**2)

    iso = IsotonicRegression(
        increasing=False,
        y_min=float(settings.depth_bounds[0]),
        y_max=float(settings.depth_bounds[1]),
        out_of_bounds="clip",
    )
    mono_depths = iso.fit_transform(ratios, depths, sample_weight=weights)

    cleaned: list[HybridCrossingRecord] = []
    for c, d in zip(ordered, mono_depths):
        method = c.method
        if abs(float(d) - float(c.crossing_depth)) > 1e-9:
            method = method + "+monotone_isotonic"
        cleaned.append(
            HybridCrossingRecord(
                ratio=float(c.ratio),
                crossing_depth=float(d),
                crossing_depth_std=float(c.crossing_depth_std),
                gp_pred_mean=float(c.gp_pred_mean),
                gp_pred_std=float(c.gp_pred_std),
                gp_depth_slope=float(c.gp_depth_slope),
                method=method,
            )
        )
    return cleaned


def plot_hybrid_gp_boundary(
    *,
    data: RMBData,
    crossings: list[HybridCrossingRecord],
    settings: HybridContourGPConfig,
) -> Path:
    out_dir = Path(settings.hybrid_output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "hybrid_gp_fid_0p5_boundary.png"

    fig, ax = plt.subplots(figsize=(8.0, 5.6))

    X, y, y_std = training_arrays_from_data(data)
    n_qubits = representative_n_qubits(data, settings)
    if len(y) > 0:
        sizes = 24 + 240 * np.clip(y_std, 0.0, 0.25)
        measured_two_qubit_gates = np.asarray(
            [expected_two_qubit_gates_from_depth_ratio(depth=depth, ratio=ratio, n_qubits=n_qubits) for ratio, depth in X],
            dtype=float,
        )
        scatter = ax.scatter(
            X[:, 0],
            measured_two_qubit_gates,
            c=y,
            s=sizes,
            edgecolors="black",
            linewidths=0.3,
            alpha=0.75,
            cmap="viridis",
            vmin=0.0,
            vmax=1.0,
            label="measured evaluations",
        )
        colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
        colorbar.set_label("observed fidelity")

    ratios = np.asarray([c.ratio for c in crossings], dtype=float)
    depths = np.asarray([c.crossing_depth for c in crossings], dtype=float)
    depth_stds = np.asarray([c.crossing_depth_std for c in crossings], dtype=float)

    measured_ratios = X[:, 0] if len(y) > 0 else np.asarray([], dtype=float)
    if settings.plot_only_supported_ratios and measured_ratios.size > 0:
        measured_ratio_arr = np.asarray(measured_ratios, dtype=float)
        supported = np.asarray([
            np.min(np.abs(measured_ratio_arr - r)) <= settings.plot_support_ratio_radius
            for r in ratios
        ])
        if np.any(supported):
            ratios = ratios[supported]
            depths = depths[supported]
            depth_stds = depth_stds[supported]

    two_qubit_gates = np.asarray(
        [expected_two_qubit_gates_from_depth_ratio(depth=d, ratio=r, n_qubits=n_qubits) for r, d in zip(ratios, depths)],
        dtype=float,
    )
    two_qubit_gate_stds = np.asarray([r * n_qubits * ds / 2.0 for r, ds in zip(ratios, depth_stds)], dtype=float)

    ax.plot(ratios, two_qubit_gates, color="black", linewidth=2.0, label=f"fid = {TARGET} boundary")
    ax.fill_between(
        ratios,
        two_qubit_gates - two_qubit_gate_stds,
        two_qubit_gates + two_qubit_gate_stds,
        color="gray",
        alpha=0.25,
        linewidth=0,
        label="GP two-qubit-gate uncertainty",
    )

    ax.set_xlabel("two_qubit_ratio")
    ax.set_ylabel("expected two-qubit gate count at fidelity = 0.5")
    ax.set_title("Hybrid contour-first + GP uncertainty boundary estimate")
    ax.set_xlim(settings.ratio_bounds[0] - 0.02, settings.ratio_bounds[1] + 0.02)
    y_upper = max(float(np.max(two_qubit_gates + two_qubit_gate_stds)) if len(two_qubit_gates) else 1.0, 1.0)
    ax.set_ylim(0, y_upper * 1.08)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", framealpha=0.92, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def plot_hybrid_gp_boundary_depth_ratio(
    *,
    data: RMBData,
    crossings: list[HybridCrossingRecord],
    settings: HybridContourGPConfig,
    trials: list[MeasurementTrialRecord] | None = None,
) -> Path:
    """Plot the same fid=TARGET boundary in ratio-vs-depth coordinates.

    This is the diagnostic companion to the two-qubit-gate-count plot.
    It shows whether the actual measured evaluations cover the ratio range,
    and it makes it easier to see whether the crossing depth is monotone.
    """
    out_dir = Path(settings.hybrid_output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "hybrid_gp_fid_0p5_boundary_depth_ratio.png"

    fig, ax = plt.subplots(figsize=(8.0, 5.6))

    X, y, y_std = training_arrays_from_data(data)
    if len(y) > 0:
        sizes = 24 + 240 * np.clip(y_std, 0.0, 0.25)
        process_by_point: dict[tuple[int, float], str] = {}
        for trial in trials or []:
            if trial.spent_shots <= 0:
                continue
            key = (int(trial.depth), round(float(trial.ratio), 12))
            process = "gp_refinement" if trial.phase == "gp_uncertainty_refinement" else "contour"
            if process == "gp_refinement" or key not in process_by_point:
                process_by_point[key] = process

        processes = np.asarray([
            process_by_point.get((int(round(depth)), round(float(ratio), 12)), "contour")
            for ratio, depth in X
        ])
        scatter = None
        for process, marker, label in (
            ("contour", "o", "contour-stage evaluations"),
            ("gp_refinement", "s", "GP-refinement evaluations"),
        ):
            mask = processes == process
            if not np.any(mask):
                continue
            scatter = ax.scatter(
                X[mask, 1],
                X[mask, 0],
                c=y[mask],
                s=sizes[mask],
                marker=marker,
                edgecolors="none",
                alpha=0.78,
                cmap="viridis",
                vmin=0.0,
                vmax=1.0,
                label=label,
            )
        near_target = np.abs(y - TARGET) <= 0.1
        if np.any(near_target):
            ax.scatter(
                X[near_target, 1],
                X[near_target, 0],
                s=sizes[near_target] + 70,
                marker="o",
                facecolors="none",
                edgecolors="black",
                linewidths=1.35,
                label=f"observed fidelity within +/-0.1 of {TARGET}",
            )
        if scatter is not None:
            colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
            colorbar.set_label("observed fidelity")

    ratios = np.asarray([c.ratio for c in crossings], dtype=float)
    depths = np.asarray([c.crossing_depth for c in crossings], dtype=float)
    depth_stds = np.asarray([c.crossing_depth_std for c in crossings], dtype=float)

    measured_ratios = X[:, 0] if len(y) > 0 else np.asarray([], dtype=float)
    if settings.plot_only_supported_ratios and measured_ratios.size > 0:
        measured_ratio_arr = np.asarray(measured_ratios, dtype=float)
        supported = np.asarray([
            np.min(np.abs(measured_ratio_arr - r)) <= settings.plot_support_ratio_radius
            for r in ratios
        ])
        if np.any(supported):
            ratios = ratios[supported]
            depths = depths[supported]
            depth_stds = depth_stds[supported]

    ax.plot(depths, ratios, color="black", linewidth=2.0, label=f"fid = {TARGET} boundary")
    ax.fill_betweenx(
        ratios,
        depths - depth_stds,
        depths + depth_stds,
        color="gray",
        alpha=0.25,
        linewidth=0,
        label="GP depth uncertainty",
    )

    ax.set_xlabel("depth at fidelity = 0.5")
    ax.set_ylabel("two_qubit_ratio")
    ax.set_title("Hybrid contour-first + GP uncertainty boundary estimate: depth view")
    ax.set_xlim(settings.depth_bounds[0] - 4, settings.depth_bounds[1] + 4)
    ax.set_ylim(settings.ratio_bounds[0] - 0.02, settings.ratio_bounds[1] + 0.02)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", framealpha=0.92, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def estimate_boundary_hybrid(settings: HybridContourGPConfig) -> RMB:
    """Run contour-first search, then GP uncertainty-aware refinement."""
    rng = default_rng(settings.rng_seed)
    backend = make_backend()
    rmb = RMB.default(rng).with_backend(backend)
    data: RMBData = rmb._data

    total_hqc_budget = (
        float(settings.hqc_budget)
        if settings.hqc_budget is not None
        else float(settings.measurement_budget)
    )
    reserve_fraction = float(np.clip(settings.gp_reserved_hqc_fraction, 0.0, 0.95))
    contour_hqc_budget = total_hqc_budget * (1.0 - reserve_fraction)

    # First stage: contour search sees only the contour budget.
    # After tracing, we restore the reserved HQC for GP uncertainty refinement.
    budget = BudgetState(
        remaining_measurements=settings.measurement_budget,
        remaining_hqc=contour_hqc_budget,
    )

    stop_reason = "budget not exhausted"
    batch_count = 0
    all_anchors: list[RMBConfig] = []
    refinement_history: list[HybridRefinementRecord] = []

    for n_qubits in settings.n_qubits_values:
        if not budget.can_spend():
            break

        anchor = find_initial_anchor(
            backend=backend,
            rng=rng,
            data=data,
            n_qubits=n_qubits,
            settings=settings,
            budget=budget,
        )
        if anchor is None:
            stop_reason = "no initial contour crossing found"
            continue

        anchor = confirm_initial_anchor(
            backend=backend,
            rng=rng,
            data=data,
            anchor=anchor,
            settings=settings,
            budget=budget,
        )

        if settings.verbose:
            print(
                "\nInitial contour anchor: "
                f"n_qubits={anchor.n_qubits}, depth={anchor.depth}, "
                f"ratio={anchor.min_two_qubit_gate_ratio:.4f}, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

        anchors = trace_from_anchor(
            backend=backend,
            rng=rng,
            data=data,
            anchor=anchor,
            settings=settings,
            budget=budget,
        )
        all_anchors.extend(anchors)
        batch_count += max(0, len(anchors) - 1)

        if settings.verbose:
            print(
                f"\nContour trace complete for n_qubits={n_qubits}: "
                f"{len(anchors)} anchors, measurements {total_measurements(data)} / "
                f"{settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

    if settings.force_ratio_coverage_before_gp and budget.can_spend():
        coverage_anchors = force_ratio_coverage_anchors(
            backend=backend,
            rng=rng,
            data=data,
            anchors=all_anchors,
            settings=settings,
            budget=budget,
        )
        if coverage_anchors:
            all_anchors.extend(coverage_anchors)
            batch_count += len(coverage_anchors)
            if settings.verbose:
                print(
                    f"\nForced ratio coverage complete: added {len(coverage_anchors)} anchors, "
                    f"measurements {total_measurements(data)} / {settings.measurement_budget}."
                )
                print_fit_reports(data, settings)

    if settings.hqc_budget is not None:
        contour_hqc_spent = contour_hqc_budget - budget.remaining_hqc
        budget.remaining_hqc = max(0.0, total_hqc_budget - contour_hqc_spent)
        if settings.verbose:
            print(
                "\nHQC reserve activated for GP refinement: "
                f"contour used {contour_hqc_spent:.3f} / {total_hqc_budget:.3f}; "
                f"available for refinement {budget.remaining_hqc:.3f}."
            )

    if settings.gp_refine_after_trace and budget.can_spend() and all_anchors:
        refinement_history = gp_uncertainty_refinement(
            backend=backend,
            rng=rng,
            data=data,
            anchors=all_anchors,
            settings=settings,
            budget=budget,
        )
        batch_count += len(refinement_history)
        if settings.verbose:
            print(
                f"\nGP uncertainty refinement complete: {len(refinement_history)} executions, "
                f"measurements {total_measurements(data)} / {settings.measurement_budget}."
            )
            print_fit_reports(data, settings)

    if budget.remaining_measurements <= 0:
        stop_reason = "measurement budget hit"
    elif settings.hqc_budget is not None and budget.remaining_hqc <= 0.0:
        stop_reason = "HQC budget hit"
    elif settings.hqc_budget is not None and budget.remaining_hqc <= settings.hqc_base_cost:
        stop_reason = "HQC budget effectively hit"
    elif stop_reason == "budget not exhausted":
        stop_reason = "contour traced/refined until no useful GP candidate or bounds reached"

    if settings.verbose:
        print_experiment_summary(
            data=data,
            settings=settings,
            stop_reason=stop_reason,
            batch_count=batch_count,
            circuit_executions=budget.circuit_executions,
            max_execution_repeats=budget.max_execution_repeats,
            remaining_measurements=budget.remaining_measurements,
            remaining_hqc=budget.remaining_hqc,
        )

    out_dir = Path(settings.hybrid_output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    trials = ensure_trial_log(budget)
    if settings.verbose:
        print_trial_log(trials)
    save_trial_log(trials, out_dir)
    if settings.verbose:
        print(f"Saved trial log: {out_dir / 'hybrid_trial_log.json'}")

    # Final GP crossing extraction and hybrid plot.
    try:
        gp = build_gp_from_data(data, settings)
        raw_crossings = extract_gp_crossings(gp=gp, settings=settings)
        write_json(out_dir / "hybrid_gp_crossings_raw.json", [asdict(c) for c in raw_crossings])
        crossings = enforce_monotone_crossing_depths(raw_crossings, settings)
        write_json(out_dir / "hybrid_gp_crossings.json", [asdict(c) for c in crossings])
        plot_path = plot_hybrid_gp_boundary(data=data, crossings=crossings, settings=settings)
        depth_plot_path = plot_hybrid_gp_boundary_depth_ratio(data=data, crossings=crossings, settings=settings, trials=trials)
        if settings.verbose:
            print(f"\nSaved hybrid GP two-qubit-gate boundary plot: {plot_path}")
            print(f"Saved hybrid GP depth-ratio boundary plot: {depth_plot_path}")
    except RuntimeError as exc:
        if settings.verbose:
            print(f"\nSkipping final GP crossing extraction: {exc}")

    write_json(out_dir / "hybrid_gp_refinement_history.json", [asdict(r) for r in refinement_history])
    write_json(
        out_dir / "hybrid_run_summary.json",
        {
            "target": TARGET,
            "measurements": total_measurements(data),
            "remaining_measurements": budget.remaining_measurements,
            "remaining_hqc": budget.remaining_hqc,
            "num_anchors": len(all_anchors),
            "num_gp_refinements": len(refinement_history),
            "stop_reason": stop_reason,
            "settings": asdict(settings),
        },
    )

    if settings.hybrid_save_path is not None:
        rmb.save(settings.hybrid_save_path)
        if settings.verbose:
            print(f"\nSaved hybrid boundary data to {settings.hybrid_save_path}")

    # Optional: keep your existing monotone bootstrap/confidence plot as a comparison.
    try:
        plot_monotone_fidelity_surface_with_confidence(
            data,
            settings,
            n_bootstrap=100,
            seed=settings.rng_seed,
        )
    except Exception as exc:  # keep final save robust even if plotting fails
        if settings.verbose:
            print(f"\nSkipping monotone confidence plot: {exc}")

    return rmb


if __name__ == "__main__":
    settings = HybridContourGPConfig(
        measurement_budget=10_000_000,
        hqc_budget=300.0,
        n_qubits_values=(20,),
        depth_bounds=(10, 80),
        ratio_bounds=(0.08, 0.8),
        random_elimination=0.3,
        scrambling_probability=0.85,
        min_adaptive_shots_per_config=1,
        max_adaptive_shots_per_config=5,
        max_shots_per_config=5,
        candidate_grid_size=(80, 80),
        boundary_width=0.08,
        shot_boundary_width=0.2,
        surface_smoothing=0.1,
        min_fit_points=8,
        monotone_l2=1e-3,
        contour_ready_probability_width=0.10,
        contour_min_anchors=4,
        contour_step_fraction=0.08,
        contour_projection_fraction=0.12,
        contour_gradient_fraction=0.01,
        ray_ratio_count=5,
        ray_probe_shots=2,
        ray_bisection_steps=5,
        ray_bisection_shots=2,
        trace_ratio_step_fraction=0.06,
        trace_depth_search_fraction=0.06,
        trace_correction_steps=2,
        trace_shots=2,
        trace_accept_probability_width=0.12,
        trace_directions=(-1, 1),
        refine_after_trace=False,  # replaced by GP uncertainty-aware refinement
        gp_refine_after_trace=True,
        gp_reserved_hqc_fraction=0.25,
        gp_refinement_shots=10,
        gp_max_refinement_executions=None,
        gp_candidate_grid_size=(80, 80),
        gp_local_depth_radius_fraction=0.10,
        gp_focus_on_measured_target_band=True,
        gp_focus_target_band_width=0.20,
        gp_focus_ratio_radius=0.06,
        gp_focus_depth_radius_fraction=0.12,
        force_ratio_coverage_before_gp=True,
        forced_ratio_coverage_count=6,
        forced_ratio_min_separation=0.04,
        contour_skip_existing_configs=True,
        contour_trace_min_ratio_step=0.01,
        max_contour_executions_per_ratio=None,
        contour_ratio_bin_width=0.01,
        plot_only_supported_ratios=True,
        plot_support_ratio_radius=0.08,
        gp_target_width_floor=0.03,
        gp_cost_power=1.0,
        hqc_cost_informed_acquisition=True,
        hqc_cost_power=1.0,
        hybrid_save_path=SCRIPT_DIR / "viarregio6_hybrid_boundary.json",
        hybrid_output_dir=SCRIPT_DIR / "viarregio6_hybrid_outputs",
        verbose=True,
    )

    rmb = estimate_boundary_hybrid(settings)
    print_fit_reports(rmb._data, settings)
