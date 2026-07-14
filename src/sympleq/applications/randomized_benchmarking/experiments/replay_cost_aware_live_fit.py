"""Replay saved COST_AWARE batches without submitting new backend jobs.

This script is intentionally load-only.  It reads the per-iteration
``batch_history`` written by ``cost_aware_surface_method.py``, rebuilds the RMB
data store one completed batch at a time, and runs the same live surface/volume
fit hooks used during a real run.  It never calls the acquisition code and never
calls ``backend.fidelity_estimation``.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    start_run,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_settings import (
    CostAwareSettings,
    control_panel_settings_kwargs,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_design import (
    _maybe_update_live_plots,
    _record_volume_point_from_posterior,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_plots import (
    _initial_rate_estimates,
    _initial_success_side_log_volume,
    _lindblad_reference_rate_estimates,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import (
    _I_L1,
    _I_L2,
    _I_M1,
    _I_M2,
    _I_N1,
    _I_N2,
    _I_V,
    _analytic_lindblad_gates,
    _config_from_checkpoint_record,
    _measured_arrays,
    _measured_config_count,
    _prior_moments,
    _stateless_posterior,
    _surface_plot_mesh,
    _wquantile,
)


DEFAULT_SEED = 2026
GENERATED_PERSONAL_ROOT = (
    Path("scripts")
    / "personal"
    / "randomized_benchmarking_personal"
    / "Personal"
)
DEFAULT_ROOT = GENERATED_PERSONAL_ROOT / "CostAware"

# --------------------------------------------------------------------------- #
# User controls
# --------------------------------------------------------------------------- #
# Edit these values and run this file directly. By default, command-line
# arguments are ignored so the replay is controlled from this block.
USE_COMMAND_LINE_ARGUMENTS = False

# ``run_numbers_to_load`` controls comparison replay. Use a single entry for
# the old one-run behaviour, or e.g. (1, 2) to compare runs side by side.
run_number_to_load = 2
run_numbers_to_load = (1, 2)

# Known run presets. Add future runs here once they exist. ``checkpoint=None``
# means "find the newest matching checkpoint under DEFAULT_ROOT/seed_*".
RUNS_TO_LOAD = {
    1: {
        "seed": 2026,
        "checkpoint": None,
        "allow_any_backend": False,
        "filename_contains": "",
        "backend_contains": "H2-",
        "match_number": 1,
    },
    2: {
        "seed": 20261,
        "checkpoint": None,
        "allow_any_backend": False,
        "filename_contains": "",
        "backend_contains": "H2-",
        "match_number": 1,
    },
}

# Generic search used when ``run_number_to_load`` is not in RUNS_TO_LOAD.
SEARCH_SEED = 2026
SEARCH_ALLOW_ANY_BACKEND = False
SEARCH_FILENAME_CONTAINS = ""
SEARCH_BACKEND_CONTAINS = "H2-"
SEARCH_MATCH_NUMBER = 1

# Replay/plot settings.
SHOW_PLOTS = True
PLOT_EVERY = 1
MAX_ITERATION = None
PAUSE_SECONDS = 0.1
NO_SURFACE_PLOT = False
NO_VOLUME_PLOT = False
# Tuple order:
#   (L1, L2, m1, m2, nu1, nu2, V)
# where L1/L2 are the base one-/two-qubit rates, m1/m2 are their linear
# Q-slopes, nu1/nu2 are their quadratic Q-curvatures, and V is visibility.
GRID_RESOLUTION_OVERRIDE = (15, 15, 9, 9, 1, 7, 1)
BOUNDARY_FIT_RESOLUTION_OVERRIDE = (17, 17, 9, 9, 1, 11, 1)
SURFACE_PNG = None
VOLUME_PNG = None


@dataclass
class ReplayResult:
    run_number: int
    label: str
    checkpoint: Path
    settings: CostAwareSettings
    rmb: object
    budget: Budget
    volume_history: list[tuple]
    surface_mesh: tuple[np.ndarray, np.ndarray, np.ndarray] | None
    total_shots: int
    replayed_batches: int


@dataclass
class ReplayState:
    run_number: int
    label: str
    checkpoint: Path
    settings: CostAwareSettings
    rmb: object
    budget: Budget
    batch_history: list[dict]
    volume_history: list[tuple]
    total_shots: int = 0
    replayed_batches: int = 0


def _parse_resolution(value: str | tuple[int, ...] | list[int] | None) -> tuple[int, ...] | None:
    if value is None:
        return None
    if isinstance(value, (tuple, list)):
        values = tuple(int(part) for part in value)
    else:
        values = tuple(int(part.strip()) for part in value.split(",") if part.strip())
    if not values:
        raise ValueError("resolution must contain at least one integer")
    return values


def _load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _checkpoint_backend(meta: dict) -> str | None:
    settings = meta.get("settings")
    if isinstance(settings, dict):
        value = settings.get("backend_model")
        if value is not None:
            return str(value)
    return None


def _checkpoint_seed(meta: dict) -> int | None:
    settings = meta.get("settings")
    if isinstance(settings, dict):
        value = settings.get("rng_seed")
        if value is not None:
            return int(value)
    return None


def _iter_checkpoint_candidates(
    seed: int,
    *,
    allow_any_backend: bool,
    root: Path = DEFAULT_ROOT,
) -> Iterable[Path]:
    seed_dir = root / f"seed_{seed}"
    if not seed_dir.exists():
        return []

    candidates = []
    for path in seed_dir.glob("*checkpoint.json"):
        try:
            meta = _load_json(path)
        except Exception:
            continue
        if _checkpoint_seed(meta) not in (None, seed):
            continue
        backend = _checkpoint_backend(meta)
        if not allow_any_backend and not (backend and backend.startswith("H2-")):
            continue
        if not isinstance(meta.get("batch_history"), list):
            continue
        candidates.append(path)
    return sorted(candidates, key=lambda path: path.stat().st_mtime, reverse=True)


def _latest_checkpoint(seed: int, *, allow_any_backend: bool) -> Path:
    candidates = list(
        _iter_checkpoint_candidates(seed, allow_any_backend=allow_any_backend)
    )
    if not candidates:
        suffix = "" if allow_any_backend else " H2"
        raise FileNotFoundError(
            f"No{suffix} COST_AWARE checkpoint found under "
            f"{DEFAULT_ROOT / f'seed_{seed}'}"
        )
    return candidates[0]


def _matching_checkpoints(
    *,
    seed: int,
    allow_any_backend: bool,
    filename_contains: str = "",
    backend_contains: str = "",
    root: Path = DEFAULT_ROOT,
) -> list[Path]:
    candidates = list(
        _iter_checkpoint_candidates(
            seed,
            allow_any_backend=allow_any_backend,
            root=root,
        )
    )
    if filename_contains:
        candidates = [path for path in candidates if filename_contains in path.name]
    if backend_contains:
        kept = []
        for path in candidates:
            try:
                backend = _checkpoint_backend(_load_json(path))
            except Exception:
                backend = None
            if backend is not None and backend_contains in backend:
                kept.append(path)
        candidates = kept
    return candidates


def _settings_from_args(args: argparse.Namespace, meta: dict) -> CostAwareSettings:
    kwargs = control_panel_settings_kwargs()
    kwargs["rng_seed"] = args.seed
    kwargs["plot"] = False
    kwargs["verbose"] = True
    kwargs["checkpoint_after_batch"] = False
    kwargs["resume_from_save"] = False
    kwargs["live_surface_plot"] = not args.no_surface_plot
    kwargs["live_volume_plot"] = not args.no_volume_plot
    kwargs["live_surface_plot_show"] = args.show
    kwargs["live_volume_plot_show"] = args.show
    kwargs["live_surface_plot_pause"] = args.pause
    kwargs["live_volume_plot_pause"] = args.pause
    kwargs["live_surface_plot_every"] = max(1, args.plot_every)
    kwargs["live_volume_plot_every"] = max(1, args.plot_every)

    saved_settings = meta.get("settings")
    if isinstance(saved_settings, dict):
        if saved_settings.get("backend_model") is not None:
            kwargs["backend_model"] = str(saved_settings["backend_model"])
        for key in (
            "initial_one_q_pauli_error",
            "initial_two_q_pauli_error",
            "initial_error_relative_uncertainty",
            "initial_one_q_error_relative_uncertainty",
            "initial_two_q_error_relative_uncertainty",
            "initial_visibility",
            "visibility_log_std",
            "visibility_bounds",
            "initial_rate_randomization_enabled",
            "initial_rate_random_seed",
            "initial_one_q_pauli_error_base",
            "initial_two_q_pauli_error_base",
            "initial_one_q_random_log_multiplier",
            "initial_two_q_random_log_multiplier",
            "initial_one_q_random_relative_std",
            "initial_two_q_random_relative_std",
            "initial_one_q_random_relative_delta",
            "initial_two_q_random_relative_delta",
        ):
            if key in saved_settings:
                kwargs[key] = saved_settings[key]

    grid_resolution = _parse_resolution(args.grid_resolution)
    boundary_fit_resolution = _parse_resolution(args.boundary_fit_resolution)
    if grid_resolution is not None:
        kwargs["grid_resolution"] = grid_resolution
    if boundary_fit_resolution is not None:
        kwargs["boundary_fit_resolution"] = boundary_fit_resolution

    settings = CostAwareSettings(**kwargs)
    if args.surface_png is not None:
        settings = replace(settings, live_surface_plot_path=Path(args.surface_png))
    if args.volume_png is not None:
        settings = replace(settings, live_volume_plot_path=Path(args.volume_png))
    return settings


def _add_request_measurements(rmb, request_record: dict) -> int:
    config_record = request_record.get("config")
    if not isinstance(config_record, dict):
        return 0
    config = _config_from_checkpoint_record(config_record)
    estimator = rmb._data.get(config)
    if estimator is None:
        estimator = rmb.backend.default_estimator()
        rmb._data[config] = estimator

    n_added = 0
    for measurement in request_record.get("measurements", []):
        if not isinstance(measurement, dict) or "outcome" not in measurement:
            continue
        estimator.record(bool(measurement["outcome"]))
        n_added += 1
    return n_added


def _sorted_batch_history(meta: dict, checkpoint_path: Path) -> list[dict]:
    batch_history = list(meta.get("batch_history", []))
    if not batch_history:
        raise RuntimeError(f"No batch_history in {checkpoint_path}")
    return sorted(batch_history, key=lambda row: int(row.get("iteration", 0)))


def _run_label(run_number: int, checkpoint_path: Path, meta: dict) -> str:
    seed = _checkpoint_seed(meta)
    backend = _checkpoint_backend(meta)
    parts = [f"run {run_number}"]
    if seed is not None:
        parts.append(f"seed {seed}")
    if backend:
        parts.append(backend)
    return " / ".join(parts)


def _make_replay_state(
    checkpoint_path: Path,
    settings: CostAwareSettings,
    *,
    run_number: int,
) -> ReplayState:
    meta = _load_json(checkpoint_path)
    _, rmb, budget = start_run(settings)
    state = ReplayState(
        run_number=run_number,
        label=_run_label(run_number, checkpoint_path, meta),
        checkpoint=checkpoint_path,
        settings=settings,
        rmb=rmb,
        budget=budget,
        batch_history=_sorted_batch_history(meta, checkpoint_path),
        volume_history=[],
    )
    print(f"Replaying {checkpoint_path}")
    print(f"  label: {state.label}")
    print(f"  backend in checkpoint: {_checkpoint_backend(meta)}")
    print(f"  grid_resolution: {settings.grid_resolution}")
    print(f"  boundary_fit_resolution: {settings.boundary_fit_resolution}")
    print("  no backend calls or acquisition calls will be made")
    return state


def _apply_replay_batch(state: ReplayState, batch: dict) -> int:
    added = 0
    for request_record in batch.get("requests", []):
        added += _add_request_measurements(state.rmb, request_record)
    state.total_shots += added
    state.replayed_batches += 1

    cost = float(batch.get("cost_hqc", 0.0))
    requested_shots = int(batch.get("requested_shots", added))
    state.budget.spend_batch(cost, requested_shots)

    iteration = int(batch.get("iteration", state.replayed_batches))
    params, weights = _stateless_posterior(state.rmb._data, state.settings)
    if params is not None and weights is not None:
        secondary_reference_gates = (
            _analytic_lindblad_gates
            if state.settings.gp_grid_surface_path is not None
            else None
        )
        _record_volume_point_from_posterior(
            params,
            weights,
            state.settings,
            iteration,
            state.volume_history,
            secondary_reference_gates=secondary_reference_gates,
        )

    print(
        f"{state.label} iteration {iteration}: added {added} shot outcomes, "
        f"configs={_measured_config_count(state.rmb._data)}, "
        f"spent={state.budget.spent_hqc:.3f} HQC"
    )
    return added


def _state_result(state: ReplayState) -> ReplayResult:
    params, weights = _stateless_posterior(state.rmb._data, state.settings)
    surface_mesh = (
        None
        if params is None or weights is None
        else _surface_plot_mesh(params, weights, state.settings)
    )
    return ReplayResult(
        run_number=state.run_number,
        label=state.label,
        checkpoint=state.checkpoint,
        settings=state.settings,
        rmb=state.rmb,
        budget=state.budget,
        volume_history=list(state.volume_history),
        surface_mesh=surface_mesh,
        total_shots=state.total_shots,
        replayed_batches=state.replayed_batches,
    )


def _weighted_mean_std(values: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(keep):
        return float("nan"), float("nan")
    v = values[keep]
    w = weights[keep]
    w = w / np.sum(w)
    mean = float(np.sum(w * v))
    var = float(np.sum(w * (v - mean) ** 2))
    return mean, float(np.sqrt(max(var, 0.0)))


def _print_fit_parameter_summary(result: ReplayResult) -> None:
    params, weights = _stateless_posterior(result.rmb._data, result.settings)
    if params is None or weights is None:
        print(f"\nFit parameter summary for {result.label}: no posterior fit available")
        return

    names = ("L1", "L2", "m1", "m2", "nu1", "nu2", "V")
    indices = (_I_L1, _I_L2, _I_M1, _I_M2, _I_N1, _I_N2, _I_V)
    print(f"\nFit parameter summary for {result.label}")
    print("  posterior weighted natural parameters:")
    for name, index in zip(names, indices):
        values = params[:, index]
        mean, std = _weighted_mean_std(values, weights)
        median = _wquantile(values, weights, 0.5)
        lo = float(np.nanmin(values))
        hi = float(np.nanmax(values))
        print(
            f"    {name:>3}: mean={mean:.6g}, median={median:.6g}, "
            f"std={std:.3g}, grid=[{lo:.6g}, {hi:.6g}]"
        )

    _, prior_std = _prior_moments(result.settings)
    q_values = tuple(float(q) for q in result.settings.q_values)
    q_span = max(max(q_values) - min(q_values), 1.0)
    q_half = max(0.5 * q_span, 1.0)
    nu2 = params[:, _I_N2]
    nu2_mean, nu2_std = _weighted_mean_std(nu2, weights)
    nu2_median = _wquantile(nu2, weights, 0.5)
    nu2_lo = float(np.nanmin(nu2))
    nu2_hi = float(np.nanmax(nu2))
    edge_mean_delta = nu2_mean * q_half * q_half
    edge_median_delta = nu2_median * q_half * q_half
    edge_grid_limit = max(abs(nu2_lo), abs(nu2_hi)) * q_half * q_half
    l2_mean, _ = _weighted_mean_std(params[:, _I_L2], weights)
    print("  two-qubit quadratic term:")
    print(
        f"    nu2 mean={nu2_mean:.6g}, median={nu2_median:.6g}, "
        f"std={nu2_std:.3g}"
    )
    print(
        f"    allowed nu2 grid support=[{nu2_lo:.6g}, {nu2_hi:.6g}], "
        f"max |nu2|={max(abs(nu2_lo), abs(nu2_hi)):.6g}"
    )
    print(
        f"    prior sigma for nu2={prior_std[_I_N2]:.6g}; "
        f"boundary_fit_resolution nu2 axis={result.settings.boundary_fit_resolution[_I_N2]}"
    )
    print(
        f"    implied L2 edge change over half Q-span {q_half:.3g}: "
        f"mean={edge_mean_delta:.6g}, median={edge_median_delta:.6g}, "
        f"grid max |delta|={edge_grid_limit:.6g}"
    )
    if np.isfinite(l2_mean) and l2_mean != 0.0:
        print(
            f"    relative to fitted L2 mean: edge mean={edge_mean_delta / l2_mean:.3%}, "
            f"grid max={edge_grid_limit / abs(l2_mean):.3%}"
        )


def replay_collect_fit(
    checkpoint_path: Path,
    settings: CostAwareSettings,
    *,
    run_number: int,
    max_iteration: int | None = None,
) -> ReplayResult:
    """Replay a run into memory and return final surface/volume summaries."""
    meta = _load_json(checkpoint_path)
    batch_history = _sorted_batch_history(meta, checkpoint_path)

    _, rmb, budget = start_run(settings)
    volume_history: list[tuple] = []
    total_shots = 0
    replayed_batches = 0

    print(f"Replaying {checkpoint_path}")
    print(f"  backend in checkpoint: {_checkpoint_backend(meta)}")
    print(f"  grid_resolution: {settings.grid_resolution}")
    print(f"  boundary_fit_resolution: {settings.boundary_fit_resolution}")
    print("  no backend calls or acquisition calls will be made")

    for batch in batch_history:
        iteration = int(batch.get("iteration", replayed_batches + 1))
        if max_iteration is not None and iteration > max_iteration:
            break

        added = 0
        for request_record in batch.get("requests", []):
            added += _add_request_measurements(rmb, request_record)
        total_shots += added
        replayed_batches += 1

        cost = float(batch.get("cost_hqc", 0.0))
        requested_shots = int(batch.get("requested_shots", added))
        budget.spend_batch(cost, requested_shots)

        params, weights = _stateless_posterior(rmb._data, settings)
        if params is not None and weights is not None:
            secondary_reference_gates = (
                _analytic_lindblad_gates
                if settings.gp_grid_surface_path is not None
                else None
            )
            _record_volume_point_from_posterior(
                params,
                weights,
                settings,
                iteration,
                volume_history,
                secondary_reference_gates=secondary_reference_gates,
            )

        print(
            f"iteration {iteration}: added {added} shot outcomes, "
            f"configs={_measured_config_count(rmb._data)}, "
            f"spent={budget.spent_hqc:.3f} HQC"
        )

    params, weights = _stateless_posterior(rmb._data, settings)
    surface_mesh = (
        None
        if params is None or weights is None
        else _surface_plot_mesh(params, weights, settings)
    )
    label = _run_label(run_number, checkpoint_path, meta)
    print(
        f"Done {label}: replayed {replayed_batches} batch(es), "
        f"{total_shots} shot outcome(s), "
        f"{_measured_config_count(rmb._data)} distinct config(s)."
    )
    return ReplayResult(
        run_number=run_number,
        label=label,
        checkpoint=checkpoint_path,
        settings=settings,
        rmb=rmb,
        budget=budget,
        volume_history=volume_history,
        surface_mesh=surface_mesh,
        total_shots=total_shots,
        replayed_batches=replayed_batches,
    )


def replay_live_fit(
    checkpoint_path: Path,
    settings: CostAwareSettings,
    *,
    max_iteration: int | None = None,
) -> None:
    meta = _load_json(checkpoint_path)
    batch_history = list(meta.get("batch_history", []))
    if not batch_history:
        raise RuntimeError(f"No batch_history in {checkpoint_path}")

    _, rmb, budget = start_run(settings)
    volume_history = []
    total_shots = 0
    replayed_batches = 0

    print(f"Replaying {checkpoint_path}")
    print(f"  backend in checkpoint: {_checkpoint_backend(meta)}")
    print(f"  grid_resolution: {settings.grid_resolution}")
    print(f"  boundary_fit_resolution: {settings.boundary_fit_resolution}")
    print("  no backend calls or acquisition calls will be made")

    for batch in sorted(batch_history, key=lambda row: int(row.get("iteration", 0))):
        iteration = int(batch.get("iteration", replayed_batches + 1))
        if max_iteration is not None and iteration > max_iteration:
            break

        added = 0
        for request_record in batch.get("requests", []):
            added += _add_request_measurements(rmb, request_record)
        total_shots += added
        replayed_batches += 1

        cost = float(batch.get("cost_hqc", 0.0))
        requested_shots = int(batch.get("requested_shots", added))
        budget.spend_batch(cost, requested_shots)

        print(
            f"iteration {iteration}: added {added} shot outcomes, "
            f"configs={_measured_config_count(rmb._data)}, "
            f"spent={budget.spent_hqc:.3f} HQC"
        )
        _maybe_update_live_plots(
            rmb,
            settings,
            budget,
            iteration,
            volume_history,
        )

    print(
        f"Done: replayed {replayed_batches} batch(es), "
        f"{total_shots} shot outcome(s), "
        f"{_measured_config_count(rmb._data)} distinct config(s)."
    )
    result = ReplayResult(
        run_number=run_number_to_load,
        label=_run_label(run_number_to_load, checkpoint_path, meta),
        checkpoint=checkpoint_path,
        settings=settings,
        rmb=rmb,
        budget=budget,
        volume_history=volume_history,
        surface_mesh=None,
        total_shots=total_shots,
        replayed_batches=replayed_batches,
    )
    _print_fit_parameter_summary(result)
    if settings.live_surface_plot_path is not None and settings.live_surface_plot:
        print(f"  surface plot: {settings.live_surface_plot_path}")
    if settings.live_volume_plot_path is not None and settings.live_volume_plot:
        print(f"  volume plot: {settings.live_volume_plot_path}")
    if settings.live_surface_plot_show or settings.live_volume_plot_show:
        plt.show()


def _comparison_output_path(path: str | Path | None, default_name: str) -> Path:
    if path is None:
        return DEFAULT_ROOT / "comparison" / default_name
    p = Path(path)
    return p.with_name(f"{p.stem}_comparison{p.suffix}")


def _plot_surface_comparison(
    results: list[ReplayResult],
    *,
    png_path: str | Path | None,
    show: bool,
    pause: float,
) -> Path | None:
    results = [result for result in results if result.surface_mesh is not None]
    if not results:
        print("surface comparison skipped: no finite surface meshes")
        return None

    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as local_plt
    from matplotlib.colors import TwoSlopeNorm
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = local_plt.figure(
        num="cost-aware replay surface comparison",
        figsize=(7.0 * len(results), 6.4),
        clear=True,
    )
    axes = []
    scatter_for_colorbar = None
    for index, result in enumerate(results, start=1):
        rr, qq, log10_gates = result.surface_mesh
        settings = result.settings
        ax = fig.add_subplot(1, len(results), index, projection="3d")
        axes.append(ax)
        ax.plot_surface(
            rr,
            qq,
            log10_gates,
            cmap="viridis",
            alpha=0.52,
            linewidth=0,
            antialiased=True,
        )

        n_bounds = (
            settings.n_gates_bounds
            if settings.surface_plot_n_gates_bounds is None
            else settings.surface_plot_n_gates_bounds
        )
        n_lo, n_hi = n_bounds
        ratios_m, gates_m, qubits_m, succ_m, fail_m = _measured_arrays(result.rmb._data)
        total = succ_m + fail_m
        mask = total > 0
        if np.any(mask):
            p_hat = succ_m[mask] / total[mask]
            size = 16.0 + 35.0 * np.clip(total[mask] / max(total[mask].max(), 1.0), 0, 1)
            scatter = ax.scatter(
                ratios_m[mask],
                qubits_m[mask],
                np.log10(np.maximum(gates_m[mask], 1e-12)),
                c=p_hat,
                cmap="coolwarm_r",
                norm=TwoSlopeNorm(vcenter=0.5, vmin=0.0, vmax=1.0),
                s=size,
                edgecolor="k",
                linewidth=0.25,
                depthshade=True,
            )
            scatter_for_colorbar = scatter

        ax.set_xlabel("two-qubit ratio r")
        ax.set_ylabel("n_qubits Q")
        ax.set_zlabel(r"$\log_{10}(n_{\rm gates})$")
        ax.set_zlim(np.log10(max(n_lo, 1e-12)), np.log10(max(n_hi, 1e-12)))
        ax.set_title(result.label)

    fig.subplots_adjust(left=0.04, right=0.90, bottom=0.08, top=0.90, wspace=0.05)
    if scatter_for_colorbar is not None:
        cbar_axis = fig.add_axes([0.925, 0.24, 0.014, 0.52])
        cbar = fig.colorbar(scatter_for_colorbar, cax=cbar_axis)
        cbar.set_label("measured survival fraction")

    out = _comparison_output_path(png_path, "boundary_surface_3d_comparison.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    if show:
        local_plt.show(block=False)
        fig.canvas.draw()
        fig.canvas.flush_events()
        local_plt.pause(max(float(pause), 0.001))
    else:
        local_plt.close(fig)
    print(f"surface comparison plot -> {out}")
    return out


def _history_arrays(history: list[tuple]) -> dict[str, np.ndarray]:
    return {
        "iteration": np.asarray([row[0] for row in history], dtype=int),
        "fit": np.asarray([row[1] for row in history], dtype=float),
        "reference": np.asarray([row[2] for row in history], dtype=float),
        "lower": np.asarray([row[3] if len(row) > 3 else np.nan for row in history], dtype=float),
        "upper": np.asarray([row[4] if len(row) > 4 else np.nan for row in history], dtype=float),
        "secondary_reference": np.asarray([row[5] if len(row) > 5 else np.nan for row in history], dtype=float),
        "l1": np.asarray([row[6] if len(row) > 6 else np.nan for row in history], dtype=float),
        "l2": np.asarray([row[7] if len(row) > 7 else np.nan for row in history], dtype=float),
        "l1_std": np.asarray([row[8] if len(row) > 8 else np.nan for row in history], dtype=float),
        "l2_std": np.asarray([row[9] if len(row) > 9 else np.nan for row in history], dtype=float),
    }


def _plot_volume_comparison(
    results: list[ReplayResult],
    *,
    png_path: str | Path | None,
    show: bool,
    pause: float,
) -> Path | None:
    results = [result for result in results if result.volume_history]
    if not results:
        print("volume comparison skipped: no finite volume histories")
        return None

    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as local_plt

    fig, axes = local_plt.subplots(
        3,
        1,
        num="cost-aware replay volume comparison",
        figsize=(10.0, 9.0),
        clear=True,
        sharex=True,
        gridspec_kw={"height_ratios": [1.45, 1.0, 1.0]},
    )
    ax_vol, ax_l1, ax_l2 = axes
    colors = local_plt.rcParams["axes.prop_cycle"].by_key().get("color", [])
    linestyles = ["-", "--", "-.", ":"]

    for index, result in enumerate(results):
        arrays = _history_arrays(result.volume_history)
        color = colors[index % len(colors)] if colors else None
        linestyle = linestyles[index % len(linestyles)]
        iterations = arrays["iteration"]
        fitted = arrays["fit"]
        lower = arrays["lower"]
        upper = arrays["upper"]
        has_band = np.isfinite(lower) & np.isfinite(upper)
        if np.any(has_band):
            yerr = np.vstack((
                np.where(has_band, np.maximum(fitted - lower, 0.0), np.nan),
                np.where(has_band, np.maximum(upper - fitted, 0.0), np.nan),
            ))
            ax_vol.errorbar(
                iterations,
                fitted,
                yerr=yerr,
                marker="o",
                linewidth=1.7,
                elinewidth=0.9,
                capsize=2,
                linestyle=linestyle,
                color=color,
                label=result.label,
            )
        else:
            ax_vol.plot(
                iterations,
                fitted,
                marker="o",
                linewidth=1.7,
                linestyle=linestyle,
                color=color,
                label=result.label,
            )

        prior_l1, prior_l2 = _initial_rate_estimates(result.settings)
        prior_volume = _initial_success_side_log_volume(result.settings, prior_l1, prior_l2)
        if np.isfinite(prior_volume):
            ax_vol.scatter([0], [prior_volume], marker="x", s=70, color=color, zorder=5)

        for axis, key, std_key, prior in (
            (ax_l1, "l1", "l1_std", prior_l1),
            (ax_l2, "l2", "l2_std", prior_l2),
        ):
            values = arrays[key]
            std = arrays[std_key]
            if np.any(np.isfinite(values)):
                axis.errorbar(
                    iterations,
                    values,
                    yerr=np.where(np.isfinite(std), std, np.nan),
                    marker="o",
                    linewidth=1.5,
                    elinewidth=0.8,
                    capsize=2,
                    linestyle=linestyle,
                    color=color,
                    label=result.label,
                )
            axis.scatter([0], [prior], marker="x", s=60, color=color, zorder=5)

    first = results[0]
    first_arrays = _history_arrays(first.volume_history)
    reference = first_arrays["reference"]
    if np.any(np.isfinite(reference)):
        ax_vol.axhline(
            float(reference[np.isfinite(reference)][-1]),
            color="black",
            linestyle="--",
            linewidth=1.2,
            alpha=0.7,
            label="primary reference",
        )
    secondary = first_arrays["secondary_reference"]
    if np.any(np.isfinite(secondary)):
        ax_vol.axhline(
            float(secondary[np.isfinite(secondary)][-1]),
            color="darkorange",
            linestyle="--",
            linewidth=1.2,
            alpha=0.75,
            label="analytic Lindblad reference",
        )
    lindblad_l1, lindblad_l2 = _lindblad_reference_rate_estimates(first.settings)
    ax_l1.axhline(lindblad_l1, color="black", linestyle="--", linewidth=1.2, alpha=0.7)
    ax_l2.axhline(lindblad_l2, color="black", linestyle="--", linewidth=1.2, alpha=0.7)

    ax_vol.set_ylabel(r"success-side volume in $(\log_{10} n, r, Q)$")
    ax_vol.set_title("COST_AWARE replay comparison")
    ax_l1.set_ylabel(r"$L_1(Q_{\rm ref})$")
    ax_l2.set_ylabel(r"$L_2(Q_{\rm ref})$")
    ax_l2.set_xlabel("iteration")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(loc="best")

    fig.tight_layout()
    out = _comparison_output_path(png_path, "boundary_volume_comparison.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    if show:
        local_plt.show(block=False)
        fig.canvas.draw()
        fig.canvas.flush_events()
        local_plt.pause(max(float(pause), 0.001))
    else:
        local_plt.close(fig)
    print(f"volume comparison plot -> {out}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Replay saved COST_AWARE iteration batches without submitting jobs."
    )
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=None,
        help="Explicit *_checkpoint.json to replay. Defaults to latest H2 checkpoint for seed.",
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--allow-any-backend",
        action="store_true",
        help="Allow latest checkpoint from any backend instead of only H2-*.",
    )
    parser.add_argument(
        "--grid-resolution",
        default=None,
        help="Comma-separated grid_resolution override, e.g. 17,17,9,9,1,9,1.",
    )
    parser.add_argument(
        "--boundary-fit-resolution",
        default=None,
        help="Comma-separated boundary_fit_resolution override.",
    )
    parser.add_argument("--max-iteration", type=int, default=None)
    parser.add_argument("--plot-every", type=int, default=1)
    parser.add_argument("--pause", type=float, default=0.1)
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--no-surface-plot", action="store_true")
    parser.add_argument("--no-volume-plot", action="store_true")
    parser.add_argument("--surface-png", type=Path, default=None)
    parser.add_argument("--volume-png", type=Path, default=None)
    args = parser.parse_args()

    checkpoint = args.checkpoint
    if checkpoint is None:
        checkpoint = _latest_checkpoint(
            args.seed,
            allow_any_backend=args.allow_any_backend,
        )
    meta = _load_json(checkpoint)
    settings = _settings_from_args(args, meta)
    replay_live_fit(checkpoint, settings, max_iteration=args.max_iteration)


def _control_panel_args_for_run(selected_run_number: int) -> argparse.Namespace:
    selected = RUNS_TO_LOAD.get(
        selected_run_number,
        {
            "seed": SEARCH_SEED,
            "checkpoint": None,
            "allow_any_backend": SEARCH_ALLOW_ANY_BACKEND,
            "filename_contains": SEARCH_FILENAME_CONTAINS,
            "backend_contains": SEARCH_BACKEND_CONTAINS,
            "match_number": SEARCH_MATCH_NUMBER,
        },
    )
    checkpoint = selected.get("checkpoint")
    if checkpoint is not None:
        checkpoint = Path(checkpoint)
    else:
        matches = _matching_checkpoints(
            seed=int(selected["seed"]),
            allow_any_backend=bool(selected.get("allow_any_backend", False)),
            filename_contains=str(selected.get("filename_contains", "")),
            backend_contains=str(selected.get("backend_contains", "")),
        )
        match_number = max(1, int(selected.get("match_number", 1)))
        if len(matches) < match_number:
            raise FileNotFoundError(
                f"Requested match {match_number}, but found {len(matches)} "
                f"checkpoint(s) for run {selected_run_number}."
            )
        checkpoint = matches[match_number - 1]

    return argparse.Namespace(
        checkpoint=checkpoint,
        seed=int(selected["seed"]),
        allow_any_backend=bool(selected.get("allow_any_backend", False)),
        grid_resolution=GRID_RESOLUTION_OVERRIDE,
        boundary_fit_resolution=BOUNDARY_FIT_RESOLUTION_OVERRIDE,
        max_iteration=MAX_ITERATION,
        plot_every=PLOT_EVERY,
        pause=PAUSE_SECONDS,
        show=SHOW_PLOTS,
        no_surface_plot=NO_SURFACE_PLOT,
        no_volume_plot=NO_VOLUME_PLOT,
        surface_png=None if SURFACE_PNG is None else Path(SURFACE_PNG),
        volume_png=None if VOLUME_PNG is None else Path(VOLUME_PNG),
    )


def _control_panel_args() -> argparse.Namespace:
    return _control_panel_args_for_run(run_number_to_load)


def main_comparison_from_control_panel() -> None:
    selected_runs = tuple(int(run_number) for run_number in run_numbers_to_load)
    if len(selected_runs) < 2:
        raise ValueError("comparison replay needs at least two run numbers")

    print(f"[replay control] comparing runs {selected_runs}")
    states: list[ReplayState] = []
    for selected_run_number in selected_runs:
        args = _control_panel_args_for_run(selected_run_number)
        print(
            f"[replay control] run_number={selected_run_number} "
            f"checkpoint={args.checkpoint}"
        )
        meta = _load_json(args.checkpoint)
        settings = _settings_from_args(args, meta)
        state = _make_replay_state(
            args.checkpoint,
            settings,
            run_number=selected_run_number,
        )
        if args.max_iteration is not None:
            state.batch_history = [
                batch
                for batch in state.batch_history
                if int(batch.get("iteration", 0)) <= int(args.max_iteration)
            ]
        states.append(state)

    iterations = sorted(
        {
            int(batch.get("iteration", 0))
            for state in states
            for batch in state.batch_history
        }
    )
    batch_by_iteration = [
        {
            int(batch.get("iteration", 0)): batch
            for batch in state.batch_history
        }
        for state in states
    ]

    for step_index, iteration in enumerate(iterations, start=1):
        print(f"[replay control] lockstep iteration {iteration}")
        for state, batches in zip(states, batch_by_iteration):
            batch = batches.get(iteration)
            if batch is None:
                print(f"{state.label} iteration {iteration}: no saved batch")
                continue
            _apply_replay_batch(state, batch)

        should_plot = (
            step_index == 1
            or step_index == len(iterations)
            or iteration % max(1, PLOT_EVERY) == 0
        )
        if not should_plot:
            continue

        results = [_state_result(state) for state in states]
        if not NO_SURFACE_PLOT:
            _plot_surface_comparison(
                results,
                png_path=SURFACE_PNG,
                show=SHOW_PLOTS,
                pause=PAUSE_SECONDS,
            )
        if not NO_VOLUME_PLOT:
            _plot_volume_comparison(
                results,
                png_path=VOLUME_PNG,
                show=SHOW_PLOTS,
                pause=PAUSE_SECONDS,
            )

    results = [_state_result(state) for state in states]
    for result in results:
        print(
            f"Done {result.label}: replayed {result.replayed_batches} batch(es), "
            f"{result.total_shots} shot outcome(s), "
            f"{_measured_config_count(result.rmb._data)} distinct config(s)."
        )
        _print_fit_parameter_summary(result)

    if SHOW_PLOTS and not (NO_SURFACE_PLOT and NO_VOLUME_PLOT):
        plt.show()


def main_from_control_panel() -> None:
    if len(tuple(run_numbers_to_load)) > 1:
        main_comparison_from_control_panel()
        return

    args = _control_panel_args()
    print(
        f"[replay control] run_number_to_load={run_number_to_load} "
        f"checkpoint={args.checkpoint}"
    )
    meta = _load_json(args.checkpoint)
    settings = _settings_from_args(args, meta)
    replay_live_fit(args.checkpoint, settings, max_iteration=args.max_iteration)


if __name__ == "__main__":
    if USE_COMMAND_LINE_ARGUMENTS and len(sys.argv) > 1:
        main()
    else:
        main_from_control_panel()
