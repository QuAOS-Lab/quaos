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
from dataclasses import replace
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt

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
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import (
    _config_from_checkpoint_record,
    _measured_config_count,
)


DEFAULT_SEED = 2026
DEFAULT_ROOT = Path("Personal") / "CostAware"

# --------------------------------------------------------------------------- #
# User controls
# --------------------------------------------------------------------------- #
# Edit these values and run this file directly. By default, command-line
# arguments are ignored so the replay is controlled from this block.
USE_COMMAND_LINE_ARGUMENTS = False

# 1 means "the current/default H2 seed-2026 run".
run_number_to_load = 1

# Known run presets. Add future runs here once they exist. ``checkpoint=None``
# means "find the newest matching checkpoint under Personal/CostAware/seed_*".
RUNS_TO_LOAD = {
    1: {
        "seed": 2026,
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
GRID_RESOLUTION_OVERRIDE = (15, 15, 9, 9, 1, 1, 1)
BOUNDARY_FIT_RESOLUTION_OVERRIDE = (17, 17, 9, 9, 1, 1, 1)
SURFACE_PNG = None
VOLUME_PNG = None


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
    if settings.live_surface_plot_path is not None and settings.live_surface_plot:
        print(f"  surface plot: {settings.live_surface_plot_path}")
    if settings.live_volume_plot_path is not None and settings.live_volume_plot:
        print(f"  volume plot: {settings.live_volume_plot_path}")
    if settings.live_surface_plot_show or settings.live_volume_plot_show:
        plt.show()


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


def _control_panel_args() -> argparse.Namespace:
    selected = RUNS_TO_LOAD.get(
        run_number_to_load,
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
                f"checkpoint(s) for run_number_to_load={run_number_to_load}."
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


def main_from_control_panel() -> None:
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
