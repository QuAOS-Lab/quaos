"""Runnable cost-aware surface-design script.

The core method implementation lives in ``cost_aware_surface_method.py``.  This
file keeps the historical import path working and owns the live/final plotting
around that method.
"""
from __future__ import annotations

import json
import numpy as np

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.experiments import (
    cost_aware_surface_method as _method,
)
from sympleq.applications.randomized_benchmarking.experiments.common import (
    print_progress,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_plots import (
    plot_boundary_surface,
    plot_live_volume_history,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import *  # noqa: F401,F403
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import (
    EXPERIMENTS_DIR,
    CostAwareSurfaceSettings,
    _LN2,
    _analytic_lindblad_gates,
    _asymptote,
    _measured_arrays,
    _measured_config_count,
    _stateless_posterior,
    _surface_s1_s2,
    _surface_plot_mesh,
    run_with_budget,
)

LiveVolumePoint = tuple[int, float, float, float, float]


def __getattr__(name: str):
    """Delegate private legacy imports to the method module."""
    return getattr(_method, name)


def _volume_reference_gates(settings: CostAwareSurfaceSettings):
    if settings.gp_grid_surface_path is not None:
        return _method._gp_grid_surface_gates
    return _analytic_lindblad_gates


def _record_live_volume_point(
    history: list[LiveVolumePoint],
    iteration: int,
    fitted_volume: float,
    true_volume: float,
    lower_volume: float,
    upper_volume: float,
) -> None:
    point = (
        int(iteration),
        float(fitted_volume),
        float(true_volume),
        float(lower_volume),
        float(upper_volume),
    )
    for i, existing in enumerate(history):
        if int(existing[0]) == int(iteration):
            history[i] = point
            return
    history.append(point)
    history.sort(key=lambda p: p[0])


def _record_volume_point_from_posterior(
    params,
    weights,
    settings: CostAwareSurfaceSettings,
    iteration: int,
    history: list[LiveVolumePoint],
) -> None:
    scores = _surface_s1_s2(params, weights, settings, _volume_reference_gates(settings))
    fitted_volume = float(scores["surface_volume_fit"])
    true_volume = float(scores["surface_volume_reference"])
    if not (np.isfinite(fitted_volume) and np.isfinite(true_volume)):
        return
    _record_live_volume_point(
        history,
        iteration,
        fitted_volume,
        true_volume,
        float(scores["surface_volume_lower_1sigma"]),
        float(scores["surface_volume_upper_1sigma"]),
    )


def _restore_volume_history_from_meta(settings: CostAwareSurfaceSettings) -> list[LiveVolumePoint]:
    if settings.save_path is None:
        return []
    meta_path = _method._checkpoint_meta_path(settings.save_path)
    if not meta_path.exists():
        return []
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception:
        return []
    if not _method._checkpoint_q_values_match(meta, settings):
        return []

    history = []
    for row in meta.get("live_volume_history", []):
        try:
            history.append(tuple(row))
        except TypeError:
            pass
    if history:
        return [
            (int(i), float(fit), float(ref), float(lo), float(hi))
            for i, fit, ref, lo, hi in history
        ]
    return _rebuild_volume_history_from_batches(meta, settings)


def _rebuild_volume_history_from_batches(
    meta: dict,
    settings: CostAwareSurfaceSettings,
) -> list[LiveVolumePoint]:
    if settings.save_path is None:
        return []
    try:
        template = RMB.load(settings.save_path)
    except Exception:
        return []

    history: list[LiveVolumePoint] = []
    data = {}
    for batch in meta.get("batch_history", []):
        iteration = batch.get("iteration")
        if iteration is None:
            continue
        for request in batch.get("requests", []):
            config_record = request.get("config")
            if not isinstance(config_record, dict):
                continue
            config = _method._config_from_checkpoint_record(config_record)
            estimator = data.setdefault(config, template.backend.default_estimator())
            for measurement in request.get("measurements", []):
                if isinstance(measurement, dict) and "outcome" in measurement:
                    estimator.record(bool(measurement["outcome"]))
        params, weights = _stateless_posterior(data, settings)
        if params is not None and weights is not None:
            _record_volume_point_from_posterior(
                params, weights, settings, int(iteration), history
            )
    return history


def _save_volume_history_to_meta(
    settings: CostAwareSurfaceSettings,
    volume_history: list[LiveVolumePoint],
) -> None:
    if settings.save_path is None or not settings.checkpoint_after_batch:
        return
    meta_path = _method._checkpoint_meta_path(settings.save_path)
    if not meta_path.exists():
        return
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if not _method._checkpoint_q_values_match(meta, settings):
            return
        meta["live_volume_history"] = [list(point) for point in volume_history]
        _method._write_json_atomic(meta_path, meta)
    except Exception:
        return


def _live_plot_params(data, settings: CostAwareSurfaceSettings, pending_requests=None):
    params, weights = _stateless_posterior(data, settings)
    if params is not None and weights is not None:
        return params, weights
    if pending_requests is None:
        return None, None

    centre, std = _method._prior_moments(settings)
    params, log_prior, _ = _method._build_grid(
        centre,
        settings.grid_halfwidth_sigmas * std,
        settings,
        resolution=settings.boundary_fit_resolution,
    )
    keep = np.isfinite(log_prior)
    if not np.any(keep):
        return None, None
    log_prior = np.where(keep, log_prior, -np.inf)
    log_prior -= np.max(log_prior[keep])
    weights = np.where(keep, np.exp(log_prior), 0.0)
    weights /= np.sum(weights)
    return params, weights


def _pending_arrays(pending_requests, settings: CostAwareSurfaceSettings):
    if not pending_requests:
        return None
    pending_requests = _method._restrict_requests_to_qubit_window(
        list(pending_requests), settings
    )
    configs = [request.config for request in pending_requests]
    return (
        np.asarray([float(config.ratio_2_qb_gates) for config in configs]),
        np.asarray([float(config.n_gates) for config in configs]),
        np.asarray([float(config.n_qubits) for config in configs]),
    )


def _plot_boundary_surface_with_pending(
    mesh,
    measured,
    settings: CostAwareSurfaceSettings,
    pending_requests=None,
) -> None:
    figure_name = "cost-aware live surface"
    plot_boundary_surface(
        *mesh,
        measured,
        settings,
        png_path=settings.live_surface_plot_path,
        show=settings.live_surface_plot_show,
        show_block=False,
        show_pause=settings.live_surface_plot_pause,
        close=False,
        figure_name=figure_name,
    )
    pending = _pending_arrays(pending_requests, settings)
    import matplotlib.pyplot as plt

    if pending is not None:
        ratios, gates, qubits = pending
        fig = plt.figure(num=figure_name)
        ax = fig.axes[0]
        ax.scatter(
            ratios,
            qubits,
            np.log10(np.maximum(gates, 1e-12)),
            c="0.55",
            s=46,
            edgecolor="k",
            linewidth=0.5,
            depthshade=True,
            label="submitted, awaiting result",
        )
        ax.legend(loc="upper left")
        if settings.live_surface_plot_path is not None:
            fig.savefig(settings.live_surface_plot_path, dpi=150)
        if settings.live_surface_plot_show:
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(max(float(settings.live_surface_plot_pause), 0.001))
        else:
            plt.close(fig)
    elif settings.live_surface_plot_show:
        fig = plt.figure(num=figure_name)
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(max(float(settings.live_surface_plot_pause), 0.001))
    else:
        plt.close(plt.figure(num=figure_name))


def _maybe_update_live_plots(
    rmb,
    settings: CostAwareSurfaceSettings,
    budget,
    iteration: int,
    volume_history: list[LiveVolumePoint],
    *,
    pending_requests=None,
) -> None:
    if (
        pending_requests is None
        and _measured_config_count(rmb._data) < settings.posterior_min_configs
    ):
        return
    try:
        params, weights = _live_plot_params(rmb._data, settings, pending_requests)
        if params is None or weights is None:
            return
        if settings.live_surface_plot and iteration % max(1, settings.live_surface_plot_every) == 0:
            mesh = _surface_plot_mesh(params, weights, settings)
            if mesh is not None:
                _plot_boundary_surface_with_pending(
                    mesh,
                    _measured_arrays(rmb._data),
                    settings,
                    pending_requests=pending_requests,
                )
        if (
            pending_requests is None
            and settings.live_volume_plot
            and iteration % max(1, settings.live_volume_plot_every) == 0
        ):
            _record_volume_point_from_posterior(
                params, weights, settings, iteration, volume_history
            )
            plot_live_volume_history(
                volume_history,
                settings,
                png_path=settings.live_volume_plot_path,
                show=settings.live_volume_plot_show,
                show_block=False,
                show_pause=settings.live_volume_plot_pause,
                close=not settings.live_volume_plot_show,
                figure_name="cost-aware live volume",
            )
    except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
        print_progress(settings, budget, f"live plotting skipped: {exc}")


def run_with_design_plots(settings: CostAwareSurfaceSettings):
    volume_history = _restore_volume_history_from_meta(settings)
    original_save_checkpoint = _method._save_measurement_checkpoint

    def save_checkpoint_with_live_volume(*args, **kwargs):
        original_save_checkpoint(*args, **kwargs)
        _save_volume_history_to_meta(settings, volume_history)

    def progress_callback(rmb, callback_settings, budget, iteration, pending_requests=None):
        _maybe_update_live_plots(
            rmb,
            callback_settings,
            budget,
            iteration,
            volume_history,
            pending_requests=pending_requests,
        )

    _method._save_measurement_checkpoint = save_checkpoint_with_live_volume
    try:
        rmb, configs, budget = run_with_budget(settings, progress_callback=progress_callback)
    finally:
        _method._save_measurement_checkpoint = original_save_checkpoint
    _save_volume_history_to_meta(settings, volume_history)

    params, weights = _stateless_posterior(rmb._data, settings)
    if params is not None and weights is not None:
        mesh = _surface_plot_mesh(params, weights, settings)
        if mesh is not None:
            plot_boundary_surface(
                *mesh,
                _measured_arrays(rmb._data),
                settings,
                png_path=settings.surface_plot_path,
                show=settings.surface_plot_show,
                show_block=True,
                figure_name="cost-aware final surface",
            )
    return rmb, configs, budget


if __name__ == "__main__":
    q_values = tuple(range(20, 51, 1))
    settings = CostAwareSurfaceSettings(
        q_values=q_values,
        n_qubits=round(sum(q_values) / len(q_values)),
        acquisition_q_resolution=30,
        acquisition_ratio_points=10,
        backend_model="sympleq",
        hqc_budget=500.0,
        plot=False,
        surface_plot_show=True,
        surface_plot_n_gates_bounds=(100, 3000),
        gp_grid_surface_path=(
            EXPERIMENTS_DIR.parent
            / "rmb_data"
            / "FLE_20260702_063751_gp_grid_3d.npz"
        ),
        gp_grid_surface_label="Rick/Shreya Grid",
        rng_seed=1234,
        verbose=True,
        use_scrambler=True,
        live_surface_plot=True,
        live_surface_plot_show=True,
        live_surface_plot_pause=0.5,
        live_volume_plot=True,
        live_volume_plot_show=True,
        live_volume_plot_pause=0.5,
        max_qubit_window=5,
        min_distinct_q_coverage=10,
        grid_resolution=(11, 11, 7, 7, 1, 5, 1),
        boundary_fit_resolution=(17, 17, 9, 9, 1, 7, 1),
        continuous_refit=True,
    )
    rmb, configs, budget = run_with_design_plots(settings)
