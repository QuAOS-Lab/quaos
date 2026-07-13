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
    plot_boundary_uncertainty_surface,
    plot_live_volume_history,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import *  # noqa: F401,F403
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import (
    EXPERIMENTS_DIR,
    CostAwareSurfaceSettings,
    _LN2,
    _I_L1,
    _I_L2,
    _analytic_lindblad_gates,
    _asymptote,
    _measured_arrays,
    _measured_config_count,
    _stateless_posterior,
    _surface_s1_s2,
    _surface_plot_mesh,
    _logn_mean_sigma_mesh,
    run_with_budget,
)

LiveVolumePoint = tuple[float, ...]


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
    secondary_reference_volume: float | None = None,
    fitted_l1_rate: float | None = None,
    fitted_l2_rate: float | None = None,
    fitted_l1_rate_std: float | None = None,
    fitted_l2_rate_std: float | None = None,
) -> None:
    point = (
        int(iteration),
        float(fitted_volume),
        float(true_volume),
        float(lower_volume),
        float(upper_volume),
        float(secondary_reference_volume)
        if secondary_reference_volume is not None
        else float("nan"),
        float(fitted_l1_rate) if fitted_l1_rate is not None else float("nan"),
        float(fitted_l2_rate) if fitted_l2_rate is not None else float("nan"),
        float(fitted_l1_rate_std)
        if fitted_l1_rate_std is not None
        else float("nan"),
        float(fitted_l2_rate_std)
        if fitted_l2_rate_std is not None
        else float("nan"),
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
    secondary_reference_gates=None,
) -> float | None:
    """Record the fitted/reference volume point; optionally return a second
    reference volume (e.g. analytic Lindblad) computed in the SAME
    ``_surface_s1_s2`` call so the fitted posterior mesh is built only once."""
    scores = _surface_s1_s2(
        params,
        weights,
        settings,
        _volume_reference_gates(settings),
        secondary_reference_gates=secondary_reference_gates,
    )
    fitted_volume = float(scores["surface_volume_fit"])
    true_volume = float(scores["surface_volume_reference"])
    if not (np.isfinite(fitted_volume) and np.isfinite(true_volume)):
        return None
    secondary = scores.get("surface_volume_secondary_reference")
    secondary_volume = (
        float(secondary)
        if secondary is not None and np.isfinite(float(secondary))
        else None
    )
    fitted_l1_rate = float(np.sum(weights * params[:, _I_L1]))
    fitted_l2_rate = float(np.sum(weights * params[:, _I_L2]))
    fitted_l1_rate_std = float(
        np.sqrt(np.sum(weights * (params[:, _I_L1] - fitted_l1_rate) ** 2))
    )
    fitted_l2_rate_std = float(
        np.sqrt(np.sum(weights * (params[:, _I_L2] - fitted_l2_rate) ** 2))
    )
    _record_live_volume_point(
        history,
        iteration,
        fitted_volume,
        true_volume,
        float(scores["surface_volume_lower_1sigma"]),
        float(scores["surface_volume_upper_1sigma"]),
        secondary_volume,
        fitted_l1_rate,
        fitted_l2_rate,
        fitted_l1_rate_std,
        fitted_l2_rate_std,
    )
    return secondary_volume


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
    if int(meta.get("live_volume_history_version", 0) or 0) >= 4:
        for row in meta.get("live_volume_history", []):
            try:
                history.append(tuple(row))
            except TypeError:
                pass
    if history:
        restored: list[LiveVolumePoint] = []
        for row in history:
            if len(row) < 5:
                continue
            restored.append(
                (
                    int(row[0]),
                    float(row[1]),
                    float(row[2]),
                    float(row[3]),
                    float(row[4]),
                    float(row[5]) if len(row) > 5 else float("nan"),
                    float(row[6]) if len(row) > 6 else float("nan"),
                    float(row[7]) if len(row) > 7 else float("nan"),
                    float(row[8]) if len(row) > 8 else float("nan"),
                    float(row[9]) if len(row) > 9 else float("nan"),
                )
            )
        return restored
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
            secondary_reference_gates = (
                _analytic_lindblad_gates
                if settings.gp_grid_surface_path is not None
                and getattr(settings, "surface_plot_analytic", True)
                else None
            )
            _record_volume_point_from_posterior(
                params,
                weights,
                settings,
                int(iteration),
                history,
                secondary_reference_gates=secondary_reference_gates,
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
        meta["live_volume_history_version"] = 4
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
            marker="x",
            s=72,
            linewidth=1.6,
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
    except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
        print_progress(settings, budget, f"live plot posterior skipped: {exc}")
        return
    if params is None or weights is None:
        return

    # Surface plot and volume plot are guarded separately, so a failure in one
    # is reported on its own and cannot silently hide behind the other.
    if settings.live_surface_plot and iteration % max(1, settings.live_surface_plot_every) == 0:
        try:
            mesh = _surface_plot_mesh(params, weights, settings)
            if mesh is not None:
                _plot_boundary_surface_with_pending(
                    mesh,
                    _measured_arrays(rmb._data),
                    settings,
                    pending_requests=pending_requests,
                )
        except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
            print_progress(settings, budget, f"live surface plot skipped: {exc}")

    # Volume history is only advanced on the settled (non-pending) callback.
    if (
        pending_requests is None
        and settings.live_volume_plot
        and iteration % max(1, settings.live_volume_plot_every) == 0
    ):
        try:
            # A second (analytic) reference line is only distinct from the primary
            # when the primary reference is a GP grid; otherwise the crimson line
            # already IS the analytic volume.  It is computed inside the single
            # _surface_s1_s2 call in the recorder (reusing the fitted mesh), so
            # there is one volume definition and no double compute.
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
            if not volume_history:
                # The plot needs at least one finite point; if we have none yet,
                # say so instead of silently drawing nothing (the usual cause is
                # the boundary lying outside n_gates_bounds so the volume is NaN).
                print_progress(
                    settings,
                    budget,
                    "live volume plot: no finite volume point yet "
                    "(boundary outside n_gates_bounds, or fit not yet started)",
                )
            else:
                png = plot_live_volume_history(
                    volume_history,
                    settings,
                    png_path=settings.live_volume_plot_path,
                    show=settings.live_volume_plot_show,
                    show_block=False,
                    show_pause=settings.live_volume_plot_pause,
                    close=not settings.live_volume_plot_show,
                    figure_name="cost-aware live volume",
                )
                if png is not None:
                    print_progress(settings, budget, f"live volume plot -> {png}")
        except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
            print_progress(settings, budget, f"live volume plot skipped: {exc}")


def run_with_design_plots(settings: CostAwareSurfaceSettings):
    volume_history = _restore_volume_history_from_meta(settings)

    def progress_callback(rmb, callback_settings, budget, iteration, pending_requests=None):
        _maybe_update_live_plots(
            rmb,
            callback_settings,
            budget,
            iteration,
            volume_history,
            pending_requests=pending_requests,
        )

    def on_checkpoint():
        # Persist the live-volume sidecar into the same metadata file right after
        # each measurement checkpoint (via run_with_budget's on_checkpoint hook,
        # so no module-level function needs to be patched).
        _save_volume_history_to_meta(settings, volume_history)

    rmb, configs, budget = run_with_budget(
        settings,
        progress_callback=progress_callback,
        on_checkpoint=on_checkpoint,
    )
    _save_volume_history_to_meta(settings, volume_history)

    params, weights = _stateless_posterior(rmb._data, settings)
    if params is not None and weights is not None:
        # Build a fitted mean and +-sigma mesh (converted to log10) and draw
        # the uncertainty surface into the main COST_AWARE surface plot so the
        # final plot includes the +-sigma envelopes.
        try:
            r_lo, r_hi = settings.ratio_bounds
            ratios = np.linspace(r_lo, r_hi, max(settings.score_ratio_points, 2))
            qubits = np.asarray(settings.q_values, dtype=float)
            rr, qq = np.meshgrid(ratios, qubits)
            mean_log, sigma_log, _ = _logn_mean_sigma_mesh(
                params, weights, settings, rr, qq, raw=getattr(settings, "plot_raw_fidelity_boundary", False)
            )
            if mean_log is not None:
                ln10 = np.log(10.0)
                mean_log10 = mean_log / ln10
                sigma_log10 = sigma_log / ln10
                k = float(getattr(settings, "surface_uncertainty_sigma", 1.0))
                lower_log10 = mean_log10 - k * sigma_log10
                upper_log10 = mean_log10 + k * sigma_log10
                png = plot_boundary_uncertainty_surface(
                    rr,
                    qq,
                    mean_log10,
                    lower_log10,
                    upper_log10,
                    _measured_arrays(rmb._data),
                    settings,
                    png_path=settings.surface_plot_path,
                    show=settings.surface_plot_show,
                    show_block=True,
                    figure_name="cost-aware final surface",
                )
        except Exception as exc:  # noqa: BLE001 - final plotting should not kill a run
            print_progress(settings, budget, f"final surface plot skipped: {exc}")
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
