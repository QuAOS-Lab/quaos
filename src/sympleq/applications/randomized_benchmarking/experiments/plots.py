"""
Plot functions shared by the boundary experiments in this package.

All plots live in (total gates, two-qubit gate ratio) space, one subplot per
``n_qubits``, so the different experiments can be compared directly.
"""
from __future__ import annotations

from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.special import betainc

from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    CrossingSettings,
    MonotoneFidelitySurface,
    bootstrap_surfaces,
    contour_points_from_surface,
    grouped_by_n_qubits,
    measured_gate_count_bounds,
    measured_items,
    try_fit_monotone_fidelity_surface,
)


REFERENCE_LINE_LABEL = "reference line"
ANALYTIC_LINE_LABEL = "analytic Lindblad line"
PARAMETRIC_BOUNDARY_LABEL = "parametric boundary fit"
PARAMETRIC_BOUNDARY_BAND_LABEL = "parametric fit 90% band"
REFERENCE_NUMERATOR = 0.7106
REFERENCE_OFFSET = 1.91e-4
REFERENCE_SLOPE = 3.65e-3


def posterior_above(estimator: BayesianEstimator) -> float:
    """Posterior probability that the Boolean success probability is above 0.5."""
    alpha, beta = estimator.posterior_alpha_beta()
    return float(1.0 - betainc(alpha, beta, 0.5))


def reference_gate_counts(ratios: np.ndarray) -> np.ndarray:
    return REFERENCE_NUMERATOR / (REFERENCE_OFFSET + REFERENCE_SLOPE * ratios)


def analytic_gate_counts(
    ratios: np.ndarray,
    *,
    one_q_noise_scale: float = 1.0,
    two_q_noise_scale: float = 1.0,
) -> np.ndarray:
    return np.log(2.0) / (
        REFERENCE_OFFSET * one_q_noise_scale
        + REFERENCE_SLOPE * two_q_noise_scale * ratios
    )


def analytic_noise_scales(settings: CrossingSettings) -> tuple[float, float]:
    return (
        float(getattr(settings, "one_q_noise_scale", 1.0)),
        float(getattr(settings, "two_q_noise_scale", 1.0)),
    )


def plot_reference_total_ratio_line(ax, settings: CrossingSettings) -> None:
    """Overlay the reference curve on an axis with x=total gates, y=ratio."""
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 400)
    gates = reference_gate_counts(ratios)
    mask = (
        np.isfinite(gates)
        & (gates > 0.0)
        & (settings.n_gates_bounds[0] <= gates)
        & (gates <= settings.n_gates_bounds[1])
    )
    if not np.any(mask):
        return
    ax.plot(
        gates[mask],
        ratios[mask],
        color="tab:purple",
        linestyle="--",
        linewidth=1.8,
        label=REFERENCE_LINE_LABEL,
        zorder=5,
    )


def plot_reference_gate_plane_line(
    ax,
    one_q_bounds: tuple[int, int],
    two_q_bounds: tuple[int, int],
    *,
    ratio_bounds: tuple[float, float] = (0.0, 1.0),
) -> None:
    """Overlay the reference curve on an axis with x=1q gates, y=2q gates."""
    ratios = np.linspace(ratio_bounds[0], ratio_bounds[1], 400)
    total_gates = reference_gate_counts(ratios)
    one_q_gates = total_gates * (1.0 - ratios)
    two_q_gates = total_gates * ratios
    mask = (
        np.isfinite(total_gates)
        & (total_gates > 0.0)
        & (one_q_bounds[0] <= one_q_gates)
        & (one_q_gates <= one_q_bounds[1])
        & (two_q_bounds[0] <= two_q_gates)
        & (two_q_gates <= two_q_bounds[1])
    )
    if not np.any(mask):
        return
    ax.plot(
        one_q_gates[mask],
        two_q_gates[mask],
        color="tab:purple",
        linestyle="--",
        linewidth=1.8,
        label=REFERENCE_LINE_LABEL,
        zorder=7,
    )


def plot_analytic_total_ratio_line(ax, settings: CrossingSettings) -> None:
    """Overlay the analytic exponential-decay contour on total-gate axes."""
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 400)
    one_q_noise_scale, two_q_noise_scale = analytic_noise_scales(settings)
    gates = analytic_gate_counts(
        ratios,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    mask = (
        np.isfinite(gates)
        & (gates > 0.0)
        & (settings.n_gates_bounds[0] <= gates)
        & (gates <= settings.n_gates_bounds[1])
    )
    if not np.any(mask):
        return
    ax.plot(
        gates[mask],
        ratios[mask],
        color="tab:green",
        linestyle=":",
        linewidth=2.0,
        label=ANALYTIC_LINE_LABEL,
        zorder=5,
    )


def plot_analytic_gate_plane_line(
    ax,
    settings: CrossingSettings,
    one_q_bounds: tuple[int, int],
    two_q_bounds: tuple[int, int],
    *,
    ratio_bounds: tuple[float, float] = (0.0, 1.0),
) -> None:
    """Overlay the analytic exponential-decay contour on gate-plane axes."""
    ratios = np.linspace(ratio_bounds[0], ratio_bounds[1], 400)
    one_q_noise_scale, two_q_noise_scale = analytic_noise_scales(settings)
    total_gates = analytic_gate_counts(
        ratios,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    one_q_gates = total_gates * (1.0 - ratios)
    two_q_gates = total_gates * ratios
    mask = (
        np.isfinite(total_gates)
        & (total_gates > 0.0)
        & (one_q_bounds[0] <= one_q_gates)
        & (one_q_gates <= one_q_bounds[1])
        & (two_q_bounds[0] <= two_q_gates)
        & (two_q_gates <= two_q_bounds[1])
    )
    if not np.any(mask):
        return
    ax.plot(
        one_q_gates[mask],
        two_q_gates[mask],
        color="tab:green",
        linestyle=":",
        linewidth=2.0,
        label=ANALYTIC_LINE_LABEL,
        zorder=7,
    )


def _fit_parametric_boundary_arrays(
    ratios: np.ndarray,
    n_gates: np.ndarray,
    successes: np.ndarray,
    failures: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float, float] | None:
    """
    Fit ``P(success) = sigmoid(k * (boundary(ratio) - n_gates) / boundary(ratio))``.

    The boundary has the inverse form ``n_gates = 1 / (q + s * ratio)``.
    """
    from scipy.optimize import minimize
    from scipy.special import expit

    if len(ratios) < 3:
        return None

    ratios = np.asarray(ratios, dtype=float)
    n_gates = np.asarray(n_gates, dtype=float)
    successes = np.asarray(successes, dtype=float)
    failures = np.asarray(failures, dtype=float)
    weights = np.asarray(weights, dtype=float)

    q0 = REFERENCE_OFFSET / REFERENCE_NUMERATOR
    s0 = REFERENCE_SLOPE / REFERENCE_NUMERATOR
    x0 = np.log([q0, s0, 4.0])

    def loss(log_params: np.ndarray) -> float:
        q, slope, sharpness = np.exp(log_params)
        boundary = 1.0 / (q + slope * ratios)
        probabilities = expit(sharpness * (boundary - n_gates) / boundary)
        eps = 1e-12
        negative_log_likelihood = -np.sum(
            weights * (
                successes * np.log(probabilities + eps)
                + failures * np.log(1.0 - probabilities + eps)
            )
        )
        # Lightly stabilize the boundary scale without forcing the reference line.
        regularization = 0.1 * np.sum((log_params[:2] - x0[:2]) ** 2)
        return float(negative_log_likelihood + regularization)

    bounds = [
        (np.log(1e-6), np.log(1.0)),
        (np.log(1e-6), np.log(1.0)),
        (np.log(0.05), np.log(100.0)),
    ]
    result = minimize(loss, x0=x0, bounds=bounds, method="L-BFGS-B")
    if not result.success:
        return None
    q, slope, sharpness = np.exp(result.x)
    return float(q), float(slope), float(sharpness)


def parametric_boundary_fit(
    data: RMBData,
) -> tuple[float, float, float] | None:
    """Fit the inverse-form parametric boundary to measured Bernoulli outcomes."""
    measured = measured_items(data)
    if len(measured) < 3:
        return None

    ratios = np.asarray([config.ratio_2_qb_gates for config, _ in measured],
                        dtype=float)
    n_gates = np.asarray([config.n_gates for config, _ in measured], dtype=float)
    successes = []
    failures = []
    weights = []
    for _, estimator in measured:
        counts = estimator.counts()
        successes.append(float(counts.get(True, 0)))
        failures.append(float(counts.get(False, 0)))
        runs = max(1, estimator.num_runs())
        variance = float(estimator.posterior_variance())
        variance_certainty = 1.0 - np.clip(variance / (1.0 / 12.0), 0.0, 1.0)
        shot_certainty = min(1.0, runs / 3.0)
        weights.append(max(0.05, shot_certainty * variance_certainty))
    successes = np.asarray(successes, dtype=float)
    failures = np.asarray(failures, dtype=float)
    weights = np.asarray(weights, dtype=float)
    return _fit_parametric_boundary_arrays(ratios, n_gates, successes, failures, weights)


def parametric_boundary_bootstrap(
    data: RMBData,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
) -> list[tuple[float, float, float]]:
    """Bootstrap inverse-boundary fits from per-config Beta posterior draws."""
    measured = measured_items(data)
    if len(measured) < 3:
        return []

    rng = np.random.default_rng(seed)
    ratios = np.asarray([config.ratio_2_qb_gates for config, _ in measured],
                        dtype=float)
    n_gates = np.asarray([config.n_gates for config, _ in measured], dtype=float)
    posteriors = [estimator.posterior_alpha_beta() for _, estimator in measured]
    alpha = np.asarray([a for a, _ in posteriors], dtype=float)
    beta = np.asarray([b for _, b in posteriors], dtype=float)
    weights = np.asarray([max(1, estimator.num_runs()) for _, estimator in measured],
                         dtype=float)

    fits = []
    for _ in range(n_bootstrap):
        sampled = rng.beta(alpha, beta)
        fit = _fit_parametric_boundary_arrays(
            ratios,
            n_gates,
            sampled,
            1.0 - sampled,
            weights,
        )
        if fit is not None:
            fits.append(fit)
    return fits


def parametric_boundary_gate_samples(
    fits: list[tuple[float, float, float]],
    ratios: np.ndarray,
) -> np.ndarray:
    """Total-gate boundary samples for fitted inverse-boundary parameters."""
    if not fits:
        return np.empty((0, len(ratios)), dtype=float)
    samples = []
    for q, slope, _ in fits:
        gates = 1.0 / (q + slope * ratios)
        gates[~np.isfinite(gates)] = np.nan
        gates[gates <= 0.0] = np.nan
        samples.append(gates)
    return np.asarray(samples, dtype=float)


def plot_parametric_boundary_total_ratio_band(
    ax,
    data: RMBData,
    settings: CrossingSettings,
    *,
    n_bootstrap: int = 100,
) -> None:
    """Overlay bootstrap confidence bands for the inverse-form boundary."""
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 160)
    samples = parametric_boundary_gate_samples(
        parametric_boundary_bootstrap(data, n_bootstrap=n_bootstrap,
                                      seed=settings.rng_seed),
        ratios,
    )
    if len(samples) == 0:
        return
    q05, q25, q75, q95 = np.nanquantile(samples, [0.05, 0.25, 0.75, 0.95], axis=0)
    valid90 = (
        np.isfinite(q05)
        & np.isfinite(q95)
        & (settings.n_gates_bounds[0] <= q05)
        & (q95 <= settings.n_gates_bounds[1])
    )
    valid50 = (
        np.isfinite(q25)
        & np.isfinite(q75)
        & (settings.n_gates_bounds[0] <= q25)
        & (q75 <= settings.n_gates_bounds[1])
    )
    if np.any(valid90):
        ax.fill_betweenx(
            ratios[valid90],
            q05[valid90],
            q95[valid90],
            color="tab:orange",
            alpha=0.12,
            linewidth=0,
            label=PARAMETRIC_BOUNDARY_BAND_LABEL,
            zorder=3,
        )
    if np.any(valid50):
        ax.fill_betweenx(
            ratios[valid50],
            q25[valid50],
            q75[valid50],
            color="tab:orange",
            alpha=0.23,
            linewidth=0,
            label="parametric fit 50% band",
            zorder=4,
        )


def plot_parametric_boundary_gate_plane_band(
    ax,
    data: RMBData,
    settings: CrossingSettings,
    one_q_bounds: tuple[int, int],
    two_q_bounds: tuple[int, int],
    *,
    n_bootstrap: int = 100,
) -> None:
    """Overlay bootstrap confidence bands for the inverse boundary in gate space."""
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 160)
    samples = parametric_boundary_gate_samples(
        parametric_boundary_bootstrap(data, n_bootstrap=n_bootstrap,
                                      seed=settings.rng_seed),
        ratios,
    )
    if len(samples) == 0:
        return
    q05, q25, q75, q95 = np.nanquantile(samples, [0.05, 0.25, 0.75, 0.95], axis=0)

    def gate_plane(gates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return gates * (1.0 - ratios), gates * ratios

    for lower, upper, alpha, label, zorder in (
        (q05, q95, 0.12, PARAMETRIC_BOUNDARY_BAND_LABEL, 3),
        (q25, q75, 0.23, "parametric fit 50% band", 4),
    ):
        lower_x, lower_y = gate_plane(lower)
        upper_x, upper_y = gate_plane(upper)
        valid = (
            np.isfinite(lower_x)
            & np.isfinite(lower_y)
            & np.isfinite(upper_x)
            & np.isfinite(upper_y)
            & (one_q_bounds[0] <= lower_x)
            & (lower_x <= one_q_bounds[1])
            & (one_q_bounds[0] <= upper_x)
            & (upper_x <= one_q_bounds[1])
            & (two_q_bounds[0] <= lower_y)
            & (lower_y <= two_q_bounds[1])
            & (two_q_bounds[0] <= upper_y)
            & (upper_y <= two_q_bounds[1])
        )
        if not np.any(valid):
            continue
        xs = np.concatenate([lower_x[valid], upper_x[valid][::-1]])
        ys = np.concatenate([lower_y[valid], upper_y[valid][::-1]])
        ax.fill(xs, ys, color="tab:orange", alpha=alpha, linewidth=0,
                label=label, zorder=zorder)


def plot_parametric_boundary_total_ratio_line(
    ax,
    data: RMBData,
    settings: CrossingSettings,
) -> None:
    """Overlay the fitted inverse-form boundary on x=total gates, y=ratio."""
    plot_parametric_boundary_total_ratio_band(ax, data, settings)
    fit = parametric_boundary_fit(data)
    if fit is None:
        return
    q, slope, _ = fit
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 400)
    gates = 1.0 / (q + slope * ratios)
    mask = (
        np.isfinite(gates)
        & (gates > 0.0)
        & (settings.n_gates_bounds[0] <= gates)
        & (gates <= settings.n_gates_bounds[1])
    )
    if not np.any(mask):
        return
    ax.plot(
        gates[mask],
        ratios[mask],
        color="tab:orange",
        linestyle="-",
        linewidth=1.9,
        label=PARAMETRIC_BOUNDARY_LABEL,
        zorder=6,
    )


def plot_parametric_boundary_gate_plane_line(
    ax,
    data: RMBData,
    settings: CrossingSettings,
    one_q_bounds: tuple[int, int],
    two_q_bounds: tuple[int, int],
) -> None:
    """Overlay the fitted inverse-form boundary on x=1q gates, y=2q gates."""
    plot_parametric_boundary_gate_plane_band(
        ax,
        data,
        settings,
        one_q_bounds,
        two_q_bounds,
    )
    fit = parametric_boundary_fit(data)
    if fit is None:
        return
    q, slope, _ = fit
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 400)
    total_gates = 1.0 / (q + slope * ratios)
    one_q_gates = total_gates * (1.0 - ratios)
    two_q_gates = total_gates * ratios
    mask = (
        np.isfinite(total_gates)
        & (total_gates > 0.0)
        & (one_q_bounds[0] <= one_q_gates)
        & (one_q_gates <= one_q_bounds[1])
        & (two_q_bounds[0] <= two_q_gates)
        & (two_q_gates <= two_q_bounds[1])
    )
    if not np.any(mask):
        return
    ax.plot(
        one_q_gates[mask],
        two_q_gates[mask],
        color="tab:orange",
        linestyle="-",
        linewidth=1.9,
        label=PARAMETRIC_BOUNDARY_LABEL,
        zorder=8,
    )


def merge_close_configs(data: RMBData, n_1qb_gates_bin: int = 10, n_2qb_gates_bin: int = 10) -> RMBData:
    """Return a copy of ``data`` with close configs collapsed and their estimators merged.

    Parameters
    ----------
    data : RMBData
        Mapping from configurations to their Bayesian estimators.
    depth_bin : int
        Configs' ``depth`` is snapped to the nearest multiple of this value
        (with a floor of ``depth_bin`` so the dataclass ``depth >= 1``
        validation is preserved).
    ratio_digits : int
        Number of decimal places to round the two-qudit ratio bounds to.

    Returns
    -------
    RMBData
        New mapping keyed by coarsened configs.
    """
    merged: RMBData = {}
    for config, estimator in data.items():
        coarse_n_1qb_gates = max(n_1qb_gates_bin, round(config.n_1qb_gates / n_1qb_gates_bin) * n_1qb_gates_bin)
        coarse_n_2qb_gates = max(n_2qb_gates_bin, round(config.n_2qb_gates / n_2qb_gates_bin) * n_2qb_gates_bin)
        coarse = replace(
            config,
            n_1qb_gates=coarse_n_1qb_gates,
            n_2qb_gates=coarse_n_2qb_gates
        )
        target = merged.setdefault(coarse, BayesianEstimator(
            threshold=estimator.threshold,
            min_runs=estimator.min_runs,
            max_runs=estimator.max_runs,
        ))
        target.merge(estimator)

    return merged


def plot_data(data: RMBData, axes=None, show: bool = True, skip_incomplete: bool = True,
              level_line: list[RMBConfig] | None = None, log_x: bool = False,
              n_1qb_gates_bin: int | None = 20, n_2qb_gates_bin: int | None = 10) -> list:
    """Scatter fidelity per ``n_qubits``: x = total gates, y = two-qudit ratio, color = fidelity.

    Parameters
    ----------
    data : RMBData
        Mapping from configurations to their Bayesian estimators.
    axes : Sequence[matplotlib.axes.Axes] | None
        Axes to plot on, one per distinct ``n_qubits`` (sorted
        ascending). If ``None``, a new figure with one subplot per
        ``n_qubits`` is created.
    show : bool
        If ``True``, call ``plt.show()`` after building the plot.
    skip_incomplete : bool
        If ``True``, drop any estimator that has not converged.
        ``n_qubits`` groups left empty after filtering are not given a
        subplot.
    level_line : list[RMBConfig] | None
        Configs lying on the fidelity = 0.5 line. They are drawn as a
        black line on the subplot matching their ``n_qubits``, in the
        same (total gates, two-qubit gate ratio) coordinates as the data.
    log_x : bool
        If ``True``, use a logarithmic gate-count axis. Useful when the
        data spans orders of magnitude in circuit size.
    n_1qb_gates_bin, n_2qb_gates_bin : int | None
        Gate-count bin sizes used to merge close configs before plotting
        (see :func:`merge_close_configs`). Either being ``None`` skips the
        merge and plots every measured config at full resolution.

    Returns
    -------
    list[matplotlib.axes.Axes]
        The axes the plots were drawn on (sorted by ``n_qubits``).
    """
    import matplotlib.pyplot as plt
    from matplotlib import patheffects
    from matplotlib.colors import LinearSegmentedColormap

    if n_1qb_gates_bin is not None and n_2qb_gates_bin is not None:
        data = merge_close_configs(data, n_1qb_gates_bin=n_1qb_gates_bin,
                                   n_2qb_gates_bin=n_2qb_gates_bin)
    if skip_incomplete:
        data = {config: estimator for config, estimator in data.items()
                if estimator.is_converged()}
    groups = grouped_by_n_qubits(data)

    if not groups:
        return []

    sorted_groups = sorted(groups.items())
    n_groups = len(sorted_groups)

    if axes is None:
        _, axes_arr = plt.subplots(
            1, n_groups, figsize=(6.5 * n_groups, 4.8), squeeze=False,
            constrained_layout=True)
        axes_list = list(axes_arr[0])
    else:
        axes_list = list(axes)
        if len(axes_list) < n_groups:
            raise ValueError(
                f"Need at least {n_groups} axes for {n_groups} n_qubits "
                f"groups, got {len(axes_list)}.")

    cmap = LinearSegmentedColormap.from_list(
        "darkred_to_lime",
        ["darkred", "red", "orange", "lime", "green"])

    for ax, (n_qubits, group) in zip(axes_list, sorted_groups):
        n_gates = np.array([c.n_gates for c in group])
        ratios = np.array(
            [c.ratio_2_qb_gates for c in group])
        fidelities = np.array([e.probability(True) for e in group.values()])
        stds = np.array([np.sqrt(e.variance(True)) for e in group.values()])

        outer_color = np.clip(fidelities - stds, 0.0, 1.0)
        middle_color = fidelities
        inner_color = np.clip(fidelities + stds, 0.0, 1.0)

        ax.scatter(n_gates, ratios, c=outer_color, cmap=cmap,
                   vmin=0.0, vmax=1.0, s=210, edgecolors="none",
                   alpha=0.85, zorder=1)
        sc = ax.scatter(n_gates, ratios, c=middle_color, cmap=cmap,
                        vmin=0.0, vmax=1.0, s=100, edgecolors="white",
                        linewidths=0.5, zorder=2)
        ax.scatter(n_gates, ratios, c=inner_color, cmap=cmap,
                   vmin=0.0, vmax=1.0, s=28, edgecolors="none", zorder=3)

        if level_line:
            points = sorted(
                (c.n_gates, c.ratio_2_qb_gates)
                for c in level_line if c.n_qubits == n_qubits
            )
            if points:
                xs, ys = zip(*points)
                line, = ax.plot(xs, ys, color="black", marker="o",
                                markersize=5, markerfacecolor="white",
                                markeredgecolor="black", linewidth=1.8,
                                label="fidelity = 0.5", zorder=5)
                line.set_path_effects([
                    patheffects.Stroke(linewidth=3.5, foreground="white"),
                    patheffects.Normal(),
                ])
                ax.legend(loc="best", frameon=True, framealpha=0.9)

        if log_x:
            ax.set_xscale("log")
        ax.grid(True, which="major", alpha=0.3, linewidth=0.6)
        ax.grid(True, which="minor", alpha=0.12, linewidth=0.4)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("# Gates", fontsize=11)
        ax.set_ylabel("Two-qudit gate ratio", fontsize=11)
        ax.set_title(f"# Qubits = {n_qubits}", fontsize=13, fontweight="bold")
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label("Fidelity", fontsize=11)
        # matplotlib's stubs mistype Colorbar.outline as a Spines mapping.
        cbar.outline.set_visible(False)  # type: ignore

    if show:
        plt.show()

    return axes_list


def plot_level_line(
    data: RMBData,
    crossings: list[RMBConfig],
    *,
    contour: np.ndarray | None = None,
    surfaces: list[MonotoneFidelitySurface] | None = None,
    settings: CrossingSettings | None = None,
    axes=None,
    png_path: str | Path | None = None,
    show: bool = True,
    n_1qb_gates_bin: int | None = 20,
    n_2qb_gates_bin: int | None = 10,
) -> list:
    """
    Scatter the recorded data with the traced fidelity = 0.5 line on top.

    Parameters
    ----------
    data : RMBData
        Mapping from configurations to their Bayesian estimators.
    crossings : list[RMBConfig]
        Configs lying on the fidelity = 0.5 line.
    contour : np.ndarray | None
        (total gates, ratio) points of a fitted p=0.5 contour (see
        :func:`~.common.contour_points_from_surface`) to overlay as a line.
    surfaces : list[MonotoneFidelitySurface] | None
        Bootstrap monotone surfaces. When provided with ``settings``, draw
        fixed-ratio contour quantile bands.
    settings : CrossingSettings | None
        Experiment settings. When given, draw the parametric boundary fit and
        reference curve in the same (total gates, ratio) coordinates.
    png_path : str | Path | None
        Where to save the figure. ``None`` skips saving.
    show : bool
        If ``True``, call ``plt.show()`` at the end.
    n_1qb_gates_bin, n_2qb_gates_bin : int | None
        Scatter bin sizes forwarded to :func:`plot_data`; either ``None``
        plots every measured config at full resolution.
    """
    import matplotlib.pyplot as plt

    axes = plot_data(data, axes=axes, show=False, skip_incomplete=False,
                     level_line=crossings, log_x=True,
                     n_1qb_gates_bin=n_1qb_gates_bin, n_2qb_gates_bin=n_2qb_gates_bin)
    if axes and contour is not None and len(contour) > 0:
        order = np.argsort(contour[:, 1])
        axes[0].plot(contour[order, 0], contour[order, 1], linestyle="-.",
                     color="0.35", linewidth=1.6, label="monotone fit p=0.5", zorder=4)
        axes[0].legend(loc="best", frameon=True, framealpha=0.9)
    if axes and settings is not None:
        if surfaces:
            plot_bootstrap_gate_bands(
                axes[0],
                settings,
                surfaces,
                total_ratio_axes=True,
            )
        plot_parametric_boundary_total_ratio_line(axes[0], data, settings)
        plot_analytic_total_ratio_line(axes[0], settings)
        plot_reference_total_ratio_line(axes[0], settings)
        axes[0].legend(loc="best", frameon=True, framealpha=0.9)
    if axes and png_path is not None:
        axes[0].figure.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    return axes


def surface_gate_crossings_at_ratios(
    surface: MonotoneFidelitySurface,
    settings: CrossingSettings,
    ratios: np.ndarray,
    *,
    level: float = 0.5,
) -> np.ndarray:
    """
    Interpolate total-gate crossings of one surface at fixed ratio values.

    Missing crossings are returned as ``nan``. This keeps the bootstrap
    uncertainty tied to the contour position at each ratio instead of averaging
    whole surfaces, which can create a displaced bootstrap center.
    """
    ratios = np.asarray(ratios, dtype=float)
    gates_axis = np.linspace(
        settings.n_gates_bounds[0],
        settings.n_gates_bounds[1],
        settings.candidate_grid_size[0],
    )
    crossings = np.full_like(ratios, np.nan, dtype=float)

    for i, ratio in enumerate(ratios):
        points = np.column_stack([gates_axis, np.full_like(gates_axis, ratio)])
        delta = surface.probability(points) - level
        exact = np.where(delta == 0.0)[0]
        if len(exact) > 0:
            crossings[i] = float(gates_axis[int(exact[0])])
            continue

        indices = np.where(delta[:-1] * delta[1:] < 0.0)[0]
        if len(indices) == 0:
            continue
        j = int(indices[0])
        denom = abs(delta[j]) + abs(delta[j + 1])
        if denom <= 0.0:
            crossings[i] = float(gates_axis[j])
            continue
        t = abs(delta[j]) / denom
        crossings[i] = float((1.0 - t) * gates_axis[j] + t * gates_axis[j + 1])

    return crossings


def bootstrap_gate_quantiles(
    surfaces: list[MonotoneFidelitySurface],
    settings: CrossingSettings,
    ratios: np.ndarray,
    *,
    level: float = 0.5,
    quantiles: tuple[float, ...] = (0.05, 0.25, 0.5, 0.75, 0.95),
) -> tuple[dict[float, np.ndarray], np.ndarray]:
    """
    Fixed-ratio bootstrap quantiles of the fitted contour.

    Returns a mapping from quantile to total-gate count and the valid fraction
    at each ratio. The valid fraction is the fraction of bootstrap surfaces
    whose fitted probability range contains ``level`` at that ratio.
    """
    ratios = np.asarray(ratios, dtype=float)
    empty = {q: np.full_like(ratios, np.nan, dtype=float) for q in quantiles}
    if not surfaces:
        return empty, np.zeros_like(ratios, dtype=float)

    samples = np.vstack([
        surface_gate_crossings_at_ratios(surface, settings, ratios, level=level)
        for surface in surfaces
    ])
    valid = np.isfinite(samples)
    valid_fraction = np.mean(valid, axis=0)
    output = {}
    for q in quantiles:
        values = np.full_like(ratios, np.nan, dtype=float)
        for i in range(len(ratios)):
            column = samples[:, i]
            column = column[np.isfinite(column)]
            if len(column) > 0:
                values[i] = float(np.quantile(column, q))
        output[q] = values
    return output, valid_fraction


def plot_bootstrap_gate_bands(
    ax,
    settings: CrossingSettings,
    surfaces: list[MonotoneFidelitySurface],
    *,
    ratios: np.ndarray | None = None,
    level: float = 0.5,
    total_ratio_axes: bool = True,
) -> None:
    """Plot fixed-ratio bootstrap contour quantile bands."""
    if ratios is None:
        ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], 160)
    quantiles, valid_fraction = bootstrap_gate_quantiles(
        surfaces,
        settings,
        ratios,
        level=level,
    )
    q05, q25, q75, q95 = (quantiles[q] for q in (0.05, 0.25, 0.75, 0.95))
    valid90 = np.isfinite(q05) & np.isfinite(q95)
    valid50 = np.isfinite(q25) & np.isfinite(q75)
    if not np.any(valid90):
        return

    if total_ratio_axes:
        ax.fill_betweenx(
            ratios[valid90],
            q05[valid90],
            q95[valid90],
            color="tab:blue",
            alpha=0.14,
            linewidth=0,
            label="bootstrap 90% contour band",
            zorder=1,
        )
        if np.any(valid50):
            ax.fill_betweenx(
                ratios[valid50],
                q25[valid50],
                q75[valid50],
                color="tab:blue",
                alpha=0.26,
                linewidth=0,
                label="bootstrap 50% contour band",
                zorder=2,
            )
    else:
        def gate_plane(gates: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            return gates * (1.0 - ratios), gates * ratios

        lower_x, lower_y = gate_plane(q05)
        upper_x, upper_y = gate_plane(q95)
        band_valid = valid90 & np.isfinite(lower_x) & np.isfinite(upper_x)
        if np.any(band_valid):
            xs = np.concatenate([lower_x[band_valid], upper_x[band_valid][::-1]])
            ys = np.concatenate([lower_y[band_valid], upper_y[band_valid][::-1]])
            ax.fill(xs, ys, color="tab:blue", alpha=0.14, linewidth=0,
                    label="bootstrap 90% contour band", zorder=1)

        lower_x, lower_y = gate_plane(q25)
        upper_x, upper_y = gate_plane(q75)
        band_valid = valid50 & np.isfinite(lower_x) & np.isfinite(upper_x)
        if np.any(band_valid):
            xs = np.concatenate([lower_x[band_valid], upper_x[band_valid][::-1]])
            ys = np.concatenate([lower_y[band_valid], upper_y[band_valid][::-1]])
            ax.fill(xs, ys, color="tab:blue", alpha=0.26, linewidth=0,
                    label="bootstrap 50% contour band", zorder=2)

    min_valid = float(np.min(valid_fraction[valid90])) if np.any(valid90) else 0.0
    ax.plot([], [], color="none", label=f"min bootstrap crossing rate {min_valid:.0%}")


def fidelity_colormap():
    """Colormap with the fidelity=0.5 midpoint shown in green."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "fidelity_red_green_purple",
        [
            (0.0, "#b2182b"),
            (0.5, "#1a9850"),
            (1.0, "#542788"),
        ],
    )


def plot_monotone_fidelity_surface_contours(
    data: RMBData,
    settings: CrossingSettings,
    *,
    surface: MonotoneFidelitySurface | None = None,
    axes=None,
    show: bool = True,
    axis_padding: float = 0.05,
    log_x: bool = False,
    label_runs: bool = True,
    run_label_limit: int = 300,
) -> list:
    """
    Scatter the data and draw the fitted monotone p=0.5 contour per ``n_qubits``.

    With ``surface``, the given fit is drawn for every group instead of
    fitting per group (the experiments in this package run at one fixed
    ``n_qubits``, so there is one group). With ``label_runs``, each point is
    annotated with its number of recorded measurements, unless a group holds
    more than ``run_label_limit`` points.
    """
    import matplotlib.pyplot as plt

    groups = sorted(grouped_by_n_qubits(data).items())
    if not groups:
        return []

    if axes is None:
        _, axes_arr = plt.subplots(1, len(groups), figsize=(5 * len(groups), 4), squeeze=False)
        axes = list(axes_arr[0])
    else:
        axes = list(axes)
        if len(axes) < len(groups):
            raise ValueError(f"Need at least {len(groups)} axes, got {len(axes)}.")
    cmap = fidelity_colormap()

    gates_span = settings.n_gates_bounds[1] - settings.n_gates_bounds[0]
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    x_pad = axis_padding * gates_span
    y_pad = axis_padding * ratio_span

    for ax, (n_qubits, group) in zip(axes, groups):
        n_gates = np.array([config.n_gates for config in group], dtype=float)
        ratios = np.array([config.ratio_2_qb_gates for config in group], dtype=float)
        fidelities = np.array([estimator.posterior_mean() for estimator in group.values()],
                              dtype=float)
        variances = np.array([estimator.posterior_variance() for estimator in group.values()],
                             dtype=float)
        certainty = 1.0 - np.clip(variances / (1.0 / 12.0), 0.0, 1.0)
        marker_sizes = 30.0 + 160.0 * certainty

        scatter = ax.scatter(
            n_gates,
            ratios,
            c=fidelities,
            cmap=cmap,
            vmin=0.0,
            vmax=1.0,
            s=marker_sizes,
            edgecolors="black",
            linewidths=0.4,
            zorder=3,
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label("Fidelity")
        ax.scatter([], [], s=30, facecolors="none", edgecolors="black", label="uncertain point")
        ax.scatter([], [], s=190, facecolors="none", edgecolors="black", label="certain point")

        if label_runs and len(group) <= run_label_limit:
            for config, estimator in group.items():
                ax.text(
                    float(config.n_gates),
                    float(config.ratio_2_qb_gates),
                    str(estimator.num_runs()),
                    ha="center",
                    va="center",
                    fontsize=6,
                    color="white",
                    weight="bold",
                    zorder=6,
                    clip_on=True,
                )

        group_surface = (surface if surface is not None
                         else try_fit_monotone_fidelity_surface(group, settings))
        if group_surface is not None:
            gates_grid, ratio_grid, probabilities = group_surface.probability_grid(
                settings.candidate_grid_size)
            if float(np.min(probabilities)) <= 0.5 <= float(np.max(probabilities)):
                ax.contour(
                    gates_grid,
                    ratio_grid,
                    probabilities,
                    levels=[0.5],
                    colors="black",
                    linewidths=2,
                    zorder=4,
                )
                ax.plot([], [], color="black", linewidth=2, label="monotone E[fidelity] = 0.5")

        plot_parametric_boundary_total_ratio_line(ax, group, settings)
        plot_analytic_total_ratio_line(ax, settings)
        plot_reference_total_ratio_line(ax, settings)

        if log_x:
            ax.set_xscale("log")
            x_lower = max(settings.n_gates_bounds[0] - x_pad,
                          0.8 * settings.n_gates_bounds[0])
        else:
            x_lower = settings.n_gates_bounds[0] - x_pad
        ax.set_xlabel("# Gates")
        ax.set_ylabel("Two-qudit gate ratio")
        ax.set_title(f"# Qubits = {n_qubits}")
        ax.set_xlim(x_lower, settings.n_gates_bounds[1] + x_pad)
        ax.set_ylim(
            max(0.0, settings.ratio_bounds[0] - y_pad),
            min(1.0, settings.ratio_bounds[1] + y_pad),
        )
        ax.legend(loc="best")

    if show:
        plt.show()

    return axes


def plot_monotone_fidelity_surface_with_confidence(
    data: RMBData,
    settings: CrossingSettings,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
    surface: MonotoneFidelitySurface | None = None,
    surfaces: list[MonotoneFidelitySurface] | None = None,
    axes=None,
    png_path: str | Path | None = None,
    show: bool = True,
    log_x: bool = False,
    label_runs: bool = True,
    run_label_limit: int = 300,
) -> list:
    """
    Plot the monotone p=0.5 line with fixed-ratio bootstrap uncertainty.

    The point monotone fit remains the central contour. Bootstrap uncertainty
    is shown as quantiles of the crossing gate count at fixed ratios, avoiding
    the displaced-center artifact caused by averaging whole fitted surfaces.
    """
    axes = plot_monotone_fidelity_surface_contours(
        data, settings, surface=surface, axes=axes, show=False, log_x=log_x,
        label_runs=label_runs, run_label_limit=run_label_limit)
    groups = sorted(grouped_by_n_qubits(data).items())

    for ax, _ in zip(axes, groups):
        group_surfaces = surfaces
        if group_surfaces is None:
            group_surfaces = bootstrap_surfaces(
                data,
                settings,
                n_bootstrap=n_bootstrap,
                seed=seed,
            )
        plot_bootstrap_gate_bands(
            ax,
            settings,
            group_surfaces,
            total_ratio_axes=True,
        )
        ax.legend(loc="best")

    if axes and png_path is not None:
        axes[0].figure.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()

    return axes


def plot_uncertainty_diagnostics(
    data: RMBData,
    settings: CrossingSettings,
    *,
    axes=None,
    png_path: str | Path | None = None,
    show: bool = True,
    log_x: bool = True,
) -> list:
    """
    Plot posterior mean and threshold ambiguity in total-gates/ratio space.

    Ambiguity is ``1 - 2 * abs(P(fidelity > 0.5) - 0.5)``: zero is confident
    either way, one is maximally undecided at the threshold.
    """
    import matplotlib.pyplot as plt

    groups = sorted(grouped_by_n_qubits(data).items())
    if not groups:
        return []

    if axes is None:
        _, axes_arr = plt.subplots(
            len(groups), 2, figsize=(11, 4.2 * len(groups)), squeeze=False,
            constrained_layout=True)
        axes_grid = axes_arr
    else:
        axes_list = list(axes)
        if len(axes_list) < 2 * len(groups):
            raise ValueError(
                f"Need at least {2 * len(groups)} axes, got {len(axes_list)}.")
        axes_grid = np.asarray(axes_list, dtype=object).reshape(len(groups), 2)

    gates_span = settings.n_gates_bounds[1] - settings.n_gates_bounds[0]
    ratio_span = settings.ratio_bounds[1] - settings.ratio_bounds[0]
    x_pad = 0.05 * gates_span
    y_pad = 0.05 * ratio_span

    for row, (n_qubits, group) in enumerate(groups):
        n_gates = np.array([config.n_gates for config in group], dtype=float)
        ratios = np.array([config.ratio_2_qb_gates for config in group], dtype=float)
        means = np.array([estimator.posterior_mean() for estimator in group.values()],
                         dtype=float)
        above = np.array([posterior_above(estimator) for estimator in group.values()],
                         dtype=float)
        ambiguity = 1.0 - 2.0 * np.abs(above - 0.5)
        ambiguity = np.clip(ambiguity, 0.0, 1.0)
        repeats = np.array([estimator.num_runs() for estimator in group.values()],
                           dtype=float)
        sizes = 35.0 + 8.0 * np.sqrt(np.maximum(repeats, 0.0))

        panels = [
            (axes_grid[row, 0], means, "Posterior mean fidelity", "Fidelity", "viridis"),
            (axes_grid[row, 1], ambiguity, "Threshold ambiguity", "Ambiguity", "magma_r"),
        ]
        for ax, values, title, colorbar_label, cmap in panels:
            scatter = ax.scatter(
                n_gates,
                ratios,
                c=values,
                cmap=cmap,
                vmin=0.0,
                vmax=1.0,
                s=sizes,
                edgecolors="black",
                linewidths=0.35,
                zorder=3,
            )
            plot_parametric_boundary_total_ratio_line(ax, group, settings)
            plot_analytic_total_ratio_line(ax, settings)
            plot_reference_total_ratio_line(ax, settings)
            if log_x:
                ax.set_xscale("log")
                x_lower = max(settings.n_gates_bounds[0] - x_pad,
                              0.8 * settings.n_gates_bounds[0])
            else:
                x_lower = settings.n_gates_bounds[0] - x_pad
            ax.set_xlim(x_lower, settings.n_gates_bounds[1] + x_pad)
            ax.set_ylim(
                max(0.0, settings.ratio_bounds[0] - y_pad),
                min(1.0, settings.ratio_bounds[1] + y_pad),
            )
            ax.set_xlabel("# Gates")
            ax.set_ylabel("Two-qudit gate ratio")
            ax.set_title(f"{title} (# Qubits = {n_qubits})")
            ax.grid(True, which="major", alpha=0.3, linewidth=0.6)
            ax.grid(True, which="minor", alpha=0.12, linewidth=0.4)
            ax.set_axisbelow(True)
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label(colorbar_label)
            ax.legend(loc="best")

    axes_out = list(axes_grid.ravel())
    if axes_out and png_path is not None:
        axes_out[0].figure.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    return axes_out


def plot_crossing_results(
    data: RMBData,
    settings: CrossingSettings,
    crossings: list[RMBConfig],
    *,
    base_path: Path | None = None,
    show: bool = True,
) -> None:
    """
    Show the standard result plots of a crossing experiment.

    The monotone-surface confidence plot, the monotone level-set plot, and
    the level-line plot are drawn for the same data, so the experiments in
    this package produce directly comparable figures. The monotone fit and
    its bootstrap surfaces are computed once and shared by all three plots;
    the level-line plot overlays the fit's p=0.5 contour. With ``base_path``
    (the resolved save path of the run, see :func:`~.common.save_crossings`)
    the figures are saved as ``*_surface.png``, ``*_levelset.png``,
    ``*_uncertainty.png``, and ``*.png`` siblings.
    """
    import matplotlib.pyplot as plt

    surface = try_fit_monotone_fidelity_surface(data, settings)
    surfaces = bootstrap_surfaces(data, settings, n_bootstrap=100, seed=settings.rng_seed)

    def sibling_path(suffix: str) -> Path | None:
        return None if base_path is None else base_path.parent / f"{base_path.stem}{suffix}"

    plot_monotone_fidelity_surface_with_confidence(
        data, settings, surface=surface, surfaces=surfaces,
        png_path=sibling_path("_surface.png"), show=False, log_x=True)
    plot_monotone_level_set(data, settings, surface=surface, surfaces=surfaces,
                            png_path=sibling_path("_levelset.png"), show=False)
    bins = settings.scatter_merge_bins or (None, None)
    plot_level_line(data, crossings,
                    contour=monotone_fit_contour(data, settings, surface=surface),
                    n_1qb_gates_bin=bins[0], n_2qb_gates_bin=bins[1],
                    png_path=sibling_path(".png"), show=False)
    if show:
        plt.show()


def monotone_fit_contour(
    data: RMBData,
    settings: CrossingSettings,
    surface: MonotoneFidelitySurface | None = None,
) -> np.ndarray | None:
    """p=0.5 contour points of the shared monotone fit, or ``None`` without a fit."""
    if surface is None:
        surface = try_fit_monotone_fidelity_surface(data, settings)
    if surface is None:
        return None
    return contour_points_from_surface(surface, settings)


def level_set_grid(one_q_bounds: tuple[int, int], two_q_bounds: tuple[int, int],
                   n: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """
    (one-qubit, two-qubit) meshgrid over the gate-count box, ``n`` per axis.

    Evaluate level-set predictions on this grid (flattened with ``ravel``)
    and pass them, reshaped back to the grid shape, to
    :func:`plot_gp_level_set` with the same bounds.
    """
    one_q_axis = np.linspace(one_q_bounds[0], one_q_bounds[1], n)
    two_q_axis = np.linspace(two_q_bounds[0], two_q_bounds[1], n)
    return np.meshgrid(one_q_axis, two_q_axis)


def gate_plane_points(one_q_grid: np.ndarray, two_q_grid: np.ndarray) -> np.ndarray:
    """(total gates, two-qubit ratio) coordinates of gate-count grid points."""
    n_gates = one_q_grid.ravel() + two_q_grid.ravel()
    return np.column_stack([n_gates, two_q_grid.ravel() / np.maximum(n_gates, 1.0)])


def plot_gp_level_set(
    probabilities: np.ndarray,
    latent_mean: np.ndarray,
    latent_variance: np.ndarray,
    results: list[tuple[float, float, int]],
    *,
    one_q_bounds: tuple[int, int],
    two_q_bounds: tuple[int, int],
    settings: CrossingSettings | None = None,
    target: float = 0.5,
    title: str = "GP level-set estimation",
    ax=None,
    png_path: str | Path | None = None,
    show: bool = True,
):
    """
    Plot a fitted probability surface with its target level set.

    The inputs are plain arrays of model predictions evaluated on the
    :func:`level_set_grid` of the box, in the grid's shape: the success
    probability and the latent (probit-space) mean and variance, from which
    the dashed and dotted lines bound the level set at 0.69 and 1.96 latent
    standard deviations. ``results`` holds (one-qubit gates, two-qubit
    gates, outcome) points drawn as success/failure markers.
    """
    import matplotlib.pyplot as plt
    from scipy.special import ndtr

    probabilities = np.asarray(probabilities, dtype=float)
    latent_mean = np.asarray(latent_mean, dtype=float).reshape(probabilities.shape)
    latent_std = np.sqrt(np.asarray(latent_variance, dtype=float)).reshape(probabilities.shape)
    one_q_grid, two_q_grid = level_set_grid(one_q_bounds, two_q_bounds,
                                            n=probabilities.shape[0])

    bands = {
        "dashed": (ndtr(latent_mean - 0.69 * latent_std),
                   ndtr(latent_mean + 0.69 * latent_std)),
        "dotted": (ndtr(latent_mean - 1.96 * latent_std),
                   ndtr(latent_mean + 1.96 * latent_std)),
    }

    own_figure = ax is None
    if own_figure:
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    else:
        fig = ax.figure
    mesh = ax.contourf(one_q_grid, two_q_grid, probabilities, levels=20, cmap="RdYlGn")
    ax.contour(one_q_grid, two_q_grid, probabilities, levels=[target],
               colors="k", linewidths=2)
    for linestyle, (low, high) in bands.items():
        ax.contour(one_q_grid, two_q_grid, low, levels=[target],
                   colors="k", linewidths=2, linestyles=linestyle)
        ax.contour(one_q_grid, two_q_grid, high, levels=[target],
                   colors="k", linewidths=2, linestyles=linestyle)
    plt.colorbar(mesh, ax=ax)

    failures = [(one_q, two_q) for one_q, two_q, outcome in results if outcome == 0]
    successes = [(one_q, two_q) for one_q, two_q, outcome in results if outcome == 1]
    if failures:
        ax.scatter(*zip(*failures), c="black", edgecolors="white", marker="x",
                   s=50, label="Failure", zorder=6, alpha=0.9)
    if successes:
        ax.scatter(*zip(*successes), c="white", edgecolors="black", marker="o",
                   s=50, label="Success", zorder=6, alpha=0.9)

    if settings is not None:
        plot_analytic_gate_plane_line(ax, settings, one_q_bounds, two_q_bounds)
    plot_reference_gate_plane_line(ax, one_q_bounds, two_q_bounds)

    # Empty handles give the level-set line styles legend entries.
    ax.plot([], [], color="k", linewidth=2, label="mean")
    ax.plot([], [], color="k", linewidth=2, linestyle="dashed", label=r"0.69$\sigma$")
    ax.plot([], [], color="k", linewidth=2, linestyle="dotted", label=r"1.96$\sigma$")
    ax.set_xlabel("# 1-qubit gates")
    ax.set_ylabel("# 2-qubit gates")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    if own_figure:
        fig.tight_layout()
    if png_path is not None:
        fig.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    return ax


def plot_monotone_level_set(
    data: RMBData,
    settings: CrossingSettings,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
    surface: MonotoneFidelitySurface | None = None,
    surfaces: list[MonotoneFidelitySurface] | None = None,
    target: float = 0.5,
    ax=None,
    png_path: str | Path | None = None,
    show: bool = True,
):
    """
    Level-set figure of the point monotone surface fit with bootstrap bands.

    The point monotone fit is the central surface and contour. Bootstrap
    uncertainty is shown as fixed-ratio contour quantile bands converted into
    the one-qubit/two-qubit gate plane. This avoids plotting a displaced
    bootstrap-averaged surface as though it were a better central estimate.
    """
    import matplotlib.pyplot as plt

    if surface is None:
        surface = try_fit_monotone_fidelity_surface(data, settings)
    if surface is None:
        return None
    if surfaces is None:
        surfaces = bootstrap_surfaces(data, settings, n_bootstrap=n_bootstrap, seed=seed)

    measured = measured_items(data)
    one_q_bounds, two_q_bounds = measured_gate_count_bounds(data)

    one_q_grid, two_q_grid = level_set_grid(one_q_bounds, two_q_bounds)
    points = gate_plane_points(one_q_grid, two_q_grid)
    probabilities = surface.probability(points).reshape(one_q_grid.shape)

    own_figure = ax is None
    if own_figure:
        fig, axis = plt.subplots(1, 1, figsize=(6, 5))
    else:
        axis = ax
        fig = axis.figure

    mesh = axis.contourf(one_q_grid, two_q_grid, probabilities, levels=20, cmap="RdYlGn")
    plt.colorbar(mesh, ax=axis).set_label("Point monotone fit fidelity")
    if float(np.min(probabilities)) <= target <= float(np.max(probabilities)):
        axis.contour(one_q_grid, two_q_grid, probabilities, levels=[target],
                     colors="black", linewidths=2)
        axis.plot([], [], color="black", linewidth=2, label="point monotone fit")

    if surfaces:
        plot_bootstrap_gate_bands(
            axis,
            settings,
            surfaces,
            level=target,
            total_ratio_axes=False,
        )

    failures = [
        (float(config.n_1qb_gates), float(config.n_2qb_gates))
        for config, estimator in measured
        if estimator.posterior_mean() < target
    ]
    successes = [
        (float(config.n_1qb_gates), float(config.n_2qb_gates))
        for config, estimator in measured
        if estimator.posterior_mean() >= target
    ]
    if failures:
        axis.scatter(*zip(*failures), c="black", marker="x",
                     s=50, label="posterior mean < target", zorder=6, alpha=0.9)
    if successes:
        axis.scatter(*zip(*successes), c="white", edgecolors="black", marker="o",
                     s=50, label="posterior mean >= target", zorder=6, alpha=0.9)

    plot_analytic_gate_plane_line(axis, settings, one_q_bounds, two_q_bounds)
    plot_reference_gate_plane_line(axis, one_q_bounds, two_q_bounds)

    plot_parametric_boundary_gate_plane_line(axis, data, settings,
                                             one_q_bounds, two_q_bounds)
    axis.set_xlabel("# 1-qubit gates")
    axis.set_ylabel("# 2-qubit gates")
    axis.set_title("Monotone-fit level set with contour bootstrap bands")
    axis.legend(loc="upper right", fontsize=8)
    if own_figure:
        fig.tight_layout()

    if png_path is not None:
        axis.figure.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        import matplotlib.pyplot as plt
        plt.show()
    return axis
