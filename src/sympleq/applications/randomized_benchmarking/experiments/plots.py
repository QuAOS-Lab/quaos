"""
Plot functions shared by the boundary experiments in this package.

All plots live in (total gates, two-qubit gate ratio) space, one subplot per
``n_qubits``, so the different experiments can be compared directly.
"""
from __future__ import annotations

from pathlib import Path
from dataclasses import replace
import numpy as np

from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    CrossingSettings,
    MonotoneFidelitySurface,
    bootstrap_contours,
    bootstrap_surfaces,
    candidate_axes,
    contour_points_from_surface,
    grouped_by_n_qubits,
    measured_gate_count_bounds,
    measured_items,
    try_fit_monotone_fidelity_surface,
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
              level_line: list[RMBConfig] | None = None, log_x: bool = False) -> list:
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

    Returns
    -------
    list[matplotlib.axes.Axes]
        The axes the plots were drawn on (sorted by ``n_qubits``).
    """
    import matplotlib.pyplot as plt
    from matplotlib import patheffects
    from matplotlib.colors import LinearSegmentedColormap

    data = merge_close_configs(data, n_1qb_gates_bin=20, n_2qb_gates_bin=10)
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
    axes=None,
    png_path: str | Path | None = None,
    show: bool = True,
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
    png_path : str | Path | None
        Where to save the figure. ``None`` skips saving.
    show : bool
        If ``True``, call ``plt.show()`` at the end.
    """
    import matplotlib.pyplot as plt

    axes = plot_data(data, axes=axes, show=False, skip_incomplete=False,
                     level_line=crossings, log_x=True)
    if axes and contour is not None and len(contour) > 0:
        order = np.argsort(contour[:, 1])
        axes[0].plot(contour[order, 0], contour[order, 1], linestyle="-.",
                     color="0.35", linewidth=1.6, label="monotone fit p=0.5", zorder=4)
        axes[0].legend(loc="best", frameon=True, framealpha=0.9)
    if axes and png_path is not None:
        axes[0].figure.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    return axes


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
    Plot the monotone p=0.5 line with monotone bootstrap uncertainty.

    The uncertainty is shown as the occupancy of bootstrap contours over the
    candidate grid cells. ``surface`` and ``surfaces`` take precomputed fits
    (see :func:`plot_monotone_fidelity_surface_contours` and
    :func:`~.common.bootstrap_contours`), so one fit and one bootstrap can be
    shared between the sibling plots of a run.
    """
    import matplotlib.pyplot as plt

    axes = plot_monotone_fidelity_surface_contours(
        data, settings, surface=surface, axes=axes, show=False, log_x=log_x,
        label_runs=label_runs, run_label_limit=run_label_limit)
    groups = sorted(grouped_by_n_qubits(data).items())

    for ax, (_, group) in zip(axes, groups):
        contours = bootstrap_contours(
            group,
            settings,
            n_bootstrap=n_bootstrap,
            seed=seed,
            surfaces=surfaces,
        )
        if not contours:
            continue

        gates_bins, ratio_bins = candidate_axes(settings)
        occupancy = np.zeros((len(ratio_bins) - 1, len(gates_bins) - 1), dtype=float)

        for contour in contours:
            gates_idx = np.searchsorted(gates_bins, contour[:, 0], side="right") - 1
            ratio_idx = np.searchsorted(ratio_bins, contour[:, 1], side="right") - 1
            valid = (
                (0 <= gates_idx)
                & (gates_idx < occupancy.shape[1])
                & (0 <= ratio_idx)
                & (ratio_idx < occupancy.shape[0])
            )
            cells = set(zip(ratio_idx[valid], gates_idx[valid]))
            for row, col in cells:
                occupancy[row, col] += 1.0

        occupancy /= len(contours)
        mesh = ax.pcolormesh(
            gates_bins,
            ratio_bins,
            occupancy,
            cmap="Blues",
            vmin=0.0,
            vmax=max(0.05, float(np.max(occupancy))),
            shading="auto",
            alpha=0.55,
            zorder=1,
        )
        cbar = plt.colorbar(mesh, ax=ax)
        cbar.set_label("Monotone bootstrap contour occupancy")
        ax.plot([], [], color="tab:blue", alpha=0.6, linewidth=6,
                label=f"monotone contour uncertainty ({len(contours)}/{n_bootstrap})")
        ax.legend(loc="best")

    if axes and png_path is not None:
        axes[0].figure.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()

    return axes


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
    the figures are saved as ``*_surface.png``, ``*_levelset.png``, and
    ``*.png`` siblings.
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
    plot_level_line(data, crossings,
                    contour=monotone_fit_contour(data, settings, surface=surface),
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
    Level-set figure of the monotone surface fit with bootstrap bands.

    The monotone fit and its bootstrap surfaces (precomputable via
    ``surface``/``surfaces``) are evaluated on the (one-qubit, two-qubit)
    gate-count box spanned by the measured configs; the latent statistics
    for the credible bands are the probit-transformed bootstrap
    probabilities, so the figure matches the GP level-set plot of the
    aepsych experiment. The plotted surface and solid mean line are the
    bootstrap band center, so the line sits between its credible bands by
    construction; the point fit's contour is overlaid separately, making
    any bootstrap location bias visible as the gap between the two lines.
    Each measured config is drawn as a success or failure marker by the
    side of ``target`` its posterior mean lies on. Returns ``None`` when
    the data cannot support a fit.
    """
    from scipy.special import ndtr, ndtri

    if surface is None:
        surface = try_fit_monotone_fidelity_surface(data, settings)
    if surface is None:
        return None
    if surfaces is None:
        surfaces = bootstrap_surfaces(data, settings, n_bootstrap=n_bootstrap, seed=seed)
    if not surfaces:
        return None

    measured = measured_items(data)
    one_q_bounds, two_q_bounds = measured_gate_count_bounds(data)

    one_q_grid, two_q_grid = level_set_grid(one_q_bounds, two_q_bounds)
    points = gate_plane_points(one_q_grid, two_q_grid)

    eps = 1e-9
    latents = np.stack([
        ndtri(np.clip(boot.probability(points), eps, 1.0 - eps))
        for boot in surfaces
    ])
    latent_mean = latents.mean(axis=0).reshape(one_q_grid.shape)
    latent_variance = latents.var(axis=0).reshape(one_q_grid.shape)
    band_center = ndtr(latent_mean)

    results = [
        (float(config.n_1qb_gates), float(config.n_2qb_gates),
         1 if estimator.posterior_mean() >= target else 0)
        for config, estimator in measured
    ]
    axis = plot_gp_level_set(
        band_center, latent_mean, latent_variance, results,
        one_q_bounds=one_q_bounds, two_q_bounds=two_q_bounds, target=target,
        title="Monotone-fit level set (bootstrap bands)",
        ax=ax, png_path=None, show=False)

    point_fit = surface.probability(points).reshape(one_q_grid.shape)
    if float(np.min(point_fit)) <= target <= float(np.max(point_fit)):
        axis.contour(one_q_grid, two_q_grid, point_fit, levels=[target],
                     colors="tab:blue", linewidths=1.6)
        axis.plot([], [], color="tab:blue", linewidth=1.6, label="point fit")
        axis.legend(loc="upper right", fontsize=8)

    if png_path is not None:
        axis.figure.savefig(png_path, dpi=200, bbox_inches="tight")
    if show:
        import matplotlib.pyplot as plt
        plt.show()
    return axis
