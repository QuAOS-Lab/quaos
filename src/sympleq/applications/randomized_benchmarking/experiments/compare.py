"""
Run all three crossing experiments on the same problem and compare them in
one figure.

Each method runs with its own settings (sharing the problem definition,
budget, cost model, and the common-random-numbers seed, so any config
measured by several methods records identical outcomes), and the comparison
figure holds one row per method with the three standard plots as columns:
the monotone-surface confidence plot, the monotone level-set plot, and the
level-line plot. :func:`replot_from_data` rebuilds the same figure from the
saved run files without re-running the simulations.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from sympleq.applications.randomized_benchmarking.RMB import RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.experiments import (
    baseline_crossing,
    charlie_crossing,
    level_crossing,
    rick_crossing,
)
from sympleq.applications.randomized_benchmarking.experiments.baseline_crossing import (
    BaselineCrossingSettings,
)
from sympleq.applications.randomized_benchmarking.experiments.charlie_crossing import (
    CharlieCrossingSettings,
)
from sympleq.applications.randomized_benchmarking.experiments.common import (
    bootstrap_surfaces,
    load_crossings,
    try_fit_monotone_fidelity_surface,
)
from sympleq.applications.randomized_benchmarking.experiments.level_crossing import (
    LevelCrossingSettings,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    monotone_fit_contour,
    plot_level_line,
    plot_monotone_fidelity_surface_with_confidence,
    plot_monotone_level_set,
)
from sympleq.applications.randomized_benchmarking.experiments.rick_crossing import (
    RickCrossingSettings,
)

COLUMN_TITLES = (
    "Data + monotone fit (bootstrap occupancy)",
    "Monotone level set (gate-count plane)",
    "Data + crossing line",
)

# (result key, row label, experiment module, settings class) of every
# compared method; both entry points below iterate this table. The baseline
# leads as the high-confidence reference the other methods are read against.
METHODS = (
    ("baseline", "Baseline (iterative decay fit)", baseline_crossing, BaselineCrossingSettings),
    # ("level", "Vanilla level crossing", level_crossing, LevelCrossingSettings),
    # ("charlie", "Charlie contour trace", charlie_crossing, CharlieCrossingSettings),
    # ("rick", "Rick GP level set", rick_crossing, RickCrossingSettings),
)


def comparison_figure(rows, *, png_path: str | Path | None = None, show: bool = True):
    """
    Draw the comparison figure: one method per row, the three plots as columns.

    Each row is a ``(name, data, settings, crossings)`` tuple. Every row's
    crossing-line plot overlays the monotone-fit p=0.5 contour, the fit
    shared by all methods. Column titles appear on the top row only; the
    method name labels each row on the left.

    Returns
    -------
    tuple
        The figure and its axes array.
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(rows), 3, figsize=(19, 4.4 * len(rows)),
                             constrained_layout=True, squeeze=False)
    for row, (name, data, settings, crossings) in enumerate(rows):
        surface = try_fit_monotone_fidelity_surface(data, settings)
        surfaces = bootstrap_surfaces(data, settings, seed=settings.rng_seed)
        plot_monotone_fidelity_surface_with_confidence(
            data, settings, surface=surface, surfaces=surfaces,
            axes=[axes[row][0]], show=False, log_x=True)
        plot_monotone_level_set(
            data, settings, surface=surface, surfaces=surfaces,
            ax=axes[row][1], show=False)
        bins = settings.scatter_merge_bins or (None, None)
        plot_level_line(data, crossings,
                        contour=monotone_fit_contour(data, settings, surface=surface),
                        axes=[axes[row][2]], show=False,
                        n_1qb_gates_bin=bins[0], n_2qb_gates_bin=bins[1])

        axes[row][0].annotate(
            name, xy=(0, 0.5), xycoords="axes fraction", xytext=(-0.42, 0.5),
            rotation=90, ha="center", va="center", fontsize=13, fontweight="bold")
        for column in range(3):
            axes[row][column].set_title(COLUMN_TITLES[column] if row == 0 else "")

    if png_path is not None:
        fig.savefig(resolve_data_path(png_path), dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    return fig, axes


def run_and_compare(
    level_settings: LevelCrossingSettings | None = None,
    charlie_settings: CharlieCrossingSettings | None = None,
    rick_settings: RickCrossingSettings | None = None,
    baseline_settings: BaselineCrossingSettings | None = None,
    *,
    png_path: str | Path | None = "crossing_comparison.png",
    show: bool = True,
) -> dict:
    """
    Run every crossing experiment and draw the comparison figure.

    Each experiment runs with its individual plotting disabled; everything
    else (verbose output, data saving) follows its settings. Bare ``png_path``
    file names are resolved inside the package's ``rmb_data`` directory. The
    baseline ignores the HQC budget, so it is the slowest row to run.

    Returns
    -------
    dict
        ``{"baseline": (rmb, crossings), "level": (rmb, crossings),
        "charlie": (rmb, crossings), "rick": (rmb, crossings)}``.
    """
    overrides = {"level": level_settings, "charlie": charlie_settings,
                 "rick": rick_settings, "baseline": baseline_settings}

    results = {}
    rows = []
    for key, name, module, settings_class in METHODS:
        settings = replace(overrides[key] or settings_class(), plot=False)
        rmb, crossings = module.run(settings)
        results[key] = (rmb, crossings)
        rows.append((name, rmb._data, settings, crossings))

    comparison_figure(rows, png_path=png_path, show=show)
    return results


def replot_from_data(
    level_path: str | Path | None = None,
    charlie_path: str | Path | None = None,
    rick_path: str | Path | None = None,
    baseline_path: str | Path | None = None,
    *,
    level_settings: LevelCrossingSettings | None = None,
    charlie_settings: CharlieCrossingSettings | None = None,
    rick_settings: RickCrossingSettings | None = None,
    baseline_settings: BaselineCrossingSettings | None = None,
    png_path: str | Path | None = "crossing_comparison.png",
    show: bool = True,
):
    """
    Rebuild the comparison figure from saved runs, without re-running them.

    The paths are the ``save_path`` values of previous runs (bare names
    resolve inside the package's ``rmb_data`` directory); ``None`` uses the
    default ``save_path`` of the method's settings. The data and the
    crossings are loaded from the files written by those runs. The settings
    only parameterize the fits and plots.

    Returns
    -------
    tuple
        The figure and its axes array.
    """
    overrides = {"level": (level_path, level_settings),
                 "charlie": (charlie_path, charlie_settings),
                 "rick": (rick_path, rick_settings),
                 "baseline": (baseline_path, baseline_settings)}

    rows = []
    for key, name, _, settings_class in METHODS:
        path, settings = overrides[key]
        settings = settings or settings_class()
        path = path if path is not None else settings.save_path
        rows.append((name, RMB.load(path)._data, settings, load_crossings(path)))
    return comparison_figure(rows, png_path=png_path, show=show)


if __name__ == "__main__":
    replot_from_data()
