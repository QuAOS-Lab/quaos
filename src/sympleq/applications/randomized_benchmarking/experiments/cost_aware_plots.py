"""Plot helpers specific to the cost-aware surface experiment."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from sympleq.applications.randomized_benchmarking.experiments.cost_aware_references import (
    analytic_lindblad_gates,
    gp_grid_surface_gates,
    gp_grid_surface_label,
    reference_surface_volume,
)


def plot_boundary_surface(
    ratios: np.ndarray,
    qubits: np.ndarray,
    log10_gates: np.ndarray,
    measured: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    settings,
    png_path: str | Path | None = None,
    *,
    show: bool = False,
    show_block: bool = True,
    show_pause: float = 0.001,
    close: bool = True,
    figure_name: str | None = None,
):
    """Plot the fitted fidelity-0.5 boundary surface."""
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

    fig = plt.figure(num=figure_name, figsize=(9.5, 7.5), clear=True)
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(
        ratios,
        qubits,
        log10_gates,
        cmap="viridis",
        alpha=0.48,
        linewidth=0,
        antialiased=True,
    )

    n_bounds = getattr(settings, "surface_plot_n_gates_bounds", None)
    if n_bounds is None:
        n_bounds = settings.n_gates_bounds
    n_lo, n_hi = n_bounds

    if getattr(settings, "surface_plot_analytic", True):
        analytic = analytic_lindblad_gates(ratios, qubits, settings)
        analytic_z = np.where(
            (analytic >= n_lo) & (analytic <= n_hi),
            np.log10(analytic),
            np.nan,
        )
        ax.plot_wireframe(
            ratios,
            qubits,
            analytic_z,
            color="crimson",
            linewidth=1.1,
            alpha=0.75,
            rstride=4,
            cstride=4,
            label="analytic Lindblad",
        )

    ratios_m, gates_m, qubits_m, succ_m, fail_m = measured
    total = succ_m + fail_m
    mask = total > 0
    if np.any(mask):
        p_hat = succ_m[mask] / total[mask]
        size = 18.0 + 40.0 * np.clip(total[mask] / max(total[mask].max(), 1.0), 0, 1)
        scatter = ax.scatter(
            ratios_m[mask],
            qubits_m[mask],
            np.log10(np.maximum(gates_m[mask], 1e-12)),
            c=p_hat,
            cmap="coolwarm_r",
            norm=TwoSlopeNorm(vcenter=0.5, vmin=0.0, vmax=1.0),
            s=size,
            edgecolor="k",
            linewidth=0.3,
            depthshade=True,
        )
        cbar = fig.colorbar(scatter, ax=ax, shrink=0.6, pad=0.10)
        cbar.set_label("measured survival fraction")

    ax.set_xlabel("two-qubit ratio  r")
    ax.set_ylabel("n_qubits  Q")
    ax.set_zlabel(r"$\log_{10}(\mathrm{gate\ count}\ n)$")
    ax.set_zlim(np.log10(max(n_lo, 1e-12)), np.log10(max(n_hi, 1e-12)))
    ax.set_title(r"Fidelity-0.5 boundary surface  $\log_{10} n_*(r, Q)$")
    if getattr(settings, "surface_plot_analytic", True):
        ax.legend(loc="upper left")

    fig.tight_layout()
    if png_path is not None:
        png_path = Path(png_path)
        png_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show(block=show_block)
        if not show_block:
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(max(float(show_pause), 0.001))
    if close:
        plt.close(fig)
    return png_path


def plot_live_volume_history(
    history: list[tuple],
    settings,
    png_path: str | Path | None = None,
    *,
    show: bool = False,
    show_block: bool = True,
    show_pause: float = 0.001,
    close: bool = True,
    figure_name: str | None = None,
):
    """Plot fitted log10-gate surface volume history with the reference line."""
    if not history:
        return None
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    iterations = np.asarray([row[0] for row in history], dtype=int)
    fitted = np.asarray([row[1] for row in history], dtype=float)
    reference = np.asarray([row[2] for row in history], dtype=float)
    lower = np.asarray(
        [row[3] if len(row) > 3 else np.nan for row in history],
        dtype=float,
    )
    upper = np.asarray(
        [row[4] if len(row) > 4 else np.nan for row in history],
        dtype=float,
    )
    ref = float(reference[-1])
    ratio = fitted / ref if ref else np.full_like(fitted, np.nan)

    fig = plt.figure(num=figure_name, figsize=(8.5, 5.2), clear=True)
    ax = fig.add_subplot(111)
    has_band = np.isfinite(lower) & np.isfinite(upper)
    if np.any(has_band):
        yerr = np.vstack((
            np.where(has_band, np.maximum(fitted - lower, 0.0), np.nan),
            np.where(has_band, np.maximum(upper - fitted, 0.0), np.nan),
        ))
        ax.errorbar(
            iterations,
            fitted,
            yerr=yerr,
            marker="o",
            linewidth=1.8,
            elinewidth=1.1,
            capsize=3,
            label=r"posterior fit $\pm 1\sigma$ volume",
        )
    else:
        ax.plot(iterations, fitted, marker="o", linewidth=1.8, label="posterior fit")
    ax.axhline(
        ref,
        color="crimson",
        linestyle="--",
        linewidth=1.6,
        label="known-rate exponential reference",
    )
    if getattr(settings, "gp_grid_surface_path", None) is not None:
        gp_ref = reference_surface_volume(settings, gp_grid_surface_gates)
        if np.isfinite(gp_ref):
            ax.axhline(
                gp_ref,
                color="black",
                linestyle="--",
                linewidth=1.6,
                label=gp_grid_surface_label(settings),
            )
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\int \log_{10} n_*(r,Q)\,dr\,dQ$")
    ax.set_title(r"Fidelity-0.5 surface volume in $\log_{10}(n)$")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax_ratio = ax.twinx()
    ax_ratio.plot(
        iterations,
        ratio,
        color="0.25",
        linestyle=":",
        marker=".",
        alpha=0.8,
        label="fit/reference",
    )
    ax_ratio.set_ylabel("fit / reference")
    ax_ratio.axhline(1.0, color="0.25", linestyle=":", linewidth=1.0, alpha=0.6)
    if ref:
        left_lo, left_hi = ax.get_ylim()
        ax_ratio.set_ylim(left_lo / ref, left_hi / ref)

    fig.tight_layout()
    if png_path is not None:
        png_path = Path(png_path)
        png_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show(block=show_block)
        if not show_block:
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(max(float(show_pause), 0.001))
    if close:
        plt.close(fig)
    return png_path
