"""Plot and score one saved FLE GP grid."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtri

from sympleq.applications.randomized_benchmarking.experiments.scores import (
    gp_grid_scores,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.fantasy_levelset_settings import (
    HQC_BUDGET,
)


# -------------------------------------------------------------------------
# PATH / SCORE HANDLES
# -------------------------------------------------------------------------

RMB_JSON_PATH = Path(
    r"Personal\seed_2025\rick_fantasy_gpu_crossing_20260623_184228.json"
)
GP_GRID_PATH: Path | None = None
PNG_PATH: Path | None = None

ONE_Q_NOISE_SCALE = 1.0
TWO_Q_NOISE_SCALE = 1.0
LAST_BACKEND_BATCH_SIZE: int | None = None


ANALYTIC_LINE_LABEL = "analytic Lindblad line"
REFERENCE_CURVES = [
    {
        "label": "reference line nq=5",
        "numerator": 0.7106,
        "offset": 1.91e-4,
        "slope": 3.65e-3,
        "color": "tab:purple",
        "linestyle": "--",
    },
    {
        "label": "reference line nq=20",
        "numerator": 0.7106,
        "offset": 1.91e-4,
        "slope": 3.65e-3,
        "color": "tab:blue",
        "linestyle": "-.",
    },
]
REFERENCE_OFFSET = REFERENCE_CURVES[0]["offset"]
REFERENCE_SLOPE = REFERENCE_CURVES[0]["slope"]


def reference_gate_counts(
    ratios: np.ndarray,
    *,
    numerator: float,
    offset: float,
    slope: float,
) -> np.ndarray:
    return numerator / (offset + slope * ratios)


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


def sibling_grid_path(json_path: Path) -> Path:
    return json_path.parent / f"{json_path.stem}_gp_grid.npz"


def sibling_png_path(json_path: Path) -> Path:
    return json_path.parent / f"{json_path.stem}_fle_levelset_score.png"


def load_points(
    json_path: Path,
    *,
    last_backend_batch_size: int | None = LAST_BACKEND_BATCH_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    gates: list[float] = []
    ratios: list[float] = []
    outcomes: list[int] = []
    records = payload.get("data", [])
    if last_backend_batch_size is not None:
        records = records[-last_backend_batch_size:]

    for record in records:
        n_1q = int(record["n_1qb_gates"])
        n_2q = int(record["n_2qb_gates"])
        total = n_1q + n_2q
        if total <= 0:
            continue
        counts = {int(outcome): int(count) for outcome, count in record["results"]}
        successes = counts.get(1, 0)
        failures = counts.get(0, 0)
        if successes + failures <= 0:
            continue
        gates.append(float(total))
        ratios.append(float(n_2q / total))
        outcomes.append(int(successes >= failures))

    return (
        np.asarray(gates, dtype=float),
        np.asarray(ratios, dtype=float),
        np.asarray(outcomes, dtype=bool),
    )


def plot_fle_grid(
    json_path: str | Path = RMB_JSON_PATH,
    grid_path: str | Path | None = GP_GRID_PATH,
    png_path: str | Path | None = PNG_PATH,
) -> dict[str, float]:
    json_path = Path(json_path)
    grid_path = sibling_grid_path(json_path) if grid_path is None else Path(grid_path)
    png_path = sibling_png_path(json_path) if png_path is None else Path(png_path)

    grid = np.load(grid_path)
    gates_grid = np.asarray(grid["gates_grid"], dtype=float)
    ratio_grid = np.asarray(grid["ratio_grid"], dtype=float)
    probabilities = np.asarray(grid["probabilities"], dtype=float)
    latent_mean = np.asarray(grid["latent_mean"], dtype=float)
    latent_variance = np.maximum(np.asarray(grid["latent_variance"], dtype=float), 0.0)
    latent_std = np.sqrt(latent_variance)
    target = float(np.asarray(grid["target"]).item())
    latent_target = float(ndtri(target))

    scores = gp_grid_scores(
        grid_path,
        one_q_noise_scale=ONE_Q_NOISE_SCALE,
        two_q_noise_scale=TWO_Q_NOISE_SCALE,
    )

    point_gates, point_ratios, point_outcomes = load_points(json_path)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.6))
    surface = ax.contourf(
        gates_grid,
        ratio_grid,
        probabilities,
        levels=np.linspace(0.05, 0.95, 19),
        cmap="RdYlGn",
        alpha=0.85,
    )
    fig.colorbar(surface, ax=ax, label="Predicted fidelity")

    ax.contour(
        gates_grid,
        ratio_grid,
        latent_mean,
        levels=[latent_target],
        colors="black",
        linewidths=2.2,
        zorder=6,
    )
    ax.contour(
        gates_grid,
        ratio_grid,
        latent_mean + latent_std,
        levels=[latent_target],
        colors="black",
        linestyles="--",
        linewidths=1.6,
        zorder=6,
    )
    ax.contour(
        gates_grid,
        ratio_grid,
        latent_mean - latent_std,
        levels=[latent_target],
        colors="black",
        linestyles="--",
        linewidths=1.6,
        zorder=6,
    )

    if len(point_gates) > 0:
        ax.scatter(
            point_gates[~point_outcomes],
            point_ratios[~point_outcomes],
            marker="x",
            s=36,
            color="black",
            linewidths=1.5,
            label="Failure",
            zorder=8,
        )
        ax.scatter(
            point_gates[point_outcomes],
            point_ratios[point_outcomes],
            marker="o",
            s=42,
            facecolors="white",
            edgecolors="black",
            linewidths=1.2,
            label="Success",
            zorder=9,
        )

    ratios = np.linspace(float(np.min(ratio_grid)), float(np.max(ratio_grid)), 400)
    analytic = analytic_gate_counts(
        ratios,
        one_q_noise_scale=ONE_Q_NOISE_SCALE,
        two_q_noise_scale=TWO_Q_NOISE_SCALE,
    )
    x_min = float(np.min(gates_grid[gates_grid > 0.0]))
    x_max = float(np.max(gates_grid))

    analytic_mask = np.isfinite(analytic) & (x_min <= analytic) & (analytic <= x_max)
    ax.plot(
        analytic[analytic_mask],
        ratios[analytic_mask],
        color="tab:green",
        linestyle=":",
        linewidth=2.0,
        label=ANALYTIC_LINE_LABEL,
        zorder=7,
    )
    for curve in REFERENCE_CURVES:
        reference = reference_gate_counts(
            ratios,
            numerator=float(curve["numerator"]),
            offset=float(curve["offset"]),
            slope=float(curve["slope"]),
        )
        reference_mask = np.isfinite(reference) & (x_min <= reference) & (reference <= x_max)
        ax.plot(
            reference[reference_mask],
            ratios[reference_mask],
            color=str(curve["color"]),
            linestyle=str(curve["linestyle"]),
            linewidth=1.8,
            label=str(curve["label"]),
            zorder=7,
        )

    ax.plot([], [], color="black", linewidth=2.2, label=f"GP mean p={target:g}")
    ax.plot([], [], color="black", linestyle="--", linewidth=1.6, label=r"$\mu \pm 1\sigma$")
    ax.set_xscale("log")
    ax.set_xlim(left=max(1.0, x_min), right=x_max)
    ax.set_ylim(0.07, float(np.max(ratio_grid)))
    ax.set_xlabel("Total gates")
    ax.set_ylabel("Two-qubit gate ratio")
    ax.set_title(
        f"FLE GP level set | HQC={HQC_BUDGET:g} | S1={scores['S1']:.4g}, "
        f"S2={scores['S2']:.4g}, A_gp={scores['A_gp']:.4g}"
    )
    ax.legend(loc="best", fontsize=8, frameon=True, framealpha=0.9)
    fig.tight_layout()

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    print(f"[score] S1={scores['S1']:.6g}")
    print(f"[score] S2={scores['S2']:.6g}")
    print(f"[score] A_gp={scores['A_gp']:.6g}")
    print(f"[saved] plot: {png_path}")
    plt.show()
    return scores


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else RMB_JSON_PATH
    plot_fle_grid(path)
