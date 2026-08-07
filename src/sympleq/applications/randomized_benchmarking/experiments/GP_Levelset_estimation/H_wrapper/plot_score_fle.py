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
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.fantasy_levelset_settings import (
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

GATES_AXIS_LIMITS: tuple[float, float] | None = (290.0, 2010.0)
RATIO_AXIS_LIMITS: tuple[float, float] | None = (0.4999, 0.9001)

SHOW_PREDICTED_FIDELITY_HUE = False
SHOW_SOBOL_POINTS = True
SOBOL_SUCCESS_COLOR = "tab:blue"
SOBOL_FAILURE_COLOR = "tab:orange"

SHOW_Q56_REFERENCE_LINE = False
SHOW_Q56_REFERENCE_POINTS = False
Q56_REFERENCE_JSON_PATH: Path | None = Path(
    r"Personal\FLE\H2_2\qband_4\accumulation"
    r"\accumulated_final_H2_2_qband_4_20260804_090823_actual_gates_q_slices_gap_3"
    r"\accumulated_final_H2_2_qband_4_20260804_090823_actual_gates_q56.json"
)
Q56_REFERENCE_GRID_PATH: Path | None = Path(
    r"Personal\FLE\H2_2\qband_4\accumulation"
    r"\accumulated_final_H2_2_qband_4_20260804_090823_actual_gates_q_slices_gap_3"
    r"\accumulated_final_H2_2_qband_4_20260804_090823_actual_gates_q56_gp_grid.npz"
)
Q56_REFERENCE_LABEL = "accumulated q=56"
Q56_REFERENCE_COLOR = "tab:blue"
Q56_REFERENCE_LINESTYLE = "-."
Q56_REFERENCE_POINTS_LABEL = "accumulated q=56"
Q56_REFERENCE_SUCCESS_COLOR = "tab:green"
Q56_REFERENCE_FAILURE_COLOR = "tab:red"


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


def record_config_key(record: dict) -> tuple[int, int, int, float, bool]:
    return (
        int(record["n_1qb_gates"]),
        int(record["n_2qb_gates"]),
        int(record["n_qubits"]),
        float(record.get("random_elimination", 0.0)),
        bool(record.get("use_scrambler", True)),
    )


def sobol_config_keys(json_path: Path) -> set[tuple[int, int, int, float, bool]]:
    keys = set()
    for path in sorted(json_path.parent.glob("measurement_*_sobol*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        for record in payload.get("data", []):
            keys.add(record_config_key(record))
    return keys


def load_points(
    json_path: Path,
    *,
    last_backend_batch_size: int | None = LAST_BACKEND_BATCH_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    sobol_keys = sobol_config_keys(json_path)
    gates: list[float] = []
    ratios: list[float] = []
    outcomes: list[int] = []
    is_sobol: list[bool] = []
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
        is_sobol.append(record_config_key(record) in sobol_keys)

    return (
        np.asarray(gates, dtype=float),
        np.asarray(ratios, dtype=float),
        np.asarray(outcomes, dtype=bool),
        np.asarray(is_sobol, dtype=bool),
    )


def qubit_title_part(json_path: Path) -> str:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    experiment_q = payload.get("experiment", {}).get("n_qubits")
    if experiment_q is not None:
        return f" | q={int(experiment_q)}"

    qubits = sorted({
        int(record["n_qubits"])
        for record in payload.get("data", [])
        if "n_qubits" in record
    })
    if len(qubits) == 1:
        return f" | q={qubits[0]}"
    return ""


def add_reference_contour(
    ax,
    *,
    grid_path: Path | None,
    label: str,
    color: str,
    linestyle: str,
) -> None:
    if grid_path is None:
        return
    if not grid_path.exists():
        print(f"[warning] skipped reference line; missing grid: {grid_path}")
        return

    grid = np.load(grid_path)
    target = float(np.asarray(grid["target"]).item())
    latent_target = float(ndtri(target))
    ax.contour(
        np.asarray(grid["ratio_grid"], dtype=float),
        np.asarray(grid["gates_grid"], dtype=float),
        np.asarray(grid["latent_mean"], dtype=float),
        levels=[latent_target],
        colors=color,
        linestyles=linestyle,
        linewidths=2.0,
        zorder=7,
    )
    ax.plot([], [], color=color, linestyle=linestyle, linewidth=2.0, label=label)


def add_reference_points(
    ax,
    *,
    json_path: Path | None,
    label: str,
) -> None:
    if json_path is None:
        return
    if not json_path.exists():
        print(f"[warning] skipped reference data; missing JSON: {json_path}")
        return

    point_gates, point_ratios, point_outcomes, _ = load_points(json_path)
    if len(point_gates) == 0:
        print(f"[warning] skipped reference data; no points in: {json_path}")
        return

    failure_mask = ~point_outcomes
    success_mask = point_outcomes
    if np.any(failure_mask):
        ax.scatter(
            point_ratios[failure_mask],
            point_gates[failure_mask],
            marker="s",
            s=24,
            color=Q56_REFERENCE_FAILURE_COLOR,
            alpha=0.75,
            linewidths=0.0,
            label=f"{label} failure",
            zorder=7,
        )
    if np.any(success_mask):
        ax.scatter(
            point_ratios[success_mask],
            point_gates[success_mask],
            marker="s",
            s=24,
            color=Q56_REFERENCE_SUCCESS_COLOR,
            alpha=0.75,
            linewidths=0.0,
            label=f"{label} success",
            zorder=7,
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

    point_gates, point_ratios, point_outcomes, point_is_sobol = load_points(json_path)

    fig, ax = plt.subplots(1, 1, figsize=(8.0, 5.6))
    if SHOW_PREDICTED_FIDELITY_HUE:
        surface = ax.contourf(
            ratio_grid,
            gates_grid,
            probabilities,
            levels=np.linspace(0.0, 1.0, 21),
            cmap="RdYlGn",
            alpha=0.85,
        )
        fig.colorbar(surface, ax=ax, label="Predicted fidelity")

    ax.contour(
        ratio_grid,
        gates_grid,
        latent_mean,
        levels=[latent_target],
        colors="black",
        linewidths=2.2,
        zorder=6,
    )
    ax.contour(
        ratio_grid,
        gates_grid,
        latent_mean + latent_std,
        levels=[latent_target],
        colors="black",
        linestyles="--",
        linewidths=1.6,
        zorder=6,
    )
    ax.contour(
        ratio_grid,
        gates_grid,
        latent_mean - latent_std,
        levels=[latent_target],
        colors="black",
        linestyles="--",
        linewidths=1.6,
        zorder=6,
    )

    if SHOW_Q56_REFERENCE_LINE:
        add_reference_contour(
            ax,
            grid_path=Q56_REFERENCE_GRID_PATH,
            label=Q56_REFERENCE_LABEL,
            color=Q56_REFERENCE_COLOR,
            linestyle=Q56_REFERENCE_LINESTYLE,
        )

    if SHOW_Q56_REFERENCE_POINTS:
        add_reference_points(
            ax,
            json_path=Q56_REFERENCE_JSON_PATH,
            label=Q56_REFERENCE_POINTS_LABEL,
        )

    if len(point_gates) > 0:
        regular_mask = ~point_is_sobol
        regular_failure = regular_mask & ~point_outcomes
        regular_success = regular_mask & point_outcomes
        sobol_failure = point_is_sobol & ~point_outcomes
        sobol_success = point_is_sobol & point_outcomes

        ax.scatter(
            point_ratios[regular_failure],
            point_gates[regular_failure],
            marker="x",
            s=36,
            color="black",
            linewidths=1.5,
            label="Failure",
            zorder=8,
        )
        ax.scatter(
            point_ratios[regular_success],
            point_gates[regular_success],
            marker="o",
            s=42,
            facecolors="white",
            edgecolors="black",
            linewidths=1.2,
            label="Success",
            zorder=9,
        )
        if SHOW_SOBOL_POINTS and np.any(sobol_failure):
            ax.scatter(
                point_ratios[sobol_failure],
                point_gates[sobol_failure],
                marker="x",
                s=42,
                color=SOBOL_FAILURE_COLOR,
                linewidths=1.7,
                label="Sobol failure",
                zorder=10,
            )
        if SHOW_SOBOL_POINTS and np.any(sobol_success):
            ax.scatter(
                point_ratios[sobol_success],
                point_gates[sobol_success],
                marker="o",
                s=48,
                facecolors=SOBOL_SUCCESS_COLOR,
                edgecolors="black",
                linewidths=1.0,
                label="Sobol success",
                zorder=11,
            )

    x_min = float(np.min(gates_grid[gates_grid > 0.0]))
    x_max = float(np.max(gates_grid))

    ax.plot([], [], color="black", linewidth=2.2, label=f"GP mean p={target:g}")
    ax.plot([], [], color="black", linestyle="--", linewidth=1.6, label=r"$\mu \pm 1\sigma$")
    ax.set_yscale("log")
    if RATIO_AXIS_LIMITS is None:
        ax.set_xlim(float(np.min(ratio_grid)), float(np.max(ratio_grid)))
    else:
        ax.set_xlim(float(RATIO_AXIS_LIMITS[0]), float(RATIO_AXIS_LIMITS[1]))
    if GATES_AXIS_LIMITS is None:
        ax.set_ylim(bottom=max(1.0, x_min), top=x_max)
    else:
        ax.set_ylim(bottom=max(1.0, float(GATES_AXIS_LIMITS[0])), top=float(GATES_AXIS_LIMITS[1]))
    ax.set_xlabel("Two-qubit gate ratio")
    ax.set_ylabel("Total gates")
    ax.set_title(
        f"FLE GP level set{qubit_title_part(json_path)} | S1={scores['S1']:.4g}, "
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
