"""Plot and score one saved FLE GP grid."""

from __future__ import annotations
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import ndtri


# -------------------------------------------------------------------------
# PATH / SCORE HANDLES
# -------------------------------------------------------------------------

# out_path = Path(r"Personal\Data\Paper_plots\plot_2\plot2.pdf")
out_path = Path(r"Personal\Data\Paper_plots\plot_2\plot2.png")

# -------------------------------------------------------------------------
# Panel-(a)-with_Fake_H
# -------------------------------------------------------------------------

Q26_fakeH_JSON_PATH = Path(
    r"Personal\FLE\H_wrapper\H2_1"
    r"\q26\seed_42\FLE_H_wrapper_20260827_151619\measurement_017_globalsur_20260831_225325_588812_actual_gates.json"
)

Q26_fakeH_GP_GRID_PATH = Path(
    r"Personal\FLE\H_wrapper\H2_1\q26"
    r"\seed_42\FLE_H_wrapper_20260827_151619\measurement_017_globalsur_20260831_225325_588812_actual_gates_gp_grid.npz"
)


Q26_fakeH_H21E_JSON_PATH = Path(
    r"Personal\FLE\Fancy_emulator\V_warp\H2_1E"
    r"\q26\seed_42\FLE_V_warp_20260827_121746"
    r"\reconstructed_native_gateset_measurement_019_globalsur_20260831_075831_968266_actual_gates.json"
)

Q26_fakeH_H21E_GP_GRID_PATH = Path(
    r"Personal\FLE\Fancy_emulator\V_warp\H2_1E\q26\seed_42\FLE_V_warp_20260827_121746"
    r"\reconstructed_native_gateset_measurement_019_globalsur_20260831_075831_968266_actual_gates_gp_grid.npz"
)


# -------------------------------------------------------------------------
# Panel-(b)-without_Fake_H
# -------------------------------------------------------------------------


Q26_accumulated_ordinary_JSON_PATH: Path | None = Path(
    r"Personal\Data\accumulated\H2-1"
    r"\accumulated_actualgr_H2-1_fle_costaware_20260813_175523_q_slices"
    r"\accumulated_actualgr_H2-1_fle_costaware_20260813_175523_q26.json"
)
Q26_accumulated_ordinary_GP_GRID_PATH: Path | None = Path(
    r"Personal\Data\accumulated\H2-1"
    r"\accumulated_actualgr_H2-1_fle_costaware_20260813_175523_q_slices"
    r"\accumulated_actualgr_H2-1_fle_costaware_20260813_175523_q26_gp_grid.npz"
)


GATES_AXIS_LIMITS: tuple[float, float] | None = (200.0, 1700.0)
RATIO_AXIS_LIMITS: tuple[float, float] | None = (0.1, 0.75)

SHOW_PREDICTED_FIDELITY_HUE = False
SHOW_CURRENT_MEASURED_POINTS_h21 = False
SHOW_CURRENT_MEASURED_POINTS_h21E = False
SHOW_SOBOL_POINTS = False
SOBOL_SUCCESS_COLOR = "tab:blue"
SOBOL_FAILURE_COLOR = "tab:orange"

CURRENT_SIGMA_BAND_ALPHA = 0.55
CURRENT_DATA_LABEL_PREFIX = "measured data"


def plot_target_from_grid(grid) -> float:
    return float(np.asarray(grid["target"]).item())


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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    sobol_keys = sobol_config_keys(json_path)
    gates: list[float] = []
    ratios: list[float] = []
    outcomes: list[int] = []
    is_sobol: list[bool] = []
    records = payload.get("data", [])

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


def measurement_data_counts(json_path: Path) -> tuple[int, int]:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    point_count = 0
    observation_count = 0
    for record in payload.get("data", []):
        n_1q = int(record["n_1qb_gates"])
        n_2q = int(record["n_2qb_gates"])
        if n_1q + n_2q <= 0:
            continue
        total = sum(int(count) for _, count in record.get("results", []))
        if total <= 0:
            continue
        point_count += 1
        observation_count += total
    return point_count, observation_count


def measurement_count_legend_label(json_path: Path, *, prefix: str) -> str:
    _, observation_count = measurement_data_counts(json_path)
    return f"{prefix}: {observation_count} obs"


def add_measurement_count_legend(
    ax,
    *,
    json_path: Path | None,
    prefix: str,
) -> None:
    if json_path is None or not json_path.exists():
        return
    ax.plot(
        [],
        [],
        color="none",
        label=measurement_count_legend_label(json_path, prefix=prefix),
    )


def add_sigma_band(
    ax,
    *,
    ratio_grid: np.ndarray,
    gates_grid: np.ndarray,
    latent_mean: np.ndarray,
    latent_std: np.ndarray,
    latent_target: float,
    color: str,
    alpha: float,
    label: str | None,
    zorder: int,
) -> None:
    band = np.abs(latent_mean - latent_target) - latent_std
    finite_band = band[np.isfinite(band)]
    if finite_band.size == 0:
        return
    band_min = float(np.min(finite_band))
    if band_min >= 0.0:
        return
    ax.contourf(
        ratio_grid,
        gates_grid,
        band,
        levels=[band_min, 0.0],
        colors=[color],
        alpha=alpha,
        zorder=zorder,
    )
    if label is not None:
        ax.plot([], [], color=color, linewidth=8.0, alpha=alpha, label=label)


def plot_fle_grid(
    json_paths: list[str | Path] = [Q26_fakeH_JSON_PATH, Q26_fakeH_H21E_JSON_PATH],
    grid_paths: list[str | Path] = [Q26_fakeH_GP_GRID_PATH, Q26_fakeH_H21E_GP_GRID_PATH],
    out_path: str | Path | None = out_path,
) -> dict[str, float]:

    # Data on H2-1
    grid_path = Path(grid_paths[0])
    json_path = Path(json_paths[0])

    grid = np.load(grid_path)
    gates_grid = np.asarray(grid["gates_grid"], dtype=float)
    ratio_grid = np.asarray(grid["ratio_grid"], dtype=float)
    probabilities = np.asarray(grid["probabilities"], dtype=float)
    latent_mean = np.asarray(grid["latent_mean"], dtype=float)
    latent_variance = np.maximum(np.asarray(grid["latent_variance"], dtype=float), 0.0)
    latent_std = np.sqrt(latent_variance)
    target = plot_target_from_grid(grid)
    latent_target = float(ndtri(target))
    print(f"[plot] grid: {grid_path}")
    print(f"[plot] contour target: p={target:g}")

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

    add_sigma_band(
        ax,
        ratio_grid=ratio_grid,
        gates_grid=gates_grid,
        latent_mean=latent_mean,
        latent_std=latent_std,
        latent_target=latent_target,
        color="navy",
        alpha=CURRENT_SIGMA_BAND_ALPHA,
        label=r"H2-1 Uncertainty",
        zorder=3,
    )

    ax.contour(
        ratio_grid,
        gates_grid,
        latent_mean,
        levels=[latent_target],
        colors="black",
        linestyles="-",
        linewidths=2.2,
        zorder=6,
    )

    if SHOW_CURRENT_MEASURED_POINTS_h21 and len(point_gates) > 0:
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

        ax.plot(
            [],
            [],
            color="none",
            label=measurement_count_legend_label(
                json_path,
                prefix=CURRENT_DATA_LABEL_PREFIX,
            ),
        )
    ax.plot([], [], color="black", linestyle="-", linewidth=2.2, label=f"H2-1")

    # Data on H2-1E
    grid_path = Path(grid_paths[1])
    json_path = Path(json_paths[1])

    grid = np.load(grid_path)
    gates_grid = np.asarray(grid["gates_grid"], dtype=float)
    ratio_grid = np.asarray(grid["ratio_grid"], dtype=float)
    probabilities = np.asarray(grid["probabilities"], dtype=float)
    latent_mean = np.asarray(grid["latent_mean"], dtype=float)
    latent_variance = np.maximum(np.asarray(grid["latent_variance"], dtype=float), 0.0)
    latent_std = np.sqrt(latent_variance)
    target = plot_target_from_grid(grid)
    latent_target = float(ndtri(target))
    print(f"[plot] grid: {grid_path}")
    print(f"[plot] contour target: p={target:g}")

    point_gates, point_ratios, point_outcomes, point_is_sobol = load_points(json_path)

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

    add_sigma_band(
        ax,
        ratio_grid=ratio_grid,
        gates_grid=gates_grid,
        latent_mean=latent_mean,
        latent_std=latent_std,
        latent_target=latent_target,
        color="cornflowerblue",
        alpha=CURRENT_SIGMA_BAND_ALPHA,
        label=r"H2-1E Uncertainty",
        zorder=3,
    )

    ax.contour(
        ratio_grid,
        gates_grid,
        latent_mean,
        levels=[latent_target],
        colors="black",
        linestyles="--",
        linewidths=2.2,
        zorder=6,
    )

    if SHOW_CURRENT_MEASURED_POINTS_h21E and len(point_gates) > 0:
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

        ax.plot(
            [],
            [],
            color="none",
            label=measurement_count_legend_label(
                json_path,
                prefix=CURRENT_DATA_LABEL_PREFIX,
            ),
        )
    ax.plot([], [], color="black", linestyle="--", linewidth=2.2, label=f"H2-1E")

    ax.set_yscale("log")
    ax.set_xlim(float(RATIO_AXIS_LIMITS[0]), float(RATIO_AXIS_LIMITS[1]))
    ax.set_ylim(bottom=max(1.0, float(GATES_AXIS_LIMITS[0])), top=float(GATES_AXIS_LIMITS[1]))
    ax.set_xlabel("Two-qubit gate ratio")
    ax.set_ylabel("Total gates")
    ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.9)
    ax.set_title("With Fake Hadamard; Add SympleQ", fontsize=10, fontweight="bold")
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"[saved] plot: {out_path}")
    plt.show()
    # return scores


if __name__ == "__main__":

    plot_fle_grid(
        # json_path=Q26_fakeH_JSON_PATH,
        # grid_path=Q26_fakeH_GP_GRID_PATH,
        # out_path=out_path
    )
