"""Static Matplotlib 3D isosurface plot for one saved 3D FLE GP grid."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.measure import marching_cubes

from fantasy_levelset_settings_3d import HQC_BUDGET


# -------------------------------------------------------------------------
# PATH / PLOT HANDLES
# -------------------------------------------------------------------------

RMB_JSON_PATH = Path(
    r"Personal\seed_2025\rick_fantasy_gpu_crossing_20260623_184228.json"
)

GP_GRID_PATH: Path | None = None
PNG_PATH: Path | None = None

SHOW_MEASURED_POINTS = True
LAST_BACKEND_BATCH_SIZE: int | None = None

ISOSURFACE_ALPHA = 0.1
SHOW_ONE_SIGMA_SURFACES = True
ONE_SIGMA_ALPHA = 0.1
FIGSIZE = (8.8, 7.0)
DPI = 220


def sibling_grid_path(json_path: Path) -> Path:
    return json_path.parent / f"{json_path.stem}_gp_grid_3d.npz"


def sibling_png_path(json_path: Path) -> Path:
    return json_path.parent / f"{json_path.stem}_fle_isosurface_3d_matplotlib.png"


def infer_n_qubits(record: dict) -> int | None:
    for key in ("n_qubits", "qubits", "n_qb", "num_qubits"):
        if key in record and record[key] is not None:
            value = int(record[key])
            if value > 0:
                return value
    return None


def load_points(
    json_path: Path,
    *,
    last_backend_batch_size: int | None = LAST_BACKEND_BATCH_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    payload = json.loads(json_path.read_text(encoding="utf-8"))

    gates: list[float] = []
    ratios: list[float] = []
    qubits: list[float] = []
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

        n_qubits = infer_n_qubits(record)
        if n_qubits is None:
            continue

        counts = {int(outcome): int(count) for outcome, count in record["results"]}
        successes = counts.get(1, 0)
        failures = counts.get(0, 0)

        if successes + failures <= 0:
            continue

        gates.append(float(total))
        ratios.append(float(n_2q / total))
        qubits.append(float(n_qubits))
        outcomes.append(int(successes >= failures))

    return (
        np.asarray(gates, dtype=float),
        np.asarray(ratios, dtype=float),
        np.asarray(qubits, dtype=float),
        np.asarray(outcomes, dtype=bool),
    )


def axis_from_grid(grid_array: np.ndarray, axis: int) -> np.ndarray:
    """
    Extract a 1D coordinate axis from a meshgrid-like 3D array.

    Expected saved grid shape:
        (n_qubits_grid, n_ratio_grid, n_gates_grid)

    axes:
        axis=0 -> qubits
        axis=1 -> ratio
        axis=2 -> gates
    """

    if axis == 0:
        return grid_array[:, 0, 0]
    if axis == 1:
        return grid_array[0, :, 0]
    if axis == 2:
        return grid_array[0, 0, :]
    raise ValueError(f"axis must be 0, 1, or 2, got {axis}")


def interp_axis(axis_values: np.ndarray, indices: np.ndarray) -> np.ndarray:
    """
    Convert marching-cubes fractional voxel indices into physical coordinates.
    """

    base = np.arange(len(axis_values), dtype=float)
    return np.interp(indices, base, axis_values)


def plot_fle_isosurface_3d_matplotlib(
    json_path: str | Path = RMB_JSON_PATH,
    grid_path: str | Path | None = GP_GRID_PATH,
    png_path: str | Path | None = PNG_PATH,
) -> None:
    json_path = Path(json_path)
    grid_path = sibling_grid_path(json_path) if grid_path is None else Path(grid_path)
    png_path = sibling_png_path(json_path) if png_path is None else Path(png_path)

    grid = np.load(grid_path)

    gates_grid = np.asarray(grid["gates_grid"], dtype=float)
    ratio_grid = np.asarray(grid["ratio_grid"], dtype=float)
    qubits_grid = np.asarray(grid["qubits_grid"], dtype=float)
    probabilities = np.asarray(grid["probabilities"], dtype=float)
    latent_mean = np.asarray(grid["latent_mean"], dtype=float)
    latent_std = np.sqrt(
        np.clip(np.asarray(grid["latent_variance"], dtype=float), 0.0, None)
    )
    target = float(np.asarray(grid["target"]).item())

    if probabilities.ndim != 3:
        raise ValueError(
            f"Expected a 3D probability grid, got ndim={probabilities.ndim}. "
            "You probably loaded an old 2D slice grid."
        )

    if not (np.nanmin(probabilities) <= target <= np.nanmax(probabilities)):
        raise ValueError(
            "Target is outside the predicted probability range. "
            f"target={target}, range=({np.nanmin(probabilities)}, {np.nanmax(probabilities)})"
        )

    # Saved grid shape should be:
    #   probabilities[q_index, ratio_index, gates_index]
    qubits_axis = axis_from_grid(qubits_grid, axis=0)
    ratio_axis = axis_from_grid(ratio_grid, axis=1)
    gates_axis = axis_from_grid(gates_grid, axis=2)
    log_gates_axis = np.log10(gates_axis)

    # marching_cubes works in array-index coordinates:
    #   vertex[:, 0] = q index
    #   vertex[:, 1] = ratio index
    #   vertex[:, 2] = gates index
    verts, faces, normals, values = marching_cubes(
        probabilities,
        level=target,
    )

    q_coords = interp_axis(qubits_axis, verts[:, 0])
    ratio_coords = interp_axis(ratio_axis, verts[:, 1])
    log_gate_coords = interp_axis(log_gates_axis, verts[:, 2])

    # Matplotlib wants vertices as (x, y, z).
    # We choose:
    #   x = log10(total gates)
    #   y = two-qubit gate ratio
    #   z = n_qubits
    surface_vertices = np.column_stack(
        [
            log_gate_coords,
            ratio_coords,
            q_coords,
        ]
    )

    mesh = Poly3DCollection(
        surface_vertices[faces],
        alpha=ISOSURFACE_ALPHA,
        linewidths=0.15,
    )
    mesh.set_facecolor("0.25")
    mesh.set_edgecolor("0.12")

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_subplot(111, projection="3d")
    ax.computed_zorder = False
    ax.add_collection3d(mesh)
    mesh.set_zorder(1)

    if SHOW_ONE_SIGMA_SURFACES:
        latent_target = NormalDist().inv_cdf(target)
        for label, level_set, color in [
            ("GP latent mean - 1 sigma", latent_mean - latent_std, "tab:blue"),
            ("GP latent mean + 1 sigma", latent_mean + latent_std, "tab:red"),
        ]:
            if not (
                np.nanmin(level_set) <= latent_target <= np.nanmax(level_set)
            ):
                print(f"[warning] skipped {label}: it does not cross the grid")
                continue

            sigma_verts, sigma_faces, _, _ = marching_cubes(
                level_set,
                level=latent_target,
            )
            sigma_vertices = np.column_stack(
                [
                    interp_axis(log_gates_axis, sigma_verts[:, 2]),
                    interp_axis(ratio_axis, sigma_verts[:, 1]),
                    interp_axis(qubits_axis, sigma_verts[:, 0]),
                ]
            )
            sigma_mesh = Poly3DCollection(
                sigma_vertices[sigma_faces],
                alpha=ONE_SIGMA_ALPHA,
                linewidths=0.0,
            )
            sigma_mesh.set_facecolor(color)
            sigma_mesh.set_edgecolor("none")
            ax.add_collection3d(sigma_mesh)
            sigma_mesh.set_zorder(2)
            ax.plot([], [], [], color=color, linewidth=5, label=label)

    if SHOW_MEASURED_POINTS:
        point_gates, point_ratios, point_qubits, point_outcomes = load_points(json_path)

        if len(point_gates) > 0:
            ax.scatter(
                np.log10(point_gates[~point_outcomes]),
                point_ratios[~point_outcomes],
                point_qubits[~point_outcomes],
                marker="x",
                color="black",
                s=36,
                linewidths=1.4,
                depthshade=False,
                zorder=20,
                label="Failure",
            )
            ax.scatter(
                np.log10(point_gates[point_outcomes]),
                point_ratios[point_outcomes],
                point_qubits[point_outcomes],
                marker="o",
                s=34,
                facecolors="white",
                edgecolors="black",
                linewidths=1.0,
                depthshade=False,
                zorder=20,
                label="Success",
            )

    ax.set_xlim(float(np.min(log_gates_axis)), float(np.max(log_gates_axis)))
    ax.set_ylim(float(np.min(ratio_axis)), float(np.max(ratio_axis)))
    ax.set_zlim(float(np.min(qubits_axis)), float(np.max(qubits_axis)))

    ax.set_xlabel("log10(total gates)")
    ax.set_ylabel("Two-qubit gate ratio")
    ax.set_zlabel("n_qubits")

    ax.set_title(
        f"3D FLE GP level set | HQC={1323} | "
        f"P(success)={target:g}"
    )

    # Useful view angle similar to your screenshot.
    ax.view_init(elev=24, azim=38)

    # Dummy legend entry for the surface.
    ax.plot([], [], [], color="0.25", linewidth=6, label=f"GP P(success)={target:g}")
    if SHOW_MEASURED_POINTS:
        ax.legend(loc="best", fontsize=8)
    else:
        ax.legend(loc="best", fontsize=8)

    fig.tight_layout()

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=DPI, bbox_inches="tight")

    print(f"[grid] loaded: {grid_path}")
    print(f"[grid] probabilities shape: {probabilities.shape}")
    print(
        "[grid] probability range: "
        f"{float(np.nanmin(probabilities)):.6g} to "
        f"{float(np.nanmax(probabilities)):.6g}"
    )
    print(f"[grid] target: {target}")
    print(f"[saved] Matplotlib 3D isosurface: {png_path}")

    plt.show()


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else RMB_JSON_PATH
    plot_fle_isosurface_3d_matplotlib(path)
