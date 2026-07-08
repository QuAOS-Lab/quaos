"""Interactive 3D isosurface plot for one saved 3D FLE GP grid."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np
import plotly.graph_objects as go

# from fantasy_levelset_settings_3d import HQC_BUDGET


# -------------------------------------------------------------------------
# PATH / PLOT HANDLES
# -------------------------------------------------------------------------

RMB_JSON_PATH = Path(
    r"Personal\seed_2025\rick_fantasy_gpu_crossing_20260623_184228.json"
)

GP_GRID_PATH: Path | None = None
HTML_PATH: Path | None = None

LAST_BACKEND_BATCH_SIZE: int | None = None
# Turn this off for a cleaner diagnostic.
SHOW_VOLUME_BACKGROUND = False

# Turn this off if you only want the GP level-set surface.
SHOW_MEASURED_POINTS = True
SHOW_FAILURE_POINTS = True
SHOW_SUCCESS_POINTS = True
SHOW_ONE_SIGMA_SURFACES = False

# Probability volume settings, only used if SHOW_VOLUME_BACKGROUND = True.
VOLUME_OPACITY = 0.1
VOLUME_SURFACE_COUNT = 12

# Main level-set surface.
ISOSURFACE_OPACITY = 0.1
ISOSURFACE_WIDTH = 1e-3

ONE_SIGMA_SURFACE_OPACITY = 0.28

# -------------------------------------------------------------------------
# PATH HELPERS
# -------------------------------------------------------------------------


def sibling_grid_path(json_path: Path) -> Path:
    """Return the expected sibling 3D GP grid path."""

    return json_path.parent / f"{json_path.stem}_gp_grid_3d.npz"


def sibling_html_path(json_path: Path) -> Path:
    """Return the expected output HTML path."""

    return json_path.parent / f"{json_path.stem}_fle_isosurface_3d.html"


def json_path_from_grid_path(grid_path: Path) -> Path | None:
    """Return the matching JSON path for a saved grid, if one exists."""

    name = grid_path.name
    for suffix in ("_gp_grid_3d.npz", "_gp_grid.npz"):
        if name.endswith(suffix):
            candidate = grid_path.parent / f"{name[:-len(suffix)]}.json"
            return candidate if candidate.exists() else None
    return None


def html_path_from_grid_path(grid_path: Path) -> Path:
    """Return the output HTML path for a direct grid input."""

    return grid_path.parent / f"{grid_path.stem}_fle_isosurface_3d.html"


def slice_json_paths_for_grid(grid_path: Path) -> list[Path]:
    """Return fixed-qubit slice JSON files next to a stacked slice grid."""

    return sorted(grid_path.parent.glob("*_q*.json"))


# -------------------------------------------------------------------------
# LOAD MEASURED POINTS
# -------------------------------------------------------------------------


def infer_n_qubits(record: dict) -> int | None:
    """
    Try to infer n_qubits from one RMB JSON record.

    Different saved formats may use different keys.
    """

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
    """
    Load measured RMB points from the saved JSON.

    Returns:
        point_gates
        point_ratios
        point_qubits
        point_outcomes

    point_outcomes:
        True  = success
        False = failure
    """

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

        counts = {
            int(outcome): int(count)
            for outcome, count in record["results"]
        }

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


def load_points_from_jsons(json_paths: list[Path]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load measured points from several fixed-qubit slice JSON files."""

    chunks = [load_points(path) for path in json_paths if path.exists()]
    if not chunks:
        return (
            np.asarray([], dtype=float),
            np.asarray([], dtype=float),
            np.asarray([], dtype=float),
            np.asarray([], dtype=bool),
        )
    return tuple(np.concatenate(parts) for parts in zip(*chunks))


# -------------------------------------------------------------------------
# GRID LOAD / VALIDATION
# -------------------------------------------------------------------------


def load_3d_grid(grid_path: Path) -> dict[str, np.ndarray | float]:
    """
    Load and validate a saved 3D GP grid.

    Expected arrays:
        gates_grid
        ratio_grid
        qubits_grid
        probabilities
        latent_mean
        latent_variance
        target
    """

    grid = np.load(grid_path)

    required_keys = [
        "gates_grid",
        "ratio_grid",
        "qubits_grid",
        "probabilities",
        "target",
    ]

    missing = [key for key in required_keys if key not in grid.files]
    if missing:
        raise KeyError(
            f"Missing keys in {grid_path}: {missing}. "
            f"Available keys: {grid.files}"
        )

    gates_grid = np.asarray(grid["gates_grid"], dtype=float)
    ratio_grid = np.asarray(grid["ratio_grid"], dtype=float)
    qubits_grid = np.asarray(grid["qubits_grid"], dtype=float)
    probabilities = np.asarray(grid["probabilities"], dtype=float)
    latent_mean = (
        np.asarray(grid["latent_mean"], dtype=float)
        if "latent_mean" in grid.files
        else None
    )
    latent_variance = (
        np.asarray(grid["latent_variance"], dtype=float)
        if "latent_variance" in grid.files
        else None
    )
    target = float(np.asarray(grid["target"]).item())
    print(target)

    if probabilities.ndim != 3:
        raise ValueError(
            "This plotter expects a 3D grid. "
            f"Got probabilities.ndim={probabilities.ndim}. "
            "You probably loaded an old 2D slice grid."
        )

    expected_shape = probabilities.shape
    for name, arr in [
        ("gates_grid", gates_grid),
        ("ratio_grid", ratio_grid),
        ("qubits_grid", qubits_grid),
        ("latent_mean", latent_mean),
        ("latent_variance", latent_variance),
    ]:
        if arr is None:
            continue
        if arr.shape != expected_shape:
            raise ValueError(
                f"Shape mismatch: {name}.shape={arr.shape}, "
                f"probabilities.shape={expected_shape}."
            )

    if not np.nanmin(probabilities) <= target <= np.nanmax(probabilities):
        print(
            "[warning] target is outside predicted probability range: "
            f"target={target}, "
            f"range=({np.nanmin(probabilities)}, {np.nanmax(probabilities)})"
        )

    return {
        "gates_grid": gates_grid,
        "ratio_grid": ratio_grid,
        "qubits_grid": qubits_grid,
        "probabilities": probabilities,
        "latent_mean": latent_mean,
        "latent_variance": latent_variance,
        "target": target,
    }


# -------------------------------------------------------------------------
# MAIN PLOTTER
# -------------------------------------------------------------------------


def plot_fle_isosurface_3d(
    json_path: str | Path = RMB_JSON_PATH,
    grid_path: str | Path | None = GP_GRID_PATH,
    html_path: str | Path | None = HTML_PATH,
    *,
    show: bool = True,
) -> None:
    """
    Plot the 3D GP level set:

        P(success | total gates, ratio, n_qubits) = target

    Saves an interactive HTML file and opens the plot.
    """

    input_path = Path(json_path)
    measured_json_paths: list[Path] = []

    if input_path.suffix.lower() == ".npz":
        grid_path = input_path if grid_path is None else Path(grid_path)
        json_candidate = json_path_from_grid_path(grid_path)
        if json_candidate is not None:
            measured_json_paths = [json_candidate]
        else:
            measured_json_paths = slice_json_paths_for_grid(grid_path)
        html_path = html_path_from_grid_path(grid_path) if html_path is None else Path(html_path)
    else:
        json_path = input_path
        grid_path = sibling_grid_path(json_path) if grid_path is None else Path(grid_path)
        html_path = sibling_html_path(json_path) if html_path is None else Path(html_path)
        measured_json_paths = [json_path]

    loaded = load_3d_grid(grid_path)

    gates_grid = loaded["gates_grid"]
    ratio_grid = loaded["ratio_grid"]
    qubits_grid = loaded["qubits_grid"]
    probabilities = loaded["probabilities"]
    latent_mean = loaded["latent_mean"]
    latent_variance = loaded["latent_variance"]
    target = float(loaded["target"])

    # Plotly coordinates.
    # Use log10(total gates), because total gates spans decades.
    x = ratio_grid.ravel()
    y = qubits_grid.ravel()
    z = np.log10(gates_grid).ravel()
    p = probabilities.ravel()
    latent_mean_flat = None if latent_mean is None else latent_mean.ravel()
    latent_std_flat = (
        None
        if latent_variance is None
        else np.sqrt(np.clip(latent_variance, 0.0, None)).ravel()
    )

    finite_mask = (
        np.isfinite(x)
        & np.isfinite(y)
        & np.isfinite(z)
        & np.isfinite(p)
    )

    x = x[finite_mask]
    y = y[finite_mask]
    z = z[finite_mask]
    p = p[finite_mask]
    if latent_mean_flat is not None:
        latent_mean_flat = latent_mean_flat[finite_mask]
    if latent_std_flat is not None:
        latent_std_flat = latent_std_flat[finite_mask]

    fig = go.Figure()

    # ------------------------------------------------------------------
    # Optional background probability volume.
    # This is useful for context, but can make the level set harder to see.
    # ------------------------------------------------------------------

    if SHOW_VOLUME_BACKGROUND:
        fig.add_trace(
            go.Volume(
                x=x,
                y=y,
                z=z,
                value=p,
                isomin=0.05,
                isomax=0.95,
                opacity=VOLUME_OPACITY,
                surface_count=VOLUME_SURFACE_COUNT,
                colorscale="RdYlGn",
                colorbar=dict(title="P(success)"),
                name="GP probability volume",
                showscale=True,
            )
        )

    # ------------------------------------------------------------------
    # Main object: the 3D level-set / isosurface.
    # ------------------------------------------------------------------

    fig.add_trace(
        go.Isosurface(
            x=x,
            y=y,
            z=z,
            value=p,
            isomin=target - ISOSURFACE_WIDTH,
            isomax=target + ISOSURFACE_WIDTH,
            surface_count=1,
            opacity=ISOSURFACE_OPACITY,
            colorscale=[
                [0.0, "black"],
                [1.0, "black"],
            ],
            caps=dict(
                x_show=False,
                y_show=False,
                z_show=False,
            ),
            name=f"GP P(success) = {target:g}",
            showscale=False,
        )
    )

    if SHOW_ONE_SIGMA_SURFACES:
        if latent_mean_flat is None or latent_std_flat is None:
            print(
                "[warning] skipped +/-1 sigma surfaces: "
                "latent_mean/latent_variance not found in grid."
            )
        else:
            latent_target = NormalDist().inv_cdf(target)
            for label, values, color in [
                (
                    "GP latent mean - 1 sigma",
                    latent_mean_flat - latent_std_flat,
                    "royalblue",
                ),
                (
                    "GP latent mean + 1 sigma",
                    latent_mean_flat + latent_std_flat,
                    "firebrick",
                ),
            ]:
                if not np.nanmin(values) <= latent_target <= np.nanmax(values):
                    print(f"[warning] skipped {label}: it does not cross the grid")
                    continue

                fig.add_trace(
                    go.Isosurface(
                        x=x,
                        y=y,
                        z=z,
                        value=values,
                        isomin=latent_target - ISOSURFACE_WIDTH,
                        isomax=latent_target + ISOSURFACE_WIDTH,
                        surface_count=1,
                        opacity=ONE_SIGMA_SURFACE_OPACITY,
                        colorscale=[
                            [0.0, color],
                            [1.0, color],
                        ],
                        caps=dict(
                            x_show=False,
                            y_show=False,
                            z_show=False,
                        ),
                        name=label,
                        showscale=False,
                    )
                )

    # ------------------------------------------------------------------
    # Measured training points.
    # ------------------------------------------------------------------

    if SHOW_MEASURED_POINTS:
        point_gates, point_ratios, point_qubits, point_outcomes = load_points_from_jsons(
            measured_json_paths
        )

        if len(point_gates) == 0:
            print(
                "[warning] no measured points were loaded. "
                "Check whether the JSON records contain n_qubits."
            )
        else:
            failure_mask = ~point_outcomes
            success_mask = point_outcomes

            if SHOW_FAILURE_POINTS and np.any(failure_mask):
                fig.add_trace(
                    go.Scatter3d(
                        x=point_ratios[failure_mask],
                        y=point_qubits[failure_mask],
                        z=np.log10(point_gates[failure_mask]),
                        mode="markers",
                        marker=dict(
                            size=3,
                            symbol="x",
                            color="seagreen",
                            line=dict(width=2),
                        ),
                        name="Failure",
                    )
                )

            if SHOW_SUCCESS_POINTS and np.any(success_mask):
                fig.add_trace(
                    go.Scatter3d(
                        x=point_ratios[success_mask],
                        y=point_qubits[success_mask],
                        z=np.log10(point_gates[success_mask]),
                        mode="markers",
                        marker=dict(
                            size=5,
                            symbol="circle-open",
                            color="mediumpurple",
                            line=dict(width=1.5),
                        ),
                        name="Success",
                    )
                )

    # ------------------------------------------------------------------
    # Layout.
    # ------------------------------------------------------------------

    fig.update_layout(
        title=(
            f"3D FLE GP level set | HQC={1323} | "
            f"P(success)={target:g}"
        ),
        scene=dict(
            xaxis=dict(
                title="Two-qubit gate ratio",
                backgroundcolor="rgba(245,245,245,0.95)",
                gridcolor="lightgray",
            ),
            yaxis=dict(
                title="n_qubits",
                backgroundcolor="rgba(245,245,245,0.95)",
                gridcolor="lightgray",
            ),
            zaxis=dict(
                title="log10(total gates)",
                backgroundcolor="rgba(245,245,245,0.95)",
                gridcolor="lightgray",
            ),
        ),
        legend=dict(
            x=0.02,
            y=0.98,
            bgcolor="rgba(255,255,255,0.75)",
        ),
        margin=dict(l=0, r=0, t=55, b=0),
    )

    html_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(html_path)

    print(f"[grid] loaded: {grid_path}")
    print(f"[grid] probabilities shape: {probabilities.shape}")
    print(
        "[grid] probability range: "
        f"{float(np.nanmin(probabilities)):.6g} to "
        f"{float(np.nanmax(probabilities)):.6g}"
    )
    print(f"[grid] target: {target}")
    print(f"[saved] interactive 3D plot: {html_path}")

    if show:
        fig.show()


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else RMB_JSON_PATH
    plot_fle_isosurface_3d(path)
