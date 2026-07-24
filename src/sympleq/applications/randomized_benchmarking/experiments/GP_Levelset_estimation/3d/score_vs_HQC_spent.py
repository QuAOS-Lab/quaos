"""Plot mean 3D success-side volume by actual HQC spent for qband runs."""

from __future__ import annotations

import json
import re
import sys
import argparse
from pathlib import Path
from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np

from score3d import axis_spacing, score_grid


ROOT_FOLDER = Path(r"Personal\FLE")
QBAND_VALUES = range(5)
SEED_VALUES = range(2026, 2030)
SAVE_FIG_PATH = Path(
    r"Personal\RMB_results_figs\Volume\score_vs_hqc_spent_H2.pdf"
)
LABELS: list[str] = ['H2-1-seed-2026', 'H2-1-seed-2027', 'H2-1-seed-2028', 'H2-1-seed-2029']

Curve = tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]


def sibling_grid(
    json_path: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> Path:
    if grid_root is None:
        return json_path.parent / f"{json_path.stem}_gp_grid_3d.npz"
    if source_root is None:
        raise ValueError("source_root must be set when grid_root is set.")
    try:
        relative_json = json_path.relative_to(source_root)
    except ValueError:
        relative_json = json_path.resolve().relative_to(source_root.resolve())
    return grid_root / relative_json.parent / f"{json_path.stem}_gp_grid_3d.npz"


def infer_source_root(input_paths: list[Path], grid_root: Path | None) -> Path | None:
    if grid_root is None:
        return None
    if grid_root.name == "uniform_10_5000_grids":
        return grid_root.parent
    return input_paths[0] if input_paths else ROOT_FOLDER


def existing_grid_path(
    json_path: Path,
    *,
    grid_root: Path | None,
    source_root: Path | None,
) -> Path:
    grid_path = sibling_grid(json_path, grid_root=grid_root, source_root=source_root)
    return grid_path


def measurement_step(json_path: Path) -> int:
    match = re.match(r"measurement_(\d+)_", json_path.name)
    return int(match.group(1)) if match else 0


def sobol_step(json_path: Path) -> int | None:
    match = re.match(r"measurement_\d+_sobol_(\d+)_", json_path.name)
    return int(match.group(1)) if match else None


def plot_order_key(json_path: Path) -> tuple[int, int, str]:
    sobol = sobol_step(json_path)
    if sobol is not None:
        return 0, sobol, json_path.name
    return 1, measurement_step(json_path), json_path.name


def hqc_spent(json_path: Path) -> float:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    experiment = payload.get("experiment", {})
    spent = experiment.get("spent_hqc")
    return float(spent) if spent is not None else float(measurement_step(json_path))


def score_grid_with_sigma_bands(grid_path: Path) -> dict[str, float]:
    """Score the central, one-sigma, and two-sigma latent level-set volumes."""

    central = score_grid(grid_path)
    grid = np.load(grid_path)

    if "latent_mean" not in grid.files or "latent_variance" not in grid.files:
        score = central["success_side_volume"]
        return {
            "score": score,
            "lower_1sigma": score,
            "upper_1sigma": score,
            "lower_2sigma": score,
            "upper_2sigma": score,
        }

    latent_mean = np.asarray(grid["latent_mean"], dtype=float)
    latent_std = np.sqrt(
        np.clip(np.asarray(grid["latent_variance"], dtype=float), 0.0, None)
    )
    target = float(np.asarray(grid["target"]).item())
    latent_target = NormalDist().inv_cdf(target)
    voxel_volume = (
        axis_spacing(np.asarray(grid["gates_grid"], dtype=float), log=True)
        * axis_spacing(np.asarray(grid["ratio_grid"], dtype=float))
        * axis_spacing(np.asarray(grid["qubits_grid"], dtype=float))
    )

    valid = np.isfinite(latent_mean) & np.isfinite(latent_std)

    def volume_at_sigma(n_sigma: float) -> float:
        surface = latent_mean + n_sigma * latent_std
        return float(np.sum(valid & (surface >= latent_target)) * voxel_volume)

    return {
        "score": central["success_side_volume"],
        "lower_1sigma": volume_at_sigma(-1.0),
        "upper_1sigma": volume_at_sigma(1.0),
        "lower_2sigma": volume_at_sigma(-2.0),
        "upper_2sigma": volume_at_sigma(2.0),
    }


def measurement_pairs(
    run_folder: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> list[tuple[Path, Path]]:
    pairs = []
    for json_path in run_folder.glob("measurement_*.json"):
        grid_path = existing_grid_path(
            json_path,
            grid_root=grid_root,
            source_root=source_root,
        )
        if grid_path.exists():
            pairs.append((json_path, grid_path))
    return sorted(pairs, key=lambda pair: plot_order_key(pair[0]))


def sobol_done_hqc(
    run_folder: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> float | None:
    sobol_hqc = [
        hqc_spent(json_path)
        for json_path, _ in measurement_pairs(
            run_folder,
            grid_root=grid_root,
            source_root=source_root,
        )
        if sobol_step(json_path) is not None
    ]
    return sobol_hqc[-1] if sobol_hqc else None


def latest_seed_run(
    seed_folder: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> Path | None:
    candidates = [
        run_folder
        for run_folder in seed_folder.iterdir()
        if run_folder.is_dir()
        and measurement_pairs(
            run_folder,
            grid_root=grid_root,
            source_root=source_root,
        )
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def run_curve(
    run_folder: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> Curve | None:
    xs: list[int] = []
    hqc: list[float] = []
    scores: list[float] = []
    lower_1sigma: list[float] = []
    upper_1sigma: list[float] = []
    lower_2sigma: list[float] = []
    upper_2sigma: list[float] = []

    pairs = measurement_pairs(
        run_folder,
        grid_root=grid_root,
        source_root=source_root,
    )
    if not pairs:
        print(f"[skip] no measurement/grid pairs in {run_folder}", flush=True)
        return None

    for submission_index, (json_path, grid_path) in enumerate(pairs, start=1):
        xs.append(submission_index)
        hqc.append(hqc_spent(json_path))
        scored = score_grid_with_sigma_bands(grid_path)
        scores.append(scored["score"])
        lower_1sigma.append(scored["lower_1sigma"])
        upper_1sigma.append(scored["upper_1sigma"])
        lower_2sigma.append(scored["lower_2sigma"])
        upper_2sigma.append(scored["upper_2sigma"])
        print(
            f"  {json_path.name}: stitched={submission_index}, "
            f"grid={grid_path.name}, "
            f"hqc={hqc[-1]:.6g}, score={scores[-1]:.6g}, "
            f"1sigma=({lower_1sigma[-1]:.6g}, {upper_1sigma[-1]:.6g}), "
            f"2sigma=({lower_2sigma[-1]:.6g}, {upper_2sigma[-1]:.6g})",
            flush=True,
        )

    return (
        np.asarray(xs, dtype=float),
        np.asarray(hqc, dtype=float),
        np.asarray(scores, dtype=float),
        np.asarray(lower_1sigma, dtype=float),
        np.asarray(upper_1sigma, dtype=float),
        np.asarray(lower_2sigma, dtype=float),
        np.asarray(upper_2sigma, dtype=float),
    )


def seed_curve(
    seed_folder: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> Curve | None:
    run_folder = latest_seed_run(
        seed_folder,
        grid_root=grid_root,
        source_root=source_root,
    )
    if run_folder is None:
        print(f"[skip] no usable run folder in {seed_folder}", flush=True)
        return None

    print(f"[seed] {seed_folder.name}: {run_folder}", flush=True)
    return run_curve(run_folder, grid_root=grid_root, source_root=source_root)


def run_label(run_folder: Path) -> str:
    seed_folder = run_folder.parent.name
    qband_folder = run_folder.parent.parent.name
    return f"{qband_folder} {seed_folder} {run_folder.name}"


def qband_folder_from_arg(path: Path) -> list[Path]:
    if path.name.startswith("qband_"):
        return [path]

    folders = [path / f"qband_{qband}" for qband in QBAND_VALUES]
    return [folder for folder in folders if folder.exists()]


def qband_label(qband_folder: Path) -> str:
    match = re.match(r"qband_(\d+)$", qband_folder.name)
    if match is None:
        return qband_folder.name.replace("_", " ")

    band_length = int(match.group(1)) + 1
    unit = "qubit" if band_length == 1 else "qubits"
    return f"band length = {band_length} {unit}"


def plot_run_folder(
    run_folder: Path,
    label: str | None = None,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> bool:
    curve = run_curve(run_folder, grid_root=grid_root, source_root=source_root)
    if curve is None:
        return False

    _, hqc, scores, lower_1sigma, upper_1sigma, lower_2sigma, upper_2sigma = curve
    line = plt.plot(hqc, scores, "o-", label=label or run_label(run_folder))[0]
    color = line.get_color()
    # plt.fill_between(
    #     hqc,
    #     lower_2sigma,
    #     upper_2sigma,
    #     color=color,
    #     alpha=0.08,
    #     linewidth=0,
    #     label=f"{label or run_label(run_folder)} 2 sigma",
    # )
    plt.fill_between(
        hqc,
        lower_1sigma,
        upper_1sigma,
        color=color,
        alpha=0.16,
        linewidth=0,
        label=f"{label or run_label(run_folder)} 1 sigma",
    )
    sobol_hqc = sobol_done_hqc(
        run_folder,
        grid_root=grid_root,
        source_root=source_root,
    )
    if sobol_hqc is not None:
        plt.axvline(
            sobol_hqc,
            color=color,
            linestyle="--",
            linewidth=1.0,
        )
    print(
        f"[run] {run_folder}: points={len(scores)}, "
        f"hqc_range=({hqc[0]:.6g}, {hqc[-1]:.6g}), "
        f"final_score={scores[-1]:.6g}",
        flush=True,
    )
    return True


def qband_curves(
    qband_folder: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> list[Curve]:
    curves = []
    for seed in SEED_VALUES:
        seed_folder = qband_folder / f"seed_{seed}"
        if not seed_folder.exists():
            print(f"[skip] missing {seed_folder}", flush=True)
            continue
        curve = seed_curve(
            seed_folder,
            grid_root=grid_root,
            source_root=source_root,
        )
        if curve is not None:
            curves.append(curve)
    if not curves:
        return curves

    return curves


def mean_bands_by_hqc(
    curves: list[Curve],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    max_len = max(len(scores) for _, _, scores, *_ in curves)
    hqc_min = max(float(np.min(hqc)) for _, hqc, *_ in curves)
    hqc_max = min(float(np.max(hqc)) for _, hqc, *_ in curves)
    hqc_axis = np.linspace(hqc_min, hqc_max, max_len)
    matrices = [
        np.full((len(curves), max_len), np.nan, dtype=float)
        for _ in range(5)
    ]

    for row, (_, hqc, scores, lower_1sigma, upper_1sigma, lower_2sigma, upper_2sigma) in enumerate(curves):
        order = np.argsort(hqc)
        unique_hqc, unique_indices = np.unique(hqc[order], return_index=True)
        for matrix, values in zip(
            matrices,
            [scores, lower_1sigma, upper_1sigma, lower_2sigma, upper_2sigma],
        ):
            unique_values = values[order][unique_indices]
            matrix[row] = np.interp(hqc_axis, unique_hqc, unique_values)

    mean, lower_1, upper_1, lower_2, upper_2 = [
        np.nanmean(matrix, axis=0)
        for matrix in matrices
    ]
    return hqc_axis, mean, lower_1, upper_1, lower_2, upper_2


def qband_sobol_done_hqc(
    qband_folder: Path,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> float | None:
    sobol_hqc = []
    for seed in SEED_VALUES:
        seed_folder = qband_folder / f"seed_{seed}"
        if not seed_folder.exists():
            continue
        run_folder = latest_seed_run(
            seed_folder,
            grid_root=grid_root,
            source_root=source_root,
        )
        if run_folder is None:
            continue
        value = sobol_done_hqc(
            run_folder,
            grid_root=grid_root,
            source_root=source_root,
        )
        if value is not None:
            sobol_hqc.append(value)
    return float(np.mean(sobol_hqc)) if sobol_hqc else None


def plot_qband(
    qband_folder: Path,
    label: str | None = None,
    *,
    grid_root: Path | None = None,
    source_root: Path | None = None,
) -> bool:
    curves = qband_curves(
        qband_folder,
        grid_root=grid_root,
        source_root=source_root,
    )
    if not curves:
        print(f"[skip] no seed curves for {qband_folder}", flush=True)
        return False

    hqc_axis, mean, lower_1sigma, upper_1sigma, lower_2sigma, upper_2sigma = mean_bands_by_hqc(curves)
    label = label or qband_label(qband_folder)

    line = plt.plot(hqc_axis, mean, "o-", label=f"{label}")[0]
    color = line.get_color()
    sobol_hqc = qband_sobol_done_hqc(
        qband_folder,
        grid_root=grid_root,
        source_root=source_root,
    )
    if sobol_hqc is not None:
        plt.axvline(
            sobol_hqc,
            color=color,
            linestyle="--",
            linewidth=1.0,
        )
    # plt.fill_between(
    #     hqc_axis,
    #     lower_2sigma,
    #     upper_2sigma,
    #     color=color,
    #     alpha=0.08,
    #     linewidth=0,
    #     label=f"{label} 2 sigma",
    # )
    plt.fill_between(
        hqc_axis,
        lower_1sigma,
        upper_1sigma,
        color=color,
        alpha=0.16,
        linewidth=0,
        label=f"{label} 1 sigma",
    )

    print(
        f"[qband] {qband_folder.name}: seeds={len(curves)}, "
        f"hqc_range=({hqc_axis[0]:.6g}, {hqc_axis[-1]:.6g})",
        flush=True,
    )
    print(
        f"[summary] {qband_folder.name}: "
        f"final_mean_score={mean[-1]:.6g}, "
        f"final_1sigma=({lower_1sigma[-1]:.6g}, {upper_1sigma[-1]:.6g}), "
        f"final_2sigma=({lower_2sigma[-1]:.6g}, {upper_2sigma[-1]:.6g}), "
        f"final_hqc={hqc_axis[-1]:.6g}",
        flush=True,
    )
    return True


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument(
        "--grid-root",
        type=Path,
        default=None,
        help=(
            "Optional mirrored grid root. JSONs are read from the input paths, "
            "but grids are read from this root using the JSON path relative to "
            "--source-root."
        ),
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        default=None,
        help=(
            "Source root used to mirror JSON paths under --grid-root. "
            "Defaults to the parent of a uniform_10_5000_grids folder."
        ),
    )
    return parser.parse_args(argv)


def main(paths: list[str | Path] | None = None) -> None:
    args = parse_args([str(path) for path in paths] if paths is not None else sys.argv[1:])
    input_paths = [ROOT_FOLDER] if not args.paths else args.paths
    grid_root = args.grid_root
    source_root = args.source_root or infer_source_root(input_paths, grid_root)

    if grid_root is not None:
        print(f"[grid-root] source_root={source_root}", flush=True)
        print(f"[grid-root] grid_root={grid_root}", flush=True)

    qband_folders: list[Path] = []
    run_folders: list[Path] = []

    for path in input_paths:
        if path.is_dir() and measurement_pairs(
            path,
            grid_root=grid_root,
            source_root=source_root,
        ):
            run_folders.append(path)
        else:
            qband_folders.extend(qband_folder_from_arg(path))

    qband_folders = sorted(set(qband_folders), key=lambda path: path.name)
    run_folders = sorted(set(run_folders), key=lambda path: str(path))
    plotted = False

    for run_folder in run_folders:
        plotted = plot_run_folder(
            run_folder,
            label=None,
            grid_root=grid_root,
            source_root=source_root,
        ) or plotted

    for qband_folder in qband_folders:
        plotted = plot_qband(
            qband_folder,
            label=None,
            grid_root=grid_root,
            source_root=source_root,
        ) or plotted

    if not plotted:
        return

    plt.xlabel("HQC spent")
    plt.ylabel("Success-side volume")
    plt.title("FLE")
    plt.legend(fontsize=8)
    plt.tight_layout()
    SAVE_FIG_PATH.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(SAVE_FIG_PATH, dpi=300, bbox_inches="tight")
    print(f"[saved] {SAVE_FIG_PATH}", flush=True)
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
