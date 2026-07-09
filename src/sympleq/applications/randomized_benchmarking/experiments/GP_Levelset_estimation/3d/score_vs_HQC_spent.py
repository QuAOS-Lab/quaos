"""Plot mean 3D success-side volume by actual HQC spent for qband runs."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from score3d import score_grid


ROOT_FOLDER = Path(r"Personal\FLE")
QBAND_VALUES = range(5)
SEED_VALUES = range(2026, 2027)
SAVE_FIG_PATH = Path(
    r"Personal\RMB_results_figs\Volume\score_vs_hqc_spent_H2.pdf"
)
LABELS: list[str] = ['H2-2', 'H2-1']


def sibling_grid(json_path: Path) -> Path:
    return json_path.parent / f"{json_path.stem}_gp_grid_3d.npz"


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


def measurement_pairs(run_folder: Path) -> list[tuple[Path, Path]]:
    pairs = []
    for json_path in run_folder.glob("measurement_*.json"):
        grid_path = sibling_grid(json_path)
        if grid_path.exists():
            pairs.append((json_path, grid_path))
    return sorted(pairs, key=lambda pair: plot_order_key(pair[0]))


def sobol_done_hqc(run_folder: Path) -> float | None:
    sobol_hqc = [
        hqc_spent(json_path)
        for json_path, _ in measurement_pairs(run_folder)
        if sobol_step(json_path) is not None
    ]
    return sobol_hqc[-1] if sobol_hqc else None


def latest_seed_run(seed_folder: Path) -> Path | None:
    candidates = [
        run_folder
        for run_folder in seed_folder.iterdir()
        if run_folder.is_dir() and measurement_pairs(run_folder)
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def run_curve(run_folder: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    xs: list[int] = []
    hqc: list[float] = []
    scores: list[float] = []

    pairs = measurement_pairs(run_folder)
    if not pairs:
        print(f"[skip] no measurement/grid pairs in {run_folder}", flush=True)
        return None

    for submission_index, (json_path, grid_path) in enumerate(pairs, start=1):
        xs.append(submission_index)
        hqc.append(hqc_spent(json_path))
        score = score_grid(grid_path)["success_side_volume"]
        scores.append(score)
        print(
            f"  {json_path.name}: stitched={submission_index}, "
            f"hqc={hqc[-1]:.6g}, score={score:.6g}",
            flush=True,
        )

    return (
        np.asarray(xs, dtype=float),
        np.asarray(hqc, dtype=float),
        np.asarray(scores, dtype=float),
    )


def seed_curve(seed_folder: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    run_folder = latest_seed_run(seed_folder)
    if run_folder is None:
        print(f"[skip] no usable run folder in {seed_folder}", flush=True)
        return None

    print(f"[seed] {seed_folder.name}: {run_folder}", flush=True)
    return run_curve(run_folder)


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


def plot_run_folder(run_folder: Path, label: str | None = None) -> bool:
    curve = run_curve(run_folder)
    if curve is None:
        return False

    _, hqc, scores = curve
    line = plt.plot(hqc, scores, "o-", label=label or run_label(run_folder))[0]
    sobol_hqc = sobol_done_hqc(run_folder)
    if sobol_hqc is not None:
        plt.axvline(
            sobol_hqc,
            color=line.get_color(),
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


def qband_curves(qband_folder: Path) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    curves = []
    for seed in SEED_VALUES:
        seed_folder = qband_folder / f"seed_{seed}"
        if not seed_folder.exists():
            print(f"[skip] missing {seed_folder}", flush=True)
            continue
        curve = seed_curve(seed_folder)
        if curve is not None:
            curves.append(curve)
    if not curves:
        return curves

    max_hqc = max(float(np.max(hqc)) for _, hqc, _ in curves)
    min_allowed_hqc = max_hqc - 20.0
    filtered_curves = []
    for curve in curves:
        _, hqc, _ = curve
        final_hqc = float(np.max(hqc))
        if final_hqc < min_allowed_hqc:
            print(
                f"[exclude] {qband_folder.name}: seed final_hqc={final_hqc:.6g} "
                f"is below max_hqc-20={min_allowed_hqc:.6g}",
                flush=True,
            )
            continue
        filtered_curves.append(curve)
    return filtered_curves


def mean_std_by_hqc(
    curves: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    max_len = max(len(scores) for _, _, scores in curves)
    hqc_min = max(float(np.min(hqc)) for _, hqc, _ in curves)
    hqc_max = min(float(np.max(hqc)) for _, hqc, _ in curves)
    hqc_axis = np.linspace(hqc_min, hqc_max, max_len)
    score_matrix = np.full((len(curves), max_len), np.nan, dtype=float)

    for row, (_, hqc, scores) in enumerate(curves):
        order = np.argsort(hqc)
        unique_hqc, unique_indices = np.unique(hqc[order], return_index=True)
        unique_scores = scores[order][unique_indices]
        score_matrix[row] = np.interp(hqc_axis, unique_hqc, unique_scores)

    mean = np.nanmean(score_matrix, axis=0)
    std = np.nanstd(score_matrix, axis=0, ddof=1)
    std[np.isnan(std)] = 0.0
    return hqc_axis, mean, std


def qband_sobol_done_hqc(qband_folder: Path) -> float | None:
    sobol_hqc = []
    for seed in SEED_VALUES:
        run_folder = latest_seed_run(qband_folder / f"seed_{seed}")
        if run_folder is None:
            continue
        value = sobol_done_hqc(run_folder)
        if value is not None:
            sobol_hqc.append(value)
    return float(np.mean(sobol_hqc)) if sobol_hqc else None


def plot_qband(qband_folder: Path, label: str | None = None) -> bool:
    curves = qband_curves(qband_folder)
    if not curves:
        print(f"[skip] no seed curves for {qband_folder}", flush=True)
        return False

    hqc_axis, mean, std = mean_std_by_hqc(curves)
    label = label or qband_label(qband_folder)

    line = plt.plot(hqc_axis, mean, "o-", label=f"{label}")[0]
    color = line.get_color()
    sobol_hqc = qband_sobol_done_hqc(qband_folder)
    if sobol_hqc is not None:
        plt.axvline(
            sobol_hqc,
            color=color,
            linestyle="--",
            linewidth=1.0,
        )
    plt.fill_between(
        hqc_axis,
        mean - std,
        mean + std,
        color=color,
        alpha=0.18,
        linewidth=0,
    )

    print(
        f"[qband] {qband_folder.name}: seeds={len(curves)}, "
        f"hqc_range=({hqc_axis[0]:.6g}, {hqc_axis[-1]:.6g})",
        flush=True,
    )
    print(
        f"[summary] {qband_folder.name}: "
        f"final_mean_score={mean[-1]:.6g}, "
        f"final_std={std[-1]:.6g}, "
        f"final_hqc={hqc_axis[-1]:.6g}",
        flush=True,
    )
    return True


def main(paths: list[str | Path] | None = None) -> None:
    input_paths = [ROOT_FOLDER] if not paths else [Path(path) for path in paths]
    qband_folders: list[Path] = []
    run_folders: list[Path] = []

    for path in input_paths:
        if path.is_dir() and measurement_pairs(path):
            run_folders.append(path)
        else:
            qband_folders.extend(qband_folder_from_arg(path))

    qband_folders = sorted(set(qband_folders), key=lambda path: path.name)
    run_folders = sorted(set(run_folders), key=lambda path: str(path))
    plotted = False
    label_index = 0

    for run_folder in run_folders:
        label = LABELS[label_index] if label_index < len(LABELS) else None
        plotted = plot_run_folder(run_folder, label=label) or plotted
        label_index += 1

    for qband_folder in qband_folders:
        label = LABELS[label_index] if label_index < len(LABELS) else None
        plotted = plot_qband(qband_folder, label=label) or plotted
        label_index += 1

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
