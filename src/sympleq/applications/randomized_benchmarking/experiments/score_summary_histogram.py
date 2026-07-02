"""Print average/std of contour scores and plot score summaries."""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

if os.environ.get("SCORE_SUMMARY_USE_AGG", "0") == "1":
    matplotlib.use("Agg", force=False)
import matplotlib.pyplot as plt
import numpy as np


RESULTS_JSON = Path(
    os.environ.get(
        "SCORE_SUMMARY_RESULTS_JSON",
        "src/sympleq/applications/randomized_benchmarking/experiments/"
        "figs/benchmarkbenchmark_scores/results.json",
    )
)


def _resolve_repo_path(path: Path) -> Path:
    repo_root = Path(__file__).resolve().parents[4]
    if path.is_absolute():
        return path
    return repo_root / path


def _finite_array(rows: list[dict], key: str) -> np.ndarray:
    values = []
    for row in rows:
        value = row.get(key)
        if value is None:
            continue
        value = float(value)
        if np.isfinite(value):
            values.append(value)
    return np.asarray(values, dtype=float)


def _finite_rows(rows: list[dict], keys: tuple[str, ...]) -> dict[str, np.ndarray]:
    values = {key: [] for key in keys}
    for row in rows:
        parsed = {}
        for key in keys:
            value = row.get(key)
            if value is None:
                break
            value = float(value)
            if not np.isfinite(value):
                break
            parsed[key] = value
        else:
            for key, value in parsed.items():
                values[key].append(value)
    return {key: np.asarray(items, dtype=float) for key, items in values.items()}


def _approach(row: dict) -> str:
    return str(row.get("approach", "method"))


def _print_score_summary(runs: list[dict], approaches: list[str]) -> None:
    print(f"runs: {len(runs)}")
    for approach in approaches:
        rows = [row for row in runs if _approach(row) == approach]
        print(f"  {approach}")
        for key in ("S1", "S2", "S_total"):
            values = _finite_array(rows, key)
            if len(values) == 0:
                print(f"    {key}: unavailable")
                continue
            std = np.std(values, ddof=1) if len(values) > 1 else 0.0
            print(
                f"    {key} mean: {np.mean(values):.6g}, "
                f"std: {std:.6g}, best: {np.min(values):.6g}, "
                f"worst: {np.max(values):.6g}"
            )


def _plot_histograms(
    runs: list[dict],
    approaches: list[str],
    path: Path,
    *,
    show: bool = False,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    colors = plt.get_cmap("tab10")
    for index, approach in enumerate(approaches):
        rows = [row for row in runs if _approach(row) == approach]
        color = colors(index % 10)
        for ax, key in zip(axes, ("S1", "S2")):
            values = _finite_array(rows, key)
            if len(values) == 0:
                continue
            label = f"{approach} mean={np.mean(values):.3g}"
            ax.hist(
                values,
                bins="auto",
                alpha=0.48,
                color=color,
                edgecolor="black",
                linewidth=0.5,
                label=label,
            )
            ax.axvline(np.mean(values), color=color, linestyle="--", linewidth=1.6)

    for ax, title in zip(axes, ("S1", "S2")):
        ax.set_title(title)
        ax.set_xlabel("Score")
        ax.set_ylabel("Runs")
        ax.legend(frameon=True, fontsize=8)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    if show:
        fig.show()
    else:
        plt.close(fig)


def _plot_noise_planes(
    runs: list[dict],
    approaches: list[str],
    path: Path,
    *,
    show: bool = False,
) -> None:
    fig, axes = plt.subplots(
        len(approaches),
        2,
        figsize=(11, max(4.2, 3.8 * len(approaches))),
        constrained_layout=True,
        squeeze=False,
    )
    for row_index, approach in enumerate(approaches):
        rows = [row for row in runs if _approach(row) == approach]
        for col_index, key in enumerate(("S1", "S2")):
            ax = axes[row_index, col_index]
            values = _finite_rows(rows, ("one_q_noise_scale", "two_q_noise_scale", key))
            if len(values[key]) == 0:
                ax.set_axis_off()
                continue
            scatter = ax.scatter(
                values["one_q_noise_scale"],
                values["two_q_noise_scale"],
                c=values[key],
                cmap="viridis",
                s=75,
                edgecolors="black",
                linewidths=0.6,
            )
            ax.axvline(1.0, color="0.45", linestyle="--", linewidth=1)
            ax.axhline(1.0, color="0.45", linestyle="--", linewidth=1)
            ax.set_xlabel("1Q noise multiplier")
            ax.set_ylabel("2Q noise multiplier")
            ax.set_title(f"{approach}: {key} on noise plane")
            fig.colorbar(scatter, ax=ax, label=key)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    if show:
        fig.show()
    else:
        plt.close(fig)


def _plot_s1_s2_plane(
    runs: list[dict],
    approaches: list[str],
    path: Path,
    *,
    show: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 5.4), constrained_layout=True)
    colors = plt.get_cmap("tab10")
    all_s1 = []
    all_s2 = []
    for index, approach in enumerate(approaches):
        rows = [row for row in runs if _approach(row) == approach]
        values = _finite_rows(rows, ("S1", "S2", "one_q_noise_scale", "two_q_noise_scale"))
        if len(values["S1"]) == 0:
            continue
        noise_multiplier = 0.5 * (
            values["one_q_noise_scale"] + values["two_q_noise_scale"]
        )
        scatter = ax.scatter(
            values["S1"],
            values["S2"],
            c=noise_multiplier,
            s=90,
            cmap="viridis",
            marker=["o", "s", "^", "D", "P", "X"][index % 6],
            edgecolors=colors(index % 10),
            linewidths=1.0,
            alpha=0.9,
            label=approach,
        )
        all_s1.extend(values["S1"])
        all_s2.extend(values["S2"])
    ax.set_xlabel("S1")
    ax.set_ylabel("S2")
    ax.set_title("S1-S2 score plane")
    ax.legend(frameon=True, fontsize=8)
    if all_s1:
        ax.set_xlim(left=0.0, right=max(all_s1) * 1.08 if max(all_s1) > 0 else 1.0)
    if all_s2:
        ax.set_ylim(bottom=0.0, top=max(all_s2) * 1.08 if max(all_s2) > 0 else 1.0)
    if all_s1 and all_s2:
        fig.colorbar(scatter, ax=ax, label="Mean noise multiplier")

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    if show:
        fig.show()
    else:
        plt.close(fig)


def plot_score_summaries(
    results_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    show: bool = False,
) -> dict[str, Path]:
    """Create S1/S2 summary plots for a benchmark ``results.json`` file."""
    results_path = _resolve_repo_path(Path(results_path))

    payload = json.loads(results_path.read_text(encoding="utf-8"))
    runs = payload["runs"]
    approaches = list(dict.fromkeys(_approach(row) for row in runs))
    if output_dir is None:
        output_dir = results_path.parent
    output_dir = _resolve_repo_path(Path(output_dir))

    _print_score_summary(runs, approaches)

    hist_path = output_dir / "score_summary_histogram.png"
    noise_path = output_dir / "score_noise_plane.png"
    score_plane_path = output_dir / "score_s1_s2_plane.png"
    _plot_histograms(runs, approaches, hist_path, show=show)
    _plot_noise_planes(runs, approaches, noise_path, show=show)
    _plot_s1_s2_plane(runs, approaches, score_plane_path, show=show)
    if show:
        plt.show()

    print(f"histogram: {hist_path}")
    print(f"noise plane: {noise_path}")
    print(f"S1-S2 plane: {score_plane_path}")
    return {
        "histogram": hist_path,
        "noise_plane": noise_path,
        "s1_s2_plane": score_plane_path,
    }


def main() -> None:
    plot_score_summaries(RESULTS_JSON, show=True)


if __name__ == "__main__":
    main()
