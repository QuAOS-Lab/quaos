"""Print average/std of modified-crossing scores and plot score summaries."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Change this to the benchmark output you want to summarize.
RESULTS_JSON = Path(
    r"C:\Users\sb1580\OneDrive - University of Exeter\SimpleQ\src\sympleq\applications\randomized_benchmarking\experiments\figs\rick_modified_scores\20260619_192426\results.json"
)

OUT_HIST_PNG = Path(
    "src/sympleq/applications/randomized_benchmarking/experiments/"
    "figs/rick_modified_scores/score_summary_histogram.png"
)
OUT_NOISE_PLANE_PNG = Path(
    "src/sympleq/applications/randomized_benchmarking/experiments/"
    "figs/rick_modified_scores/score_noise_plane.png"
)
OUT_S1_S2_PLANE_PNG = Path(
    "src/sympleq/applications/randomized_benchmarking/experiments/"
    "figs/rick_modified_scores/score_s1_s2_plane.png"
)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[4]
    results_path = RESULTS_JSON
    if not results_path.is_absolute():
        results_path = repo_root / results_path

    payload = json.loads(results_path.read_text(encoding="utf-8"))
    runs = payload["runs"]
    s1 = np.asarray([float(row["S1"]) for row in runs], dtype=float)
    s2 = np.asarray([float(row["S2"]) for row in runs], dtype=float)
    one_q = np.asarray([float(row["one_q_noise_scale"]) for row in runs], dtype=float)
    two_q = np.asarray([float(row["two_q_noise_scale"]) for row in runs], dtype=float)

    print(f"runs: {len(runs)}")
    print(f"S1 mean: {np.mean(s1):.6g}")
    print(f"S1 std:  {np.std(s1, ddof=1) if len(s1) > 1 else 0.0:.6g}")
    print(f"S2 mean: {np.mean(s2):.6g}")
    print(f"S2 std:  {np.std(s2, ddof=1) if len(s2) > 1 else 0.0:.6g}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axes[0].hist(s1, bins="auto", color="tab:blue", alpha=0.8, edgecolor="black")
    axes[0].axvline(np.mean(s1), color="black", linestyle="--", linewidth=2, label="mean")
    axes[0].set_title("S1")
    axes[0].set_xlabel("Score")
    axes[0].set_ylabel("Runs")
    axes[0].legend(frameon=True)

    axes[1].hist(s2, bins="auto", color="tab:green", alpha=0.8, edgecolor="black")
    axes[1].axvline(np.mean(s2), color="black", linestyle="--", linewidth=2, label="mean")
    axes[1].set_title("S2")
    axes[1].set_xlabel("Score")
    axes[1].set_ylabel("Runs")
    axes[1].legend(frameon=True)

    hist_path = OUT_HIST_PNG
    if not hist_path.is_absolute():
        hist_path = repo_root / hist_path
    hist_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(hist_path, dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for ax, values, title in [
        (axes[0], s1, "S1 on noise plane"),
        (axes[1], s2, "S2 on noise plane"),
    ]:
        scatter = ax.scatter(
            one_q,
            two_q,
            c=values,
            cmap="viridis",
            s=90,
            edgecolors="black",
            linewidths=0.6,
        )
        ax.axvline(1.0, color="0.45", linestyle="--", linewidth=1)
        ax.axhline(1.0, color="0.45", linestyle="--", linewidth=1)
        ax.set_xlabel("1Q noise multiplier")
        ax.set_ylabel("2Q noise multiplier")
        ax.set_title(title)
        fig.colorbar(scatter, ax=ax, label="Score")

    noise_path = OUT_NOISE_PLANE_PNG
    if not noise_path.is_absolute():
        noise_path = repo_root / noise_path
    noise_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(noise_path, dpi=180)
    plt.close(fig)

    noise_multiplier = 0.5 * (one_q + two_q)

    fig, ax = plt.subplots(figsize=(6.4, 5.2), constrained_layout=True)
    scatter = ax.scatter(
        s1,
        s2,
        c=noise_multiplier,
        s=110,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.6,
        alpha=0.9,
    )
    ax.set_xlabel("S1")
    ax.set_ylabel("S2")
    ax.set_xlim(0.0, 0.1)
    fig.colorbar(scatter, ax=ax, label="Mean noise multiplier")

    score_plane_path = OUT_S1_S2_PLANE_PNG
    if not score_plane_path.is_absolute():
        score_plane_path = repo_root / score_plane_path
    score_plane_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(score_plane_path, dpi=180)
    plt.close(fig)

    print(f"histogram: {hist_path}")
    print(f"noise plane: {noise_path}")
    print(f"S1-S2 plane: {score_plane_path}")


if __name__ == "__main__":
    main()
