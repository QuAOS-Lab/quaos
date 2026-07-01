"""Minimal 3D FLE grid volume score."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


GRID_PATH = Path(
    r"Personal\seed_2025\rick_fantasy_gpu_crossing_20260623_184228_gp_grid_3d.npz"
)


def axis_spacing(values: np.ndarray, *, log: bool = False) -> float:
    axis = np.unique(values)
    if log:
        axis = np.log10(axis)
    if len(axis) < 2:
        return 1.0
    return float(np.mean(np.diff(np.sort(axis))))


def score_grid(grid_path: str | Path = GRID_PATH) -> dict[str, float]:
    grid = np.load(Path(grid_path))

    probabilities = np.asarray(grid["probabilities"], dtype=float)
    target = float(np.asarray(grid["target"]).item())

    voxel_volume = (
        axis_spacing(np.asarray(grid["gates_grid"], dtype=float), log=True)
        * axis_spacing(np.asarray(grid["ratio_grid"], dtype=float))
        * axis_spacing(np.asarray(grid["qubits_grid"], dtype=float))
    )

    valid = np.isfinite(probabilities)
    success_volume = float(np.sum(valid & (probabilities >= target)) * voxel_volume)
    failure_volume = float(np.sum(valid & (probabilities < target)) * voxel_volume)

    return {
        "target": target,
        "voxel_volume": voxel_volume,
        "success_side_volume": success_volume,
        "failure_side_volume": failure_volume,
        "total_grid_volume": success_volume + failure_volume,
    }


if __name__ == "__main__":
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else GRID_PATH
    paths = sorted(path.glob("*.npz")) if path.is_dir() else [path]
    for grid_path in paths:
        print(grid_path)
        for key, value in score_grid(grid_path).items():
            print(f"  {key}: {value:.8g}")
