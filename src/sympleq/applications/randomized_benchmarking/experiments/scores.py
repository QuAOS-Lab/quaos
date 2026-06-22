"""Small score helpers for modified-crossing GP grids."""
from __future__ import annotations

from pathlib import Path

import numpy as np


REFERENCE_OFFSET = 1.91e-4
REFERENCE_SLOPE = 3.65e-3
SIGMA_FLOOR = 1e-8


def true_log_gates(
    ratio: np.ndarray,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> np.ndarray:
    """Analytic Lindblad contour in log(total gates)."""
    gates = np.log(2.0) / (
        REFERENCE_OFFSET * one_q_noise_scale
        + REFERENCE_SLOPE * two_q_noise_scale * ratio
    )
    return np.log(gates)


def gp_mean_contour_from_grid(grid_path: str | Path) -> dict[str, np.ndarray]:
    """Extract the GP mean latent_mean=0 contour and propagated uncertainty."""
    data = np.load(grid_path)
    gates = np.asarray(data["gates_grid"], dtype=float)
    ratio = np.asarray(data["ratio_grid"], dtype=float)
    latent_mean = np.asarray(data["latent_mean"], dtype=float)
    latent_variance = np.maximum(np.asarray(data["latent_variance"], dtype=float), 0.0)
    latent_std = np.sqrt(latent_variance)

    log_gates_axis = np.log(gates[0, :])
    ratio_axis = ratio[:, 0]
    dmu_dratio, dmu_dlogg = np.gradient(
        latent_mean,
        ratio_axis,
        log_gates_axis,
        edge_order=1,
    )
    grad_norm = np.sqrt(dmu_dratio**2 + dmu_dlogg**2)
    sigma_surface = latent_std / np.maximum(grad_norm, SIGMA_FLOOR)

    contour_ratio: list[float] = []
    contour_log_gates: list[float] = []
    contour_sigma: list[float] = []
    for row, r in enumerate(ratio_axis):
        values = latent_mean[row, :]
        crossings = np.where(values[:-1] * values[1:] <= 0.0)[0]
        if len(crossings) == 0:
            continue
        col = int(crossings[0])
        denom = abs(values[col]) + abs(values[col + 1])
        t = 0.0 if denom <= 0.0 else abs(values[col]) / denom
        log_g = (1.0 - t) * log_gates_axis[col] + t * log_gates_axis[col + 1]
        sigma = (1.0 - t) * sigma_surface[row, col] + t * sigma_surface[row, col + 1]
        if np.isfinite(log_g) and np.isfinite(sigma) and sigma > 0.0:
            contour_ratio.append(float(r))
            contour_log_gates.append(float(log_g))
            contour_sigma.append(float(max(sigma, SIGMA_FLOOR)))

    return {
        "ratio": np.asarray(contour_ratio, dtype=float),
        "gp_log_gates": np.asarray(contour_log_gates, dtype=float),
        "sigma_contour": np.asarray(contour_sigma, dtype=float),
    }


def gp_grid_scores(
    grid_path: str | Path,
    *,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> dict[str, float]:
    """Return S1 and S2 for a saved GP grid NPZ."""
    contour = gp_mean_contour_from_grid(grid_path)
    ratio = contour["ratio"]
    gp_log_gates = contour["gp_log_gates"]
    sigma = contour["sigma_contour"]
    if len(ratio) < 2:
        raise ValueError(f"Not enough GP contour points in {grid_path}")

    analytic_log_gates = true_log_gates(
        ratio,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    delta = np.abs(gp_log_gates - analytic_log_gates)
    a_gp = float(np.trapz(gp_log_gates, ratio))

    return {
        "S1": float(np.trapz(delta, ratio) / a_gp),
        "S2": float(np.trapz(sigma, ratio) / a_gp),
        "A_gp": a_gp,
        "mean_delta_log_gates": float(np.mean(delta)),
        "mean_sigma_contour": float(np.mean(sigma)),
        "n_contour_points": int(len(ratio)),
    }
