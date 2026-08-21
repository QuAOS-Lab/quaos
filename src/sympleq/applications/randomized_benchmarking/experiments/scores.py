"""Score helpers for randomized-benchmarking boundary fits."""
from __future__ import annotations

from pathlib import Path

import numpy as np


REFERENCE_OFFSET = 1.91e-4
REFERENCE_SLOPE = 3.65e-3
SIGMA_FLOOR = 1e-8


def axis_spacing(values: np.ndarray, *, log: bool = False) -> float:
    axis = np.unique(np.asarray(values, dtype=float))
    if log:
        axis = np.log10(axis)
    if len(axis) < 2:
        return 1.0
    return float(np.mean(np.diff(np.sort(axis))))


def integrate_trapezoid(
    y: np.ndarray,
    x: np.ndarray,
    *,
    axis: int = -1,
) -> np.ndarray:
    """NumPy-version compatible trapezoid integration."""
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is not None:
        return trapezoid(y, x, axis=axis)
    return np.trapz(y, x, axis=axis)


def surface_average(
    values: np.ndarray,
    ratio_axis: np.ndarray,
    qubit_axis: np.ndarray,
    *,
    eps: float = 1e-12,
) -> float:
    """Average a surface over ratio and qubit axes with trapezoid integration."""
    values = np.nan_to_num(np.asarray(values, dtype=float))
    ratio_axis = np.asarray(ratio_axis, dtype=float)
    qubit_axis = np.asarray(qubit_axis, dtype=float)
    r_span = max(float(np.nanmax(ratio_axis) - np.nanmin(ratio_axis)), eps)
    if len(qubit_axis) == 1:
        return float(integrate_trapezoid(values[0], ratio_axis) / r_span)
    q_span = max(float(np.nanmax(qubit_axis) - np.nanmin(qubit_axis)), eps)
    return float(
        integrate_trapezoid(
            integrate_trapezoid(values, ratio_axis, axis=1),
            qubit_axis,
        )
        / max(r_span * q_span, eps)
    )


def surface_log_scores(
    fitted_log_gates: np.ndarray,
    reference_log_gates: np.ndarray,
    sigma_log_gates: np.ndarray,
    ratio_axis: np.ndarray,
    qubit_axis: np.ndarray,
    *,
    eps: float = 1e-12,
) -> dict[str, float | int]:
    """Return S1/S2 scores for fitted/reference log-gate surfaces."""
    fitted_log_gates = np.asarray(fitted_log_gates, dtype=float)
    reference_log_gates = np.asarray(reference_log_gates, dtype=float)
    sigma_log_gates = np.asarray(sigma_log_gates, dtype=float)
    valid = (
        np.isfinite(fitted_log_gates)
        & np.isfinite(reference_log_gates)
        & np.isfinite(sigma_log_gates)
    )
    n_points = int(np.count_nonzero(valid))
    if n_points < 2:
        return {
            "S1": float("nan"),
            "S2": float("nan"),
            "mean_delta_log_gates": float("nan"),
            "mean_sigma_log_gates": float("nan"),
            "score_points": n_points,
        }

    fitted = np.where(valid, fitted_log_gates, np.nan)
    reference = np.where(valid, reference_log_gates, np.nan)
    sigma = np.where(valid, sigma_log_gates, np.nan)
    delta = np.where(valid, np.abs(fitted - reference), np.nan)

    normaliser = abs(surface_average(fitted, ratio_axis, qubit_axis, eps=eps))
    if not np.isfinite(normaliser) or normaliser <= eps:
        normaliser = float(np.nanmean(np.abs(fitted)))
    normaliser = max(float(normaliser), eps)
    return {
        "S1": float(
            surface_average(delta, ratio_axis, qubit_axis, eps=eps) / normaliser
        ),
        "S2": float(
            surface_average(sigma, ratio_axis, qubit_axis, eps=eps) / normaliser
        ),
        "mean_delta_log_gates": float(np.nanmean(delta)),
        "mean_sigma_log_gates": float(np.nanmean(sigma)),
        "score_points": n_points,
    }


def success_side_log_volume(
    boundary_gates: np.ndarray,
    ratio_axis: np.ndarray,
    qubit_axis: np.ndarray,
    min_gates: float,
    max_gates: float,
    *,
    eps: float = 1e-12,
) -> float:
    """Log10-gate success-side volume under a boundary surface."""
    boundary_gates = np.asarray(boundary_gates, dtype=float)
    ratio_axis = np.asarray(ratio_axis, dtype=float)
    qubit_axis = np.asarray(qubit_axis, dtype=float)
    valid_boundary = np.isfinite(boundary_gates) & (boundary_gates > 0.0)
    if not np.any(valid_boundary):
        return float("nan")
    min_gates = max(float(min_gates), eps)
    max_gates = max(float(max_gates), min_gates)
    clipped = np.clip(boundary_gates, min_gates, max_gates)
    height = np.where(
        valid_boundary,
        np.maximum(np.log10(clipped) - np.log10(min_gates), 0.0),
        np.nan,
    )
    height = np.nan_to_num(height)
    if len(qubit_axis) == 1:
        return float(integrate_trapezoid(height[0], ratio_axis))
    return float(
        integrate_trapezoid(
            integrate_trapezoid(height, ratio_axis, axis=1),
            qubit_axis,
        )
    )


def surface_boundary_scores_2d(
    *,
    fitted_log_gates: np.ndarray,
    reference_log_gates: np.ndarray,
    sigma_log_gates: np.ndarray,
    ratio_axis: np.ndarray,
    qubit_axis: np.ndarray,
    volume_fit_gates: np.ndarray,
    volume_sigma_log_gates: np.ndarray,
    volume_reference_gates: np.ndarray,
    volume_ratio_axis: np.ndarray,
    volume_qubit_axis: np.ndarray,
    min_gates: float,
    max_gates: float,
    eps: float = 1e-12,
) -> dict[str, float | int]:
    """S1/S2 and success-side volume scores for 2-D boundary surfaces."""
    scores = surface_log_scores(
        fitted_log_gates,
        reference_log_gates,
        sigma_log_gates,
        ratio_axis,
        qubit_axis,
        eps=eps,
    )
    if int(scores["score_points"]) < 2:
        return {
            "surface_S1": float("nan"),
            "surface_S2": float("nan"),
            "surface_volume_fit": float("nan"),
            "surface_volume_lower_1sigma": float("nan"),
            "surface_volume_upper_1sigma": float("nan"),
            "surface_volume_reference": float("nan"),
            "surface_volume_ratio": float("nan"),
            "surface_mean_delta_log_gates": float("nan"),
            "surface_mean_sigma_log_gates": float("nan"),
            "surface_score_points": int(scores["score_points"]),
        }

    volume_fit_gates = np.asarray(volume_fit_gates, dtype=float)
    volume_sigma_log_gates = np.asarray(volume_sigma_log_gates, dtype=float)
    volume_reference_gates = np.asarray(volume_reference_gates, dtype=float)
    valid_volume = (
        np.isfinite(volume_fit_gates)
        & np.isfinite(volume_reference_gates)
        & np.isfinite(volume_sigma_log_gates)
        & (volume_fit_gates > 0.0)
        & (volume_reference_gates > 0.0)
    )
    volume_fit_gates = np.where(valid_volume, volume_fit_gates, np.nan)
    volume_sigma_log_gates = np.where(valid_volume, volume_sigma_log_gates, np.nan)
    volume_reference_gates = np.where(valid_volume, volume_reference_gates, np.nan)
    volume_fit = success_side_log_volume(
        volume_fit_gates,
        volume_ratio_axis,
        volume_qubit_axis,
        min_gates,
        max_gates,
        eps=eps,
    )
    volume_lower = success_side_log_volume(
        volume_fit_gates * np.exp(-volume_sigma_log_gates),
        volume_ratio_axis,
        volume_qubit_axis,
        min_gates,
        max_gates,
        eps=eps,
    )
    volume_upper = success_side_log_volume(
        volume_fit_gates * np.exp(volume_sigma_log_gates),
        volume_ratio_axis,
        volume_qubit_axis,
        min_gates,
        max_gates,
        eps=eps,
    )
    volume_reference = success_side_log_volume(
        volume_reference_gates,
        volume_ratio_axis,
        volume_qubit_axis,
        min_gates,
        max_gates,
        eps=eps,
    )
    volume_ratio = (
        float(volume_fit / volume_reference)
        if (
            np.isfinite(volume_fit)
            and np.isfinite(volume_reference)
            and abs(volume_reference) > eps
        )
        else float("nan")
    )
    return {
        "surface_S1": float(scores["S1"]),
        "surface_S2": float(scores["S2"]),
        "surface_volume_fit": float(volume_fit),
        "surface_volume_lower_1sigma": float(volume_lower),
        "surface_volume_upper_1sigma": float(volume_upper),
        "surface_volume_reference": float(volume_reference),
        "surface_volume_ratio": volume_ratio,
        "surface_mean_delta_log_gates": float(scores["mean_delta_log_gates"]),
        "surface_mean_sigma_log_gates": float(scores["mean_sigma_log_gates"]),
        "surface_score_points": int(scores["score_points"]),
    }


def surface_log_ratio_rows(
    fitted_gates: np.ndarray,
    reference_gates: np.ndarray,
    qubit_axis: np.ndarray,
) -> list[tuple[int, list[float]]]:
    """Rows of log(fitted/reference) values for printing or saving."""
    fitted_gates = np.asarray(fitted_gates, dtype=float)
    reference_gates = np.asarray(reference_gates, dtype=float)
    rows: list[tuple[int, list[float]]] = []
    for q, fitted_row, reference_row in zip(qubit_axis, fitted_gates, reference_gates):
        valid = (
            np.isfinite(fitted_row)
            & np.isfinite(reference_row)
            & (fitted_row > 0.0)
            & (reference_row > 0.0)
        )
        row = np.full_like(fitted_row, np.nan, dtype=float)
        row[valid] = np.log(fitted_row[valid] / reference_row[valid])
        rows.append((int(q), [float(value) for value in row]))
    return rows


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
    print(grid_path)
    ratio = contour["ratio"]
    gp_log_gates = contour["gp_log_gates"]
    sigma = contour["sigma_contour"]
    if len(ratio) < 2:
        print(contour["ratio"])
        print(contour["gp_log_gates"])
        print(contour["sigma_contour"])
        raise ValueError(f"Not enough GP contour points in {grid_path}")

    analytic_log_gates = true_log_gates(
        ratio,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    delta = np.abs(gp_log_gates - analytic_log_gates)
    a_gp = float(integrate_trapezoid(gp_log_gates, ratio))

    return {
        "S1": float(integrate_trapezoid(delta, ratio) / a_gp),
        "S2": float(integrate_trapezoid(sigma, ratio) / a_gp),
        "A_gp": a_gp,
        "mean_delta_log_gates": float(np.mean(delta)),
        "mean_sigma_contour": float(np.mean(sigma)),
        "n_contour_points": int(len(ratio)),
    }
