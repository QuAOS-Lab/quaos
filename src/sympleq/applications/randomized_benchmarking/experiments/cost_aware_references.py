"""Reference-surface helpers for cost-aware randomized benchmarking runs."""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable

import numpy as np

from sympleq.applications.randomized_benchmarking.backends.exponential import (
    asymptote_value,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    analytic_gate_counts,
)
from sympleq.applications.randomized_benchmarking.experiments.scores import (
    integrate_trapezoid,
)

_EPS = 1e-12
_LN2 = np.log(2.0)


def raw_fidelity_gate_factor(
    settings: Any,
    qubits: np.ndarray,
    visibility: float | np.ndarray = 1.0,
) -> np.ndarray:
    """Factor g(Q) mapping the renormalised boundary n* onto the raw p=0.5 depth.

    The renormalised boundary n* sits where the *renormalised* fidelity
    2^{-nD} = 0.5.  The raw survival p = V(1-B)2^{-nD} + B crosses 0.5 at
    n_raw = n* * g(Q), with

        g(Q) = -log2[ (0.5 - B(Q)) / (V (1 - B(Q))) ],

    where B(Q) is the survival asymptote (``settings.asymptote_model``) and V the
    visibility.  g depends only on B(Q) and V (not on r); it is >= 1 wherever the
    crossing exists and -> 1 as B -> 0 (large Q).  Returns nan where no raw p=0.5
    depth exists (B(Q) >= 0.5, or the un-decayed top V(1-B)+B < 0.5).  This is the
    single definition shared by the surface method (fitted boundary/volume) and
    the plot module (analytic reference overlays), so raw-fidelity remaps stay
    consistent across every figure and score.
    """
    v = np.asarray(visibility, dtype=float)
    q = np.asarray(qubits, dtype=float)
    model = getattr(settings, "asymptote_model", "depolarizing")
    unique_q = np.unique(q)
    b_of_q = {float(x): float(asymptote_value(model, float(x))) for x in unique_q}
    b = np.vectorize(b_of_q.get, otypes=[float])(q)
    denom = v * (1.0 - b)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(denom > 0.0, (0.5 - b) / denom, np.nan)
        return np.where((x > 0.0) & (x < 1.0), -np.log2(x), np.nan)


def analytic_lindblad_gates(
    ratios: np.ndarray,
    qubits: np.ndarray,
    settings: Any,
) -> np.ndarray:
    """Analytic Lindblad fidelity-0.5 surface.

    The current analytic reference is independent of Q because the reference
    model is expressed in total gate count and two-qubit ratio only.
    """
    one_q = float(getattr(settings, "one_q_noise_scale", 1.0))
    two_q = float(getattr(settings, "two_q_noise_scale", 1.0))
    ratio_grid, qubit_grid = np.broadcast_arrays(ratios, qubits)
    gates = analytic_gate_counts(
        ratio_grid,
        one_q_noise_scale=one_q,
        two_q_noise_scale=two_q,
    )
    return np.full_like(qubit_grid, 1.0, dtype=float) * gates


def reference_surface_volume(
    settings: Any,
    reference_gates: Callable[[np.ndarray, np.ndarray, Any], np.ndarray] = (
        analytic_lindblad_gates
    ),
) -> float:
    """Integral of reference log10(n*(r,Q)) over the scored rectangle."""
    r_lo, r_hi = settings.ratio_bounds
    ratios = np.linspace(r_lo, r_hi, max(settings.score_ratio_points, 2))
    qubits = np.asarray(settings.q_values, dtype=float)
    rr, qq = np.meshgrid(ratios, qubits)
    values = reference_gates(rr, qq, settings)
    finite = np.isfinite(values)
    if not np.any(finite):
        return float("nan")
    # Some optional reference grids cover only part of the scored rectangle.
    # Integrate over the supported rectangular subgrid instead of silently
    # replacing unsupported cells by zero.
    if not np.all(finite):
        row_keep = np.any(finite, axis=1)
        col_keep = np.any(finite, axis=0)
        if np.count_nonzero(row_keep) == 0 or np.count_nonzero(col_keep) < 2:
            return float("nan")
        qubits = qubits[row_keep]
        ratios = ratios[col_keep]
        values = values[np.ix_(row_keep, col_keep)]
        if not np.all(np.isfinite(values)):
            return float("nan")
    if np.any(values <= 0.0):
        return float("nan")
    values = np.log10(values)
    if len(qubits) == 1:
        return float(integrate_trapezoid(values[0], ratios))
    return float(
        integrate_trapezoid(integrate_trapezoid(values, ratios, axis=1), qubits)
    )


@lru_cache(maxsize=8)
def load_calibrated_surface(path: str) -> dict:
    with open(path) as f:
        payload = json.load(f)
    for key in ("lambda1_poly_coefficients", "lambda2_poly_coefficients"):
        if key not in payload:
            raise ValueError(f"Calibrated surface {path!r} is missing {key!r}.")
    return payload


def calibrated_surface_label(settings: Any) -> str:
    if settings.calibrated_surface_path is None:
        return "calibrated"
    payload = load_calibrated_surface(str(Path(settings.calibrated_surface_path)))
    seed = payload.get("settings", {}).get("rng_seed")
    suffix = f" seed={seed}" if seed is not None else ""
    return f"{payload.get('label', 'calibrated SympleQ')}{suffix}"


def calibrated_surface_gates(
    ratios: np.ndarray,
    qubits: np.ndarray,
    settings: Any,
) -> np.ndarray:
    """Evaluate the saved calibrated effective-rate reference surface.

    The calibration file stores natural-log rates lambda_i(Q).  The boundary is

        n_ref(r,Q) = ln 2 / [(1-r) lambda_1(Q) + r lambda_2(Q)].
    """
    if settings.calibrated_surface_path is None:
        raise ValueError("No calibrated_surface_path configured.")
    payload = load_calibrated_surface(str(Path(settings.calibrated_surface_path)))
    ratio_grid, qubit_grid = np.broadcast_arrays(ratios, qubits)
    variable = payload.get("polynomial_variable", "q_minus_reference")
    if variable == "q_minus_reference":
        x = qubit_grid.astype(float) - float(payload.get("q_reference", 0.0))
    elif variable == "q":
        x = qubit_grid.astype(float)
    else:
        raise ValueError(f"Unknown calibrated polynomial variable {variable!r}.")
    lam1 = np.polyval(np.asarray(payload["lambda1_poly_coefficients"], dtype=float), x)
    lam2 = np.polyval(np.asarray(payload["lambda2_poly_coefficients"], dtype=float), x)
    denominator = (1.0 - ratio_grid) * lam1 + ratio_grid * lam2
    gates = np.full_like(ratio_grid, np.nan, dtype=float)
    ok = denominator > 0.0
    gates[ok] = _LN2 / denominator[ok]
    return gates


def gp_grid_surface_label(settings: Any) -> str:
    if settings.gp_grid_surface_label is not None:
        return settings.gp_grid_surface_label
    if settings.gp_grid_surface_path is None:
        return "external GP grid"
    return Path(settings.gp_grid_surface_path).stem


def gate_crossing_from_probability(
    gates: np.ndarray,
    probabilities: np.ndarray,
    target: float,
) -> float:
    """Linear p=target crossing along increasing gate counts."""
    gates = np.asarray(gates, dtype=float)
    probabilities = np.asarray(probabilities, dtype=float)
    ok = np.isfinite(gates) & np.isfinite(probabilities)
    if np.count_nonzero(ok) < 2:
        return float("nan")
    gates = gates[ok]
    probabilities = probabilities[ok]
    diff = probabilities - target
    exact = np.flatnonzero(np.isclose(diff, 0.0, atol=1e-12))
    if len(exact):
        return float(gates[int(exact[0])])
    sign_change = np.flatnonzero(diff[:-1] * diff[1:] < 0.0)
    if not len(sign_change):
        return float("nan")
    # Prefer the physical high-to-low survival transition if it exists.
    down = [i for i in sign_change if diff[i] > 0.0 and diff[i + 1] < 0.0]
    i = int(down[0] if down else sign_change[0])
    p0, p1 = probabilities[i], probabilities[i + 1]
    if abs(p1 - p0) <= _EPS:
        return float(0.5 * (gates[i] + gates[i + 1]))
    frac = (target - p0) / (p1 - p0)
    return float(gates[i] + frac * (gates[i + 1] - gates[i]))


@lru_cache(maxsize=8)
def load_gp_grid_surface(path: str) -> dict:
    with np.load(path, allow_pickle=True) as payload:
        qubits = np.asarray(payload["qubits_axis"], dtype=float)
        ratios = np.asarray(payload["ratio_axis"], dtype=float)
        gates = np.asarray(payload["gates_axis"], dtype=float)
        probabilities = np.asarray(payload["probabilities"], dtype=float)
        target = float(payload["target"]) if "target" in payload else 0.5

    expected = (len(qubits), len(ratios), len(gates))
    if probabilities.shape != expected:
        raise ValueError(
            f"GP grid probabilities shape {probabilities.shape} does not match "
            f"(qubits, ratios, gates) = {expected}."
        )
    contour = np.full((len(qubits), len(ratios)), np.nan, dtype=float)
    for qi in range(len(qubits)):
        for ri in range(len(ratios)):
            contour[qi, ri] = gate_crossing_from_probability(
                gates, probabilities[qi, ri, :], target
            )
    return {
        "qubits": qubits,
        "ratios": ratios,
        "gates": gates,
        "contour": contour,
        "target": target,
    }


def interp_grid2d_nan(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    values: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Bilinear interpolation on a regular grid, returning nan outside/near gaps."""
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)
    values = np.asarray(values, dtype=float)
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    out = np.full_like(x, np.nan, dtype=float)
    inside = (
        (x >= x_axis[0])
        & (x <= x_axis[-1])
        & (y >= y_axis[0])
        & (y <= y_axis[-1])
    )
    if not np.any(inside):
        return out
    xi = np.searchsorted(x_axis, x[inside], side="right") - 1
    yi = np.searchsorted(y_axis, y[inside], side="right") - 1
    xi = np.clip(xi, 0, len(x_axis) - 2)
    yi = np.clip(yi, 0, len(y_axis) - 2)
    x0, x1 = x_axis[xi], x_axis[xi + 1]
    y0, y1 = y_axis[yi], y_axis[yi + 1]
    tx = np.divide(x[inside] - x0, x1 - x0, out=np.zeros_like(x0), where=x1 != x0)
    ty = np.divide(y[inside] - y0, y1 - y0, out=np.zeros_like(y0), where=y1 != y0)
    v00 = values[xi, yi]
    v10 = values[xi + 1, yi]
    v01 = values[xi, yi + 1]
    v11 = values[xi + 1, yi + 1]
    valid = np.isfinite(v00) & np.isfinite(v10) & np.isfinite(v01) & np.isfinite(v11)
    interp = (
        (1.0 - tx) * (1.0 - ty) * v00
        + tx * (1.0 - ty) * v10
        + (1.0 - tx) * ty * v01
        + tx * ty * v11
    )
    inside_idx = np.flatnonzero(inside)
    out.flat[inside_idx[valid]] = interp[valid]
    return out


def gp_grid_surface_gates(
    ratios: np.ndarray,
    qubits: np.ndarray,
    settings: Any,
) -> np.ndarray:
    """Evaluate the p=target contour extracted from an external GP grid."""
    if settings.gp_grid_surface_path is None:
        raise ValueError("No gp_grid_surface_path configured.")
    surface = load_gp_grid_surface(str(Path(settings.gp_grid_surface_path)))
    return interp_grid2d_nan(
        surface["qubits"],
        surface["ratios"],
        surface["contour"],
        np.asarray(qubits, dtype=float),
        np.asarray(ratios, dtype=float),
    )
