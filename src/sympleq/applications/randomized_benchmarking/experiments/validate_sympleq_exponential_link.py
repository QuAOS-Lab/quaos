"""
Validate whether local SympleQ RB data follow the assumed exponential link.

For fixed two-qubit ratio ``r`` and register size ``Q``, the cost-aware surface
model assumes

    p(n,r,Q) = V (1 - B(Q)) 2^{-n D(r,Q)} + B(Q),    B(Q)=2^{-Q}.

With V=1, the renormalized survival

    F(n) = (p(n)-B)/(1-B)

should satisfy

    log2 F(n) = -n D(r,Q).

This script measures a fixed grid of depths around the analytic Lindblad
crossing, fits ``log2 F_hat`` versus total gates, and reports whether the
empirical decay rate agrees with the analytic rate.  It is intentionally
independent of the adaptive acquisition in ``cost_aware_surface_design.py``.
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

_SRC_ROOT = Path(__file__).resolve().parents[4]
sys.path = [path for path in sys.path if path != str(_SRC_ROOT)]
sys.path.insert(0, str(_SRC_ROOT))

import matplotlib

SHOW = os.environ.get("SYMPLEQ_LINK_SHOW", "0") != "0"
if not SHOW:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng

from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementRequest,
)
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBData
from sympleq.applications.randomized_benchmarking.experiments.common import (
    spend_request_batch,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_design import (
    BASE_1Q_PAULI_ERROR,
    BASE_2Q_PAULI_ERROR,
    CostAwareSurfaceSettings,
    make_config_q,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    REFERENCE_OFFSET,
    REFERENCE_SLOPE,
)
from sympleq.core.noise.noise_model import GenericNoise


def parse_float_tuple(env_name: str, default: tuple[float, ...]) -> tuple[float, ...]:
    value = os.environ.get(env_name)
    if not value:
        return default
    return tuple(float(part.strip()) for part in value.split(",") if part.strip())


def parse_int_tuple(env_name: str, default: tuple[int, ...]) -> tuple[int, ...]:
    value = os.environ.get(env_name)
    if not value:
        return default
    return tuple(int(part.strip()) for part in value.split(",") if part.strip())


def parse_bool(env_name: str, default: bool) -> bool:
    value = os.environ.get(env_name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off"}


Q_VALUES = parse_int_tuple("SYMPLEQ_LINK_Q_VALUES", (5, 10, 20, 30, 40, 50))
RATIOS = parse_float_tuple("SYMPLEQ_LINK_RATIOS", (0.05, 0.2, 0.5, 0.8, 0.9))
DEPTH_FACTORS = parse_float_tuple(
    "SYMPLEQ_LINK_DEPTH_FACTORS",
    (0.3, 0.5, 0.8, 1.0, 1.3, 1.8),
)
SHOTS_PER_POINT = int(os.environ.get("SYMPLEQ_LINK_SHOTS", "500"))
RNG_SEED = int(os.environ.get("SYMPLEQ_LINK_SEED", "20260629"))
USE_SCRAMBLER = parse_bool("SYMPLEQ_LINK_USE_SCRAMBLER", True)
ONE_Q_NOISE_SCALE = float(os.environ.get("SYMPLEQ_LINK_1Q_SCALE", "1.0"))
TWO_Q_NOISE_SCALE = float(os.environ.get("SYMPLEQ_LINK_2Q_SCALE", "1.0"))
OUT_DIR = (
    Path(__file__).resolve().parent
    / "figs"
    / "sympleq_exponential_link_validation"
)
ROWS_CSV = OUT_DIR / "measurements.csv"
FITS_CSV = OUT_DIR / "fits.csv"
RESULTS_JSON = OUT_DIR / "results.json"


def analytic_rate(ratio: float) -> float:
    """Analytic ``D`` in ``2^{-nD}`` units."""
    return (
        REFERENCE_OFFSET * ONE_Q_NOISE_SCALE
        + REFERENCE_SLOPE * TWO_Q_NOISE_SCALE * ratio
    ) / np.log(2.0)


def analytic_gate_count(ratio: float) -> float:
    return 1.0 / analytic_rate(ratio)


def build_backend(rng) -> SympleqBackend:
    noise_model = GenericNoise.from_paulis(
        [BASE_1Q_PAULI_ERROR * ONE_Q_NOISE_SCALE] * 3,
        rng,
    )
    two_qubit_noise_model = GenericNoise.from_paulis(
        [BASE_2Q_PAULI_ERROR * TWO_Q_NOISE_SCALE] * 3,
        rng,
    )
    return SympleqBackend(
        noise_model=noise_model,
        two_qubit_noise_model=two_qubit_noise_model,
    )


def unique_even_depths(settings: CostAwareSurfaceSettings, ratio: float) -> list[int]:
    n_lo, n_hi = settings.n_gates_bounds
    n_star = analytic_gate_count(ratio)
    depths = []
    for factor in DEPTH_FACTORS:
        depth = int(np.clip(2 * round((factor * n_star) / 2), n_lo, n_hi))
        depths.append(depth)
    return sorted(set(depths))


def config_key(config) -> tuple[int, int, int, int]:
    return (
        int(config.n_qubits),
        int(config.n_gates),
        int(config.n_1qb_gates),
        int(config.n_2qb_gates),
    )


def measure_grid() -> tuple[RMBData, list[dict], dict[tuple[int, int, int, int], float]]:
    rng = default_rng(RNG_SEED)
    backend = build_backend(rng)
    settings = CostAwareSurfaceSettings(
        q_values=Q_VALUES,
        n_qubits=int(round(np.median(Q_VALUES))),
        rng_seed=RNG_SEED,
        use_scrambler=USE_SCRAMBLER,
        initial_one_q_pauli_error=BASE_1Q_PAULI_ERROR * ONE_Q_NOISE_SCALE,
        initial_two_q_pauli_error=BASE_2Q_PAULI_ERROR * TWO_Q_NOISE_SCALE,
    )
    requests = []
    request_rows = []
    target_ratio_by_config = {}
    for q in Q_VALUES:
        for ratio in RATIOS:
            for n_gates in unique_even_depths(settings, ratio):
                config = make_config_q(settings, n_gates, ratio, q)
                requests.append(MeasurementRequest(config, SHOTS_PER_POINT))
                target_ratio_by_config[config_key(config)] = float(ratio)
                request_rows.append({
                    "q": int(q),
                    "target_ratio": float(ratio),
                    "realized_ratio": float(config.ratio_2_qb_gates),
                    "n_gates": int(config.n_gates),
                    "n_1qb_gates": int(config.n_1qb_gates),
                    "n_2qb_gates": int(config.n_2qb_gates),
                    "shots_requested": SHOTS_PER_POINT,
                })

    print(
        f"Measuring {len(requests)} configs, "
        f"{sum(req.shots for req in requests)} total local SympleQ shots"
    )
    data: RMBData = {}
    spend_request_batch(backend, rng, data, requests, seed=RNG_SEED)
    return data, request_rows, target_ratio_by_config


def measurement_rows(
    data: RMBData,
    target_ratio_by_config: dict[tuple[int, int, int, int], float],
) -> list[dict]:
    rows = []
    for config, estimator in data.items():
        realized_ratio = float(config.ratio_2_qb_gates)
        target_ratio = target_ratio_by_config.get(config_key(config), realized_ratio)
        counts = estimator.counts()
        success = int(counts.get(True, 0))
        failure = int(counts.get(False, 0))
        shots = success + failure
        p_hat = success / shots if shots else float("nan")
        b = 2.0 ** (-config.n_qubits)
        f_renorm = (p_hat - b) / (1.0 - b)
        f_clipped = float(np.clip(f_renorm, 1e-9, 1.0))
        rows.append({
            "q": int(config.n_qubits),
            "ratio": float(target_ratio),
            "target_ratio": float(target_ratio),
            "realized_ratio": realized_ratio,
            "n_gates": int(config.n_gates),
            "success": success,
            "failure": failure,
            "shots": shots,
            "p_hat": float(p_hat),
            "asymptote": float(b),
            "renormalized_fidelity": float(f_renorm),
            "log2_renormalized_fidelity": float(np.log2(f_clipped)),
        })
    return sorted(rows, key=lambda row: (row["q"], row["ratio"], row["n_gates"]))


def fit_decay(rows: list[dict]) -> list[dict]:
    fits = []
    grouped: dict[tuple[int, float], list[dict]] = {}
    for row in rows:
        grouped.setdefault((row["q"], row["ratio"]), []).append(row)

    for (q, ratio), group in sorted(grouped.items()):
        xs = np.asarray([row["n_gates"] for row in group], dtype=float)
        ys = np.asarray([row["log2_renormalized_fidelity"] for row in group], dtype=float)
        valid = np.isfinite(xs) & np.isfinite(ys)
        if np.count_nonzero(valid) < 2:
            fits.append(unavailable_fit_row(q, ratio, len(group)))
            continue

        x = xs[valid]
        y = ys[valid]
        slope, intercept = np.polyfit(x, y, deg=1)
        y_fit = slope * x + intercept
        ss_res = float(np.sum((y - y_fit) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
        d_emp = float(-slope)
        d_ref = float(analytic_rate(ratio))
        fits.append({
            "q": int(q),
            "ratio": float(ratio),
            "target_ratio": float(ratio),
            "n_points": int(np.count_nonzero(valid)),
            "realized_ratio_min": float(np.min([row["realized_ratio"] for row in group])),
            "realized_ratio_max": float(np.max([row["realized_ratio"] for row in group])),
            "intercept": float(intercept),
            "slope": float(slope),
            "D_empirical": d_emp,
            "D_analytic": d_ref,
            "D_ratio_empirical_over_analytic": (
                float(d_emp / d_ref) if d_ref > 0.0 else float("nan")
            ),
            "r_squared": float(r_squared),
            "rmse_log2": float(np.sqrt(np.mean((y - y_fit) ** 2))),
        })
    return fits


def unavailable_fit_row(q: int, ratio: float, n_points: int) -> dict:
    return {
        "q": int(q),
        "ratio": float(ratio),
        "target_ratio": float(ratio),
        "n_points": int(n_points),
        "realized_ratio_min": float("nan"),
        "realized_ratio_max": float("nan"),
        "intercept": float("nan"),
        "slope": float("nan"),
        "D_empirical": float("nan"),
        "D_analytic": float(analytic_rate(ratio)),
        "D_ratio_empirical_over_analytic": float("nan"),
        "r_squared": float("nan"),
        "rmse_log2": float("nan"),
    }


def fit_matrix(fits: list[dict], key: str) -> np.ndarray:
    matrix = np.full((len(Q_VALUES), len(RATIOS)), np.nan, dtype=float)
    q_index = {q: i for i, q in enumerate(Q_VALUES)}
    ratio_index = {float(r): i for i, r in enumerate(RATIOS)}
    for row in fits:
        i = q_index.get(int(row["q"]))
        j = ratio_index.get(float(row["ratio"]))
        if i is not None and j is not None:
            matrix[i, j] = float(row[key])
    return matrix


def plot_heatmap(matrix: np.ndarray, title: str, colorbar_label: str, path: Path, *,
                 vmin=None, vmax=None, cmap="viridis") -> None:
    fig, ax = plt.subplots(figsize=(1.35 * len(RATIOS) + 3.2, 0.55 * len(Q_VALUES) + 3.0))
    image = ax.imshow(matrix, origin="lower", aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(np.arange(len(RATIOS)), [f"{r:.2g}" for r in RATIOS])
    ax.set_yticks(np.arange(len(Q_VALUES)), [str(q) for q in Q_VALUES])
    ax.set_xlabel("two-qubit ratio r")
    ax.set_ylabel("n_qubits Q")
    ax.set_title(title)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(image, ax=ax, label=colorbar_label)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    if not SHOW:
        plt.close(fig)


def plot_decay_lines(rows: list[dict], fits: list[dict], path: Path) -> None:
    n_rows = len(Q_VALUES)
    n_cols = len(RATIOS)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(3.4 * n_cols, 2.45 * n_rows),
        squeeze=False,
        sharex=False,
        sharey=False,
        constrained_layout=True,
    )
    by_fit = {(int(row["q"]), float(row["ratio"])): row for row in fits}
    by_data: dict[tuple[int, float], list[dict]] = {}
    for row in rows:
        by_data.setdefault((int(row["q"]), float(row["ratio"])), []).append(row)

    for row_index, q in enumerate(Q_VALUES):
        for col_index, ratio in enumerate(RATIOS):
            ax = axes[row_index, col_index]
            group = by_data.get((q, float(ratio)), [])
            fit = by_fit.get((q, float(ratio)))
            if not group or fit is None:
                ax.set_axis_off()
                continue
            x = np.asarray([row["n_gates"] for row in group], dtype=float)
            y = np.asarray([row["log2_renormalized_fidelity"] for row in group], dtype=float)
            ax.scatter(x, y, s=24, color="tab:blue", label="SympleQ")
            if np.isfinite(fit["slope"]):
                xx = np.linspace(float(np.min(x)), float(np.max(x)), 100)
                yy = fit["slope"] * xx + fit["intercept"]
                analytic = -analytic_rate(ratio) * xx
                ax.plot(xx, yy, color="black", linewidth=1.4, label="fit")
                ax.plot(xx, analytic, color="crimson", linestyle="--", linewidth=1.2, label="analytic")
            ax.axhline(-1.0, color="0.75", linewidth=0.8)
            ax.set_title(
                f"Q={q}, r={ratio:.2g}\n"
                f"D/Dref={fit['D_ratio_empirical_over_analytic']:.2f}, "
                f"R2={fit['r_squared']:.2f}",
                fontsize=9,
            )
            if row_index == n_rows - 1:
                ax.set_xlabel("total gates n")
            if col_index == 0:
                ax.set_ylabel(r"$\log_2 F_{\rm renorm}$")
            if row_index == 0 and col_index == 0:
                ax.legend(fontsize=8)
    fig.savefig(path, dpi=180)
    if not SHOW:
        plt.close(fig)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data, request_rows, target_ratio_by_config = measure_grid()
    rows = measurement_rows(data, target_ratio_by_config)
    fits = fit_decay(rows)
    write_csv(ROWS_CSV, rows)
    write_csv(FITS_CSV, fits)

    d_ratio = fit_matrix(fits, "D_ratio_empirical_over_analytic")
    intercept = fit_matrix(fits, "intercept")
    r_squared = fit_matrix(fits, "r_squared")

    plot_heatmap(
        d_ratio,
        "Empirical / analytic decay rate",
        "D_emp / D_analytic",
        OUT_DIR / "decay_rate_ratio_heatmap.png",
        vmin=0.75,
        vmax=1.25,
        cmap="coolwarm",
    )
    plot_heatmap(
        intercept,
        "Fitted intercept in log2-renormalized survival",
        "intercept",
        OUT_DIR / "intercept_heatmap.png",
        cmap="coolwarm",
    )
    plot_heatmap(
        r_squared,
        "Linearity of log2-renormalized survival",
        "R^2",
        OUT_DIR / "r_squared_heatmap.png",
        vmin=0.0,
        vmax=1.0,
        cmap="viridis",
    )
    plot_decay_lines(rows, fits, OUT_DIR / "decay_line_fits.png")

    finite_ratio = d_ratio[np.isfinite(d_ratio)]
    summary = {
        "q_values": list(Q_VALUES),
        "ratios": list(RATIOS),
        "depth_factors": list(DEPTH_FACTORS),
        "shots_per_point": SHOTS_PER_POINT,
        "rng_seed": RNG_SEED,
        "use_scrambler": USE_SCRAMBLER,
        "one_q_noise_scale": ONE_Q_NOISE_SCALE,
        "two_q_noise_scale": TWO_Q_NOISE_SCALE,
        "n_configs": len(rows),
        "n_requested_configs": len(request_rows),
        "mean_D_ratio": float(np.mean(finite_ratio)) if len(finite_ratio) else None,
        "min_D_ratio": float(np.min(finite_ratio)) if len(finite_ratio) else None,
        "max_D_ratio": float(np.max(finite_ratio)) if len(finite_ratio) else None,
        "outputs": {
            "measurements_csv": str(ROWS_CSV),
            "fits_csv": str(FITS_CSV),
            "decay_rate_ratio_heatmap": str(OUT_DIR / "decay_rate_ratio_heatmap.png"),
            "intercept_heatmap": str(OUT_DIR / "intercept_heatmap.png"),
            "r_squared_heatmap": str(OUT_DIR / "r_squared_heatmap.png"),
            "decay_line_fits": str(OUT_DIR / "decay_line_fits.png"),
        },
    }
    RESULTS_JSON.write_text(json.dumps({"summary": summary, "fits": fits}, indent=2), encoding="utf-8")

    print("\nSummary")
    print(f"  configs measured: {len(rows)}")
    print(f"  shots per point: {SHOTS_PER_POINT}")
    print(f"  use scrambler: {USE_SCRAMBLER}")
    if len(finite_ratio):
        print(f"  mean D_emp/D_analytic: {np.mean(finite_ratio):.3f}")
        print(f"  range D_emp/D_analytic: {np.min(finite_ratio):.3f} to {np.max(finite_ratio):.3f}")
    print(f"  wrote outputs to {OUT_DIR}")
    if SHOW:
        plt.show()


if __name__ == "__main__":
    main()
