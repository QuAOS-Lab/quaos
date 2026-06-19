"""Clean score benchmark for ``rick_modified_crossing.py``.

Edit only the handles near the top for ordinary runs:

* ``N_NOISE_VALUES`` and ``NOISE_SEED`` choose the noise multipliers.
* ``REALISATION_SEEDS`` chooses the Rick/AEPsych random seeds per noise value.

Each Rick run writes its own RMB JSON, GP-grid NPZ, level-set plot, and scores.
This benchmark only collects those saved scores into one CSV/JSON summary.
"""
from __future__ import annotations

import csv
import json
import sys
from datetime import datetime
from pathlib import Path
from statistics import mean

from numpy.random import Generator as RNGGenerator, default_rng

_SRC_ROOT = Path(__file__).resolve().parents[4]
sys.path = [path for path in sys.path if path != str(_SRC_ROOT)]
sys.path.insert(0, str(_SRC_ROOT))

from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.experiments.common import CrossingSettings
from sympleq.applications.randomized_benchmarking.experiments.monotone_tracing import (
    BASE_1Q_PAULI_ERROR,
    BASE_2Q_PAULI_ERROR,
)
from sympleq.applications.randomized_benchmarking.experiments.rick_modified_crossing import (
    RickFantasyGPUSettings,
    run_with_budget,
)
from sympleq.core.noise.noise_model import GenericNoise


# ---------------------------------------------------------------------------
# BENCHMARK HANDLES
# ---------------------------------------------------------------------------

N_NOISE_VALUES = 10
NOISE_SEED = 20260618
NOISE_SPREAD = 0.20

# These are the realisation ids inside each noise value. The actual Rick seed is
# ``10_000 * noise_index + realisation_seed``.
REALISATION_SEEDS = list(range(10))

HQC_BUDGET = 200.0
MAX_COST_PER_RUN = 15.0

FIG_DIR = (
    Path(__file__).resolve().parent
    / "figs"
    / "rick_modified_scores"
    / datetime.now().strftime("%Y%m%d_%H%M%S")
)
RESULTS_JSON = FIG_DIR / "results.json"
RUNS_CSV = FIG_DIR / "runs.csv"


def noise_values() -> list[tuple[float, float]]:
    """Deterministic one-qubit and two-qubit noise multipliers."""
    rng = default_rng(NOISE_SEED)
    return [
        tuple(float(x) for x in rng.uniform(
            1.0 - NOISE_SPREAD,
            1.0 + NOISE_SPREAD,
            size=2,
        ))
        for _ in range(N_NOISE_VALUES)
    ]


def scaled_noise_backend_factory(
    one_q_noise_scale: float,
    two_q_noise_scale: float,
):
    """Build a SympleQ backend with scaled Lindblad noise."""

    def factory(settings: CrossingSettings, rng: RNGGenerator) -> RMBBackend:
        noise_model = GenericNoise.from_paulis(
            [BASE_1Q_PAULI_ERROR * one_q_noise_scale] * 3,
            rng,
        )
        two_qubit_noise_model = GenericNoise.from_paulis(
            [BASE_2Q_PAULI_ERROR * two_q_noise_scale] * 3,
            rng,
        )
        return SympleqBackend(
            noise_model=noise_model,
            two_qubit_noise_model=two_qubit_noise_model,
        )

    return factory


def make_settings(
    *,
    seed: int,
    save_path: Path,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> RickFantasyGPUSettings:
    """Settings for one Rick modified run."""
    settings = RickFantasyGPUSettings(
        rng_seed=seed,
        save_path=save_path,
        plot=True,
        verbose_fantasies=False,
        print_diagnostics=False,
        backend_factory=scaled_noise_backend_factory(
            one_q_noise_scale,
            two_q_noise_scale,
        ),
        hqc_budget=HQC_BUDGET,
        max_cost_per_run=MAX_COST_PER_RUN,
        n_gates_bounds=(10, 5000),
        ratio_bounds=(0.1, 0.9),
        use_fake_corners=True,
        easy_corner_outcome=1,
        hard_corner_outcome=0,
        extra_fake_anchors=[],
        initial_sobol_samples=10,
        initial_sobol_max_cost_per_run=30.0,
        sobol_scramble=True,
        validate_gp_contour=False,
        reserve_validation_budget=False,
        validation_hqc_budget=0.0,
        plot_gp_level_set_result=True,
        save_gp_prediction_grid=True,
        plots_module="plots_1",
        bayesian_estimation_module="sympleq.core.bayesian_estimation",
        estimator_threshold=0.0,
        estimator_min_runs=1,
        optimization_steps=500,
        inducing_size=150,
        acquisition_function="GlobalSUR",
        acquisition_restarts=2,
        acquisition_samples=300,
        batching=True,
        max_batch_size=80,
        fantasy_batching=True,
        use_gpu=False,
        force_default_device_during_aepsych=True,
    )
    object.__setattr__(settings, "one_q_noise_scale", one_q_noise_scale)
    object.__setattr__(settings, "two_q_noise_scale", two_q_noise_scale)
    return settings


def read_scores(path: Path) -> dict[str, float]:
    """Read the scores written by ``rick_modified_crossing.py``."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    return dict(payload.get("scores", {}))


def run_one(
    *,
    noise_index: int,
    realisation_seed: int,
    one_q_noise_scale: float,
    two_q_noise_scale: float,
) -> dict:
    seed = 10_000 * noise_index + int(realisation_seed)
    run_dir = FIG_DIR / f"noise_{noise_index:02d}" / f"seed_{seed:06d}"
    run_dir.mkdir(parents=True, exist_ok=True)
    save_path = run_dir / "rick_modified_crossing.json"

    settings = make_settings(
        seed=seed,
        save_path=save_path,
        one_q_noise_scale=one_q_noise_scale,
        two_q_noise_scale=two_q_noise_scale,
    )
    rmb, crossings, budget, validation_rmb = run_with_budget(settings)
    scores = read_scores(save_path)

    return {
        "noise_index": noise_index,
        "realisation_seed": int(realisation_seed),
        "seed": seed,
        "one_q_noise_scale": one_q_noise_scale,
        "two_q_noise_scale": two_q_noise_scale,
        "S1": scores.get("S1"),
        "S2": scores.get("S2"),
        "A_gp": scores.get("A_gp"),
        "mean_delta_log_gates": scores.get("mean_delta_log_gates"),
        "mean_sigma_contour": scores.get("mean_sigma_contour"),
        "n_contour_points": scores.get("n_contour_points"),
        "spent_hqc": budget.spent_hqc,
        "n_crossings": len(crossings),
        "n_training_configs": len(rmb._data),
        "n_validation_configs": 0 if validation_rmb is None else len(validation_rmb._data),
        "json_path": str(save_path),
        "grid_path": str(save_path.with_name("rick_modified_crossing_gp_grid.npz")),
        "level_set_figure": str(save_path.with_name("rick_modified_crossing_level_set.png")),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []

    for noise_index, (one_q, two_q) in enumerate(noise_values()):
        for realisation_seed in REALISATION_SEEDS:
            row = run_one(
                noise_index=noise_index,
                realisation_seed=realisation_seed,
                one_q_noise_scale=one_q,
                two_q_noise_scale=two_q,
            )
            rows.append(row)
            print(
                f"noise={noise_index:02d} seed={row['seed']:06d} "
                f"S1={row['S1']:.6g} S2={row['S2']:.6g} "
                f"spent={row['spent_hqc']:.1f} HQC"
            )

    summary = {
        "n_noise_values": N_NOISE_VALUES,
        "realisation_seeds": REALISATION_SEEDS,
        "n_runs": len(rows),
        "noise_seed": NOISE_SEED,
        "noise_spread": NOISE_SPREAD,
        "mean_S1": mean(row["S1"] for row in rows),
        "mean_S2": mean(row["S2"] for row in rows),
        "mean_spent_hqc": mean(row["spent_hqc"] for row in rows),
    }
    RESULTS_JSON.write_text(
        json.dumps({"summary": summary, "runs": rows}, indent=2),
        encoding="utf-8",
    )
    write_csv(RUNS_CSV, rows)

    print("\nSummary")
    print(f"  runs: {summary['n_runs']}")
    print(f"  mean S1: {summary['mean_S1']:.6g}")
    print(f"  mean S2: {summary['mean_S2']:.6g}")
    print(f"  mean spent: {summary['mean_spent_hqc']:.1f} HQC")
    print(f"  output: {FIG_DIR}")


if __name__ == "__main__":
    main()
