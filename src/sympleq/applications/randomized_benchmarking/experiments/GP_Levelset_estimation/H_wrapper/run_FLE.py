"""
Algorithm
---------
1. Add fake anchor points:
       easy corner -> assumed success
       hard corner -> assumed failure

2. Generate explicit Sobol warm-up points.

3. Measure Sobol batch using the RMB backend (SympleQ/Quantinuum).

4. Fit/update AEPsych GP classifier.

5. Build GlobalSUR batches using fantasy outcomes:
       copy real strategy
       ask copied strategy for GlobalSUR point
       predict P(success)
       sample virtual success/failure
       add virtual outcome to copied strategy only
       repeat until batch/budget is full

6. Measure the selected batch for real on CPU RMB backend.

7. Add only real outcomes to the real strategy.

8. Repeat until HQC budget is exhausted.

9. Extract the GP-predicted P(success)=target_threshold contour.

"""

from __future__ import annotations


import logging
import json
import re
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import torch

from aepsych.strategy import SequentialStrategy

from sympleq.applications.randomized_benchmarking.RMB import RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.experiments.common import (
    default_backend_factory,
    measurement_rng,
    print_experiment_summary,
    print_progress,
    start_run,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.fantasy_levelset_estimation import (
    Observation,
    add_observation_to_strategy,
    build_strategy,
    choose_gp_device,
    contour_target,
    level_set_configs,
    measure_batch_and_update_real_strategy,
    move_strategy_models_to_device,
    point_from_config,
    repair_stats,
    refreshed_strategy_for_prediction,
    reset_repair_stats,
    seed_fake_corners,
    selected_config_metadata,
    select_affordable_prefix,
    select_fantasy_globalsur_batch,
    select_plain_aepsych_batch,
    sobol_initial_candidates,
    temporary_torch_default_device,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.fantasy_levelset_settings import (
    FantasySettings,
    RNG_SEEDS,
    control_panel_settings_kwargs,
)


# Storing

def backend_folder_from_kwargs(seed_kwargs: dict) -> str:
    """Return the run-folder backend name from the backend's device_name."""

    settings = FantasySettings(**seed_kwargs)
    backend_factory = seed_kwargs.get("backend_factory", default_backend_factory)
    backend = backend_factory(settings, np.random.default_rng(settings.rng_seed))
    device_name = getattr(backend, "device_name", type(backend).__name__)
    return str(device_name).replace("-", "_")


def timestamped_personal_save_path(
    seed: int | None = None,
    *,
    timestamp: str | None = None,
    n_qubits: int | None = None,
    multi_qubit: bool = False,
    backend_folder: str = "backend_unknown",
) -> Path:
    """Timestamped 2D FLE JSON path under ``Personal/FLE``."""

    timestamp = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    seed_folder = "seed_unseeded" if seed is None else f"seed_{seed}"
    qubit_folder = "q_unknown" if n_qubits is None else f"q{int(n_qubits)}"
    suffix = f"_q{int(n_qubits)}" if multi_qubit and n_qubits is not None else ""
    run_folder = (
        Path("Personal")
        / "FLE"
        / "H_wrapper"
        / backend_folder
        / qubit_folder
        / seed_folder
        / f"FLE_H_wrapper_{timestamp}"
    )
    return run_folder / f"FLE_{timestamp}{suffix}.json"


def timestamped_recovery_save_path(recovery_folder: str | Path) -> Path:
    """Timestamped final JSON path inside an existing recovery folder."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path(recovery_folder) / f"FLE_recovery_{timestamp}.json"


def qubit_slices(value) -> list[int]:
    if isinstance(value, (list, tuple)):
        return [int(q) for q in value]
    return [int(value)]


# Prediction
def predict_native_level_set(
    strategy: SequentialStrategy,
    settings: FantasySettings,
    *,
    device: torch.device,
    n_grid: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate fitted GP on its native total-gates/ratio mesh."""
    gates_axis = np.geomspace(
        max(1.0, float(settings.n_gates_bounds[0])),
        float(settings.n_gates_bounds[1]),
        n_grid,
    )
    ratio_axis = np.linspace(
        float(settings.ratio_bounds[0]),
        float(settings.ratio_bounds[1]),
        n_grid,
    )
    gates_grid, ratio_grid = np.meshgrid(gates_axis, ratio_axis)
    grid = torch.tensor(
        np.column_stack([gates_grid.ravel(), ratio_grid.ravel()]),
        dtype=torch.double,
        device=device,
    )

    move_strategy_models_to_device(strategy, device)

    with temporary_torch_default_device(
        device,
        enabled=settings.force_default_device_during_aepsych,
    ):
        with torch.no_grad():
            probabilities, _ = strategy.model.predict(
                grid,
                probability_space=True,
            )
            latent_mean, latent_variance = strategy.model.predict(grid)

    return (
        probabilities.detach().cpu().numpy().reshape(gates_grid.shape),
        latent_mean.detach().cpu().numpy().reshape(gates_grid.shape),
        latent_variance.detach().cpu().numpy().reshape(gates_grid.shape),
        gates_grid,
        ratio_grid,
    )


def save_gp_prediction_grid(
    strategy: SequentialStrategy,
    settings: FantasySettings,
    *,
    device: torch.device,
    json_path: Path,
) -> Path:
    (
        probabilities,
        latent_mean,
        latent_variance,
        gates_grid,
        ratio_grid,
    ) = predict_native_level_set(
        strategy,
        settings,
        device=device,
    )

    grid_path = json_path.parent / f"{json_path.stem}_gp_grid.npz"
    np.savez_compressed(
        grid_path,
        gates_grid=gates_grid,
        ratio_grid=ratio_grid,
        x_grid=gates_grid,
        y_grid=ratio_grid,
        probabilities=probabilities,
        latent_mean=latent_mean,
        latent_variance=latent_variance,
        target=np.asarray(contour_target(settings), dtype=float),
        coordinate_system=np.asarray("total_ratio"),
        rng_seed=np.asarray(
            -1 if settings.rng_seed is None else settings.rng_seed,
            dtype=int,
        ),
    )
    print(f"[saved] GP grid: {grid_path}")
    return grid_path


def sent_config_payload(config) -> dict:
    payload = {
        "n_1qb_gates": int(config.n_1qb_gates),
        "n_2qb_gates": int(config.n_2qb_gates),
        "n_qubits": int(config.n_qubits),
        "n_gates": int(config.n_gates),
        "ratio_2_qb_gates": float(config.ratio_2_qb_gates),
        "random_elimination": float(config.random_elimination),
        "use_scrambler": bool(config.use_scrambler),
    }
    if hasattr(config, "circuit_metadata"):
        payload.update(config.circuit_metadata())
    return payload


def actual_after_elimination_payload(
    config,
    *,
    rng_seed: int | None,
    first_shot_index: int | None,
) -> dict:
    if rng_seed is None or first_shot_index is None:
        return {
            "n_1qb_gates": None,
            "n_2qb_gates": None,
            "n_qubits": int(config.n_qubits),
            "n_gates": None,
            "ratio_2_qb_gates": None,
            "first_shot_index": first_shot_index,
            "available": False,
        }

    circuit = config.random_circuit(
        rng=measurement_rng(int(rng_seed), config, int(first_shot_index))
    )
    n_1q = sum(
        1
        for gate in circuit.gates
        if gate.n_qudits == 1 and gate.name != "Id"
    )
    n_2q = sum(
        1
        for gate in circuit.gates
        if gate.n_qudits == 2 and gate.name != "Id"
    )
    n_gates = n_1q + n_2q
    payload = {
        "n_1qb_gates": int(n_1q),
        "n_2qb_gates": int(n_2q),
        "n_qubits": int(config.n_qubits),
        "n_gates": int(n_gates),
        "ratio_2_qb_gates": float(n_2q / n_gates) if n_gates else 0.0,
        "first_shot_index": int(first_shot_index),
        "available": True,
    }
    if hasattr(config, "circuit_metadata"):
        payload.update(config.circuit_metadata())
    return payload


def sent_batch_metadata(
    settings: FantasySettings,
    sent_configs,
    data=None,
) -> list[dict]:
    metadata = []
    batch_offsets = {}
    for config in sent_configs:
        first_shot_index = None
        if data is not None:
            estimator = data.get(config) if hasattr(data, "get") else None
            previous_runs = (
                int(estimator.num_runs())
                if estimator is not None and hasattr(estimator, "num_runs")
                else 0
            )
            offset = batch_offsets.get(config, 0)
            batch_offsets[config] = offset + 1
            first_shot_index = previous_runs + offset

        metadata.append(
            {
                "sent_config": sent_config_payload(config),
                "proposal": selected_config_metadata(config),
                "requested_shots": 1,
                "first_shot_index": first_shot_index,
                "actual_after_elimination": actual_after_elimination_payload(
                    config,
                    rng_seed=settings.rng_seed,
                    first_shot_index=first_shot_index,
                ),
            }
        )
    return metadata


def write_hqc_metadata(
    json_path: Path,
    settings: FantasySettings,
    budget,
    phase: str,
    step: int | None,
    sent_configs=None,
    batch_metadata: list[dict] | None = None,
) -> None:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    payload["experiment"] = {
        "phase": phase,
        "step": step,
        "spent_hqc": float(settings.hqc_budget - budget.remaining_hqc),
        "remaining_hqc": float(budget.remaining_hqc),
        "hqc_budget": float(settings.hqc_budget),
        "repair_stats": repair_stats(),
        "circuit_variant": getattr(settings, "circuit_variant", "unknown"),
        "protected_h_wrapper": bool(getattr(settings, "protected_h_wrapper", False)),
        "h_wrapper_1q_gates_per_qubit": int(
            getattr(settings, "h_wrapper_1q_gates_per_qubit", 0)
        ),
    }
    if sent_configs is not None:
        if batch_metadata is None:
            batch_metadata = sent_batch_metadata(settings, sent_configs)
        payload["experiment"]["sent_configs"] = [
            item["sent_config"]
            for item in batch_metadata
        ]
        payload["experiment"]["proposal_configs"] = [
            item["proposal"]
            for item in batch_metadata
        ]
        payload["experiment"]["actual_after_elimination"] = [
            item["actual_after_elimination"]
            for item in batch_metadata
        ]
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _settings_run_folder(settings: FantasySettings) -> Path | None:
    if settings.save_path is not None:
        return Path(settings.save_path).parent
    if settings.recovery_folder is not None:
        return Path(settings.recovery_folder)
    return None


def pending_backend_batch_path(settings: FantasySettings) -> Path | None:
    folder = _settings_run_folder(settings)
    if folder is None:
        return None
    return folder / "pending_backend_batch.json"


def pending_settings_snapshot(settings: FantasySettings) -> dict:
    backend_factory_name = getattr(
        settings.backend_factory,
        "__name__",
        type(settings.backend_factory).__name__,
    )
    return {
        "rng_seed": settings.rng_seed,
        "n_qubits": int(settings.n_qubits),
        "n_gates_bounds": [float(value) for value in settings.n_gates_bounds],
        "ratio_bounds": [float(value) for value in settings.ratio_bounds],
        "hqc_budget": float(settings.hqc_budget),
        "max_cost_per_run": float(settings.max_cost_per_run),
        "backend_factory": backend_factory_name,
        "save_path": None if settings.save_path is None else str(settings.save_path),
        "circuit_variant": getattr(settings, "circuit_variant", "unknown"),
        "protected_h_wrapper": bool(getattr(settings, "protected_h_wrapper", False)),
        "h_wrapper_1q_gates_per_qubit": int(
            getattr(settings, "h_wrapper_1q_gates_per_qubit", 0)
        ),
    }


def write_pending_backend_batch(
    *,
    settings: FantasySettings,
    phase: str,
    step: int,
    sent_configs,
    data=None,
    batch_metadata: list[dict] | None = None,
) -> None:
    """Persist the exact FLE batch before submitting stitched circuits."""

    if not settings.save_real_checkpoints:
        return

    path = pending_backend_batch_path(settings)
    if path is None:
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    requests = []
    if batch_metadata is None:
        batch_metadata = sent_batch_metadata(settings, sent_configs, data=data)
    for item in batch_metadata:
        requests.append(
            {
                "config": item["sent_config"],
                "proposal": item["proposal"],
                "requested_shots": item["requested_shots"],
                "first_shot_index": item["first_shot_index"],
                "actual_after_elimination": item["actual_after_elimination"],
            }
        )

    payload = {
        "reason": "before_backend_submit",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "model": "FLE_H_wrapper",
        "phase": phase,
        "step": int(step),
        "settings": pending_settings_snapshot(settings),
        "pending_batch": {
            "requests": requests,
        },
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[saved] pending backend batch: {path}")


def clear_pending_backend_batch(settings: FantasySettings) -> None:
    path = pending_backend_batch_path(settings)
    if path is not None and path.exists():
        path.unlink()


def latest_recovery_json(settings: FantasySettings) -> Path | None:
    """Return the latest cumulative RMB JSON checkpoint for recovery."""

    if settings.recovery_folder is not None:
        folder = Path(settings.recovery_folder)
    elif settings.save_path is not None:
        folder = Path(settings.save_path).parent
    else:
        return None

    candidates = sorted(folder.glob("measurement_*.json"))
    if not candidates:
        candidates = sorted(folder.glob("FLE_*.json"))

    return candidates[-1] if candidates else None


def measurement_checkpoint_step(json_path: Path) -> int:
    """Return the numeric step from a measurement checkpoint filename."""

    match = re.match(r"measurement_(\d+)_", json_path.name)
    return int(match.group(1)) if match else 0


def latest_measurement_checkpoint_step(folder: Path) -> int:
    """Return the largest measurement checkpoint step already in a folder."""

    return max(
        (measurement_checkpoint_step(path) for path in folder.glob("measurement_*.json")),
        default=0,
    )


def recovery_hqc_spent(json_path: Path) -> float | None:
    """Read spent HQC metadata from a recovery checkpoint if present."""

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    experiment = payload.get("experiment", {})
    spent_hqc = experiment.get("spent_hqc")
    return None if spent_hqc is None else float(spent_hqc)


def recover_observations_from_json(
    json_path: Path,
    *,
    strategy: SequentialStrategy,
    rmb,
    observations: list[Observation],
    results_for_plot: list[tuple[float, float, int]],
    device: torch.device,
) -> int:
    """Load saved RMB data and replay its real observations into AEPsych."""

    recovered_rmb = RMB.load(json_path)
    rmb._data.update(recovered_rmb._data)

    n_observations = 0
    for config, estimator in recovered_rmb._data.items():
        for outcome, count in estimator.counts().items():
            for _ in range(int(count)):
                obs = Observation(
                    x_cpu=point_from_config(config, device=torch.device("cpu")),
                    y=int(outcome),
                    source="recovery",
                )
                add_observation_to_strategy(strategy, obs, device=device)
                observations.append(obs)
                results_for_plot.append(
                    (
                        float(config.n_gates),
                        float(config.ratio_2_qb_gates),
                        int(outcome),
                    )
                )
                n_observations += 1

    return n_observations


def save_real_checkpoint(
    *,
    rmb,
    strategy: SequentialStrategy,
    settings: FantasySettings,
    observations: list[Observation],
    device: torch.device,
    phase: str,
    step: int,
    budget,
    sent_configs=None,
    batch_metadata: list[dict] | None = None,
) -> None:
    if not settings.save_real_checkpoints or settings.save_path is None:
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    checkpoint_path = (
        Path(settings.save_path).parent
        / f"measurement_{step:03d}_{phase}_{timestamp}.json"
    )
    rmb.save(checkpoint_path)
    checkpoint_path = resolve_data_path(checkpoint_path)
    write_hqc_metadata(
        checkpoint_path,
        settings,
        budget,
        phase,
        step,
        sent_configs,
        batch_metadata=batch_metadata,
    )
    print(f"[saved] real RMB checkpoint: {checkpoint_path}")

    if not settings.save_gp_prediction_grid:
        return

    try:
        checkpoint_strategy = refreshed_strategy_for_prediction(
            strategy,
            settings,
            observations,
            device=device,
        )
        save_gp_prediction_grid(
            checkpoint_strategy,
            settings,
            device=device,
            json_path=checkpoint_path,
        )
    except Exception as exc:
        print(f"[checkpoint] GP grid was not saved: {exc}")


def print_run_handles(
    settings: FantasySettings,
    *,
    gp_device: torch.device,
) -> None:
    """
    Print the tunable handles for reproducibility.
    """

    print("\n================ RUN HANDLES ================")

    print("[target]")
    print(f"  target_threshold                     = {settings.target_threshold}")

    print("[qubits]")
    print(f"  qubits                     = {settings.n_qubits}")

    print("[fake anchors]")
    print(f"  use_fake_corners                     = {settings.use_fake_corners}")
    print(f"  easy_corner_outcome                  = {settings.easy_corner_outcome}")
    print(f"  hard_corner_outcome                  = {settings.hard_corner_outcome}")
    print(f"  extra_fake_anchors                   = {settings.extra_fake_anchors}")

    print("[sobol warm-up]")
    print(f"  initial_sobol_samples                = {settings.initial_sobol_samples}")
    print(f"  initial_sobol_max cost_per_run       = "
          f"{settings.initial_sobol_max_cost_per_run}")
    print(f"  sobol_scramble                       = {settings.sobol_scramble}")

    print("[recovery]")
    print(f"  recovery_mode                        = {settings.recovery_mode}")
    print(f"  recovery_folder                      = {settings.recovery_folder}")

    print("[GP / AEPsych]")
    print(f"  acquisition_function                 = {settings.acquisition_function}")
    print(f"  optimization_steps                   = {settings.optimization_steps}")
    print(f"  inducing_size                        = {settings.inducing_size}")
    print(f"  acquisition_restarts                 = {settings.acquisition_restarts}")
    print(f"  acquisition_samples                  = {settings.acquisition_samples}")

    print("[batching / budget]")
    print(f"  batching                             = {settings.batching}")
    print(f"  max_batch_size                       = {settings.max_batch_size}")
    print(f"  fantasy_batching                     = {settings.fantasy_batching}")
    print(f"  max_cost_per_run                     = {settings.max_cost_per_run}")
    print(f"  hqc_budget                           = {settings.hqc_budget}")

    print("[search box]")
    print(f"  n_gates_bounds                       = {settings.n_gates_bounds}")
    print(f"  ratio_bounds                         = {settings.ratio_bounds}")

    print("[device]")
    print(f"  use_gpu                              = {settings.use_gpu}")
    print(f"  selected GP device                   = {gp_device}")
    print(f"  force_default_device_during_aepsych  = "
          f"{settings.force_default_device_during_aepsych}")
    backend_factory_name = getattr(
        settings.backend_factory,
        "__name__",
        type(settings.backend_factory).__name__,
    )
    print(f"  backend_factory                      = {backend_factory_name}")

    print("[debug]")
    print(f"  rng_seed                             = {settings.rng_seed}")
    print(f"  verbose_fantasies                    = {settings.verbose_fantasies}")

    print("=============================================\n")


# =============================================================================
# Main run logic
# =============================================================================


def run(
    settings: FantasySettings,
    *,
    return_budget: bool = False,
):
    """
    Run the experiment.
    """

    logging.getLogger().setLevel(logging.WARNING)
    warnings.filterwarnings("ignore")

    torch.set_default_dtype(torch.float64)

    gp_device = choose_gp_device(settings)

    if settings.rng_seed is not None:
        torch.manual_seed(settings.rng_seed)
        if gp_device.type == "cuda":
            torch.cuda.manual_seed_all(settings.rng_seed)

    print_run_handles(settings, gp_device=gp_device)
    reset_repair_stats()

    rng, rmb, budget = start_run(settings)
    backend_details = [
        f"{name}={value}"
        for name in ("device_name", "project_name")
        if (value := getattr(rmb.backend, name, None)) is not None
    ]
    backend_text = type(rmb.backend).__name__
    if backend_details:
        backend_text = f"{backend_text} ({', '.join(backend_details)})"
    # print("[RMB backend]")
    print(f"  actual rmb.backend                   = {backend_text}")

    data = rmb._data

    fantasy_seed = None if settings.rng_seed is None else settings.rng_seed + 99173
    fantasy_rng = np.random.default_rng(fantasy_seed)

    strategy = build_strategy(settings)

    observations: list[Observation] = []
    results_for_plot: list[tuple[float, float, int]] = []
    checkpoint_step = 0

    # -------------------------------------------------------------------------
    # 1. Fake anchors
    # -------------------------------------------------------------------------

    seed_fake_corners(
        strategy,
        settings,
        observations,
        results_for_plot,
        device=gp_device,
    )

    recovered_observations = 0
    if settings.recovery_mode:
        recovery_json = latest_recovery_json(settings)
        if recovery_json is None:
            print("[recovery] enabled, but no checkpoint JSON found; running Sobol.")
        else:
            output_folder = (
                Path(settings.save_path).parent
                if settings.save_path is not None
                else recovery_json.parent
            )
            checkpoint_step = latest_measurement_checkpoint_step(output_folder)
            print(f"[recovery] next measurement step starts after {checkpoint_step}")

            recovered_hqc_spent = recovery_hqc_spent(recovery_json)
            if recovered_hqc_spent is None:
                print(f"[recovery] HQC spent unavailable in {recovery_json}")
            else:
                budget.spent_hqc = recovered_hqc_spent
                budget.remaining_hqc = max(
                    0.0,
                    float(settings.hqc_budget) - recovered_hqc_spent,
                )
                print(f"[recovery] recovered HQC spent = {recovered_hqc_spent:.6g}")
                print(f"[recovery] remaining HQC budget = {budget.remaining_hqc:.6g}")

            recovered_observations = recover_observations_from_json(
                recovery_json,
                strategy=strategy,
                rmb=rmb,
                observations=observations,
                results_for_plot=results_for_plot,
                device=gp_device,
            )
            if recovered_observations:
                print(
                    f"[recovery] loaded {recovered_observations} observations "
                    f"from {recovery_json}; skipping Sobol warm-up."
                )
            else:
                print(
                    f"[recovery] found {recovery_json}, but it contained no "
                    "observations; running Sobol."
                )

    # -------------------------------------------------------------------------
    # 2. Sobol warm-up batch
    # -------------------------------------------------------------------------

    exhausted = False

    if recovered_observations:
        print_progress(settings, budget, "sobol: skipped after recovery")
    else:
        sobol_candidates = sobol_initial_candidates(settings)
        sobol_batch, exhausted = select_affordable_prefix(
            sobol_candidates,
            settings,
            budget,
            max_cost_per_run=settings.initial_sobol_max_cost_per_run,
        )

    if not recovered_observations and sobol_batch:
        pending_step = checkpoint_step + 1
        sobol_batch_metadata = sent_batch_metadata(settings, sobol_batch, data=data)
        write_pending_backend_batch(
            settings=settings,
            phase="sobol",
            step=pending_step,
            sent_configs=sobol_batch,
            data=data,
            batch_metadata=sobol_batch_metadata,
        )

        measure_batch_and_update_real_strategy(
            phase="sobol",
            selected=sobol_batch,
            strategy=strategy,
            rmb=rmb,
            rng=rng,
            data=data,
            budget=budget,
            settings=settings,
            observations=observations,
            results_for_plot=results_for_plot,
            device=gp_device,
        )
        checkpoint_step += 1
        save_real_checkpoint(
            rmb=rmb,
            strategy=strategy,
            settings=settings,
            observations=observations,
            device=gp_device,
            phase="sobol",
            step=checkpoint_step,
            budget=budget,
            sent_configs=sobol_batch,
            batch_metadata=sobol_batch_metadata,
        )
        clear_pending_backend_batch(settings)
    elif not recovered_observations:
        print_progress(settings, budget, "sobol: no affordable initial batch")

    # -------------------------------------------------------------------------
    # 3. GlobalSUR / GP loop
    # -------------------------------------------------------------------------

    while not exhausted and budget.remaining_hqc > 0:
        use_fantasy_globalsur = (
            settings.fantasy_batching
            and settings.batching
            and settings.acquisition_function == "GlobalSUR"
        )

        if use_fantasy_globalsur:
            selected, exhausted = select_fantasy_globalsur_batch(
                strategy,
                settings,
                budget,
                observations,
                fantasy_rng,
                device=gp_device,
            )
        else:
            selected, exhausted = select_plain_aepsych_batch(
                strategy,
                settings,
                budget,
                device=gp_device,
            )

        if not selected:
            break

        stitched_total_gates = sum(config.n_gates for config in selected)

        print()
        print(
            f"[real batch] phase=globalsur "
            f"configs={len(selected)} "
            f"stitched_total_gates={stitched_total_gates}"
        )
        print()

        pending_step = checkpoint_step + 1
        selected_batch_metadata = sent_batch_metadata(settings, selected, data=data)
        write_pending_backend_batch(
            settings=settings,
            phase="globalsur",
            step=pending_step,
            sent_configs=selected,
            data=data,
            batch_metadata=selected_batch_metadata,
        )

        measure_batch_and_update_real_strategy(
            phase="globalsur",
            selected=selected,
            strategy=strategy,
            rmb=rmb,
            rng=rng,
            data=data,
            budget=budget,
            settings=settings,
            observations=observations,
            results_for_plot=results_for_plot,
            device=gp_device,
        )
        checkpoint_step += 1
        save_real_checkpoint(
            rmb=rmb,
            strategy=strategy,
            settings=settings,
            observations=observations,
            device=gp_device,
            phase="globalsur",
            step=checkpoint_step,
            budget=budget,
            sent_configs=selected,
            batch_metadata=selected_batch_metadata,
        )
        clear_pending_backend_batch(settings)

    # -------------------------------------------------------------------------
    # 4. Summary and contour extraction
    # -------------------------------------------------------------------------

    stop_reason = (
        "HQC budget exhausted"
        if exhausted or budget.remaining_hqc <= 0
        else "no affordable batch left"
    )

    print_experiment_summary(
        data,
        settings,
        budget,
        stop_reason=stop_reason,
    )

    plot_strategy = refreshed_strategy_for_prediction(
        strategy,
        settings,
        observations,
        device=gp_device,
    )

    crossings = level_set_configs(
        plot_strategy,
        settings,
        device=gp_device,
        measured_data=data,
        debug=False,
    )

    base_path = None
    if settings.save_path is not None:
        save_path = Path(settings.save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        rmb.save(save_path)
        base_path = resolve_data_path(save_path)
        write_hqc_metadata(base_path, settings, budget, "final", None)
        print(f"[saved] RMB data: {base_path}")

    if settings.save_gp_prediction_grid and base_path is not None:
        save_gp_prediction_grid(
            plot_strategy,
            settings,
            device=gp_device,
            json_path=base_path,
        )

    if return_budget:
        return rmb, crossings, budget
    return rmb, crossings


def run_with_budget(settings: FantasySettings):
    """Run the modified crossing experiment and return the final budget."""
    return run(settings, return_budget=True)


def main() -> None:
    kwargs = control_panel_settings_kwargs()
    qubits = qubit_slices(kwargs.pop("n_qubits"))
    multi_qubit = len(qubits) > 1

    for run_index, rng_seed in enumerate(RNG_SEEDS, start=1):
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        for qubit_index, n_qubits in enumerate(qubits, start=1):
            seed_kwargs = dict(kwargs)
            seed_kwargs["rng_seed"] = rng_seed
            seed_kwargs["n_qubits"] = n_qubits
            if seed_kwargs.get("recovery_mode"):
                if seed_kwargs.get("recovery_folder") is None:
                    raise ValueError(
                        "RECOVERY_FOLDER must be set when RECOVERY_MODE=True "
                        "so resumed outputs are written into that folder."
                    )
                seed_kwargs["save_path"] = timestamped_recovery_save_path(
                    seed_kwargs["recovery_folder"]
                )
            else:
                backend_folder = backend_folder_from_kwargs(seed_kwargs)
                seed_kwargs["save_path"] = timestamped_personal_save_path(
                    seed=rng_seed,
                    timestamp=timestamp,
                    n_qubits=n_qubits,
                    multi_qubit=multi_qubit,
                    backend_folder=backend_folder,
                )
            print(
                f"\n[seed run] {run_index}/{len(RNG_SEEDS)} "
                f"qubit slice {qubit_index}/{len(qubits)} "
                f"rng_seed={rng_seed} n_qubits={n_qubits} "
                f"save_path={seed_kwargs['save_path']}\n"
            )
            run(FantasySettings(**seed_kwargs))


if __name__ == "__main__":
    main()
