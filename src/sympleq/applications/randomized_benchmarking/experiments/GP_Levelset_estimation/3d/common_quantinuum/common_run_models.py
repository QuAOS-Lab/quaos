"""
Common run sript for FLE and Cost_Aware models.
Run with python path/to/common_run_models.py "FLE" or "COST_AWARE" to run the corresponding model.
FLE is implemened.

For Cost_aware:
Should define a settings class as FLE_settings and definition as FLE_3d_fix_qubit_band
       use write_hqc_metadata for storing intermediate results in json


The GP grid saving is not implemented for Cost_Aware
The save_real_checkpoint saves only the real RMB data in a json file and not the GP grid for Cost_Aware
Saves in "Path("Personal") / model_folder / seed_folder / f"FLE_{timestamp}"

run_FLE runs only FLE; the storing of configs and data/grid is done through this after each *real* measurement

IMPORTANT: The 'main' function checks if the settings allows for more than 7000 gates for the emulator,
However, THE ACTUAL CHECK IS PASSED IN THE DEFINITION FILE (FLE_3d_fix_qubit_band.py)

(SEE SELECT_AFFORDABLE_PREFIX() AND SELECT_FANTASY_GLOBALSUR_BATCH();
if settings.backend_factory is quantinuum_emulator_backend_factory:
    if stitched_total_gates > gate_budget:
        break
This should be done as soon as the stitching is done.

"""

import logging
import json
import re
import warnings
from datetime import datetime
from pathlib import Path
import numpy as np


from sympleq.applications.randomized_benchmarking.RMB import RMB, RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.experiments.common import (
    batch_hqc_cost,
    default_backend_factory,
    quantinuum_H2_backend_factory,
    quantinuum_H21E_backend_factory,
    quantinuum_emulator_backend_factory,
    print_experiment_summary,
    print_progress,
    print_progress,
    start_run,
)


MODEL = "FLE"


if MODEL == "FLE":
    import torch
    from FLE_settings import (
        FantasySettings as SettingsClass,
        RNG_SEEDS,
        QUBIT_BAND_LENGTHS,
        QUBIT_BAND_LENGTHS,
        control_panel_settings_kwargs as settings_kwargs
    )

    from FLE_3d_fix_qubit_band import (
        Observation,
        add_observation_to_strategy,
        add_observation_to_strategy,
        build_strategy,
        choose_gp_device,
        contour_target,
        measure_batch_and_update_real_strategy,
        move_strategy_models_to_device,
        one_shot_requests,
        point_from_config,
        repair_stats,
        refreshed_strategy_for_prediction,
        reset_repair_stats,
        seed_fake_corners,
        select_affordable_prefix,
        select_fantasy_globalsur_batch,
        sobol_initial_candidates,
        temporary_torch_default_device,
        valid_config,
    )


elif MODEL == "COST_AWARE":
    from cost_aware_surface_settings import CostAwareSettings as SettingsClass

else:
    raise ValueError(f"Unknown model: {MODEL}")


# Storing
def timestamped_personal_save_path(seed: int | None = None, qubit_band_length: int | None = None) -> Path:
    """Timestamped run folder and final JSON path under ``Personal``."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    seed_folder = "seed_unseeded" if seed is None else f"seed_{seed}"
    band_folder = f"qband_{qubit_band_length}" if qubit_band_length is not None else "qband_unseeded"
    backend_folder = "H2_1"
    model_folder = "CostAware" if MODEL == "COST_AWARE" else "FLE"

    run_folder = Path("Personal") / model_folder / backend_folder / band_folder / seed_folder / f"FLE_{timestamp}"

    return run_folder / f"FLE_{timestamp}.json"


def timestamped_recovery_save_path(recovery_folder: str | Path) -> Path:
    """Timestamped final JSON path inside an existing recovery folder."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path(recovery_folder) / f"FLE_recovery_{timestamp}.json"


def timestamped_recovery_save_path(recovery_folder: str | Path) -> Path:
    """Timestamped final JSON path inside an existing recovery folder."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path(recovery_folder) / f"FLE_recovery_{timestamp}.json"


# Bookkeeping_
def write_hqc_metadata(
    json_path: Path,
    settings,
    budget,
    phase: str,
    step: int | None,
    sent_configs=None,
) -> None:
    """Store experiment metadata inside an RMB JSON checkpoint."""

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    payload["experiment"] = {
        "phase": phase,
        "step": step,
        "spent_hqc": float(settings.hqc_budget - budget.remaining_hqc),
        "remaining_hqc": float(budget.remaining_hqc),
        "hqc_budget": float(settings.hqc_budget),
        "repair_stats": repair_stats(),
    }
    if sent_configs is not None:
        sent_configs = sorted(
            sent_configs,
            key=lambda config: config.n_qubits,
            reverse=True,
        )
        sent_configs = sorted(
            sent_configs,
            key=lambda config: config.n_qubits,
            reverse=True,
        )
        payload["experiment"]["sent_configs"] = [
            {
                "n_1qb_gates": int(config.n_1qb_gates),
                "n_2qb_gates": int(config.n_2qb_gates),
                "n_qubits": int(config.n_qubits),
                "n_gates": int(config.n_gates),
                "ratio_2_qb_gates": float(config.ratio_2_qb_gates),
            }
            for config in sent_configs
        ]
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _settings_run_folder(settings) -> Path | None:
    if settings.save_path is not None:
        return Path(settings.save_path).parent
    if settings.recovery_folder is not None:
        return Path(settings.recovery_folder)
    return None


def pending_backend_batch_path(settings) -> Path | None:
    folder = _settings_run_folder(settings)
    if folder is None:
        return None
    return folder / "pending_backend_batch.json"


def write_pending_backend_batch(
    *,
    settings,
    phase: str,
    step: int,
    sent_configs,
    data=None,
) -> None:
    """Persist the exact FLE batch before submitting stitched circuits."""

    if not settings.save_real_checkpoints:
        return

    path = pending_backend_batch_path(settings)
    if path is None:
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    requests = []
    for config in sent_configs:
        first_shot_index = None
        if data is not None:
            estimator = data.get(config) if hasattr(data, "get") else None
            if estimator is not None and hasattr(estimator, "num_runs"):
                first_shot_index = int(estimator.num_runs())
        requests.append(
            {
                "config": {
                    "n_1qb_gates": int(config.n_1qb_gates),
                    "n_2qb_gates": int(config.n_2qb_gates),
                    "n_qubits": int(config.n_qubits),
                    "n_gates": int(config.n_gates),
                    "ratio_2_qb_gates": float(config.ratio_2_qb_gates),
                    "random_elimination": float(config.random_elimination),
                    "use_scrambler": bool(config.use_scrambler),
                },
                "requested_shots": 1,
                "first_shot_index": first_shot_index,
            }
        )

    payload = {
        "reason": "before_backend_submit",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "model": MODEL,
        "phase": phase,
        "step": int(step),
        "pending_batch": {
            "requests": requests,
        },
    }
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[saved] pending backend batch: {path}")


def clear_pending_backend_batch(settings) -> None:
    path = pending_backend_batch_path(settings)
    if path is not None and path.exists():
        path.unlink()


def predict_native_level_set(
    strategy,
    settings,
    *,
    device,
    n_gates_grid: int = 80,
    n_ratio_grid: int = 80,
    n_qubits_grid: int = 16,
):
    gates_axis = np.geomspace(
        max(1.0, float(settings.n_gates_bounds[0])),
        float(settings.n_gates_bounds[1]),
        n_gates_grid,
    )
    ratio_axis = np.linspace(
        float(settings.ratio_bounds[0]),
        float(settings.ratio_bounds[1]),
        n_ratio_grid,
    )
    qubits_axis = np.linspace(
        float(settings.n_qubits_bounds[0]),
        float(settings.n_qubits_bounds[1]),
        n_qubits_grid,
    )

    qubits_grid, ratio_grid, gates_grid = np.meshgrid(
        qubits_axis,
        ratio_axis,
        gates_axis,
        indexing="ij",
    )

    grid = torch.tensor(
        np.column_stack(
            [
                gates_grid.ravel(),
                ratio_grid.ravel(),
                qubits_grid.ravel(),
            ]
        ),
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

    shape = gates_grid.shape

    return (
        probabilities.detach().cpu().numpy().reshape(shape),
        latent_mean.detach().cpu().numpy().reshape(shape),
        latent_variance.detach().cpu().numpy().reshape(shape),
        gates_grid,
        ratio_grid,
        qubits_grid,
        gates_axis,
        ratio_axis,
        qubits_axis,
    )


def is_h2_quantinuum_backend(settings) -> bool:
    return settings.backend_factory in {
        quantinuum_H2_backend_factory,
        quantinuum_H21E_backend_factory,
    }


def select_h2_sobol_batch(candidates, settings, budget) -> tuple[list, bool]:
    """Select one cross-qubit Sobol batch for H2, capped by Sobol HQC."""

    selected = []
    cost_cap = (
        settings.max_cost_per_run
        if settings.initial_sobol_max_cost_per_run is None
        else settings.initial_sobol_max_cost_per_run
    )

    for candidate in candidates:
        trial_batch = selected + [candidate]
        next_cost = batch_hqc_cost(one_shot_requests(trial_batch))

        if next_cost > cost_cap:
            break
        if next_cost > budget.remaining_hqc:
            return selected, True

        selected.append(candidate)
        print(f"Stitched total gates: {sum(config.n_gates for config in trial_batch)}")

    return selected, False


def save_gp_prediction_grid(strategy,
                            *,
                            model: str,
                            settings: SettingsClass,
                            device: torch.device,
                            json_path: Path,
                            ) -> Path:
    """Save the current 3D GP prediction grid beside one RMB JSON snapshot."""

    (
        probabilities,
        latent_mean,
        latent_variance,
        gates_grid,
        ratio_grid,
        qubits_grid,
        gates_axis,
        ratio_axis,
        qubits_axis,
    ) = predict_native_level_set(strategy, settings, device=device)

    grid_path = json_path.parent / f"{json_path.stem}_gp_grid_3d.npz"
    np.savez_compressed(
        grid_path,
        gates_grid=gates_grid,
        ratio_grid=ratio_grid,
        qubits_grid=qubits_grid,
        gates_axis=gates_axis,
        ratio_axis=ratio_axis,
        qubits_axis=qubits_axis,
        x_grid=gates_grid,
        y_grid=ratio_grid,
        z_grid=qubits_grid,
        probabilities=probabilities,
        latent_mean=latent_mean,
        latent_variance=latent_variance,
        target=np.asarray(contour_target(settings), dtype=float),
        coordinate_system=np.asarray("total_ratio_n_qubits_3d"),
        rng_seed=np.asarray(
            -1 if settings.rng_seed is None else settings.rng_seed,
            dtype=int,
        ),
    )
    print(f"[saved] 3D GP grid: {grid_path}")
    return grid_path


def save_real_checkpoint(
    *,
    rmb,
    model: str,
    settings,
    phase: str,
    step: int,
    budget,
    sent_configs=None,
    strategy=None,
    observations=None,
    device=None,
) -> None:
    """Save real RMB data and a GP grid after one real measurement batch."""

    if not settings.save_real_checkpoints or settings.save_path is None:
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    checkpoint_path = (
        Path(settings.save_path).parent / f"measurement_{step:03d}_{phase}_{timestamp}.json"
    )
    rmb.save(checkpoint_path)
    checkpoint_path = resolve_data_path(checkpoint_path)
    write_hqc_metadata(checkpoint_path, settings, budget, phase, step, sent_configs)
    print(f"[saved] real RMB checkpoint: {checkpoint_path}")

    if not settings.save_gp_prediction_grid:
        return

    if model == "FLE":
        checkpoint_strategy = refreshed_strategy_for_prediction(
            strategy,
            settings,
            observations,
            device=device,
        )
        save_gp_prediction_grid(checkpoint_strategy,
                                settings=settings,
                                model=model,
                                device=device,
                                json_path=checkpoint_path,
                                )
    elif model == "COST_AWARE":  # TODO add grid saving for Cost-Aware
        warnings.warn(
            "GP grid saving not implemented for Cost-Aware model. "
            "Skipping GP grid save."
        )


def latest_recovery_json(settings) -> Path | None:
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


def recovery_phase(json_path: Path) -> str | None:
    """Read the checkpoint phase from recovery metadata if present."""

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    experiment = payload.get("experiment", {})
    phase = experiment.get("phase")
    return None if phase is None else str(phase)


def read_recovered_sobol_submissions(json_path: Path, max_sobol_submissions: int) -> int:
    """Return how many Sobol submissions are already in the recovery checkpoint."""

    phase = recovery_phase(json_path)
    if phase is None:
        return 0
    if phase.lower().startswith("globalsur"):
        return max_sobol_submissions
    match = re.fullmatch(r"sobol_(\d+)", phase)
    if match:
        return min(int(match.group(1)), max_sobol_submissions)
    return max_sobol_submissions


def recover_observations_from_json(
    json_path: Path,
    *,
    strategy,
    rmb,
    observations,
    results_for_plot,
    device,
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


def latest_recovery_json(settings) -> Path | None:
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


def recovery_phase(json_path: Path) -> str | None:
    """Read the checkpoint phase from recovery metadata if present."""

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    experiment = payload.get("experiment", {})
    phase = experiment.get("phase")
    return None if phase is None else str(phase)


def read_recovered_sobol_submissions(json_path: Path, max_sobol_submissions: int) -> int:
    """Return how many Sobol submissions are already in the recovery checkpoint."""

    phase = recovery_phase(json_path)
    if phase is None:
        return 0
    if phase.lower().startswith("globalsur"):
        return max_sobol_submissions
    match = re.fullmatch(r"sobol_(\d+)", phase)
    if match:
        return min(int(match.group(1)), max_sobol_submissions)
    return max_sobol_submissions


def recover_observations_from_json(
    json_path: Path,
    *,
    strategy,
    rmb,
    observations,
    results_for_plot,
    device,
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


def print_run_handles(
    settings: SettingsClass,
    *,
    gp_device: None,
) -> None:
    """
    Print the tunable handles for reproducibility.
    """

    print("\n================ RUN HANDLES ================")
    print(f"  model                                = {MODEL}")

    if MODEL == 'FLE':
        #gp_device = torch.device
        print("[target]")
        print(f"  target_threshold                     = {settings.target_threshold}")

        print("[qubits]")
        print(f"  qubits                     = {settings.n_qubits}")
        print(f"  qubit_band_length                     = {settings.qubit_band_length}")

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

    else:
        raise NotImplementedError(f"Run handles not implemented for model: {MODEL}")


def run_FLE(
    settings: SettingsClass,
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
    recovered_sobol_count = 0
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
                budget.remaining_hqc = max(0.0, float(settings.hqc_budget) - recovered_hqc_spent)
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
                recovered_sobol_count = read_recovered_sobol_submissions(
                    recovery_json,
                    settings.initial_sobol_submissions,
                )
                print(
                    f"[recovery] loaded {recovered_observations} observations "
                    f"from {recovery_json}; recovered "
                    f"{recovered_sobol_count} Sobol submissions."
                )
            else:
                print(
                    f"[recovery] found {recovery_json}, but it contained no "
                    "observations; running Sobol."
                )

    recovered_observations = 0
    recovered_sobol_count = 0
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
                budget.remaining_hqc = max(0.0, float(settings.hqc_budget) - recovered_hqc_spent)
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
                recovered_sobol_count = read_recovered_sobol_submissions(
                    recovery_json,
                    settings.initial_sobol_submissions,
                )
                print(
                    f"[recovery] loaded {recovered_observations} observations "
                    f"from {recovery_json}; recovered "
                    f"{recovered_sobol_count} Sobol submissions."
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

    sobol_candidates = [] if recovered_observations else sobol_initial_candidates(settings)
    sobol_submissions = 0

    max_sobol_submissions = settings.initial_sobol_submissions
    sobol_candidates = sobol_initial_candidates(settings)
    if recovered_observations and recovered_sobol_count < max_sobol_submissions:
        recovered_sobol_qubits = {
            config.n_qubits for config in rmb._data
        }
        print(
            "[recovery] skipping recovered Sobol qubits: "
            f"{sorted(recovered_sobol_qubits)}"
        )
        sobol_candidates = [
            config for config in sobol_candidates
            if config.n_qubits not in recovered_sobol_qubits
        ]
    sobol_submissions = recovered_sobol_count
    remaining_sobol_candidates = list(sobol_candidates)

    while (remaining_sobol_candidates and sobol_submissions < max_sobol_submissions and not exhausted
           and budget.remaining_hqc > 0
           ):
        print(f"[sobol configs] submission={sobol_submissions + 1}")

        if is_h2_quantinuum_backend(settings):
            sobol_batch, exhausted = select_h2_sobol_batch(
                remaining_sobol_candidates,
                settings,
                budget,
            )
        else:
            sobol_batch, exhausted = select_affordable_prefix(
                remaining_sobol_candidates,
                settings,
                budget,
                max_cost_per_run=settings.initial_sobol_max_cost_per_run,
            )

        if not sobol_batch:
            break

        # Remove the configs that are about to be measured.
        remaining_sobol_candidates = remaining_sobol_candidates[len(sobol_batch):]

        for index, config in enumerate(sobol_batch, start=1):
            is_valid = valid_config(config)
            print(
                f"  config={index:03d}: "
                f"n_gates={config.n_gates} "
                f"n_1q={config.n_1qb_gates} "
                f"n_2q={config.n_2qb_gates} "
                f"ratio={config.ratio_2_qb_gates:.4f} "
                f"n_qubits={config.n_qubits} "
                f"valid={is_valid}"
            )

            if not is_valid:
                print(
                    "[sobol invalid sent config] "
                    f"n_1q={config.n_1qb_gates} "
                    f"n_2q={config.n_2qb_gates} "
                    f"n_gates={config.n_gates} "
                    f"ratio={config.ratio_2_qb_gates:.4f} "
                    f"n_qubits={config.n_qubits}"
                )

        pending_step = checkpoint_step + 1
        pending_phase = f"sobol_{sobol_submissions + 1:03d}"
        write_pending_backend_batch(
            settings=settings,
            phase=pending_phase,
            step=pending_step,
            sent_configs=sobol_batch,
            data=data,
        )

        measure_batch_and_update_real_strategy(
            phase=f"sobol_{sobol_submissions + 1}",
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

        sobol_submissions += 1
        checkpoint_step += 1

        save_real_checkpoint(model=model,
                             rmb=rmb,
                             strategy=strategy,
                             settings=settings,
                             observations=observations,
                             device=gp_device,
                             phase=f"sobol_{sobol_submissions:03d}",
                             step=checkpoint_step,
                             budget=budget,
                             sent_configs=sobol_batch,
                             )
        clear_pending_backend_batch(settings)

    new_sobol_submissions = sobol_submissions - recovered_sobol_count
    if new_sobol_submissions:
        print(
            f"sobol: completed {new_sobol_submissions} submissions after recovery, "
            f"{len(sobol_candidates) - len(remaining_sobol_candidates)} configs measured"
        )
    else:
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
        write_pending_backend_batch(
            settings=settings,
            phase="globalsur",
            step=pending_step,
            sent_configs=selected,
            data=data,
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
            model=model,
            rmb=rmb,
            strategy=strategy,
            settings=settings,
            observations=observations,
            device=gp_device,
            phase="globalsur",
            step=checkpoint_step,
            budget=budget,
            sent_configs=selected,
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
            model=model,
            device=gp_device,
            json_path=base_path,
        )

    if return_budget:
        return rmb, budget
    return rmb


def main(model, SettingsClass, settings_kwargs) -> None:
    kwargs = settings_kwargs()
    backend_factory = kwargs.get("backend_factory", default_backend_factory)

    # THIS JUST CHECKS IF THE SETTINGS ARE VALID;
    # THE ACTUAL CHECK IS PASSED IN THE DEFINITION FILE (FLE_3d_fix_qubit_band.py)
    # SEE SELECT_AFFORDABLE_PREFIX() AND SELECT_FANTASY_GLOBALSUR_BATCH();
    # if settings.backend_factory is quantinuum_emulator_backend_factory:
        # if stitched_total_gates > gate_budget:
        #     break

    if backend_factory == quantinuum_emulator_backend_factory and kwargs.get("gate_budget") > 7000:
        raise ValueError(
            "Gate budget too high for Quantinuum emulator. Please set gate_budget <= 7000."
        )
   # TODO if backend_factory == ACTUAL DEVICE and kwargs.get("max_cost_per_run") > 35:
        # raise ValueError(
        #     "HQC per stitched circuit too high  for Quantinuum. Please set max_cost_per_run <= 35."
        # )

    if model == "FLE":
        for qubit_band_length in QUBIT_BAND_LENGTHS:
            for run_index, rng_seed in enumerate(RNG_SEEDS, start=1):
                seed_kwargs = dict(kwargs)
                seed_kwargs["rng_seed"] = rng_seed
                seed_kwargs["qubit_band_length"] = qubit_band_length
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
                    seed_kwargs["save_path"] = timestamped_personal_save_path(
                        seed=rng_seed, qubit_band_length=qubit_band_length)
                print(
                    f"\n[seed run] {run_index}/{len(RNG_SEEDS)} "
                    f"rng_seed={rng_seed} save_path={seed_kwargs['save_path']}\n"
                )
                run_FLE(SettingsClass(**seed_kwargs))

    elif model == "COST_AWARE":
        raise NotImplementedError("Cost-Aware model is not implemented in this script.")


if __name__ == "__main__":
    model = MODEL
    main(model, SettingsClass, settings_kwargs)
