""" 
Common run sript for FLE and Cost_Aware models.
Run with python path/to/common_run_models.py "FLE" or "COST_AWARE" to run the corresponding model.
FLE is implemened.

For Cost_aware:
Should define a settings class as FLE_settings and definition as FLE_3d_fix_qubit_band
       use write_hqc_metadata for storing intermediate results in json
       

The GP grid saving is not implemented for Cost_Aware
The save_real_checkpoint saves only the real RMB data in a json file and not the GP grid for Cost_Aware
Saves generated data under scripts/personal/randomized_benchmarking_personal/Personal.

run_FLE runs only FLE; the storing of configs and data/grid is done through this after each *real* measurement

IMPORTANT: The 'main' function checks if the settings allows for more than 7000 gates for the emulator, 
However, THE ACTUAL CHECK IS PASSED IN THE DEFINITION FILE (FLE_3d_fix_qubit_band.py)

(SEE SELECT_AFFORDABLE_PREFIX() AND SELECT_FANTASY_GLOBALSUR_BATCH();
if settings.backend_factory is quantinuum_emulator_backend_factory:
    if stitched_total_gates > gate_budget:
        break 
This should be done as soon as the stitching is done.

"""
from __future__ import annotations

import logging
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path
import numpy as np


from sympleq.applications.randomized_benchmarking.RMB import resolve_data_path
from sympleq.applications.randomized_benchmarking.experiments.common import (
    default_backend_factory,
    quantinuum_emulator_backend_factory,
    print_experiment_summary,
    start_run,
)


_SUPPORTED_MODELS = ("FLE", "COST_AWARE")
_DEFAULT_MODEL = "COST_AWARE"
GENERATED_PERSONAL_ROOT = (
    Path("scripts")
    / "personal"
    / "randomized_benchmarking_personal"
    / "Personal"
)


def _resolve_model(argv: list[str]) -> str:
    """Pick the model from the CLI, e.g. ``python common_run_models.py FLE``.

    Falls back to ``_DEFAULT_MODEL`` when no argument is given (so imports still
    work), and errors clearly on an unrecognised value rather than silently
    running the default.
    """
    if len(argv) <= 1:
        return _DEFAULT_MODEL
    requested = argv[1].strip().upper()
    if requested not in _SUPPORTED_MODELS:
        raise SystemExit(
            f"Unknown model {argv[1]!r}; choose one of {', '.join(_SUPPORTED_MODELS)}."
        )
    return requested


# Model is taken from the command line (the conditional imports below depend on
# it, so it must be resolved at import time).
MODEL = _resolve_model(sys.argv)



if MODEL == "FLE":
    import torch
    from FLE_settings import (
        FantasySettings as SettingsClass,
        RNG_SEEDS,
        control_panel_settings_kwargs as settings_kwargs
    )

    from FLE_3d_fix_qubit_band import (
        Observation,
        build_strategy,
        choose_gp_device,
        contour_target,
        measure_batch_and_update_real_strategy,
        move_strategy_models_to_device,
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
    from cost_aware_settings import (
        CostAwareSettings as SettingsClass,
        RNG_SEEDS,
        control_panel_settings_kwargs as settings_kwargs,
        with_randomized_initial_rate_guesses,
    )
    from cost_aware_surface_design import run_with_design_plots

else:
    raise ValueError(f"Unknown model: {MODEL}")


# Storing
def timestamped_personal_save_path(seed: int | None = None) -> Path:
    """Timestamped generated run folder and final JSON path."""

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    seed_folder = "seed_unseeded" if seed is None else f"seed_{seed}"
    model_folder = "CostAware" if MODEL == "COST_AWARE" else "FLE"
    prefix = "CostAware" if MODEL == "COST_AWARE" else "FLE"
    run_folder = GENERATED_PERSONAL_ROOT / model_folder / seed_folder / f"{prefix}_{timestamp}"

    return run_folder / f"{prefix}_{timestamp}.json"


def restartable_cost_aware_save_path(seed: int | None = None) -> Path:
    seed_folder = "seed_unseeded" if seed is None else f"seed_{seed}"
    return (
        GENERATED_PERSONAL_ROOT
        / "CostAware"
        / seed_folder
        / "CostAware_restartable.json"
    )


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
    }
    if "repair_stats" in globals():
        payload["experiment"]["repair_stats"] = repair_stats()
    if sent_configs is not None:
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

    elif MODEL == "COST_AWARE":
        print("[surface]")
        print(f"  q_values                             = {settings.q_values}")
        print(f"  n_gates_bounds                       = {settings.n_gates_bounds}")
        print(f"  ratio_bounds                         = {settings.ratio_bounds}")

        print("[batching / budget]")
        print(f"  hqc_budget                           = {settings.hqc_budget}")
        print(f"  max_cost_per_run                     = {settings.max_cost_per_run}")
        print(f"  gate_budget                          = {settings.gate_budget}")
        if hasattr(settings, "single_batch_gate_budget"):
            print(f"  single_batch_gate_budget             = {settings.single_batch_gate_budget}")
        if hasattr(settings, "emulator_max_batch_cost"):
            print(f"  emulator_max_batch_cost              = {settings.emulator_max_batch_cost}")
        print(f"  max_qubit_window                     = {settings.max_qubit_window}")

        print("[acquisition]")
        print(f"  acquisition_q_resolution             = {settings.acquisition_q_resolution}")
        print(f"  acquisition_ratio_points             = {settings.acquisition_ratio_points}")
        print(f"  grid_resolution                      = {settings.grid_resolution}")
        print(f"  boundary_fit_resolution              = {settings.boundary_fit_resolution}")

        print("[prior]")
        print(f"  initial_one_q_pauli_error            = {settings.initial_one_q_pauli_error:.6g}")
        print(f"  initial_two_q_pauli_error            = {settings.initial_two_q_pauli_error:.6g}")
        one_q_prefactor = (
            settings.initial_one_q_pauli_error
            / settings.initial_one_q_pauli_error_base
            if settings.initial_one_q_pauli_error_base
            else float("nan")
        )
        two_q_prefactor = (
            settings.initial_two_q_pauli_error
            / settings.initial_two_q_pauli_error_base
            if settings.initial_two_q_pauli_error_base
            else float("nan")
        )
        print(
            "  initial_rate_random_prefactors       = "
            f"({one_q_prefactor:.6g}, {two_q_prefactor:.6g})"
        )
        print(f"  initial_error_relative_uncertainty   = {settings.initial_error_relative_uncertainty}")
        if settings.initial_rate_randomization_enabled:
            print(f"  initial_rate_random_seed             = {settings.initial_rate_random_seed}")
            print(
                "  initial_rate_relative_stds           = "
                f"({settings.initial_one_q_random_relative_std:.4g}, "
                f"{settings.initial_two_q_random_relative_std:.4g})"
            )
            print(
                "  initial_rate_relative_deltas         = "
                f"({settings.initial_one_q_random_relative_delta:+.4g}, "
                f"{settings.initial_two_q_random_relative_delta:+.4g})"
            )

        backend_factory_name = getattr(
            settings.backend_factory,
            "__name__",
            type(settings.backend_factory).__name__,
        )
        print("[backend]")
        print(f"  backend_model                        = {settings.backend_model}")
        print(f"  backend_factory                      = {backend_factory_name}")

        print("[debug]")
        print(f"  rng_seed                             = {settings.rng_seed}")
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

    # -------------------------------------------------------------------------
    # 2. Sobol warm-up batch
    # -------------------------------------------------------------------------

    exhausted = False

    sobol_candidates = sobol_initial_candidates(settings)
    sobol_submissions = 0

    max_sobol_submissions = settings.initial_sobol_submissions
    remaining_sobol_candidates = list(sobol_candidates)

    while (remaining_sobol_candidates and sobol_submissions < max_sobol_submissions and not exhausted
           and budget.remaining_hqc > 0
           ):
        print(f"[sobol configs] submission={sobol_submissions + 1}")

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

    if sobol_submissions:
        print(
            f"sobol: completed {sobol_submissions} submissions, "
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
            device=gp_device,
            json_path=base_path,
        )

    if return_budget:
        return rmb, budget
    return rmb


def run_COST_AWARE(
    settings: SettingsClass,
    *,
    return_budget: bool = False,
):
    logging.getLogger().setLevel(logging.WARNING)
    print_run_handles(settings, gp_device=None)
    rmb, configs, budget = run_with_design_plots(settings)
    backend_details = [
        f"{name}={value}"
        for name in ("device_name", "project_name")
        if (value := getattr(rmb.backend, name, None)) is not None
    ]
    backend_text = type(rmb.backend).__name__
    if backend_details:
        backend_text = f"{backend_text} ({', '.join(backend_details)})"
    print(f"  actual rmb.backend                   = {backend_text}")
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

    gate_budget = kwargs.get("gate_budget")
    if (
        backend_factory == quantinuum_emulator_backend_factory
        and gate_budget is not None
        and int(gate_budget) > 7000
    ):
        raise ValueError(
            "Gate budget too high for Quantinuum emulator. Please set gate_budget <= 7000."
        )
   # TODO if backend_factory == ACTUAL DEVICE and kwargs.get("max_cost_per_run") > 35:
        # raise ValueError(
        #     "HQC per stitched circuit too high  for Quantinuum. Please set max_cost_per_run <= 35."
        # )

    if model == "FLE":
        for run_index, rng_seed in enumerate(RNG_SEEDS, start=1):
            seed_kwargs = dict(kwargs)
            seed_kwargs["rng_seed"] = rng_seed
            seed_kwargs["save_path"] = timestamped_personal_save_path(seed=rng_seed)
            print(
                f"\n[seed run] {run_index}/{len(RNG_SEEDS)} "
                f"rng_seed={rng_seed} save_path={seed_kwargs['save_path']}\n"
            )
            run_FLE(SettingsClass(**seed_kwargs))

    elif model == "COST_AWARE":
        for run_index, rng_seed in enumerate(RNG_SEEDS, start=1):
            seed_kwargs = dict(kwargs)
            seed_kwargs["rng_seed"] = rng_seed
            seed_kwargs = with_randomized_initial_rate_guesses(seed_kwargs, rng_seed)
            if seed_kwargs.get("save_path") is None:
                seed_kwargs["save_path"] = restartable_cost_aware_save_path(seed=rng_seed)
            print(
                f"\n[seed run] {run_index}/{len(RNG_SEEDS)} "
                f"rng_seed={rng_seed} save_path={seed_kwargs['save_path']}\n"
            )
            run_COST_AWARE(SettingsClass(**seed_kwargs))


if __name__ == "__main__":
    model = MODEL
    main(model, SettingsClass, settings_kwargs)
