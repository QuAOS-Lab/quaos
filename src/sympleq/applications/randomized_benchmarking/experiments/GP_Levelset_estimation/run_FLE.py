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
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import torch

from aepsych.strategy import SequentialStrategy

from sympleq.applications.randomized_benchmarking.RMB import resolve_data_path
from sympleq.applications.randomized_benchmarking.experiments.common import (
    print_experiment_summary,
    print_progress,
    start_run,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.fantasy_levelset_estimation import (
    Observation,
    build_strategy,
    choose_gp_device,
    contour_target,
    level_set_configs,
    measure_batch_and_update_real_strategy,
    move_strategy_models_to_device,
    refreshed_strategy_for_prediction,
    seed_fake_corners,
    select_affordable_prefix,
    select_fantasy_globalsur_batch,
    select_plain_aepsych_batch,
    sobol_initial_candidates,
    temporary_torch_default_device,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.fantasy_levelset_settings import (
    FantasySettings,
    RNG_SEEDS,
    control_panel_settings_kwargs,
)


# Storing

def timestamped_personal_save_path(seed: int | None = None) -> Path:
    """Timestamped JSON output path under the repository's Personal folder."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if seed is None:
        return Path("Personal") / f"FLE_{timestamp}.json"
    return Path("Personal") / f"seed_{seed}" / f"FLE_{timestamp}.json"


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

    print("[fake anchors]")
    print(f"  use_fake_corners                     = {settings.use_fake_corners}")
    print(f"  easy_corner_outcome                  = {settings.easy_corner_outcome}")
    print(f"  hard_corner_outcome                  = {settings.hard_corner_outcome}")
    print(f"  extra_fake_anchors                   = {settings.extra_fake_anchors}")

    print("[sobol warm-up]")
    print(f"  initial_sobol_samples                = {settings.initial_sobol_samples}")
    print(f"  initial_sobol_max_cost_per_run       = "
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

    rng, rmb, budget = start_run(settings)
    backend_details = [
        f"{name}={value}"
        for name in ("device_name", "project_name")
        if (value := getattr(rmb.backend, name, None)) is not None
    ]
    backend_text = type(rmb.backend).__name__
    if backend_details:
        backend_text = f"{backend_text} ({', '.join(backend_details)})"
    print("[RMB backend]")
    print(f"  actual rmb.backend                   = {backend_text}")

    data = rmb._data

    fantasy_seed = None if settings.rng_seed is None else settings.rng_seed + 99173
    fantasy_rng = np.random.default_rng(fantasy_seed)

    strategy = build_strategy(settings)

    observations: list[Observation] = []
    results_for_plot: list[tuple[float, float, int]] = []

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
    sobol_batch, exhausted = select_affordable_prefix(
        sobol_candidates,
        settings,
        budget,
        max_cost_per_run=settings.initial_sobol_max_cost_per_run,
    )

    if sobol_batch:
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
        print(f"[saved] RMB data: {base_path}")

    if settings.save_gp_prediction_grid and base_path is not None:
        (
            probabilities,
            latent_mean,
            latent_variance,
            gates_grid,
            ratio_grid,
        ) = predict_native_level_set(
            plot_strategy,
            settings,
            device=gp_device,
        )
        if settings.save_gp_prediction_grid:
            grid_path = base_path.parent / f"{base_path.stem}_gp_grid.npz"
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

    if return_budget:
        return rmb, crossings, budget
    return rmb, crossings


def run_with_budget(settings: FantasySettings):
    """Run the modified crossing experiment and return the final budget."""
    return run(settings, return_budget=True)


def main() -> None:
    kwargs = control_panel_settings_kwargs()

    for run_index, rng_seed in enumerate(RNG_SEEDS, start=1):
        seed_kwargs = dict(kwargs)
        seed_kwargs["rng_seed"] = rng_seed
        seed_kwargs["save_path"] = timestamped_personal_save_path(seed=rng_seed)
        print(
            f"\n[seed run] {run_index}/{len(RNG_SEEDS)} "
            f"rng_seed={rng_seed} save_path={seed_kwargs['save_path']}\n"
        )
        run(FantasySettings(**seed_kwargs))


if __name__ == "__main__":
    main()
