"""
Fantasy-batched GPU-assisted GP level-set estimation for RMB.

What this file does
-------------------
Estimate the P(success) = target_threshold contour in:

    (total gates, two-qubit gate ratio)

space.

Important split
---------------
RMB backend:
    CPU only. This file does not try to GPU-enable RMB circuit execution.

AEPsych / GP / GlobalSUR:
    Uses CUDA if available and if USE_GPU=True in main().

Algorithm
---------
1. Add fake anchor points:
       easy corner -> assumed success
       hard corner -> assumed failure

2. Generate explicit Sobol warm-up points.

3. Measure Sobol batch using the CPU RMB backend.

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

import copy
import importlib
import json
import logging
import re
import sys
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import NamedTuple

import numpy as np
import torch
from aepsych.config import Config
from aepsych.strategy import SequentialStrategy
import aepsych.transforms.parameters as aepsych_parameter_transforms

from sympleq.applications.randomized_benchmarking.RMB import RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    MeasurementRequest,
    batch_hqc_cost,
    candidate_axes,
    config_points,
    print_experiment_summary,
    print_progress,
    spend_request_batch,
    start_run,
)
from sympleq.applications.randomized_benchmarking.experiments.scores import (
    gp_grid_scores,
)


_ORIGINAL_AEPSYCH_TRANSFORM_OPTIONS = aepsych_parameter_transforms.transform_options
_ORIGINAL_AEPSYCH_STR_TO_LIST = Config._str_to_list
_ORIGINAL_AEPSYCH_STR_TO_ARRAY = Config._str_to_array


def timestamped_personal_save_path(seed: int | None = None) -> Path:
    """Timestamped JSON output path under the repository's Personal folder."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if seed is None:
        return Path("Personal") / f"rick_fantasy_gpu_crossing_{timestamp}.json"
    return Path("Personal") / f"seed_{seed}" / f"rick_fantasy_gpu_crossing_{timestamp}.json"


def import_experiment_module(module_name: str):
    """Import an experiment helper module by short or fully qualified name."""
    if "." not in module_name:
        module_name = (
            "sympleq.applications.randomized_benchmarking.experiments."
            f"{module_name}"
        )
    return importlib.import_module(module_name)


def apply_bayesian_estimation_module(
    module_name: str,
    *,
    default_threshold: float,
    default_min_runs: int,
) -> None:
    """Point SympleQ modules at Rick's local estimator default.

    The estimator implementation comes from ``module_name``.  Only
    ``BayesianEstimator.default()`` is made local to this experiment, so shared
    code such as Charlie keeps its own default behavior.
    """
    module = importlib.import_module(module_name)
    base_estimator = module.BayesianEstimator

    class RickBayesianEstimator(base_estimator):
        @classmethod
        def default(cls):
            return cls(
                threshold=default_threshold,
                min_runs=default_min_runs,
            )

    RickBayesianEstimator.__name__ = base_estimator.__name__
    RickBayesianEstimator.__qualname__ = base_estimator.__qualname__
    RickBayesianEstimator.__module__ = base_estimator.__module__

    for loaded in list(sys.modules.values()):
        loaded_name = getattr(loaded, "__name__", "")
        if loaded_name.startswith("sympleq.") and hasattr(loaded, "BayesianEstimator"):
            setattr(loaded, "BayesianEstimator", RickBayesianEstimator)


# =============================================================================
# Data containers
# =============================================================================


class Observation(NamedTuple):
    """
    One binary observation given to AEPsych.

    x_cpu:
        Shape (1, 2). Stored on CPU so it can safely be replayed into a
        copied/rebuilt strategy.

    y:
        Binary outcome. 1 = success, 0 = failure.

    source:
        "fake", "sobol", "globalsur", or "fantasy".

    Important:
        fantasy observations must only be added to copied/fantasy strategies.
        They must never be added to the real strategy.
    """

    x_cpu: torch.Tensor
    y: int
    source: str


@dataclass(frozen=True)
class RickFantasyGPUSettings(CrossingSettings):
    """
    Settings container.

    Most user-tunable values are set in main() at the bottom of the file.

    This inherits the usual CrossingSettings handles, including things like:
        n_gates_bounds
        ratio_bounds
        hqc_budget
        max_cost_per_run
        rng_seed
        plot
        make_config(...)
    """

    # Output
    save_path: str | Path | None = field(default_factory=timestamped_personal_save_path)

    # Target contour
    target_threshold: float = 0.5

    # Fake anchors
    use_fake_corners: bool = True
    easy_corner_outcome: int = 1
    hard_corner_outcome: int = 0
    extra_fake_anchors: list[tuple[float, float, int]] = field(default_factory=list)

    # Explicit Sobol warm-up
    initial_sobol_samples: int = 10
    initial_sobol_max_cost_per_run: float | None = None
    sobol_scramble: bool = False

    # Validation pass on the trained-GP contour
    validate_gp_contour: bool = False
    reserve_validation_budget: bool = True
    validation_hqc_budget: float = 0.0
    validation_max_cost_per_run: float | None = None
    validation_target_configs: int = 30
    validation_min_configs: int = 15
    validation_max_restarts: int = 20
    validation_shots_per_config: int = 2
    validation_max_attempts_multiplier: int = 10
    validation_probability_band: tuple[float, float] = (0.4, 0.6)
    validation_seed_offset: int = 271828
    validation_save_path: str | Path | None = None
    plot_gp_level_set_result: bool = False
    save_gp_prediction_grid: bool = False
    plots_module: str = "plots_1"
    bayesian_estimation_module: str = "sympleq.core.bayesian_estimation"
    estimator_threshold: float = 0.0
    estimator_min_runs: int = 1

    # AEPsych / GP / acquisition
    optimization_steps: int = 10000
    inducing_size: int = 150
    acquisition_function: str = "GlobalSUR"
    acquisition_restarts: int = 8
    acquisition_samples: int = 30000

    # Batching
    batching: bool = True
    max_batch_size: int = 80
    fantasy_batching: bool = True

    # GPU
    use_gpu: bool = True
    force_default_device_during_aepsych: bool = True

    # Debugging
    verbose_fantasies: bool = True
    print_diagnostics: bool = False


# =============================================================================
# Device handling
# =============================================================================


def choose_gp_device(settings: RickFantasyGPUSettings) -> torch.device:
    """
    Choose device for AEPsych/GP/acquisition.

    RMB measurement remains CPU regardless of this.
    """

    if settings.use_gpu and torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


@contextmanager
def temporary_torch_default_device(device: torch.device, enabled: bool):
    """
    Temporarily make PyTorch create tensors on the chosen device.

    This helps AEPsych/BoTorch internals create acquisition tensors on CUDA.

    Turn off with:
        force_default_device_during_aepsych=False

    if your AEPsych version has device mismatch problems.
    """

    if not enabled or device.type == "cpu" or not hasattr(torch, "set_default_device"):
        yield
        return

    get_default_device = getattr(torch, "get_default_device", None)
    old_device = get_default_device() if get_default_device is not None else "cpu"

    torch.set_default_device(device)
    try:
        yield
    finally:
        torch.set_default_device(old_device)


def move_strategy_models_to_device(
    strategy: SequentialStrategy,
    device: torch.device,
) -> None:
    """
    Best-effort move of AEPsych model objects to CUDA.

    This is defensive because AEPsych internals differ slightly across versions.
    """

    if device.type == "cpu":
        return

    possible_objects = [strategy]

    for attr in ("_strat", "strat", "current_strategy"):
        obj = getattr(strategy, attr, None)
        if obj is not None:
            possible_objects.append(obj)

    for attr in ("strat_list", "_strat_list", "strategies"):
        objs = getattr(strategy, attr, None)
        if objs is not None:
            possible_objects.extend(list(objs))

    for obj in possible_objects:
        model = getattr(obj, "model", None)
        if model is not None and hasattr(model, "to"):
            try:
                model.to(device)
            except Exception:
                pass


# =============================================================================
# Basic coordinate helpers
# =============================================================================


def contour_target(settings: RickFantasyGPUSettings) -> float:
    """Return and validate the target contour value."""

    target = float(settings.target_threshold)
    if not 0.0 < target < 1.0:
        raise ValueError(f"target_threshold must be in (0, 1), got {target}.")
    return target


def point_from_config(config: RMBConfig, *, device: torch.device) -> torch.Tensor:
    """
    Convert RMBConfig to AEPsych coordinate.

    Returns shape:
        (1, 2)

    Coordinates:
        x[0, 0] = total gates
        x[0, 1] = two-qubit gate ratio
    """

    return torch.tensor(
        config_points([config]),
        dtype=torch.double,
        device=device,
    )


def raw_point(
    n_gates: float,
    ratio: float,
    *,
    device: torch.device,
) -> torch.Tensor:
    """Create an AEPsych point directly from raw coordinates."""

    return torch.tensor(
        [[float(n_gates), float(ratio)]],
        dtype=torch.double,
        device=device,
    )


def config_from_aepsych_x(
    x: torch.Tensor,
    settings: RickFantasyGPUSettings,
) -> RMBConfig:
    """Convert AEPsych's generated point into an RMBConfig."""

    x_cpu = x.detach().cpu()
    return settings.make_config(
        float(x_cpu[0, 0].item()),
        float(x_cpu[0, 1].item()),
    )


def one_shot_requests(configs: list[RMBConfig]) -> list[MeasurementRequest]:
    """
    Convert configs into one-shot measurement requests.

    Each selected config gives one binary success/failure outcome.
    """

    return [MeasurementRequest(config, 1) for config in configs]


def validation_requests(
    configs: list[RMBConfig],
    settings: RickFantasyGPUSettings,
) -> list[MeasurementRequest]:
    """Validation requests can use more shots per config than training."""
    shots = max(1, int(settings.validation_shots_per_config))
    return [MeasurementRequest(config, shots) for config in configs]


# =============================================================================
# AEPsych strategy construction
# =============================================================================


def aepsych_config_string(settings: RickFantasyGPUSettings) -> str:
    """
    Build AEPsych config for the GP/GlobalSUR phase only.

    Note:
        Sobol is handled explicitly by this file, not by AEPsych's init_strat.
    """

    target = float(contour_target(settings))
    n_gates_lower = float(settings.n_gates_bounds[0])
    n_gates_upper = float(settings.n_gates_bounds[1])
    ratio_lower = float(settings.ratio_bounds[0])
    ratio_upper = float(settings.ratio_bounds[1])
    optimization_steps = int(settings.optimization_steps)
    inducing_size = int(settings.inducing_size)
    acquisition_restarts = int(settings.acquisition_restarts)
    acquisition_samples = int(settings.acquisition_samples)

    return f"""
    [common]
    parnames = [n_gates, ratio]
    outcome_types = [binary]
    strategy_names = [opt_strat]

    [n_gates]
    par_type = continuous
    lower_bound = {n_gates_lower}
    upper_bound = {n_gates_upper}

    [ratio]
    par_type = continuous
    lower_bound = {ratio_lower}
    upper_bound = {ratio_upper}

    [opt_strat]
    min_asks = {optimization_steps}
    model = GPClassificationModel
    generator = OptimizeAcqfGenerator

    [GPClassificationModel]
    inducing_size = {inducing_size}
    likelihood = BernoulliLikelihood

    [OptimizeAcqfGenerator]
    acqf = {settings.acquisition_function}
    restarts = {acquisition_restarts}
    samps = {acquisition_samples}

    [{settings.acquisition_function}]
    target = {target}
    """.strip()


def sanitize_aepsych_config_numbers(config: Config) -> None:
    """Rewrite NumPy scalar reprs in AEPsych config values to plain numbers."""

    for section in config.sections():
        for option, value in config.items(section):
            cleaned = clean_aepsych_numeric_string(value)
            if cleaned != value:
                config.set(section, option, cleaned)


def clean_aepsych_numeric_string(value: str) -> str:
    """Convert NumPy scalar reprs in config strings to plain numeric literals."""

    numpy_scalar_pattern = re.compile(r"(?:np\.)?(?:float|int)\d*\(([^()]*)\)")
    return numpy_scalar_pattern.sub(r"\1", value)


def sanitized_aepsych_transform_options(config: Config, transforms=None) -> Config:
    """AEPsych transform_options wrapper that cleans NumPy scalar reprs."""

    transformed_config = _ORIGINAL_AEPSYCH_TRANSFORM_OPTIONS(config, transforms)
    sanitize_aepsych_config_numbers(transformed_config)
    return transformed_config


def sanitized_aepsych_str_to_list(self: Config, value: str, element_type=float):
    """AEPsych parser wrapper that accepts NumPy scalar reprs."""

    cleaned = clean_aepsych_numeric_string(value)
    return _ORIGINAL_AEPSYCH_STR_TO_LIST(self, cleaned, element_type)


def sanitized_aepsych_str_to_array(self: Config, value: str) -> np.ndarray:
    """AEPsych parser wrapper that accepts NumPy scalar reprs."""

    cleaned = clean_aepsych_numeric_string(value)
    return _ORIGINAL_AEPSYCH_STR_TO_ARRAY(self, cleaned)


def build_strategy(settings: RickFantasyGPUSettings) -> SequentialStrategy:
    """Build the AEPsych GP/acquisition strategy."""

    config = Config()
    config.update(config_str=aepsych_config_string(settings))
    sanitize_aepsych_config_numbers(config)
    aepsych_parameter_transforms.transform_options = sanitized_aepsych_transform_options
    Config._str_to_list = sanitized_aepsych_str_to_list
    Config._str_to_array = sanitized_aepsych_str_to_array
    return SequentialStrategy.from_config(config)


def add_observation_to_strategy(
    strategy: SequentialStrategy,
    observation: Observation,
    *,
    device: torch.device,
) -> None:
    """Add one observation to an AEPsych strategy."""

    x = observation.x_cpu.to(device=device, dtype=torch.double)
    strategy.add_data(x, [int(observation.y)])


def replay_observations(
    strategy: SequentialStrategy,
    observations: list[Observation],
    *,
    device: torch.device,
) -> None:
    """Replay stored fake/real observations into a rebuilt strategy."""

    for obs in observations:
        add_observation_to_strategy(strategy, obs, device=device)


# =============================================================================
# Fake anchors
# =============================================================================


def seed_fake_corners(
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    observations: list[Observation],
    results_for_plot: list[tuple[float, float, int]],
    *,
    device: torch.device,
) -> None:
    """
    Add fake easy/hard corner anchors.

    These are not measured.
    These are not charged.
    """

    if not settings.use_fake_corners:
        return

    corners = [
        (
            float(settings.n_gates_bounds[0]),
            float(settings.ratio_bounds[0]),
            int(settings.easy_corner_outcome),
        ),
        (
            float(settings.n_gates_bounds[1]),
            float(settings.ratio_bounds[1]),
            int(settings.hard_corner_outcome),
        ),
    ]
    corners.extend(settings.extra_fake_anchors)

    for n_gates, ratio, outcome in corners:
        n_gates = float(np.clip(n_gates, *settings.n_gates_bounds))
        ratio = float(np.clip(ratio, *settings.ratio_bounds))
        outcome = int(outcome)
        obs = Observation(
            x_cpu=raw_point(n_gates, ratio, device=torch.device("cpu")),
            y=outcome,
            source="fake",
        )

        add_observation_to_strategy(strategy, obs, device=device)
        observations.append(obs)

        # Plotting uses native AEPsych coordinates: total gates and ratio.
        config = settings.make_config(n_gates, ratio)
        results_for_plot.append(
            (
                float(config.n_gates),
                float(config.ratio_2_qb_gates),
                outcome,
            )
        )


# =============================================================================
# Explicit Sobol warm-up
# =============================================================================


def sobol_initial_candidates(
    settings: RickFantasyGPUSettings,
) -> list[RMBConfig]:
    """
    Generate explicit Sobol warm-up candidates.

    Sobol gives points in [0, 1]^2, then we rescale to:
        [n_gates_min, n_gates_max] x [ratio_min, ratio_max]
    """

    if settings.initial_sobol_samples <= 0:
        return []

    engine = torch.quasirandom.SobolEngine(
        dimension=2,
        scramble=bool(settings.sobol_scramble),
        seed=settings.rng_seed if settings.sobol_scramble else None,
    )

    unit_points = engine.draw(settings.initial_sobol_samples).cpu().numpy()

    n_min, n_max = settings.n_gates_bounds
    r_min, r_max = settings.ratio_bounds

    candidates: list[RMBConfig] = []
    seen: set[RMBConfig] = set()

    for u_n, u_r in unit_points:
        n_gates = float(n_min + u_n * (n_max - n_min))
        ratio = float(r_min + u_r * (r_max - r_min))

        config = settings.make_config(n_gates, ratio)

        # Realized gate-count rounding can create duplicates.
        if config not in seen:
            candidates.append(config)
            seen.add(config)

    return candidates


def select_affordable_prefix(
    candidates: list[RMBConfig],
    settings: RickFantasyGPUSettings,
    budget: Budget,
    *,
    max_cost_per_run: float | None = None,
) -> tuple[list[RMBConfig], bool]:
    """
    Select as many candidates as fit in one stitched batch.

    Returns:
        selected, exhausted
    """

    selected: list[RMBConfig] = []
    cost_cap = settings.max_cost_per_run if max_cost_per_run is None else max_cost_per_run

    for candidate in candidates:
        next_cost = batch_hqc_cost(one_shot_requests(selected + [candidate]))

        if next_cost > cost_cap:
            break

        if next_cost > budget.remaining_hqc:
            return selected, True

        selected.append(candidate)

    return selected, False


# =============================================================================
# Real measurement
# =============================================================================


def measure_batch_and_update_real_strategy(
    *,
    phase: str,
    selected: list[RMBConfig],
    strategy: SequentialStrategy,
    rmb: RMB,
    rng: np.random.Generator,
    data,
    budget: Budget,
    settings: RickFantasyGPUSettings,
    observations: list[Observation],
    results_for_plot: list[tuple[float, float, int]],
    device: torch.device,
) -> None:
    """
    Measure selected configs using the CPU RMB backend.

    Only real outcomes from here are added to the real strategy.
    """

    if not selected:
        return

    # CPU-side RMB measurement.
    outcomes_by_config = spend_request_batch(
        rmb.backend,
        rng,
        data,
        one_shot_requests(selected),
        seed=settings.rng_seed,
    )

    # Add real measured data to the real GP strategy.
    for config, outcomes in outcomes_by_config.items():
        for outcome in outcomes:
            obs = Observation(
                x_cpu=point_from_config(config, device=torch.device("cpu")),
                y=int(outcome),
                source=phase,
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

    cost = batch_hqc_cost(one_shot_requests(selected))
    budget.spend_batch(cost, len(selected))

    print_progress(
        settings,
        budget,
        f"{phase}: measured {len(selected)} configs",
    )


# =============================================================================
# Fantasy GlobalSUR batching
# =============================================================================


def copy_strategy_for_fantasies(
    real_strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    observations: list[Observation],
    *,
    device: torch.device,
) -> SequentialStrategy:
    """
    Copy the real strategy for fantasy batching.

    Preferred path:
        deepcopy(real_strategy)

    Fallback path:
        build a fresh strategy and replay fake + real observations.
    """

    try:
        fantasy_strategy = copy.deepcopy(real_strategy)
        move_strategy_models_to_device(fantasy_strategy, device)
        return fantasy_strategy
    except Exception:
        fantasy_strategy = build_strategy(settings)
        replay_observations(fantasy_strategy, observations, device=device)
        move_strategy_models_to_device(fantasy_strategy, device)
        return fantasy_strategy


def refreshed_strategy_for_prediction(
    real_strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    observations: list[Observation],
    *,
    device: torch.device,
) -> SequentialStrategy:
    """
    Copy the real strategy and force AEPsych to refresh its fitted model.

    AEPsych updates stored data with add_data(...), but the fitted model can
    remain stale until gen() runs. Calling gen() on a copy gives plotting and
    diagnostics the same freshly-fit model path used by GlobalSUR, without
    mutating the real training strategy.
    """
    prediction_strategy = copy_strategy_for_fantasies(
        real_strategy,
        settings,
        observations,
        device=device,
    )
    move_strategy_models_to_device(prediction_strategy, device)
    with temporary_torch_default_device(
        device,
        enabled=settings.force_default_device_during_aepsych,
    ):
        try:
            prediction_strategy.gen()
        except Exception as exc:
            print(f"[plot model refresh] warning: gen() refresh failed: {exc}")
    return prediction_strategy


def predict_success_probability(
    strategy: SequentialStrategy,
    config: RMBConfig,
    settings: RickFantasyGPUSettings,
    *,
    device: torch.device,
) -> float:
    """
    Predict P(success | config) from the GP classifier.
    """

    if strategy.model is None:
        return 0.5

    x = point_from_config(config, device=device)
    move_strategy_models_to_device(strategy, device)

    with temporary_torch_default_device(
        device,
        enabled=settings.force_default_device_during_aepsych,
    ):
        with torch.no_grad():
            try:
                p, _ = strategy.model.predict(x, probability_space=True)
            except RuntimeError:
                # Fallback for device mismatch.
                p, _ = strategy.model.predict(
                    x.detach().cpu(),
                    probability_space=True,
                )

    p_float = float(p.detach().cpu().reshape(-1)[0])
    return min(max(p_float, 0.0), 1.0)


def print_real_strategy_prediction_diagnostic(
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    results_for_plot: list[tuple[float, float, int]],
    *,
    device: torch.device,
    limit: int = 12,
) -> None:
    """Print real-strategy predictions at actual plotted measurement points."""
    if strategy.model is None:
        print("[real strategy prediction diagnostic] model is None")
        return

    unique_points: list[tuple[float, float, int]] = []
    seen: set[tuple[int, float]] = set()
    for n_gates, ratio, outcome in reversed(results_for_plot):
        key = (int(round(n_gates)), round(float(ratio), 6))
        if key in seen:
            continue
        seen.add(key)
        unique_points.append((float(n_gates), float(ratio), int(outcome)))
        if len(unique_points) >= limit:
            break
    unique_points.reverse()

    predictions: list[float] = []
    print("[real strategy prediction diagnostic]")
    print(f"  checked_points                       = {len(unique_points)}")
    for index, (n_gates, ratio, outcome) in enumerate(unique_points, start=1):
        config = settings.make_config(n_gates, ratio)
        p_success = predict_success_probability(
            strategy,
            config,
            settings,
            device=device,
        )
        predictions.append(p_success)
        print(
            f"  {index:03d}: "
            f"n_gates={config.n_gates} "
            f"ratio={config.ratio_2_qb_gates:.4f} "
            f"outcome={outcome} "
            f"real_model_p_success={p_success:.6f}"
        )
    if predictions:
        print(
            "  real_model_p_success_range           = "
            f"{min(predictions):.6f} to {max(predictions):.6f}"
        )


def add_fantasy_outcome(
    fantasy_strategy: SequentialStrategy,
    config: RMBConfig,
    virtual_outcome: int,
    *,
    device: torch.device,
) -> None:
    """
    Add virtual success/failure to copied strategy only.
    """

    obs = Observation(
        x_cpu=point_from_config(config, device=torch.device("cpu")),
        y=int(virtual_outcome),
        source="fantasy",
    )
    add_observation_to_strategy(fantasy_strategy, obs, device=device)


def select_fantasy_globalsur_batch(
    real_strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    budget: Budget,
    observations: list[Observation],
    fantasy_rng: np.random.Generator,
    *,
    device: torch.device,
) -> tuple[list[RMBConfig], bool]:
    """
    Select one batch using copied-strategy fantasy GlobalSUR.

    This is the key modified algorithm.
    """

    selected: list[RMBConfig] = []
    seen: set[RMBConfig] = set()

    fantasy_strategy = copy_strategy_for_fantasies(
        real_strategy,
        settings,
        observations,
        device=device,
    )

    batch_limit = settings.max_batch_size if settings.batching else 1
    max_attempts = max(5 * batch_limit, batch_limit + 10)

    attempts = 0

    while len(selected) < batch_limit and attempts < max_attempts:
        attempts += 1

        move_strategy_models_to_device(fantasy_strategy, device)

        with temporary_torch_default_device(
            device,
            enabled=settings.force_default_device_during_aepsych,
        ):
            x = fantasy_strategy.gen()

        candidate = config_from_aepsych_x(x, settings)

        # Rounding can cause duplicate realized configs.
        if candidate in seen:
            continue

        next_cost = batch_hqc_cost(one_shot_requests(selected + [candidate]))

        if next_cost > settings.max_cost_per_run:
            break

        if next_cost > budget.remaining_hqc:
            return selected, True

        p_success = predict_success_probability(
            fantasy_strategy,
            candidate,
            settings,
            device=device,
        )

        virtual_outcome = int(fantasy_rng.random() < p_success)

        add_fantasy_outcome(
            fantasy_strategy,
            candidate,
            virtual_outcome,
            device=device,
        )

        selected.append(candidate)
        seen.add(candidate)

        if settings.verbose_fantasies:
            print(
                "[fantasy] "
                f"batch_index={len(selected):03d} "
                f"n_gates={candidate.n_gates} "
                f"ratio={candidate.ratio_2_qb_gates:.4f} "
                f"p_success={p_success:.4f} "
                f"virtual={virtual_outcome} "
                f"batch_cost={next_cost:.3f}"
            )

    return selected, False


def select_plain_aepsych_batch(
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    budget: Budget,
    *,
    device: torch.device,
) -> tuple[list[RMBConfig], bool]:
    """
    Non-fantasy fallback batch selection.

    Useful for comparison/debugging.
    """

    selected: list[RMBConfig] = []
    seen: set[RMBConfig] = set()

    batch_limit = settings.max_batch_size if settings.batching else 1
    max_attempts = max(5 * batch_limit, batch_limit + 10)

    attempts = 0

    while len(selected) < batch_limit and attempts < max_attempts:
        attempts += 1

        move_strategy_models_to_device(strategy, device)

        with temporary_torch_default_device(
            device,
            enabled=settings.force_default_device_during_aepsych,
        ):
            x = strategy.gen()

        candidate = config_from_aepsych_x(x, settings)

        if candidate in seen:
            continue

        next_cost = batch_hqc_cost(one_shot_requests(selected + [candidate]))

        if next_cost > settings.max_cost_per_run:
            break

        if next_cost > budget.remaining_hqc:
            return selected, True

        selected.append(candidate)
        seen.add(candidate)

    return selected, False


# =============================================================================
# Final GP contour extraction
# =============================================================================


def predict_level_set(
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    *,
    one_q_bounds: tuple[int, int],
    two_q_bounds: tuple[int, int],
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluate fitted GP on a grid for plotting.
    """

    plots = import_experiment_module(settings.plots_module)
    gate_plane_points = plots.gate_plane_points
    level_set_grid = plots.level_set_grid

    one_q_grid, two_q_grid = level_set_grid(one_q_bounds, two_q_bounds)
    points = gate_plane_points(one_q_grid, two_q_grid)

    grid = torch.tensor(
        np.column_stack(
            [
                np.clip(
                    points[:, 0],
                    settings.n_gates_bounds[0],
                    settings.n_gates_bounds[1],
                ),
                np.clip(
                    points[:, 1],
                    settings.ratio_bounds[0],
                    settings.ratio_bounds[1],
                ),
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

    return (
        probabilities.detach().cpu().numpy().reshape(one_q_grid.shape),
        latent_mean.detach().cpu().numpy().reshape(one_q_grid.shape),
        latent_variance.detach().cpu().numpy().reshape(one_q_grid.shape),
    )


def predict_native_level_set(
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
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


def level_set_configs(
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    *,
    device: torch.device,
    measured_data=None,
    debug: bool = False,
    debug_limit: int = 12,
) -> list[RMBConfig]:
    """
    Extract configs on the GP-predicted P(success)=target contour.
    """

    if strategy.model is None:
        return []

    gates_axis, ratio_axis = candidate_axes(settings)
    gates_grid, ratio_grid = np.meshgrid(gates_axis, ratio_axis)

    grid = torch.tensor(
        np.stack([gates_grid.ravel(), ratio_grid.ravel()], axis=1),
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

    probability_grid = probabilities.detach().cpu().numpy().reshape(gates_grid.shape)
    target = contour_target(settings)
    delta = probability_grid - target

    if debug:
        flat_order = np.argsort(np.abs(delta.ravel()))
        measured_points = None
        if measured_data:
            measured_points = config_points(list(measured_data.keys()))
            lower = np.array([settings.n_gates_bounds[0], settings.ratio_bounds[0]], dtype=float)
            upper = np.array([settings.n_gates_bounds[1], settings.ratio_bounds[1]], dtype=float)
            span = np.maximum(upper - lower, 1e-12)
            measured_points = (measured_points - lower) / span

        def nearest_measured_distance(n_gates: float, ratio: float) -> float | None:
            if measured_points is None or len(measured_points) == 0:
                return None
            point = (
                np.array([float(n_gates), float(ratio)], dtype=float)
                - lower
            ) / span
            distances = np.linalg.norm(measured_points - point, axis=1)
            return float(np.min(distances))

        print("[GP contour grid debug]")
        print(f"  target                               = {target}")
        print(f"  probability range on grid            = "
              f"{float(np.min(probability_grid)):.6g} to {float(np.max(probability_grid)):.6g}")
        band_lo, band_hi = settings.validation_probability_band
        in_band = (band_lo <= probability_grid) & (probability_grid <= band_hi)
        print(f"  validation probability band          = [{band_lo}, {band_hi}]")
        print(f"  grid points in validation band       = {int(np.count_nonzero(in_band))}")
        print(f"  grid points closest to target        = {min(debug_limit, len(flat_order))}")
        for flat_index in flat_order[:debug_limit]:
            row, col = np.unravel_index(int(flat_index), delta.shape)
            config = settings.make_config(float(gates_grid[row, col]), float(ratio_grid[row, col]))
            nearest_distance = nearest_measured_distance(
                float(gates_grid[row, col]),
                float(ratio_grid[row, col]),
            )
            nearest_text = "n/a" if nearest_distance is None else f"{nearest_distance:.4f}"
            print(
                "    "
                f"p={probability_grid[row, col]:.6g}, "
                f"delta={delta[row, col]:+.6g}, "
                f"grid_gates={float(gates_grid[row, col]):.6g}, "
                f"grid_ratio={float(ratio_grid[row, col]):.6g}, "
                f"nearest_measured_norm={nearest_text}, "
                f"config=(gates={config.n_gates}, ratio={config.ratio_2_qb_gates:.4f})"
            )

    configs: list[RMBConfig] = []
    seen: set[RMBConfig] = set()

    def add_config(n_gates: float, ratio: float) -> None:
        config = settings.make_config(n_gates, ratio)
        if config not in seen:
            seen.add(config)
            configs.append(config)

    band_lo, band_hi = settings.validation_probability_band
    band_rows, band_cols = np.where((band_lo <= probability_grid) & (probability_grid <= band_hi))
    for row, col in zip(band_rows, band_cols):
        add_config(float(gates_grid[row, col]), float(ratio_grid[row, col]))

    for row in range(delta.shape[0]):
        crossing = np.where(delta[row, :-1] * delta[row, 1:] < 0)[0]

        if len(crossing) == 0:
            continue

        i = int(crossing[0])

        left = abs(delta[row, i])
        right = abs(delta[row, i + 1])
        t = left / (left + right)

        n_gates = float((1.0 - t) * gates_axis[i] + t * gates_axis[i + 1])
        ratio = float(ratio_axis[row])
        add_config(n_gates, ratio)

    return sorted(
        configs,
        key=lambda c: (c.ratio_2_qb_gates, c.n_gates),
    )


def validation_seed(settings: RickFantasyGPUSettings) -> int | None:
    if settings.rng_seed is None:
        return None
    return int(settings.rng_seed) + int(settings.validation_seed_offset)


def validation_output_path(
    settings: RickFantasyGPUSettings,
    base_path: Path | None,
) -> Path | None:
    if settings.validation_save_path is not None:
        return Path(settings.validation_save_path)
    if base_path is None:
        return None
    return base_path.parent / f"{base_path.stem}_validation.json"


def write_scores_to_json(base_path: Path, scores: dict[str, float]) -> None:
    """Attach modified-crossing scores to the main RMB JSON file."""
    payload = json.loads(base_path.read_text(encoding="utf-8"))
    payload["scores"] = scores
    base_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def effective_training_hqc_budget(settings: RickFantasyGPUSettings) -> float:
    """Training/acquisition budget after optional validation reservation."""
    if not settings.validate_gp_contour or not settings.reserve_validation_budget:
        return settings.hqc_budget
    return max(0.0, settings.hqc_budget - settings.validation_hqc_budget)


def select_aepsych_validation_configs(
    *,
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    budget: Budget,
    observations: list[Observation],
    fantasy_rng: np.random.Generator,
    device: torch.device,
    cost_cap: float | None = None,
) -> tuple[list[RMBConfig], bool]:
    """Ask the trained AEPsych strategy for validation configs."""
    if strategy.model is None:
        return [], False

    selected: list[RMBConfig] = []
    seen: set[RMBConfig] = set()
    fantasy_strategy = copy_strategy_for_fantasies(
        strategy,
        settings,
        observations,
        device=device,
    )

    target_count = max(0, int(settings.validation_target_configs))
    max_attempts = max(
        target_count + 10,
        int(settings.validation_max_attempts_multiplier) * max(1, target_count),
    )
    if cost_cap is None:
        cost_cap = (
            settings.max_cost_per_run
            if settings.validation_max_cost_per_run is None
            else settings.validation_max_cost_per_run
        )

    attempts = 0
    exhausted = False
    while len(selected) < target_count and attempts < max_attempts:
        attempts += 1
        move_strategy_models_to_device(fantasy_strategy, device)

        with temporary_torch_default_device(
            device,
            enabled=settings.force_default_device_during_aepsych,
        ):
            x = fantasy_strategy.gen()

        candidate = config_from_aepsych_x(x, settings)
        if candidate in seen:
            continue

        next_cost = batch_hqc_cost(validation_requests(selected + [candidate], settings))
        if next_cost > cost_cap:
            break
        if next_cost > budget.remaining_hqc:
            exhausted = True
            break

        p_success = predict_success_probability(
            fantasy_strategy,
            candidate,
            settings,
            device=device,
        )
        virtual_outcome = int(fantasy_rng.random() < p_success)
        add_fantasy_outcome(
            fantasy_strategy,
            candidate,
            virtual_outcome,
            device=device,
        )

        selected.append(candidate)
        seen.add(candidate)
        print(
            "[validation candidate] "
            f"index={len(selected):03d} "
            f"n_gates={candidate.n_gates} "
            f"ratio={candidate.ratio_2_qb_gates:.4f} "
            f"p_success={p_success:.4f} "
            f"virtual={virtual_outcome} "
            f"batch_cost={next_cost:.3f}"
        )

    return selected, exhausted


def validate_gp_contour_configs(
    *,
    strategy: SequentialStrategy,
    observations: list[Observation],
    rmb: RMB,
    rng: np.random.Generator,
    fantasy_rng: np.random.Generator,
    settings: RickFantasyGPUSettings,
    base_path: Path | None,
    device: torch.device,
    validation_hqc_budget: float | None = None,
) -> tuple[RMB | None, Path | None, list[RMBConfig], list[tuple[float, float, int]], Budget | None]:
    """Measure trained-AEPsych validation suggestions into a separate RMBData."""
    if not settings.validate_gp_contour:
        return None, None, [], [], None
    validation_hqc_budget = (
        settings.validation_hqc_budget
        if validation_hqc_budget is None
        else validation_hqc_budget
    )
    if validation_hqc_budget <= 0.0:
        print("[validation] skipped: validation_hqc_budget <= 0")
        return None, None, [], [], None

    validation_budget = Budget(remaining_hqc=validation_hqc_budget)
    selected: list[RMBConfig] = []
    exhausted = False
    restart_index = 0
    while len(selected) < settings.validation_min_configs:
        selected, exhausted = select_aepsych_validation_configs(
            strategy=strategy,
            settings=settings,
            budget=validation_budget,
            observations=observations,
            fantasy_rng=fantasy_rng,
            device=device,
            cost_cap=validation_budget.remaining_hqc,
        )
        if len(selected) >= settings.validation_min_configs:
            break
        restart_index += 1
        print(
            "[validation restart] "
            f"selected={len(selected)} "
            f"minimum={settings.validation_min_configs} "
            f"restart={restart_index}/{settings.validation_max_restarts}"
        )
        if restart_index >= settings.validation_max_restarts:
            raise RuntimeError(
                "Validation candidate generation could not reach "
                f"{settings.validation_min_configs} configs after "
                f"{settings.validation_max_restarts} restarts."
            )

    print("[validation configs]")
    for index, config in enumerate(selected, start=1):
        print(
            f"  {index:03d}: "
            f"n_gates={config.n_gates}, "
            f"ratio={config.ratio_2_qb_gates:.4f}, "
            f"n_1qb={config.n_1qb_gates}, "
            f"n_2qb={config.n_2qb_gates}"
        )

    validation_rmb = RMB.default(rng).with_backend(rmb.backend)
    requests = validation_requests(selected, settings)
    outcomes_by_config = spend_request_batch(
        validation_rmb.backend,
        rng,
        validation_rmb._data,
        requests,
        seed=validation_seed(settings),
    )
    validation_results_for_plot: list[tuple[float, float, int]] = []
    print("[validation measurement outcomes]")
    for index, (config, outcomes) in enumerate(outcomes_by_config.items(), start=1):
        estimator = validation_rmb._data[config]
        validation_results_for_plot.extend(
            (
                float(config.n_1qb_gates),
                float(config.n_2qb_gates),
                int(outcome),
            )
            for outcome in outcomes
        )
        print(
            f"  {index:03d}: "
            f"n_gates={config.n_gates}, "
            f"ratio={config.ratio_2_qb_gates:.4f}, "
            f"outcomes={list(map(int, outcomes))}, "
            f"probability_true={estimator.probability(True):.4f}, "
            f"posterior_mean={estimator.posterior_mean():.4f}"
        )
    cost = batch_hqc_cost(requests)
    validation_budget.spend_batch(cost, len(selected))

    output_path = validation_output_path(settings, base_path)
    resolved_path = None
    if output_path is not None:
        validation_rmb.save(output_path)
        resolved_path = output_path if output_path.is_absolute() else Path(output_path)

    print("[validation]")
    print(f"  validation_available_hqc             = {validation_hqc_budget:.6g}")
    print(f"  requested validation configs         = {settings.validation_target_configs}")
    print(f"  selected validation configs          = {len(selected)}")
    print(f"  measured validation configs          = {len(outcomes_by_config)}")
    print(f"  validation_spent_hqc                 = {validation_budget.spent_hqc:.6g}")
    print(f"  validation_remaining_hqc             = {validation_budget.remaining_hqc:.6g}")
    print(f"  validation_exhausted                 = {exhausted}")
    if resolved_path is not None:
        print(f"  validation_save_path                 = {resolved_path}")

    return validation_rmb, resolved_path, selected, validation_results_for_plot, validation_budget


# =============================================================================
# Logging handles
# =============================================================================


def print_run_handles(
    settings: RickFantasyGPUSettings,
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

    if settings.validate_gp_contour:
        print("[validation]")
        print(f"  reserve_validation_budget            = {settings.reserve_validation_budget}")
        print(f"  validation_hqc_budget                = {settings.validation_hqc_budget}")
        print(f"  validation_max_cost_per_run          = {settings.validation_max_cost_per_run}")
        print(f"  validation_target_configs            = {settings.validation_target_configs}")
        print(f"  validation_min_configs               = {settings.validation_min_configs}")
        print(f"  validation_max_restarts              = {settings.validation_max_restarts}")
        print(f"  validation_shots_per_config          = {settings.validation_shots_per_config}")
        print(f"  validation_max_attempts_multiplier   = "
              f"{settings.validation_max_attempts_multiplier}")
        print(f"  validation_probability_band          = {settings.validation_probability_band}")
        print(f"  validation_seed_offset               = {settings.validation_seed_offset}")
        print(f"  validation_save_path                 = {settings.validation_save_path}")
    print(f"  plot_gp_level_set_result             = {settings.plot_gp_level_set_result}")
    print(f"  plots_module                         = {settings.plots_module}")
    print(f"  bayesian_estimation_module           = {settings.bayesian_estimation_module}")
    print(f"  estimator_threshold                  = {settings.estimator_threshold}")
    print(f"  estimator_min_runs                   = {settings.estimator_min_runs}")

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
    print(f"  effective_training_hqc_budget        = "
          f"{effective_training_hqc_budget(settings)}")

    print("[search box]")
    print(f"  n_gates_bounds                       = {settings.n_gates_bounds}")
    print(f"  ratio_bounds                         = {settings.ratio_bounds}")

    print("[device]")
    print(f"  use_gpu                              = {settings.use_gpu}")
    print(f"  selected GP device                   = {gp_device}")
    print(f"  force_default_device_during_aepsych  = "
          f"{settings.force_default_device_during_aepsych}")
    print("  RMB backend                          = SympleQ")

    print("[debug]")
    print(f"  rng_seed                             = {settings.rng_seed}")
    print(f"  verbose_fantasies                    = {settings.verbose_fantasies}")
    print(f"  print_diagnostics                    = {settings.print_diagnostics}")

    print("=============================================\n")


# =============================================================================
# Main run logic
# =============================================================================


def run(
    settings: RickFantasyGPUSettings,
    *,
    return_budget: bool = False,
):
    """
    Run the experiment.
    """

    logging.getLogger().setLevel(logging.WARNING)
    warnings.filterwarnings("ignore")

    apply_bayesian_estimation_module(
        settings.bayesian_estimation_module,
        default_threshold=settings.estimator_threshold,
        default_min_runs=settings.estimator_min_runs,
    )

    torch.set_default_dtype(torch.float64)

    gp_device = choose_gp_device(settings)

    if settings.rng_seed is not None:
        torch.manual_seed(settings.rng_seed)
        if gp_device.type == "cuda":
            torch.cuda.manual_seed_all(settings.rng_seed)

    print_run_handles(settings, gp_device=gp_device)

    rng, rmb, budget = start_run(settings)
    budget.remaining_hqc = effective_training_hqc_budget(settings)
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

    crossings = level_set_configs(
        strategy,
        settings,
        device=gp_device,
        measured_data=data,
        debug=False,
    )

    base_path = None
    if settings.save_path is not None:
        rmb.save(settings.save_path)
        base_path = resolve_data_path(settings.save_path)

    gp_posterior_mean_for_validation_plot = None
    if settings.plot and settings.plot_gp_level_set_result and base_path is not None:
        plots = import_experiment_module(settings.plots_module)
        plot_gp_level_set = plots.plot_gp_level_set

        plot_strategy = refreshed_strategy_for_prediction(
            strategy,
            settings,
            observations,
            device=gp_device,
        )
        if settings.print_diagnostics:
            print_real_strategy_prediction_diagnostic(
                plot_strategy,
                settings,
                results_for_plot,
                device=gp_device,
            )

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
            score_summary = gp_grid_scores(
                grid_path,
                one_q_noise_scale=float(getattr(settings, "one_q_noise_scale", 1.0)),
                two_q_noise_scale=float(getattr(settings, "two_q_noise_scale", 1.0)),
            )
            write_scores_to_json(base_path, score_summary)
            print("[scores]")
            print(f"  S1                                  = {score_summary['S1']:.6g}")
            print(f"  S2                                  = {score_summary['S2']:.6g}")
            print(f"  A_gp                                = {score_summary['A_gp']:.6g}")
        gp_posterior_mean_for_validation_plot = (gates_grid, ratio_grid, probabilities)
        plot_gp_level_set(
            probabilities,
            latent_mean,
            latent_variance,
            results_for_plot,
            one_q_bounds=(0, 1),
            two_q_bounds=(0, 1),
            settings=settings,
            target=contour_target(settings),
            title=f"AEPsych level-set result | seed={settings.rng_seed}",
            log_x=True,
            x_grid=gates_grid,
            y_grid=ratio_grid,
            x_label="# Gates",
            y_label="Two-qubit gate ratio",
            coordinate_system="total_ratio",
            png_path=base_path.parent / f"{base_path.stem}_level_set.png",
            show=False,
        )

    (
        validation_rmb,
        validation_base_path,
        validation_crossings,
        validation_results_for_plot,
        validation_budget,
    ) = validate_gp_contour_configs(
        strategy=strategy,
        observations=observations,
        rmb=rmb,
        rng=rng,
        fantasy_rng=fantasy_rng,
        settings=settings,
        base_path=base_path,
        device=gp_device,
        validation_hqc_budget=(
            settings.validation_hqc_budget
            + max(0.0, budget.remaining_hqc)
        ),
    )
    # -------------------------------------------------------------------------
    # 5. Optional plotting
    # -------------------------------------------------------------------------

    if settings.plot:
        plots = import_experiment_module(settings.plots_module)
        plot_monotone_fidelity_surface_with_confidence = (
            plots.plot_monotone_fidelity_surface_with_confidence
        )

        if validation_rmb is not None and validation_base_path is not None:
            axes = plot_monotone_fidelity_surface_with_confidence(
                validation_rmb._data,
                settings,
                png_path=None,
                show=False,
                log_x=True,
                bootstrap_kind="parametric",
            )
            if axes and gp_posterior_mean_for_validation_plot is not None:
                gates_grid, ratio_grid, probabilities = gp_posterior_mean_for_validation_plot
                ax = axes[0]
                ax.set_title(f"{ax.get_title()} | seed={settings.rng_seed}")
                if float(np.min(probabilities)) <= contour_target(settings) <= float(np.max(probabilities)):
                    ax.contour(
                        gates_grid,
                        ratio_grid,
                        probabilities,
                        levels=[contour_target(settings)],
                        colors="tab:blue",
                        linewidths=2.2,
                        zorder=8,
                    )
                    ax.plot([], [], color="tab:blue", linewidth=2.2, label="GP mean p=0.5")
                    ax.legend(
                        loc="upper left",
                        bbox_to_anchor=(1.38, 1.0),
                        borderaxespad=0.0,
                        frameon=True,
                        framealpha=0.9,
                    )
            if axes:
                axes[0].figure.savefig(
                    validation_base_path.parent / f"{validation_base_path.stem}_surface.png",
                    dpi=200,
                    bbox_inches="tight",
                )

    if return_budget:
        total_budget = Budget(
            remaining_hqc=(
                validation_budget.remaining_hqc
                if validation_budget is not None
                else max(
                    0.0,
                    settings.hqc_budget - budget.spent_hqc,
                )
            ),
            spent_hqc=(
                budget.spent_hqc
                + (validation_budget.spent_hqc if validation_budget is not None else 0.0)
            ),
            jobs=(
                budget.jobs
                + (validation_budget.jobs if validation_budget is not None else 0)
            ),
            max_job_circuits=max(
                budget.max_job_circuits,
                validation_budget.max_job_circuits if validation_budget is not None else 0,
            ),
        )
        return rmb, crossings, total_budget, validation_rmb
    return rmb, crossings


def run_with_budget(settings: RickFantasyGPUSettings):
    """Run the modified crossing experiment and return the final budget."""
    return run(settings, return_budget=True)


# =============================================================================
# MAIN CONTROL PANEL
# =============================================================================


def main() -> None:
    """
    Main function.

    This is the only place you should need to edit for ordinary experiments.
    Treat this as the control panel / handle section.
    """

    # -------------------------------------------------------------------------
    # TARGET HANDLE
    # -------------------------------------------------------------------------

    TARGET_THRESHOLD = 0.5

    # -------------------------------------------------------------------------
    # SEARCH-BOX HANDLES
    # -------------------------------------------------------------------------
    # Keep these consistent with your RMB/CrossingSettings defaults.
    # Uncomment and edit if you want to override the parent defaults.

    N_GATES_BOUNDS = (10, 5000)
    RATIO_BOUNDS = (0.08, 0.98)

    # Example:
    # N_GATES_BOUNDS = (10, 180)
    # RATIO_BOUNDS = (0.2, 0.7)

    # -------------------------------------------------------------------------
    # HQC BUDGET HANDLES
    # -------------------------------------------------------------------------
    # If None, use the defaults inherited from CrossingSettings.

    HQC_BUDGET = 200.0
    MAX_COST_PER_RUN = 15.0

    # Example:
    # HQC_BUDGET = 500.0
    # MAX_COST_PER_RUN = 100.0

    # -------------------------------------------------------------------------
    # FAKE-ANCHOR HANDLES
    # -------------------------------------------------------------------------

    USE_FAKE_CORNERS = True
    EASY_CORNER_OUTCOME = 1
    HARD_CORNER_OUTCOME = 0
    EXTRA_FAKE_ANCHORS = []

    # -------------------------------------------------------------------------
    # SOBOL WARM-UP HANDLES
    # -------------------------------------------------------------------------

    INITIAL_SOBOL_SAMPLES = 10
    INITIAL_SOBOL_MAX_COST_PER_RUN = 30.0
    SOBOL_SCRAMBLE = True

    # -------------------------------------------------------------------------
    # VALIDATION HANDLES
    # -------------------------------------------------------------------------
    # After the GP is trained, measure only configs on its predicted contour
    # into a separate validation RMBData and plot/save that validation data.

    VALIDATE_GP_CONTOUR = False
    RESERVE_VALIDATION_BUDGET = False
    VALIDATION_HQC_BUDGET = 0.0
    VALIDATION_MAX_COST_PER_RUN = 30.0
    VALIDATION_TARGET_CONFIGS = 30
    VALIDATION_MIN_CONFIGS = 15
    VALIDATION_MAX_RESTARTS = 20
    VALIDATION_SHOTS_PER_CONFIG = 1
    VALIDATION_MAX_ATTEMPTS_MULTIPLIER = 10
    VALIDATION_PROBABILITY_BAND = (0.4, 0.6)
    VALIDATION_SEED_OFFSET = 271828
    VALIDATION_SAVE_PATH = None
    PLOT_GP_LEVEL_SET_RESULT = True
    SAVE_GP_PREDICTION_GRID = True
    PLOTS_MODULE = "plots_1"
    BAYESIAN_ESTIMATION_MODULE = "sympleq.core.bayesian_estimation"

    # -------------------------------------------------------------------------
    # ESTIMATOR HANDLES
    # -------------------------------------------------------------------------
    # Local to Rick modified; does not alter shared BayesianEstimator.default().

    ESTIMATOR_THRESHOLD = 0.0
    ESTIMATOR_MIN_RUNS = 1

    # -------------------------------------------------------------------------
    # GP / AEPSYCH HANDLES
    # -------------------------------------------------------------------------

    OPTIMIZATION_STEPS = 500
    INDUCING_SIZE = 150
    ACQUISITION_FUNCTION = "GlobalSUR"

    # These dominate suggestion time.
    # Reduce for quick tests.
    ACQUISITION_RESTARTS = 2 # Local optimization restarts for acquisition function optimization.
    ACQUISITION_SAMPLES = 300

    # -------------------------------------------------------------------------
    # BATCHING HANDLES
    # -------------------------------------------------------------------------

    BATCHING = True
    MAX_BATCH_SIZE = 80
    FANTASY_BATCHING = True

    # -------------------------------------------------------------------------
    # GPU HANDLES
    # -------------------------------------------------------------------------
    # RMB backend remains CPU.
    # This only controls AEPsych/GP/acquisition.

    USE_GPU = False

    # If CUDA device mismatch errors happen, set this to False.
    FORCE_DEFAULT_DEVICE_DURING_AEPSYCH = True

    # -------------------------------------------------------------------------
    # REPRODUCIBILITY / DEBUG HANDLES
    # -------------------------------------------------------------------------

    RNG_SEEDS = [2025, 2026, 2027, 2028, 2029, 2030]
    VERBOSE_FANTASIES = False
    PRINT_DIAGNOSTICS = False
    PLOT = True

    # -------------------------------------------------------------------------
    # Build settings.
    # -------------------------------------------------------------------------
    # We only pass optional inherited handles when the user explicitly sets them,
    # so this file remains compatible with your existing CrossingSettings defaults.

    kwargs = dict(
        target_threshold=TARGET_THRESHOLD,
        use_fake_corners=USE_FAKE_CORNERS,
        easy_corner_outcome=EASY_CORNER_OUTCOME,
        hard_corner_outcome=HARD_CORNER_OUTCOME,
        extra_fake_anchors=EXTRA_FAKE_ANCHORS,
        initial_sobol_samples=INITIAL_SOBOL_SAMPLES,
        initial_sobol_max_cost_per_run=INITIAL_SOBOL_MAX_COST_PER_RUN,
        sobol_scramble=SOBOL_SCRAMBLE,
        validate_gp_contour=VALIDATE_GP_CONTOUR,
        reserve_validation_budget=RESERVE_VALIDATION_BUDGET,
        validation_hqc_budget=VALIDATION_HQC_BUDGET,
        validation_max_cost_per_run=VALIDATION_MAX_COST_PER_RUN,
        validation_target_configs=VALIDATION_TARGET_CONFIGS,
        validation_min_configs=VALIDATION_MIN_CONFIGS,
        validation_max_restarts=VALIDATION_MAX_RESTARTS,
        validation_shots_per_config=VALIDATION_SHOTS_PER_CONFIG,
        validation_max_attempts_multiplier=VALIDATION_MAX_ATTEMPTS_MULTIPLIER,
        validation_probability_band=VALIDATION_PROBABILITY_BAND,
        validation_seed_offset=VALIDATION_SEED_OFFSET,
        validation_save_path=VALIDATION_SAVE_PATH,
        plot_gp_level_set_result=PLOT_GP_LEVEL_SET_RESULT,
        save_gp_prediction_grid=SAVE_GP_PREDICTION_GRID,
        plots_module=PLOTS_MODULE,
        bayesian_estimation_module=BAYESIAN_ESTIMATION_MODULE,
        estimator_threshold=ESTIMATOR_THRESHOLD,
        estimator_min_runs=ESTIMATOR_MIN_RUNS,
        optimization_steps=OPTIMIZATION_STEPS,
        inducing_size=INDUCING_SIZE,
        acquisition_function=ACQUISITION_FUNCTION,
        acquisition_restarts=ACQUISITION_RESTARTS,
        acquisition_samples=ACQUISITION_SAMPLES,
        batching=BATCHING,
        max_batch_size=MAX_BATCH_SIZE,
        fantasy_batching=FANTASY_BATCHING,
        use_gpu=USE_GPU,
        force_default_device_during_aepsych=FORCE_DEFAULT_DEVICE_DURING_AEPSYCH,
        verbose_fantasies=VERBOSE_FANTASIES,
        print_diagnostics=PRINT_DIAGNOSTICS,
        plot=PLOT,
    )

    if N_GATES_BOUNDS is not None:
        kwargs["n_gates_bounds"] = N_GATES_BOUNDS

    if RATIO_BOUNDS is not None:
        kwargs["ratio_bounds"] = RATIO_BOUNDS

    if HQC_BUDGET is not None:
        kwargs["hqc_budget"] = HQC_BUDGET

    if MAX_COST_PER_RUN is not None:
        kwargs["max_cost_per_run"] = MAX_COST_PER_RUN

    for run_index, rng_seed in enumerate(RNG_SEEDS, start=1):
        seed_kwargs = dict(kwargs)
        seed_kwargs["rng_seed"] = rng_seed
        seed_kwargs["save_path"] = timestamped_personal_save_path(seed=rng_seed)
        print(
            f"\n[seed run] {run_index}/{len(RNG_SEEDS)} "
            f"rng_seed={rng_seed} save_path={seed_kwargs['save_path']}\n"
        )
        settings = RickFantasyGPUSettings(**seed_kwargs)
        run(settings)

    if PLOT:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
