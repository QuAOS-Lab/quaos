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
import logging
import re
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

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    MeasurementRequest,
    batch_hqc_cost,
    candidate_axes,
    config_points,
    measured_gate_count_bounds,
    print_experiment_summary,
    print_fit_reports,
    print_progress,
    save_crossings,
    spend_request_batch,
    start_run,
)


_ORIGINAL_AEPSYCH_TRANSFORM_OPTIONS = aepsych_parameter_transforms.transform_options
_ORIGINAL_AEPSYCH_STR_TO_LIST = Config._str_to_list
_ORIGINAL_AEPSYCH_STR_TO_ARRAY = Config._str_to_array


def timestamped_personal_save_path() -> Path:
    """Timestamped JSON output path under the repository's Personal folder."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return Path("Personal") / f"rick_fantasy_gpu_crossing_{timestamp}.json"


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

    # Explicit Sobol warm-up
    initial_sobol_samples: int = 10
    sobol_scramble: bool = False

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

    for n_gates, ratio, outcome in corners:
        obs = Observation(
            x_cpu=raw_point(n_gates, ratio, device=torch.device("cpu")),
            y=outcome,
            source="fake",
        )

        add_observation_to_strategy(strategy, obs, device=device)
        observations.append(obs)

        # Plotting uses realized one/two-qubit counts.
        config = settings.make_config(n_gates, ratio)
        results_for_plot.append(
            (
                float(config.n_1qb_gates),
                float(config.n_2qb_gates),
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
) -> tuple[list[RMBConfig], bool]:
    """
    Select as many candidates as fit in one stitched batch.

    Returns:
        selected, exhausted
    """

    selected: list[RMBConfig] = []

    for candidate in candidates:
        next_cost = batch_hqc_cost(one_shot_requests(selected + [candidate]))

        if next_cost > settings.max_cost_per_run:
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
                    float(config.n_1qb_gates),
                    float(config.n_2qb_gates),
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

    from sympleq.applications.randomized_benchmarking.experiments.plots import (
        gate_plane_points,
        level_set_grid,
    )

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


def level_set_configs(
    strategy: SequentialStrategy,
    settings: RickFantasyGPUSettings,
    *,
    device: torch.device,
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

    delta = (
        probabilities.detach().cpu().numpy().reshape(gates_grid.shape)
        - contour_target(settings)
    )

    configs: list[RMBConfig] = []
    seen: set[RMBConfig] = set()

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

        config = settings.make_config(n_gates, ratio)

        if config not in seen:
            seen.add(config)
            configs.append(config)

    return sorted(
        configs,
        key=lambda c: (c.ratio_2_qb_gates, c.n_gates),
    )


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

    print("[sobol warm-up]")
    print(f"  initial_sobol_samples                = {settings.initial_sobol_samples}")
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
    print("  RMB backend                          = CPU")

    print("[debug]")
    print(f"  rng_seed                             = {settings.rng_seed}")
    print(f"  verbose_fantasies                    = {settings.verbose_fantasies}")

    print("=============================================\n")


# =============================================================================
# Main run logic
# =============================================================================


def run(settings: RickFantasyGPUSettings) -> tuple[RMB, list[RMBConfig]]:
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
    print_fit_reports(data, settings)

    crossings = level_set_configs(
        strategy,
        settings,
        device=gp_device,
    )

    base_path = None
    if settings.save_path is not None:
        base_path = save_crossings(
            rmb,
            settings,
            budget,
            crossings,
        )

    # -------------------------------------------------------------------------
    # 5. Optional plotting
    # -------------------------------------------------------------------------

    if settings.plot:
        import matplotlib.pyplot as plt

        from sympleq.applications.randomized_benchmarking.experiments.plots import (
            plot_crossing_results,
            plot_gp_level_set,
        )

        plot_crossing_results(
            data,
            settings,
            crossings,
            base_path=base_path,
            show=False,
        )

        if strategy.model is not None:
            gp_png_path = None
            if base_path is not None:
                gp_png_path = base_path.parent / f"{base_path.stem}_gp.png"

            one_q_bounds, two_q_bounds = measured_gate_count_bounds(data)

            probabilities, latent_mean, latent_variance = predict_level_set(
                strategy,
                settings,
                one_q_bounds=one_q_bounds,
                two_q_bounds=two_q_bounds,
                device=gp_device,
            )

            plot_gp_level_set(
                probabilities,
                latent_mean,
                latent_variance,
                results_for_plot,
                one_q_bounds=one_q_bounds,
                two_q_bounds=two_q_bounds,
                target=contour_target(settings),
                png_path=gp_png_path,
                show=False,
            )

        plt.show()

    return rmb, crossings


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
    RATIO_BOUNDS = (0.08, 1.0)

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

    # -------------------------------------------------------------------------
    # SOBOL WARM-UP HANDLES
    # -------------------------------------------------------------------------

    INITIAL_SOBOL_SAMPLES = 30
    SOBOL_SCRAMBLE = False

    # -------------------------------------------------------------------------
    # GP / AEPSYCH HANDLES
    # -------------------------------------------------------------------------

    OPTIMIZATION_STEPS = 500
    INDUCING_SIZE = 150
    ACQUISITION_FUNCTION = "GlobalSUR"

    # These dominate suggestion time.
    # Reduce for quick tests.
    ACQUISITION_RESTARTS = 2 # Local optimization restarts for acquisition function optimization.
    ACQUISITION_SAMPLES = 300 #

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

    RNG_SEED = 2026
    VERBOSE_FANTASIES = True
    PLOT = True
    SAVE_PATH = timestamped_personal_save_path()

    # -------------------------------------------------------------------------
    # Build settings.
    # -------------------------------------------------------------------------
    # We only pass optional inherited handles when the user explicitly sets them,
    # so this file remains compatible with your existing CrossingSettings defaults.

    kwargs = dict(
        save_path=SAVE_PATH,
        target_threshold=TARGET_THRESHOLD,
        use_fake_corners=USE_FAKE_CORNERS,
        easy_corner_outcome=EASY_CORNER_OUTCOME,
        hard_corner_outcome=HARD_CORNER_OUTCOME,
        initial_sobol_samples=INITIAL_SOBOL_SAMPLES,
        sobol_scramble=SOBOL_SCRAMBLE,
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
        rng_seed=RNG_SEED,
        verbose_fantasies=VERBOSE_FANTASIES,
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

    settings = RickFantasyGPUSettings(**kwargs)

    run(settings)


if __name__ == "__main__":
    main()
