"""
Fantasy-batched GP level-set definitions for RMB.

What this file does
-------------------
Stores the definitions used for fantasy-batched GP level-set estimation of the
P(success) = target_threshold contour in:

    (total gates, two-qubit gate ratio)

space.

Important split
---------------
RMB backend:
    SympleQ / Quantinuum.
    For SympleQ, circuit execution remains CPU-only. This file does not try to
    GPU-enable RMB circuit execution.

AEPsych / GP / GlobalSUR:
    Uses CUDA if available and if settings.use_gpu=True.
"""


from __future__ import annotations

import copy
import re
import warnings
from contextlib import contextmanager
from typing import NamedTuple

import numpy as np
import torch
from aepsych.acquisition.objective import ProbitObjective
from aepsych.config import Config
from aepsych.generators.monotonic_thompson_sampler_generator import (
    MonotonicThompsonSamplerGenerator,
)
from aepsych.factory.monotonic import monotonic_mean_covar_factory
from aepsych.models.monotonic_rejection_gp import MonotonicRejectionGP
from aepsych.strategy import SequentialStrategy
from aepsych.utils import get_dims, get_optimizer_options
import aepsych.transforms.parameters as aepsych_parameter_transforms
from botorch.utils.sampling import draw_sobol_samples

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    MeasurementRequest,
    batch_hqc_cost,
    candidate_axes,
    config_points,
    print_progress,
    spend_request_batch,
    default_backend_factory,
    quantinuum_emulator_backend_factory,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.fantasy_levelset_settings import (
    FantasySettings,
)


_ORIGINAL_AEPSYCH_TRANSFORM_OPTIONS = aepsych_parameter_transforms.transform_options
_ORIGINAL_AEPSYCH_STR_TO_LIST = Config._str_to_list
_ORIGINAL_AEPSYCH_STR_TO_ARRAY = Config._str_to_array
_ORIGINAL_MONOTONIC_GP_FROM_CONFIG = MonotonicRejectionGP.from_config
_ORIGINAL_MONOTONIC_TS_FROM_CONFIG = MonotonicThompsonSamplerGenerator.from_config


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


## =============================================================================
# Basic helpers
# =============================================================================


def contour_target(settings: FantasySettings) -> float:
    """Return and validate the target contour value."""

    target = float(settings.target_threshold)
    if not 0.0 < target < 1.0:
        raise ValueError(f"target_threshold must be in (0, 1), got {target}.")
    return target


def gp_target(settings: FantasySettings) -> float:
    """
    Return the target value used by the monotone GP.

    Stored/user-facing contour is P(success) = target_threshold.
    The monotone GP is trained on P(failure), so its target is 1 - target.
    """

    return 1.0 - contour_target(settings)


def success_outcome_to_gp_outcome(y_success: int) -> int:
    """
    Convert stored RMB outcome to the outcome seen by the monotone GP.

    Stored RMB convention:
        1 = success
        0 = failure

    Monotone GP convention in this file:
        1 = failure
        0 = success

    This makes the GP's learned probability increase with total gates and ratio.
    """

    return 1 - int(y_success)


def gp_probability_to_success_probability(p_failure: float) -> float:
    """Convert GP-predicted P(failure) back to P(success)."""

    return 1.0 - float(p_failure)


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
    settings: FantasySettings,
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


def select_affordable_prefix(
    candidates: list[RMBConfig],
    settings: FantasySettings,
    budget: Budget,
    *,
    max_cost_per_run: float | None = None,
) -> tuple[list[RMBConfig], bool]:
    """
    Select as many candidates as fit in one stitched batch.
    This is for the initial SOBOL warm-up.
    Returns:
        selected, exhausted
    """

    selected: list[RMBConfig] = []
    cost_cap = settings.max_cost_per_run if max_cost_per_run is None else max_cost_per_run

    for candidate in candidates:
        trial_batch = selected + [candidate]
        next_cost = batch_hqc_cost(one_shot_requests(trial_batch))
        stitched_total_gates = sum(config.n_gates for config in trial_batch)

        if next_cost > cost_cap:
            break

        gate_budget = getattr(settings, "gate_budget", None)

        if settings.backend_factory is quantinuum_emulator_backend_factory:
            if stitched_total_gates > gate_budget:
                break

        if next_cost > budget.remaining_hqc:
            return selected, True

        selected.append(candidate)
        print(f"Stitched total gates: {stitched_total_gates}")

    return selected, False


# =============================================================================
# Fake anchors
# =============================================================================


def seed_fake_corners(
    strategy: SequentialStrategy,
    settings: FantasySettings,
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
    settings: FantasySettings,
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


# =============================================================================
# AEPsych strategy construction
# =============================================================================


def aepsych_config_string(settings: FantasySettings) -> str:
    """
    Build AEPsych config for the monotone GP phase only.

    Note:
        Sobol is handled explicitly by this file, not by AEPsych's init_strat.

    Important:
        AEPsych's monotone GP is monotone increasing. Physically, P(success)
        should decrease with total gates and two-qubit ratio. Therefore this
        file trains the GP on P(failure), then converts predictions back to
        P(success) wherever the rest of the RMB code expects it.
    """

    target = float(gp_target(settings))
    n_gates_lower = float(settings.n_gates_bounds[0])
    n_gates_upper = float(settings.n_gates_bounds[1])
    ratio_lower = float(settings.ratio_bounds[0])
    ratio_upper = float(settings.ratio_bounds[1])
    optimization_steps = int(settings.optimization_steps)
    inducing_size = int(settings.inducing_size)
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
    model = MonotonicRejectionGP
    generator = MonotonicThompsonSamplerGenerator

    [MonotonicRejectionGP]
    lb = [{n_gates_lower}, {ratio_lower}]
    ub = [{n_gates_upper}, {ratio_upper}]
    dim = 2
    monotonic_idxs = [0, 1]
    num_induc = {inducing_size}
    num_samples = 128
    num_rejection_samples = 40

    [MonotonicThompsonSamplerGenerator]
    target = {target}
    dim = 2
    num_samples = 1
    num_rejection_samples = 40
    num_ts_points = {acquisition_samples}
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


def patched_monotonic_thompson_from_config(
    cls,
    config: Config,
) -> MonotonicThompsonSamplerGenerator:
    """AEPsych version compatibility patch: pass the required ``dim`` argument."""

    classname = cls.__name__
    n_samples = config.getint(classname, "num_samples", fallback=1)
    n_rejection_samples = config.getint(
        classname,
        "num_rejection_samples",
        fallback=500,
    )
    num_ts_points = config.getint(classname, "num_ts_points", fallback=1000)
    target = config.getfloat(classname, "target", fallback=0.75)
    objective = config.getobj(classname, "objective", fallback=ProbitObjective)
    dim = config.getint(classname, "dim", fallback=2)
    explore_features = config.getlist(
        classname,
        "explore_idxs",
        element_type=int,
        fallback=None,
    )

    return cls(
        n_samples=n_samples,
        n_rejection_samples=n_rejection_samples,
        num_ts_points=num_ts_points,
        target_value=target,
        objective=objective,
        dim=dim,
        explore_features=explore_features,
    )


def patched_monotonic_thompson_gen(
    self: MonotonicThompsonSamplerGenerator,
    num_points: int,
    model: MonotonicRejectionGP,
    fixed_features=None,
    **kwargs,
) -> torch.Tensor:
    """AEPsych version compatibility patch: sample from the wrapped base model."""

    if fixed_features is not None:
        warnings.warn(
            "Cannot fix features when generating from MonotonicRejectionGenerator"
        )

    sample_model = getattr(model, "_base_obj", model)
    x = draw_sobol_samples(
        bounds=sample_model.bounds_,
        n=self.num_ts_points,
        q=1,
    ).squeeze(1)

    if self.explore_features is not None:
        for idx in self.explore_features:
            val = (
                sample_model.bounds_[0, idx]
                + torch.rand(1, device=x.device)
                * (sample_model.bounds_[1, idx] - sample_model.bounds_[0, idx])
            ).item()
            x[:, idx] = val

    f_samp = sample_model.sample(
        x,
        num_samples=self.n_samples,
        num_rejection_samples=self.n_rejection_samples,
    )
    dist = torch.abs(self.objective(f_samp) - self.target_value)
    if dist.ndim == 1:
        best_idx = torch.argmin(dist).reshape(1)
    else:
        best_idx = torch.argmin(dist, dim=1)
    return torch.Tensor(x[best_idx])


def patched_monotonic_gp_from_config(cls, config: Config) -> MonotonicRejectionGP:
    """AEPsych version compatibility patch: parse scalar counts as ints."""

    classname = cls.__name__
    num_induc = config.getint(classname, "num_induc", fallback=25)
    num_samples = config.getint(classname, "num_samples", fallback=250)
    num_rejection_samples = config.getint(
        classname,
        "num_rejection_samples",
        fallback=5000,
    )

    lb = config.gettensor(classname, "lb")
    ub = config.gettensor(classname, "ub")
    dim = get_dims(config)
    mean_covar_factory = config.getobj(
        classname,
        "mean_covar_factory",
        fallback=monotonic_mean_covar_factory,
    )
    mean, covar = mean_covar_factory(config)
    monotonic_idxs = config.getlist(classname, "monotonic_idxs", fallback=[-1])
    optimizer_options = get_optimizer_options(config, classname)

    return cls(
        monotonic_idxs=monotonic_idxs,
        lb=lb,
        ub=ub,
        dim=dim,
        num_induc=num_induc,
        num_samples=num_samples,
        num_rejection_samples=num_rejection_samples,
        mean_module=mean,
        covar_module=covar,
        optimizer_options=optimizer_options,
    )


def build_strategy(settings: FantasySettings) -> SequentialStrategy:
    """Build the AEPsych monotone GP/acquisition strategy."""

    # Make sure the config parser can resolve these class names from strings.
    for obj in (MonotonicRejectionGP, MonotonicThompsonSamplerGenerator):
        registered_names = getattr(Config, "registered_names", {})
        if obj.__name__ not in registered_names:
            Config.register_object(obj)

    config = Config()
    config.update(config_str=aepsych_config_string(settings))
    sanitize_aepsych_config_numbers(config)
    aepsych_parameter_transforms.transform_options = sanitized_aepsych_transform_options
    Config._str_to_list = sanitized_aepsych_str_to_list
    Config._str_to_array = sanitized_aepsych_str_to_array
    MonotonicRejectionGP.from_config = classmethod(patched_monotonic_gp_from_config)
    MonotonicThompsonSamplerGenerator.from_config = classmethod(
        patched_monotonic_thompson_from_config
    )
    MonotonicThompsonSamplerGenerator.gen = patched_monotonic_thompson_gen
    return SequentialStrategy.from_config(config)


def add_observation_to_strategy(
    strategy: SequentialStrategy,
    observation: Observation,
    *,
    device: torch.device,
) -> None:
    """
    Add one observation to an AEPsych strategy.

    The Observation object keeps the RMB convention:
        1 = success
        0 = failure

    The monotone GP receives the inverted convention:
        1 = failure
        0 = success
    """

    x = observation.x_cpu.to(device=device, dtype=torch.double)
    y_gp = success_outcome_to_gp_outcome(observation.y)
    strategy.add_data(x, [int(y_gp)])


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
# Device handling for AEPsych/GP/GlobalSUR (Kinda legacy; if GPU enabled ever)
# =============================================================================


def choose_gp_device(settings: FantasySettings) -> torch.device:
    """
    Choose device for AEPsych/GP/acquisition.

    RMB measurement on SympleQ remains CPU regardless of this.
    """

    if settings.use_gpu and torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


@contextmanager
def temporary_torch_default_device(device: torch.device, enabled: bool):
    """
    Temporarily make PyTorch create tensors on the chosen device.

    This helps AEPsych/BoTorch internals create acquisition tensors on CUDA.
    Redundant if run with CPU only. Does nothing if CUDA is unavailable.

    Turn off with in FantasySettings:
        force_default_device_during_aepsych=False
        It means: Do not temporarily change PyTorch’s default tensor device.

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
    settings: FantasySettings,
    observations: list[Observation],
    results_for_plot: list[tuple[float, float, int]],
    device: torch.device,
) -> None:
    """
    Measure selected configs using the RMB backend (SimpleQ/Emulator/Quantinuum Device).

    Only real outcomes from here are added to the real strategy.
    """

    if not selected:
        return

    # RMB measurement.
    outcomes_by_config = spend_request_batch(
        rmb.backend,  # To be set later
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
    settings: FantasySettings,
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
    settings: FantasySettings,
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
    settings: FantasySettings,
    *,
    device: torch.device,
) -> float:
    """
    Predict P(success | config).

    Internally the monotone GP predicts P(failure), because P(failure) should
    increase with total gates and two-qubit gate ratio. This function converts
    the result back to P(success) for the rest of the code.
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
                p_failure, _ = strategy.model.predict(x, probability_space=True)
            except RuntimeError:
                # Fallback for device mismatch.
                p_failure, _ = strategy.model.predict(
                    x.detach().cpu(),
                    probability_space=True,
                )

    p_failure_float = float(p_failure.detach().cpu().reshape(-1)[0])
    p_success = gp_probability_to_success_probability(p_failure_float)
    return min(max(p_success, 0.0), 1.0)


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
    settings: FantasySettings,
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
            if settings.verbose_fantasies:
                generator = fantasy_strategy._strat.generator
                base_generator = getattr(generator, "_base_obj", generator)
                rejection_samples = getattr(base_generator, "n_rejection_samples", None)
                ts_points = getattr(base_generator, "num_ts_points", None)
                print(
                    "[fantasy-gen] "
                    f"attempt={attempts:03d} "
                    f"accepted={len(selected):03d}/{batch_limit} "
                    f"rejection_samples={rejection_samples} "
                    f"ts_points={ts_points}"
                )
            x = fantasy_strategy.gen()
            if settings.verbose_fantasies:
                print("[fantasy-gen] generated candidate")

        candidate = config_from_aepsych_x(x, settings)

        # Rounding can cause duplicate realized configs.
        if candidate in seen:
            continue

        trial_batch = selected + [candidate]
        next_cost = batch_hqc_cost(one_shot_requests(trial_batch))
        stitched_total_gates = sum(config.n_gates for config in trial_batch)
        # print(f"Stitched total gates: {stitched_total_gates}")

        if next_cost > settings.max_cost_per_run:
            break

        gate_budget = getattr(settings, "gate_budget", None)

        if settings.backend_factory is quantinuum_emulator_backend_factory:
            if stitched_total_gates > gate_budget:
                break

        if next_cost > budget.remaining_hqc:
            return selected, True

        if settings.verbose_fantasies:
            print(
                "[fantasy-predict] "
                f"n_gates={candidate.n_gates} "
                f"ratio={candidate.ratio_2_qb_gates:.4f}"
            )
        p_success = predict_success_probability(
            fantasy_strategy,
            candidate,
            settings,
            device=device,
        )
        if settings.verbose_fantasies:
            print(f"[fantasy-predict] p_success={p_success:.4f}")

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
    settings: FantasySettings,
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


def level_set_configs(
    strategy: SequentialStrategy,
    settings: FantasySettings,
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

    failure_probability_grid = (
        probabilities.detach().cpu().numpy().reshape(gates_grid.shape)
    )

    # Convert the GP's monotone P(failure) surface back to P(success).
    probability_grid = 1.0 - failure_probability_grid

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
        print(
            f"  probability range on grid            = "
            f"{float(np.min(probability_grid)):.6g} to {float(np.max(probability_grid)):.6g}"
        )
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


# # =============================================================================
# # Final GP contour extraction
# # =============================================================================


# def predict_level_set(
#     strategy: SequentialStrategy,
#     settings: FantasySettings,
#     *,
#     one_q_bounds: tuple[int, int],
#     two_q_bounds: tuple[int, int],
#     device: torch.device,
# ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
#     """
#     Evaluate fitted GP on a grid for plotting.
#     """

#     plots = import_experiment_module(settings.plots_module)
#     gate_plane_points = plots.gate_plane_points
#     level_set_grid = plots.level_set_grid

#     one_q_grid, two_q_grid = level_set_grid(one_q_bounds, two_q_bounds)
#     points = gate_plane_points(one_q_grid, two_q_grid)

#     grid = torch.tensor(
#         np.column_stack(
#             [
#                 np.clip(
#                     points[:, 0],
#                     settings.n_gates_bounds[0],
#                     settings.n_gates_bounds[1],
#                 ),
#                 np.clip(
#                     points[:, 1],
#                     settings.ratio_bounds[0],
#                     settings.ratio_bounds[1],
#                 ),
#             ]
#         ),
#         dtype=torch.double,
#         device=device,
#     )

#     move_strategy_models_to_device(strategy, device)

#     with temporary_torch_default_device(
#         device,
#         enabled=settings.force_default_device_during_aepsych,
#     ):
#         with torch.no_grad():
#             probabilities, _ = strategy.model.predict(
#                 grid,
#                 probability_space=True,
#             )
#             latent_mean, latent_variance = strategy.model.predict(grid)

#     return (
#         probabilities.detach().cpu().numpy().reshape(one_q_grid.shape),
#         latent_mean.detach().cpu().numpy().reshape(one_q_grid.shape),
#         latent_variance.detach().cpu().numpy().reshape(one_q_grid.shape),
#     )


# def predict_native_level_set(
#     strategy: SequentialStrategy,
#     settings: FantasySettings,
#     *,
#     device: torch.device,
#     n_grid: int = 100,
# ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
#     """Evaluate fitted GP on its native total-gates/ratio mesh."""
#     gates_axis = np.geomspace(
#         max(1.0, float(settings.n_gates_bounds[0])),
#         float(settings.n_gates_bounds[1]),
#         n_grid,
#     )
#     ratio_axis = np.linspace(
#         float(settings.ratio_bounds[0]),
#         float(settings.ratio_bounds[1]),
#         n_grid,
#     )
#     gates_grid, ratio_grid = np.meshgrid(gates_axis, ratio_axis)
#     grid = torch.tensor(
#         np.column_stack([gates_grid.ravel(), ratio_grid.ravel()]),
#         dtype=torch.double,
#         device=device,
#     )

#     move_strategy_models_to_device(strategy, device)

#     with temporary_torch_default_device(
#         device,
#         enabled=settings.force_default_device_during_aepsych,
#     ):
#         with torch.no_grad():
#             probabilities, _ = strategy.model.predict(
#                 grid,
#                 probability_space=True,
#             )
#             latent_mean, latent_variance = strategy.model.predict(grid)

#     return (
#         probabilities.detach().cpu().numpy().reshape(gates_grid.shape),
#         latent_mean.detach().cpu().numpy().reshape(gates_grid.shape),
#         latent_variance.detach().cpu().numpy().reshape(gates_grid.shape),
#         gates_grid,
#         ratio_grid,
#     )
