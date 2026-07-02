"""
Fantasy-batched GP level-set definitions for RMB.

What this file does
-------------------
Stores the definitions used for fantasy-batched GP level-set estimation of the
P(success) = target_threshold contour in:

    (total gates, two-qubit gate ratio, n_qubits)

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
from contextlib import contextmanager
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
    MeasurementRequest,
    batch_hqc_cost,
    print_progress,
    spend_request_batch,
    default_backend_factory,
    quantinuum_emulator_backend_factory,
)
from fantasy_levelset_settings_3d import FantasySettings
from sympleq.integrations.quantinuum.utils import NATIVE_GATES_SET


_ORIGINAL_AEPSYCH_TRANSFORM_OPTIONS = aepsych_parameter_transforms.transform_options
_ORIGINAL_AEPSYCH_STR_TO_LIST = Config._str_to_list
_ORIGINAL_AEPSYCH_STR_TO_ARRAY = Config._str_to_array
_REPAIR_STATS = {"accepted_configs": 0, "repaired_configs": 0}


class Observation(NamedTuple):
    """
    One binary observation given to AEPsych.

    x_cpu:
        Shape (1, 3). Stored on CPU so it can safely be replayed into a
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
# Basic helpers h
# =============================================================================


def contour_target(settings: FantasySettings) -> float:
    """Return and validate the target contour value."""

    target = float(settings.target_threshold)
    if not 0.0 < target < 1.0:
        raise ValueError(f"target_threshold must be in (0, 1), got {target}.")
    return target


def point_from_config(config: RMBConfig, *, device: torch.device) -> torch.Tensor:
    """
    Convert RMBConfig to AEPsych coordinate.

    Returns shape:
        (1, 3)

    Coordinates:
        x[0, 0] = total gates
        x[0, 1] = two-qubit gate ratio
        x[0, 2] = number of qubits
    """

    return torch.tensor(
        [[
            float(config.n_gates),
            float(config.ratio_2_qb_gates),
            float(config.n_qubits),
        ]],
        dtype=torch.double,
        device=device,
    )


def n_qubits_bounds(settings: FantasySettings) -> tuple[int, int]:
    """Return the AEPsych qubit-coordinate bounds."""

    bounds = getattr(settings, "n_qubits_bounds", (5, 20))
    return int(bounds[0]), int(bounds[1])

def make_config_3d(
    settings: FantasySettings,
    n_gates: float,
    ratio: float,
    n_qubits: float,
) -> RMBConfig:
    """Build an RMBConfig from the 3D AEPsych coordinate."""

    q_min, q_max = n_qubits_bounds(settings)
    q = int(np.clip(round(float(n_qubits)), q_min, q_max))
    n_2qb_gates = max(0, 2 * round(float(ratio) * float(n_gates) / 2))
    n_1qb_gates = max(2 * q, 2 * round((float(n_gates) - n_2qb_gates) / 2))
    return (
        RMBConfig.default()
        .with_n_qubits(q)
        .with_n_1qb_gates(n_1qb_gates)
        .with_n_2qb_gates(n_2qb_gates)
        .with_random_elimination(settings.random_elimination)
        .with_gates_set(tuple(NATIVE_GATES_SET))
    )


def valid_config(config: RMBConfig) -> bool:
    """Return whether this config has enough 1Q gates for scrambler + inverse."""

    return int(config.n_1qb_gates) >= 2 * int(config.n_qubits)


def raw_validity(
    settings: FantasySettings,
    n_gates: float,
    ratio: float,
    n_qubits: float,
) -> tuple[bool, int, int, int]:
    """Check the raw AEPsych point before make_config_3d repairs it."""

    q_min, q_max = n_qubits_bounds(settings)
    q = int(np.clip(round(float(n_qubits)), q_min, q_max))
    n_2q = max(0, 2 * round(float(ratio) * float(n_gates) / 2))
    n_1q = 2 * round((float(n_gates) - n_2q) / 2)
    return n_1q >= 2 * q, n_1q, n_2q, q


def reset_repair_stats() -> None:
    _REPAIR_STATS["accepted_configs"] = 0
    _REPAIR_STATS["repaired_configs"] = 0


def record_repair(raw_is_valid: bool) -> None:
    _REPAIR_STATS["accepted_configs"] += 1
    if not raw_is_valid:
        _REPAIR_STATS["repaired_configs"] += 1


def repair_stats() -> dict[str, float | int]:
    accepted = int(_REPAIR_STATS["accepted_configs"])
    repaired = int(_REPAIR_STATS["repaired_configs"])
    return {
        "accepted_configs": accepted,
        "repaired_configs": repaired,
        "repair_fraction": float(repaired / accepted) if accepted else 0.0,
    }


def raw_point(
    n_gates: float,
    ratio: float,
    n_qubits: float,
    *,
    device: torch.device,
) -> torch.Tensor:
    """Create an AEPsych point directly from raw coordinates."""

    return torch.tensor(
        [[float(n_gates), float(ratio), float(n_qubits)]],
        dtype=torch.double,
        device=device,
    )

def qubit_band_from_first_selection(
    settings: FantasySettings,
    first_n_qubits: float,
) -> tuple[int, int]:
    """Build the per-batch qubit band from the first fantasy proposal."""

    q_min, q_max = n_qubits_bounds(settings)
    band_length = max(0, int(settings.qubit_band_length))
    band_start_max = max(q_min, q_max - band_length)
    band_start = int(np.clip(round(float(first_n_qubits)), q_min, band_start_max))
    band_end = min(q_max, band_start + band_length)
    return band_start, band_end


def config_from_aepsych_x_qubit_band(
    x: torch.Tensor,
    settings: FantasySettings,
    qubit_band: tuple[int, int],
) -> RMBConfig:
    """Convert AEPsych's point into an RMBConfig inside one qubit band."""

    x_cpu = x.detach().cpu()
    band_start, band_end = qubit_band
    n_qubits = int(np.clip(round(float(x_cpu[0, 2].item())), band_start, band_end))
    return make_config_3d(
        settings,
        float(x_cpu[0, 0].item()),
        float(x_cpu[0, 1].item()),
        float(n_qubits),
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
    This is for the initial SOBOL warm-up.
    Returns:
        selected, exhausted
    """

    selected: list[RMBConfig] = []
    qubit_band: tuple[int, int] | None = None
    cost_cap = settings.max_cost_per_run if max_cost_per_run is None else max_cost_per_run

    for candidate in candidates:
        if qubit_band is None:
            qubit_band = qubit_band_from_first_selection(settings, candidate.n_qubits)
        elif not qubit_band[0] <= candidate.n_qubits <= qubit_band[1]:
            continue

        trial_batch = selected + [candidate]
        next_cost = batch_hqc_cost(one_shot_requests(trial_batch))
        stitched_total_gates = sum(config.n_gates for config in trial_batch)

        if next_cost > cost_cap:
            break

        gate_budget = getattr(settings, "gate_budget", None)

        # if settings.backend_factory is default_backend_factory or settings.backend_factory is quantinuum_emulator_backend_factory:
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

    q_min, q_max = n_qubits_bounds(settings)
    corners = []
    fake_qubit_slices = max(1, int(settings.fake_qubit_slices))
    qubit_slices = np.rint(np.linspace(q_min, q_max, fake_qubit_slices)).astype(int)
    for q_slice in dict.fromkeys(int(q) for q in qubit_slices):
        corners.extend(
            [
                (
                    float(settings.n_gates_bounds[0]),
                    float(settings.ratio_bounds[0]),
                    float(q_slice),
                    int(settings.easy_corner_outcome),
                ),
                (
                    float(settings.n_gates_bounds[1]),
                    float(settings.ratio_bounds[1]),
                    float(q_slice),
                    int(settings.hard_corner_outcome),
                ),
            ]
        )
    for anchor in settings.extra_fake_anchors:
        if len(anchor) == 3:
            n_gates, ratio, outcome = anchor
            corners.append((n_gates, ratio, float(settings.n_qubits), outcome))
        else:
            n_gates, ratio, n_qubits, outcome = anchor
            corners.append((n_gates, ratio, n_qubits, outcome))

    for n_gates, ratio, n_qubits, outcome in corners:
        n_gates = float(np.clip(n_gates, *settings.n_gates_bounds))
        ratio = float(np.clip(ratio, *settings.ratio_bounds))
        n_qubits = float(np.clip(n_qubits, q_min, q_max))
        outcome = int(outcome)
        obs = Observation(
            x_cpu=raw_point(n_gates, ratio, n_qubits, device=torch.device("cpu")),
            y=outcome,
            source="fake",
        )

        add_observation_to_strategy(strategy, obs, device=device)
        observations.append(obs)

        # Plotting uses native AEPsych coordinates: total gates and ratio.
        config = make_config_3d(settings, n_gates, ratio, n_qubits)
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

    Sobol gives points in [0, 1]^2 for:
        [n_gates_min, n_gates_max] x [ratio_min, ratio_max]

    The qubit coordinate is assigned as distinct integer slices so Sobol can be
    measured as separate same-qubit stitched submissions.
    """

    if settings.initial_sobol_samples <= 0:
        return []

    q_min, q_max = n_qubits_bounds(settings)
    band_length = max(0, int(settings.qubit_band_length))
    band_start_max = max(q_min, q_max - band_length)

    n_bands = max(1, int(settings.sobol_band_batches))
    band_starts = np.rint(np.linspace(q_min, band_start_max, n_bands)).astype(int)
    qubit_values = list(dict.fromkeys(int(q) for q in band_starts))
    n_bands = len(qubit_values)

    total_samples = int(settings.initial_sobol_samples)
    samples_per_band = max(1, int(np.ceil(total_samples / n_bands)))

    engine = torch.quasirandom.SobolEngine(
        dimension=2,
        scramble=bool(settings.sobol_scramble),
        seed=settings.rng_seed if settings.sobol_scramble else None,
    )

    unit_points = engine.draw(samples_per_band).cpu().numpy()

    n_min, n_max = settings.n_gates_bounds
    r_min, r_max = settings.ratio_bounds

    candidates: list[RMBConfig] = []
    seen: set[RMBConfig] = set()

    for n_qubits in qubit_values:
        for u_n, u_r in unit_points:
            n_gates = float(n_min + u_n * (n_max - n_min))
            ratio = float(r_min + u_r * (r_max - r_min))

            raw_is_valid, _, _, _ = raw_validity(
                settings,
                n_gates,
                ratio,
                n_qubits,
            )
            config = make_config_3d(settings, n_gates, ratio, n_qubits)

            # Realized gate-count rounding can create duplicates.
            if config not in seen:
                record_repair(raw_is_valid)
                candidates.append(config)
                seen.add(config)

            if len(candidates) >= total_samples:
                return candidates

    return candidates


# =============================================================================
# AEPsych strategy construction
# =============================================================================


def aepsych_config_string(settings: FantasySettings) -> str:
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
    qubits_lower, qubits_upper = n_qubits_bounds(settings)
    optimization_steps = int(settings.optimization_steps)
    inducing_size = int(settings.inducing_size)
    acquisition_restarts = int(settings.acquisition_restarts)
    acquisition_samples = int(settings.acquisition_samples)

    return f"""
    [common]
    parnames = [n_gates, ratio, n_qubits]
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

    [n_qubits]
    par_type = continuous
    lower_bound = {float(qubits_lower)}
    upper_bound = {float(qubits_upper)}

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
    stimuli_per_trial = 1

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


def build_strategy(settings: FantasySettings) -> SequentialStrategy:
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
    qubit_band: tuple[int, int] | None = None

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

        x_cpu = x.detach().cpu()
        raw_n_gates = float(x_cpu[0, 0].item())
        raw_ratio = float(x_cpu[0, 1].item())
        raw_n_qubits = float(x_cpu[0, 2].item())

        if qubit_band is None:
            qubit_band = qubit_band_from_first_selection(settings, raw_n_qubits)

        candidate = config_from_aepsych_x_qubit_band(
            x,
            settings,
            qubit_band,
        )
        raw_n_qubits = float(candidate.n_qubits)

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
        raw_is_valid, raw_n_1q, raw_n_2q, raw_q = raw_validity(
            settings,
            raw_n_gates,
            raw_ratio,
            raw_n_qubits,
        )
        record_repair(raw_is_valid)

        if settings.verbose_fantasies:
            is_valid = valid_config(candidate)
            print(
                "[fantasy] "
                f"batch_index={len(selected):03d} "
                f"raw_gates={raw_n_gates:.4g} "
                f"raw_ratio={raw_ratio:.4f} "
                f"raw_q={raw_q} "
                f"raw_valid={raw_is_valid} "
                f"n_gates={candidate.n_gates} "
                f"n_qubits={candidate.n_qubits} "
                f"ratio={candidate.ratio_2_qb_gates:.4f} "
                f"valid={is_valid} "
                f"p_success={p_success:.4f} "
                f"virtual={virtual_outcome} "
                f"batch_cost={next_cost:.3f}"
            )
            if not raw_is_valid:
                print(
                    "[fantasy repaired config] "
                    f"raw_n_1q={raw_n_1q} "
                    f"raw_n_2q={raw_n_2q} "
                    f"chosen_n_1q={candidate.n_1qb_gates} "
                    f"chosen_n_2q={candidate.n_2qb_gates} "
                    f"chosen_gates={candidate.n_gates} "
                    f"chosen_ratio={candidate.ratio_2_qb_gates:.4f} "
                    f"chosen_q={candidate.n_qubits}"
                )
            if not is_valid:
                print(
                    "[fantasy invalid sent config] "
                    f"n_1q={candidate.n_1qb_gates} "
                    f"n_2q={candidate.n_2qb_gates} "
                    f"n_gates={candidate.n_gates} "
                    f"ratio={candidate.ratio_2_qb_gates:.4f} "
                    f"n_qubits={candidate.n_qubits}"
                )

    return selected, False
