from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from scipy.optimize import minimize
from scipy.special import expit

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData

from viarregio2 import (
    config_from_parameters,
    fidelity_mean,
    fidelity_variance,
    make_backend,
    total_measurements,
)
from viarregio4 import fidelity_colormap
from viarregio7 import (
    BudgetState,
    ContourFirstExperimentConfig,
    MeasurementRequest,
    batch_config_key,
    batched_hqc_cost,
    max_cost_per_batch,
    measured_probability,
    native_quantinuum_hqc_cost,
    script_default_settings as script_default_settings_v7,
    spend_configs_batch,
    template_config,
)


@dataclass(frozen=True)
class VolumeExperimentConfig(ContourFirstExperimentConfig):
    """
    Coupled 3D HQC-aware fidelity-volume estimation settings.
    """
    n_qubits_values: tuple[int, ...] = (10, 15, 20, 25, 30, 35, 40, 45, 50)
    seed_n_qubits_values: tuple[int, ...] = (10, 30, 50)
    volume_grid_size: tuple[int, int, int] = (45, 37, 25)
    acquisition_ratio_count: int = 8
    acquisition_qubit_count: int = 9
    acquisition_passes: int = 8
    acquisition_candidates_per_pass: int = 27
    initial_seed_ratio_count: int = 6
    initial_seed_depth_grid_count: int = 13
    initial_low_ratio_depth_fractions: tuple[float, ...] = (0.0, 0.45, 0.6, 0.75, 0.9, 1.0)
    initial_seed_shots: int = 2
    qubit_depth_cap_enabled: bool = True
    qubit_depth_cap_plateau: int = 30
    qubit_depth_cap_power: float = 1.0
    qubit_depth_cap_margin_fraction: float = 0.0
    qubit_depth_cap_acquisition: bool = False
    min_volume_fit_points: int = 18
    volume_sparsity_radius: float = 0.12
    volume_uncertainty_weight: float = 1.0
    volume_sparsity_weight: float = 0.45
    volume_cost_power: float = 1.0
    volume_guardrails_enabled: bool = True
    volume_guardrail_min_successes_per_qubit: int = 1
    volume_guardrail_min_failures_per_qubit: int = 2
    volume_guardrail_shots: int = 2
    volume_guardrail_max_configs_per_pass: int = 14
    volume_guardrail_failure_ratios: tuple[float, ...] = (1.0, 0.72, 0.45)
    volume_guardrail_success_ratios: tuple[float, ...] = (0.0, 0.28, 0.55)
    bracket_completion_enabled: bool = True
    bracket_completion_max_groups_per_pass: int = 18
    bracket_completion_ratio_decimals: int = 2
    bracket_completion_probability_width: float = 0.32
    bracket_completion_depth_fractions: tuple[float, ...] = (0.08, 0.16, 0.28)
    volume_save_path: str | Path | None = "viarregio8_volume_grid.json"
    plot_save_path: str | Path | None = "viarregio8_volume.png"
    diagnostics_path: str | Path | None = "viarregio8_diagnostics.json"
    save_path: str | Path | None = "viarregio8_boundary.json"


def script_default_settings(**overrides) -> VolumeExperimentConfig:
    base = script_default_settings_v7(
        measurement_budget=1_000_000_000,
        hqc_budget=500.0,
        n_qubits_values=(10, 15, 20, 25, 30, 35, 40, 45, 50),
        depth_bounds=(4, 150),
        ratio_bounds=(0.08, 0.8),
        max_cost_per_batch=50.0,
        batch_max_configs=56,
        batch_fill_repeats=True,
        batch_fill_max_shots_per_config=8,
        batch_discovery_fill_max_shots_per_config=4,
        batch_fill_all_stages=False,
        contour_bracket_probe_shots=2,
        contour_bracket_depth_fractions=(0.04, 0.08, 0.12),
        contour_bracket_max_relative_depth=0.20,
        save_path=None,
        diagnostics_path=None,
        verbose=True,
    )
    params = dict(base.__dict__)
    params.update(
        n_qubits_values=(10, 15, 20, 25, 30, 35, 40, 45, 50),
        seed_n_qubits_values=(10, 30, 50),
        volume_save_path="viarregio8_volume_grid.json",
        plot_save_path="viarregio8_volume.png",
        diagnostics_path="viarregio8_diagnostics.json",
        save_path="viarregio8_boundary.json",
    )
    params.update(overrides)
    return VolumeExperimentConfig(**params)


@dataclass
class MonotoneFidelityVolume:
    alpha: float
    coefficients: np.ndarray
    depth_bounds: tuple[int, int]
    ratio_bounds: tuple[float, float]
    n_qubits_bounds: tuple[int, int]
    n_points: int

    def _scale(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        lower = np.array(
            [self.depth_bounds[0], self.ratio_bounds[0], self.n_qubits_bounds[0]],
            dtype=float,
        )
        upper = np.array(
            [self.depth_bounds[1], self.ratio_bounds[1], self.n_qubits_bounds[1]],
            dtype=float,
        )
        span = np.maximum(upper - lower, 1e-12)
        return np.clip((points - lower) / span, 0.0, 1.0)

    def features(self, points: np.ndarray) -> np.ndarray:
        scaled = self._scale(points)
        depth = scaled[:, 0]
        ratio = scaled[:, 1]
        qubits = scaled[:, 2]
        return np.column_stack([
            depth,
            ratio,
            qubits,
            depth**2,
            ratio**2,
            qubits**2,
            depth * ratio,
            depth * qubits,
            ratio * qubits,
            depth**3,
            ratio**3,
            qubits**3,
            depth * ratio * qubits,
        ])

    def probability(self, points: np.ndarray) -> np.ndarray:
        return expit(self.alpha - self.features(points) @ self.coefficients)


def volume_points_from_data(data: RMBData) -> tuple[list[RMBConfig], np.ndarray, np.ndarray, np.ndarray]:
    configs = [config for config, estimator in data.items() if estimator.num_runs() > 0]
    points = np.array(
        [
            [
                float(config.depth),
                float(config.min_two_qubit_gate_ratio),
                float(config.n_qubits),
            ]
            for config in configs
        ],
        dtype=float,
    )
    successes = []
    failures = []
    for config in configs:
        counts = data[config].counts()
        successes.append(float(counts.get(True, 0)))
        failures.append(float(counts.get(False, 0)))
    return configs, points, np.asarray(successes, dtype=float), np.asarray(failures, dtype=float)


def n_qubits_bounds(settings: VolumeExperimentConfig) -> tuple[int, int]:
    values = tuple(int(value) for value in settings.n_qubits_values)
    return min(values), max(values)


def estimated_max_depth_for_qubits(
    n_qubits: int,
    settings: VolumeExperimentConfig,
) -> float:
    """
    Qubit-dependent depth domain for the 3D benchmark.

    Defaults encode the observed scale: about depth 150 at 10 qubits, depth
    50 at 30 qubits, and about depth 30 at 50 qubits.
    """
    d_min, d_max = settings.depth_bounds
    if not settings.qubit_depth_cap_enabled:
        return float(d_max)
    reference_qubits = max(1.0, float(settings.qubit_depth_cap_plateau))
    q = max(1.0, float(n_qubits))
    reference_depth = min(float(d_max), 50.0)
    scaled_depth = reference_depth * (reference_qubits / q) ** max(0.0, settings.qubit_depth_cap_power)
    margin = max(0.0, settings.qubit_depth_cap_margin_fraction) * float(d_max - d_min)
    return float(np.clip(scaled_depth + margin, d_min, d_max))


def depth_domain_for_qubits(
    n_qubits: float,
    settings: VolumeExperimentConfig,
) -> tuple[float, float]:
    return float(settings.depth_bounds[0]), estimated_max_depth_for_qubits(int(round(n_qubits)), settings)


def volume_contour_bracket_depths(
    *,
    depth: float,
    n_qubits: int,
    settings: VolumeExperimentConfig,
) -> list[float]:
    d_min, d_max = depth_domain_for_qubits(n_qubits, settings)
    depth_span = max(1.0, float(d_max - d_min))
    candidate_depths = [float(depth)]
    for fraction in settings.contour_bracket_depth_fractions:
        span_offset = abs(float(fraction)) * depth_span
        relative_offset = (
            max(0.0, settings.contour_bracket_max_relative_depth)
            * max(1.0, float(depth))
        )
        offset = max(1.0, min(span_offset, max(1.0, relative_offset)))
        candidate_depths.extend([
            float(depth) - offset,
            float(depth) + offset,
        ])
    rounded_unique = sorted({round(float(value), 6) for value in candidate_depths})
    return [
        float(candidate_depth)
        for candidate_depth in rounded_unique
        if float(d_min) <= candidate_depth <= float(d_max)
    ]


def capped_contour_bracket_requests(
    *,
    config: RMBConfig,
    data: RMBData,
    settings: VolumeExperimentConfig,
    template: RMBConfig | None = None,
    seen: set[RMBConfig] | None = None,
    force_cap: bool = False,
) -> list[MeasurementRequest]:
    if seen is None:
        seen = set()
    if template is None:
        template = template_config(settings, config.n_qubits)

    ratio = float(config.min_two_qubit_gate_ratio)
    requests: list[MeasurementRequest] = []
    for depth in volume_contour_bracket_depths(
        depth=float(config.depth),
        n_qubits=config.n_qubits,
        settings=settings,
    ):
        bracket_config = config_from_parameters(
            template=template,
            depth=depth,
            ratio=ratio,
        )
        if bracket_config in seen:
            continue
        seen.add(bracket_config)
        current_runs = data[bracket_config].num_runs() if bracket_config in data else 0
        missing = max(0, settings.max_shots_per_config - current_runs)
        shots = min(max(1, settings.contour_bracket_probe_shots), missing)
        if shots > 0:
            requests.append(MeasurementRequest(bracket_config, shots))
    return requests


def fit_monotone_fidelity_volume(
    data: RMBData,
    settings: VolumeExperimentConfig,
    initial_surface: MonotoneFidelityVolume | None = None,
) -> MonotoneFidelityVolume:
    configs, points, successes, failures = volume_points_from_data(data)
    min_points = max(settings.min_volume_fit_points, settings.min_fit_points)
    if len(configs) < min_points:
        raise ValueError(f"Need at least {min_points} configs to fit 3D volume, got {len(configs)}.")

    probe = MonotoneFidelityVolume(
        alpha=0.0,
        coefficients=np.zeros(13, dtype=float),
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_qubits_bounds=n_qubits_bounds(settings),
        n_points=len(configs),
    )
    features = probe.features(points)

    def loss_and_grad(params: np.ndarray) -> tuple[float, np.ndarray]:
        alpha = params[0]
        coefficients = params[1:]
        z = alpha - features @ coefficients
        p = expit(z)
        eps = 1e-12
        loss = -np.sum(successes * np.log(p + eps) + failures * np.log(1.0 - p + eps))
        loss += 0.5 * settings.monotone_l2 * np.sum(coefficients**2)
        residual = p * (successes + failures) - successes
        grad_alpha = np.sum(residual)
        grad_coefficients = -features.T @ residual + settings.monotone_l2 * coefficients
        return loss, np.concatenate([[grad_alpha], grad_coefficients])

    if (
        initial_surface is not None
        and len(initial_surface.coefficients) == features.shape[1]
        and initial_surface.depth_bounds == settings.depth_bounds
        and initial_surface.ratio_bounds == settings.ratio_bounds
        and initial_surface.n_qubits_bounds == n_qubits_bounds(settings)
    ):
        x0 = np.concatenate([[initial_surface.alpha], initial_surface.coefficients])
    else:
        x0 = np.zeros(features.shape[1] + 1, dtype=float)
    bounds = [(None, None)] + [(0.0, None)] * features.shape[1]
    result = minimize(
        fun=lambda params: loss_and_grad(params)[0],
        x0=x0,
        jac=lambda params: loss_and_grad(params)[1],
        bounds=bounds,
        method="L-BFGS-B",
    )
    if not result.success:
        raise RuntimeError(result.message)

    return MonotoneFidelityVolume(
        alpha=float(result.x[0]),
        coefficients=result.x[1:],
        depth_bounds=settings.depth_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_qubits_bounds=n_qubits_bounds(settings),
        n_points=len(configs),
    )


def depth50_at(
    surface: MonotoneFidelityVolume,
    *,
    ratio: float,
    n_qubits: float,
    settings: VolumeExperimentConfig,
    steps: int = 28,
) -> tuple[float, str]:
    low, high = depth_domain_for_qubits(n_qubits, settings)
    p_low = float(surface.probability(np.array([[low, ratio, n_qubits]], dtype=float))[0])
    p_high = float(surface.probability(np.array([[high, ratio, n_qubits]], dtype=float))[0])
    if p_low < 0.5:
        return low, "all_below"
    if p_high >= 0.5:
        return high, "all_above"
    for _ in range(max(1, steps)):
        mid = 0.5 * (low + high)
        p_mid = float(surface.probability(np.array([[mid, ratio, n_qubits]], dtype=float))[0])
        if p_mid >= 0.5:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high), "crossing"


def extract_volume_grid(
    surface: MonotoneFidelityVolume,
    settings: VolumeExperimentConfig,
) -> dict:
    n_depths, n_ratios, n_qubits_count = settings.volume_grid_size
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], max(2, n_ratios))
    qubits = np.linspace(n_qubits_bounds(settings)[0], n_qubits_bounds(settings)[1], max(2, n_qubits_count))
    depth50 = np.zeros((len(qubits), len(ratios)), dtype=float)
    status = np.empty(depth50.shape, dtype=object)
    for q_index, q in enumerate(qubits):
        for r_index, ratio in enumerate(ratios):
            depth, state = depth50_at(surface, ratio=float(ratio), n_qubits=float(q), settings=settings)
            depth50[q_index, r_index] = depth
            status[q_index, r_index] = state

    d_min = float(settings.depth_bounds[0])
    domain_max = np.array(
        [depth_domain_for_qubits(float(q), settings)[1] for q in qubits],
        dtype=float,
    )[:, None]
    normalized_heights = np.clip(
        (depth50 - d_min) / np.maximum(1e-12, domain_max - d_min),
        0.0,
        1.0,
    )
    normalized_volume = float(np.mean(normalized_heights))
    return {
        "ratios": ratios.tolist(),
        "n_qubits": qubits.tolist(),
        "depth50": depth50.tolist(),
        "depth_max": domain_max[:, 0].tolist(),
        "status": status.tolist(),
        "normalized_volume": normalized_volume,
    }


def save_volume_grid(path: str | Path | None, grid: dict) -> None:
    if path is None:
        return
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(grid, indent=2), encoding="utf-8")


def existing_scaled_points(data: RMBData, settings: VolumeExperimentConfig) -> np.ndarray:
    configs, points, _, _ = volume_points_from_data(data)
    if not configs:
        return np.empty((0, 3), dtype=float)
    return np.asarray(
        [
            scaled_volume_point(
                float(config.depth),
                float(config.min_two_qubit_gate_ratio),
                float(config.n_qubits),
                settings,
            )
            for config in configs
        ],
        dtype=float,
    )


def scaled_volume_point(depth: float, ratio: float, n_qubits: float, settings: VolumeExperimentConfig) -> np.ndarray:
    d_min, d_max = depth_domain_for_qubits(n_qubits, settings)
    q_min, q_max = n_qubits_bounds(settings)
    return np.array([
        (float(depth) - d_min) / max(1e-12, d_max - d_min),
        (float(ratio) - settings.ratio_bounds[0]) / max(1e-12, settings.ratio_bounds[1] - settings.ratio_bounds[0]),
        (float(n_qubits) - q_min) / max(1e-12, q_max - q_min),
    ], dtype=float)


def initial_seed_qubits(settings: VolumeExperimentConfig) -> tuple[int, ...]:
    available = set(int(value) for value in settings.n_qubits_values)
    seeds = [int(value) for value in settings.seed_n_qubits_values if int(value) in available]
    if seeds:
        return tuple(dict.fromkeys(seeds))
    values = sorted(available)
    return tuple(dict.fromkeys([values[0], values[len(values) // 2], values[-1]]))


def low_ratio_anchor_depth(
    *,
    data: RMBData,
    n_qubits: int,
    settings: VolumeExperimentConfig,
) -> float:
    d_min, d_max = depth_domain_for_qubits(n_qubits, settings)
    ratio = float(settings.ratio_bounds[0])
    measured = [
        (config, measured_probability(data, config))
        for config in data
        if config.n_qubits == n_qubits
        and abs(float(config.min_two_qubit_gate_ratio) - ratio) < 1e-9
    ]
    measured = sorted(
        [(config, p) for config, p in measured if p is not None],
        key=lambda item: float(item[0].depth),
    )
    for (low_config, low_p), (high_config, high_p) in zip(measured, measured[1:]):
        if float(low_p) >= 0.5 and float(high_p) <= 0.5:
            low_depth = float(low_config.depth)
            high_depth = float(high_config.depth)
            span = max(1e-12, high_p - low_p)
            if abs(span) > 1e-12:
                fraction = np.clip((0.5 - low_p) / span, 0.0, 1.0)
                return float(low_depth + fraction * (high_depth - low_depth))
            return 0.5 * (low_depth + high_depth)
    if measured:
        return float(min(measured, key=lambda item: abs(float(item[1]) - 0.5))[0].depth)
    return float(d_min + 0.75 * (d_max - d_min))


def request_initial_low_ratio_seed(
    *,
    data: RMBData,
    settings: VolumeExperimentConfig,
) -> list[MeasurementRequest]:
    requests: list[MeasurementRequest] = []
    seen: set[RMBConfig] = set()
    for n_qubits in initial_seed_qubits(settings):
        template = template_config(settings, n_qubits)
        d_min, d_max = depth_domain_for_qubits(n_qubits, settings)
        fractions = tuple(settings.initial_low_ratio_depth_fractions)
        if not fractions:
            fractions = tuple(np.linspace(0.0, 1.0, max(2, settings.initial_seed_depth_grid_count)))
        low_ratio_depths = [
            d_min + np.clip(float(fraction), 0.0, 1.0) * (d_max - d_min)
            for fraction in fractions
        ]
        for depth in low_ratio_depths:
            config = config_from_parameters(template=template, depth=float(depth), ratio=settings.ratio_bounds[0])
            if config not in seen:
                seen.add(config)
                requests.append(MeasurementRequest(config, settings.initial_seed_shots))
    return requests


def request_initial_contour_seed(
    *,
    data: RMBData,
    settings: VolumeExperimentConfig,
) -> list[MeasurementRequest]:
    requests: list[MeasurementRequest] = []
    seen: set[RMBConfig] = set()
    ratios = np.linspace(
        settings.ratio_bounds[0],
        settings.ratio_bounds[1],
        max(2, settings.initial_seed_ratio_count),
    )
    for n_qubits in initial_seed_qubits(settings):
        template = template_config(settings, n_qubits)
        d_min, seed_depth_max = depth_domain_for_qubits(n_qubits, settings)

        anchor_depth = low_ratio_anchor_depth(
            data=data,
            n_qubits=n_qubits,
            settings=settings,
        )

        for ratio in ratios[1:]:
            predicted_depth = anchor_depth * (settings.ratio_bounds[0] / max(float(ratio), 1e-12)) ** 0.55
            predicted_depth = float(np.clip(predicted_depth, d_min, seed_depth_max))
            center = config_from_parameters(template=template, depth=predicted_depth, ratio=float(ratio))
            requests.extend(
                capped_contour_bracket_requests(
                    config=center,
                    data=data,
                    settings=settings,
                    template=template,
                    seen=seen,
                    force_cap=True,
                )
            )
    return requests


def request_initial_seed(
    *,
    data: RMBData,
    settings: VolumeExperimentConfig,
) -> list[MeasurementRequest]:
    return request_initial_low_ratio_seed(data=data, settings=settings) + request_initial_contour_seed(
        data=data,
        settings=settings,
    )


def candidate_score(
    *,
    config: RMBConfig,
    surface: MonotoneFidelityVolume,
    existing_points: np.ndarray,
    selected_points: list[np.ndarray],
    data: RMBData,
    settings: VolumeExperimentConfig,
) -> float:
    if config in data and data[config].num_runs() >= settings.max_shots_per_config:
        return 0.0
    point = np.array([[float(config.depth), float(config.min_two_qubit_gate_ratio), float(config.n_qubits)]])
    p = float(surface.probability(point)[0])
    uncertainty = max(1e-6, p * (1.0 - p))
    scaled = scaled_volume_point(
        float(config.depth),
        float(config.min_two_qubit_gate_ratio),
        float(config.n_qubits),
        settings,
    )
    if existing_points.size:
        sparsity = min(1.0, float(np.min(np.linalg.norm(existing_points - scaled, axis=1))) / settings.volume_sparsity_radius)
    else:
        sparsity = 1.0
    if selected_points:
        diversity = min(1.0, float(np.min(np.linalg.norm(np.asarray(selected_points) - scaled, axis=1))) / settings.volume_sparsity_radius)
    else:
        diversity = 1.0
    request = MeasurementRequest(config, max(1, settings.contour_bracket_probe_shots))
    cost = max(1e-12, batched_hqc_cost([request], settings))
    return float(
        (settings.volume_uncertainty_weight * uncertainty)
        * (1.0 + settings.volume_sparsity_weight * sparsity)
        * max(0.15, diversity)
        / (cost ** max(0.0, settings.volume_cost_power))
    )


def request_volume_acquisition(
    *,
    surface: MonotoneFidelityVolume,
    data: RMBData,
    settings: VolumeExperimentConfig,
) -> list[MeasurementRequest]:
    ratio_count = max(2, settings.acquisition_ratio_count)
    qubit_count = max(2, settings.acquisition_qubit_count)
    ratios = np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1], ratio_count)
    valid_qubits = np.asarray(settings.n_qubits_values, dtype=float)
    if qubit_count >= len(valid_qubits):
        qubits = valid_qubits
    else:
        qubits = np.linspace(n_qubits_bounds(settings)[0], n_qubits_bounds(settings)[1], qubit_count)
    existing_points = existing_scaled_points(data, settings)
    scored: list[tuple[float, RMBConfig]] = []
    for q_value in qubits:
        n_qubits = int(valid_qubits[np.argmin(np.abs(valid_qubits - q_value))])
        template = template_config(settings, n_qubits)
        for ratio in ratios:
            depth, state = depth50_at(surface, ratio=float(ratio), n_qubits=float(n_qubits), settings=settings)
            if state == "all_below":
                continue
            if settings.qubit_depth_cap_acquisition:
                depth = min(float(depth), estimated_max_depth_for_qubits(n_qubits, settings))
            config = config_from_parameters(template=template, depth=depth, ratio=float(ratio))
            score = candidate_score(
                config=config,
                surface=surface,
                existing_points=existing_points,
                selected_points=[],
                data=data,
                settings=settings,
            )
            if score > 0.0:
                scored.append((score, config))

    scored.sort(key=lambda item: item[0], reverse=True)
    selected: list[RMBConfig] = []
    selected_points: list[np.ndarray] = []
    seen_configs: set[RMBConfig] = set()
    for _, config in scored:
        if config in seen_configs:
            continue
        score = candidate_score(
            config=config,
            surface=surface,
            existing_points=existing_points,
            selected_points=selected_points,
            data=data,
            settings=settings,
        )
        if score <= 0.0:
            continue
        selected.append(config)
        seen_configs.add(config)
        selected_points.append(
            scaled_volume_point(
                float(config.depth),
                float(config.min_two_qubit_gate_ratio),
                float(config.n_qubits),
                settings,
            )
        )
        if len(selected) >= max(1, settings.acquisition_candidates_per_pass):
            break

    requests: list[MeasurementRequest] = []
    seen_requests: set[RMBConfig] = set()
    templates = {n: template_config(settings, n) for n in settings.n_qubits_values}
    for config in selected:
        requests.extend(
            capped_contour_bracket_requests(
                config=config,
                data=data,
                settings=settings,
                template=templates[config.n_qubits],
                seen=seen_requests,
            )
        )
    return requests


def request_bracket_completion(
    *,
    surface: MonotoneFidelityVolume,
    data: RMBData,
    settings: VolumeExperimentConfig,
) -> list[MeasurementRequest]:
    """
    Repair one-sided same-ratio groups by adding the missing side.

    Good contour learning needs, for many fixed (ratio, n_qubits) locations,
    measured depths on both sides of p=0.5. This pass looks at already sampled
    ratio/qubit groups and asks for a deeper point after all-success groups or
    a shallower point after all-failure groups, prioritizing groups near the
    fitted contour.
    """
    if not settings.bracket_completion_enabled:
        return []

    groups: dict[tuple[int, float], list[tuple[RMBConfig, float]]] = {}
    decimals = max(0, int(settings.bracket_completion_ratio_decimals))
    for config, estimator in data.items():
        if estimator.num_runs() <= 0:
            continue
        key = (
            int(config.n_qubits),
            round(float(config.min_two_qubit_gate_ratio), decimals),
        )
        groups.setdefault(key, []).append((config, fidelity_mean(estimator)))

    candidates: list[tuple[float, RMBConfig]] = []
    for (n_qubits, ratio_key), values in groups.items():
        if n_qubits not in settings.n_qubits_values:
            continue
        pass_depths = [float(config.depth) for config, p in values if p >= 0.5]
        fail_depths = [float(config.depth) for config, p in values if p < 0.5]
        if pass_depths and fail_depths:
            continue

        ratio_values = [float(config.min_two_qubit_gate_ratio) for config, _ in values]
        ratio = float(np.median(ratio_values))
        d_min, d_max = depth_domain_for_qubits(n_qubits, settings)
        predicted_depth, state = depth50_at(surface, ratio=ratio, n_qubits=float(n_qubits), settings=settings)
        if state == "all_below":
            continue
        predicted_p = float(surface.probability(np.array([[predicted_depth, ratio, n_qubits]], dtype=float))[0])
        if abs(predicted_p - 0.5) > settings.bracket_completion_probability_width:
            continue

        domain_span = max(1.0, d_max - d_min)
        measured_depths = [float(config.depth) for config, _ in values]
        template = template_config(settings, n_qubits)
        proposed_depths: list[float] = []
        if pass_depths and not fail_depths:
            base_depth = max(max(pass_depths), predicted_depth)
            for fraction in settings.bracket_completion_depth_fractions:
                proposed_depths.append(base_depth + abs(float(fraction)) * domain_span)
        elif fail_depths and not pass_depths:
            base_depth = min(min(fail_depths), predicted_depth)
            for fraction in settings.bracket_completion_depth_fractions:
                proposed_depths.append(base_depth - abs(float(fraction)) * domain_span)
        else:
            continue

        for proposed_depth in proposed_depths:
            depth = float(np.clip(proposed_depth, d_min, d_max))
            if any(abs(depth - existing_depth) < 0.5 for existing_depth in measured_depths):
                continue
            config = config_from_parameters(template=template, depth=depth, ratio=ratio)
            if config in data and data[config].num_runs() >= settings.max_shots_per_config:
                continue
            distance_to_contour = min(
                abs(float(config.depth) - predicted_depth) / domain_span,
                1.0,
            )
            group_weight = min(1.0, len(values) / 3.0)
            score = (1.0 - distance_to_contour) + 0.25 * group_weight
            candidates.append((score, config))
            break

    candidates.sort(key=lambda item: item[0], reverse=True)
    requests: list[MeasurementRequest] = []
    seen: set[RMBConfig] = set()
    for _, config in candidates:
        if config in seen:
            continue
        seen.add(config)
        current_runs = data[config].num_runs() if config in data else 0
        missing = max(0, settings.max_shots_per_config - current_runs)
        shots = min(max(1, settings.contour_bracket_probe_shots), missing)
        if shots > 0:
            requests.append(MeasurementRequest(config, shots))
        if len(requests) >= max(0, settings.bracket_completion_max_groups_per_pass):
            break
    return requests


def qubit_outcome_summary(data: RMBData, n_qubits: int) -> tuple[int, int, int]:
    successes = 0
    failures = 0
    configs = 0
    for config, estimator in data.items():
        if config.n_qubits != n_qubits or estimator.num_runs() <= 0:
            continue
        configs += 1
        p = fidelity_mean(estimator)
        if p >= 0.5:
            successes += 1
        else:
            failures += 1
    return configs, successes, failures


def ratio_from_fraction(fraction: float, settings: VolumeExperimentConfig) -> float:
    r_min, r_max = settings.ratio_bounds
    return float(r_min + np.clip(float(fraction), 0.0, 1.0) * (r_max - r_min))


def request_volume_guardrails(
    *,
    data: RMBData,
    settings: VolumeExperimentConfig,
) -> list[MeasurementRequest]:
    """
    Add sparse easy/hard anchors so the coupled 3D fit cannot drift from
    mostly one-sided evidence in a qubit slice.
    """
    if not settings.volume_guardrails_enabled:
        return []

    requests: list[tuple[float, MeasurementRequest]] = []
    seen: set[RMBConfig] = set()
    for n_qubits in settings.n_qubits_values:
        configs, successes, failures = qubit_outcome_summary(data, int(n_qubits))
        template = template_config(settings, int(n_qubits))
        d_min, d_max = depth_domain_for_qubits(int(n_qubits), settings)

        failure_deficit = max(0, settings.volume_guardrail_min_failures_per_qubit - failures)
        for offset, ratio_fraction in enumerate(settings.volume_guardrail_failure_ratios[:failure_deficit]):
            ratio = ratio_from_fraction(ratio_fraction, settings)
            config = config_from_parameters(template=template, depth=d_max, ratio=ratio)
            if config in seen or (
                config in data and data[config].num_runs() >= settings.max_shots_per_config
            ):
                continue
            seen.add(config)
            priority = 3.0 + failure_deficit + 0.1 * offset - 0.01 * configs
            requests.append((priority, MeasurementRequest(config, settings.volume_guardrail_shots)))

        success_deficit = max(0, settings.volume_guardrail_min_successes_per_qubit - successes)
        for offset, ratio_fraction in enumerate(settings.volume_guardrail_success_ratios[:success_deficit]):
            ratio = ratio_from_fraction(ratio_fraction, settings)
            config = config_from_parameters(template=template, depth=d_min, ratio=ratio)
            if config in seen or (
                config in data and data[config].num_runs() >= settings.max_shots_per_config
            ):
                continue
            seen.add(config)
            priority = 2.0 + success_deficit + 0.1 * offset - 0.01 * configs
            requests.append((priority, MeasurementRequest(config, settings.volume_guardrail_shots)))

    requests.sort(key=lambda item: item[0], reverse=True)
    return [
        request
        for _, request in requests[: max(0, settings.volume_guardrail_max_configs_per_pass)]
    ]


def add_diagnostic(diagnostics: list[dict] | None, event: str, **values) -> None:
    if diagnostics is None:
        return
    clean = {}
    for key, value in values.items():
        if isinstance(value, np.generic):
            value = value.item()
        clean[key] = value
    diagnostics.append({"event": event, **clean})


def estimate_boundary(settings: VolumeExperimentConfig) -> RMB:
    rng = default_rng(settings.rng_seed)
    backend = make_backend()
    rmb = RMB.default(rng).with_backend(backend)
    data: RMBData = rmb._data
    budget = BudgetState(
        remaining_measurements=settings.measurement_budget,
        remaining_hqc=float(settings.hqc_budget) if settings.hqc_budget is not None else float(settings.measurement_budget),
    )
    diagnostics: list[dict] = []

    seed_requests = request_initial_low_ratio_seed(data=data, settings=settings)
    spend_configs_batch(
        backend=backend,
        rng=rng,
        data=data,
        requests=seed_requests,
        budget=budget,
        settings=settings,
        diagnostics=diagnostics,
        stage="initial_grid",
    )
    seed_trace_requests = request_initial_contour_seed(data=data, settings=settings)
    spend_configs_batch(
        backend=backend,
        rng=rng,
        data=data,
        requests=seed_trace_requests,
        budget=budget,
        settings=settings,
        diagnostics=diagnostics,
        stage="initial",
    )

    fitted_surface: MonotoneFidelityVolume | None = None
    for pass_index in range(max(0, settings.acquisition_passes)):
        if not budget.can_spend():
            break
        guardrail_requests = request_volume_guardrails(data=data, settings=settings)
        completion_requests: list[MeasurementRequest] = []
        try:
            fitted_surface = fit_monotone_fidelity_volume(
                data,
                settings,
                initial_surface=fitted_surface,
            )
        except (RuntimeError, ValueError) as exc:
            if not guardrail_requests:
                add_diagnostic(diagnostics, "volume_fit_failed", pass_index=pass_index, reason=str(exc))
                break
            requests = guardrail_requests
        else:
            completion_requests = request_bracket_completion(
                surface=fitted_surface,
                data=data,
                settings=settings,
            )
            acquisition_requests = request_volume_acquisition(surface=fitted_surface, data=data, settings=settings)
            requests = guardrail_requests + completion_requests + acquisition_requests
        if not requests:
            add_diagnostic(diagnostics, "volume_acquisition_empty", pass_index=pass_index)
            break
        before = total_measurements(data)
        spends = spend_configs_batch(
            backend=backend,
            rng=rng,
            data=data,
            requests=requests,
            budget=budget,
            settings=settings,
            diagnostics=diagnostics,
            stage="boundary_acquisition",
        )
        after = total_measurements(data)
        add_diagnostic(
            diagnostics,
            "volume_acquisition",
            pass_index=pass_index,
            requested_configs=len(requests),
            guardrail_configs=len(guardrail_requests),
            completion_configs=len(completion_requests),
            spent=sum(spend.spent_shots for spend in spends),
            measurements=after,
            remaining_hqc=round(float(budget.remaining_hqc), 6),
        )
        if after == before:
            break

    if fitted_surface is None:
        try:
            fitted_surface = fit_monotone_fidelity_volume(
                data,
                settings,
                initial_surface=fitted_surface,
            )
        except (RuntimeError, ValueError):
            fitted_surface = None

    batch_saving = budget.native_hqc_estimate - budget.stitched_hqc_spent
    batch_saving_fraction = (
        batch_saving / budget.native_hqc_estimate
        if budget.native_hqc_estimate > 0.0
        else 0.0
    )
    rmb._batch_cost_summary = {
        "stitched_hqc_spent": budget.stitched_hqc_spent,
        "native_hqc_estimate": budget.native_hqc_estimate,
        "estimated_batching_saving_hqc": batch_saving,
        "estimated_batching_saving_fraction": batch_saving_fraction,
        "batched_jobs": budget.batched_jobs,
        "max_batched_job_size": budget.max_batched_job_size,
        "max_cost_per_batch": max_cost_per_batch(settings),
    }
    rmb._batch_config_ids = budget.batch_config_ids

    if fitted_surface is not None:
        grid = extract_volume_grid(fitted_surface, settings)
        rmb._volume_grid = grid
        save_volume_grid(settings.volume_save_path, grid)
    else:
        grid = None

    if settings.verbose:
        print(
            "3D RMB estimate: "
            f"{total_measurements(data)} outcomes, {len(data)} configs, "
            f"stitched HQC {budget.stitched_hqc_spent:.3f}"
        )
        if grid is not None:
            print(f"  normalized fidelity>=0.5 volume: {grid['normalized_volume']:.4f}")
        print(
            f"  batching saving: {batch_saving:.3f} "
            f"({batch_saving_fraction:.1%})"
        )

    if settings.save_path is not None:
        rmb.save(settings.save_path)
    if settings.diagnostics_path is not None:
        payload = {
            "settings": {
                "measurement_budget": settings.measurement_budget,
                "hqc_budget": settings.hqc_budget,
                "n_qubits_values": settings.n_qubits_values,
                "seed_n_qubits_values": settings.seed_n_qubits_values,
                "depth_bounds": settings.depth_bounds,
                "ratio_bounds": settings.ratio_bounds,
                "volume_grid_size": settings.volume_grid_size,
                "qubit_depth_cap_enabled": settings.qubit_depth_cap_enabled,
                "qubit_depth_cap_plateau": settings.qubit_depth_cap_plateau,
                "qubit_depth_cap_power": settings.qubit_depth_cap_power,
                "qubit_depth_cap_margin_fraction": settings.qubit_depth_cap_margin_fraction,
                "qubit_depth_cap_acquisition": settings.qubit_depth_cap_acquisition,
                "volume_guardrails_enabled": settings.volume_guardrails_enabled,
                "volume_guardrail_min_successes_per_qubit": settings.volume_guardrail_min_successes_per_qubit,
                "volume_guardrail_min_failures_per_qubit": settings.volume_guardrail_min_failures_per_qubit,
                "volume_guardrail_shots": settings.volume_guardrail_shots,
                "volume_guardrail_max_configs_per_pass": settings.volume_guardrail_max_configs_per_pass,
                "bracket_completion_enabled": settings.bracket_completion_enabled,
                "bracket_completion_max_groups_per_pass": settings.bracket_completion_max_groups_per_pass,
                "bracket_completion_ratio_decimals": settings.bracket_completion_ratio_decimals,
                "bracket_completion_probability_width": settings.bracket_completion_probability_width,
                "bracket_completion_depth_fractions": settings.bracket_completion_depth_fractions,
                "execution_backend": "local_sympleq",
                "cost_model": "quantinuum_stitched_estimate",
            },
            "total_measurements": total_measurements(data),
            "n_configs": len(data),
            "volume_grid": grid,
            "batch_config_ids": budget.batch_config_ids,
            **rmb._batch_cost_summary,
            "events": diagnostics,
        }
        path = Path(settings.diagnostics_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return rmb


def volume_metrics(
    estimated_grid: dict,
    reference_grid: dict | None = None,
    reference_surface: MonotoneFidelityVolume | None = None,
) -> dict:
    metrics = {
        "normalized_volume": float(estimated_grid["normalized_volume"]),
    }
    if reference_grid is not None:
        estimate_depths = np.asarray(estimated_grid["depth50"], dtype=float)
        reference_depths = np.asarray(reference_grid["depth50"], dtype=float)
        if estimate_depths.shape == reference_depths.shape:
            metrics["reference_normalized_volume"] = float(reference_grid["normalized_volume"])
            metrics["absolute_volume_error"] = abs(
                float(estimated_grid["normalized_volume"])
                - float(reference_grid["normalized_volume"])
            )
            metrics["mean_abs_depth50_error"] = float(np.mean(np.abs(estimate_depths - reference_depths)))
    if reference_surface is not None:
        ratios = np.asarray(estimated_grid["ratios"], dtype=float)
        qubits = np.asarray(estimated_grid["n_qubits"], dtype=float)
        depth50 = np.asarray(estimated_grid["depth50"], dtype=float)
        points = []
        for q_index, q in enumerate(qubits):
            for r_index, ratio in enumerate(ratios):
                points.append([depth50[q_index, r_index], ratio, q])
        probabilities = reference_surface.probability(np.asarray(points, dtype=float))
        metrics["calibration_error"] = float(np.mean(np.abs(probabilities - 0.5)))
    return metrics


def plot_volume_estimate_3d(
    rmb: RMB,
    settings: VolumeExperimentConfig,
    *,
    show: bool = True,
    save_path: str | Path | None = None,
):
    """
    Plot the v8 fidelity=0.5 contour as one 3D height-field.
    """
    cache_root = Path(tempfile.gettempdir()) / "sympleq_viarregio8_plot_cache"
    mpl_cache = cache_root / "matplotlib"
    xdg_cache = cache_root / "xdg"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    xdg_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_cache))

    import matplotlib.pyplot as plt

    grid = getattr(rmb, "_volume_grid", None)
    if grid is None:
        surface = fit_monotone_fidelity_volume(rmb._data, settings)
        grid = extract_volume_grid(surface, settings)
        rmb._volume_grid = grid

    ratios = np.asarray(grid["ratios"], dtype=float)
    qubits = np.asarray(grid["n_qubits"], dtype=float)
    depth50 = np.asarray(grid["depth50"], dtype=float)
    ratio_grid, qubit_grid = np.meshgrid(ratios, qubits)

    fig = plt.figure(figsize=(9.5, 7.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(
        depth50,
        ratio_grid,
        qubit_grid,
        color="#4c78a8",
        alpha=0.42,
        linewidth=0,
        antialiased=True,
    )

    configs = [
        config
        for config, estimator in rmb._data.items()
        if estimator.num_runs() > 0
    ]
    if configs:
        depths = np.asarray([config.depth for config in configs], dtype=float)
        point_ratios = np.asarray([config.min_two_qubit_gate_ratio for config in configs], dtype=float)
        point_qubits = np.asarray([config.n_qubits for config in configs], dtype=float)
        fidelities = np.asarray([fidelity_mean(rmb._data[config]) for config in configs], dtype=float)
        scatter = ax.scatter(
            depths,
            point_ratios,
            point_qubits,
            c=fidelities,
            cmap=fidelity_colormap(),
            vmin=0.0,
            vmax=1.0,
            s=22,
            edgecolors="black",
            linewidths=0.3,
        )
        fig.colorbar(
            scatter,
            ax=ax,
            shrink=0.72,
            pad=0.08,
            label="Measured mean fidelity",
        )

    ax.set_xlabel("Depth")
    ax.set_ylabel("Two-qubit gate ratio")
    ax.set_zlabel("n_qubits")
    ax.set_xlim(settings.depth_bounds)
    ax.set_ylim(settings.ratio_bounds)
    ax.set_zlim(n_qubits_bounds(settings))
    ax.set_title(
        "v8 fidelity = 0.5 contour\n"
        f"normalized fidelity>=0.5 volume = {float(grid['normalized_volume']):.4f}"
    )

    output = settings.plot_save_path if save_path is None else save_path
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=160)
        if settings.verbose:
            print(f"Saved 3D volume plot to {path}")
    if show:
        plt.show()
    return fig, ax


if __name__ == "__main__":
    settings = script_default_settings()
    rmb = estimate_boundary(settings)
    grid = getattr(rmb, "_volume_grid", None)
    if grid is not None:
        print(json.dumps(volume_metrics(grid), indent=2))
        plot_volume_estimate_3d(rmb, settings)
