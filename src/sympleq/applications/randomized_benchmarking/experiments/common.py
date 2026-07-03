"""
Utilities shared by the boundary experiments in this package.

Configs live in (total gates, two-qubit gate ratio) space at fixed
``n_qubits``. Every probe is metered in Quantinuum credits (HQC) with the
base submission cost paid once per stitched batch. By default circuits are
simulated locally with SympleQ; ``CrossingSettings.backend_factory``
switches a run to another backend, e.g. real Quantinuum hardware.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable

import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng
from scipy.optimize import minimize
from scipy.special import expit

from sympleq.applications.randomized_benchmarking.RMB import RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.backends.base import MeasurementRequest, RMBBackend
from sympleq.applications.randomized_benchmarking.backends.quantinuum import QuantinuumBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.integrations.quantinuum.utils import (
    BASE_SIMULATION_COST,
    NATIVE_GATES_SET,
    pytket_bare_simulation_cost,
    to_pytket_circuit,
)


def default_backend_factory(settings: CrossingSettings, rng: RNGGenerator) -> RMBBackend:
    """SympleQ emulation of Quantinuum hardware; nothing is submitted."""
    return QuantinuumBackend.default_sympleq_backend(rng)


def quantinuum_emulator_backend_factory(settings: CrossingSettings, rng: RNGGenerator) -> RMBBackend:
    """SympleQ emulation of Quantinuum hardware; nothing is submitted."""
    return QuantinuumBackend(
        device_name="H2-Emulator",
        project_name="level-benchmark",
        batch_size=1,
        max_cost_per_run=settings.max_cost_per_run
    )

def quantinuum_H2_backend_factory(settings: CrossingSettings, rng: RNGGenerator) -> RMBBackend:
    """SympleQ emulation of Quantinuum hardware; nothing is submitted."""
    return QuantinuumBackend(
        device_name="H2-1",
        project_name="level-benchmark",
        batch_size=1,
        max_cost_per_run=settings.max_cost_per_run
    )

@dataclass(frozen=True)
class CrossingSettings:
    """
    Problem definition shared by the fidelity = 0.5 crossing experiments.

    Each experiment subclasses this with its own method knobs, so different
    strategies can be run and compared on the same problem with the same
    budget, data handling, plots, and status reporting.

    Parameters
    ----------
    n_qubits : int
        Fixed number of qubits for the whole run.
    random_elimination : float
        ``RMBConfig.random_elimination`` used for every probed config.
    use_scrambler : bool
        ``RMBConfig.use_scrambler`` used for every probed config. Disabling it
        removes the ``2 * n_qubits`` one-qubit gate floor from config creation.
    n_gates_bounds : tuple[int, int]
        Search range for the total gate count.
    ratio_bounds : tuple[float, float]
        Search range for the two-qubit gate ratio.
    hqc_budget : float
        Total Quantinuum credits the run may spend.
    max_cost_per_run : float
        HQC cap of one stitched submission, mirroring
        ``QuantinuumBackend.max_cost_per_run``. It sets how many probe
        circuits share one base submission cost.
    monotone_l2 : float
        L2 regularization of the monotone fidelity-surface fit.
    min_fit_points : int
        Minimum measured configs required to fit a surface.
    candidate_grid_size : tuple[int, int]
        (gates, ratio) grid resolution used by surface evaluation and plots.
    rng_seed : int | None
        Seed for the run RNG. Each measurement additionally derives its own
        stream from (seed, config, shot index), so two experiments run with
        the same seed record identical outcomes at any config they both
        measure. ``None`` uses fresh entropy.
    backend_factory : Callable[[CrossingSettings, RNGGenerator], RMBBackend]
        Builds the backend the run measures on, from the settings and the
        run rng. The default simulates locally with the SympleQ emulation of
        Quantinuum hardware; pass a factory returning a
        :class:`~..backends.quantinuum.QuantinuumBackend` to submit to real
        hardware (align its ``max_cost_per_run`` with the settings' so the
        planned batches match what the backend stitches).
    save_path : str | Path | None
        Where to save the RMB data (see :meth:`RMB.save`). ``None`` skips
        saving. Crossings go to a sibling ``*_crossings.json`` file and
        plots to sibling ``*.png`` files.
    plot : bool
        Show the data and the traced line at the end of the run.
    verbose : bool
        Print progress while running.
    scatter_merge_bins : tuple[int, int] | None
        ``(one-qubit, two-qubit)`` gate-count bin sizes used to coarsen the
        scatter plots, merging the estimators of configs that fall in the same
        bin. ``None`` plots every measured config at full resolution; use it
        for densely-sampled runs whose points the default bins would collapse.
    """
    n_qubits: int = 5
    random_elimination: float = 0.1
    use_scrambler: bool = True
    n_gates_bounds: tuple[int, int] = (10, 3000)
    ratio_bounds: tuple[float, float] = (0.1, 0.9)
    hqc_budget: float = 250.0
    max_cost_per_run: float = 15.0
    monotone_l2: float = 1e-3
    min_fit_points: int = 16
    candidate_grid_size: tuple[int, int] = (50, 50)
    rng_seed: int | None = 1234
    backend_factory: Callable[[CrossingSettings, RNGGenerator], RMBBackend] = \
        default_backend_factory
    save_path: str | Path | None = None
    plot: bool = True
    verbose: bool = True
    scatter_merge_bins: tuple[int, int] | None = (20, 10)

    def make_config(self, n_gates: float, ratio: float) -> RMBConfig:
        """
        Build a valid native-gate RMBConfig from total gates and two-qubit ratio.

        Gate counts are rounded to the even values the circuit construction
        actually realizes.  With the scrambler enabled there must be at least
        the 2 * n_qubits scrambler one-qubit gates, so the config matches the
        generated circuits and the Quantinuum cost is computed for what would
        really run.
        """
        n_2qb_gates = max(0, 2 * round(ratio * n_gates / 2))
        min_1qb_gates = 2 * self.n_qubits if self.use_scrambler else 0
        n_1qb_gates = max(min_1qb_gates, 2 * round((n_gates - n_2qb_gates) / 2))
        return (
            RMBConfig.default()
            .with_n_qubits(self.n_qubits)
            .with_n_1qb_gates(n_1qb_gates)
            .with_n_2qb_gates(n_2qb_gates)
            .with_random_elimination(self.random_elimination)
            .with_use_scrambler(self.use_scrambler)
            .with_gates_set(tuple(NATIVE_GATES_SET))
        )


@dataclass
class Budget:
    """Remaining and spent Quantinuum credits, counting stitched submissions."""
    remaining_hqc: float
    spent_hqc: float = 0.0
    jobs: int = 0
    max_job_circuits: int = 0

    def spend(self, cost: float) -> None:
        self.remaining_hqc -= cost
        self.spent_hqc += cost

    def spend_batch(self, cost: float, n_circuits: int) -> None:
        """Charge one stitched submission of ``n_circuits`` circuits."""
        self.spend(cost)
        self.jobs += 1
        self.max_job_circuits = max(self.max_job_circuits, n_circuits)

    def can_afford(self, cost: float) -> bool:
        return self.remaining_hqc >= cost


def start_run(settings: CrossingSettings) -> tuple[RNGGenerator, RMB, Budget]:
    """Seeded RNG, RMB on the settings' backend, and budget for one run."""
    rng = default_rng(settings.rng_seed)
    rmb = RMB.default(rng).with_backend(settings.backend_factory(settings, rng))
    return rng, rmb, Budget(remaining_hqc=settings.hqc_budget)


@lru_cache(maxsize=None)
def single_circuit_bare_hqc(config: RMBConfig) -> float:
    """
    Bare Quantinuum credits for one circuit of this config at one shot,
    excluding the per-submission base cost but including the qubit resets
    needed to stitch it after another circuit.

    The cost depends only on the gate counts (the pricing circuit is drawn
    with a fixed rng), so the result is memoized per config.
    """
    circuit = to_pytket_circuit(config.random_circuit(rng=default_rng(0)))
    return pytket_bare_simulation_cost(circuit) + config.n_qubits / 5000


def stitched_batch_hqc(total_bare_hqc: float) -> float:
    """
    Quantinuum credits for one stitched submission with this total bare cost.

    Mirrors ``QuantinuumBackend.fidelity_estimation``, which stitches circuits
    into a single program, so the base submission cost is paid once per batch.
    """
    return BASE_SIMULATION_COST + total_bare_hqc


def batch_hqc_cost(requests: list[MeasurementRequest]) -> float:
    """HQC cost of one stitched submission holding all the requested circuits."""
    requests = [request for request in requests if request.shots > 0]
    if not requests:
        return 0.0
    return stitched_batch_hqc(sum(request.shots * single_circuit_bare_hqc(request.config)
                                  for request in requests))


def bare_hqc_at_register_width(config: RMBConfig, register_width: int) -> float:
    """Bare HQC for the real circuit, reset/register-priced at width W.

    The circuit gates stay those of ``config``; only the per-shot
    register/reset term is lifted from q/5000 to W/5000.  This is useful when
    a stitched batch runs on a physical register whose width is set by the
    largest circuit in the batch, while smaller circuits leave qubits idle.
    """
    width = max(int(register_width), int(config.n_qubits))
    return single_circuit_bare_hqc(config) + (width - int(config.n_qubits)) / 5000


def batch_hqc_cost_at_physical_width(requests: list[MeasurementRequest]) -> float:
    """HQC cost with every circuit priced at the batch's physical register width."""
    requests = [request for request in requests if request.shots > 0]
    if not requests:
        return 0.0
    width = max(int(request.config.n_qubits) for request in requests)
    total_bare = sum(
        request.shots * bare_hqc_at_register_width(request.config, width)
        for request in requests
    )
    return stitched_batch_hqc(total_bare)


def marginal_hqc_cost(config: RMBConfig, shots: int) -> float:
    """Marginal stitched cost of one config's shots inside a non-empty batch."""
    return max(1e-12, shots * single_circuit_bare_hqc(config))


def stitch_batch_size(bare_hqc: float, max_cost_per_run: float, max_circuits: int) -> int:
    """Number of circuits at this bare cost that fit in one stitched submission."""
    if bare_hqc <= 0:
        return max_circuits
    affordable = int((max_cost_per_run - BASE_SIMULATION_COST) // bare_hqc)
    return max(1, min(affordable, max_circuits))


def measurement_rng(seed: int, config: RMBConfig, shot_index: int) -> RNGGenerator:
    """RNG of one measurement, derived from the seed, the config, and the shot index."""
    return default_rng([seed, config.n_qubits, config.n_1qb_gates,
                        config.n_2qb_gates, shot_index])


def spend_request_batch(backend, rng: RNGGenerator, data: RMBData,
                        requests: list[MeasurementRequest],
                        seed: int | None = None) -> dict[RMBConfig, list[bool]]:
    """
    Record the requested Boolean fidelity outcomes into ``data``.

    The whole batch goes to the backend in one call, so backends that stitch
    circuits into device submissions see the full batch. With ``seed``, every
    shot draws from :func:`measurement_rng` keyed by the config's total
    recorded shots, so a config's outcome stream does not depend on how its
    measurements are grouped into batches or on the order in which configs
    are visited: different experiments run with the same seed see identical
    outcomes at shared configs (common random numbers).
    """
    requests = [request for request in requests if request.shots > 0]
    if not requests:
        return {}
    offsets: dict[RMBConfig, int] = {}
    for request in requests:
        estimator = data.setdefault(request.config, BayesianEstimator.default())
        offsets.setdefault(request.config, estimator.num_runs())

    shot_rng = None
    if seed is not None:
        def shot_rng(config: RMBConfig, index: int) -> RNGGenerator:
            return measurement_rng(seed, config, offsets[config] + index)

    outcomes = backend.fidelity_estimation(requests, rng, shot_rng=shot_rng).outcomes
    for config, results in outcomes.items():
        for outcome in results:
            data[config].record(bool(outcome))
    return outcomes


def measured_items(data: RMBData) -> list[tuple[RMBConfig, BayesianEstimator]]:
    """(config, estimator) pairs with at least one recorded outcome."""
    return [(config, estimator) for config, estimator in data.items()
            if estimator.num_runs() > 0]


def measured_gate_count_bounds(
    data: RMBData,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """
    (one-qubit, two-qubit) gate-count box spanned by the measured configs.

    Degenerate axes are widened to at least 2 gates so the box can hold a
    plotting grid.
    """
    measured = measured_items(data)
    one_q_values = [config.n_1qb_gates for config, _ in measured]
    two_q_values = [config.n_2qb_gates for config, _ in measured]
    return (
        (min(one_q_values), max(max(one_q_values), min(one_q_values) + 2)),
        (min(two_q_values), max(max(two_q_values), min(two_q_values) + 2)),
    )


def grouped_by_n_qubits(data: RMBData) -> dict[int, RMBData]:
    """Group measured configs by ``n_qubits``, dropping empty estimators."""
    groups: dict[int, RMBData] = {}
    for config, estimator in data.items():
        if estimator.num_runs() == 0:
            continue
        groups.setdefault(config.n_qubits, {})[config] = estimator
    return groups


def config_points(configs: list[RMBConfig]) -> np.ndarray:
    """(total gates, two-qubit ratio) coordinates of configs as a float array."""
    return np.array(
        [[float(config.n_gates), float(config.ratio_2_qb_gates)] for config in configs],
        dtype=float,
    )


def candidate_axes(settings: CrossingSettings) -> tuple[np.ndarray, np.ndarray]:
    """(total gates, ratio) axes of the candidate grid over the search box."""
    return (
        np.linspace(settings.n_gates_bounds[0], settings.n_gates_bounds[1],
                    settings.candidate_grid_size[0]),
        np.linspace(settings.ratio_bounds[0], settings.ratio_bounds[1],
                    settings.candidate_grid_size[1]),
    )


def _surface_features(points: np.ndarray,
                      n_gates_bounds: tuple[int, int],
                      ratio_bounds: tuple[float, float]) -> np.ndarray:
    """Cubic polynomial features of points scaled into the unit square."""
    points = np.asarray(points, dtype=float)
    lower = np.array([n_gates_bounds[0], ratio_bounds[0]], dtype=float)
    upper = np.array([n_gates_bounds[1], ratio_bounds[1]], dtype=float)
    scaled = np.clip((points - lower) / (upper - lower), 0.0, 1.0)
    gates = scaled[:, 0]
    ratio = scaled[:, 1]
    return np.column_stack([
        gates,
        ratio,
        gates**2,
        ratio**2,
        gates * ratio,
        gates**3,
        ratio**3,
        gates**2 * ratio,
        gates * ratio**2,
    ])


@dataclass
class MonotoneFidelitySurface:
    """
    Fitted surface for E[fidelity | total gates, two-qubit ratio].

    The logit of the fidelity is a cubic polynomial in the scaled coordinates
    with non-negative coefficients, so the expected fidelity cannot increase
    as either the gate count or the ratio increases.
    """
    alpha: float
    coefficients: np.ndarray
    n_gates_bounds: tuple[int, int]
    ratio_bounds: tuple[float, float]
    n_points: int

    def probability(self, points: np.ndarray) -> np.ndarray:
        features = _surface_features(points, self.n_gates_bounds, self.ratio_bounds)
        return expit(self.alpha - features @ self.coefficients)

    def probability_grid(
        self,
        grid_size: tuple[int, int],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n_gates_axis, n_ratios = grid_size
        gates = np.linspace(self.n_gates_bounds[0], self.n_gates_bounds[1], n_gates_axis)
        ratios = np.linspace(self.ratio_bounds[0], self.ratio_bounds[1], n_ratios)
        gates_grid, ratio_grid = np.meshgrid(gates, ratios)
        points = np.column_stack([gates_grid.ravel(), ratio_grid.ravel()])
        probabilities = self.probability(points).reshape(gates_grid.shape)
        return gates_grid, ratio_grid, probabilities

    def report(self, grid_size: tuple[int, int]) -> str:
        _, _, probabilities = self.probability_grid(grid_size)
        min_probability = float(np.min(probabilities))
        max_probability = float(np.max(probabilities))
        has_contour = min_probability <= 0.5 <= max_probability
        return (
            "Fitted monotone fidelity surface:\n"
            f"  fit points: {self.n_points}\n"
            f"  predicted fidelity range on grid: "
            f"{min_probability:.4f} to {max_probability:.4f}\n"
            f"  contains fidelity=0.5 contour: {has_contour}"
        )


def fit_monotone_surface_from_counts(
    points: np.ndarray,
    successes: np.ndarray,
    failures: np.ndarray,
    settings: CrossingSettings,
    n_points: int,
) -> MonotoneFidelitySurface:
    """Fit the monotone surface with the binomial likelihood for Boolean outcomes."""
    points = np.asarray(points, dtype=float)
    successes = np.asarray(successes, dtype=float)
    failures = np.asarray(failures, dtype=float)
    if np.any(successes < 0.0) or np.any(failures < 0.0):
        raise ValueError("Success and failure counts must be non-negative.")

    features = _surface_features(points, settings.n_gates_bounds, settings.ratio_bounds)

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

    x0 = np.zeros(features.shape[1] + 1, dtype=float)
    bounds = [(None, None)] + [(0.0, None)] * features.shape[1]
    result = minimize(fun=loss_and_grad, jac=True, x0=x0, bounds=bounds, method="L-BFGS-B")
    if not result.success:
        raise RuntimeError(result.message)

    return MonotoneFidelitySurface(
        alpha=float(result.x[0]),
        coefficients=result.x[1:],
        n_gates_bounds=settings.n_gates_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=n_points,
    )


def fit_monotone_fidelity_surface(data: RMBData,
                                  settings: CrossingSettings) -> MonotoneFidelitySurface:
    """Fit the monotone surface to the recorded Boolean outcomes in ``data``."""
    measured = measured_items(data)
    if len(measured) < settings.min_fit_points:
        raise ValueError(
            f"Need at least {settings.min_fit_points} data points to fit a surface, "
            f"got {len(measured)}."
        )

    successes = []
    failures = []
    for _, estimator in measured:
        counts = estimator.counts()
        successes.append(float(counts.get(True, 0)))
        failures.append(float(counts.get(False, 0)))

    return fit_monotone_surface_from_counts(
        points=config_points([config for config, _ in measured]),
        successes=np.asarray(successes, dtype=float),
        failures=np.asarray(failures, dtype=float),
        settings=settings,
        n_points=len(measured),
    )


def try_fit_monotone_fidelity_surface(
    data: RMBData,
    settings: CrossingSettings,
) -> MonotoneFidelitySurface | None:
    """Fit the monotone surface, or ``None`` when the data cannot support a fit."""
    try:
        return fit_monotone_fidelity_surface(data, settings)
    except (RuntimeError, ValueError):
        return None


def fit_monotone_surface_from_values(
    configs: list[RMBConfig],
    fidelities: np.ndarray,
    weights: np.ndarray,
    settings: CrossingSettings,
) -> MonotoneFidelitySurface:
    """
    Fit the monotone surface to fidelity values with per-config weights.

    The weighted Bernoulli likelihood equals the binomial likelihood with
    ``weights * fidelities`` successes and ``weights * (1 - fidelities)``
    failures, so this delegates to :func:`fit_monotone_surface_from_counts`.
    """
    fidelities = np.asarray(fidelities, dtype=float)
    weights = np.asarray(weights, dtype=float)
    return fit_monotone_surface_from_counts(
        points=config_points(configs),
        successes=weights * fidelities,
        failures=weights * (1.0 - fidelities),
        settings=settings,
        n_points=len(configs),
    )


@dataclass
class PhysicalDecaySurface:
    """
    Fitted randomized-benchmarking decay E[fidelity | one-qubit, two-qubit gates].

    Errors of each gate type compound multiplicatively, so the expected
    fidelity decays exponentially in the one- and two-qubit gate counts,

        p(N1, N2) = B + A * exp(-(gamma_1 * N1 + gamma_2 * N2)),

    between an amplitude ``A`` and a depolarized floor ``B``, with non-negative
    per-gate decay rates ``gamma_1``, ``gamma_2``. The fidelity = level contour
    is closed form, ``gamma_1 * N1 + gamma_2 * N2 = ln(A / (level - B))``, so in
    (total gates, ratio) coordinates the inverse crossing size is linear in the
    ratio.
    """
    gamma_1: float
    gamma_2: float
    amplitude: float
    floor: float
    n_gates_bounds: tuple[int, int]
    ratio_bounds: tuple[float, float]
    n_points: int

    def probability(self, points: np.ndarray) -> np.ndarray:
        """Expected fidelity at (total gates, two-qubit ratio) points."""
        points = np.asarray(points, dtype=float)
        n_gates = points[:, 0]
        ratio = points[:, 1]
        n_2 = ratio * n_gates
        n_1 = (1.0 - ratio) * n_gates
        return self.floor + self.amplitude * np.exp(
            -(self.gamma_1 * n_1 + self.gamma_2 * n_2))

    def crossing_n_gates(self, ratio: float, level: float = 0.5) -> float | None:
        """
        Total gate count where the fidelity crosses ``level`` at ``ratio``.

        Returns ``None`` when the level lies outside the surface's range or the
        decay rate at this ratio is non-positive, so the contour does not exist.
        """
        if not self.floor < level < self.floor + self.amplitude:
            return None
        denominator = self.gamma_1 + (self.gamma_2 - self.gamma_1) * ratio
        if denominator <= 0.0:
            return None
        return float(np.log(self.amplitude / (level - self.floor)) / denominator)

    def probability_grid(
        self,
        grid_size: tuple[int, int],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n_gates_axis, n_ratios = grid_size
        gates = np.linspace(self.n_gates_bounds[0], self.n_gates_bounds[1], n_gates_axis)
        ratios = np.linspace(self.ratio_bounds[0], self.ratio_bounds[1], n_ratios)
        gates_grid, ratio_grid = np.meshgrid(gates, ratios)
        points = np.column_stack([gates_grid.ravel(), ratio_grid.ravel()])
        probabilities = self.probability(points).reshape(gates_grid.shape)
        return gates_grid, ratio_grid, probabilities

    def params(self) -> np.ndarray:
        """``(gamma_1, gamma_2, amplitude, floor)`` as a float array."""
        return np.array([self.gamma_1, self.gamma_2, self.amplitude, self.floor],
                        dtype=float)

    def crossing_line(self, level: float = 0.5) -> tuple[float, float, float] | None:
        """
        Coefficients of the closed-form fidelity = ``level`` line.

        The contour is ``gamma_1 N1 + gamma_2 N2 = K`` with
        ``K = ln(A / (level - B))``, so in (total gates, ratio) coordinates the
        inverse crossing size is linear in the ratio,
        ``1 / n_gates(ratio) = a + b * ratio``.

        Returns
        -------
        tuple[float, float, float] | None
            ``(a, b, K)``, or ``None`` when ``level`` lies outside the
            surface's range so the line does not exist.
        """
        if not self.floor < level < self.floor + self.amplitude:
            return None
        constant = float(np.log(self.amplitude / (level - self.floor)))
        if constant <= 0.0:
            return None
        return self.gamma_1 / constant, (self.gamma_2 - self.gamma_1) / constant, constant

    def report(self) -> str:
        lines = [
            "Fitted physical decay surface:",
            f"  fit points: {self.n_points}",
            f"  fidelity(N1, N2) = {self.floor:.4f} + {self.amplitude:.4f} * "
            f"exp(-({self.gamma_1:.3e} * N1 + {self.gamma_2:.3e} * N2))",
            f"  fidelity(n_gates, ratio) = {self.floor:.4f} + {self.amplitude:.4f} * "
            f"exp(-n_gates * ({self.gamma_1:.3e} + {self.gamma_2 - self.gamma_1:.3e} * ratio))",
            f"  gamma_1 (one-qubit): {self.gamma_1:.3e}",
            f"  gamma_2 (two-qubit): {self.gamma_2:.3e}",
            f"  amplitude A: {self.amplitude:.4f}  floor B: {self.floor:.4f}",
        ]
        line = self.crossing_line()
        if line is not None:
            a, b, constant = line
            lines.append(
                f"  fidelity=0.5 line: n_gates(ratio) = 1 / "
                f"({self.gamma_1 / constant:.3e} + {(self.gamma_2 - self.gamma_1) / constant:.3e} * ratio)")
            lines.append(
                f"                     1/n_gates(ratio) = {a:.4e} + {b:.4e} * ratio")
        return "\n".join(lines)


def fit_physical_decay_from_counts(
    n_1_gates: np.ndarray,
    n_2_gates: np.ndarray,
    successes: np.ndarray,
    failures: np.ndarray,
    settings: CrossingSettings,
    n_points: int,
    *,
    n_restarts: int = 6,
) -> PhysicalDecaySurface:
    """
    Fit ``B + A exp(-(g1 N1 + g2 N2))`` by binomial maximum likelihood.

    The amplitude is reparametrized ``A = (1 - B) * c`` with ``c, B`` in
    ``[0, 1]`` so the predicted fidelity stays in ``[0, 1]`` under plain box
    constraints. Several gamma-scale restarts guard against local minima.
    """
    n_1 = np.asarray(n_1_gates, dtype=float)
    n_2 = np.asarray(n_2_gates, dtype=float)
    successes = np.asarray(successes, dtype=float)
    failures = np.asarray(failures, dtype=float)
    totals = successes + failures

    def loss_and_grad(theta: np.ndarray) -> tuple[float, np.ndarray]:
        gamma_1, gamma_2, c, floor = theta
        amplitude = (1.0 - floor) * c
        decay = np.exp(-(gamma_1 * n_1 + gamma_2 * n_2))
        eps = 1e-12
        p = np.clip(floor + amplitude * decay, eps, 1.0 - eps)
        loss = -np.sum(successes * np.log(p) + failures * np.log(1.0 - p))
        residual = (successes - totals * p) / (p * (1.0 - p))
        amplitude_decay = amplitude * decay
        grad_gamma_1 = float(np.sum(residual * amplitude_decay * n_1))
        grad_gamma_2 = float(np.sum(residual * amplitude_decay * n_2))
        grad_c = -float(np.sum(residual * (1.0 - floor) * decay))
        grad_floor = -float(np.sum(residual * (1.0 - c * decay)))
        return loss, np.array([grad_gamma_1, grad_gamma_2, grad_c, grad_floor])

    bounds = [(0.0, None), (0.0, None), (0.0, 1.0), (0.0, 1.0)]
    floor_guess = 1.0 / (2.0 ** settings.n_qubits)
    best = None
    for restart in range(max(1, n_restarts)):
        gamma_guess = 10.0 ** (-4.0 + restart)
        x0 = np.array([gamma_guess, gamma_guess, 0.9, floor_guess])
        result = minimize(fun=loss_and_grad, jac=True, x0=x0, bounds=bounds,
                          method="L-BFGS-B")
        if result.success and (best is None or result.fun < best.fun):
            best = result
    if best is None:
        raise RuntimeError("physical decay fit did not converge")

    gamma_1, gamma_2, c, floor = best.x
    return PhysicalDecaySurface(
        gamma_1=float(gamma_1),
        gamma_2=float(gamma_2),
        amplitude=float((1.0 - floor) * c),
        floor=float(floor),
        n_gates_bounds=settings.n_gates_bounds,
        ratio_bounds=settings.ratio_bounds,
        n_points=n_points,
    )


def fit_physical_decay(data: RMBData, settings: CrossingSettings) -> PhysicalDecaySurface:
    """Fit the physical decay surface to the recorded Boolean outcomes in ``data``."""
    measured = measured_items(data)
    if len(measured) < settings.min_fit_points:
        raise ValueError(
            f"Need at least {settings.min_fit_points} data points to fit a surface, "
            f"got {len(measured)}."
        )
    n_1_gates = [float(config.n_1qb_gates) for config, _ in measured]
    n_2_gates = [float(config.n_2qb_gates) for config, _ in measured]
    successes = []
    failures = []
    for _, estimator in measured:
        counts = estimator.counts()
        successes.append(float(counts.get(True, 0)))
        failures.append(float(counts.get(False, 0)))
    return fit_physical_decay_from_counts(
        np.asarray(n_1_gates), np.asarray(n_2_gates),
        np.asarray(successes), np.asarray(failures),
        settings, len(measured),
    )


def try_fit_physical_decay(data: RMBData,
                           settings: CrossingSettings) -> PhysicalDecaySurface | None:
    """Fit the physical decay surface, or ``None`` when the data cannot support a fit."""
    try:
        return fit_physical_decay(data, settings)
    except (RuntimeError, ValueError):
        return None


def super_level_fraction(surface: MonotoneFidelitySurface,
                         settings: CrossingSettings,
                         level: float = 0.5) -> float:
    """
    Fraction of the search box with fidelity above ``level``.

    The fraction is measured by area in the physical (one-qubit, two-qubit)
    gate-count plane: the (total gates, ratio) -> gate-count map has Jacobian
    equal to the total gate count, so the candidate-grid cells are weighted
    by it. Use with :func:`bootstrap_surfaces` for a credible interval.
    """
    gates_grid, _, probabilities = surface.probability_grid(settings.candidate_grid_size)
    return float(np.sum((probabilities > level) * gates_grid) / np.sum(gates_grid))


def contour_points_from_surface(
    surface: MonotoneFidelitySurface,
    settings: CrossingSettings,
    level: float = 0.5,
) -> np.ndarray:
    """Approximate contour points by linearly interpolating grid-edge crossings."""
    gates_grid, ratio_grid, probabilities = surface.probability_grid(settings.candidate_grid_size)
    delta = probabilities - level

    exact = np.zeros_like(delta, dtype=bool)
    exact[:, :-1] |= delta[:, :-1] == 0.0
    exact[:-1, :] |= delta[:-1, :] == 0.0
    points = [np.column_stack([gates_grid[exact], ratio_grid[exact]])]

    # Edges along the gates axis: interpolate the gates coordinate.
    p0, p1 = delta[:, :-1], delta[:, 1:]
    crossing = p0 * p1 < 0.0
    t = np.abs(p0[crossing]) / (np.abs(p0[crossing]) + np.abs(p1[crossing]))
    points.append(np.column_stack([
        (1.0 - t) * gates_grid[:, :-1][crossing] + t * gates_grid[:, 1:][crossing],
        ratio_grid[:, :-1][crossing],
    ]))

    # Edges along the ratio axis: interpolate the ratio coordinate.
    p0, p1 = delta[:-1, :], delta[1:, :]
    crossing = p0 * p1 < 0.0
    t = np.abs(p0[crossing]) / (np.abs(p0[crossing]) + np.abs(p1[crossing]))
    points.append(np.column_stack([
        gates_grid[:-1, :][crossing],
        (1.0 - t) * ratio_grid[:-1, :][crossing] + t * ratio_grid[1:, :][crossing],
    ]))

    stacked = np.concatenate(points)
    if len(stacked) == 0:
        return np.empty((0, 2), dtype=float)
    return np.unique(stacked, axis=0)


def bootstrap_surfaces(
    data: RMBData,
    settings: CrossingSettings,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
) -> list[MonotoneFidelitySurface]:
    """Fit monotone surfaces to posterior draws of the per-config fidelities."""
    rng = default_rng(seed)
    measured = measured_items(data)
    if len(measured) < settings.min_fit_points:
        return []

    configs = [config for config, _ in measured]
    posteriors = [estimator.posterior_alpha_beta() for _, estimator in measured]
    alpha_array = np.asarray([alpha for alpha, _ in posteriors], dtype=float)
    beta_array = np.asarray([beta for _, beta in posteriors], dtype=float)
    weights_array = np.asarray(
        [max(1, estimator.num_runs()) for _, estimator in measured], dtype=float)
    surfaces = []

    for _ in range(n_bootstrap):
        sampled_fidelities = rng.beta(alpha_array, beta_array)
        try:
            surfaces.append(fit_monotone_surface_from_values(
                configs,
                sampled_fidelities,
                weights_array,
                settings,
            ))
        except (RuntimeError, ValueError):
            continue

    return surfaces


def bootstrap_contours(
    data: RMBData,
    settings: CrossingSettings,
    *,
    n_bootstrap: int = 100,
    seed: int | None = None,
    surfaces: list[MonotoneFidelitySurface] | None = None,
) -> list[np.ndarray]:
    """
    Draw posterior monotone surfaces and return their p=0.5 contours.

    With ``surfaces``, the given surfaces are used instead of fitting new
    bootstrap draws, so a caller that needs both can bootstrap once.
    """
    if surfaces is None:
        surfaces = bootstrap_surfaces(data, settings, n_bootstrap=n_bootstrap, seed=seed)
    contours = []
    for surface in surfaces:
        contour = contour_points_from_surface(surface, settings)
        if len(contour) > 0:
            contours.append(contour)
    return contours


def save_crossings(rmb: RMB, settings: CrossingSettings, budget: Budget,
                   crossings: list[RMBConfig]) -> Path:
    """
    Save the RMB data to ``settings.save_path`` and the crossings to a sibling
    ``*_crossings.json``; return the resolved base path.
    """
    if settings.save_path is None:
        raise ValueError("settings.save_path is None; nothing to save.")
    rmb.save(settings.save_path)
    base_path = resolve_data_path(settings.save_path)
    payload = {
        "n_qubits": settings.n_qubits,
        "spent_hqc": budget.spent_hqc,
        "crossings": [
            {
                "n_gates": c.n_gates,
                "ratio_2qb_gates": round(c.ratio_2_qb_gates, 4),
                "n_1qb_gates": c.n_1qb_gates,
                "n_2qb_gates": c.n_2qb_gates,
                "use_scrambler": c.use_scrambler,
            }
            for c in crossings
        ],
    }
    crossings_path = base_path.parent / f"{base_path.stem}_crossings.json"
    crossings_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return base_path


def load_crossings(path: str | Path) -> list[RMBConfig]:
    """Load the crossing configs saved by :func:`save_crossings` next to ``path``."""
    base_path = resolve_data_path(path)
    crossings_path = base_path.parent / f"{base_path.stem}_crossings.json"
    payload = json.loads(crossings_path.read_text(encoding="utf-8"))
    return [
        RMBConfig(n_1qb_gates=record["n_1qb_gates"],
                  n_2qb_gates=record["n_2qb_gates"],
                  n_qubits=payload["n_qubits"],
                  use_scrambler=record.get("use_scrambler", True))
        for record in payload["crossings"]
    ]


def print_progress(settings: CrossingSettings, budget: Budget, message: str) -> None:
    """Print a progress line with the spent budget when ``settings.verbose``."""
    if settings.verbose:
        print(f"{message} (spent {budget.spent_hqc:.1f} / {settings.hqc_budget} HQC)")


def print_crossing(settings: CrossingSettings, budget: Budget,
                   label: str, config: RMBConfig) -> None:
    """Print a labelled crossing/anchor point when ``settings.verbose``."""
    print_progress(settings, budget,
                   f"{label}: ratio={config.ratio_2_qb_gates:.3f} n_gates={config.n_gates}")


def print_experiment_summary(data: RMBData, settings: CrossingSettings, budget: Budget,
                             stop_reason: str | None = None) -> None:
    """Print measured-data and spent-budget statistics when ``settings.verbose``."""
    if not settings.verbose:
        return
    measured_configs = [config for config, _ in measured_items(data)]
    repeats = [data[config].num_runs() for config in measured_configs]

    print("\nExperiment summary")
    if stop_reason is not None:
        print(f"  stop reason: {stop_reason}")
    print(f"  stitched submissions: {budget.jobs}")
    print(f"  max circuits in one submission: {budget.max_job_circuits}")
    print(f"  distinct configs measured: {len(measured_configs)}")
    print(f"  total circuits: {sum(repeats)}")
    print(f"  max repeats of one config: {max(repeats, default=0)}")
    print(f"  HQC spent: {budget.spent_hqc:.1f} / {settings.hqc_budget}")
    if measured_configs:
        n_gates = [config.n_gates for config in measured_configs]
        ratios = [config.ratio_2_qb_gates for config in measured_configs]
        print(f"  sampled gates range: {min(n_gates)} to {max(n_gates)}")
        print(f"  sampled two-qubit ratio range: {min(ratios):.2f} to {max(ratios):.2f}")


def print_fit_reports(data: RMBData, settings: CrossingSettings) -> None:
    """Print a monotone surface-fit report per ``n_qubits`` when ``settings.verbose``."""
    if not settings.verbose:
        return
    for n_qubits, group in sorted(grouped_by_n_qubits(data).items()):
        surface = try_fit_monotone_fidelity_surface(group, settings)
        if surface is None:
            print(f"n_qubits={n_qubits}: not enough data for a monotone surface fit yet.")
            continue
        print(f"\nn_qubits={n_qubits}")
        print(surface.report(settings.candidate_grid_size))
