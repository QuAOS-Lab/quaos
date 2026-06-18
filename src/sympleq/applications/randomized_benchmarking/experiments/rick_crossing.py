"""
GP level-set estimation of the fidelity = 0.5 line in (total gates,
two-qubit gate ratio) space at fixed n_qubits, driven by aepsych.

Stripped-down port of ``LSE_manual_integration.ipynb`` (rick_benchmark
branch). A Gaussian-process classification model (aepsych/BoTorch) is fitted
to the Boolean fidelity outcomes; an acquisition function (``GlobalSUR`` or
``MCLevelSetEstimation``) proposes the next configs to measure, batched
greedily into stitched submissions, until the HQC budget is exhausted. The
GP is seeded with one assumed success at the easy corner and one assumed
failure at the hard corner of the search box.

Measurements are priced with the same HQC model as ``level_crossing.py``
(one base submission cost per stitched batch plus the pytket bare cost of
every circuit; the notebook's hand-rolled weights were replaced by the
shared model), and outcomes are drawn with the same common-random-numbers
scheme, so runs are directly comparable with the other experiments.

Ported deviations from the notebook: ``scrambling_probability`` no longer
exists on this branch and was dropped; gate counts are rounded to the even
values the circuit construction realizes; ``LSE_budget``/``batch_budget``
map to the shared ``hqc_budget``/``max_cost_per_run``; the notebook's raw
(one-qubit gates, two-qubit gates) search box was replaced by the (total
gates, two-qubit gate ratio) box all experiments share.
"""
from __future__ import annotations

import logging
import warnings
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
from aepsych.config import Config
from aepsych.strategy import SequentialStrategy

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


@dataclass(frozen=True)
class RickCrossingSettings(CrossingSettings):
    """
    Settings for a GP level-set estimation crossing experiment.

    The problem definition, including the search box, lives in
    :class:`~.common.CrossingSettings`. The fields here are the aepsych
    knobs with the notebook's run profile values.
    """
    save_path: str | Path | None = "rick_crossing.json"

    target_threshold: float = 0.5

    # aepsych strategy. The acquisition optimizer dominates the runtime at
    # ~restarts * samples acquisition evaluations per suggested point; these
    # are the notebook profile values, and cutting them (e.g. 2 restarts,
    # 2000 samples) speeds suggestions up ~15x for little gain on a smooth
    # two-parameter box.
    initial_random_samples: int = 10
    optimization_steps: int = 10000
    inducing_size: int = 150
    acquisition_function: str = "GlobalSUR"
    acquisition_restarts: int = 4
    acquisition_samples: int = 6000

    # Stitched batch planning.
    batching: bool = True
    max_batch_size: int = 80


def point_from_config(config: RMBConfig) -> torch.Tensor:
    """aepsych observation point of a config: its (total gates, ratio) coordinates."""
    return torch.tensor(config_points([config]), dtype=torch.double)


def one_shot_requests(configs: list[RMBConfig]) -> list[MeasurementRequest]:
    return [MeasurementRequest(config, 1) for config in configs]


def contour_target(settings: RickCrossingSettings) -> float:
    target = float(settings.target_threshold)
    if not 0.0 < target < 1.0:
        raise ValueError(f"target_threshold must be within (0, 1), got {target}.")
    return target


def aepsych_config_string(settings: RickCrossingSettings) -> str:
    target = contour_target(settings)
    return f"""
    [common]
    parnames = [n_gates, ratio]
    outcome_types = [binary]
    strategy_names = [init_strat, opt_strat]

    [n_gates]
    par_type = continuous
    lower_bound = {settings.n_gates_bounds[0]}
    upper_bound = {settings.n_gates_bounds[1]}

    [ratio]
    par_type = continuous
    lower_bound = {settings.ratio_bounds[0]}
    upper_bound = {settings.ratio_bounds[1]}

    [init_strat]
    min_asks = {settings.initial_random_samples}
    generator = SobolGenerator

    [opt_strat]
    min_asks = {settings.optimization_steps}
    model = GPClassificationModel
    generator = OptimizeAcqfGenerator

    [GPClassificationModel]
    inducing_size = {settings.inducing_size}
    likelihood = BernoulliLikelihood

    [OptimizeAcqfGenerator]
    acqf = {settings.acquisition_function}
    restarts = {settings.acquisition_restarts}
    samps = {settings.acquisition_samples}

    [{settings.acquisition_function}]
    target = {target}
    """.strip()


def build_strategy(settings: RickCrossingSettings) -> SequentialStrategy:
    config = Config()
    config.update(config_str=aepsych_config_string(settings))
    return SequentialStrategy.from_config(config)


def batch_candidates(
    strategy: SequentialStrategy,
    settings: RickCrossingSettings,
) -> Iterator[RMBConfig]:
    """
    Candidate configs for one stitched batch, in acquisition order.

    ``MCLevelSetEstimation`` scores a whole batch in one call; the other
    acquisition functions are asked one point at a time with the batch
    selected so far pending.
    """
    if settings.acquisition_function == "MCLevelSetEstimation":
        x_batch = strategy.gen(num_points=settings.max_batch_size)
        for i in range(settings.max_batch_size):
            yield settings.make_config(
                float(x_batch[i, 0].item()), float(x_batch[i, 1].item()))
        return
    pending_points: list[torch.Tensor] = []
    while True:
        x_pending = torch.cat(pending_points, dim=0) if pending_points else None
        x = strategy.gen(X_pending=x_pending)
        candidate = settings.make_config(float(x[0, 0].item()), float(x[0, 1].item()))
        pending_points.append(point_from_config(candidate))
        yield candidate


def select_batch(
    strategy: SequentialStrategy,
    settings: RickCrossingSettings,
    budget: Budget,
) -> tuple[list[RMBConfig], bool]:
    """
    Ask the acquisition function for the next stitched batch of configs.

    Returns the selected configs and whether the run budget was found to be
    exhausted during selection. Batches are grown greedily until the next
    candidate would exceed the per-submission cost cap or the budget.
    """
    if not settings.batching:
        x = strategy.gen()
        return [settings.make_config(
            float(x[0, 0].item()), float(x[0, 1].item()))], False

    selected: list[RMBConfig] = []
    for candidate in batch_candidates(strategy, settings):
        cost = batch_hqc_cost(one_shot_requests(selected + [candidate]))
        if cost > settings.max_cost_per_run:
            break
        if cost > budget.remaining_hqc:
            return selected, True
        selected.append(candidate)
    return selected, False


def run(settings: RickCrossingSettings) -> tuple[RMB, list[RMBConfig]]:
    """
    Run the GP level-set estimation experiment.

    Returns
    -------
    tuple[RMB, list[RMBConfig]]
        The RMB holding all recorded data and the configs on the fitted
        GP's fidelity = target contour, in increasing ratio order.
    """
    logging.getLogger().setLevel(logging.WARNING)
    warnings.filterwarnings("ignore")
    torch.set_default_dtype(torch.float64)
    if settings.rng_seed is not None:
        # The GP fit and the acquisition optimization draw from torch's
        # global RNG; seeding it makes the suggested configs reproducible.
        torch.manual_seed(settings.rng_seed)

    rng, rmb, budget = start_run(settings)
    data = rmb._data
    strategy = build_strategy(settings)
    results: list[tuple[float, float, int]] = []

    # Seed the GP with an assumed success at the easy corner and an assumed
    # failure at the hard corner; neither is measured nor charged.
    corners = [
        (float(settings.n_gates_bounds[0]), float(settings.ratio_bounds[0]), 1),
        (float(settings.n_gates_bounds[1]), float(settings.ratio_bounds[1]), 0),
    ]
    for n_gates, ratio, outcome in corners:
        strategy.add_data(torch.tensor([[n_gates, ratio]], dtype=torch.double), [outcome])
        corner = settings.make_config(n_gates, ratio)
        results.append((float(corner.n_1qb_gates), float(corner.n_2qb_gates), outcome))

    exhausted = False
    while not exhausted and budget.remaining_hqc > 0:
        selected, exhausted = select_batch(strategy, settings, budget)
        if not selected:
            # Nothing affordable fits in a submission; stop instead of asking
            # the acquisition function for the same point forever.
            break
        outcomes_by_config = spend_request_batch(
            rmb.backend, rng, data, one_shot_requests(selected), seed=settings.rng_seed)
        for config, outcomes in outcomes_by_config.items():
            for outcome in outcomes:
                strategy.add_data(point_from_config(config), [int(outcome)])
                results.append((float(config.n_1qb_gates), float(config.n_2qb_gates),
                                int(outcome)))
        budget.spend_batch(batch_hqc_cost(one_shot_requests(selected)), len(selected))
        print_progress(settings, budget, f"Measured batch of {len(selected)} configs")

    stop_reason = ("HQC budget exhausted" if exhausted or budget.remaining_hqc <= 0
                   else "no affordable batch left")
    print_experiment_summary(data, settings, budget, stop_reason=stop_reason)
    print_fit_reports(data, settings)

    crossings = level_set_configs(strategy, settings)

    base_path = None
    if settings.save_path is not None:
        base_path = save_crossings(rmb, settings, budget, crossings)

    if settings.plot:
        import matplotlib.pyplot as plt
        from sympleq.applications.randomized_benchmarking.experiments.plots import (
            plot_crossing_results,
            plot_gp_level_set,
        )
        plot_crossing_results(data, settings, crossings, base_path=base_path, show=False)
        if strategy.model is not None:
            gp_png_path = None
            if base_path is not None:
                gp_png_path = base_path.parent / f"{base_path.stem}_gp.png"
            one_q_bounds, two_q_bounds = measured_gate_count_bounds(data)
            probabilities, latent_mean, latent_variance = predict_level_set(
                strategy, settings, one_q_bounds=one_q_bounds, two_q_bounds=two_q_bounds)
            plot_gp_level_set(probabilities, latent_mean, latent_variance, results,
                              one_q_bounds=one_q_bounds,
                              two_q_bounds=two_q_bounds,
                              target=contour_target(settings),
                              png_path=gp_png_path, show=False)
        plt.show()

    return rmb, crossings


def predict_level_set(
    strategy: SequentialStrategy,
    settings: RickCrossingSettings,
    *,
    one_q_bounds: tuple[int, int],
    two_q_bounds: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluate the fitted GP on the level-set grid of a gate-count box.

    The (one-qubit, two-qubit) grid points are mapped into the (total gates,
    two-qubit gate ratio) coordinates the GP is fitted in, clipped to the
    search box. Returns the success probability and the latent mean and
    variance as plain grid-shaped arrays, ready for
    :func:`~.plots.plot_gp_level_set`.
    """
    from sympleq.applications.randomized_benchmarking.experiments.plots import (
        gate_plane_points,
        level_set_grid,
    )

    one_q_grid, two_q_grid = level_set_grid(one_q_bounds, two_q_bounds)
    points = gate_plane_points(one_q_grid, two_q_grid)
    grid = torch.tensor(np.column_stack([
        np.clip(points[:, 0], settings.n_gates_bounds[0], settings.n_gates_bounds[1]),
        np.clip(points[:, 1], settings.ratio_bounds[0], settings.ratio_bounds[1]),
    ]), dtype=torch.double)
    with torch.no_grad():
        probabilities, _ = strategy.model.predict(grid, probability_space=True)
        latent_mean, latent_variance = strategy.model.predict(grid)
    return (probabilities.numpy().reshape(one_q_grid.shape),
            latent_mean.numpy().reshape(one_q_grid.shape),
            latent_variance.numpy().reshape(one_q_grid.shape))


def level_set_configs(strategy: SequentialStrategy,
                      settings: RickCrossingSettings) -> list[RMBConfig]:
    """
    Configs on the GP's fidelity = target contour, in increasing ratio order.

    At each ratio of the candidate grid, the total gate count where the
    predicted success probability first crosses the target is linearly
    interpolated; ratios the GP keeps entirely on one side of the target
    contribute no config. Returns no configs when the budget died before
    the strategy could fit a model.
    """
    if strategy.model is None:
        return []
    gates_axis, ratio_axis = candidate_axes(settings)
    gates_grid, ratio_grid = np.meshgrid(gates_axis, ratio_axis)
    grid = torch.tensor(np.stack([gates_grid.ravel(), ratio_grid.ravel()], axis=1),
                        dtype=torch.double)
    with torch.no_grad():
        probabilities, _ = strategy.model.predict(grid, probability_space=True)
    delta = probabilities.numpy().reshape(gates_grid.shape) - contour_target(settings)

    configs: list[RMBConfig] = []
    seen: set[RMBConfig] = set()
    for row in range(delta.shape[0]):
        crossing = np.where(delta[row, :-1] * delta[row, 1:] < 0)[0]
        if len(crossing) == 0:
            continue
        i = int(crossing[0])
        t = abs(delta[row, i]) / (abs(delta[row, i]) + abs(delta[row, i + 1]))
        config = settings.make_config(
            float((1.0 - t) * gates_axis[i] + t * gates_axis[i + 1]),
            float(ratio_axis[row]))
        if config not in seen:
            seen.add(config)
            configs.append(config)
    # Gate-count rounding can locally reorder the realized ratios.
    return sorted(configs, key=lambda c: (c.ratio_2_qb_gates, c.n_gates))


if __name__ == "__main__":
    rmb, crossings = run(RickCrossingSettings())
