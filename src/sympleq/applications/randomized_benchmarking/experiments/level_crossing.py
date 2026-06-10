"""
Trace the fidelity = 0.5 line in (total gates, two-qubit gate ratio) space at
fixed n_qubits.

The experiment finds a first crossing by bisecting in the total gate count at
a small fixed two-qubit ratio, then tracks the line by stepping the ratio up
and re-bisecting inside a bracket predicted from the previous crossings.
Configs are built from the transformed parameters on demand. Every probe is
metered in Quantinuum credits (HQC) via ``pytket_bare_simulation_cost`` with
the base submission cost paid once per stitched batch, and the whole run stops
when the HQC budget is exhausted.

Circuits are only ever simulated with the SympleQ emulation of Quantinuum
hardware (``QuantinuumBackend.default_sympleq_backend``); nothing is submitted
to Quantinuum.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.random import default_rng
from scipy.stats import beta

from sympleq.applications.randomized_benchmarking.RMB import RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.backends.quantinuum import QuantinuumBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.integrations.quantinuum.utils import (
    BASE_SIMULATION_COST,
    NATIVE_GATES_SET,
    pytket_bare_simulation_cost,
    to_pytket_circuit,
)


@dataclass(frozen=True)
class LevelCrossingSettings:
    """
    Settings for a fidelity = 0.5 level-crossing experiment.

    Parameters
    ----------
    n_qubits : int
        Fixed number of qubits for the whole run.
    random_elimination : float
        ``RMBConfig.random_elimination`` used for every probed config.
    n_gates_bounds : tuple[int, int]
        Search range for the total-gates crossing at each ratio value.
    ratio_start : float
        Two-qubit gate ratio of the first (cheap) crossing search.
    ratio_step : float
        Increment of the ratio between consecutive traced crossings.
    ratio_max : float
        Largest ratio to attempt.
    n_gates_resolution : float
        Bisection stops once the bracket is narrower than this fraction of
        its upper edge.
    max_shots_per_probe : int
        Hard cap on circuits spent on a single config probe.
    decision_confidence : float
        Posterior mass required on one side of 0.5 to end a probe early.
        Values at or below 0.875 let one or two unanimous shots decide a
        side, which makes bracket checks unreliable; keep it above that.
    max_cost_per_run : float
        HQC cap of one stitched submission, mirroring
        ``QuantinuumBackend.max_cost_per_run``. It sets how many probe
        circuits share one base submission cost.
    hqc_budget : float
        Total Quantinuum credits the run may spend.
    rng_seed : int | None
        Seed for the run RNG. ``None`` uses fresh entropy.
    save_path : str | Path | None
        Where to save the RMB data (see :meth:`RMB.save`). ``None`` skips
        saving. Crossings go to a sibling ``*_crossings.json`` file and the
        plot to a sibling ``*.png``.
    plot : bool
        Show the data and the traced line at the end of the run.
    verbose : bool
        Print progress while running.
    """
    n_qubits: int = 4
    random_elimination: float = 0.1
    n_gates_bounds: tuple[int, int] = (10, 30000)
    ratio_start: float = 0.0
    ratio_step: float = 0.05
    ratio_max: float = 1.0
    n_gates_resolution: float = 0.1
    max_shots_per_probe: int = 16
    decision_confidence: float = 0.876
    max_cost_per_run: float = 10.0
    hqc_budget: float = 2000.0
    rng_seed: int | None = 1234
    save_path: str | Path | None = "level_crossing.json"
    plot: bool = True
    verbose: bool = True


@dataclass
class Budget:
    """Remaining and spent Quantinuum credits."""
    remaining_hqc: float
    spent_hqc: float = 0.0

    def spend(self, cost: float) -> None:
        self.remaining_hqc -= cost
        self.spent_hqc += cost

    def can_afford(self, cost: float) -> bool:
        return self.remaining_hqc >= cost


def make_config(settings: LevelCrossingSettings, n_gates: int, ratio: float) -> RMBConfig:
    """
    Build a valid native-gate RMBConfig from total gates and two-qubit ratio.

    Gate counts are rounded to the even values the circuit construction
    actually realizes, with at least the 2 * n_qubits scrambler one-qubit
    gates, so the config matches the generated circuits and the Quantinuum
    cost is computed for what would really run.
    """
    n_2qb_gates = max(0, 2 * round(ratio * n_gates / 2))
    n_1qb_gates = max(2 * settings.n_qubits, 2 * round((n_gates - n_2qb_gates) / 2))
    return (
        RMBConfig.default()
        .with_n_qubits(settings.n_qubits)
        .with_n_1qb_gates(n_1qb_gates)
        .with_n_2qb_gates(n_2qb_gates)
        .with_random_elimination(settings.random_elimination)
        .with_gates_set(tuple(NATIVE_GATES_SET))
    )


def single_circuit_bare_hqc(config: RMBConfig) -> float:
    """
    Bare Quantinuum credits for one circuit of this config at one shot,
    excluding the per-submission base cost but including the qubit resets
    needed to stitch it after another circuit.

    The cost depends only on the gate counts, so any circuit drawn from the
    config prices all of them.
    """
    circuit = to_pytket_circuit(config.random_circuit(rng=default_rng(0)))
    return pytket_bare_simulation_cost(circuit) + config.n_qubits / 5000


def stitched_batch_hqc(bare_hqc: float, n_circuits: int) -> float:
    """
    Quantinuum credits for one stitched submission of ``n_circuits`` circuits.

    Mirrors ``QuantinuumBackend.fidelity_estimation``, which stitches circuits
    into a single program, so the base submission cost is paid once per batch.
    """
    return BASE_SIMULATION_COST + n_circuits * bare_hqc


def stitch_batch_size(bare_hqc: float, settings: LevelCrossingSettings) -> int:
    """Number of circuits at this bare cost that fit in one stitched submission."""
    if bare_hqc <= 0:
        return settings.max_shots_per_probe
    affordable = int((settings.max_cost_per_run - BASE_SIMULATION_COST) // bare_hqc)
    return max(1, min(affordable, settings.max_shots_per_probe))


def posterior_above(estimator: BayesianEstimator) -> float:
    """Posterior probability that the fidelity of ``estimator`` is above 0.5."""
    counts = estimator.counts()
    return float(beta.sf(0.5, counts.get(True, 0) + 1, counts.get(False, 0) + 1))


def posterior_mean(estimator: BayesianEstimator) -> float:
    """Beta(1, 1) posterior mean fidelity of ``estimator``."""
    counts = estimator.counts()
    return (counts.get(True, 0) + 1.0) / (estimator.num_runs() + 2.0)


def implied_above(rmb: RMB, config: RMBConfig, settings: LevelCrossingSettings) -> float | None:
    """
    Side of 0.5 implied for ``config`` by monotonicity of already-decided data.

    Fidelity decreases when gates of either kind are added, so a config with
    at least as many one- and two-qubit gates as a confidently-below config
    is below, and one with at most as many of each as a confidently-above
    config is above. Such configs need no measurements at all.

    Returns
    -------
    float | None
        1.0 (above), 0.0 (below), or ``None`` when existing data does not
        determine the side.
    """
    for other, estimator in rmb._data.items():
        if other.n_qubits != config.n_qubits:
            continue
        above = posterior_above(estimator)
        if (above >= settings.decision_confidence
                and config.n_1qb_gates <= other.n_1qb_gates
                and config.n_2qb_gates <= other.n_2qb_gates):
            return 1.0
        if (1.0 - above >= settings.decision_confidence
                and config.n_1qb_gates >= other.n_1qb_gates
                and config.n_2qb_gates >= other.n_2qb_gates):
            return 0.0
    return None


def probe_fidelity(rmb: RMB, config: RMBConfig, budget: Budget,
                   settings: LevelCrossingSettings) -> tuple[float | None, float]:
    """
    Estimate the fidelity of ``config`` with as few circuits as possible.

    Records circuit outcomes into ``rmb._data``, in stitched batches priced
    like the Quantinuum backend, until the Beta posterior places
    ``decision_confidence`` mass on one side of 0.5, the shot cap is reached,
    or the budget runs out. Configs whose side is already implied by
    monotonicity from existing data are not measured at all.

    Returns
    -------
    tuple[float | None, float]
        Posterior mean fidelity and posterior probability that the fidelity
        is above 0.5. The mean is ``None`` when the budget stopped the probe
        before it could either decide a side or reach the shot cap.
    """
    estimator = rmb._data.setdefault(config, BayesianEstimator(threshold=0.0, min_runs=0))
    above = posterior_above(estimator)
    if max(above, 1.0 - above) < settings.decision_confidence:
        implied = implied_above(rmb, config, settings)
        if implied is not None:
            return posterior_mean(estimator), implied

    bare_hqc = single_circuit_bare_hqc(config)

    while estimator.num_runs() < settings.max_shots_per_probe:
        above = posterior_above(estimator)
        if max(above, 1.0 - above) >= settings.decision_confidence:
            break
        n_circuits = min(stitch_batch_size(bare_hqc, settings),
                         settings.max_shots_per_probe - estimator.num_runs())
        while n_circuits > 0 and not budget.can_afford(stitched_batch_hqc(bare_hqc, n_circuits)):
            n_circuits -= 1
        if n_circuits <= 0:
            break
        for _ in range(n_circuits):
            for outcome in rmb.backend.fidelity_estimation(config, rmb.rng):
                estimator.record(bool(outcome))
        budget.spend(stitched_batch_hqc(bare_hqc, n_circuits))

    above = posterior_above(estimator)
    decided = max(above, 1.0 - above) >= settings.decision_confidence
    if not decided and estimator.num_runs() < settings.max_shots_per_probe:
        return None, above
    return posterior_mean(estimator), above


def predict_crossing(crossings: list[RMBConfig], ratio: float) -> int | None:
    """
    Predict the crossing total gates at ``ratio`` from previous crossings.

    At the crossing the accumulated noise is roughly constant, so the inverse
    crossing size is approximately linear in the ratio: 1/n*(r) = a + b*r.

    Returns
    -------
    int | None
        Predicted total gates, or ``None`` with fewer than two crossings or
        when the fit does not cross at this ratio.
    """
    if len(crossings) < 2:
        return None
    ratios = [c.ratio_2_qb_gates for c in crossings]
    inverse_sizes = [1.0 / c.n_gates for c in crossings]
    slope, intercept = np.polyfit(ratios, inverse_sizes, 1)
    inverse = intercept + slope * ratio
    if inverse <= 0:
        return None
    return round(1.0 / inverse)


def find_crossing(rmb: RMB, budget: Budget, settings: LevelCrossingSettings,
                  ratio: float, lo: int, hi: int) -> RMBConfig | None:
    """
    Bisect in total gates for the fidelity = 0.5 crossing at fixed ratio.

    Fidelity decreases monotonically with gate count, so the bracket needs
    fidelity confidently above 0.5 at ``lo`` and confidently below at ``hi``.
    The bisection only branches on confident side decisions; any point whose
    shot cap cannot tell it from 0.5, endpoints included, is statistically on
    the line and is returned as the crossing directly.

    Returns
    -------
    RMBConfig | None
        The crossing config, or ``None`` when the line lies confidently
        outside the bracket or the budget runs out.
    """
    config_lo = make_config(settings, lo, ratio)
    p_lo, above_lo = probe_fidelity(rmb, config_lo, budget, settings)
    if p_lo is None:
        return None
    if above_lo < settings.decision_confidence:
        if 1.0 - above_lo < settings.decision_confidence:
            return config_lo
        return None
    config_hi = make_config(settings, hi, ratio)
    p_hi, above_hi = probe_fidelity(rmb, config_hi, budget, settings)
    if p_hi is None:
        return None
    if 1.0 - above_hi < settings.decision_confidence:
        if above_hi < settings.decision_confidence:
            return config_hi
        return None

    while hi - lo > max(2, settings.n_gates_resolution * hi):
        mid = 2 * round((lo + hi) / 4)
        if mid in (lo, hi):
            break
        p_mid, above_mid = probe_fidelity(rmb, make_config(settings, mid, ratio), budget, settings)
        if p_mid is None:
            break
        if settings.verbose:
            print(f"  bisect ratio={ratio:.3f}: n_gates={mid} p={p_mid:.2f} "
                  f"(spent {budget.spent_hqc:.1f} HQC)")
        if above_mid >= settings.decision_confidence:
            lo = mid
        elif 1.0 - above_mid >= settings.decision_confidence:
            hi = mid
        else:
            return make_config(settings, mid, ratio)

    return make_config(settings, (lo + hi) // 2, ratio)


def trace_level_line(settings: LevelCrossingSettings) -> tuple[RMB, list[RMBConfig]]:
    """
    Find the fidelity = 0.5 line and trace it up in the two-qubit gate ratio.

    Returns
    -------
    tuple[RMB, list[RMBConfig]]
        The RMB holding all recorded data and the crossing configs in
        increasing ratio order.
    """
    rng = default_rng(settings.rng_seed)
    rmb = RMB.default(rng).with_backend(QuantinuumBackend.default_sympleq_backend())
    budget = Budget(remaining_hqc=settings.hqc_budget)
    crossings: list[RMBConfig] = []

    n_gates_min, n_gates_max = settings.n_gates_bounds
    lo, hi = n_gates_min, n_gates_max
    step_index = 0
    ratio = settings.ratio_start

    while ratio <= settings.ratio_max + 1e-9 and budget.remaining_hqc > 0:
        crossing = find_crossing(rmb, budget, settings, ratio, lo, hi)
        if crossing is None and (lo, hi) != (n_gates_min, n_gates_max):
            # The predicted bracket missed the line; retry with full bounds.
            lo, hi = n_gates_min, n_gates_max
            crossing = find_crossing(rmb, budget, settings, ratio, lo, hi)

        if crossing is not None:
            crossings.append(crossing)
            if settings.verbose:
                print(f"Crossing: ratio={crossing.ratio_2_qb_gates:.3f} n_gates={crossing.n_gates} "
                      f"(spent {budget.spent_hqc:.1f} / {settings.hqc_budget} HQC)")
        elif settings.verbose:
            # A crossing in total gates exists at every ratio, so a miss means
            # either the budget ran dry or unlucky endpoint reads; the next
            # ratio is still worth trying with the remaining budget.
            reason = "budget exhausted" if budget.remaining_hqc <= 0 else "endpoints unresolved"
            print(f"No crossing found at ratio={ratio:.3f} ({reason}).")

        step_index += 1
        ratio = settings.ratio_start + step_index * settings.ratio_step

        # Bracket the next crossing around the inverse-linear fit prediction;
        # before the fit is possible, search below the previous crossing.
        predicted = predict_crossing(crossings, ratio)
        if predicted is not None:
            lo = max(n_gates_min, predicted // 2)
            hi = min(n_gates_max, 2 * predicted)
        elif crossings:
            lo, hi = n_gates_min, crossings[-1].n_gates
        else:
            lo, hi = n_gates_min, n_gates_max
        if hi <= lo:
            lo, hi = n_gates_min, n_gates_max

    if settings.verbose:
        print(f"\nTraced {len(crossings)} crossings, "
              f"spent {budget.spent_hqc:.1f} / {settings.hqc_budget} HQC.")
        for config in crossings:
            print(f"  ratio={config.ratio_2_qb_gates:.3f} n_gates={config.n_gates:>6}")

    if settings.save_path is not None:
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
                }
                for c in crossings
            ],
        }
        crossings_path = base_path.parent / f"{base_path.stem}_crossings.json"
        crossings_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if settings.plot:
        import matplotlib.pyplot as plt
        axes = RMB.plot_data(rmb._data, show=False, skip_incomplete=False, level_line=crossings)
        if axes and settings.save_path is not None:
            axes[0].figure.savefig(base_path.parent / f"{base_path.stem}.png",
                                   dpi=150, bbox_inches="tight")
        plt.show()

    return rmb, crossings


if __name__ == "__main__":
    rmb, crossings = trace_level_line(LevelCrossingSettings())
