"""
Cost-aware global active design for the 3-D fidelity = 0.5 boundary surface.

This module estimates the randomized-benchmarking survival-fidelity 0.5 boundary
jointly over

    (r, Q) = (two-qubit gate ratio, number of qubits),

i.e. the total gate count n*(r, Q) at which the survival probability crosses
0.5, as one smooth low-dimensional *surface* rather than one 1-D contour per Q.

How this differs from bootstrap_batched_monotone_tracing.py
-----------------------------------------------------------
1.  Global parametric surface.  The inverse boundary is the two-rate RB
    exponential, parameterised by the per-gate log-error rates themselves:

        1 / n*(r, Q) = [ L1(Q) (1 - r) + L2(Q) * r ] / ln 2,
        L1(Q) = L10 + m1 * (Q - Qref),   L2(Q) = L20 + m2 * (Q - Qref),

    where L1, L2 are the one- and two-qubit per-gate rates (so the boundary
    slope (L2 - L1)/ln 2 correctly carries the -L1 term), and m1, m2 let each
    rate vary linearly with register size Q.  The survival probability uses the
    physically-motivated randomized-benchmarking exponential link

        p(n, r, Q) = V * (1 - B(Q)) * 2^{-n / n*(r,Q)} + B(Q),

    where B(Q) is the survival asymptote (default 2^{-Q}, the depolarizing
    value) and V is the global visibility (SPAM amplitude, V = 1 with no SPAM).
    The transition sharpness is the physical constant ln 2 in the exponent, not
    a free parameter; the boundary n* = 1/D sits at n/n* = 1 (renormalised
    fidelity 0.5) independently of V and B.  Pooling every measurement across all
    r and Q lets the whole surface be inferred from far fewer circuits than
    tracing each Q independently.  Expanded in (r, Q) the model is bilinear, so
    when m1 = m2 = 0 (rates independent of Q) it reduces to the single-register
    inverse boundary 1/(q + s r) of the original method.

2.  Deterministic grid posterior.  The 5-parameter posterior (L1, L2,
    V) is represented on a moderate grid that is re-centred on the running
    posterior each update.  By default the visibility axis is pinned at V = 1
    (resolution 1), leaving the four rate parameters free; free V to fit SPAM on
    hardware.  There is no importance sampling from a fixed prior, so the
    effective sample size cannot collapse as data accumulate, and the fit /
    acquisition are reproducible with no per-step Monte-Carlo jitter.

3.  No top-down trace.  Each iteration assembles ONE stitched submission
    (cost <= max_cost_per_run) of (n, r, Q) probes chosen globally by expected
    information gain per HQC.  Probes AND per-probe shot counts are packed
    greedily to fill the submission, amortising the per-submission base cost.
    The acquisition score is BALD -- the mutual information between the next
    Bernoulli outcome and the surface parameters -- divided by the exact
    stitched HQC cost from batch_hqc_cost.  Because the scored surface is a
    deterministic function of the parameters, reducing parameter uncertainty
    reduces surface uncertainty, and BALD concentrates probes where the model
    is least certain about outcomes (i.e. near the boundary across all r, Q).

4.  Cheap misspecification diagnostics.  Standardised residuals of every
    measured configuration are reported each update (zero extra cost); an
    optional deep-probe check can be enabled to test the link in the tails.

Scoring target
--------------
The benchmark score is the log-ratio surface error

    exp( -sqrt( mean_{r,Q} [ log( n_hat(r,Q) / n_truth(r,Q) ) ]^2 ) ),

so the acquisition minimises posterior variance of log n*(r,Q) integrated over
the scored (r, Q) region.  If a ground-truth surface callable is supplied via
``settings.truth_boundary`` the score is reported; otherwise the integrated
posterior log-boundary uncertainty is reported as the stopping diagnostic.
"""
from __future__ import annotations

import json
from collections.abc import Callable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.RMB import RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementRequest,
    RMBBackend,
)
from sympleq.applications.randomized_benchmarking.backends.exponential import (
    ExponentialBackend,
    asymptote_value,
)
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    bare_hqc_at_register_width,
    batch_hqc_cost_at_physical_width,
    measured_items,
    print_experiment_summary,
    print_progress,
    spend_request_batch,
    start_run,
    stitched_batch_hqc,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_plots import (
    plot_boundary_uncertainty_surface,
    plot_boundary_surface,
    plot_live_volume_history,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_references import (
    analytic_lindblad_gates as _analytic_lindblad_gates,
    calibrated_surface_gates as _calibrated_surface_gates,
    calibrated_surface_label as _calibrated_surface_label,
    gp_grid_surface_gates as _gp_grid_surface_gates,
    gp_grid_surface_label as _gp_grid_surface_label,
    load_gp_grid_surface as _load_gp_grid_surface,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    REFERENCE_OFFSET,
    REFERENCE_SLOPE,
    parametric_boundary_fit_score,
    parametric_bootstrap_analytic_coverage,
)
from sympleq.applications.randomized_benchmarking.experiments.scores import (
    integrate_trapezoid as _trapezoid,
)
from sympleq.core.noise.noise_model import GenericNoise


BASE_1Q_PAULI_ERROR = 0.000025
BASE_2Q_PAULI_ERROR = 0.00079
_EPS = 1e-12
_LN2 = np.log(2.0)

# Natural parameter-vector layout.  L1, L2 are the one-/two-qubit per-gate
# log-error rates expanded linearly about the reference register size Qref:
#     Li(Q) = Li0 + mi (Q - Qref).
_I_L1, _I_L2 = 0, 1   # rate intercepts at Qref  (positive; log-transformed)
_I_M1, _I_M2 = 2, 3   # linear Q-slopes  dLi/dQ
_I_N1, _I_N2 = 4, 5   # quadratic Q-curvatures  d^2 Li/dQ^2 (pinned to 0 by default)
_I_V = 6              # visibility V in (0, 1]   (log-transformed)
_N_PARAMS = 7
# Logistic-equivalent sharpness of the RB exponential at the 0.5 crossing
# (slope-matched: the universal curve 2^{-n/n*} has slope -ln2/2 there, equal to
# a logistic of sharpness 2 ln 2).  Reported to the shared single-Q plotting
# code, which expects a sharpness; the boundary location itself uses only the
# rate parameters and does not depend on it.
_PHYSICAL_SHARPNESS = 2.0 * np.log(2.0)
_CHECKPOINT_SCHEMA_VERSION = 4
LiveVolumePoint = tuple[int, float, float, float, float]


# --------------------------------------------------------------------------- #
# Backend and configuration helpers
# --------------------------------------------------------------------------- #
def sympleq_backend_factory(
    settings: CrossingSettings,
    rng: RNGGenerator,
) -> RMBBackend:
    """Local SympleQ backend with the default benchmarking noise model.

    The same backend instance serves every number of qubits; SympleqBackend
    builds each circuit from ``request.config`` and the noise models are not
    dimension-specific.
    """
    noise_model = GenericNoise.from_paulis([BASE_1Q_PAULI_ERROR] * 3, rng)
    two_qubit_noise_model = GenericNoise.from_paulis([BASE_2Q_PAULI_ERROR] * 3, rng)
    return SympleqBackend(
        noise_model=noise_model,
        two_qubit_noise_model=two_qubit_noise_model,
    )


def exponential_backend_factory(
    settings: CrossingSettings,
    rng: RNGGenerator,
) -> RMBBackend:
    """Direct positive-control backend for the assumed RB exponential link."""
    one_q_error = getattr(settings, "initial_one_q_pauli_error", BASE_1Q_PAULI_ERROR)
    two_q_error = getattr(settings, "initial_two_q_pauli_error", BASE_2Q_PAULI_ERROR)
    one_q_scale = float(one_q_error) / BASE_1Q_PAULI_ERROR
    two_q_scale = float(two_q_error) / BASE_2Q_PAULI_ERROR
    one_q_rate = REFERENCE_OFFSET * one_q_scale
    two_q_rate = one_q_rate + REFERENCE_SLOPE * two_q_scale
    q_values = tuple(getattr(settings, "q_values", (settings.n_qubits,)))
    q_reference = float(np.median(q_values))
    return ExponentialBackend(
        one_qubit_rate=one_q_rate,
        two_qubit_rate=two_q_rate,
        one_qubit_q_slope=float(getattr(settings, "exponential_one_q_q_slope", 0.0)),
        two_qubit_q_slope=float(getattr(settings, "exponential_two_q_q_slope", 0.0)),
        q_reference=q_reference,
        visibility=float(getattr(settings, "initial_visibility", 1.0)),
        asymptote_model=getattr(settings, "asymptote_model", "depolarizing"),
    )


@lru_cache(maxsize=None)
def make_config_q(
    settings: "CostAwareSurfaceSettings",
    n_gates: float,
    ratio: float,
    n_qubits: int,
) -> RMBConfig:
    """Map (total gates, two-qubit ratio, n_qubits) to a realizable RMBConfig.

    Memoized because the greedy batch assembler evaluates batch_hqc_cost over
    many candidate (n, r, Q) probes repeatedly; rebuilding configs each time
    would dominate runtime.
    """
    n_2qb_gates = max(0, 2 * round(ratio * n_gates / 2))
    min_1qb_gates = 2 * n_qubits if settings.use_scrambler else 0
    n_1qb_gates = max(min_1qb_gates, 2 * round((n_gates - n_2qb_gates) / 2))
    return (
        RMBConfig.default()
        .with_n_qubits(int(n_qubits))
        .with_n_1qb_gates(int(n_1qb_gates))
        .with_n_2qb_gates(int(n_2qb_gates))
        .with_random_elimination(settings.random_elimination)
        .with_use_scrambler(settings.use_scrambler)
    )


# --------------------------------------------------------------------------- #
# Settings
# --------------------------------------------------------------------------- #
EXPERIMENTS_DIR = Path(__file__).resolve().parent
DEFAULT_SURFACE_PLOT_PATH = EXPERIMENTS_DIR / "figs" / "boundary_surface_3d.png"
DEFAULT_LIVE_SURFACE_PLOT_PATH = (
    EXPERIMENTS_DIR / "figs" / "boundary_surface_3d_live.png"
)
DEFAULT_LIVE_VOLUME_PLOT_PATH = (
    EXPERIMENTS_DIR / "figs" / "boundary_volume_live.png"
)
DEFAULT_UNCERTAINTY_PLOT_PATH = (
    EXPERIMENTS_DIR / "figs" / "boundary_surface_3d_uncertainty.png"
)
DEFAULT_GRID_VS_CONTINUOUS_PLOT_PATH = (
    EXPERIMENTS_DIR / "figs" / "boundary_surface_3d_grid_vs_continuous.png"
)


def _with_run_suffix(
    path: str | Path | None,
    *,
    backend_model: str,
    rng_seed: int | None,
    max_qubit_window: int | None,
    q_values: tuple[int, ...],
) -> str | Path | None:
    """Append Q-grid/window/backend/seed to output file names for local runs."""
    if path is None or rng_seed is None:
        return path
    if backend_model not in {"sympleq", "exponential"}:
        return path
    path_type = Path if isinstance(path, Path) else str
    p = Path(path)
    qwin = 0 if max_qubit_window is None else int(max_qubit_window)
    q_suffix = _q_values_suffix(q_values)
    suffix = f"_{q_suffix}_maxqwin{qwin}_{backend_model}_seed{rng_seed}"
    if p.stem.endswith(suffix):
        return path
    new_path = p.with_name(f"{p.stem}{suffix}{p.suffix}")
    return new_path if path_type is Path else str(new_path)


def _q_values_suffix(q_values: tuple[int, ...]) -> str:
    """Compact filename component for the exact requested scoring Q grid."""
    qs = tuple(int(q) for q in q_values)
    if not qs:
        return "qnone"
    if len(qs) == 1:
        return f"q{qs[0]}"
    steps = [b - a for a, b in zip(qs[:-1], qs[1:])]
    if steps and len(set(steps)) == 1:
        return f"q{qs[0]}-{qs[-1]}s{steps[0]}"
    return "q" + "-".join(str(q) for q in qs)


@dataclass(frozen=True)
class CostAwareSurfaceSettings(CrossingSettings):
    """Settings for cost-aware global surface design.

    The survival probability uses the physically-motivated RB exponential link
    p(n, r, Q) = V * (1 - B(Q)) * 2^{-n D(r,Q)} + B(Q), where D = 1/n* is the
    inverse-boundary rate, B(Q) is the asymptote, and V is the global
    visibility.  The boundary n* = 1/D sits at n D = 1 (renormalised fidelity
    0.5) independently of V and B.  The prior centre for (L1, L2) is the analytic
    line implied by the configured initial Pauli errors; the rate Q-slopes
    (m1, m2) are centred at zero (the agnostic "no Q dependence" prior) with a
    weakly informative width.
    """

    save_path: str | Path | None = "physics_informed_cost_aware_surface_design.json"
    # Save an RMB-format measurement JSON after every submitted batch and resume
    # from it on restart.  A small sibling checkpoint file stores budget state;
    # the main JSON remains the same measurement-outcome format as FLE_*.json.
    checkpoint_after_batch: bool = True
    resume_from_save: bool = True
    # Where the 3-D boundary-surface plot is written (multi-Q runs with plot=True).
    surface_plot_path: str | Path | None = DEFAULT_SURFACE_PLOT_PATH
    surface_plot_show: bool = False
    surface_plot_analytic: bool = True
    # Plot-only z-axis/gate-count bounds.  Leave None to use n_gates_bounds;
    # set e.g. (100, 3000) to zoom the 3-D surface without changing acquisition.
    surface_plot_n_gates_bounds: tuple[int, int] | None = None
    # Optional second comparison/reference surface produced by
    # calibrate_sympleq_effective_surface.py.  Analytic Lindblad remains the
    # primary reference; this adds an extra calibrated grid/score/wireframe.
    calibrated_surface_path: str | Path | None = None
    # Optional probability grid from another method.  Expected arrays are
    # qubits_axis, ratio_axis, gates_axis, probabilities, and target.  The
    # p=target contour is extracted along gates and shown as an extra surface.
    gp_grid_surface_path: str | Path | None = None
    gp_grid_surface_label: str | None = None
    # If enabled, refresh the 3-D boundary-surface plot after batches during
    # the run. The fixed-Q threshold-ambiguity diagnostic is written only at
    # run end by the final plot block.
    # The PNG path is updated in place, so it is safe for long unattended runs.
    # ``live_surface_plot_show`` additionally opens a non-blocking Matplotlib
    # window when the local plotting backend supports it.
    live_surface_plot: bool = False
    live_surface_plot_path: str | Path | None = DEFAULT_LIVE_SURFACE_PLOT_PATH
    live_surface_plot_every: int = 1
    live_surface_plot_show: bool = False
    live_surface_plot_pause: float = 0.25
    # Live history of the integrated volume under the fitted n*(r,Q) surface.
    # The dashed reference line is the analytic exponential surface built from
    # the known simulation noise scales (one_q_noise_scale/two_q_noise_scale),
    # not from the current posterior estimate.
    live_volume_plot: bool = False
    live_volume_plot_path: str | Path | None = DEFAULT_LIVE_VOLUME_PLOT_PATH
    live_volume_plot_every: int = 1
    live_volume_plot_show: bool = False
    live_volume_plot_pause: float = 0.25
    # Static +-1 sigma boundary-uncertainty surface plot (written at run end).
    # Draws the posterior-mean log-boundary surface enveloped by the +-k sigma
    # surfaces of log n*(r,Q), where k = uncertainty_plot_sigma.
    surface_uncertainty_plot: bool = True
    surface_uncertainty_plot_path: str | Path | None = DEFAULT_UNCERTAINTY_PLOT_PATH
    surface_uncertainty_plot_show: bool = False
    surface_uncertainty_sigma: float = 1.0
    # Plot the *raw* survival p=0.5 boundary (real fidelity) instead of the
    # renormalised-fidelity boundary n* = 1/D.  The renormalised boundary sits
    # where 2^{-nD} = 0.5; the raw survival p = V(1-B)2^{-nD}+B crosses 0.5
    # deeper, at
    #     n_raw(r,Q) = n*(r,Q) * ( -log2[ (0.5 - B(Q)) / (V (1 - B(Q))) ] ),
    # a per-Q rescale depending only on the asymptote B(Q) and visibility V (not
    # on r).  Where the raw curve never reaches 0.5 (B(Q) >= 0.5, or the
    # un-decayed top V(1-B)+B < 0.5) the crossing does not exist and those (r,Q)
    # are masked.  This remaps the three 3-D boundary-surface plots (main,
    # +-sigma uncertainty, grid-vs-continuous).  The S1/S2/volume *scores* and
    # any analytic/calibrated/GP reference overlays drawn by the shared plotting
    # module stay in the renormalised convention (they cannot be remapped from
    # here); at large Q (B ~ 0) the factor -> 1, so the two nearly coincide.
    plot_raw_fidelity_boundary: bool = True
    # Continuous MAP re-fit (grid-initialised) run once at the end on all
    # gathered data.  It optimises the *same* free parameters as the grid (the
    # axes with grid_resolution > 1), holding the pinned axes at their prior
    # centre, but over a continuous parameter space -- so the boundary and the
    # curvature nu2 land between grid ticks and the reported goodness-of-fit is
    # evaluated at the true optimum rather than at a quantised grid node.  A
    # separate comparison plot overlays the grid-median and continuous surfaces.
    continuous_refit: bool = True
    grid_vs_continuous_plot_path: str | Path | None = (
        DEFAULT_GRID_VS_CONTINUOUS_PLOT_PATH
    )
    backend_factory: Callable[[CrossingSettings, RNGGenerator], RMBBackend] = (
        sympleq_backend_factory
    )
    # Convenience selector for built-in local backends.  Custom code may still
    # pass ``backend_factory`` directly; setting this to "exponential" selects
    # the exact Bernoulli simulator for the assumed RB exponential link.
    backend_model: str = "sympleq"

    # --- region scored / measured -----------------------------------------
    q_values: tuple[int, ...] = (3, 4, 5, 6, 7, 8)
    score_ratio_points: int = 6
    # Optional ground-truth surface for scoring: callable(r, Q) -> n*(r, Q).
    truth_boundary: Callable[[float, int], float] | None = None

    # --- budgets -----------------------------------------------------------
    target_log_rms: float = 0.04       # stop when sqrt(mean Var[log n*]) <= this
    max_iterations: int = 400

    # --- prior on the boundary surface -------------------------------------
    initial_one_q_pauli_error: float = BASE_1Q_PAULI_ERROR # * 0.6
    initial_two_q_pauli_error: float = BASE_2Q_PAULI_ERROR # * 0.6
    initial_error_relative_uncertainty: float = 0.30
    initial_one_q_error_relative_uncertainty: float | None = None
    initial_two_q_error_relative_uncertainty: float | None = None
    # Visibility V in [0, 1] is the global SPAM amplitude of the RB exponential
    # link  p = V*(1-B)*2^{-n D} + B.  V = 1 is the no-SPAM / ideal-simulator
    # case (and is pinned by default via grid_resolution[6] = 1); free it (set
    # grid_resolution[6] > 1) to fit SPAM on hardware.
    initial_visibility: float = 1.0
    visibility_log_std: float = 0.20
    # Asymptote B(Q) of the survival curve.  "depolarizing" -> B = 2^{-Q} (the
    # physically correct value for raw register-survival RB); "zero" -> B = 0
    # (use only if the stored fidelity is already asymptote-subtracted); a float
    # overrides with a constant.  See the deep-probe diagnostic to verify.
    asymptote_model: str | float = "depolarizing"
    # Width of the rate Q-slope prior as a fraction of the rate over the Q span.
    q_slope_prior_fraction: float = 0.30
    # Width of the rate Q-curvature (quadratic) prior, as a fraction of the rate
    # over the squared half-span.  Only active when the nu axes are unpinned
    # (grid_resolution entries 4, 5 > 1); held in reserve for hardware whose
    # Q-dependence may be nonlinear.
    q_curve_prior_fraction: float = 0.30
    # Optional Q-dependence of the direct ExponentialBackend truth rates.  These
    # are natural-log rate slopes dL_i/dQ, not posterior-prior widths.
    exponential_one_q_q_slope: float = 0.0
    exponential_two_q_q_slope: float = 0.0

    def __post_init__(self) -> None:
        if self.live_surface_plot_show and not self.live_surface_plot:
            object.__setattr__(self, "live_surface_plot", True)
        if self.live_volume_plot_show and not self.live_volume_plot:
            object.__setattr__(self, "live_volume_plot", True)
        if self.surface_uncertainty_plot_show and not self.surface_uncertainty_plot:
            object.__setattr__(self, "surface_uncertainty_plot", True)
        if self.backend_model in {"sympleq", "exponential"}:
            for attr in (
                "save_path",
                "surface_plot_path",
                "live_surface_plot_path",
                "live_volume_plot_path",
                "surface_uncertainty_plot_path",
                "grid_vs_continuous_plot_path",
            ):
                object.__setattr__(
                    self,
                    attr,
                    _with_run_suffix(
                        getattr(self, attr),
                        backend_model=self.backend_model,
                        rng_seed=self.rng_seed,
                        max_qubit_window=getattr(self, "max_qubit_window", 0),
                        q_values=tuple(int(q) for q in self.q_values),
                    ),
                )
        if self.backend_model == "sympleq":
            return
        if self.backend_model == "exponential":
            object.__setattr__(self, "backend_factory", exponential_backend_factory)
            return
        raise ValueError(
            "backend_model must be 'sympleq' or 'exponential' "
            f"(got {self.backend_model!r})"
        )

    # --- grid posterior ----------------------------------------------------
    # Axes are (log L1, log L2, m1, m2, nu1, nu2, log V).  By default the two
    # curvature axes (nu1, nu2) and the visibility axis are pinned (resolution
    # 1): the curvatures to 0 (linear-in-Q rates) and V to initial_visibility
    # (no SPAM).  Free the curvature axes (resolution > 1 on entries 4, 5) only
    # for a device whose Q-dependence is expected to be nonlinear -- e.g. on
    # hardware where crosstalk or spectator errors can make lambda(Q) curve.
    # Free the visibility axis (entry 6 > 1) only where SPAM is appreciable.
    grid_resolution: tuple[int, ...] = (11, 11, 7, 7, 1, 1, 1)
    grid_halfwidth_sigmas: float = 3.0
    grid_halfwidth_floor_fraction: float = 0.35
    posterior_min_configs: int = 4
    denom_floor: float = 1e-6
    # >1 widens the posterior.  With the physically-correct exponential link the
    # posterior should be close to calibrated, so this defaults to 1; raise it
    # only if the 90% bootstrap coverage still comes out below nominal.
    likelihood_temperature: float = 1.0
    # Finer grid used only by the (out-of-loop) scoring/plotting hooks so the
    # bootstrap bands are not quantised by the coarse design-loop grid.
    boundary_fit_resolution: tuple[int, ...] = (17, 17, 9, 9, 1, 1, 1)

    # --- stopping ----------------------------------------------------------
    # Stop after `stop_patience` consecutive batches that each reduce the
    # integrated log-boundary uncertainty by less than `min_rel_improvement`
    # (relative).  This makes spend track marginal value instead of always
    # exhausting the budget on an unreachable absolute target.
    min_rel_improvement: float = 0.01
    stop_patience: int = 3

    # --- acquisition (BALD: parameter mutual information per HQC) -----------
    boundary_probe_quantiles: tuple[float, ...] = (0.25, 0.5, 0.75)
    # Multiples of the median predicted n* at which to also place probes; under
    # the RB exponential link a measurement is most informative about the rate
    # somewhat deeper than the crossing, so these reach a little beyond n*.  The
    # deepest multiple is kept modest so the large-n* / low-r / high-Q corner
    # (where the global linear fit already extrapolates the boundary) does not
    # attract wastefully deep, expensive probes.
    probe_depth_fractions: tuple[float, ...] = (0.7, 1.0, 1.3)
    # Hard cap on any probe depth, as a fraction of the deepest measurable gate
    # count, so a single high-n* slice cannot spawn runaway-depth circuits.
    max_probe_depth_fraction: float = 0.6
    # Local depth cap: no probe deeper than this multiple of the slice's own
    # predicted median n*.  Sits above the deepest probe_depth_fraction (1.3) so
    # genuine near-boundary and information-peak probes are untouched, but clips
    # the broad-posterior boundary quantiles that would otherwise reach 2-4x the
    # crossing in the low-r / high-Q corner -- past where survival has saturated
    # and a measurement is uninformative.  This is the cap that actually bounds
    # the corner, since it scales with the local boundary rather than the global
    # gate bound.  Set to 0.0 to use only the global gate-bound cap above.
    max_probe_depth_nstar_multiple: float = 1.6
    # Drop acquisition probes whose requested (n, r, Q) cannot be realised by
    # RMBConfig after even-gate rounding and the 2Q scrambler one-qubit floor.
    # If the scrambler floor forces the config above the requested total-gate
    # count, the point is unreachable at that depth.  If rounding shifts the
    # realised two-qubit ratio by more than this tolerance, the point belongs to
    # a different slice and should not be scored as the requested r.
    acquisition_realized_ratio_tolerance: float = 0.02
    # Acquisition floor on the model-predicted survival probability.  Candidate
    # depths whose posterior-mean predicted survival falls below this are dropped
    # before scoring: deep in the tail the outcome is near-certain failure
    # (p -> B(Q)), so such shots carry little information about the boundary
    # location (the Fisher information peaks near the crossing, not the tail) yet
    # are the most sensitive to the asymptote B(Q).  Capping in p-space (not in
    # n) is self-adjusting across (r, Q): it bites only where the model already
    # expects failure, so it cannot waste credits and behaves identically on the
    # simulator and on hardware.  Set to 0.0 to disable.  Keep modest (~0.15) so
    # enough near-asymptote coverage remains to constrain B(Q) if it is freed.
    min_predicted_survival: float = 0.15
    # Acquisition CAP on the model-predicted renormalised fidelity F = 2^{-nD}.
    # Candidate depths whose posterior-mean predicted F exceeds this are dropped
    # before scoring: at small n the survival sits near the un-decayed top of the
    # curve (F -> 1), which is exactly where the single-exponential link is least
    # valid -- the random circuit is not yet twirled, so an early-depth point
    # carries a finite-depth transient (an apparent F(0) > 1 intercept) that
    # biases the fitted slope rather than informing it.  Capping in F-space
    # (renormalised, so independent of B(Q) and V) is self-adjusting across
    # (r, Q): it always trims the same portion of the *decay* regardless of how
    # fast a slice decays.  Together with min_predicted_survival this keeps probes
    # in the informative mid-fidelity band (roughly F in [0.15, 0.85]).  Set to
    # 1.0 to disable (no upper cap).
    max_predicted_fidelity: float = 0.85
    # Predicted-survival band (p_bar below this) treated as the near-asymptote
    # "deep tail" by the passive residual diagnostic.  Reported but never acted
    # on: a systematically one-signed mean residual in this band indicates the
    # asymptote B(Q) is mis-modelled (the shared rates are being dragged to
    # reconcile a tail the link cannot hit); a centred-but-scattered tail means
    # the points are merely over-leveraged.  Truth-free, so it reads identically
    # on the simulator and on hardware.
    tail_diagnostic_pmax: float = 0.25
    shot_diminish_rho: float = 0.6     # per-(r,Q)-slice diminishing returns
    max_shots_per_probe: int = 8
    warmup_iterations: int = 2
    coverage_bonus: float = 3.0        # early bonus for unmeasured (r,Q) slices
    # --- joint (Q, r, n) acquisition ---------------------------------------
    # The acquisition chooses (Q, r, n) per probe by BALD-per-HQC.  By default
    # the candidate Q set is settings.q_values and the candidate r set is the
    # score grid, i.e. the historical behaviour.  Set acquisition_q_resolution > 0
    # to let the design choose Q from a denser set spanning the Q range (so it can
    # load probes where the rate Q-slope/curvature is least constrained -- the Q
    # extremes for the slope, the interior for curvature), and
    # acquisition_ratio_points to use a denser r set than the scoring grid.
    acquisition_q_resolution: int = 0
    acquisition_ratio_points: int = 0
    # Distinct-Q coverage floor.  Until at least this many distinct Q values have
    # been measured, unmeasured Q values receive the coverage bonus, so a free-Q
    # design cannot collapse onto the two Q extremes (the optimal-but-blind design
    # for a strictly linear-in-Q rate) before the Q axis is spanned.  This is the
    # guard that keeps interior Q sampled so nonlinearity can be detected; pair a
    # high value with freed nu axes on hardware.  0 disables (historical
    # behaviour: only the warmup (r,Q)-slice bonus applies).
    min_distinct_q_coverage: int = 0
    # --- hardware Q-window constraint (e.g. Quantinuum H2) ------------------
    # On a device where one stitched submission runs on a *single physical
    # register* of constant width, a circuit that uses fewer qubits than the
    # register leaves the rest idle.  A few idle qubits are fine, but too many
    # inject crosstalk/memory error into the active qubits (biasing the measured
    # q-qubit decay) and waste the wider register.  ``max_qubit_window`` caps the
    # spread of Q within a single batch: all probes in one submission must have
    # their Q inside a contiguous window of this width, so the physical width is
    # W = max(Q in batch) and no config leaves more than max_qubit_window - 1
    # idle qubits.  Example: max_qubit_window = 3 admits a batch spanning
    # {17, 18, 19} (W = 19, at most 2 idle).  The acquisition still scores the
    # whole (Q, r, n) grid; it then packs the *best feasible window*: each
    # candidate window is greedily packed and the window whose packed batch
    # delivers the most information (realized BALD, after diminishing returns and
    # the shared base cost) is submitted.  Over successive batches the myopic
    # per-window choice tiles Q, since a measured window loses its uncertainty.
    # 0 disables the constraint (a single batch may span any Q); a value >= the
    # full Q span is equivalent to disabled.  NOTE: for a window to hold several
    # Q values the acquisition Q-candidates must be spaced finely enough to fall
    # inside it -- set acquisition_q_resolution so the candidate Q set is dense
    # (e.g. integer-spaced) rather than only q_values.  The batch cost is always
    # charged at the physical register width W = max Q in the batch, independent
    # of this cap (see batch_hqc_cost_at_physical_width), matching how H2 prices a
    # submission; the window only bounds the idle-qubit *bias*, not the price.
    max_qubit_window: int = 0
    # When the Q-window constraint is active, exact-pack only this many
    # top-ranked candidate windows (ranked by a cheap per-slice pre-score) before
    # choosing the best by realized packed value.  This bounds the extra cost of
    # the per-window packing; 0 exact-packs every window (most thorough, slowest).
    qubit_window_shortlist: int = 3

    # ---- plotting / scoring compatibility (single-Q slice at reference Q) --
    def boundary_fit(self, data) -> tuple[float, float, float] | None:
        params, weights = _stateless_posterior(data, self)
        if params is None:
            return None
        q = max(int(round(self.n_qubits)), 1)
        dq = q - _q_reference(self)
        lam1, lam2 = _lambdas(params, dq)
        q_eff = lam1 / _LN2
        s_eff = (lam2 - lam1) / _LN2
        med = (
            float(_wquantile(q_eff, weights, 0.5)),
            float(_wquantile(s_eff, weights, 0.5)),
            _PHYSICAL_SHARPNESS,
        )
        return med

    def boundary_bootstrap(
        self,
        data,
        *,
        n_bootstrap: int = 100,
        seed: int | None = None,
    ) -> list[tuple[float, float, float]]:
        params, weights = _stateless_posterior(data, self)
        if params is None:
            return []
        q = max(int(round(self.n_qubits)), 1)
        dq = q - _q_reference(self)
        lam1, lam2 = _lambdas(params, dq)
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(params), size=n_bootstrap, replace=True, p=weights)
        return [
            (
                float(lam1[i] / _LN2),
                float((lam2[i] - lam1[i]) / _LN2),
                _PHYSICAL_SHARPNESS,
            )
            for i in idx
        ]


# --------------------------------------------------------------------------- #
# Prior, transform, grid
# --------------------------------------------------------------------------- #
def _initial_relative_uncertainties(
    settings: CostAwareSurfaceSettings,
) -> tuple[float, float]:
    one_q = (
        settings.initial_error_relative_uncertainty
        if settings.initial_one_q_error_relative_uncertainty is None
        else settings.initial_one_q_error_relative_uncertainty
    )
    two_q = (
        settings.initial_error_relative_uncertainty
        if settings.initial_two_q_error_relative_uncertainty is None
        else settings.initial_two_q_error_relative_uncertainty
    )
    return max(one_q, 0.0), max(two_q, 0.0)


def _prior_moments(
    settings: CostAwareSurfaceSettings,
) -> tuple[np.ndarray, np.ndarray]:
    """Prior centre and std in transformed space
    t = (log L1, log L2, m1, m2, nu1, nu2, log V).

    The model is the two-rate RB exponential

        1/n*(r,Q) = [ L1(Q) (1-r) + L2(Q) r ] / ln 2,
        Li(Q) = Li0 + mi (Q - Qref) + nu_i (Q - Qref)^2,

    so L10, L20 are the one-/two-qubit per-gate log-error rates at the reference
    register size Qref, m1, m2 their linear Q-dependence, and nu1, nu2 the
    optional quadratic curvature (pinned to zero by default).  The rate centres
    come from the configured initial Pauli-error guesses; crucially the
    two-qubit rate centre includes the one-qubit rate, so the boundary slope
    (L2 - L1)/ln 2 carries the physically-correct -L1 term.  The Q-slopes and
    curvatures are centred at zero (agnostic) with weakly informative widths.
    """
    one_q_scale = settings.initial_one_q_pauli_error / BASE_1Q_PAULI_ERROR
    two_q_scale = settings.initial_two_q_pauli_error / BASE_2Q_PAULI_ERROR
    lam1_prior = REFERENCE_OFFSET * one_q_scale
    lam2_prior = lam1_prior + REFERENCE_SLOPE * two_q_scale
    u1, u2 = _initial_relative_uncertainties(settings)

    q_span = max(max(settings.q_values) - min(settings.q_values), 1)
    m1_std = settings.q_slope_prior_fraction * lam1_prior / q_span
    m2_std = settings.q_slope_prior_fraction * lam2_prior / q_span
    # Quadratic curvature prior width: fraction of the rate over the squared
    # half-span, so a nu (Q-Qref)^2 excursion at the Q-range edge is comparable
    # in size to the configured fraction of the rate.  Only active if the nu
    # axes are unpinned (grid_resolution[4], [5] > 1).
    q_half = max(0.5 * q_span, 1.0)
    n1_std = settings.q_curve_prior_fraction * lam1_prior / (q_half * q_half)
    n2_std = settings.q_curve_prior_fraction * lam2_prior / (q_half * q_half)

    centre = np.array(
        [
            np.log(max(lam1_prior, 1e-12)),
            np.log(max(lam2_prior, 1e-12)),
            0.0,
            0.0,
            0.0,
            0.0,
            np.log(min(max(settings.initial_visibility, 1e-6), 1.0)),
        ],
        dtype=float,
    )
    std = np.array(
        [
            max(0.05, np.log1p(u1)),
            max(0.05, np.log1p(u2)),
            max(1e-9, m1_std),
            max(1e-9, m2_std),
            max(1e-12, n1_std),
            max(1e-12, n2_std),
            max(0.01, settings.visibility_log_std),
        ],
        dtype=float,
    )
    return centre, std


def _params_from_t(t: np.ndarray) -> np.ndarray:
    """Map transformed coordinates to natural parameters (L10, L20, m1, m2, V).

    The two rate intercepts and the visibility are exp of their coordinates
    (positive); the Q-slopes m are linear; V is clipped to (0, 1] so the RB link
    stays a valid probability.
    """
    params = np.empty_like(t)
    params[:, _I_L1] = np.exp(t[:, _I_L1])
    params[:, _I_L2] = np.exp(t[:, _I_L2])
    params[:, _I_M1] = t[:, _I_M1]
    params[:, _I_M2] = t[:, _I_M2]
    params[:, _I_N1] = t[:, _I_N1]
    params[:, _I_N2] = t[:, _I_N2]
    params[:, _I_V] = np.clip(np.exp(t[:, _I_V]), 1e-6, 1.0)
    return params


def _build_grid(
    centre_grid: np.ndarray,
    halfwidth_grid: np.ndarray,
    settings: CostAwareSurfaceSettings,
    resolution: tuple[int, ...] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build an axis-aligned tensor grid in transformed space.

    The prior is fixed (``_prior_moments``); only the grid *layout* moves so we
    spend resolution where the posterior actually is.  ``resolution`` overrides
    ``settings.grid_resolution`` (used by the finer scoring/plotting hooks).
    """
    prior_centre, prior_std = _prior_moments(settings)
    res = settings.grid_resolution if resolution is None else resolution
    axes = [
        (
            np.array([centre_grid[d]])
            if n <= 1
            else np.linspace(
                centre_grid[d] - halfwidth_grid[d],
                centre_grid[d] + halfwidth_grid[d],
                n,
            )
        )
        for d, n in enumerate(res)
    ]
    mesh = np.meshgrid(*axes, indexing="ij")
    t = np.stack([m.ravel() for m in mesh], axis=1)
    z = (t - prior_centre) / prior_std
    log_prior = -0.5 * np.sum(z * z, axis=1)
    params = _params_from_t(t)
    return params, log_prior, t


def _axis_layout(
    t: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    prior_std: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Axis-aligned next layout: posterior mean centre and per-axis half-width."""
    mean_t = np.sum(weights[:, None] * t, axis=0)
    var_t = np.sum(weights[:, None] * (t - mean_t) ** 2, axis=0)
    std_t = np.sqrt(np.maximum(var_t, 0.0))
    k = settings.grid_halfwidth_sigmas
    floor = settings.grid_halfwidth_floor_fraction * k * prior_std
    halfwidth = np.maximum(k * std_t, floor)
    return mean_t, halfwidth


def _region_corners(settings: CostAwareSurfaceSettings) -> np.ndarray:
    """The four (r, Q) corners of the scored rectangle.

    The inverse boundary is bilinear in (r, Q), so its minimum over the
    rectangle is attained at a corner; checking the four corners therefore
    certifies positivity (admissibility) everywhere in the region exactly.
    """
    r_lo, r_hi = settings.ratio_bounds
    q_lo, q_hi = min(settings.q_values), max(settings.q_values)
    return np.array(
        [[r_lo, q_lo], [r_lo, q_hi], [r_hi, q_lo], [r_hi, q_hi]],
        dtype=float,
    )


def _q_reference(settings: CostAwareSurfaceSettings) -> float:
    """Reference register size at which the rate intercepts are defined.

    The midpoint of the Q range decorrelates each rate intercept from its slope.
    """
    return 0.5 * (min(settings.q_values) + max(settings.q_values))


def _lambdas(params: np.ndarray, dq: float) -> tuple[np.ndarray, np.ndarray]:
    """Per-gate rates at register offset dq = Q - Qref, for every parameter row.

        Li(Q) = Li0 + mi (Q - Qref) + nu_i (Q - Qref)^2,   i = 1, 2.

    The quadratic coefficients nu_i (columns _I_N1, _I_N2) are pinned to zero by
    default (their grid axes have resolution 1), so this reduces to the linear
    law unless the curvature axes are explicitly freed for a device whose
    Q-dependence is expected to be nonlinear (e.g. crosstalk on hardware).
    """
    dq2 = dq * dq
    lam1 = params[:, _I_L1] + params[:, _I_M1] * dq + params[:, _I_N1] * dq2
    lam2 = params[:, _I_L2] + params[:, _I_M2] * dq + params[:, _I_N2] * dq2
    return lam1, lam2


def _denominator(
    params: np.ndarray,
    r: float,
    q: float,
    settings: CostAwareSurfaceSettings,
) -> np.ndarray:
    """Inverse boundary 1/n*(r, Q) for every parameter row.

        1/n* = [ L1(Q) (1-r) + L2(Q) r ] / ln 2,
        Li(Q) = Li0 + mi (Q - Qref) + nu_i (Q - Qref)^2.
    """
    dq = q - _q_reference(settings)
    lam1, lam2 = _lambdas(params, dq)
    return (lam1 * (1.0 - r) + lam2 * r) / _LN2


# --------------------------------------------------------------------------- #
# Data, likelihood, posterior
# --------------------------------------------------------------------------- #
def _measured_arrays(data) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Measured configs as (ratio, n_gates, n_qubits, successes, failures)."""
    measured = measured_items(data)
    ratios, gates, qubits, succ, fail = [], [], [], [], []
    for config, estimator in measured:
        counts = estimator.counts()
        ratios.append(float(config.ratio_2_qb_gates))
        gates.append(float(config.n_gates))
        qubits.append(float(config.n_qubits))
        succ.append(float(counts.get(True, 0)))
        fail.append(float(counts.get(False, 0)))
    return (
        np.asarray(ratios), np.asarray(gates), np.asarray(qubits),
        np.asarray(succ), np.asarray(fail),
    )


def _asymptote(settings: CostAwareSurfaceSettings, q: float) -> float:
    """Survival asymptote B(Q) of the RB decay."""
    return asymptote_value(settings.asymptote_model, q)


def _link_probability(
    params: np.ndarray,
    n: float,
    r: float,
    q: float,
    settings: CostAwareSurfaceSettings,
) -> np.ndarray:
    """RB exponential survival probability at (n, r, Q) for every parameter row.

        p = V * (1 - B(Q)) * 2^{-n D(r,Q)} + B(Q),     D = 1 / n*.

    The transition sharpness is the physical constant ln 2 inside the exponent
    (not a free parameter); V is the visibility (params[:, _I_V]) and B(Q) the
    asymptote.  At n D = 1 the renormalised fidelity (p - B)/(V(1-B)) equals 0.5,
    so the boundary n* = 1/D is independent of V and B.
    """
    d = _denominator(params, r, q, settings)
    b = _asymptote(settings, q)
    visibility = params[:, _I_V]
    fidelity = np.power(2.0, -n * d)  # 2^{-n D} in (0, 1]
    p = visibility * (1.0 - b) * fidelity + b
    return np.clip(p, _EPS, 1.0 - _EPS)


def _grid_posterior(
    data,
    params: np.ndarray,
    log_prior: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> tuple[np.ndarray, float]:
    """Return normalised posterior weights over the grid and the ESS.

    Parameter rows whose denominator is non-positive anywhere in the scored
    (r, Q) box are excluded (the boundary must be finite and positive over the
    region we report).  Because the denominator is bilinear, its minimum over a
    rectangle is attained at a corner, so checking the four corners is exact.
    """
    log_post = log_prior.copy()

    corners = _region_corners(settings)
    valid = np.ones(len(params), dtype=bool)
    for r, q in corners:
        valid &= _denominator(params, r, q, settings) > settings.denom_floor

    ratios, gates, qubits, succ, fail = _measured_arrays(data)
    log_lik = np.zeros(len(params))
    for i in range(len(ratios)):
        pij = _link_probability(params, gates[i], ratios[i], qubits[i], settings)
        log_lik += succ[i] * np.log(pij + _EPS) + fail[i] * np.log(1.0 - pij + _EPS)

    temperature = max(settings.likelihood_temperature, 1e-6)
    log_post = log_prior + log_lik / temperature
    log_post = np.where(valid, log_post, -np.inf)
    if not np.any(np.isfinite(log_post)):
        weights = np.full(len(params), 1.0 / len(params))
        return weights, float(len(params))

    log_post -= np.max(log_post[np.isfinite(log_post)])
    weights = np.where(np.isfinite(log_post), np.exp(log_post), 0.0)
    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0.0:
        weights = np.full(len(params), 1.0 / len(params))
    else:
        weights /= total
    ess = 1.0 / float(np.sum(weights**2))
    return weights, ess


def _stateless_posterior(
    data,
    settings: CostAwareSurfaceSettings,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    """Posterior on a prior-centred grid (used by plotting/scoring hooks)."""
    if _measured_config_count(data) < settings.posterior_min_configs:
        return None, None
    centre, std = _prior_moments(settings)
    halfwidth = settings.grid_halfwidth_sigmas * std
    params, log_prior, _ = _build_grid(
        centre, halfwidth, settings, resolution=settings.boundary_fit_resolution
    )
    weights, _ = _grid_posterior(data, params, log_prior, settings)
    return params, weights


def _measured_config_count(data) -> int:
    return len(measured_items(data))


# --------------------------------------------------------------------------- #
# Weighted statistics
# --------------------------------------------------------------------------- #
def _wquantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    keep = np.isfinite(values) & (weights > 0.0)
    if not np.any(keep):
        return float("nan")
    v, w = values[keep], weights[keep]
    order = np.argsort(v)
    v, w = v[order], w[order]
    cum = np.cumsum(w)
    cum /= cum[-1]
    return float(np.interp(q, cum, v))


def _binary_entropy(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, _EPS, 1.0 - _EPS)
    return -p * np.log(p) - (1.0 - p) * np.log(1.0 - p)


def _effective_free_params(settings: CostAwareSurfaceSettings) -> int:
    """Number of unpinned posterior axes (grid resolution > 1).

    Used as the nominal parameter count k when forming chi^2 / (N - k).  This is
    only a nominal dof: the prior constrains the parameters, so the *effective*
    number of parameters is somewhat smaller, but with a few hundred configs the
    difference is negligible for a health check.
    """
    return int(sum(1 for res in settings.grid_resolution if int(res) > 1))


# --------------------------------------------------------------------------- #
# Surface uncertainty (the scoring-aligned objective)
# --------------------------------------------------------------------------- #
def _score_grid(settings: CostAwareSurfaceSettings) -> list[tuple[float, int]]:
    r_lo, r_hi = settings.ratio_bounds
    ratios = np.linspace(r_lo, r_hi, settings.score_ratio_points)
    return [(float(r), int(q)) for q in settings.q_values for r in ratios]


def _acquisition_slices(settings: CostAwareSurfaceSettings) -> list[tuple[float, int]]:
    """Candidate (r, Q) slices the acquisition may probe.

    Decoupled from the *scoring* grid (``_score_grid``): the surface is always
    scored over ``settings.q_values`` x ``score_ratio_points``, but the design is
    free to *measure* at a denser set if requested, so it can place probes where
    the rate Q-slope/curvature is least constrained rather than only at the
    scored nodes.  Defaults reproduce the historical behaviour (candidate Q set =
    ``q_values``, candidate r set = the scoring grid).
    """
    r_lo, r_hi = settings.ratio_bounds
    n_r = settings.acquisition_ratio_points or settings.score_ratio_points
    ratios = np.linspace(r_lo, r_hi, max(int(n_r), 2))
    q_res = int(settings.acquisition_q_resolution)
    if q_res > 0:
        q_lo, q_hi = min(settings.q_values), max(settings.q_values)
        qs = sorted({int(round(v)) for v in np.linspace(q_lo, q_hi, q_res)})
    else:
        qs = [int(v) for v in settings.q_values]
    return [(float(r), int(q)) for q in qs for r in ratios]


def _surface_log_uncertainty(
    params: np.ndarray,
    weights: np.ndarray,
    grid: list[tuple[float, int]],
    settings: CostAwareSurfaceSettings,
) -> tuple[float, dict[tuple[float, int], float]]:
    """Mean over (r,Q) of posterior Var[log n*]; also the per-node variance.

    Shares the vectorised log-n* machinery used by the acquisition objective.
    Over the admissible region every weighted node has a positive denominator,
    so this matches the previous per-node loop exactly.
    """
    if not grid:
        return 0.0, {}
    log_n, log_n_sq = _contour_log_matrices(params, grid, settings)
    node_var = _node_log_variances(log_n, log_n_sq, weights)
    per_node = {node: float(v) for node, v in zip(grid, node_var)}
    return float(np.mean(node_var)), per_node


# --------------------------------------------------------------------------- #
# Candidate probes and BALD acquisition
# --------------------------------------------------------------------------- #
def _candidate_n_values(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    r: float,
    q: int,
) -> list[int]:
    """Even gate counts spanning the informative depths at (r, Q).

    Candidates combine quantiles of the predicted boundary *location* (covering
    posterior spread) with fixed multiples of the median predicted n* (covering
    the depth where a measurement is most informative about the rate).  For the
    RB exponential link this information peak sits deeper than n* (around
    1/ln 2 times n*, pushed further by Bernoulli weighting), so probing only at
    n* would miss it once the boundary location is well localised.
    """
    d = _denominator(params, r, q, settings)
    ok = (d > 0.0) & (weights > 0.0)
    if not np.any(ok):
        return []
    boundary = 1.0 / d[ok]
    w = weights[ok]
    n_lo, n_hi = settings.n_gates_bounds
    median_nstar = _wquantile(boundary, w, 0.5)
    global_cap = min(n_hi, max(n_lo, int(settings.max_probe_depth_fraction * n_hi)))
    mult = settings.max_probe_depth_nstar_multiple
    if mult > 0.0 and np.isfinite(median_nstar) and median_nstar > 0:
        local_cap = int(mult * median_nstar)
        depth_cap = max(n_lo, min(global_cap, local_cap))
    else:
        depth_cap = global_cap
    gates = []
    for quantile in settings.boundary_probe_quantiles:
        val = _wquantile(boundary, w, quantile)
        if np.isfinite(val):
            gates.append(int(np.clip(2 * round(val / 2), n_lo, depth_cap)))
    if np.isfinite(median_nstar):
        for fraction in settings.probe_depth_fractions:
            val = fraction * median_nstar
            gates.append(int(np.clip(2 * round(val / 2), n_lo, depth_cap)))
    candidates = sorted(set(g for g in gates if g > 0))

    floor = settings.min_predicted_survival
    if floor > 0.0 and candidates:
        kept = []
        for g in candidates:
            pij = _link_probability(params, g, r, q, settings)
            if float(np.sum(weights * pij)) >= floor:
                kept.append(g)
        candidates = kept if kept else candidates[:1]

    fmax = settings.max_predicted_fidelity
    if 0.0 < fmax < 1.0 and candidates:
        kept = []
        for g in candidates:
            fid = np.power(2.0, -g * np.maximum(d, 0.0))   # renormalised F per row
            if float(np.sum(weights * fid)) <= fmax:
                kept.append(g)
        candidates = kept if kept else candidates[-1:]
    return candidates


def _realizable_probe_config(
    settings: CostAwareSurfaceSettings,
    n: int,
    r: float,
    q: int,
) -> RMBConfig | None:
    """Return the submitted config if it realises requested (n,r,Q), else None."""
    config = make_config_q(settings, n, r, q)
    if int(config.n_gates) != int(n):
        return None
    tol = max(float(settings.acquisition_realized_ratio_tolerance), 0.0)
    if abs(float(config.ratio_2_qb_gates) - float(r)) > tol:
        return None
    return config


def _bald_per_shot(
    params: np.ndarray,
    weights: np.ndarray,
    n: float,
    r: float,
    q: int,
    settings: CostAwareSurfaceSettings,
) -> float:
    """Mutual information between one Bernoulli outcome at (n,r,Q) and the parameters."""
    pij = _link_probability(params, n, r, q, settings)
    p_bar = float(np.sum(weights * pij))
    expected_cond_entropy = float(np.sum(weights * _binary_entropy(pij)))
    return float(_binary_entropy(np.array([p_bar]))[0] - expected_cond_entropy)


def _contour_log_matrices(
    params: np.ndarray,
    score_grid: list[tuple[float, int]],
    settings: CostAwareSurfaceSettings,
) -> tuple[np.ndarray, np.ndarray]:
    """Return L and L^2 where L[g, k] = log n*(node g) for parameter row k."""
    log_n = np.zeros((len(score_grid), len(params)))
    for g, (r, q) in enumerate(score_grid):
        d = _denominator(params, r, q, settings)
        ok = d > 0.0
        log_n[g, ok] = -np.log(d[ok])
    return log_n, log_n * log_n


def _node_log_variances(
    log_n: np.ndarray,
    log_n_sq: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    """Posterior variance of log n* at each scored (r,Q) node (vectorised)."""
    mean = log_n @ weights
    mean_sq = log_n_sq @ weights
    return np.maximum(mean_sq - mean * mean, 0.0)


def _pack_candidates(
    candidates: list[dict],
    settings: CostAwareSurfaceSettings,
    affordable_cap: float,
) -> tuple[list[MeasurementRequest], float, float]:
    """Greedily pack one stitched submission from a candidate pool.

    The batch is priced at the *physical register width* W = max Q in the batch
    (the H2 model): a smaller-Q circuit sharing the submission pays for the idle
    qubits of the wider register.  Because widening the batch reprices every
    circuit already in it, the greedy density (value per marginal HQC) sees a
    high-Q probe as expensive when the rest of the batch is low-Q, so the packer
    prefers Q-homogeneous submissions on its own.
    """
    if not candidates:
        return [], 0.0, 0.0

    rho = settings.shot_diminish_rho
    shots: dict[tuple[int, float, int], int] = {}   # (n, r, q) -> shots
    slice_shots: dict[tuple[float, int], int] = {}
    config_by_key = {
        (cand["n"], cand["r"], cand["q"]): cand["config"] for cand in candidates
    }

    def requests_from(shots_map: dict) -> list[MeasurementRequest]:
        return [
            MeasurementRequest(config_by_key[(n, r, q)], s)
            for (n, r, q), s in shots_map.items()
            if s > 0
        ]

    def cost_of(shots_map: dict) -> float:
        items = [(n, r, q, s) for (n, r, q), s in shots_map.items() if s > 0]
        if not items:
            return 0.0
        width = max(q for (_, _, q, _) in items)
        total_bare = sum(
            s * bare_hqc_at_register_width(config_by_key[(n, r, q)], width)
            for (n, r, q, s) in items
        )
        return stitched_batch_hqc(total_bare)

    current_cost = 0.0
    total_value = 0.0
    while True:
        best = None
        best_density = 0.0
        best_value_gain = 0.0
        for cand in candidates:
            key = (cand["n"], cand["r"], cand["q"])
            if shots.get(key, 0) >= settings.max_shots_per_probe:
                continue
            marginal_value = cand["value"] * (rho ** slice_shots.get(cand["slice"], 0))
            trial = dict(shots)
            trial[key] = trial.get(key, 0) + 1
            trial_cost = cost_of(trial)
            if trial_cost > affordable_cap:
                continue
            marginal_cost = trial_cost - current_cost
            if marginal_cost <= 0.0:
                marginal_cost = _EPS
            density = marginal_value / marginal_cost
            if density > best_density:
                best_density = density
                best = (key, cand["slice"], trial_cost)
                best_value_gain = marginal_value
        if best is None:
            break
        key, slc, trial_cost = best
        shots[key] = shots.get(key, 0) + 1
        slice_shots[slc] = slice_shots.get(slc, 0) + 1
        current_cost = trial_cost
        total_value += best_value_gain

    return requests_from(shots), total_value, current_cost


def _qubit_windows(candidate_qs, width: int) -> list[tuple[int, ...]]:
    """Maximal contiguous Q windows of the given width over the candidate Q set."""
    qs = sorted({int(q) for q in candidate_qs})
    raw: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for top in qs:
        lo = top - int(width) + 1
        grp = tuple(q for q in qs if lo <= q <= top)
        if grp and grp not in seen:
            seen.add(grp)
            raw.append(grp)
    maximal: list[tuple[int, ...]] = []
    for grp in raw:
        gset = set(grp)
        if not any(gset < set(other) for other in raw):
            maximal.append(grp)
    return maximal


def _assemble_batch(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    affordable_cap: float,
    measured_slices: set[tuple[float, int]],
    early: bool,
) -> list[MeasurementRequest]:
    """Assemble one stitched submission, optionally within a single Q-window."""
    acq_slices = _acquisition_slices(settings)
    measured_q = {q for (_, q) in measured_slices}
    need_q_coverage = len(measured_q) < int(settings.min_distinct_q_coverage)

    candidates: list[dict] = []
    for r, q in acq_slices:
        for n in _candidate_n_values(params, weights, settings, r, q):
            config = _realizable_probe_config(settings, n, r, q)
            if config is None:
                continue
            value = _bald_per_shot(params, weights, n, r, q, settings)
            if value <= 0.0:
                continue
            slice_new = (r, q) not in measured_slices
            q_new = q not in measured_q
            bonus = 1.0
            if early and slice_new:
                bonus = settings.coverage_bonus
            if need_q_coverage and q_new:
                bonus = max(bonus, settings.coverage_bonus)
            candidates.append(
                {
                    "n": n,
                    "r": r,
                    "q": q,
                    "config": config,
                    "slice": (r, q),
                    "value": value * bonus,
                }
            )
    if not candidates:
        return []

    width = int(getattr(settings, "max_qubit_window", 0) or 0)
    if width <= 0:
        requests, _, _ = _pack_candidates(candidates, settings, affordable_cap)
        return requests

    windows = _qubit_windows({c["q"] for c in candidates}, width)
    if len(windows) <= 1:
        requests, _, _ = _pack_candidates(candidates, settings, affordable_cap)
        return requests

    def _prescore(group: tuple[int, ...]) -> float:
        gset = set(group)
        best_per_slice: dict[tuple[float, int], float] = {}
        for c in candidates:
            if c["q"] in gset:
                s = c["slice"]
                if c["value"] > best_per_slice.get(s, 0.0):
                    best_per_slice[s] = c["value"]
        return float(sum(best_per_slice.values()))

    windows_ranked = sorted(windows, key=_prescore, reverse=True)
    shortlist = int(getattr(settings, "qubit_window_shortlist", 0) or 0)
    if shortlist > 0:
        windows_ranked = windows_ranked[:shortlist]

    best_requests: list[MeasurementRequest] = []
    best_value = -1.0
    for group in windows_ranked:
        gset = set(group)
        group_candidates = [c for c in candidates if c["q"] in gset]
        requests, value, _ = _pack_candidates(
            group_candidates, settings, affordable_cap
        )
        if requests and value > best_value:
            best_value = value
            best_requests = requests
    return best_requests


# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #
def _residual_diagnostics(
    data,
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> dict:
    """Standardised-residual + goodness-of-fit diagnostics over measured points.

    ``chi2_dof`` is the Pearson chi^2 / (N - k); ``dev_dof`` the Bernoulli
    deviance / (N - k).  Evaluated at the posterior-mean predicted p (single-row
    ``params`` with weight 1 gives the value at a point estimate, e.g. the
    continuous MAP fit).
    """
    ratios, gates, qubits, succ, fail = _measured_arrays(data)
    r_lo, r_hi = settings.ratio_bounds
    r_mid = 0.5 * (r_lo + r_hi)
    pmax = settings.tail_diagnostic_pmax

    max_abs = 0.0
    n_flag = 0
    lowr_sum = highr_sum = 0.0
    lowr_n = highr_n = 0
    tail_n = 0
    chi2 = 0.0
    deviance = 0.0
    n_points = 0
    for i in range(len(ratios)):
        total = succ[i] + fail[i]
        if total <= 0:
            continue
        n_points += 1
        pij = _link_probability(params, gates[i], ratios[i], qubits[i], settings)
        p_bar = float(np.sum(weights * pij))
        p_bar = min(max(p_bar, _EPS), 1.0 - _EPS)
        p_hat = succ[i] / total
        sd = np.sqrt(max(p_bar * (1.0 - p_bar) / total, _EPS))
        z = (p_hat - p_bar) / sd          # signed
        max_abs = max(max_abs, abs(z))
        n_flag += int(abs(z) > 3.0)
        if ratios[i] <= r_mid:
            lowr_sum += z
            lowr_n += 1
        else:
            highr_sum += z
            highr_n += 1
        tail_n += int(p_hat < pmax)
        chi2 += z * z
        s, f = succ[i], fail[i]
        if s > 0:
            deviance += 2.0 * s * np.log(max(p_hat, _EPS) / p_bar)
        if f > 0:
            deviance += 2.0 * f * np.log(max(1.0 - p_hat, _EPS) / (1.0 - p_bar))

    k = _effective_free_params(settings)
    dof = max(n_points - k, 1)
    return {
        "max_abs": max_abs,
        "n_flag": n_flag,
        "lowr_mean": lowr_sum / lowr_n if lowr_n else 0.0,
        "lowr_n": lowr_n,
        "highr_mean": highr_sum / highr_n if highr_n else 0.0,
        "highr_n": highr_n,
        "tail_n": tail_n,
        "chi2": float(chi2),
        "deviance": float(deviance),
        "dof": int(dof),
        "n_points": int(n_points),
        "chi2_dof": float(chi2 / dof),
        "dev_dof": float(deviance / dof),
    }


def _surface_plot_mesh(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Posterior-median boundary surface as log10 gates on a plot mesh."""
    med = np.asarray(
        [float(_wquantile(params[:, j], weights, 0.5)) for j in range(_N_PARAMS)]
    )
    return _boundary_log10_from_param_vector(med, settings)


def _surface_uncertainty_plot_mesh(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Posterior mean boundary and +-1 sigma surfaces as log10 gates."""
    r_lo, r_hi = settings.ratio_bounds
    q_lo, q_hi = min(settings.q_values), max(settings.q_values)
    n_bounds = (
        settings.n_gates_bounds
        if settings.surface_plot_n_gates_bounds is None
        else settings.surface_plot_n_gates_bounds
    )
    n_lo, n_hi = n_bounds
    grid_r = np.linspace(r_lo, r_hi, 40)
    grid_q = np.linspace(q_lo, q_hi, 40)
    rr, qq = np.meshgrid(grid_r, grid_q)

    mean_log = np.full(rr.shape, np.nan, dtype=float)
    sigma_log = np.full(rr.shape, np.nan, dtype=float)
    flat_r = rr.ravel()
    flat_q = qq.ravel()
    mean_flat = mean_log.ravel()
    sigma_flat = sigma_log.ravel()
    raw = bool(getattr(settings, "plot_raw_fidelity_boundary", False))
    for idx in range(flat_r.size):
        q = float(flat_q[idx])
        d = _denominator(params, float(flat_r[idx]), q, settings)
        ok = (d > 0.0) & (weights > 0.0)
        if not np.any(ok):
            continue
        w = weights[ok]
        log_n = -np.log(d[ok])                    # log n* (renormalised) per row
        if raw:
            # log n_raw = log n* + log g(Q, V_row); this folds any V posterior
            # spread into the band and masks rows with no raw p=0.5 crossing.
            b = _asymptote(settings, q)
            vrow = params[ok, _I_V]
            denom = vrow * (1.0 - b)
            with np.errstate(divide="ignore", invalid="ignore"):
                x = np.where(denom > 0.0, (0.5 - b) / denom, np.nan)
                g = np.where((x > 0.0) & (x < 1.0), -np.log2(x), np.nan)
            good = np.isfinite(g) & (g > 0.0)
            if not np.any(good):
                continue
            log_n = log_n[good] + np.log(g[good])
            w = w[good]
        w = w / np.sum(w)
        mean = float(np.sum(w * log_n))
        var = float(np.sum(w * log_n * log_n) - mean * mean)
        mean_flat[idx] = mean
        sigma_flat[idx] = float(np.sqrt(max(var, 0.0)))

    k = float(settings.surface_uncertainty_sigma)
    ln10 = np.log(10.0)

    def _clip_log10(natural_log_n: np.ndarray) -> np.ndarray:
        gates = np.exp(natural_log_n)
        z = natural_log_n / ln10
        return np.where((gates >= n_lo) & (gates <= n_hi), z, np.nan)

    mean_z = _clip_log10(mean_log)
    lower_z = _clip_log10(mean_log - k * sigma_log)
    upper_z = _clip_log10(mean_log + k * sigma_log)
    if not np.any(np.isfinite(mean_z)):
        return None
    return rr, qq, mean_z, lower_z, upper_z


# --------------------------------------------------------------------------- #
# Continuous MAP re-fit (grid-initialised) and grid-vs-continuous comparison
# --------------------------------------------------------------------------- #
def _boundary_log10_from_param_vector(
    param_vector: np.ndarray,
    settings: CostAwareSurfaceSettings,
    *,
    n_r: int = 48,
    n_q: int = 48,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """log10 n*(r, Q) mesh for a single natural parameter vector.

    Shared by the grid-median surface plot and the continuous-fit surface plot so
    both are drawn with identical meshing and gate-count clipping.
    """
    p = np.asarray(param_vector, dtype=float)
    r_lo, r_hi = settings.ratio_bounds
    q_lo, q_hi = min(settings.q_values), max(settings.q_values)
    n_bounds = (
        settings.n_gates_bounds
        if settings.surface_plot_n_gates_bounds is None
        else settings.surface_plot_n_gates_bounds
    )
    n_lo, n_hi = n_bounds
    grid_r = np.linspace(r_lo, r_hi, n_r)
    grid_q = np.linspace(q_lo, q_hi, n_q)
    rr, qq = np.meshgrid(grid_r, grid_q)
    dq = qq - _q_reference(settings)
    lam1 = p[_I_L1] + p[_I_M1] * dq + p[_I_N1] * dq * dq
    lam2 = p[_I_L2] + p[_I_M2] * dq + p[_I_N2] * dq * dq
    d = (lam1 * (1.0 - rr) + lam2 * rr) / _LN2
    gates = np.where(d > 0.0, 1.0 / np.maximum(d, _EPS), np.nan)
    if getattr(settings, "plot_raw_fidelity_boundary", False):
        # Remap n* (renormalised 2^{-nD}=0.5 crossing) onto the raw survival
        # p=0.5 depth n_raw = n* * g(Q), with visibility V from this parameter
        # vector.  g is nan where the raw 0.5 crossing does not exist (masked).
        gates = gates * _raw_fidelity_gate_factor(settings, qq, float(p[_I_V]))
    log10_gates = np.where(
        np.isfinite(gates) & (gates >= n_lo) & (gates <= n_hi),
        np.log10(np.where(np.isfinite(gates) & (gates > 0.0), gates, np.nan)),
        np.nan,
    )
    if not np.any(np.isfinite(log10_gates)):
        return None
    return rr, qq, log10_gates


def _raw_fidelity_gate_factor(
    settings: CostAwareSurfaceSettings,
    qq: np.ndarray,
    visibility: float,
) -> np.ndarray:
    """Per-(r,Q) factor g(Q) mapping the renormalised boundary onto raw p=0.5.

    The renormalised boundary n* sits where 2^{-nD} = 0.5.  The raw survival
    p = V(1-B)2^{-nD} + B crosses 0.5 at n_raw = n* * g(Q), with

        g(Q) = -log2[ (0.5 - B(Q)) / (V (1 - B(Q))) ].

    g depends only on the asymptote B(Q) and visibility V (not on r); it is >= 1
    wherever the crossing exists and -> 1 as B -> 0 (large Q), so the raw and
    renormalised boundaries nearly coincide at high Q.  Returns nan where no raw
    p=0.5 depth exists (B(Q) >= 0.5, or the un-decayed top V(1-B)+B < 0.5), i.e.
    where 0 < (0.5-B)/(V(1-B)) < 1 fails.
    """
    v = float(visibility)
    uq = np.unique(qq)
    b_of_q = {float(q): float(_asymptote(settings, float(q))) for q in uq}
    b = np.vectorize(b_of_q.get, otypes=[float])(qq)
    denom = v * (1.0 - b)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(denom > 0.0, (0.5 - b) / denom, np.nan)
        g = np.where((x > 0.0) & (x < 1.0), -np.log2(x), np.nan)
    return g


def _free_axis_indices(settings: CostAwareSurfaceSettings) -> list[int]:
    """Transformed-space axes the continuous fit optimises (grid resolution > 1).

    The continuous re-fit deliberately fits the *same* model the grid did: the
    axes pinned in the grid (resolution 1) stay pinned at their prior centre, so
    e.g. nu2 is only fitted when the grid freed it (nesting the linear model
    inside the quadratic one).
    """
    return [d for d, res in enumerate(settings.grid_resolution) if int(res) > 1]


def _continuous_map_fit(
    data,
    settings: CostAwareSurfaceSettings,
    t_init: np.ndarray,
) -> tuple[np.ndarray | None, object]:
    """Continuous MAP fit over the free transformed axes, initialised at t_init.

    Maximises log_likelihood + log_prior (the same objective the grid scores),
    but over a continuous parameter vector rather than a lattice, so the boundary
    and any freed curvature nu2 land between grid ticks.  Pinned axes are held at
    the prior centre.  Optimisation runs in prior-standardised coordinates so the
    wildly different natural scales (log-rates ~ 1, slopes ~ 1e-6, curvatures ~
    1e-8) are all O(1) and Nelder-Mead is well conditioned.  Returns
    (t_hat, optimiser_result); t_hat is None if there is nothing to fit.
    """
    free = _free_axis_indices(settings)
    if not free:
        return None, None
    ratios, gates, qubits, succ, fail = _measured_arrays(data)
    if ratios.size == 0:
        return None, None

    prior_centre, prior_std = _prior_moments(settings)
    b = np.array([_asymptote(settings, q) for q in qubits], dtype=float)
    qref = _q_reference(settings)
    dq = qubits - qref
    corners = _region_corners(settings)
    corner_dq = corners[:, 1] - qref
    corner_r = corners[:, 0]

    # Base vector: pinned axes fixed at their prior centre; free axes overwritten.
    t_base = t_init.astype(float).copy()
    for d in range(_N_PARAMS):
        if d not in free:
            t_base[d] = prior_centre[d]

    def t_from_z(z_free: np.ndarray) -> np.ndarray:
        t = t_base.copy()
        for j, d in enumerate(free):
            t[d] = prior_centre[d] + z_free[j] * prior_std[d]
        return t

    def neg_log_post(z_free: np.ndarray) -> float:
        t = t_from_z(z_free)
        p = _params_from_t(t[None, :])[0]
        # Admissibility over the scored rectangle (bilinear -> corners suffice).
        cl1 = p[_I_L1] + p[_I_M1] * corner_dq + p[_I_N1] * corner_dq * corner_dq
        cl2 = p[_I_L2] + p[_I_M2] * corner_dq + p[_I_N2] * corner_dq * corner_dq
        cd = (cl1 * (1.0 - corner_r) + cl2 * corner_r) / _LN2
        if not np.all(cd > settings.denom_floor):
            return 1e18
        lam1 = p[_I_L1] + p[_I_M1] * dq + p[_I_N1] * dq * dq
        lam2 = p[_I_L2] + p[_I_M2] * dq + p[_I_N2] * dq * dq
        d = (lam1 * (1.0 - ratios) + lam2 * ratios) / _LN2
        if not np.all(d > 0.0):
            return 1e18
        fid = np.power(2.0, -gates * d)
        pr = np.clip(p[_I_V] * (1.0 - b) * fid + b, _EPS, 1.0 - _EPS)
        loglik = float(np.sum(succ * np.log(pr) + fail * np.log(1.0 - pr)))
        zvec = (t - prior_centre) / prior_std
        logprior = -0.5 * float(np.sum(zvec * zvec))
        return -(loglik + logprior)

    z0 = np.array(
        [(t_init[d] - prior_centre[d]) / prior_std[d] for d in free], dtype=float
    )
    try:
        from scipy.optimize import minimize
    except Exception as exc:  # noqa: BLE001 - scipy is the only external need here
        raise RuntimeError(
            "continuous re-fit requires scipy.optimize (install scipy)"
        ) from exc
    result = minimize(
        neg_log_post,
        z0,
        method="Nelder-Mead",
        options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 8000, "maxfev": 16000},
    )
    t_hat = t_from_z(np.asarray(result.x, dtype=float))
    return t_hat, result


def _report_continuous_fit(
    data,
    settings: CostAwareSurfaceSettings,
    t_hat: np.ndarray,
    grid_params: np.ndarray,
    grid_weights: np.ndarray,
    result: object,
) -> None:
    """Print the continuous MAP parameters, its goodness-of-fit, and grid deltas."""
    params_hat = _params_from_t(np.asarray(t_hat, dtype=float)[None, :])
    a = params_hat[0]
    diag = _residual_diagnostics(data, params_hat, np.array([1.0]), settings)
    k = _effective_free_params(settings)
    free = _free_axis_indices(settings)
    names = ["L1", "L2", "m1", "m2", "nu1", "nu2", "V"]
    converged = bool(getattr(result, "success", True)) if result is not None else True
    print("\nContinuous MAP re-fit (grid-initialised; free axes = grid_resolution>1)")
    print(
        "  free parameters: "
        + ", ".join(names[d] for d in free)
        + f"  ({'converged' if converged else 'did NOT fully converge'})"
    )
    print(
        f"  MAP rates @Qref (L1, L2; m1, m2; nu1, nu2; V) = "
        f"({a[_I_L1]:.3g}, {a[_I_L2]:.3g}; {a[_I_M1]:.3g}, {a[_I_M2]:.3g}; "
        f"{a[_I_N1]:.3g}, {a[_I_N2]:.3g}; {a[_I_V]:.3g})"
    )
    print(
        "  goodness of fit @ MAP: "
        f"chi2/dof={diag['chi2_dof']:.3f} deviance/dof={diag['dev_dof']:.3f} "
        f"(N={diag['n_points']}, k={k}, dof={diag['dof']}, "
        f"max|z|={diag['max_abs']:.1f}, flagged={diag['n_flag']})"
    )
    gmed = [
        float(_wquantile(grid_params[:, j], grid_weights, 0.5))
        for j in range(_N_PARAMS)
    ]
    print("  grid median vs continuous MAP (natural params):")
    for j, name in enumerate(names):
        marker = "  <- freed" if j in free else ""
        print(f"    {name:>3s}: grid={gmed[j]:+.4g}  cont={a[j]:+.4g}{marker}")


def plot_grid_vs_continuous_surface(
    grid_param_vector: np.ndarray,
    continuous_param_vector: np.ndarray,
    measured_arrays: tuple[np.ndarray, ...],
    settings: CostAwareSurfaceSettings,
    png_path: str | Path | None = None,
    *,
    show: bool = False,
    show_block: bool = True,
    show_pause: float = 0.001,
    close: bool = True,
    figure_name: str | None = None,
):
    """Overlay the grid-median and continuous-MAP boundary surfaces in log10 n.

    Blue is the grid posterior-median surface (quantised by the design lattice);
    orange is the continuous MAP surface (between ticks).  Where they separate is
    where grid quantisation moved the estimate -- typically along a freed axis
    such as the curvature nu2.  Measured configs are scattered for context.
    """
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

    grid_mesh = _boundary_log10_from_param_vector(grid_param_vector, settings)
    cont_mesh = _boundary_log10_from_param_vector(continuous_param_vector, settings)
    if grid_mesh is None and cont_mesh is None:
        return None

    fig = plt.figure(num=figure_name, figsize=(9.5, 7.5), clear=True)
    ax = fig.add_subplot(projection="3d")

    if grid_mesh is not None:
        rr, qq, zg = grid_mesh
        ax.plot_surface(
            rr, qq, zg, color="tab:blue", alpha=0.30,
            linewidth=0, antialiased=True, rstride=1, cstride=1, shade=False,
        )
        ax.plot_wireframe(
            rr, qq, zg, color="tab:blue", alpha=0.55,
            linewidth=0.6, rstride=4, cstride=4,
        )
    if cont_mesh is not None:
        rr, qq, zc = cont_mesh
        ax.plot_surface(
            rr, qq, zc, color="tab:orange", alpha=0.30,
            linewidth=0, antialiased=True, rstride=1, cstride=1, shade=False,
        )
        ax.plot_wireframe(
            rr, qq, zc, color="tab:orange", alpha=0.65,
            linewidth=0.6, rstride=4, cstride=4,
        )

    ratios, gates, qubits, succ, fail = measured_arrays
    total = succ + fail
    mask = total > 0
    if np.any(mask):
        ax.scatter(
            ratios[mask], qubits[mask], np.log10(np.maximum(gates[mask], _EPS)),
            c="0.25", s=10, alpha=0.5, depthshade=True,
        )

    handles = [
        Line2D([0], [0], color="tab:blue", lw=2, label="grid posterior median"),
        Line2D([0], [0], color="tab:orange", lw=2, label="continuous MAP fit"),
    ]
    ax.legend(handles=handles, loc="upper left")
    ax.set_xlabel("two-qubit ratio  r")
    ax.set_ylabel("n_qubits  Q")
    ax.set_zlabel(r"$\log_{10}(\mathrm{gate\ count}\ n)$")
    if getattr(settings, "plot_raw_fidelity_boundary", False):
        ax.set_title(
            "Grid-median vs continuous-MAP boundary surface\n"
            r"(raw survival $p=0.5$ depth $\log_{10} n_{\mathrm{raw}}(r, Q)$)"
        )
    else:
        ax.set_title(
            r"Grid-median vs continuous-MAP boundary surface  $\log_{10} n_*(r, Q)$"
        )
    fig.tight_layout()
    if png_path is not None:
        png_path = Path(png_path)
        png_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_path, dpi=150)
    if show:
        plt.show(block=show_block)
        if not show_block:
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(max(float(show_pause), 0.001))
    if close:
        plt.close(fig)
    return png_path


def _maybe_update_live_surface_plot(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    iteration: int,
) -> None:
    """Best-effort in-run refresh of the boundary surface plot."""
    if not settings.live_surface_plot:
        return
    every = max(1, int(settings.live_surface_plot_every))
    if iteration % every != 0:
        return
    if settings.live_surface_plot_path is None and not settings.live_surface_plot_show:
        return
    if _measured_config_count(rmb._data) < settings.posterior_min_configs:
        return

    try:
        params, weights = _stateless_posterior(rmb._data, settings)
        if params is None or weights is None:
            return
        mesh = _surface_plot_mesh(params, weights, settings)
        if mesh is None:
            return
        png_path = plot_boundary_surface(
            *mesh,
            _measured_arrays(rmb._data),
            settings,
            png_path=settings.live_surface_plot_path,
            show=settings.live_surface_plot_show,
            show_block=False,
            show_pause=settings.live_surface_plot_pause,
            close=not settings.live_surface_plot_show,
            figure_name="cost-aware live surface",
        )
        if png_path is not None:
            print_progress(
                settings,
                budget,
                f"live surface plot -> {png_path}",
            )
    except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
        print_progress(settings, budget, f"live surface plot skipped: {exc}")


def _maybe_update_live_volume_plot(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    iteration: int,
    history: list[LiveVolumePoint],
) -> None:
    """Best-effort live plot of fitted surface volume against iteration."""
    if not settings.live_volume_plot:
        return
    every = max(1, int(settings.live_volume_plot_every))
    if iteration % every != 0:
        return
    if settings.live_volume_plot_path is None and not settings.live_volume_plot_show:
        return
    if _measured_config_count(rmb._data) < settings.posterior_min_configs:
        return

    try:
        params, weights = _stateless_posterior(rmb._data, settings)
        if params is None or weights is None:
            return
        point = _record_volume_point_from_posterior(
            params,
            weights,
            settings,
            iteration,
            history,
        )
        if point is None:
            return
        png_path = plot_live_volume_history(
            history,
            settings,
            png_path=settings.live_volume_plot_path,
            show=settings.live_volume_plot_show,
            show_block=False,
            show_pause=settings.live_volume_plot_pause,
            close=not settings.live_volume_plot_show,
            figure_name="cost-aware live volume",
        )
        if png_path is not None:
            _, fitted_volume, true_volume, _, _ = point
            ratio = fitted_volume / true_volume if true_volume else float("nan")
            print_progress(
                settings,
                budget,
                f"live volume plot -> {png_path} (ratio {ratio:.4g})",
            )
    except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
        print_progress(settings, budget, f"live volume plot skipped: {exc}")


def _record_volume_point_from_posterior(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    iteration: int,
    history: list[LiveVolumePoint],
) -> LiveVolumePoint | None:
    """Record fitted/reference/+-1 sigma volume for this posterior state."""
    scores = _surface_s1_s2(
        params,
        weights,
        settings,
        _volume_reference_gates(settings),
    )
    fitted_volume = float(scores["surface_volume_fit"])
    lower_volume = float(scores["surface_volume_lower_1sigma"])
    upper_volume = float(scores["surface_volume_upper_1sigma"])
    true_volume = float(scores["surface_volume_reference"])
    if not (np.isfinite(fitted_volume) and np.isfinite(true_volume)):
        return None
    _record_live_volume_point(
        history,
        iteration,
        fitted_volume,
        true_volume,
        lower_volume,
        upper_volume,
    )
    return (
        int(iteration),
        fitted_volume,
        true_volume,
        lower_volume,
        upper_volume,
    )


def _add_prior_live_volume_point(
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    grid_centre: np.ndarray,
    grid_halfwidth: np.ndarray,
    history: list[LiveVolumePoint],
) -> None:
    """Seed the live volume curve with the prior predictive volume at iteration 0."""
    if not settings.live_volume_plot:
        return
    if settings.live_volume_plot_path is None and not settings.live_volume_plot_show:
        return
    try:
        params, log_prior, _ = _build_grid(grid_centre, grid_halfwidth, settings)
        log_w = log_prior - np.max(log_prior)
        weights = np.exp(log_w)
        weights = weights / np.sum(weights)
        scores = _surface_s1_s2(
            params,
            weights,
            settings,
            _volume_reference_gates(settings),
        )
        fitted_volume = float(scores["surface_volume_fit"])
        lower_volume = float(scores["surface_volume_lower_1sigma"])
        upper_volume = float(scores["surface_volume_upper_1sigma"])
        true_volume = float(scores["surface_volume_reference"])
        if not (np.isfinite(fitted_volume) and np.isfinite(true_volume)):
            return
        _record_live_volume_point(
            history,
            0,
            fitted_volume,
            true_volume,
        )
        png_path = plot_live_volume_history(
            history,
            settings,
            png_path=settings.live_volume_plot_path,
            show=settings.live_volume_plot_show,
            show_block=False,
            show_pause=settings.live_volume_plot_pause,
            close=not settings.live_volume_plot_show,
            figure_name="cost-aware live volume",
        )
        if png_path is not None:
            ratio = fitted_volume / true_volume if true_volume else float("nan")
            print_progress(
                settings,
                budget,
                f"live volume prior -> {png_path} (ratio {ratio:.4g})",
            )
    except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
        print_progress(settings, budget, f"live volume prior skipped: {exc}")


def _volume_reference_gates(settings: CostAwareSurfaceSettings):
    """Reference surface used by the volume-vs-iteration history."""
    if settings.gp_grid_surface_path is not None:
        return _gp_grid_surface_gates
    return _analytic_lindblad_gates


def _record_live_volume_point(
    history: list[LiveVolumePoint],
    iteration: int,
    fitted_volume: float,
    true_volume: float,
    lower_volume: float = float("nan"),
    upper_volume: float = float("nan"),
) -> None:
    """Append or replace the live-volume point for an iteration."""
    point = (
        int(iteration),
        float(fitted_volume),
        float(true_volume),
        float(lower_volume),
        float(upper_volume),
    )
    for i, existing in enumerate(history):
        existing_iteration = existing[0]
        if int(existing_iteration) == int(iteration):
            history[i] = point
            return
    history.append(point)
    history.sort(key=lambda p: p[0])


def _config_checkpoint_record(config: RMBConfig) -> dict:
    """JSON-serializable description of an RMB config."""
    return {
        "n_gates": int(config.n_gates),
        "ratio_2qb_gates": float(config.ratio_2_qb_gates),
        "n_1qb_gates": int(config.n_1qb_gates),
        "n_2qb_gates": int(config.n_2qb_gates),
        "n_qubits": int(config.n_qubits),
        "random_elimination": float(config.random_elimination),
        "use_scrambler": bool(config.use_scrambler),
        "gates_set": [g.name for g in config.gates_set],
    }


def _config_from_checkpoint_record(record: dict) -> RMBConfig:
    """Rebuild the config fields needed by plotting/saving from checkpoint JSON."""
    return RMBConfig(
        n_1qb_gates=int(record["n_1qb_gates"]),
        n_2qb_gates=int(record["n_2qb_gates"]),
        n_qubits=int(record["n_qubits"]),
        random_elimination=float(record.get("random_elimination", 0.0)),
        use_scrambler=bool(record.get("use_scrambler", True)),
    )


def _batch_checkpoint_record(
    *,
    batch_number: int,
    iteration: int,
    cost_hqc: float,
    requests: list[MeasurementRequest],
    outcomes: dict[RMBConfig, list[bool]],
    first_shot_indices: dict[RMBConfig, int],
) -> dict:
    """Record every shot outcome in the submitted batch."""
    request_records = []
    for request in requests:
        config = request.config
        result_list = [bool(outcome) for outcome in outcomes.get(config, [])]
        first_shot = int(first_shot_indices.get(config, 0))
        request_records.append(
            {
                "config": _config_checkpoint_record(config),
                "requested_shots": int(request.shots),
                "first_shot_index": first_shot,
                "measurements": [
                    {"shot_index": first_shot + i, "outcome": bool(outcome)}
                    for i, outcome in enumerate(result_list)
                ],
            }
        )
    return {
        "batch": int(batch_number),
        "iteration": int(iteration),
        "cost_hqc": float(cost_hqc),
        "requested_shots": int(sum(req.shots for req in requests)),
        "requests": request_records,
    }


def _submitted_from_batch_history(batch_history: list[dict]) -> list[RMBConfig]:
    """Recover submitted configs in submission order from checkpoint metadata."""
    submitted: list[RMBConfig] = []
    for batch in batch_history:
        for request in batch.get("requests", []):
            config_record = request.get("config")
            if isinstance(config_record, dict):
                submitted.append(_config_from_checkpoint_record(config_record))
    return submitted


def _legacy_batch_history_from_data(data) -> list[dict]:
    """Represent pre-ledger aggregate data as one batch with unknown provenance."""
    requests = []
    total_shots = 0
    for config, estimator in measured_items(data):
        outcomes = []
        for outcome, count in estimator.counts().items():
            outcomes.extend([bool(outcome)] * int(count))
        total_shots += len(outcomes)
        requests.append(
            {
                "config": _config_checkpoint_record(config),
                "requested_shots": len(outcomes),
                "first_shot_index": 0,
                "measurements": [
                    {"shot_index": i, "outcome": bool(outcome)}
                    for i, outcome in enumerate(outcomes)
                ],
            }
        )
    if not requests:
        return []
    return [
        {
            "batch": None,
            "iteration": None,
            "source": "legacy_aggregate_import",
            "note": "Batch assignment was not present in this older checkpoint.",
            "cost_hqc": None,
            "requested_shots": total_shots,
            "requests": requests,
        }
    ]


def _live_volume_history_from_meta(meta: dict) -> list[LiveVolumePoint]:
    """Read persisted live-volume points, tolerating old checkpoint files."""
    history = []
    for point in meta.get("live_volume_history", []):
        try:
            history.append(
                (
                    int(point["iteration"]),
                    float(point["fitted_volume"]),
                    float(point["reference_volume"]),
                    float(point.get("lower_volume", float("nan"))),
                    float(point.get("upper_volume", float("nan"))),
                )
            )
        except (KeyError, TypeError, ValueError):
            continue
    history.sort(key=lambda p: p[0])
    return history


def _write_json_atomic(path: Path, payload: dict) -> None:
    """Write JSON through a sibling temp file to avoid partial metadata files."""
    tmp_path = path.with_name(f"{path.name}.tmp")
    tmp_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp_path.replace(path)


def _checkpoint_meta_path(save_path: str | Path) -> Path:
    base_path = resolve_data_path(save_path)
    return base_path.parent / f"{base_path.stem}_checkpoint.json"


def _settings_q_values(settings: CostAwareSurfaceSettings) -> list[int]:
    return [int(q) for q in settings.q_values]


def _checkpoint_q_values_match(meta: dict, settings: CostAwareSurfaceSettings) -> bool:
    try:
        saved_q_values = [int(q) for q in meta["settings"]["q_values"]]
    except (KeyError, TypeError, ValueError):
        return False
    return saved_q_values == _settings_q_values(settings)


def _save_measurement_checkpoint(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    submitted: list[RMBConfig],
    batch_history: list[dict],
    live_volume_history: list[LiveVolumePoint],
    *,
    iteration: int | None = None,
    reason: str = "checkpoint",
) -> None:
    """Save RMB-format measurements plus resumable batch/plot metadata."""
    if settings.save_path is None or not settings.checkpoint_after_batch:
        return
    try:
        rmb.save(settings.save_path)
        base_path = resolve_data_path(settings.save_path)
        meta = {
            "schema_version": _CHECKPOINT_SCHEMA_VERSION,
            "reason": reason,
            "iteration": iteration,
            "settings": {
                "q_values": _settings_q_values(settings),
                "backend_model": settings.backend_model,
                "rng_seed": settings.rng_seed,
                "max_qubit_window": int(settings.max_qubit_window),
            },
            "hqc_budget": settings.hqc_budget,
            "spent_hqc": budget.spent_hqc,
            "remaining_hqc": budget.remaining_hqc,
            "jobs": budget.jobs,
            "max_job_circuits": budget.max_job_circuits,
            "submitted_configs": len(submitted),
            "measured_configs": _measured_config_count(rmb._data),
            "batch_history": batch_history,
            "live_volume_history": [
                {
                    "iteration": int(i),
                    "fitted_volume": float(fitted),
                    "reference_volume": float(reference),
                    "lower_volume": float(lower),
                    "upper_volume": float(upper),
                }
                for i, fitted, reference, lower, upper in live_volume_history
            ],
        }
        _write_json_atomic(_checkpoint_meta_path(settings.save_path), meta)
        print_progress(settings, budget, f"saved measurements -> {base_path}")
    except Exception as exc:  # noqa: BLE001 - checkpointing should not kill a run
        print_progress(settings, budget, f"checkpoint save skipped: {exc}")


def _resume_measurement_checkpoint(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
) -> tuple[RMB, Budget, list[RMBConfig], list[dict], list[LiveVolumePoint]]:
    """Load saved measurements into the current backend and restore budget state."""
    if settings.save_path is None or not settings.resume_from_save:
        return rmb, budget, [], [], []
    base_path = resolve_data_path(settings.save_path)
    if not base_path.exists():
        return rmb, budget, [], [], []

    meta: dict | None = None
    meta_path = _checkpoint_meta_path(settings.save_path)
    if not meta_path.exists():
        print_progress(
            settings,
            budget,
            "resume skipped: checkpoint metadata missing, so q_values cannot be verified",
        )
        return rmb, budget, [], [], []
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - metadata is required for safe resume
        print_progress(settings, budget, f"resume skipped: checkpoint metadata unreadable: {exc}")
        return rmb, budget, [], [], []
    if not _checkpoint_q_values_match(meta, settings):
        saved_q_values = meta.get("settings", {}).get("q_values")
        print_progress(
            settings,
            budget,
            "resume skipped: checkpoint q_values "
            f"{saved_q_values!r} do not match current q_values {_settings_q_values(settings)!r}",
        )
        return rmb, budget, [], [], []

    try:
        loaded = RMB.load(base_path, rng=rmb.rng)
        rmb._data = loaded._data
    except Exception as exc:  # noqa: BLE001 - resume is best-effort
        print_progress(settings, budget, f"resume skipped: {exc}")
        return rmb, budget, [], [], []

    submitted = [config for config, _ in measured_items(rmb._data)]
    batch_history: list[dict] = []
    live_volume_history: list[LiveVolumePoint] = []
    try:
        spent = float(meta.get("spent_hqc", 0.0))
        budget = Budget(
            remaining_hqc=max(settings.hqc_budget - spent, 0.0),
            spent_hqc=spent,
            jobs=int(meta.get("jobs", 0)),
            max_job_circuits=int(meta.get("max_job_circuits", 0)),
        )
        raw_batch_history = meta.get("batch_history", [])
        if isinstance(raw_batch_history, list):
            batch_history = raw_batch_history
            if not batch_history:
                batch_history = _legacy_batch_history_from_data(rmb._data)
            restored_submitted = _submitted_from_batch_history(batch_history)
            if restored_submitted:
                submitted = restored_submitted
        live_volume_history = _live_volume_history_from_meta(meta)
    except Exception as exc:  # noqa: BLE001 - metadata passed q check but may be partial
        print_progress(settings, budget, f"checkpoint metadata ignored: {exc}")

    print_progress(
        settings,
        budget,
        f"resumed {len(submitted)} measured configs from {base_path}",
    )
    return rmb, budget, submitted, batch_history, live_volume_history


# --------------------------------------------------------------------------- #
# Main run loop
# --------------------------------------------------------------------------- #
def run_with_budget(
    settings: CostAwareSurfaceSettings,
) -> tuple[RMB, list[RMBConfig], Budget]:
    """Run the cost-aware global surface designer."""
    rng, rmb, budget = start_run(settings)
    (
        rmb,
        budget,
        submitted,
        batch_history,
        volume_history,
    ) = _resume_measurement_checkpoint(rmb, settings, budget)
    data = rmb._data
    backend = rmb.backend

    prior_centre, prior_std = _prior_moments(settings)
    grid_centre = prior_centre.copy()
    grid_halfwidth = settings.grid_halfwidth_sigmas * prior_std

    score_grid = _score_grid(settings)
    _add_prior_live_volume_point(
        settings, budget, grid_centre, grid_halfwidth, volume_history
    )
    prev_log_rms = float("inf")
    stalled = 0

    start_iteration = min(max(budget.jobs, 0), settings.max_iterations)
    if _measured_config_count(data) >= settings.posterior_min_configs:
        _maybe_update_live_surface_plot(rmb, settings, budget, start_iteration)
        _maybe_update_live_volume_plot(
            rmb, settings, budget, start_iteration, volume_history
        )

    for iteration in range(start_iteration, settings.max_iterations):
        remaining_global = budget.remaining_hqc
        if remaining_global <= 0.0:
            break

        params, log_prior, t = _build_grid(grid_centre, grid_halfwidth, settings)
        weights, ess = _grid_posterior(data, params, log_prior, settings)

        n_measured = _measured_config_count(data)
        log_rms = float("inf")
        if n_measured >= settings.posterior_min_configs:
            mean_var, _ = _surface_log_uncertainty(params, weights, score_grid, settings)
            log_rms = float(np.sqrt(mean_var))
            diag = _residual_diagnostics(data, params, weights, settings)
            max_z, n_flag = diag["max_abs"], diag["n_flag"]
            if settings.verbose:
                print_progress(
                    settings,
                    budget,
                    f"iter {iteration}: configs={n_measured} ess={ess:.0f}/{len(params)} "
                    f"log-rms={log_rms:.4f} chi2/dof={diag['chi2_dof']:.2f} "
                    f"dev/dof={diag['dev_dof']:.2f} max|z|={max_z:.1f} flagged={n_flag} "
                    f"resid<z> low-r={diag['lowr_mean']:+.2f}(n={diag['lowr_n']}) "
                    f"high-r={diag['highr_mean']:+.2f}(n={diag['highr_n']}) "
                    f"obs-tail={diag['tail_n']}",
                )
            if log_rms <= settings.target_log_rms:
                print_progress(settings, budget, "target surface accuracy reached")
                break
            if np.isfinite(prev_log_rms):
                rel_improvement = (prev_log_rms - log_rms) / max(prev_log_rms, _EPS)
                stalled = stalled + 1 if rel_improvement < settings.min_rel_improvement else 0
            prev_log_rms = log_rms
            if stalled >= settings.stop_patience:
                print_progress(settings, budget, "stopping: information plateau (marginal value low)")
                break

        affordable_cap = min(settings.max_cost_per_run, remaining_global)
        if affordable_cap <= 0.0:
            break

        early = iteration < settings.warmup_iterations
        measured_slices = {
            (float(c.ratio_2_qb_gates), int(c.n_qubits)) for c, _ in measured_items(data)
        }
        requests = _assemble_batch(
            params, weights, settings, affordable_cap, measured_slices, early
        )

        requests = [req for req in requests if req.shots > 0]
        if not requests:
            print_progress(settings, budget, "no affordable informative probe; stopping")
            break

        if settings.verbose and int(settings.max_qubit_window or 0) > 0:
            batch_qs = sorted({int(req.config.n_qubits) for req in requests})
            if batch_qs:
                print_progress(
                    settings,
                    budget,
                    f"batch Q-window {batch_qs[0]}..{batch_qs[-1]} "
                    f"(register W={batch_qs[-1]}, max idle {batch_qs[-1] - batch_qs[0]}) "
                    f"over Q={batch_qs}",
                )

        # Charge the batch at the physical register width W = max Q in the batch
        # (H2 runs one submission on one register; idle qubits are paid for).
        cost = batch_hqc_cost_at_physical_width(requests)
        if cost > remaining_global + 1e-9:
            print_progress(settings, budget, "next batch exceeds remaining budget; stopping")
            break

        first_shot_indices = {
            req.config: data.get(req.config, backend.default_estimator()).num_runs()
            for req in requests
        }
        outcomes = spend_request_batch(
            backend, rng, data, requests, seed=settings.rng_seed
        )
        budget.spend_batch(cost, sum(req.shots for req in requests))
        submitted.extend(req.config for req in requests)
        batch_history.append(
            _batch_checkpoint_record(
                batch_number=budget.jobs,
                iteration=iteration + 1,
                cost_hqc=cost,
                requests=requests,
                outcomes=outcomes,
                first_shot_indices=first_shot_indices,
            )
        )
        _maybe_update_live_surface_plot(rmb, settings, budget, iteration + 1)
        _maybe_update_live_volume_plot(
            rmb, settings, budget, iteration + 1, volume_history
        )
        _save_measurement_checkpoint(
            rmb,
            settings,
            budget,
            submitted,
            batch_history,
            volume_history,
            iteration=iteration + 1,
            reason="after_batch",
        )

        # Re-centre the grid on the running posterior for the next iteration.
        grid_centre, grid_halfwidth = _axis_layout(t, weights, settings, prior_std)

    # --- final report ------------------------------------------------------
    params, log_prior, t_final = _build_grid(grid_centre, grid_halfwidth, settings)
    weights, ess = _grid_posterior(data, params, log_prior, settings)
    _report(rmb, settings, budget, params, weights, score_grid)

    # --- continuous MAP re-fit on all gathered data (grid-initialised) -----
    # Runs once at the end.  Optimises the same free parameters as the grid over
    # a continuous space (so nu2 / the boundary land between grid ticks) and
    # reports the goodness-of-fit at the true optimum -- this is the chi2/dof to
    # trust in preference to the grid's quantised value.
    continuous_t_hat = None
    continuous_result = None
    if settings.continuous_refit:
        try:
            t_init = np.sum(weights[:, None] * t_final, axis=0)
            continuous_t_hat, continuous_result = _continuous_map_fit(
                data, settings, t_init
            )
            if continuous_t_hat is not None:
                _report_continuous_fit(
                    data, settings, continuous_t_hat, params, weights, continuous_result
                )
            else:
                print_progress(
                    settings,
                    budget,
                    "continuous re-fit skipped: no free axes or no measured data",
                )
        except Exception as exc:  # noqa: BLE001 - re-fit is a best-effort readout
            print_progress(settings, budget, f"continuous re-fit skipped: {exc}")
            continuous_t_hat = None

    _record_volume_point_from_posterior(
        params,
        weights,
        settings,
        budget.jobs,
        volume_history,
    )
    _save_measurement_checkpoint(
        rmb,
        settings,
        budget,
        submitted,
        batch_history,
        volume_history,
        iteration=budget.jobs,
        reason="final",
    )

    if settings.save_path is not None:
        try:
            from sympleq.applications.randomized_benchmarking.experiments.common import (
                save_crossings,
            )
            base_path = save_crossings(rmb, settings, budget, submitted)
        except Exception as exc:  # noqa: BLE001 - saving is best-effort
            print_progress(settings, budget, f"save skipped: {exc}")
            base_path = None
    else:
        base_path = None

    if settings.plot:
        try:
            if (
                settings.surface_uncertainty_plot
                and (
                    settings.surface_uncertainty_plot_path is not None
                    or settings.surface_uncertainty_plot_show
                    or settings.surface_plot_show
                )
            ):
                mesh = _surface_uncertainty_plot_mesh(params, weights, settings)
                if mesh is None:
                    print_progress(
                        settings,
                        budget,
                        "uncertainty surface plot skipped: no finite surface in plot bounds",
                    )
                else:
                    show_uncertainty = (
                        settings.surface_uncertainty_plot_show
                        or settings.surface_plot_show
                    )
                    png_path = plot_boundary_uncertainty_surface(
                        *mesh,
                        _measured_arrays(data),
                        settings,
                        png_path=settings.surface_uncertainty_plot_path,
                        show=show_uncertainty,
                        show_block=False,
                        show_pause=settings.live_surface_plot_pause,
                        close=not show_uncertainty,
                        figure_name="cost-aware final uncertainty surface",
                    )
                    if png_path is not None:
                        print_progress(
                            settings,
                            budget,
                            f"uncertainty surface plot -> {png_path}",
                        )
        except Exception as exc:  # noqa: BLE001 - plotting is best-effort
            print_progress(settings, budget, f"uncertainty surface plot skipped: {exc}")
        # Grid-median vs continuous-MAP comparison surface.
        try:
            if (
                settings.continuous_refit
                and continuous_t_hat is not None
                and (
                    settings.grid_vs_continuous_plot_path is not None
                    or settings.surface_plot_show
                    or settings.surface_uncertainty_plot_show
                )
            ):
                grid_med = np.asarray(
                    [
                        float(_wquantile(params[:, j], weights, 0.5))
                        for j in range(_N_PARAMS)
                    ]
                )
                cont_params = _params_from_t(
                    np.asarray(continuous_t_hat, dtype=float)[None, :]
                )[0]
                show_cmp = (
                    settings.surface_plot_show
                    or settings.surface_uncertainty_plot_show
                )
                png_path = plot_grid_vs_continuous_surface(
                    grid_med,
                    cont_params,
                    _measured_arrays(data),
                    settings,
                    png_path=settings.grid_vs_continuous_plot_path,
                    show=show_cmp,
                    show_block=False,
                    show_pause=settings.live_surface_plot_pause,
                    close=not show_cmp,
                    figure_name="grid vs continuous surface",
                )
                if png_path is not None:
                    print_progress(
                        settings,
                        budget,
                        f"grid-vs-continuous surface plot -> {png_path}",
                    )
        except Exception as exc:  # noqa: BLE001 - plotting is best-effort
            print_progress(
                settings, budget, f"grid-vs-continuous plot skipped: {exc}"
            )
        try:
            if not settings.surface_plot_show:
                import matplotlib
                matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from sympleq.applications.randomized_benchmarking.experiments.plots import (
                plot_uncertainty_diagnostics,
            )
            fixed_q_path = (
                None
                if base_path is None
                else base_path.parent / f"{base_path.stem}_uncertainty.png"
            )
            axes = plot_uncertainty_diagnostics(
                data,
                settings,
                png_path=fixed_q_path,
                show=settings.surface_plot_show,
            )
            if axes and fixed_q_path is not None:
                print_progress(
                    settings,
                    budget,
                    f"fixed-Q threshold plot -> {fixed_q_path}",
                )
            if axes and not settings.surface_plot_show:
                plt.close(axes[0].figure)
        except Exception as exc:  # noqa: BLE001 - plotting is best-effort
            print_progress(settings, budget, f"plot skipped: {exc}")
        try:
            if settings.live_volume_plot_path is not None or settings.live_volume_plot_show:
                png_path = plot_live_volume_history(
                    volume_history,
                    settings,
                    png_path=settings.live_volume_plot_path,
                    show=settings.live_volume_plot_show,
                    show_block=False,
                    show_pause=settings.live_volume_plot_pause,
                    close=not settings.live_volume_plot_show,
                    figure_name="cost-aware final volume",
                )
                if png_path is not None:
                    print_progress(settings, budget, f"volume plot -> {png_path}")
        except Exception as exc:  # noqa: BLE001 - plotting is best-effort
            print_progress(settings, budget, f"volume plot skipped: {exc}")

    return rmb, submitted, budget


def _report(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    params: np.ndarray,
    weights: np.ndarray,
    score_grid: list[tuple[float, int]],
) -> None:
    print_experiment_summary(rmb._data, settings, budget)
    mean_var, per_node = _surface_log_uncertainty(params, weights, score_grid, settings)
    print("\nSurface posterior")
    print(f"  measured configs: {_measured_config_count(rmb._data)}")
    print(f"  log-rms uncertainty over scored (r,Q): {np.sqrt(mean_var):.4f}")
    diag = _residual_diagnostics(rmb._data, params, weights, settings)
    print(
        "  goodness of fit (grid): "
        f"chi2/dof={diag['chi2_dof']:.3f} deviance/dof={diag['dev_dof']:.3f} "
        f"(N={diag['n_points']}, k={_effective_free_params(settings)}, "
        f"dof={diag['dof']}, max|z|={diag['max_abs']:.1f}, flagged={diag['n_flag']})"
    )
    a = [float(_wquantile(params[:, j], weights, 0.5)) for j in range(_N_PARAMS)]
    print(f"  median rates @Qref (L1, L2; dL/dQ: m1, m2; d2L/dQ2: nu1, nu2; V) = "
          f"({a[_I_L1]:.3g}, {a[_I_L2]:.3g}; {a[_I_M1]:.3g}, {a[_I_M2]:.3g}; "
          f"{a[_I_N1]:.3g}, {a[_I_N2]:.3g}; {a[_I_V]:.3g})")
    surface_scores = _surface_s1_s2(params, weights, settings, _analytic_lindblad_gates)
    print(
        "  analytic Lindblad surface S1/S2: "
        f"{surface_scores['surface_S1']:.5g} / {surface_scores['surface_S2']:.5g}"
    )
    print(
        "  mean |delta log gates| / sigma(log gates): "
        f"{surface_scores['surface_mean_delta_log_gates']:.5g} / "
        f"{surface_scores['surface_mean_sigma_log_gates']:.5g}"
    )
    print(
        "  analytic Lindblad success-side log-volume fit/reference "
        "(log10 gates x ratio x Q): "
        f"{surface_scores['surface_volume_fit']:.5g} / "
        f"{surface_scores['surface_volume_reference']:.5g} "
        f"(ratio {surface_scores['surface_volume_ratio']:.5g})"
    )
    _print_surface_log_ratio_grid(
        params,
        weights,
        settings,
        _analytic_lindblad_gates,
        "analytic Lindblad",
    )

    if settings.calibrated_surface_path is not None:
        calibrated_label = _calibrated_surface_label(settings)
        calibrated_scores = _surface_s1_s2(
            params,
            weights,
            settings,
            _calibrated_surface_gates,
        )
        print(
            f"  {calibrated_label} surface S1/S2: "
            f"{calibrated_scores['surface_S1']:.5g} / "
            f"{calibrated_scores['surface_S2']:.5g}"
        )
        print(
            f"  {calibrated_label} success-side log-volume fit/reference "
            "(log10 gates x ratio x Q): "
            f"{calibrated_scores['surface_volume_fit']:.5g} / "
            f"{calibrated_scores['surface_volume_reference']:.5g} "
            f"(ratio {calibrated_scores['surface_volume_ratio']:.5g})"
        )
        _print_surface_log_ratio_grid(
            params,
            weights,
            settings,
            _calibrated_surface_gates,
            calibrated_label,
        )

    if settings.gp_grid_surface_path is not None:
        gp_label = _gp_grid_surface_label(settings)
        gp_scores = _surface_s1_s2(
            params,
            weights,
            settings,
            _gp_grid_surface_gates,
        )
        print(
            f"  {gp_label} surface S1/S2: "
            f"{gp_scores['surface_S1']:.5g} / "
            f"{gp_scores['surface_S2']:.5g}"
        )
        print(
            f"  {gp_label} success-side log-volume fit/reference "
            "(log10 gates x ratio x Q): "
            f"{gp_scores['surface_volume_fit']:.5g} / "
            f"{gp_scores['surface_volume_reference']:.5g} "
            f"(ratio {gp_scores['surface_volume_ratio']:.5g})"
        )
        _print_surface_log_ratio_grid(
            params,
            weights,
            settings,
            _gp_grid_surface_gates,
            gp_label,
        )

    if settings.truth_boundary is not None:
        sq_log = []
        for r, q in score_grid:
            d = _denominator(params, r, q, settings)
            ok = (d > 0.0) & (weights > 0.0)
            if not np.any(ok):
                continue
            n_hat = float(np.sum((weights[ok] / np.sum(weights[ok])) * (1.0 / d[ok])))
            n_true = float(settings.truth_boundary(r, q))
            if n_hat > 0 and n_true > 0:
                sq_log.append((np.log(n_hat / n_true)) ** 2)
        if sq_log:
            score = float(np.exp(-np.sqrt(np.mean(sq_log))))
            print(f"  surface score (vs truth): {score:.3f}")


def run(settings: CostAwareSurfaceSettings) -> tuple[RMB, list[RMBConfig]]:
    rmb, submitted, _ = run_with_budget(settings)
    return rmb, submitted


def _surface_log_ratio_grid(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    reference_gates,
) -> tuple[np.ndarray, list[tuple[int, list[float]]]]:
    """Rows of log(n_fit/n_reference) using the volume convention."""
    r_lo, r_hi = settings.ratio_bounds
    ratios = np.linspace(r_lo, r_hi, max(settings.score_ratio_points, 2))
    rows: list[tuple[int, list[float]]] = []
    for q_raw in settings.q_values:
        q = int(q_raw)
        row = []
        reference = reference_gates(
            ratios,
            np.full_like(ratios, float(q)),
            settings,
        )
        for r, n_true in zip(ratios, reference):
            d = _denominator(params, float(r), float(q), settings)
            ok = (d > 0.0) & (weights > 0.0)
            value = float("nan")
            if np.any(ok) and np.isfinite(n_true) and n_true > 0.0:
                w = weights[ok] / np.sum(weights[ok])
                mean_d = float(np.sum(w * d[ok]))
                if mean_d > 0.0:
                    n_fit = 1.0 / mean_d
                    value = float(np.log(n_fit / n_true))
            row.append(value)
        rows.append((q, row))
    return ratios, rows


def _print_surface_log_ratio_grid(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    reference_gates,
    label: str,
) -> None:
    ratios, rows = _surface_log_ratio_grid(params, weights, settings, reference_gates)
    print(f"  log(n_fit/n_reference) grid ({label}; fit uses 1/mean D)")
    print("    r      " + " ".join(f"{r:>+6.2f}" for r in ratios))
    for q, row in rows:
        print(f"    Q={q:<3d} " + " ".join(f"{x:+.2f}" for x in row))


def _surface_s1_s2(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    reference_gates,
) -> dict[str, float | int]:
    """3-D extension of S1/S2 over the scored (r,Q) surface."""
    r_lo, r_hi = settings.ratio_bounds
    ratios = np.linspace(r_lo, r_hi, max(settings.score_ratio_points, 2))
    qubits = np.asarray(settings.q_values, dtype=float)
    rr, qq = np.meshgrid(ratios, qubits)
    med = np.asarray([float(_wquantile(params[:, j], weights, 0.5)) for j in range(_N_PARAMS)])
    _qref = _q_reference(settings)
    _dq = qq - _qref
    _lam1 = med[_I_L1] + med[_I_M1] * _dq + med[_I_N1] * _dq * _dq
    _lam2 = med[_I_L2] + med[_I_M2] * _dq + med[_I_N2] * _dq * _dq
    d_med = (_lam1 * (1.0 - rr) + _lam2 * rr) / _LN2
    fitted_log = np.where(d_med > 0.0, -np.log(d_med), np.nan)
    reference = reference_gates(rr, qq, settings)
    reference_log = np.where(reference > 0.0, np.log(reference), np.nan)

    sigma = np.full_like(fitted_log, np.nan, dtype=float)
    mean_d_gates = np.full_like(fitted_log, np.nan, dtype=float)
    for q_index, q in enumerate(qubits):
        for r_index, r in enumerate(ratios):
            d = _denominator(params, r, q, settings)
            ok = (d > 0.0) & (weights > 0.0)
            if not np.any(ok):
                continue
            w = weights[ok] / np.sum(weights[ok])
            mean_d = float(np.sum(w * d[ok]))
            if mean_d > 0.0:
                mean_d_gates[q_index, r_index] = 1.0 / mean_d
            log_n = -np.log(d[ok])
            mean = float(np.sum(w * log_n))
            var = float(np.sum(w * log_n * log_n) - mean * mean)
            sigma[q_index, r_index] = float(np.sqrt(max(var, 0.0)))

    valid = (
        np.isfinite(fitted_log)
        & np.isfinite(reference_log)
        & np.isfinite(sigma)
    )
    if np.count_nonzero(valid) < 2:
        return {
            "surface_S1": float("nan"),
            "surface_S2": float("nan"),
            "surface_volume_fit": float("nan"),
            "surface_volume_lower_1sigma": float("nan"),
            "surface_volume_upper_1sigma": float("nan"),
            "surface_volume_reference": float("nan"),
            "surface_volume_ratio": float("nan"),
            "surface_mean_delta_log_gates": float("nan"),
            "surface_mean_sigma_log_gates": float("nan"),
            "surface_score_points": int(np.count_nonzero(valid)),
        }

    fitted_log = np.where(valid, fitted_log, np.nan)
    delta = np.where(valid, np.abs(fitted_log - reference_log), np.nan)
    sigma = np.where(valid, sigma, np.nan)
    r_span = max(r_hi - r_lo, _EPS)
    q_span = float(max(qubits) - min(qubits))

    def _surface_average(values: np.ndarray) -> float:
        values = np.nan_to_num(values)
        if len(qubits) == 1:
            return float(_trapezoid(values[0], ratios) / r_span)
        return float(
            _trapezoid(_trapezoid(values, ratios, axis=1), qubits)
            / max(r_span * q_span, _EPS)
        )

    normaliser = abs(_surface_average(fitted_log))
    if not np.isfinite(normaliser) or normaliser <= _EPS:
        normaliser = float(np.nanmean(np.abs(fitted_log)))
    normaliser = max(float(normaliser), _EPS)
    s1 = float(_surface_average(delta) / normaliser)
    s2 = float(_surface_average(sigma) / normaliser)

    def _volume_axes() -> tuple[np.ndarray, np.ndarray, float, float]:
        if (
            reference_gates is _gp_grid_surface_gates
            and settings.gp_grid_surface_path is not None
        ):
            gp_surface = _load_gp_grid_surface(str(Path(settings.gp_grid_surface_path)))
            gate_axis = np.asarray(gp_surface["gates"], dtype=float)
            return (
                np.asarray(gp_surface["ratios"], dtype=float),
                np.asarray(gp_surface["qubits"], dtype=float),
                float(np.nanmin(gate_axis)),
                float(np.nanmax(gate_axis)),
            )
        n_lo, n_hi = settings.n_gates_bounds
        return ratios, qubits, float(n_lo), float(n_hi)

    def _success_side_volume(
        boundary_gates: np.ndarray,
        ratio_axis: np.ndarray,
        qubit_axis: np.ndarray,
        min_gates: float,
        max_gates: float,
    ) -> float:
        valid_boundary = np.isfinite(boundary_gates) & (boundary_gates > 0.0)
        if not np.any(valid_boundary):
            return float("nan")
        min_gates = max(float(min_gates), _EPS)
        max_gates = max(float(max_gates), min_gates)
        clipped = np.clip(boundary_gates, min_gates, max_gates)
        height = np.where(
            valid_boundary,
            np.maximum(np.log10(clipped) - np.log10(min_gates), 0.0),
            np.nan,
        )
        height = np.nan_to_num(height)
        if len(qubit_axis) == 1:
            return float(_trapezoid(height[0], ratio_axis))
        return float(_trapezoid(_trapezoid(height, ratio_axis, axis=1), qubit_axis))

    volume_ratios, volume_qubits, volume_min_gates, volume_max_gates = _volume_axes()
    v_rr, v_qq = np.meshgrid(volume_ratios, volume_qubits)
    volume_mean_d_gates = np.full(v_rr.shape, np.nan, dtype=float)
    volume_sigma = np.full(v_rr.shape, np.nan, dtype=float)
    for q_index, q in enumerate(volume_qubits):
        for r_index, r in enumerate(volume_ratios):
            d = _denominator(params, float(r), float(q), settings)
            ok = (d > 0.0) & (weights > 0.0)
            if not np.any(ok):
                continue
            w = weights[ok] / np.sum(weights[ok])
            mean_d = float(np.sum(w * d[ok]))
            if mean_d > 0.0:
                volume_mean_d_gates[q_index, r_index] = 1.0 / mean_d
            log_n = -np.log(d[ok])
            mean = float(np.sum(w * log_n))
            var = float(np.sum(w * log_n * log_n) - mean * mean)
            volume_sigma[q_index, r_index] = float(np.sqrt(max(var, 0.0)))
    volume_reference_gates = reference_gates(v_rr, v_qq, settings)

    valid_volume = (
        np.isfinite(volume_mean_d_gates)
        & np.isfinite(volume_reference_gates)
        & np.isfinite(volume_sigma)
        & (volume_mean_d_gates > 0.0)
        & (volume_reference_gates > 0.0)
    )
    volume_mean_d_gates = np.where(valid_volume, volume_mean_d_gates, np.nan)
    volume_sigma = np.where(valid_volume, volume_sigma, np.nan)
    volume_reference_gates = np.where(valid_volume, volume_reference_gates, np.nan)
    volume_fit = _success_side_volume(
        volume_mean_d_gates,
        volume_ratios,
        volume_qubits,
        volume_min_gates,
        volume_max_gates,
    )
    volume_lower = _success_side_volume(
        volume_mean_d_gates * np.exp(-volume_sigma),
        volume_ratios,
        volume_qubits,
        volume_min_gates,
        volume_max_gates,
    )
    volume_upper = _success_side_volume(
        volume_mean_d_gates * np.exp(volume_sigma),
        volume_ratios,
        volume_qubits,
        volume_min_gates,
        volume_max_gates,
    )
    volume_reference = _success_side_volume(
        volume_reference_gates,
        volume_ratios,
        volume_qubits,
        volume_min_gates,
        volume_max_gates,
    )
    volume_ratio = (
        float(volume_fit / volume_reference)
        if (
            np.isfinite(volume_fit)
            and np.isfinite(volume_reference)
            and abs(volume_reference) > _EPS
        )
        else float("nan")
    )
    return {
        "surface_S1": s1,
        "surface_S2": s2,
        "surface_volume_fit": float(volume_fit),
        "surface_volume_lower_1sigma": float(volume_lower),
        "surface_volume_upper_1sigma": float(volume_upper),
        "surface_volume_reference": float(volume_reference),
        "surface_volume_ratio": volume_ratio,
        "surface_mean_delta_log_gates": float(np.nanmean(delta)),
        "surface_mean_sigma_log_gates": float(np.nanmean(sigma)),
        "surface_score_points": int(np.count_nonzero(valid)),
    }


def run_surface(
    *,
    q_values: tuple[int, ...] = tuple(range(5, 51, 5)),
    hqc_budget: float = 1500.0,
    rng_seed: int | None = 0,
    backend_model: str = "sympleq",
    plot: bool = True,
    surface_plot_path: str | Path | None = DEFAULT_SURFACE_PLOT_PATH,
    surface_plot_show: bool = False,
    surface_plot_analytic: bool = True,
    surface_plot_n_gates_bounds: tuple[int, int] | None = None,
    surface_uncertainty_plot: bool = True,
    surface_uncertainty_plot_path: str | Path | None = DEFAULT_UNCERTAINTY_PLOT_PATH,
    surface_uncertainty_plot_show: bool = False,
    surface_uncertainty_sigma: float = 1.0,
    calibrated_surface_path: str | Path | None = None,
    gp_grid_surface_path: str | Path | None = None,
    gp_grid_surface_label: str | None = None,
    live_surface_plot: bool = False,
    live_surface_plot_path: str | Path | None = DEFAULT_LIVE_SURFACE_PLOT_PATH,
    live_surface_plot_every: int = 1,
    live_surface_plot_show: bool = False,
    live_surface_plot_pause: float = 0.25,
    live_volume_plot: bool = False,
    live_volume_plot_path: str | Path | None = DEFAULT_LIVE_VOLUME_PLOT_PATH,
    live_volume_plot_every: int = 1,
    live_volume_plot_show: bool = False,
    live_volume_plot_pause: float = 0.25,
    **settings_overrides,
) -> tuple[RMB, list[RMBConfig], Budget]:
    """Run a full multi-Q surface fit and (optionally) write the 3-D plot."""
    hqc_budget = float(settings_overrides.pop("hqc_budget", hqc_budget))
    backend_model = str(settings_overrides.pop("backend_model", backend_model))
    settings = CostAwareSurfaceSettings(
        q_values=tuple(q_values),
        n_qubits=int(round(float(np.median(q_values)))),
        hqc_budget=hqc_budget,
        rng_seed=rng_seed,
        backend_model=backend_model,
        plot=plot,
        surface_plot_path=surface_plot_path,
        surface_plot_show=surface_plot_show,
        surface_plot_analytic=surface_plot_analytic,
        surface_plot_n_gates_bounds=surface_plot_n_gates_bounds,
        surface_uncertainty_plot=surface_uncertainty_plot or surface_uncertainty_plot_show,
        surface_uncertainty_plot_path=surface_uncertainty_plot_path,
        surface_uncertainty_plot_show=surface_uncertainty_plot_show,
        surface_uncertainty_sigma=surface_uncertainty_sigma,
        calibrated_surface_path=calibrated_surface_path,
        gp_grid_surface_path=gp_grid_surface_path,
        gp_grid_surface_label=gp_grid_surface_label,
        live_surface_plot=live_surface_plot or live_surface_plot_show,
        live_surface_plot_path=live_surface_plot_path,
        live_surface_plot_every=live_surface_plot_every,
        live_surface_plot_show=live_surface_plot_show,
        live_surface_plot_pause=live_surface_plot_pause,
        live_volume_plot=live_volume_plot or live_volume_plot_show,
        live_volume_plot_path=live_volume_plot_path,
        live_volume_plot_every=live_volume_plot_every,
        live_volume_plot_show=live_volume_plot_show,
        live_volume_plot_pause=live_volume_plot_pause,
        **settings_overrides,
    )
    return run_with_budget(settings)


def print_run_score(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
) -> None:
    """Original-style fit-to-analytic-line score for the reference-Q slice."""
    print(f"\nRun score (n_qubits = {settings.n_qubits})")
    try:
        score = parametric_boundary_fit_score(rmb._data, settings)
        coverage50, coverage90 = parametric_bootstrap_analytic_coverage(
            rmb._data,
            settings,
            seed=None if settings.rng_seed is None else settings.rng_seed + 500_000,
        )
        print(f"  score: {score:.3f}")
        print(f"  analytic line inside bootstrap 50% band: {100.0 * coverage50:.1f}%")
        print(f"  analytic line inside bootstrap 90% band: {100.0 * coverage90:.1f}%")
    except Exception as exc:  # noqa: BLE001 - report rather than crash the summary
        print(f"  score unavailable: {exc}")
    print(f"  spent: {budget.spent_hqc:.1f} HQC")
    print(f"  measured configs: {_measured_config_count(rmb._data)}")


if __name__ == "__main__":
    # --- surface run with quadratic (nu2) Q-dependence freed ----------------
    # grid_resolution / boundary_fit_resolution free the two-qubit curvature axis
    # nu2 (index 5): the linear model is nested inside the quadratic one, so with
    # the mean-zero shrinkage prior nu2 stays near 0 if the rate is linear and
    # moves off 0 only if there is real curvature.  nu1 (index 4) and V (index 6)
    # stay pinned.  The continuous MAP re-fit then reports nu2 and the
    # goodness-of-fit at the true (unquantised) optimum, and a grid-vs-continuous
    # comparison surface is written alongside the usual plots.
    rmb, configs, budget = run_surface(
        q_values=tuple(range(20, 51, 1)),   # score over Q = 20, 21, ..., 50
        acquisition_q_resolution=12,     # Q candidates across the full Q range
        acquisition_ratio_points=15,     # r candidates across ratio_bounds
        backend_model="sympleq",
        hqc_budget=500.0,
        plot=True,
        surface_plot_show=True,
        surface_plot_n_gates_bounds=(100, 3000),
        gp_grid_surface_path=(
            EXPERIMENTS_DIR.parent
            / "rmb_data"
            / "FLE_20260702_063751_gp_grid_3d.npz"
        ),
        gp_grid_surface_label="Rick/Shreya Grid",
        rng_seed=123,
        verbose=True,
        use_scrambler=True,
        live_surface_plot_show=True,
        live_surface_plot_pause=0.5,
        live_volume_plot=True,
        live_volume_plot_show=True,
        live_volume_plot_pause=0.5,
        max_qubit_window=1,
        min_distinct_q_coverage=10,
        # Free the two-qubit quadratic curvature nu2 (grid index 5) in both the
        # design-loop grid and the finer scoring grid; nu1 and V stay pinned.
        grid_resolution=(11, 11, 7, 7, 1, 5, 1),
        boundary_fit_resolution=(17, 17, 9, 9, 1, 7, 1),
        # Continuous MAP re-fit + grid-vs-continuous comparison plot (defaults on).
        continuous_refit=True,
    )
