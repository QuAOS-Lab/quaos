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
from sympleq.applications.randomized_benchmarking.backends.exponential import ExponentialBackend
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    Budget,
    CrossingSettings,
    batch_hqc_cost,
    measured_items,
    print_experiment_summary,
    print_progress,
    spend_request_batch,
    start_run,
)
from sympleq.applications.randomized_benchmarking.experiments.plots import (
    REFERENCE_OFFSET,
    REFERENCE_SLOPE,
    parametric_boundary_fit_score,
    parametric_bootstrap_analytic_coverage,
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


def _with_backend_seed_suffix(
    path: str | Path | None,
    *,
    backend_model: str,
    rng_seed: int | None,
) -> str | Path | None:
    """Append backend/seed to output file names for reproducible local runs."""
    if path is None or rng_seed is None:
        return path
    if backend_model not in {"sympleq", "exponential"}:
        return path
    path_type = Path if isinstance(path, Path) else str
    p = Path(path)
    suffix = f"_{backend_model}_seed{rng_seed}"
    if p.stem.endswith(suffix):
        return path
    new_path = p.with_name(f"{p.stem}{suffix}{p.suffix}")
    return new_path if path_type is Path else str(new_path)


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
    # If enabled, refresh the 3-D surface plot after batches during the run.
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
    surface_uncertainty_plot: bool = False
    surface_uncertainty_plot_path: str | Path | None = DEFAULT_UNCERTAINTY_PLOT_PATH
    surface_uncertainty_plot_show: bool = False
    surface_uncertainty_sigma: float = 1.0
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
    initial_one_q_pauli_error: float = BASE_1Q_PAULI_ERROR * 0.6
    initial_two_q_pauli_error: float = BASE_2Q_PAULI_ERROR * 0.6
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
            ):
                object.__setattr__(
                    self,
                    attr,
                    _with_backend_seed_suffix(
                        getattr(self, attr),
                        backend_model=self.backend_model,
                        rng_seed=self.rng_seed,
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
    model = settings.asymptote_model
    if isinstance(model, str):
        if model == "depolarizing":
            return 2.0 ** (-float(q))
        if model == "zero":
            return 0.0
        raise ValueError(f"Unknown asymptote_model {model!r}.")
    return float(model)


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
    # Cap probe depth at the tighter of two scales:
    #   (i)  a fraction of the deepest measurable gate count (global guard), and
    #   (ii) a multiple of *this slice's* predicted boundary n* (local guard).
    # The local guard is the important one for the low-r / high-Q corner: there
    # n* is large and the broad-posterior boundary quantiles below would
    # otherwise reach far past the crossing, where survival has saturated to the
    # floor and a measurement is uninformative (flat curve, negligible Fisher
    # information) yet, because low-r gates are cheap, still wins on BALD-per-HQC.
    # Tying the cap to the local n* removes those runaway probes at their cause
    # without suppressing genuine near-boundary measurements anywhere.
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

    # Predicted-survival floor: drop candidate depths where the model already
    # expects near-certain failure (posterior-mean p below the floor).  These
    # deep-tail shots are uninformative about the crossing yet leverage the
    # asymptote, so they skew the shared rate fit (notably overestimating n* at
    # low r).  Applied in p-space so it self-adjusts across (r, Q).
    floor = settings.min_predicted_survival
    if floor > 0.0 and candidates:
        kept = []
        for g in candidates:
            pij = _link_probability(params, g, r, q, settings)
            if float(np.sum(weights * pij)) >= floor:
                kept.append(g)
        # Never return empty purely because of the floor: if every candidate is
        # below it (e.g. a very small-B(Q) high-Q slice), keep the shallowest so
        # the slice still contributes its most-informative available depth.
        candidates = kept if kept else candidates[:1]

    # Predicted-fidelity cap: drop candidate depths that sit too high on the
    # decay curve (posterior-mean renormalised F = 2^{-nD} above the cap).  These
    # shallow points are where the single-exponential link is least valid -- the
    # circuit is not yet twirled, so an early-depth measurement carries the
    # finite-depth transient (apparent F(0) > 1) that biases the fitted slope.
    # F is independent of B(Q) and V, so the cap self-adjusts across (r, Q).
    fmax = settings.max_predicted_fidelity
    if 0.0 < fmax < 1.0 and candidates:
        kept = []
        for g in candidates:
            fid = np.power(2.0, -g * np.maximum(d, 0.0))   # renormalised F per row
            if float(np.sum(weights * fid)) <= fmax:
                kept.append(g)
        # If every candidate is above the cap (all too shallow, e.g. a very
        # large-n* slice whose boundary lies beyond the depth cap), keep the
        # deepest -- the lowest-F, most-decayed available probe.
        candidates = kept if kept else candidates[-1:]
    return candidates


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
    """Return L and L^2 where L[g, k] = log n*(node g) for parameter row k.

    Entries with non-positive denominator are set to 0; the corresponding
    parameter rows always have zero posterior weight, so they never contribute.
    """
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


def _assemble_batch(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    affordable_cap: float,
    measured_slices: set[tuple[float, int]],
    early: bool,
) -> list[MeasurementRequest]:
    """Greedily pack one stitched submission to maximise value per HQC.

    The per-shot value is the parameter mutual information (BALD).  Each (r, Q)
    slice has diminishing returns (rho^slice_shots), which spreads shots across
    slices.  Marginal HQC of every increment is the exact difference of
    batch_hqc_cost, so the per-submission base cost is amortised correctly: the
    first increment pays the base, later increments pay only marginal
    bare/reset/measurement cost.

    Candidates span the joint (Q, r, n) design (``_acquisition_slices`` x the
    per-slice candidate depths), so the design chooses register size as freely as
    ratio and depth.  A distinct-Q coverage floor keeps interior Q sampled until
    the Q axis is spanned, so a free-Q design cannot collapse onto the two Q
    extremes before any nonlinearity in the rate could be detected.
    """
    acq_slices = _acquisition_slices(settings)
    measured_q = {q for (_, q) in measured_slices}
    need_q_coverage = len(measured_q) < int(settings.min_distinct_q_coverage)

    # Build the candidate pool: a few near-boundary n per candidate (r, Q) slice.
    candidates: list[dict] = []
    for r, q in acq_slices:
        for n in _candidate_n_values(params, weights, settings, r, q):
            value = _bald_per_shot(params, weights, n, r, q, settings)
            if value <= 0.0:
                continue
            # Coverage bonus: during warmup for any unmeasured (r, Q) slice, and
            # (independently) for any unmeasured Q while the distinct-Q floor is
            # unmet -- the latter is what prevents endpoint collapse under a
            # free-Q design.
            slice_new = (r, q) not in measured_slices
            q_new = q not in measured_q
            bonus = 1.0
            if early and slice_new:
                bonus = settings.coverage_bonus
            if need_q_coverage and q_new:
                bonus = max(bonus, settings.coverage_bonus)
            candidates.append(
                {"n": n, "r": r, "q": q, "slice": (r, q), "value": value * bonus}
            )
    if not candidates:
        return []

    rho = settings.shot_diminish_rho
    shots: dict[tuple[int, float, int], int] = {}   # (n, r, q) -> shots
    slice_shots: dict[tuple[float, int], int] = {}

    def requests_from(shots_map: dict) -> list[MeasurementRequest]:
        return [
            MeasurementRequest(make_config_q(settings, n, r, q), s)
            for (n, r, q), s in shots_map.items()
            if s > 0
        ]

    current_cost = 0.0
    while True:
        best = None
        best_density = 0.0
        for cand in candidates:
            key = (cand["n"], cand["r"], cand["q"])
            if shots.get(key, 0) >= settings.max_shots_per_probe:
                continue
            marginal_value = cand["value"] * (rho ** slice_shots.get(cand["slice"], 0))
            trial = dict(shots)
            trial[key] = trial.get(key, 0) + 1
            trial_cost = batch_hqc_cost(requests_from(trial))
            if trial_cost > affordable_cap:
                continue
            marginal_cost = trial_cost - current_cost
            if marginal_cost <= 0.0:
                marginal_cost = _EPS
            density = marginal_value / marginal_cost
            if density > best_density:
                best_density = density
                best = (key, cand["slice"], trial_cost)
        if best is None:
            break
        key, slc, trial_cost = best
        shots[key] = shots.get(key, 0) + 1
        slice_shots[slc] = slice_shots.get(slc, 0) + 1
        current_cost = trial_cost

    return requests_from(shots)


# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #
def _residual_diagnostics(
    data,
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> dict:
    """Standardised-residual diagnostics over measured points (truth-free).

    Returns a dict with the global ``max_abs`` |z| and ``n_flag`` (count of
    |z| > 3), plus a breakdown of the *signed* mean residual by ratio region and
    an empirical-tail count.  The signed means are the informative part:

    * ``lowr_mean`` strongly negative  ->  at low r the model predicts higher
      survival than observed, i.e. it has placed the boundary *too deep*
      (overestimated n*).  This is the fingerprint of the low-r skew, and it is
      invisible to any predicted-p test because the model believes those probes
      sit near p = 0.5.
    * ``tail_n`` counts points whose *observed* survival p_hat is near the
      asymptote (< ``tail_diagnostic_pmax``): empirically dead-tail shots,
      regardless of what the model predicted for them.

    Goodness-of-fit: two aggregate statistics summarise whether the surface model
    is *adequate* (not merely precise).  ``chi2_dof`` is the Pearson
    chi^2 = sum z_i^2 over (N - k) dof; ``dev_dof`` is the Bernoulli
    deviance / (N - k), which is the correct likelihood-ratio analogue and stays
    valid near p = 0 or 1 where the Gaussian z approximation degrades.  Both ~1
    mean the model fits within shot noise; >> 1 flags systematic misfit (the
    z-region breakdown then says *where*).

    All quantities use only measured counts and the current posterior, so they
    read identically on the simulator and on hardware.
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
        tail_n += int(p_hat < pmax)        # empirical tail, keyed on observed survival
        # Pearson chi^2 contribution (one Bernoulli cell per measured config).
        chi2 += z * z
        # Bernoulli deviance contribution: 2[ s log(p_hat/p_bar)
        #   + f log((1-p_hat)/(1-p_bar)) ], with 0 log 0 = 0.
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


def _maybe_update_live_surface_plot(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    iteration: int,
) -> None:
    """Best-effort in-run refresh of the 3-D surface diagnostic plot."""
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
        if settings.live_surface_plot_show:
            import matplotlib.pyplot as plt
            plt.ion()
        png_path = plot_surface_3d(
            rmb,
            settings,
            png_path=settings.live_surface_plot_path,
            show=settings.live_surface_plot_show,
            show_block=False,
            show_pause=settings.live_surface_plot_pause,
            close=not settings.live_surface_plot_show,
            figure_name="cost-aware live surface",
        )
        if png_path is not None:
            print_progress(settings, budget, f"live surface plot -> {png_path}")
    except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
        print_progress(settings, budget, f"live surface plot skipped: {exc}")


def _maybe_update_live_volume_plot(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    iteration: int,
    history: list[tuple[int, float, float]],
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
        scores = _surface_s1_s2(params, weights, settings, _analytic_lindblad_gates)
        fitted_volume = float(scores["surface_volume_fit"])
        true_volume = _reference_surface_volume(settings, _analytic_lindblad_gates)
        if not (np.isfinite(fitted_volume) and np.isfinite(true_volume)):
            return
        history.append((iteration, fitted_volume, true_volume))
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
                f"live volume plot -> {png_path} (ratio {ratio:.4g})",
            )
    except Exception as exc:  # noqa: BLE001 - live plotting should never kill a run
        print_progress(settings, budget, f"live volume plot skipped: {exc}")


def _add_prior_live_volume_point(
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    grid_centre: np.ndarray,
    grid_halfwidth: np.ndarray,
    history: list[tuple[int, float, float]],
) -> None:
    """Seed the live volume curve with the prior predictive volume at iteration 0."""
    if not settings.live_volume_plot:
        return
    if settings.live_volume_plot_path is None and not settings.live_volume_plot_show:
        return
    try:
        params, log_prior, _ = _build_grid(grid_centre, grid_halfwidth, settings)
        # With no measurements this is just the normalized prior over the grid.
        log_w = log_prior - np.max(log_prior)
        weights = np.exp(log_w)
        weights = weights / np.sum(weights)
        scores = _surface_s1_s2(params, weights, settings, _analytic_lindblad_gates)
        fitted_volume = float(scores["surface_volume_fit"])
        true_volume = _reference_surface_volume(settings, _analytic_lindblad_gates)
        if not (np.isfinite(fitted_volume) and np.isfinite(true_volume)):
            return
        history.append((0, fitted_volume, true_volume))
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


def _checkpoint_meta_path(save_path: str | Path) -> Path:
    base_path = resolve_data_path(save_path)
    return base_path.parent / f"{base_path.stem}_checkpoint.json"


def _save_measurement_checkpoint(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
    submitted: list[RMBConfig],
    *,
    iteration: int | None = None,
    reason: str = "checkpoint",
) -> None:
    """Save RMB-format measurements plus a small budget sidecar."""
    if settings.save_path is None or not settings.checkpoint_after_batch:
        return
    try:
        rmb.save(settings.save_path)
        base_path = resolve_data_path(settings.save_path)
        meta = {
            "reason": reason,
            "iteration": iteration,
            "hqc_budget": settings.hqc_budget,
            "spent_hqc": budget.spent_hqc,
            "remaining_hqc": budget.remaining_hqc,
            "jobs": budget.jobs,
            "max_job_circuits": budget.max_job_circuits,
            "submitted_configs": len(submitted),
            "measured_configs": _measured_config_count(rmb._data),
        }
        _checkpoint_meta_path(settings.save_path).write_text(
            json.dumps(meta, indent=2),
            encoding="utf-8",
        )
        print_progress(settings, budget, f"saved measurements -> {base_path}")
    except Exception as exc:  # noqa: BLE001 - checkpointing should not kill a run
        print_progress(settings, budget, f"checkpoint save skipped: {exc}")


def _resume_measurement_checkpoint(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    budget: Budget,
) -> tuple[RMB, Budget, list[RMBConfig]]:
    """Load saved measurements into the current backend and restore budget state."""
    if settings.save_path is None or not settings.resume_from_save:
        return rmb, budget, []
    base_path = resolve_data_path(settings.save_path)
    if not base_path.exists():
        return rmb, budget, []

    try:
        loaded = RMB.load(base_path, rng=rmb.rng)
        # Keep the currently configured backend, but reuse the saved estimators.
        rmb._data = loaded._data
    except Exception as exc:  # noqa: BLE001 - resume is best-effort
        print_progress(settings, budget, f"resume skipped: {exc}")
        return rmb, budget, []

    submitted = [config for config, _ in measured_items(rmb._data)]
    meta_path = _checkpoint_meta_path(settings.save_path)
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            spent = float(meta.get("spent_hqc", 0.0))
            budget = Budget(
                remaining_hqc=max(settings.hqc_budget - spent, 0.0),
                spent_hqc=spent,
                jobs=int(meta.get("jobs", 0)),
                max_job_circuits=int(meta.get("max_job_circuits", 0)),
            )
        except Exception as exc:  # noqa: BLE001 - metadata is useful but optional
            print_progress(settings, budget, f"checkpoint metadata ignored: {exc}")
    else:
        print_progress(
            settings,
            budget,
            "checkpoint metadata missing; saved measurements loaded as a warm start",
        )

    print_progress(
        settings,
        budget,
        f"resumed {len(submitted)} measured configs from {base_path}",
    )
    return rmb, budget, submitted


# --------------------------------------------------------------------------- #
# Main run loop
# --------------------------------------------------------------------------- #
def run_with_budget(
    settings: CostAwareSurfaceSettings,
) -> tuple[RMB, list[RMBConfig], Budget]:
    """Run the cost-aware global surface designer."""
    rng, rmb, budget = start_run(settings)
    rmb, budget, submitted = _resume_measurement_checkpoint(rmb, settings, budget)
    data = rmb._data
    backend = rmb.backend

    prior_centre, prior_std = _prior_moments(settings)
    grid_centre = prior_centre.copy()
    grid_halfwidth = settings.grid_halfwidth_sigmas * prior_std

    score_grid = _score_grid(settings)
    volume_history: list[tuple[int, float, float]] = []
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

        cost = batch_hqc_cost(requests)
        if cost > remaining_global + 1e-9:
            print_progress(settings, budget, "next batch exceeds remaining budget; stopping")
            break

        spend_request_batch(backend, rng, data, requests, seed=settings.rng_seed)
        budget.spend_batch(cost, sum(req.shots for req in requests))
        submitted.extend(req.config for req in requests)
        _maybe_update_live_surface_plot(rmb, settings, budget, iteration + 1)
        _maybe_update_live_volume_plot(
            rmb, settings, budget, iteration + 1, volume_history
        )
        _save_measurement_checkpoint(
            rmb,
            settings,
            budget,
            submitted,
            iteration=iteration + 1,
            reason="after_batch",
        )

        # Re-centre the grid on the running posterior for the next iteration.
        grid_centre, grid_halfwidth = _axis_layout(t, weights, settings, prior_std)

    # --- final report ------------------------------------------------------
    params, log_prior, _ = _build_grid(grid_centre, grid_halfwidth, settings)
    weights, ess = _grid_posterior(data, params, log_prior, settings)
    _report(rmb, settings, budget, params, weights, score_grid)
    _save_measurement_checkpoint(
        rmb,
        settings,
        budget,
        submitted,
        iteration=budget.jobs,
        reason="final",
    )

    if settings.save_path is not None:
        try:
            from sympleq.applications.randomized_benchmarking.experiments.common import (
                save_crossings,
            )
            save_crossings(rmb, settings, budget, submitted)
        except Exception as exc:  # noqa: BLE001 - saving is best-effort
            print_progress(settings, budget, f"save skipped: {exc}")

    if settings.plot:
        try:
            from sympleq.applications.randomized_benchmarking.experiments.plots import (
                plot_crossing_results,
            )
            plot_crossing_results(data, settings, submitted, base_path=None)
        except Exception as exc:  # noqa: BLE001 - plotting is best-effort
            print_progress(settings, budget, f"plot skipped: {exc}")
        if len(settings.q_values) > 1 and settings.surface_plot_path is not None:
            try:
                plot_surface_3d(
                    rmb,
                    settings,
                    png_path=settings.surface_plot_path,
                    show=settings.surface_plot_show,
                )
                print_progress(settings, budget, f"surface plot -> {settings.surface_plot_path}")
            except Exception as exc:  # noqa: BLE001 - plotting is best-effort
                print_progress(settings, budget, f"surface plot skipped: {exc}")
        if (
            len(settings.q_values) > 1
            and settings.surface_uncertainty_plot
            and (
                settings.surface_uncertainty_plot_path is not None
                or settings.surface_uncertainty_plot_show
            )
        ):
            try:
                unc_path = plot_surface_uncertainty_3d(
                    rmb,
                    settings,
                    png_path=settings.surface_uncertainty_plot_path,
                    show=settings.surface_uncertainty_plot_show,
                )
                if unc_path is not None:
                    print_progress(
                        settings, budget, f"uncertainty surface plot -> {unc_path}"
                    )
            except Exception as exc:  # noqa: BLE001 - plotting is best-effort
                print_progress(
                    settings, budget, f"uncertainty surface plot skipped: {exc}"
                )

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
        "  goodness of fit: "
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
        "  analytic Lindblad log10-gate surface volume fit/reference "
        "(fit uses log10(1/mean D)): "
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
            f"  {calibrated_label} log10-gate surface volume fit/reference "
            "(fit uses log10(1/mean D)): "
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


def _surface_median_params(
    data,
    settings: CostAwareSurfaceSettings,
) -> np.ndarray | None:
    """Weighted-median surface parameters (L1, L2, m1, m2, V), or None."""
    params, weights = _stateless_posterior(data, settings)
    if params is None:
        return None
    return np.array(
        [float(_wquantile(params[:, j], weights, 0.5)) for j in range(_N_PARAMS)]
    )


def _surface_logn_mean_sigma(
    data,
    settings: CostAwareSurfaceSettings,
    rr: np.ndarray,
    qq: np.ndarray,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Posterior mean and std of log n*(r,Q) on a (rr, qq) mesh.

    Returns (mean_logn, sigma_logn) with the same shape as ``rr``, computed from
    the full grid posterior.  Memory-light: loops mesh points and reduces over
    the posterior per point, storing only per-node scalars.  The std is the
    pointwise posterior uncertainty of the natural-log boundary -- the quantity
    the +-sigma surfaces envelope.
    """
    params, weights = _stateless_posterior(data, settings)
    if params is None or weights is None:
        return None
    mean_logn = np.full(rr.shape, np.nan, dtype=float)
    sigma_logn = np.full(rr.shape, np.nan, dtype=float)
    flat_r = rr.ravel()
    flat_q = qq.ravel()
    mean_flat = mean_logn.ravel()
    sigma_flat = sigma_logn.ravel()
    for idx in range(flat_r.size):
        d = _denominator(params, float(flat_r[idx]), float(flat_q[idx]), settings)
        ok = (d > 0.0) & (weights > 0.0)
        if not np.any(ok):
            continue
        w = weights[ok] / np.sum(weights[ok])
        log_n = -np.log(d[ok])
        mean = float(np.sum(w * log_n))
        var = float(np.sum(w * log_n * log_n) - mean * mean)
        mean_flat[idx] = mean
        sigma_flat[idx] = float(np.sqrt(max(var, 0.0)))
    return mean_flat.reshape(rr.shape), sigma_flat.reshape(rr.shape)


def _trapezoid(y: np.ndarray, x: np.ndarray, *, axis: int = -1) -> np.ndarray:
    trapezoid = getattr(np, "trapezoid", None)
    if trapezoid is not None:
        return trapezoid(y, x, axis=axis)
    return np.trapz(y, x, axis=axis)


def _noise_scales(settings: CostAwareSurfaceSettings) -> tuple[float, float]:
    one_q = float(getattr(settings, "one_q_noise_scale", 1.0))
    two_q = float(getattr(settings, "two_q_noise_scale", 1.0))
    return one_q, two_q


def _analytic_lindblad_gates(
    ratios: np.ndarray,
    qubits: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> np.ndarray:
    """Analytic Lindblad fidelity-0.5 surface.

    The current analytic reference is independent of Q because the reference
    model is expressed in total gate count and two-qubit ratio only.
    """
    one_q_scale, two_q_scale = _noise_scales(settings)
    ratio_grid, qubit_grid = np.broadcast_arrays(ratios, qubits)
    gates = np.log(2.0) / (
        REFERENCE_OFFSET * one_q_scale
        + REFERENCE_SLOPE * two_q_scale * ratio_grid
    )
    return np.full_like(qubit_grid, 1.0, dtype=float) * gates


def _reference_surface_volume(
    settings: CostAwareSurfaceSettings,
    reference_gates=_analytic_lindblad_gates,
) -> float:
    """Integral of reference log10(n*(r,Q)) over the scored rectangle."""
    r_lo, r_hi = settings.ratio_bounds
    ratios = np.linspace(r_lo, r_hi, max(settings.score_ratio_points, 2))
    qubits = np.asarray(settings.q_values, dtype=float)
    rr, qq = np.meshgrid(ratios, qubits)
    values = reference_gates(rr, qq, settings)
    finite = np.isfinite(values)
    if not np.any(finite):
        return float("nan")
    # Some optional reference grids cover only part of the scored rectangle.
    # Integrate over the supported rectangular subgrid instead of silently
    # replacing unsupported cells by zero.
    if not np.all(finite):
        row_keep = np.any(finite, axis=1)
        col_keep = np.any(finite, axis=0)
        if np.count_nonzero(row_keep) == 0 or np.count_nonzero(col_keep) < 2:
            return float("nan")
        qubits = qubits[row_keep]
        ratios = ratios[col_keep]
        values = values[np.ix_(row_keep, col_keep)]
        if not np.all(np.isfinite(values)):
            return float("nan")
    if np.any(values <= 0.0):
        return float("nan")
    values = np.log10(values)
    if len(qubits) == 1:
        return float(_trapezoid(values[0], ratios))
    return float(_trapezoid(_trapezoid(values, ratios, axis=1), qubits))


@lru_cache(maxsize=8)
def _load_calibrated_surface(path: str) -> dict:
    with open(path) as f:
        payload = json.load(f)
    for key in ("lambda1_poly_coefficients", "lambda2_poly_coefficients"):
        if key not in payload:
            raise ValueError(f"Calibrated surface {path!r} is missing {key!r}.")
    return payload


def _calibrated_surface_label(settings: CostAwareSurfaceSettings) -> str:
    if settings.calibrated_surface_path is None:
        return "calibrated"
    payload = _load_calibrated_surface(str(Path(settings.calibrated_surface_path)))
    seed = payload.get("settings", {}).get("rng_seed")
    suffix = f" seed={seed}" if seed is not None else ""
    return f"{payload.get('label', 'calibrated SympleQ')}{suffix}"


def _calibrated_surface_gates(
    ratios: np.ndarray,
    qubits: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> np.ndarray:
    """Evaluate the saved calibrated effective-rate reference surface.

    The calibration file stores natural-log rates lambda_i(Q).  The boundary is

        n_ref(r,Q) = ln 2 / [(1-r) lambda_1(Q) + r lambda_2(Q)].
    """
    if settings.calibrated_surface_path is None:
        raise ValueError("No calibrated_surface_path configured.")
    payload = _load_calibrated_surface(str(Path(settings.calibrated_surface_path)))
    ratio_grid, qubit_grid = np.broadcast_arrays(ratios, qubits)
    variable = payload.get("polynomial_variable", "q_minus_reference")
    if variable == "q_minus_reference":
        x = qubit_grid.astype(float) - float(payload.get("q_reference", 0.0))
    elif variable == "q":
        x = qubit_grid.astype(float)
    else:
        raise ValueError(f"Unknown calibrated polynomial variable {variable!r}.")
    lam1 = np.polyval(np.asarray(payload["lambda1_poly_coefficients"], dtype=float), x)
    lam2 = np.polyval(np.asarray(payload["lambda2_poly_coefficients"], dtype=float), x)
    denominator = (1.0 - ratio_grid) * lam1 + ratio_grid * lam2
    gates = np.full_like(ratio_grid, np.nan, dtype=float)
    ok = denominator > 0.0
    gates[ok] = _LN2 / denominator[ok]
    return gates


def _gp_grid_surface_label(settings: CostAwareSurfaceSettings) -> str:
    if settings.gp_grid_surface_label is not None:
        return settings.gp_grid_surface_label
    if settings.gp_grid_surface_path is None:
        return "external GP grid"
    return Path(settings.gp_grid_surface_path).stem


def _gate_crossing_from_probability(
    gates: np.ndarray,
    probabilities: np.ndarray,
    target: float,
) -> float:
    """Linear p=target crossing along increasing gate counts."""
    gates = np.asarray(gates, dtype=float)
    probabilities = np.asarray(probabilities, dtype=float)
    ok = np.isfinite(gates) & np.isfinite(probabilities)
    if np.count_nonzero(ok) < 2:
        return float("nan")
    gates = gates[ok]
    probabilities = probabilities[ok]
    diff = probabilities - target
    exact = np.flatnonzero(np.isclose(diff, 0.0, atol=1e-12))
    if len(exact):
        return float(gates[int(exact[0])])
    sign_change = np.flatnonzero(diff[:-1] * diff[1:] < 0.0)
    if not len(sign_change):
        return float("nan")
    # Prefer the physical high-to-low survival transition if it exists.
    down = [i for i in sign_change if diff[i] > 0.0 and diff[i + 1] < 0.0]
    i = int(down[0] if down else sign_change[0])
    p0, p1 = probabilities[i], probabilities[i + 1]
    if abs(p1 - p0) <= _EPS:
        return float(0.5 * (gates[i] + gates[i + 1]))
    frac = (target - p0) / (p1 - p0)
    return float(gates[i] + frac * (gates[i + 1] - gates[i]))


@lru_cache(maxsize=8)
def _load_gp_grid_surface(path: str) -> dict:
    with np.load(path, allow_pickle=True) as payload:
        qubits = np.asarray(payload["qubits_axis"], dtype=float)
        ratios = np.asarray(payload["ratio_axis"], dtype=float)
        gates = np.asarray(payload["gates_axis"], dtype=float)
        probabilities = np.asarray(payload["probabilities"], dtype=float)
        target = float(payload["target"]) if "target" in payload else 0.5

    expected = (len(qubits), len(ratios), len(gates))
    if probabilities.shape != expected:
        raise ValueError(
            f"GP grid probabilities shape {probabilities.shape} does not match "
            f"(qubits, ratios, gates) = {expected}."
        )
    contour = np.full((len(qubits), len(ratios)), np.nan, dtype=float)
    for qi in range(len(qubits)):
        for ri in range(len(ratios)):
            contour[qi, ri] = _gate_crossing_from_probability(
                gates, probabilities[qi, ri, :], target
            )
    return {
        "qubits": qubits,
        "ratios": ratios,
        "gates": gates,
        "contour": contour,
        "target": target,
    }


def _interp_grid2d_nan(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    values: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
) -> np.ndarray:
    """Bilinear interpolation on a regular grid, returning nan outside/near gaps."""
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)
    values = np.asarray(values, dtype=float)
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    out = np.full_like(x, np.nan, dtype=float)
    inside = (
        (x >= x_axis[0])
        & (x <= x_axis[-1])
        & (y >= y_axis[0])
        & (y <= y_axis[-1])
    )
    if not np.any(inside):
        return out
    xi = np.searchsorted(x_axis, x[inside], side="right") - 1
    yi = np.searchsorted(y_axis, y[inside], side="right") - 1
    xi = np.clip(xi, 0, len(x_axis) - 2)
    yi = np.clip(yi, 0, len(y_axis) - 2)
    x0, x1 = x_axis[xi], x_axis[xi + 1]
    y0, y1 = y_axis[yi], y_axis[yi + 1]
    tx = np.divide(x[inside] - x0, x1 - x0, out=np.zeros_like(x0), where=x1 != x0)
    ty = np.divide(y[inside] - y0, y1 - y0, out=np.zeros_like(y0), where=y1 != y0)
    v00 = values[xi, yi]
    v10 = values[xi + 1, yi]
    v01 = values[xi, yi + 1]
    v11 = values[xi + 1, yi + 1]
    valid = np.isfinite(v00) & np.isfinite(v10) & np.isfinite(v01) & np.isfinite(v11)
    interp = (
        (1.0 - tx) * (1.0 - ty) * v00
        + tx * (1.0 - ty) * v10
        + (1.0 - tx) * ty * v01
        + tx * ty * v11
    )
    inside_idx = np.flatnonzero(inside)
    out.flat[inside_idx[valid]] = interp[valid]
    return out


def _gp_grid_surface_gates(
    ratios: np.ndarray,
    qubits: np.ndarray,
    settings: CostAwareSurfaceSettings,
) -> np.ndarray:
    """Evaluate the p=target contour extracted from an external GP grid."""
    if settings.gp_grid_surface_path is None:
        raise ValueError("No gp_grid_surface_path configured.")
    surface = _load_gp_grid_surface(str(Path(settings.gp_grid_surface_path)))
    return _interp_grid2d_nan(
        surface["qubits"],
        surface["ratios"],
        surface["contour"],
        np.asarray(qubits, dtype=float),
        np.asarray(ratios, dtype=float),
    )


def _surface_log_ratio_grid(
    params: np.ndarray,
    weights: np.ndarray,
    settings: CostAwareSurfaceSettings,
    reference_gates,
) -> tuple[np.ndarray, list[tuple[int, list[float]]]]:
    """Rows of log(n_fit/n_reference) using the volume convention.

    The fitted contour here is exactly the same estimator used by the printed
    volume ratio: n_fit(r,Q) = 1 / E[D(r,Q)].
    """
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
    """3-D extension of S1/S2 over the scored (r,Q) surface.

    S1 integrates absolute log-boundary error over r and Q and normalises by
    the fitted log-boundary volume.  S2 integrates posterior standard deviation
    of log n*(r,Q) over the same surface and uses the same normalisation.
    """
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

    def _surface_integral(values: np.ndarray) -> float:
        values = np.nan_to_num(values)
        if len(qubits) == 1:
            return float(_trapezoid(values[0], ratios))
        return float(_trapezoid(_trapezoid(values, ratios, axis=1), qubits))

    # Fill invalid cells as zero for integration.  The score grid is regular,
    # and invalid cells are rare unless the posterior gives nonphysical rates.
    normaliser = abs(_surface_average(fitted_log))
    if not np.isfinite(normaliser) or normaliser <= _EPS:
        normaliser = float(np.nanmean(np.abs(fitted_log)))
    normaliser = max(float(normaliser), _EPS)
    s1 = float(_surface_average(delta) / normaliser)
    s2 = float(_surface_average(sigma) / normaliser)
    valid_volume = (
        np.isfinite(mean_d_gates)
        & np.isfinite(reference)
        & (mean_d_gates > 0.0)
        & (reference > 0.0)
    )
    volume_fit = _surface_integral(
        np.where(valid_volume, np.log10(mean_d_gates), np.nan)
    )
    volume_reference = _surface_integral(
        np.where(valid_volume, np.log10(reference), np.nan)
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
        "surface_volume_reference": float(volume_reference),
        "surface_volume_ratio": volume_ratio,
        "surface_mean_delta_log_gates": float(np.nanmean(delta)),
        "surface_mean_sigma_log_gates": float(np.nanmean(sigma)),
        "surface_score_points": int(np.count_nonzero(valid)),
    }


def plot_live_volume_history(
    history: list[tuple[int, float, float]],
    settings: CostAwareSurfaceSettings,
    png_path: str | Path | None = None,
    *,
    show: bool = False,
    show_block: bool = True,
    show_pause: float = 0.001,
    close: bool = True,
    figure_name: str | None = None,
):
    """Plot fitted log10-gate surface volume history with the known-rate reference line."""
    if not history:
        return None
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    iterations = np.asarray([row[0] for row in history], dtype=int)
    fitted = np.asarray([row[1] for row in history], dtype=float)
    reference = np.asarray([row[2] for row in history], dtype=float)
    ref = float(reference[-1])
    ratio = fitted / ref if ref else np.full_like(fitted, np.nan)

    fig = plt.figure(num=figure_name, figsize=(8.5, 5.2), clear=True)
    ax = fig.add_subplot(111)
    ax.plot(iterations, fitted, marker="o", linewidth=1.8, label="posterior fit")
    ax.axhline(
        ref,
        color="crimson",
        linestyle="--",
        linewidth=1.6,
        label="known-rate exponential reference",
    )
    if settings.gp_grid_surface_path is not None:
        gp_ref = _reference_surface_volume(settings, _gp_grid_surface_gates)
        if np.isfinite(gp_ref):
            ax.axhline(
                gp_ref,
                color="black",
                linestyle="--",
                linewidth=1.6,
                label=_gp_grid_surface_label(settings),
            )
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\int \log_{10} n_*(r,Q)\,dr\,dQ$")
    ax.set_title(r"Fidelity-0.5 surface volume in $\log_{10}(n)$")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")

    ax_ratio = ax.twinx()
    ax_ratio.plot(
        iterations,
        ratio,
        color="0.25",
        linestyle=":",
        marker=".",
        alpha=0.8,
        label="fit/reference",
    )
    ax_ratio.set_ylabel("fit / reference")
    ax_ratio.axhline(1.0, color="0.25", linestyle=":", linewidth=1.0, alpha=0.6)
    if ref:
        left_lo, left_hi = ax.get_ylim()
        ax_ratio.set_ylim(left_lo / ref, left_hi / ref)

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


def plot_surface_3d(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    png_path: str | Path | None = None,
    *,
    show: bool = False,
    show_block: bool = True,
    show_pause: float = 0.001,
    close: bool = True,
    figure_name: str | None = None,
    n_surface_draws: int = 0,
):
    """Plot the fitted n*(r, Q) boundary surface with the measured points.

    The translucent surface is the posterior-median boundary
    n*(r,Q) = ln 2 / (L1(Q)(1-r) + L2(Q) r).  Measured configurations are
    scattered at their (r, Q, n) with colour set by the observed survival
    fraction (diverging about 0.5) and size by the number of shots, so points
    sitting above the surface (high fidelity) and below it (low fidelity) are
    visible.  Optionally overlays ``n_surface_draws`` faint posterior-draw
    surfaces as an uncertainty band.
    """
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

    data = rmb._data
    med = _surface_median_params(data, settings)
    if med is None:
        return None

    r_lo, r_hi = settings.ratio_bounds
    q_lo, q_hi = min(settings.q_values), max(settings.q_values)
    n_lo, n_hi = (
        settings.n_gates_bounds
        if settings.surface_plot_n_gates_bounds is None
        else settings.surface_plot_n_gates_bounds
    )
    z_lo, z_hi = np.log10([max(n_lo, _EPS), max(n_hi, _EPS)])

    grid_r = np.linspace(r_lo, r_hi, 48)
    grid_q = np.linspace(q_lo, q_hi, 48)
    rr, qq = np.meshgrid(grid_r, grid_q)

    q_ref = _q_reference(settings)

    def _surface(p):
        dq = qq - q_ref
        lam1 = p[_I_L1] + p[_I_M1] * dq + p[_I_N1] * dq * dq
        lam2 = p[_I_L2] + p[_I_M2] * dq + p[_I_N2] * dq * dq
        d = (lam1 * (1.0 - rr) + lam2 * rr) / _LN2
        z = np.where(d > 0.0, 1.0 / np.maximum(d, 1e-12), np.nan)
        return np.where((z >= n_lo) & (z <= n_hi), np.log10(z), np.nan)

    fig = plt.figure(num=figure_name, figsize=(9.5, 7.5), clear=True)
    ax = fig.add_subplot(projection="3d")

    if n_surface_draws > 0:
        params, weights = _stateless_posterior(data, settings)
        if params is not None:
            seed = 0 if settings.rng_seed is None else settings.rng_seed
            idx = np.random.default_rng(seed).choice(
                len(params), size=n_surface_draws, replace=True, p=weights
            )
            for i in idx:
                ax.plot_wireframe(
                    rr, qq, _surface(params[i]), color="0.5",
                    alpha=0.12, linewidth=0.4, rstride=6, cstride=6,
                )

    ax.plot_surface(
        rr, qq, _surface(med), cmap="viridis", alpha=0.45,
        linewidth=0, antialiased=True, rstride=1, cstride=1,
    )
    if settings.surface_plot_analytic:
        analytic_surface = _analytic_lindblad_gates(rr, qq, settings)
        analytic_surface = np.where(
            (analytic_surface >= n_lo) & (analytic_surface <= n_hi),
            np.log10(analytic_surface),
            np.nan,
        )
        ax.plot_wireframe(
            rr,
            qq,
            analytic_surface,
            color="crimson",
            linewidth=1.1,
            alpha=0.75,
            rstride=4,
            cstride=4,
            label="analytic Lindblad",
        )
    if settings.calibrated_surface_path is not None:
        calibrated_surface = _calibrated_surface_gates(rr, qq, settings)
        calibrated_surface = np.where(
            (calibrated_surface >= n_lo) & (calibrated_surface <= n_hi),
            np.log10(calibrated_surface),
            np.nan,
        )
        ax.plot_wireframe(
            rr,
            qq,
            calibrated_surface,
            color="black",
            linewidth=1.1,
            alpha=0.75,
            rstride=4,
            cstride=4,
            label=_calibrated_surface_label(settings),
        )
    if settings.gp_grid_surface_path is not None:
        gp_surface = _gp_grid_surface_gates(rr, qq, settings)
        gp_surface = np.where(
            (gp_surface >= n_lo) & (gp_surface <= n_hi),
            np.log10(gp_surface),
            np.nan,
        )
        ax.plot_wireframe(
            rr,
            qq,
            gp_surface,
            color="black",
            linewidth=1.3,
            alpha=0.85,
            rstride=4,
            cstride=4,
            label=_gp_grid_surface_label(settings),
        )

    ratios, gates, qubits, succ, fail = _measured_arrays(data)
    total = succ + fail
    mask = total > 0
    if np.any(mask):
        p_hat = succ[mask] / total[mask]
        size = 18.0 + 40.0 * np.clip(total[mask] / max(total[mask].max(), 1.0), 0, 1)
        sc = ax.scatter(
            ratios[mask], qubits[mask], np.log10(np.maximum(gates[mask], _EPS)),
            c=p_hat, cmap="coolwarm_r",
            norm=TwoSlopeNorm(vcenter=0.5, vmin=0.0, vmax=1.0),
            s=size, edgecolor="k", linewidth=0.3, depthshade=True,
        )
        cbar = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.10)
        cbar.set_label("measured survival fraction")

    ax.set_xlabel("two-qubit ratio  r")
    ax.set_ylabel("n_qubits  Q")
    ax.set_zlabel(r"$\log_{10}(\mathrm{gate\ count}\ n)$")
    ax.set_zlim(z_lo, z_hi)
    ax.set_title(r"Fidelity-0.5 boundary surface  $\log_{10} n_*(r, Q)$")
    if (
        settings.surface_plot_analytic
        or settings.calibrated_surface_path is not None
        or settings.gp_grid_surface_path is not None
    ):
        ax.legend(loc="upper left")
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


def plot_surface_uncertainty_3d(
    rmb: RMB,
    settings: CostAwareSurfaceSettings,
    png_path: str | Path | None = None,
    *,
    show: bool = False,
    show_block: bool = True,
    show_pause: float = 0.001,
    close: bool = True,
    figure_name: str | None = None,
):
    """Plot the posterior-mean boundary surface enveloped by +- k sigma surfaces.

    At every (r, Q) the grid posterior gives a full distribution over
    log n*(r,Q); ``_surface_logn_mean_sigma`` returns its mean and standard
    deviation on the plot mesh.  This draws three surfaces in log10(n):
    the mean boundary and the mean +- k sigma envelope (k =
    settings.surface_uncertainty_sigma).  The vertical gap between the upper and
    lower sheets is the spatially-resolved boundary uncertainty -- it shows
    *where* on (r, Q) the boundary is well- vs poorly-determined (typically
    widest in the sparsely-sampled, extrapolated low-r / high-Q corner).

    Note: the band is the posterior's *self-reported* uncertainty, so it is only
    honest if the posterior is calibrated.  Read it together with the
    goodness-of-fit chi2/dof printed in the report -- a tight band with a large
    chi2/dof means the boundary is precisely placed at a systematically wrong
    location (the band cannot see model misspecification).
    """
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

    data = rmb._data
    r_lo, r_hi = settings.ratio_bounds
    q_lo, q_hi = min(settings.q_values), max(settings.q_values)
    n_lo, n_hi = (
        settings.n_gates_bounds
        if settings.surface_plot_n_gates_bounds is None
        else settings.surface_plot_n_gates_bounds
    )
    z_lo, z_hi = np.log10([max(n_lo, _EPS), max(n_hi, _EPS)])
    k = float(settings.surface_uncertainty_sigma)

    # A coarser mesh than plot_surface_3d: each node reduces over the full
    # posterior, so keep it modest for speed/memory.
    grid_r = np.linspace(r_lo, r_hi, 40)
    grid_q = np.linspace(q_lo, q_hi, 40)
    rr, qq = np.meshgrid(grid_r, grid_q)

    result = _surface_logn_mean_sigma(data, settings, rr, qq)
    if result is None:
        return None
    mean_logn, sigma_logn = result

    ln10 = np.log(10.0)

    def _clip_log10(natural_log_n):
        n = np.exp(natural_log_n)
        z = natural_log_n / ln10  # log10(n)
        return np.where((n >= n_lo) & (n <= n_hi), z, np.nan)

    mean_z = _clip_log10(mean_logn)
    upper_z = _clip_log10(mean_logn + k * sigma_logn)
    lower_z = _clip_log10(mean_logn - k * sigma_logn)

    fig = plt.figure(num=figure_name, figsize=(9.5, 7.5), clear=True)
    ax = fig.add_subplot(projection="3d")

    # Mean boundary surface, coloured by the local 1-sigma width (in log10 n)
    # so the colour itself encodes where uncertainty is large.
    sigma_log10 = sigma_logn / ln10
    facenorm = plt.Normalize(
        vmin=float(np.nanmin(sigma_log10)) if np.any(np.isfinite(sigma_log10)) else 0.0,
        vmax=float(np.nanmax(sigma_log10)) if np.any(np.isfinite(sigma_log10)) else 1.0,
    )
    facecolors = plt.cm.viridis(facenorm(np.nan_to_num(sigma_log10)))
    ax.plot_surface(
        rr, qq, mean_z, facecolors=facecolors, alpha=0.85,
        linewidth=0, antialiased=True, rstride=1, cstride=1, shade=False,
    )
    # +-k sigma envelope as translucent grey sheets.
    ax.plot_surface(
        rr, qq, upper_z, color="0.5", alpha=0.18,
        linewidth=0, antialiased=True, rstride=1, cstride=1, shade=False,
    )
    ax.plot_surface(
        rr, qq, lower_z, color="0.5", alpha=0.18,
        linewidth=0, antialiased=True, rstride=1, cstride=1, shade=False,
    )
    ax.plot_wireframe(
        rr, qq, upper_z, color="0.35", alpha=0.45,
        linewidth=0.5, rstride=4, cstride=4,
    )
    ax.plot_wireframe(
        rr, qq, lower_z, color="0.35", alpha=0.45,
        linewidth=0.5, rstride=4, cstride=4,
    )

    # Overlay the measured points for context (same colour convention as the
    # main surface plot).
    ratios, gates, qubits, succ, fail = _measured_arrays(data)
    total = succ + fail
    mask = total > 0
    if np.any(mask):
        p_hat = succ[mask] / total[mask]
        ax.scatter(
            ratios[mask], qubits[mask], np.log10(np.maximum(gates[mask], _EPS)),
            c=p_hat, cmap="coolwarm_r",
            norm=TwoSlopeNorm(vcenter=0.5, vmin=0.0, vmax=1.0),
            s=12, edgecolor="k", linewidth=0.2, depthshade=True, alpha=0.6,
        )

    mappable = plt.cm.ScalarMappable(norm=facenorm, cmap="viridis")
    mappable.set_array([])
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.10)
    cbar.set_label(r"$\sigma$ of $\log_{10} n_*(r,Q)$  (boundary uncertainty)")

    ax.set_xlabel("two-qubit ratio  r")
    ax.set_ylabel("n_qubits  Q")
    ax.set_zlabel(r"$\log_{10}(\mathrm{gate\ count}\ n)$")
    ax.set_zlim(z_lo, z_hi)
    ax.set_title(
        rf"Fidelity-0.5 boundary with $\pm{k:g}\sigma$ uncertainty surfaces"
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
    surface_uncertainty_plot: bool = False,
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
    """Run a full multi-Q surface fit and (optionally) write the 3-D plot.

    Leaves the rate Q-slopes m1, m2 free (the Q-dependence) via the default grid
    resolution and pins only the visibility; the acquisition measures a sparse
    subset of Q and the linear-in-Q surface fills the range.
    """
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


def run_slope_scan(
    *,
    ratios: tuple[float, ...] = (0.1, 0.3),
    q_values: tuple[int, ...] = (5, 25, 50),
    depth_multiples: tuple[float, ...] = (0.3, 0.5, 0.7, 0.9, 1.1, 1.4, 1.7, 2.0),
    shots_per_point: int = 300,
    backend_model: str = "sympleq",
    rng_seed: int | None = 0,
    fit_f_range: tuple[float, float] = (0.08, 0.90),
    realized_ratio_tol: float = 0.05,
    **settings_overrides,
) -> dict:
    """Directly measure the RB effective decay D(r,Q) = -d(log2 F)/dn.

    This is a *diagnostic*, not a fit: it bypasses the posterior, the acquisition
    loop and the analytic reference entirely.  For each (r, Q) it runs a dense
    depth sweep, computes the renormalised fidelity

        F(n) = (p_hat(n) - B(Q)) / (V (1 - B(Q))),     p_hat = survivals / shots,

    from the *raw* survival counts, and fits a straight line log2 F = a0 - D n.
    The slope magnitude D is the effective per-gate decay; n* = 1/D.
    """
    q_values = tuple(sorted({int(q) for q in q_values}))
    settings = CostAwareSurfaceSettings(
        q_values=q_values,
        n_qubits=int(round(float(np.median(q_values)))),
        hqc_budget=float(settings_overrides.pop("hqc_budget", 1e12)),
        rng_seed=rng_seed,
        backend_model=backend_model,
        plot=False,
        save_path=None,
        verbose=False,
        **settings_overrides,
    )
    rng, rmb, budget = start_run(settings)
    data = rmb._data
    backend = rmb.backend
    vis = float(settings.initial_visibility)

    requests: list[MeasurementRequest] = []
    for r in ratios:
        n_star_analytic = float(
            _analytic_lindblad_gates(
                np.array([r]), np.array([float(min(q_values))]), settings
            ).ravel()[0]
        )
        for q in q_values:
            for u in depth_multiples:
                n_t = max(2, int(2 * round(u * n_star_analytic / 2)))
                requests.append(
                    MeasurementRequest(make_config_q(settings, n_t, float(r), int(q)), shots_per_point)
                )
    spend_request_batch(backend, rng, data, requests, seed=rng_seed)

    ratios_m, gates_m, qubits_m, succ_m, fail_m = _measured_arrays(data)

    results: dict[tuple[float, int], dict | None] = {}
    print(f"slope scan [{backend_model}]: D(r,Q) = -d(log2 F)/dn   "
          "(raw survival; no fit, no reference)")
    for r in ratios:
        print(f"  r_target={r:.2f}")
        d_by_q: list[tuple[int, float]] = []
        for q in q_values:
            b = _asymptote(settings, float(q))
            pts = []
            for i in range(len(gates_m)):
                if int(qubits_m[i]) != int(q):
                    continue
                if abs(float(ratios_m[i]) - r) > realized_ratio_tol:
                    continue
                total = succ_m[i] + fail_m[i]
                if total <= 0:
                    continue
                p_hat = succ_m[i] / total
                f = (p_hat - b) / max(vis * (1.0 - b), _EPS)
                if fit_f_range[0] <= f <= fit_f_range[1]:
                    pts.append((float(gates_m[i]), float(f), float(ratios_m[i])))
            if len(pts) < 3:
                print(f"    Q={int(q):<3d}  too few usable points ({len(pts)}) "
                      "-- target r not realizable here, or F out of range")
                results[(r, int(q))] = None
                continue
            pts.sort()
            n_arr = np.array([p[0] for p in pts])
            f_arr = np.array([p[1] for p in pts])
            r_real = float(np.mean([p[2] for p in pts]))
            y = np.log2(f_arr)
            slope, intercept = np.polyfit(n_arr, y, 1)
            d_meas = -float(slope)
            n_star_meas = (1.0 / d_meas) if d_meas > 0 else float("inf")
            yhat = slope * n_arr + intercept
            ss_res = float(np.sum((y - yhat) ** 2))
            ss_tot = float(np.sum((y - y.mean()) ** 2))
            r2 = (1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
            n_star_an = float(
                _analytic_lindblad_gates(
                    np.array([r_real]), np.array([float(q)]), settings
                ).ravel()[0]
            )
            d_analytic = 1.0 / n_star_an if n_star_an > 0 else float("nan")
            ratio = d_meas / d_analytic if d_analytic > 0 else float("nan")
            print(
                f"    Q={int(q):<3d}  D={d_meas:.4e}  n*={n_star_meas:8.1f}  "
                f"a0={intercept:+.3f}  R2={r2:.3f}  pts={len(pts)}  "
                f"r_real={r_real:.2f}  D/D_analytic={ratio:.3f}"
            )
            results[(r, int(q))] = {
                "D": d_meas, "n_star": n_star_meas, "intercept": float(intercept),
                "r2": r2, "n_points": len(pts), "r_real": r_real,
                "D_over_analytic": ratio,
            }
            d_by_q.append((int(q), d_meas))
        if len(d_by_q) >= 2:
            (q_lo, d_lo), (q_hi, d_hi) = d_by_q[0], d_by_q[-1]
            trend = d_hi / d_lo if d_lo > 0 else float("nan")
            verdict = (
                "Q-INDEPENDENT (ramp is fit/reference)" if 0.95 <= trend <= 1.05
                else "Q-DEPENDENT (fit tracks real SympleQ rate; reference is idealised)"
            )
            print(f"    -> D(Q={q_hi})/D(Q={q_lo}) = {trend:.3f}  [{verdict}]")
    return results


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
    # Single n_qubits = 5 fit: pin m1, m2 = 0 (the pure single-Q inverse boundary
    # ln 2 / (L1 (1-r) + L2 r)) and pin V = 1 (no SPAM on the simulator), so
    # only the two rate parameters are free.  Directly comparable to the original
    # experiment; reports the original-style score.
    # settings = CostAwareSurfaceSettings(
    #     n_qubits=5,
    #     q_values=(5,),
    #     grid_resolution=(25, 25, 1, 1, 1, 1, 1),
    #     boundary_fit_resolution=(41, 41, 1, 1, 1, 1, 1),
    #     hqc_budget=100.0,
    #     rng_seed=0,
    #     verbose=False,
    # )
    # rmb, _configs, budget = run_with_budget(settings)
    # print_run_score(rmb, settings, budget)

    # --- NEXT CHECK: direct slope scan -------------------------------------
    # Measures D(r,Q) = -d(log2 F)/dn straight from raw survival, with no
    # posterior, no acquisition and no analytic reference.  Read the printed
    # D(Q_hi)/D(Q_lo) trend per r: ~1.0 means the effective decay is Q-flat (so
    # the Q-ramp in the fit/analytic comparison is the fit or the Q-flat
    # reference), while a clear >1 trend means the SympleQ effective decay
    # genuinely steepens with Q and the fit is tracking real physics the
    # reference omits.  Also watch a0 (intercept): a0 != 0 flags a broken
    # unit-visibility assumption (assumption 6).  Cheap; runs many shots but no
    # budget loop.
    # run_slope_scan(
    #     ratios=(0.1, 0.3),
    #     q_values=(5, 25, 50),
    #     shots_per_point=300,
    #     backend_model="sympleq",
    #     rng_seed=123,
    # )


    # --- previous surface run (uncomment to go back to the full design) -----
    rmb, configs, budget = run_surface(
        q_values=tuple(range(5, 21, 5)),   # score over Q = 5,10,...,50
        acquisition_q_resolution=7,     # Q candidates across the full Q range
        acquisition_ratio_points=12,     # r candidates across ratio_bounds
        backend_model="sympleq",
        hqc_budget=300.0,
        plot=True,
        surface_plot_path=DEFAULT_SURFACE_PLOT_PATH,
        surface_plot_show=True,
        surface_plot_n_gates_bounds=(100, 3000),
        # New: write/show the +-1 sigma boundary-uncertainty surface at run end.
        surface_uncertainty_plot_show=True,
        surface_uncertainty_sigma=1.0,
        gp_grid_surface_path=(
            EXPERIMENTS_DIR.parent
            / "rmb_data"
            / "FLE_20260630_114233_gp_grid_3d.npz"
        ),
        gp_grid_surface_label="Rick/Shreya Surface",
        rng_seed=123,
        verbose=True,
        use_scrambler=True,
        live_surface_plot_show=True,
        live_surface_plot_pause=0.5,
        live_volume_plot=True,
        live_volume_plot_show=True,
        live_volume_plot_pause=0.5,
    )
