"""Settings handles for fantasy-batched GP level-set estimation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from sympleq.applications.randomized_benchmarking.experiments.common import (
    CrossingSettings,
    default_backend_factory,
    quantinuum_emulator_backend_factory,
    quantinuum_H2_backend_factory,
    dephasing_sympleq_backend_factory,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.Fancy_emulator.unstitched_quantinuum import (
    fancy_unstitched_h21e,
    fancy_unstitched_h22e,
)


@dataclass(frozen=True)
class FantasySettings(CrossingSettings):
    """
    User-facing handles for fantasy-batched GP level-set estimation.

    Inherited CrossingSettings handles include:
        n_gates_bounds
        ratio_bounds
        hqc_budget
        max_cost_per_run
        rng_seed
        make_config(...)
        backend_factory
    """

    # Output
    save_path: str | Path | None = None

    # Target contour
    target_threshold: float = 0.5

    # Fake anchors
    use_fake_corners: bool = True
    easy_corner_outcome: int = 1
    hard_corner_outcome: int = 0
    extra_fake_anchors: list[tuple] = field(default_factory=list)

    # Explicit Sobol warm-up
    initial_sobol_samples: int = 10
    initial_sobol_max_cost_per_run: float | None = None
    sobol_scramble: bool = False

    #GP grid saving
    save_gp_prediction_grid: bool = False
    save_real_checkpoints: bool = True

    # Gate budget per stitched submission (for emulator only)
    gate_budget: int = 7000

    # Recovery
    recovery_mode: bool = False
    recovery_folder: str | Path | None = None

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
# CONTROL PANEL
# =============================================================================
# Edit this section for ordinary FLE runs.

# -------------------------------------------------------------------------
# TARGET HANDLE
# -------------------------------------------------------------------------

TARGET_THRESHOLD = 0.5

# -------------------------------------------------------------------------
# NuUMBER of Qubits
# -------------------------------------------------------------------------

N_QUBITS = 26
# For a multi-slice run, use e.g.:
# N_QUBITS = [5,20]

# -------------------------------------------------------------------------
# SEARCH-BOX HANDLES
# -------------------------------------------------------------------------
# Set to None to use the defaults inherited from CrossingSettings.

N_GATES_BOUNDS = (100, 5000)
RATIO_BOUNDS = (0.1, 0.9)

# -------------------------------------------------------------------------
# HQC BUDGET HANDLES
# -------------------------------------------------------------------------
# Set to None to use the defaults inherited from CrossingSettings.

HQC_BUDGET = 10000
MAX_COST_PER_RUN = 100

# -------------------------------------------------------------------------
# Gate BUDGET HANDLES
# -------------------------------------------------------------------------
# Set to None to use the defaults inherited from CrossingSettings.
GATE_BUDGET = 7000

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

INITIAL_SOBOL_SAMPLES = 15
INITIAL_SOBOL_MAX_COST_PER_RUN = 30
SOBOL_SCRAMBLE = True

SAVE_GP_PREDICTION_GRID = True
SAVE_REAL_CHECKPOINTS = True

# -------------------------------------------------------------------------
# RECOVERY HANDLES
# -------------------------------------------------------------------------

RECOVERY_MODE = True
RECOVERY_FOLDER = Path(r"Personal\FLE\H2_1E\q26\seed_72\FLE_20260827_120009")


# -------------------------------------------------------------------------
# GP / AEPSYCH HANDLES
# -------------------------------------------------------------------------

OPTIMIZATION_STEPS = 500
INDUCING_SIZE = 150
ACQUISITION_FUNCTION = "GlobalSUR"
ACQUISITION_RESTARTS = 2
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

USE_GPU = False
FORCE_DEFAULT_DEVICE_DURING_AEPSYCH = True

# -------------------------------------------------------------------------
# BACKEND HANDLE
# -------------------------------------------------------------------------

# BACKEND_FACTORY = default_backend_factory
# BACKEND_FACTORY = quantinuum_emulator_backend_factory
BACKEND_FACTORY = fancy_unstitched_h21e
# BACKEND_FACTORY = fancy_unstitched_h22e
# BACKEND_FACTORY = quantinuum_H2_backend_factory
# BACKEND_FACTORY = dephasing_sympleq_backend_factory

# -------------------------------------------------------------------------
# REPRODUCIBILITY / DEBUG HANDLES
# -------------------------------------------------------------------------

RNG_SEEDS = [42]
VERBOSE_FANTASIES = True
PLOT = True


def control_panel_settings_kwargs() -> dict:
    """Build FantasySettings keyword arguments from the control panel."""
    kwargs = dict(
        target_threshold=TARGET_THRESHOLD,
        use_fake_corners=USE_FAKE_CORNERS,
        easy_corner_outcome=EASY_CORNER_OUTCOME,
        hard_corner_outcome=HARD_CORNER_OUTCOME,
        extra_fake_anchors=EXTRA_FAKE_ANCHORS,
        initial_sobol_samples=INITIAL_SOBOL_SAMPLES,
        initial_sobol_max_cost_per_run=INITIAL_SOBOL_MAX_COST_PER_RUN,
        sobol_scramble=SOBOL_SCRAMBLE,
        save_gp_prediction_grid=SAVE_GP_PREDICTION_GRID,
        save_real_checkpoints=SAVE_REAL_CHECKPOINTS,
        recovery_mode=RECOVERY_MODE,
        recovery_folder=RECOVERY_FOLDER,
        optimization_steps=OPTIMIZATION_STEPS,
        inducing_size=INDUCING_SIZE,
        acquisition_function=ACQUISITION_FUNCTION,
        acquisition_restarts=ACQUISITION_RESTARTS,
        acquisition_samples=ACQUISITION_SAMPLES,
        batching=BATCHING,
        max_batch_size=MAX_BATCH_SIZE,
        fantasy_batching=FANTASY_BATCHING,
        backend_factory=BACKEND_FACTORY,
        use_gpu=USE_GPU,
        force_default_device_during_aepsych=FORCE_DEFAULT_DEVICE_DURING_AEPSYCH,
        verbose_fantasies=VERBOSE_FANTASIES,
        plot=PLOT,
        gate_budget=GATE_BUDGET,
        n_qubits=N_QUBITS,
    )

    if N_GATES_BOUNDS is not None:
        kwargs["n_gates_bounds"] = N_GATES_BOUNDS

    if RATIO_BOUNDS is not None:
        kwargs["ratio_bounds"] = RATIO_BOUNDS

    if HQC_BUDGET is not None:
        kwargs["hqc_budget"] = HQC_BUDGET

    if MAX_COST_PER_RUN is not None:
        kwargs["max_cost_per_run"] = MAX_COST_PER_RUN

    return kwargs
