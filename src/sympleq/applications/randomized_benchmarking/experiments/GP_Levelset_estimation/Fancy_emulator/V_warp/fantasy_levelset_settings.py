"""Settings handles for noisy Fancy-emulator native-V wrapper FLE runs."""

from __future__ import annotations

from pathlib import Path

from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.Fancy_emulator.V_warp.unstitched_quantinuum import (
    vwarp_fancy_h21e,
    vwarp_fancy_h22e,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.fantasy_levelset_settings import (
    FantasySettings,
)


# =============================================================================
# CONTROL PANEL
# =============================================================================
# Edit this section for noisy Fancy-emulator native-V wrapper runs.

# -------------------------------------------------------------------------
# TARGET HANDLE
# -------------------------------------------------------------------------

TARGET_THRESHOLD = 0.5

# -------------------------------------------------------------------------
# NUMBER OF QUBITS
# -------------------------------------------------------------------------

N_QUBITS = 26

# -------------------------------------------------------------------------
# SEARCH-BOX HANDLES
# -------------------------------------------------------------------------

N_GATES_BOUNDS = (200, 1500)
RATIO_BOUNDS = (0.1, 0.7)

# -------------------------------------------------------------------------
# HQC BUDGET HANDLES
# -------------------------------------------------------------------------

HQC_BUDGET = 10000
MAX_COST_PER_RUN = 100

# -------------------------------------------------------------------------
# GATE BUDGET HANDLES
# -------------------------------------------------------------------------

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
INITIAL_SOBOL_MAX_COST_PER_RUN = 100
SOBOL_SCRAMBLE = True

SAVE_GP_PREDICTION_GRID = True
SAVE_REAL_CHECKPOINTS = True

# -------------------------------------------------------------------------
# RECOVERY HANDLES
# -------------------------------------------------------------------------

RECOVERY_MODE = True
RECOVERY_FOLDER = Path(
    r"Personal\FLE\Fancy_emulator\V_warp\H2_1E\q26\seed_42\FLE_V_warp_20260827_121746"
)

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

BACKEND_FACTORY = vwarp_fancy_h21e
# BACKEND_FACTORY = vwarp_fancy_h22e

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
