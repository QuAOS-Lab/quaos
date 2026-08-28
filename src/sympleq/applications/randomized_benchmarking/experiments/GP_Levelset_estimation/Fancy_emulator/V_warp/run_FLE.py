from __future__ import annotations

from datetime import datetime
from pathlib import Path

from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper import (
    run_FLE as h_wrapper_run,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.Fancy_emulator.V_warp import (
    fantasy_levelset_settings as settings_module,
)


def timestamped_personal_save_path(
    seed: int | None = None,
    *,
    timestamp: str | None = None,
    n_qubits: int | None = None,
    multi_qubit: bool = False,
    backend_folder: str = "backend_unknown",
) -> Path:
    timestamp = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    seed_folder = "seed_unseeded" if seed is None else f"seed_{seed}"
    qubit_folder = "q_unknown" if n_qubits is None else f"q{int(n_qubits)}"
    suffix = f"_q{int(n_qubits)}" if multi_qubit and n_qubits is not None else ""
    run_folder = (
        Path("Personal")
        / "FLE"
        / "Fancy_emulator"
        / "V_warp"
        / backend_folder
        / qubit_folder
        / seed_folder
        / f"FLE_V_warp_{timestamp}"
    )
    return run_folder / f"FLE_{timestamp}{suffix}.json"


def main() -> None:
    h_wrapper_run.FantasySettings = settings_module.FantasySettings
    h_wrapper_run.control_panel_settings_kwargs = (
        settings_module.control_panel_settings_kwargs
    )
    h_wrapper_run.RNG_SEEDS = settings_module.RNG_SEEDS
    h_wrapper_run.PLOT = settings_module.PLOT
    h_wrapper_run.timestamped_personal_save_path = timestamped_personal_save_path
    h_wrapper_run.main()


if __name__ == "__main__":
    main()
