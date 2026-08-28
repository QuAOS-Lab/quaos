from __future__ import annotations

from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.Fancy_emulator.unstitched_quantinuum import (
    FancyUnstitchedQuantinuumBackend,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.h_wrapper_config import (
    H_WRAPPER_1Q_GATES_PER_QUBIT,
    H_WRAPPER_VARIANT,
)


class VwarpFancyUnstitchedQuantinuumBackend(FancyUnstitchedQuantinuumBackend):
    """Noisy Quantinuum emulator backend for native-V wrapper experiments."""

    log_prefix = "[vwarp-fancy]"
    program_name_prefix = "vwarp-fancy-unstitched"

    def to_dict(self) -> dict:
        payload = super().to_dict()
        payload.update(
            {
                "backend_family": "vwarp_fancy_emulator",
                "circuit_variant": H_WRAPPER_VARIANT,
                "protected_v_wrapper": True,
                "v_wrapper_1q_gates_per_qubit": H_WRAPPER_1Q_GATES_PER_QUBIT,
            }
        )
        return payload


def _vwarp_fancy_backend(settings, device_name: str, project_name: str):
    return VwarpFancyUnstitchedQuantinuumBackend(
        device_name=device_name,
        project_name=project_name,
        batch_size=1,
        max_cost_per_run=settings.max_cost_per_run,
    )


def vwarp_fancy_h21e(settings, rng):
    return _vwarp_fancy_backend(
        settings,
        device_name="H2-1E",
        project_name="Vwarp_Fancy_H21E",
    )


def vwarp_fancy_h22e(settings, rng):
    return _vwarp_fancy_backend(
        settings,
        device_name="H2-2E",
        project_name="Vwarp_Fancy_H22E",
    )
