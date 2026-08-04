from __future__ import annotations

import copy
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")

from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.plot_score_fle import (  # noqa: E402
    plot_fle_grid,
)


# Handles.
ACCUMULATED_JSON = Path(
    r"Personal\FLE\H2_2\qband_4\accumulation\accumulated_final_H2_2_qband_4_20260804_090823_actual_gates.json"
)
GRID_3D = Path(
    r"Personal\FLE\H2_2\uniform_10_5000_grids\qband_4\accumulation\accumulated_final_H2_2_qband_4_20260804_090823_actual_gates_gp_grid_3d.npz"
)
QUBIT_GAP = 3
QUBITS: list[int] | None = None


def spaced_qubits(records: list[dict]) -> list[int]:
    if QUBITS is not None:
        return QUBITS
    available = sorted({int(record["n_qubits"]) for record in records})
    picked: list[int] = []
    for q in available:
        if not picked or q - picked[-1] >= QUBIT_GAP:
            picked.append(q)
    if available and picked[-1] != available[-1]:
        picked.append(available[-1])
    return picked


def grid_slice(values: np.ndarray, qubits_axis: np.ndarray, q: int) -> np.ndarray:
    hi = int(np.searchsorted(qubits_axis, float(q), side="left"))
    if hi <= 0:
        return values[0]
    if hi >= len(qubits_axis):
        return values[-1]
    lo = hi - 1
    weight = (float(q) - qubits_axis[lo]) / (qubits_axis[hi] - qubits_axis[lo])
    return (1.0 - weight) * values[lo] + weight * values[hi]


payload = json.loads(ACCUMULATED_JSON.read_text(encoding="utf-8"))
records = payload.get("data", [])
grid = np.load(GRID_3D)
qubits_axis = np.asarray(grid["qubits_axis"], dtype=float)
output_folder = ACCUMULATED_JSON.parent / f"{ACCUMULATED_JSON.stem}_q_slices_gap_{QUBIT_GAP}"
output_folder.mkdir(parents=True, exist_ok=True)

for q in spaced_qubits(records):
    q_records = [record for record in records if int(record["n_qubits"]) == q]
    if not q_records:
        print(f"[skip] q={q}: no records in {ACCUMULATED_JSON}")
        continue

    slice_payload = copy.deepcopy(payload)
    slice_payload["data"] = q_records
    slice_payload["experiment"] = {
        **payload.get("experiment", {}),
        "source": "accumulated_qubit_slice",
        "source_json": str(ACCUMULATED_JSON),
        "source_grid_3d": str(GRID_3D),
        "n_qubits": q,
    }

    slice_json = output_folder / f"{ACCUMULATED_JSON.stem}_q{q:02d}.json"
    slice_grid = output_folder / f"{slice_json.stem}_gp_grid.npz"
    slice_json.write_text(json.dumps(slice_payload, indent=2), encoding="utf-8")
    np.savez_compressed(
        slice_grid,
        gates_grid=grid_slice(np.asarray(grid["gates_grid"], dtype=float), qubits_axis, q),
        ratio_grid=grid_slice(np.asarray(grid["ratio_grid"], dtype=float), qubits_axis, q),
        probabilities=grid_slice(np.asarray(grid["probabilities"], dtype=float), qubits_axis, q),
        latent_mean=grid_slice(np.asarray(grid["latent_mean"], dtype=float), qubits_axis, q),
        latent_variance=grid_slice(np.asarray(grid["latent_variance"], dtype=float), qubits_axis, q),
        target=np.asarray(grid["target"], dtype=float),
    )
    print(f"[slice] q={q} records={len(slice_payload['data'])} json={slice_json}")
    plot_fle_grid(slice_json)
