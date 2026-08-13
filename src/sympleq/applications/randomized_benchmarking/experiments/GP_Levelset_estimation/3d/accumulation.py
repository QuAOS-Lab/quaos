from __future__ import annotations

import copy
import json
import re
from datetime import datetime
from pathlib import Path

COST_AWARE = True
H2_DEVICE = "H2-2" if COST_AWARE else "H2_2"
FLE_SEEDS = [2026, 2027, 2028, 2029]
COST_AWARE_SEEDS = {
    "H2-1": [2026, 20261, 20267, 20268],
    "H2-2": [20262, 20266, 20269, 20270],
}
SEEDS = COST_AWARE_SEEDS[H2_DEVICE.replace("_", "-")] if COST_AWARE else FLE_SEEDS

# If this list is non-empty, the script accumulates these JSONs directly
# instead of searching seed folders. Use this for accumulation-of-accumulations.
EXPLICIT_SOURCE_JSONS: list[Path] = [
    Path(
        r"Personal\FLE\H2_2\qband_4\accumulation"
        r"\accumulated_final_H2_2_qband_4_20260804_090823_actual_gates.json"
    ),
    Path(
        r"Personal\Cost_aware\H2-2\accumulation"
        r"\accumulated_final_H2-2_cost_aware_20260811_160022.json"
    ),
]
EXPLICIT_OUTPUT_FOLDER: Path | None = Path(r"Personal\Cost_aware\H2-2\accumulation")
EXPLICIT_OUTPUT_TAG = "H2-2_FLE_cost_aware_combined_accumulations"


# 1. Folder creation.
if COST_AWARE:
    source_root = Path("Personal") / "Cost_aware" / H2_DEVICE.replace("_", "-")
else:
    source_root = Path("Personal") / "FLE" / H2_DEVICE.replace("-", "_") / "qband_4"
if EXPLICIT_SOURCE_JSONS:
    output_folder = EXPLICIT_OUTPUT_FOLDER or (
        Path("Personal") / ("Cost_aware" if COST_AWARE else "FLE") / "accumulation"
    )
else:
    output_folder = source_root / "accumulation"
output_folder.mkdir(parents=True, exist_ok=True)


# 2. Accumulation of final JSONs.
def measurement_step(path: Path) -> int:
    match = re.match(r"measurement_(\d+)_", path.name)
    return int(match.group(1)) if match else -1


final_paths = []
if EXPLICIT_SOURCE_JSONS:
    for path in EXPLICIT_SOURCE_JSONS:
        path = Path(path)
        if not path.exists():
            print(f"[skip] explicit source missing: {path}")
            continue
        final_paths.append(path)
        print(f"[source] {path}")
else:
    for seed in SEEDS:
        if COST_AWARE:
            paths = [source_root / f"seed_{seed}" / "CostAware_restartable.json"]
            paths = [path for path in paths if path.exists()]
        else:
            paths = list((source_root / f"seed_{seed}").rglob("measurement_*.json"))
        if not paths:
            print(f"[skip] seed_{seed}: no source JSONs")
            continue
        final_path = paths[0] if COST_AWARE else max(paths, key=measurement_step)
        final_paths.append(final_path)
        print(f"[final] seed_{seed}: {final_path}")

if not final_paths:
    raise RuntimeError(f"No final JSONs found under {source_root}")

payload = copy.deepcopy(json.loads(final_paths[0].read_text(encoding="utf-8")))
payload["data"] = []
payload["experiment"] = {
    "source": "accumulated_explicit_jsons"
    if EXPLICIT_SOURCE_JSONS
    else (
        "accumulated_cost_aware_final_jsons"
        if COST_AWARE
        else "accumulated_final_jsons"
    ),
    "cost_aware": COST_AWARE,
    "h2_device": H2_DEVICE,
    "source_jsons": [str(path) for path in final_paths],
    "created_at": datetime.now().isoformat(timespec="seconds"),
}
if EXPLICIT_SOURCE_JSONS:
    payload["experiment"]["explicit_source_jsons"] = True
    payload["experiment"]["output_tag"] = EXPLICIT_OUTPUT_TAG
else:
    payload["experiment"]["seeds"] = SEEDS
if not COST_AWARE:
    payload["experiment"]["qband"] = 4

for path in final_paths:
    payload["data"].extend(json.loads(path.read_text(encoding="utf-8")).get("data", []))

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
if EXPLICIT_SOURCE_JSONS:
    output_path = output_folder / f"accumulated_{EXPLICIT_OUTPUT_TAG}_{timestamp}.json"
elif COST_AWARE:
    output_path = output_folder / f"accumulated_final_{H2_DEVICE}_cost_aware_{timestamp}.json"
else:
    output_path = output_folder / f"accumulated_final_{H2_DEVICE}_qband_4_{timestamp}.json"
output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(f"[saved] {output_path}")
print(f"[data] total records={len(payload['data'])}")
