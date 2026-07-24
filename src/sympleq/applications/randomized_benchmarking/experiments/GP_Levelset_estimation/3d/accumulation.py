from __future__ import annotations

import copy
import json
import re
from datetime import datetime
from pathlib import Path

H2_DEVICE = "H2_2"
SEEDS = [2026, 2027, 2028,
         2029,
         ]


# 1. Folder creation.
source_root = Path("Personal") / "FLE" / H2_DEVICE / "qband_4"
output_folder = source_root / "accumulation"
output_folder.mkdir(parents=True, exist_ok=True)


# 2. Accumulation of final JSONs.
def measurement_step(path: Path) -> int:
    match = re.match(r"measurement_(\d+)_", path.name)
    return int(match.group(1)) if match else -1


final_paths = []
for seed in SEEDS:
    paths = list((source_root / f"seed_{seed}").rglob("measurement_*.json"))
    if not paths:
        print(f"[skip] seed_{seed}: no measurement JSONs")
        continue
    final_path = max(paths, key=measurement_step)
    final_paths.append(final_path)
    print(f"[final] seed_{seed}: {final_path}")

if not final_paths:
    raise RuntimeError(f"No final JSONs found under {source_root}")

payload = copy.deepcopy(json.loads(final_paths[0].read_text(encoding="utf-8")))
payload["data"] = []
payload["experiment"] = {
    "source": "accumulated_final_jsons",
    "h2_device": H2_DEVICE,
    "qband": 4,
    "seeds": SEEDS,
    "source_jsons": [str(path) for path in final_paths],
    "created_at": datetime.now().isoformat(timespec="seconds"),
}

for path in final_paths:
    payload["data"].extend(json.loads(path.read_text(encoding="utf-8")).get("data", []))

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_path = output_folder / f"accumulated_final_{H2_DEVICE}_qband_4_{timestamp}.json"
output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
print(f"[saved] {output_path}")
print(f"[data] total records={len(payload['data'])}")
