from __future__ import annotations

import copy
import json
import re
from pathlib import Path

import numpy as np

from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.circuits.gates import DEFAULT_GATES_SET
from sympleq.integrations.quantinuum.utils import NATIVE_GATES_SET


# Handles.
INPUT_JSON = Path(
    r"Personal\FLE\H2_2\qband_4\accumulation\accumulated_final_H2_2_qband_4_20260804_090823.json"
)
OUTPUT_JSON = INPUT_JSON.with_name(f"{INPUT_JSON.stem}_actual_gates.json")


GATE_BY_NAME = {
    gate.name: gate
    for gate in [*DEFAULT_GATES_SET, *NATIVE_GATES_SET]
}


def seed_from_path(path: Path) -> int:
    match = re.search(r"seed_(\d+)", str(path))
    if match is None:
        raise ValueError(f"Could not infer seed from source path: {path}")
    return int(match.group(1))


def config_from_record(record: dict) -> RMBConfig:
    gates_set = tuple(
        GATE_BY_NAME[name]
        for name in record.get("gates_set", [])
    )
    if not gates_set:
        gates_set = tuple(NATIVE_GATES_SET)
    return RMBConfig(
        n_1qb_gates=int(record["n_1qb_gates"]),
        n_2qb_gates=int(record["n_2qb_gates"]),
        n_qubits=int(record["n_qubits"]),
        gates_set=gates_set,
        random_elimination=float(record.get("random_elimination", 0.0)),
        use_scrambler=bool(record.get("use_scrambler", True)),
    )


def n_results(record: dict) -> int:
    return sum(int(count) for _, count in record.get("results", []))


def shot_rng(seed: int, config: RMBConfig, shot_index: int) -> np.random.Generator:
    return np.random.default_rng(
        [
            seed,
            config.n_qubits,
            config.n_1qb_gates,
            config.n_2qb_gates,
            shot_index,
        ]
    )


def non_identity_counts(config: RMBConfig, seed: int, shot_index: int) -> tuple[int, int]:
    circuit = config.random_circuit(rng=shot_rng(seed, config, shot_index))
    n_1q = sum(
        1
        for gate in circuit.gates
        if gate.n_qudits == 1 and gate.name != "Id"
    )
    n_2q = sum(
        1
        for gate in circuit.gates
        if gate.n_qudits == 2 and gate.name != "Id"
    )
    return n_1q, n_2q


def same_config(left: dict, right: dict) -> bool:
    keys = ["n_1qb_gates", "n_2qb_gates", "n_qubits", "random_elimination", "use_scrambler"]
    return all(left.get(key) == right.get(key) for key in keys)


def source_records(payload: dict) -> list[tuple[dict, int]]:
    source_jsons = payload.get("experiment", {}).get("source_jsons", [])
    if not source_jsons:
        raise RuntimeError(
            "Expected experiment.source_jsons in the accumulated JSON, "
            "so each record can be mapped back to its seed."
        )

    records_with_seed: list[tuple[dict, int]] = []
    for source_json in source_jsons:
        path = Path(source_json)
        seed = seed_from_path(path)
        source_payload = json.loads(path.read_text(encoding="utf-8"))
        source_data = source_payload.get("data", [])
        print(f"[source] seed={seed} records={len(source_data)} path={path}")
        records_with_seed.extend((record, seed) for record in source_data)
    return records_with_seed


payload = json.loads(INPUT_JSON.read_text(encoding="utf-8"))
rewritten = copy.deepcopy(payload)
records = rewritten.get("data", [])
source_data = source_records(payload)

if len(source_data) != len(records):
    raise RuntimeError(
        f"Accumulated/source record count mismatch: "
        f"{len(records)} accumulated vs {len(source_data)} from source_jsons."
    )

total_nominal_1q = 0
total_actual_1q = 0
total_nominal_2q = 0
total_actual_2q = 0
repeated_records = 0

for index, (record, (source_record, seed)) in enumerate(zip(records, source_data), start=1):
    if not same_config(record, source_record):
        raise RuntimeError(f"Record {index} does not match its source JSON record.")

    config = config_from_record(record)
    shots = n_results(record)
    if shots <= 0:
        print(f"[skip] record={index}: no recorded shots")
        continue
    if shots > 1:
        repeated_records += 1

    counts = [
        non_identity_counts(config, seed, shot_index)
        for shot_index in range(shots)
    ]
    actual_1q = int(round(sum(n_1q for n_1q, _ in counts) / shots))
    actual_2q = int(round(sum(n_2q for _, n_2q in counts) / shots))

    total_nominal_1q += int(record["n_1qb_gates"]) * shots
    total_actual_1q += actual_1q * shots
    total_nominal_2q += int(record["n_2qb_gates"]) * shots
    total_actual_2q += actual_2q * shots

    record["n_1qb_gates"] = actual_1q
    record["n_2qb_gates"] = actual_2q
    if "n_gates" in record:
        record["n_gates"] = actual_1q + actual_2q
    if "ratio_2_qb_gates" in record:
        total = actual_1q + actual_2q
        record["ratio_2_qb_gates"] = actual_2q / total if total else 0.0

    if index % 100 == 0 or shots > 1:
        spread_1q = sorted({n_1q for n_1q, _ in counts})
        print(
            f"[record] {index}/{len(records)} seed={seed} shots={shots} "
            f"1q {config.n_1qb_gates}->{actual_1q} spread={spread_1q} "
            f"2q {config.n_2qb_gates}->{actual_2q}"
        )

if OUTPUT_JSON.resolve() == INPUT_JSON.resolve():
    raise RuntimeError("OUTPUT_JSON must be different from INPUT_JSON.")

OUTPUT_JSON.write_text(json.dumps(rewritten, indent=2), encoding="utf-8")

print(f"[saved] {OUTPUT_JSON}")
print(f"[records] {len(records)}")
print(f"[repeated-records] {repeated_records}")
print(f"[1q] nominal={total_nominal_1q} actual={total_actual_1q}")
print(f"[2q] nominal={total_nominal_2q} actual={total_actual_2q}")
