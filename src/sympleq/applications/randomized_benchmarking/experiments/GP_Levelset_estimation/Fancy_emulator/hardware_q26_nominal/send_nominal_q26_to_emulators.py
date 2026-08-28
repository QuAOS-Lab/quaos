from __future__ import annotations

import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import qnexus as qnx
from numpy.random import default_rng
from qnexus.models.references import IncompleteJobItemRef

from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.Fancy_emulator.unstitched_quantinuum import (
    fancy_noisy_emulator_config,
)
from sympleq.core.circuits.gates import GATES
from sympleq.integrations.quantinuum.utils import (
    BASE_SIMULATION_COST,
    pytket_bare_simulation_cost,
    to_pytket_circuit,
)
from sympleq.integrations.quantinuum.workflow import (
    build_and_compile_circuits,
    setup,
)


QUBITS = 26
N_SHOTS = 1
OUTPUT_ROOT = Path(r"Personal\FLE\Fancy_emulator\hardware_q26_nominal")
REUSE_CIRCUIT_LIST_JSON: Path | None = Path(
    r"Personal\FLE\Fancy_emulator\hardware_q26_nominal"
    r"\nominal_q26_20260827_134207"
    r"\nominal_q26_circuit_list.json"
)
TARGET_DEVICES_TO_SEND = ("H2-2E",)
WAIT_FOR_RESULTS = False

SOURCE_FILES = [
    (
        "H2-1",
        2026,
        Path(
            r"Personal\Data\FLE\H2_1\qband_4\seed_2026"
            r"\FLE_20260708_175833"
            r"\measurement_018_globalsur_20260710_200745_087360.json"
        ),
    ),
    (
        "H2-1",
        2027,
        Path(
            r"Personal\Data\FLE\H2_1\qband_4\seed_2027"
            r"\FLE_20260710_173814"
            r"\measurement_018_globalsur_20260712_035438_338739.json"
        ),
    ),
    (
        "H2-1",
        2028,
        Path(
            r"Personal\Data\FLE\H2_1\qband_4\seed_2028"
            r"\FLE_20260712_102620"
            r"\measurement_018_globalsur_20260716_043515_296732.json"
        ),
    ),
    (
        "H2-1",
        2029,
        Path(
            r"Personal\Data\FLE\H2_1\qband_4\seed_2029"
            r"\FLE_20260717_083507"
            r"\measurement_018_globalsur_20260719_185930_205592.json"
        ),
    ),
    (
        "H2-2",
        2026,
        Path(
            r"Personal\Data\FLE\H2_2\qband_4\seed_2026"
            r"\FLE_20260708_113322"
            r"\measurement_019_globalsur_20260719_144303_503855.json"
        ),
    ),
    (
        "H2-2",
        2027,
        Path(
            r"Personal\Data\FLE\H2_2\qband_4\seed_2027"
            r"\FLE_20260719_172406"
            r"\measurement_018_globalsur_20260721_093645_150943.json"
        ),
    ),
    (
        "H2-2",
        2028,
        Path(
            r"Personal\Data\FLE\H2_2\qband_4\seed_2028"
            r"\FLE_20260720_145715"
            r"\measurement_018_globalsur_20260724_055405_380918.json"
        ),
    ),
    (
        "H2-2",
        2029,
        Path(
            r"Personal\Data\FLE\H2_2\qband_4\seed_2029"
            r"\FLE_20260724_102820"
            r"\measurement_018_globalsur_20260730_161735_139460.json"
        ),
    ),
]

TARGET_DEVICE = {
    "H2-1": "H2-1E",
    "H2-2": "H2-2E",
}

PROJECT_NAME = {
    "H2-1E": "Fancy_Emulator_Nominal_Q26_H21E",
    "H2-2E": "Fancy_Emulator_Nominal_Q26_H22E",
}


def measurement_rng(seed: int, config: RMBConfig, shot_index: int):
    return default_rng(
        [
            int(seed),
            int(config.n_qubits),
            int(config.n_1qb_gates),
            int(config.n_2qb_gates),
            int(shot_index),
        ]
    )


def gate_set_from_names(names: list[str]):
    aliases = {
        "ZZP": "ZZMax",
        "ZZP_inv": "ZZMax_inv",
    }
    return tuple(getattr(GATES, aliases.get(str(name), str(name))) for name in names)


def nominal_config(record: dict) -> RMBConfig:
    return (
        RMBConfig.default()
        .with_n_qubits(int(record["n_qubits"]))
        .with_n_1qb_gates(int(record["n_1qb_gates"]))
        .with_n_2qb_gates(int(record["n_2qb_gates"]))
        .with_random_elimination(float(record.get("random_elimination", 0.1)))
        .with_use_scrambler(bool(record.get("use_scrambler", True)))
        .with_gates_set(gate_set_from_names(record["gates_set"]))
    )


def nominal_circuit(config: RMBConfig, rng):
    circuit = to_pytket_circuit(config.random_circuit(rng=rng))
    if circuit.n_gates == 0:
        raise RuntimeError(f"Generated empty circuit for config {config!r}.")
    return circuit


def result_counts(record: dict) -> list[tuple[bool, int]]:
    return [(bool(outcome), int(count)) for outcome, count in record["results"]]


def config_key(config: RMBConfig) -> tuple:
    return (
        config.n_qubits,
        config.n_1qb_gates,
        config.n_2qb_gates,
        config.random_elimination,
        config.use_scrambler,
        tuple(gate.name for gate in config.gates_set),
    )


def outcome_counts(result) -> dict[str, int]:
    counts = result.get_empirical_distribution().as_counter()
    return {
        "".join(str(int(bit)) for bit in outcome): int(count)
        for outcome, count in counts.items()
    }


def success_from_result(result) -> bool:
    counts = result.get_empirical_distribution().as_counter()
    if not counts:
        return False
    top_outcome, _ = counts.most_common()[0]
    return all(bit == 0 for bit in top_outcome)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[write] {path}", flush=True)


def identifier_strings(obj) -> list[str]:
    values = []
    for attr in ("id", "uuid", "item_id", "job_id"):
        value = getattr(obj, attr, None)
        if value is not None:
            values.append(str(value))
    for attr in ("ref", "annotations"):
        nested = getattr(obj, attr, None)
        if nested is None:
            continue
        for nested_attr in ("id", "uuid", "item_id", "job_id"):
            value = getattr(nested, nested_attr, None)
            if value is not None:
                values.append(str(value))
    values.append(str(obj))
    return list(dict.fromkeys(values))


def best_identifier(obj) -> str:
    return identifier_strings(obj)[0]


def load_nominal_circuit_list() -> list[dict]:
    items = []
    shot_offsets = defaultdict(int)

    for source_device, seed, path in SOURCE_FILES:
        target_device = TARGET_DEVICE[source_device]
        if target_device not in TARGET_DEVICES_TO_SEND:
            print(
                "[source-skip] "
                f"device={source_device} target={target_device} "
                "not in TARGET_DEVICES_TO_SEND",
                flush=True,
            )
            continue

        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)

        q26_records = 0
        q26_circuits = 0
        for record_index, record in enumerate(payload["data"]):
            if int(record["n_qubits"]) != QUBITS:
                continue

            q26_records += 1
            config = nominal_config(record)
            key = (source_device, seed, config_key(config))

            for hardware_success, count in result_counts(record):
                for _ in range(count):
                    shot_index = shot_offsets[key]
                    shot_offsets[key] += 1

                    rng = measurement_rng(seed, config, shot_index)
                    circuit = nominal_circuit(config, rng)
                    item_index = len(items)
                    item = {
                        "index": item_index,
                        "source_device": source_device,
                        "target_device": target_device,
                        "seed": seed,
                        "source_json": str(path),
                        "record_index": record_index,
                        "shot_index": shot_index,
                        "hardware_success": hardware_success,
                        "nominal_n_1qb_gates": config.n_1qb_gates,
                        "nominal_n_2qb_gates": config.n_2qb_gates,
                        "nominal_n_gates": config.n_gates,
                        "nominal_ratio_2qb_gates": config.ratio_2_qb_gates,
                        "pytket_n_qubits": circuit.n_qubits,
                        "pytket_n_bits": circuit.n_bits,
                        "pytket_n_gates": circuit.n_gates,
                        "pytket_n_1qb_gates": circuit.n_1qb_gates(),
                        "pytket_n_2qb_gates": circuit.n_2qb_gates(),
                    }
                    items.append({"metadata": item, "circuit": circuit})
                    q26_circuits += 1
                    print(
                        "[circuit-list] "
                        f"{item_index:04d} "
                        f"source={source_device} target={target_device} "
                        f"seed={seed} record={record_index} shot={shot_index} "
                        f"hardware_success={hardware_success} "
                        f"q={config.n_qubits} "
                        f"nominal_n1={config.n_1qb_gates} "
                        f"nominal_n2={config.n_2qb_gates} "
                        f"nominal_total={config.n_gates} "
                        f"nominal_ratio={config.ratio_2_qb_gates:.6f} "
                        f"pytket_n1={circuit.n_1qb_gates()} "
                        f"pytket_n2={circuit.n_2qb_gates()} "
                        f"pytket_gates={circuit.n_gates}",
                        flush=True,
                    )

        print(
            "[source] "
            f"device={source_device} seed={seed} "
            f"q26_records={q26_records} q26_circuits={q26_circuits} "
            f"path={path}",
            flush=True,
        )

    return items


def config_from_manifest(record: dict) -> RMBConfig:
    return (
        RMBConfig.default()
        .with_n_qubits(int(record["pytket_n_qubits"]))
        .with_n_1qb_gates(int(record["nominal_n_1qb_gates"]))
        .with_n_2qb_gates(int(record["nominal_n_2qb_gates"]))
        .with_random_elimination(0.1)
        .with_use_scrambler(True)
        .with_gates_set(gate_set_from_names([
            "Id",
            "H",
            "H_inv",
            "S",
            "S_inv",
            "V",
            "V_inv",
            "X",
            "X_inv",
            "Y",
            "Y_inv",
            "Z",
            "Z_inv",
            "CX",
            "CX_inv",
            "SWAP",
            "CZ",
            "ZZP",
            "ZZP_inv",
        ]))
    )


def load_nominal_circuit_list_from_manifest(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    items = []
    for record in manifest["circuits"]:
        if record["target_device"] not in TARGET_DEVICES_TO_SEND:
            continue

        config = config_from_manifest(record)
        rng = measurement_rng(int(record["seed"]), config, int(record["shot_index"]))
        circuit = nominal_circuit(config, rng)
        metadata = dict(record)
        items.append({"metadata": metadata, "circuit": circuit})
        print(
            "[reuse-circuit-list] "
            f"{metadata['index']:04d} "
            f"source={metadata['source_device']} "
            f"target={metadata['target_device']} "
            f"seed={metadata['seed']} "
            f"record={metadata['record_index']} "
            f"shot={metadata['shot_index']} "
            f"pytket_gates={circuit.n_gates}",
            flush=True,
        )

    print(
        "[reuse-circuit-list] "
        f"loaded {len(items)} circuits from {path} "
        f"for targets={TARGET_DEVICES_TO_SEND}",
        flush=True,
    )
    return items


def submission_cost(items: list[dict]) -> float:
    return float(BASE_SIMULATION_COST) + sum(
        pytket_bare_simulation_cost(item["circuit"])
        for item in items
    )

def run_device(
    items: list[dict],
    target_device: str,
    *,
    manifest: dict,
    manifest_path: Path,
) -> list[dict]:
    device_items = [
        item for item in items
        if item["metadata"]["target_device"] == target_device
    ]
    if not device_items:
        return []

    backend_config = fancy_noisy_emulator_config(target_device)
    project_name = PROJECT_NAME[target_device]
    print(
        "[backend-config] "
        f"target={target_device} project={project_name} simulator=state-vector "
        "noisy_simulation=True stitching=unstitched",
        flush=True,
    )
    setup(project_name)

    completed = []
    for batch_index, batch in enumerate([device_items]):
        cost = submission_cost(batch)
        circuits = [item["circuit"] for item in batch]
        print(
            "[submit] "
            f"target={target_device} batch={batch_index:03d} "
            f"programs={len(circuits)} estimated_cost={cost:.6f}",
            flush=True,
        )
        refs = build_and_compile_circuits(
            circuits,
            backend_config=backend_config,
            name=f"nominal-q26-{target_device}-{batch_index:03d}",
        )
        execute_job = qnx.start_execute_job(
            programs=refs,
            n_shots=N_SHOTS,
            backend_config=backend_config,
            name=f"nominal-q26-{target_device}-{batch_index:03d}",
        )
        job_id = best_identifier(execute_job)
        for item in batch:
            item["metadata"]["execute_job_id"] = job_id
            item["metadata"]["execute_job_identifiers"] = identifier_strings(execute_job)
            item["metadata"]["execute_job_batch_index"] = batch_index
            item["metadata"]["execute_job_project_name"] = project_name

        batch_manifest = {
            "target_device": target_device,
            "project_name": project_name,
            "batch_index": batch_index,
            "execute_job_id": job_id,
            "execute_job_identifiers": identifier_strings(execute_job),
            "circuit_indices": [
                int(item["metadata"]["index"])
                for item in batch
            ],
            "n_programs": len(batch),
            "estimated_cost": cost,
            "status": "submitted",
        }
        manifest["batches"].append(batch_manifest)
        manifest["circuits"] = [item["metadata"] for item in items]
        write_json(manifest_path, manifest)
        print(
            "[job] "
            f"target={target_device} batch={batch_index:03d} "
            f"job_id={job_id}",
            flush=True,
        )

        if not WAIT_FOR_RESULTS:
            batch_manifest["status"] = "submitted_no_wait"
            manifest["circuits"] = [item["metadata"] for item in items]
            write_json(manifest_path, manifest)
            print(
                "[submit-only] "
                f"target={target_device} batch={batch_index:03d} "
                "execute job id saved; not waiting for emulator results",
                flush=True,
            )
            continue

        qnx.jobs.wait_for(execute_job)
        ref_results = qnx.jobs.results(execute_job)
        if isinstance(ref_results, IncompleteJobItemRef):
            raise RuntimeError(f"Qnexus execution failed: {ref_results}")
        results = [
            ref_result.download_result()
            for ref_result in ref_results
            if not isinstance(ref_result, IncompleteJobItemRef)
        ]
        batch_manifest["status"] = "completed"
        manifest["circuits"] = [item["metadata"] for item in items]
        write_json(manifest_path, manifest)
        if len(results) != len(batch):
            raise RuntimeError(
                f"{target_device}: got {len(results)} results for "
                f"{len(batch)} circuits."
            )

        for item, result in zip(batch, results):
            metadata = dict(item["metadata"])
            metadata["emulator_success"] = success_from_result(result)
            metadata["emulator_counts"] = outcome_counts(result)
            completed.append(metadata)
            print(
                "[result] "
                f"{metadata['index']:04d} target={target_device} "
                f"hardware_success={metadata['hardware_success']} "
                f"emulator_success={metadata['emulator_success']} "
                f"counts={metadata['emulator_counts']}",
                flush=True,
            )

    return completed


def summary(records: list[dict]) -> dict:
    by_target = defaultdict(int)
    by_source = defaultdict(int)
    for record in records:
        by_target[record["target_device"]] += 1
        by_source[record["source_device"]] += 1
    return {
        "total_circuits": len(records),
        "by_source_device": dict(by_source),
        "by_target_device": dict(by_target),
    }


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = OUTPUT_ROOT / f"nominal_q26_{timestamp}"

    if REUSE_CIRCUIT_LIST_JSON is None:
        items = load_nominal_circuit_list()
    else:
        items = load_nominal_circuit_list_from_manifest(REUSE_CIRCUIT_LIST_JSON)
    metadata = [item["metadata"] for item in items]
    manifest = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "description": (
            "Nominal q26 circuits from H2-1/H2-2 hardware FLE final JSONs, "
            "sent unstitched to noisy Quantinuum emulators."
        ),
        "q": QUBITS,
        "n_shots_per_circuit": N_SHOTS,
        "project_name_by_target_device": PROJECT_NAME,
        "simulator": "state-vector",
        "noisy_simulation": True,
        "stitching": "unstitched",
        "target_devices_to_send": list(TARGET_DEVICES_TO_SEND),
        "wait_for_results": WAIT_FOR_RESULTS,
        "source_files": [
            {
                "source_device": source_device,
                "seed": seed,
                "path": str(path),
                "target_device": TARGET_DEVICE[source_device],
            }
            for source_device, seed, path in SOURCE_FILES
        ],
        "summary": summary(metadata),
        "batches": [],
        "circuits": metadata,
    }
    circuit_list_path = output_dir / "nominal_q26_circuit_list.json"
    job_manifest_path = output_dir / "nominal_q26_job_manifest.json"
    write_json(circuit_list_path, manifest)
    write_json(job_manifest_path, manifest)

    results = []
    for target_device in TARGET_DEVICES_TO_SEND:
        results.extend(
            run_device(
                items,
                target_device,
                manifest=manifest,
                manifest_path=job_manifest_path,
            )
        )

    if WAIT_FOR_RESULTS:
        output = dict(manifest)
        output["completed_at"] = datetime.now().isoformat(timespec="seconds")
        output["summary"] = summary(results)
        output["circuits"] = sorted(results, key=lambda record: record["index"])
        write_json(output_dir / "nominal_q26_emulator_results.json", output)
    else:
        manifest["submitted_at"] = datetime.now().isoformat(timespec="seconds")
        manifest["summary"] = summary([item["metadata"] for item in items])
        write_json(job_manifest_path, manifest)
        print(
            "[done] submit-only mode; use report_nominal_q26_from_jobs.py "
            "with nominal_q26_job_manifest.json to fetch results later",
            flush=True,
        )


if __name__ == "__main__":
    main()
