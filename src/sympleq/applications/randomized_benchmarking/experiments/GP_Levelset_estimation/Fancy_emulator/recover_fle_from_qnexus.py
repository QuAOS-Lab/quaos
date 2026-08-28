"""Recover a pending Fancy-emulator FLE batch from completed QNexus jobs.

This is intentionally standalone: it does not resume the FLE runner, rebuild
the GP, or submit anything. It reads ``pending_backend_batch.json`` from a run
folder, rebuilds the circuits in the same order used by the Fancy unstitched
backend, fetches one or more completed QNexus execute jobs, appends the
recovered success/failure outcomes to the latest measurement checkpoint, and
writes a new timestamped measurement JSON in the same run folder.

Supported pending batches:
    * ordinary Fancy emulator: ``fancy_unstitched_emulator_backend_factory``
    * V-wrapper Fancy emulator: ``vwarp_fancy_h21e`` / ``vwarp_fancy_h22e``
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import qnexus as qnx
from numpy.random import default_rng
from pytket.backends.backendresult import BackendResult
from pytket.circuit import Circuit as PytketCircuit
from qnexus.client.auth import is_logged_in
from qnexus.models.filters import SortFilterEnum
from qnexus.models.job_status import JobStatusEnum
from qnexus.models.references import (
    CircuitRef,
    ExecutionResultRef,
    IncompleteJobItemRef,
    JobType,
)

from sympleq.applications.randomized_benchmarking.backends.quantinuum import (
    QuantinuumBackend,
)
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.h_wrapper_config import (
    HWrappedRMBConfig,
)
from sympleq.core.circuits.gates import GATES
from sympleq.integrations.quantinuum.utils import (
    BASE_SIMULATION_COST,
    NATIVE_GATES_SET,
    pytket_bare_simulation_cost,
    to_pytket_circuit,
)


KNOWN_BACKENDS = {
    "fancy_unstitched_emulator_backend_factory": {
        "device_name": "H2-1E",
        "project_name": "Fancy_Emulator_Unstitched_H21E",
        "wrapper": False,
    },
    "vwarp_fancy_h21e": {
        "device_name": "H2-1E",
        "project_name": "Vwarp_Fancy_H21E",
        "wrapper": True,
    },
    "vwarp_fancy_h22e": {
        "device_name": "H2-2E",
        "project_name": "Vwarp_Fancy_H22E",
        "wrapper": True,
    },
}

GATE_ALIASES = {
    "ZZP": "ZZMax",
    "ZZP_inv": "ZZMax_inv",
}


@dataclass(frozen=True)
class PendingCircuit:
    request_index: int
    shot_index: int
    config: RMBConfig
    circuit: PytketCircuit
    request_payload: dict


def debug(message: str) -> None:
    print(f"[fancy-recovery] {message}", flush=True)


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


@lru_cache(maxsize=None)
def single_circuit_bare_hqc(config: RMBConfig) -> float:
    circuit = to_pytket_circuit(config.random_circuit(rng=default_rng(0)))
    return pytket_bare_simulation_cost(circuit) + int(config.n_qubits) / 5000


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    debug(f"wrote recovered measurement JSON: {path}")


def pending_path(run_folder: Path) -> Path:
    return run_folder / "pending_backend_batch.json"


def load_pending(run_folder: Path) -> dict:
    path = pending_path(run_folder)
    if not path.exists():
        raise FileNotFoundError(f"No pending backend batch found: {path}")
    payload = read_json(path)
    payload["_path"] = str(path)
    return payload


def pending_requests(payload: dict) -> list[dict]:
    requests = payload.get("pending_batch", {}).get("requests", [])
    if not requests:
        raise ValueError("pending_backend_batch.json has no pending_batch.requests.")
    for index, request in enumerate(requests):
        if request.get("config") is None:
            raise ValueError(f"pending request {index} has no config.")
    return requests


def measurement_step(path: Path) -> int:
    parts = path.name.split("_")
    if len(parts) < 2 or parts[0] != "measurement":
        return 0
    try:
        return int(parts[1])
    except ValueError:
        return 0


def latest_measurement_checkpoint(run_folder: Path) -> Path:
    checkpoints = sorted(
        path
        for path in run_folder.glob("measurement_*.json")
        if "_actual_gates" not in path.stem
        and "_gp_grid" not in path.stem
        and "_combined" not in path.stem
    )
    if not checkpoints:
        raise FileNotFoundError(f"No measurement_*.json checkpoints in {run_folder}")
    return checkpoints[-1]


def output_measurement_path(run_folder: Path, pending: dict, *, partial: bool) -> Path:
    step = int(pending.get("step", 0))
    phase = str(pending.get("phase", "globalsur"))
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    suffix = "_partial_recovery" if partial else ""
    return run_folder / f"measurement_{step:03d}_{phase}{suffix}_{timestamp}.json"


def parse_job_ids(values: list[str] | None) -> list[str]:
    if not values:
        return []
    job_ids: list[str] = []
    for value in values:
        job_ids.extend(part.strip() for part in value.split(",") if part.strip())
    return job_ids


def device_from_folder(run_folder: Path) -> str | None:
    text = str(run_folder).replace("/", "\\")
    if "H2_1E" in text or "H2-1E" in text:
        return "H2-1E"
    if "H2_2E" in text or "H2-2E" in text:
        return "H2-2E"
    return None


def backend_factory_name(pending: dict) -> str:
    return str(pending.get("settings", {}).get("backend_factory", ""))


def is_wrapper_pending(pending: dict, run_folder: Path) -> bool:
    settings = pending.get("settings", {})
    backend_name = backend_factory_name(pending)
    if backend_name in KNOWN_BACKENDS:
        return bool(KNOWN_BACKENDS[backend_name]["wrapper"])
    return (
        bool(settings.get("protected_h_wrapper"))
        or "wrapper" in str(settings.get("circuit_variant", "")).lower()
        or "V_warp" in str(run_folder)
    )


def infer_device_name(pending: dict, run_folder: Path) -> str:
    backend_name = backend_factory_name(pending)
    if backend_name in KNOWN_BACKENDS:
        return str(KNOWN_BACKENDS[backend_name]["device_name"])
    device_name = device_from_folder(run_folder)
    if device_name is not None:
        return device_name
    raise ValueError(
        "Could not infer device name. Pass --device-name explicitly."
    )


def infer_project_name(pending: dict, run_folder: Path, device_name: str) -> str:
    backend_name = backend_factory_name(pending)
    if backend_name in KNOWN_BACKENDS:
        return str(KNOWN_BACKENDS[backend_name]["project_name"])
    if "V_warp" in str(run_folder):
        suffix = "H21E" if device_name == "H2-1E" else "H22E"
        return f"Vwarp_Fancy_{suffix}"
    suffix = "H21E" if device_name == "H2-1E" else "H22E"
    return f"Fancy_Emulator_Unstitched_{suffix}"


def gate_set_from_names(names: list[str] | None):
    if not names:
        return tuple(NATIVE_GATES_SET)
    gates = []
    for name in names:
        attr = GATE_ALIASES.get(str(name), str(name))
        try:
            gates.append(getattr(GATES, attr))
        except AttributeError as exc:
            raise AttributeError(f"Unknown gate name in pending config: {name!r}") from exc
    return tuple(gates)


def config_from_record(record: dict, *, wrapper: bool) -> RMBConfig:
    cls = HWrappedRMBConfig if wrapper else RMBConfig
    return (
        cls.default()
        .with_n_qubits(int(record["n_qubits"]))
        .with_n_1qb_gates(int(record["n_1qb_gates"]))
        .with_n_2qb_gates(int(record["n_2qb_gates"]))
        .with_random_elimination(float(record.get("random_elimination", 0.1)))
        .with_use_scrambler(bool(record.get("use_scrambler", True)))
        .with_gates_set(gate_set_from_names(record.get("gates_set")))
    )


def config_key(config: RMBConfig) -> tuple[int, int, int, float, bool]:
    return (
        int(config.n_1qb_gates),
        int(config.n_2qb_gates),
        int(config.n_qubits),
        float(config.random_elimination),
        bool(config.use_scrambler),
    )


def record_key(record: dict) -> tuple[int, int, int, float, bool]:
    return (
        int(record["n_1qb_gates"]),
        int(record["n_2qb_gates"]),
        int(record["n_qubits"]),
        float(record.get("random_elimination", 0.1)),
        bool(record.get("use_scrambler", True)),
    )


def existing_shot_offsets(previous_payload: dict) -> dict[tuple, int]:
    offsets: dict[tuple, int] = defaultdict(int)
    for record in previous_payload.get("data", []):
        offsets[record_key(record)] += sum(
            int(count) for _, count in record.get("results", [])
        )
    return offsets


def sent_config_payload(config: RMBConfig) -> dict:
    payload = {
        "n_1qb_gates": int(config.n_1qb_gates),
        "n_2qb_gates": int(config.n_2qb_gates),
        "n_qubits": int(config.n_qubits),
        "n_gates": int(config.n_gates),
        "ratio_2_qb_gates": float(config.ratio_2_qb_gates),
        "random_elimination": float(config.random_elimination),
        "use_scrambler": bool(config.use_scrambler),
    }
    if hasattr(config, "circuit_metadata"):
        payload.update(config.circuit_metadata())
    return payload


def actual_after_elimination_payload(
    config: RMBConfig,
    *,
    rng_seed: int | None,
    shot_index: int | None,
) -> dict:
    if rng_seed is None or shot_index is None:
        payload = {
            "n_1qb_gates": None,
            "n_2qb_gates": None,
            "n_qubits": int(config.n_qubits),
            "n_gates": None,
            "ratio_2_qb_gates": None,
            "first_shot_index": shot_index,
            "available": False,
        }
        if hasattr(config, "circuit_metadata"):
            payload.update(config.circuit_metadata())
        return payload

    circuit = config.random_circuit(
        rng=measurement_rng(int(rng_seed), config, int(shot_index))
    )
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
    n_gates = n_1q + n_2q
    payload = {
        "n_1qb_gates": int(n_1q),
        "n_2qb_gates": int(n_2q),
        "n_qubits": int(config.n_qubits),
        "n_gates": int(n_gates),
        "ratio_2_qb_gates": float(n_2q / n_gates) if n_gates else 0.0,
        "first_shot_index": int(shot_index),
        "available": True,
    }
    if hasattr(config, "circuit_metadata"):
        payload.update(config.circuit_metadata())
    return payload


def build_pending_circuits(
    *,
    pending: dict,
    requests: list[dict],
    previous_payload: dict,
    run_folder: Path,
) -> list[PendingCircuit]:
    settings = pending.get("settings", {})
    rng_seed = settings.get("rng_seed")
    if rng_seed is None:
        raise ValueError("Cannot recover Fancy emulator batch without settings.rng_seed.")

    wrapper = is_wrapper_pending(pending, run_folder)
    fallback_offsets = existing_shot_offsets(previous_payload)
    local_offsets: dict[tuple, int] = defaultdict(int)

    debug(f"config family={'V/H-wrapper' if wrapper else 'ordinary RMB'}")
    circuits: list[PendingCircuit] = []
    for request_index, request in enumerate(requests):
        config = config_from_record(request["config"], wrapper=wrapper)
        shots = int(request.get("requested_shots", 1))
        first_shot_index = request.get("first_shot_index")
        key = config_key(config)
        debug(
            "pending request "
            f"{request_index:02d}: shots={shots} "
            f"first_shot_index={first_shot_index} "
            f"n_1q={config.n_1qb_gates} n_2q={config.n_2qb_gates} "
            f"q={config.n_qubits} ratio={config.ratio_2_qb_gates:.6f}"
        )

        for shot_offset in range(shots):
            if first_shot_index is None:
                shot_index = fallback_offsets[key] + local_offsets[key]
            else:
                shot_index = int(first_shot_index) + shot_offset
            local_offsets[key] += 1
            circuit_rng = measurement_rng(int(rng_seed), config, int(shot_index))
            circuit = QuantinuumBackend.compatible_circuit(config, circuit_rng)
            circuits.append(
                PendingCircuit(
                    request_index=request_index,
                    shot_index=int(shot_index),
                    config=config,
                    circuit=circuit,
                    request_payload=request,
                )
            )
            debug(
                "  built circuit "
                f"shot_index={shot_index} qubits={circuit.n_qubits} "
                f"bits={circuit.n_bits} gates={circuit.n_gates} "
                f"n_1q={circuit.n_1qb_gates()} n_2q={circuit.n_2qb_gates()}"
            )

    return circuits


def pack_unstitched_submissions(
    circuits: list[PendingCircuit],
    *,
    max_cost_per_run: float,
) -> list[list[PendingCircuit]]:
    submissions: list[list[PendingCircuit]] = []
    current: list[PendingCircuit] = []
    current_cost = float(BASE_SIMULATION_COST)

    for item in circuits:
        append_cost = pytket_bare_simulation_cost(item.circuit)
        if current and current_cost + append_cost > max_cost_per_run:
            submissions.append(current)
            current = []
            current_cost = float(BASE_SIMULATION_COST)
        current.append(item)
        current_cost += append_cost

    if current:
        submissions.append(current)

    debug(f"packed into {len(submissions)} unstitched execute submission(s)")
    for submission_index, submission in enumerate(submissions):
        cost = unstitched_submission_cost(submission)
        shapes = [circuit_shape(item.circuit) for item in submission]
        debug(
            f"  submission {submission_index:02d}: programs={len(submission)} "
            f"estimated_cost={cost:.6g} shapes={shapes}"
        )
    return submissions


def unstitched_submission_cost(submission: list[PendingCircuit]) -> float:
    return float(BASE_SIMULATION_COST) + sum(
        pytket_bare_simulation_cost(item.circuit)
        for item in submission
    )


def runner_batch_cost(items: list[PendingCircuit]) -> float:
    if not items:
        return 0.0
    return float(BASE_SIMULATION_COST) + sum(
        single_circuit_bare_hqc(item.config)
        for item in items
    )


def circuit_shape(circuit: PytketCircuit) -> tuple[int, int, int, int, int]:
    return (
        int(circuit.n_qubits),
        int(circuit.n_bits),
        int(circuit.n_gates),
        int(circuit.n_1qb_gates()),
        int(circuit.n_2qb_gates()),
    )


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


def result_identifier_strings(ref) -> list[str]:
    values = []
    for attr in (
        "id",
        "uuid",
        "item_id",
        "job_id",
        "job_item_id",
        "job_item_integer_id",
    ):
        value = getattr(ref, attr, None)
        if value is not None:
            values.append(str(value))
    values.append(str(ref))
    return list(dict.fromkeys(values))


def ensure_qnexus_login() -> None:
    if not is_logged_in():
        debug("QNexus login not active; opening login flow")
        qnx.login()
    else:
        debug("QNexus login already active")


def find_execute_job(project_name: str, job_id: str):
    project = qnx.projects.get(name=project_name)
    try:
        job = qnx.jobs.get(
            id=job_id,
            project=project,
            job_type=[JobType.EXECUTE],
        )
        if job is not None:
            return job
    except Exception as exc:  # noqa: BLE001
        debug(f"direct qnx.jobs.get(id=...) lookup failed: {exc}")

    for job in qnx.jobs.get_all(
        project=project,
        job_type=[JobType.EXECUTE],
        sort_filters=[SortFilterEnum.CREATED_DESC],
    ):
        if any(job_id in value for value in identifier_strings(job)):
            return job
        try:
            for ref in qnx.jobs.results(job, allow_incomplete=True):
                if any(job_id in value for value in result_identifier_strings(ref)):
                    return job
        except Exception:
            continue

    raise RuntimeError(f"No execute job {job_id!r} in project {project_name!r}.")


def list_recent_execute_jobs(
    *,
    project_name: str,
    device_name: str | None,
    n_recent_jobs: int,
) -> None:
    ensure_qnexus_login()
    project = qnx.projects.get(name=project_name)
    debug(f"recent execute jobs in project={project_name!r}")
    shown = 0
    for job in qnx.jobs.get_all(
        project=project,
        job_type=[JobType.EXECUTE],
        sort_filters=[SortFilterEnum.CREATED_DESC],
    ):
        system = None if getattr(job, "system", None) is None else job.system.name
        if device_name is not None and system != device_name:
            continue
        debug(
            f"  job id candidates={identifier_strings(job)[:3]} "
            f"status={getattr(job, 'last_status', None)} system={system!r}"
        )
        shown += 1
        if shown >= n_recent_jobs:
            break
    debug(f"listed {shown} job(s)")


def job_is_completed(job) -> bool:
    status = getattr(job, "last_status", None)
    return status == JobStatusEnum.COMPLETED or "COMPLETED" in str(status)


def download_program_results(job) -> list[tuple[PytketCircuit, BackendResult]]:
    fetched: list[tuple[PytketCircuit, BackendResult]] = []
    refs = qnx.jobs.results(job)
    if isinstance(refs, IncompleteJobItemRef):
        raise RuntimeError(f"Qnexus execution failed: {refs}")
    for ref in refs:
        if isinstance(ref, IncompleteJobItemRef):
            raise RuntimeError(f"Incomplete Qnexus job item: {ref}")
        if not isinstance(ref, ExecutionResultRef):
            continue
        result = ref.download_result()
        if not isinstance(result, BackendResult):
            continue
        input_program = ref.get_input()
        if not isinstance(input_program, CircuitRef):
            raise RuntimeError(f"Execution result input is not a CircuitRef: {input_program}")
        fetched.append((input_program.download_circuit(), result))
    return fetched


def fetch_job_results(
    *,
    project_name: str,
    job_id: str,
) -> list[tuple[PytketCircuit, BackendResult]]:
    job = find_execute_job(project_name, job_id)
    system = None if getattr(job, "system", None) is None else job.system.name
    debug(
        f"using execute job id={job_id} "
        f"status={getattr(job, 'last_status', None)} system={system!r}"
    )
    if not job_is_completed(job):
        raise RuntimeError(f"Execute job {job_id!r} is not completed.")
    fetched = download_program_results(job)
    debug(f"downloaded {len(fetched)} program result(s) for job {job_id}")
    for index, (circuit, _) in enumerate(fetched):
        debug(f"  fetched {index:02d}: shape={circuit_shape(circuit)}")
    return fetched


def match_submission_index(
    *,
    fetched: list[tuple[PytketCircuit, BackendResult]],
    submissions: list[list[PendingCircuit]],
    used_indices: set[int],
    requested_index: int | None,
    allow_circuit_mismatch: bool,
) -> int:
    fetched_shapes = [circuit_shape(circuit) for circuit, _ in fetched]

    if requested_index is not None:
        if requested_index < 0 or requested_index >= len(submissions):
            raise IndexError(
                f"--submission-index {requested_index} is outside "
                f"0..{len(submissions) - 1}."
            )
        expected_len = len(submissions[requested_index])
        if expected_len != len(fetched):
            raise RuntimeError(
                f"submission {requested_index:02d} has {expected_len} program(s), "
                f"but fetched job has {len(fetched)}."
            )
        if not allow_circuit_mismatch:
            expected_shapes = [
                circuit_shape(item.circuit)
                for item in submissions[requested_index]
            ]
            if expected_shapes != fetched_shapes:
                raise RuntimeError(
                    f"submission {requested_index:02d} circuit shapes mismatch: "
                    f"expected={expected_shapes} fetched={fetched_shapes}"
                )
        return requested_index

    strict_candidates = []
    for index, submission in enumerate(submissions):
        if index in used_indices:
            continue
        expected_shapes = [circuit_shape(item.circuit) for item in submission]
        if expected_shapes == fetched_shapes:
            strict_candidates.append(index)

    if len(strict_candidates) == 1:
        return strict_candidates[0]
    if len(strict_candidates) > 1:
        raise RuntimeError(
            "Fetched job matches multiple submissions by shape: "
            f"{strict_candidates}. Pass --submission-index."
        )

    if allow_circuit_mismatch:
        length_candidates = [
            index
            for index, submission in enumerate(submissions)
            if index not in used_indices and len(submission) == len(fetched)
        ]
        if len(length_candidates) == 1:
            debug(
                "WARNING: no shape match; using the only unused submission with "
                f"the same program count: {length_candidates[0]:02d}"
            )
            return length_candidates[0]

    raise RuntimeError(
        "Could not match fetched job to a reconstructed Fancy submission. "
        f"fetched_shapes={fetched_shapes}. Pass --submission-index if you know it, "
        "or add --allow-circuit-mismatch for a count-only fallback."
    )


def outcome_counts(result: BackendResult) -> dict[str, int]:
    counts = result.get_empirical_distribution().as_counter()
    return {
        "".join(str(int(bit)) for bit in outcome): int(count)
        for outcome, count in counts.items()
    }


def success_from_result(result: BackendResult) -> bool | None:
    counts = result.get_empirical_distribution().as_counter()
    if not counts:
        return None
    top_outcome, _ = counts.most_common(1)[0]
    return all(bit == 0 for bit in top_outcome)


def append_recovered_records(
    payload: dict,
    recovered: list[tuple[PendingCircuit, bool]],
) -> int:
    records = payload.setdefault("data", [])
    record_by_key = {record_key(record): record for record in records}
    added_points = 0

    for item, outcome in recovered:
        key = config_key(item.config)
        record = record_by_key.get(key)
        if record is None:
            record = {
                "n_1qb_gates": int(item.config.n_1qb_gates),
                "n_2qb_gates": int(item.config.n_2qb_gates),
                "n_qubits": int(item.config.n_qubits),
                "random_elimination": float(item.config.random_elimination),
                "use_scrambler": bool(item.config.use_scrambler),
                "gates_set": [gate.name for gate in item.config.gates_set],
                "results": [],
            }
            records.append(record)
            record_by_key[key] = record

        existing_results = {
            bool(old_outcome): int(count)
            for old_outcome, count in record.get("results", [])
        }
        existing_results[bool(outcome)] = existing_results.get(bool(outcome), 0) + 1
        record["results"] = [
            [old_outcome, count]
            for old_outcome, count in sorted(existing_results.items())
        ]
        added_points += 1

    return added_points


def recovered_request_metadata(
    recovered_items: list[PendingCircuit],
    *,
    pending: dict,
) -> tuple[list[dict], list[dict], list[dict]]:
    settings = pending.get("settings", {})
    rng_seed = settings.get("rng_seed")
    sent_configs = []
    proposal_configs = []
    actual_after_elimination = []

    for item in recovered_items:
        request = item.request_payload
        sent_configs.append(dict(request.get("config") or sent_config_payload(item.config)))
        if "proposal" in request:
            proposal_configs.append(request["proposal"])
        if "actual_after_elimination" in request:
            actual_after_elimination.append(request["actual_after_elimination"])
        else:
            actual_after_elimination.append(
                actual_after_elimination_payload(
                    item.config,
                    rng_seed=None if rng_seed is None else int(rng_seed),
                    shot_index=item.shot_index,
                )
            )

    return sent_configs, proposal_configs, actual_after_elimination


def update_experiment_metadata(
    payload: dict,
    *,
    pending: dict,
    previous_measurement: Path,
    recovered_items: list[PendingCircuit],
    job_ids: list[str],
    matched_indices: list[int],
    project_name: str,
    device_name: str,
    complete: bool,
    recovered_points: int,
    expected_points: int,
) -> None:
    previous_experiment = dict(payload.get("experiment", {}))
    settings = pending.get("settings", {})
    hqc_budget = float(
        previous_experiment.get(
            "hqc_budget",
            settings.get("hqc_budget", 0.0),
        )
    )
    previous_spent = float(previous_experiment.get("spent_hqc", 0.0))
    recovered_batch_cost = (
        runner_batch_cost(recovered_items)
        if not complete
        else runner_batch_cost(recovered_items)
    )
    spent_hqc = previous_spent + recovered_batch_cost

    sent_configs, proposal_configs, actual_after_elimination = recovered_request_metadata(
        recovered_items,
        pending=pending,
    )

    experiment = {
        "phase": str(pending.get("phase", "globalsur")),
        "step": int(pending.get("step", previous_experiment.get("step", 0) + 1)),
        "spent_hqc": float(spent_hqc),
        "remaining_hqc": float(max(hqc_budget - spent_hqc, 0.0)),
        "hqc_budget": float(hqc_budget),
        "repair_stats": previous_experiment.get("repair_stats", {}),
        "sent_configs": sent_configs,
        "actual_after_elimination": actual_after_elimination,
        "recovery": {
            "source": "fancy_emulator_pending_backend_batch_and_qnexus_job",
            "job_ids": job_ids,
            "project_name": project_name,
            "device_name": device_name,
            "matched_submission_indices": matched_indices,
            "pending_backend_batch": pending.get("_path"),
            "previous_measurement": str(previous_measurement),
            "complete_pending_batch": bool(complete),
            "recovered_batch_cost": float(recovered_batch_cost),
            "recovered_points": int(recovered_points),
            "expected_points": int(expected_points),
        },
    }
    if proposal_configs:
        experiment["proposal_configs"] = proposal_configs
    for key in (
        "circuit_variant",
        "protected_h_wrapper",
        "h_wrapper_1q_gates_per_qubit",
    ):
        if key in settings:
            experiment[key] = settings[key]

    payload["experiment"] = experiment


def recover(args: argparse.Namespace) -> Path:
    run_folder = args.run_folder
    pending = load_pending(run_folder)
    requests = pending_requests(pending)
    previous_measurement = args.previous_measurement or latest_measurement_checkpoint(
        run_folder
    )
    previous_payload = read_json(previous_measurement)
    device_name = args.device_name or infer_device_name(pending, run_folder)
    project_name = args.project_name or infer_project_name(
        pending,
        run_folder,
        device_name,
    )
    max_cost_per_run = float(
        args.max_cost_per_run
        if args.max_cost_per_run is not None
        else pending.get("settings", {}).get("max_cost_per_run", 30.0)
    )
    job_ids = parse_job_ids(args.job_id)
    submission_indices = args.submission_index or []

    debug(f"run folder={run_folder}")
    debug(f"pending batch={pending['_path']}")
    debug(f"previous measurement={previous_measurement}")
    debug(f"pending model={pending.get('model')}")
    debug(f"pending phase={pending.get('phase')} step={pending.get('step')}")
    debug(f"backend_factory={backend_factory_name(pending)}")
    debug(f"project={project_name}")
    debug(f"device={device_name}")
    debug(f"job ids={job_ids}")
    debug(f"max_cost_per_run={max_cost_per_run}")

    if not job_ids:
        raise ValueError("No --job-id supplied.")
    if submission_indices and len(submission_indices) not in (1, len(job_ids)):
        raise ValueError(
            "Pass either one --submission-index for one job id, or one "
            "--submission-index per job id."
        )

    built_circuits = build_pending_circuits(
        pending=pending,
        requests=requests,
        previous_payload=previous_payload,
        run_folder=run_folder,
    )
    submissions = pack_unstitched_submissions(
        built_circuits,
        max_cost_per_run=max_cost_per_run,
    )

    ensure_qnexus_login()

    recovered: list[tuple[PendingCircuit, bool]] = []
    recovered_details = []
    matched_indices: list[int] = []
    used_indices: set[int] = set()

    for job_index, job_id in enumerate(job_ids):
        fetched = fetch_job_results(project_name=project_name, job_id=job_id)
        requested_index = None
        if submission_indices:
            requested_index = (
                submission_indices[0]
                if len(submission_indices) == 1
                else submission_indices[job_index]
            )
        submission_index = match_submission_index(
            fetched=fetched,
            submissions=submissions,
            used_indices=used_indices,
            requested_index=requested_index,
            allow_circuit_mismatch=args.allow_circuit_mismatch,
        )
        used_indices.add(submission_index)
        matched_indices.append(submission_index)
        submission = submissions[submission_index]
        debug(
            f"recovering job {job_id} as submission {submission_index:02d} "
            f"with {len(submission)} program(s)"
        )

        for result_index, (item, (_, result)) in enumerate(zip(submission, fetched)):
            success = success_from_result(result)
            counts = outcome_counts(result)
            if success is None:
                debug(
                    f"  result {result_index:02d}: request={item.request_index:02d} "
                    "empty counts, skipped"
                )
                continue
            recovered.append((item, bool(success)))
            recovered_details.append(
                {
                    "job_id": job_id,
                    "submission_index": int(submission_index),
                    "program_index": int(result_index),
                    "request_index": int(item.request_index),
                    "shot_index": int(item.shot_index),
                    "success": bool(success),
                    "counts": counts,
                    "n_1qb_gates": int(item.config.n_1qb_gates),
                    "n_2qb_gates": int(item.config.n_2qb_gates),
                    "n_qubits": int(item.config.n_qubits),
                }
            )
            debug(
                f"  result {result_index:02d}: request={item.request_index:02d} "
                f"shot_index={item.shot_index} q={item.config.n_qubits} "
                f"n_1q={item.config.n_1qb_gates} n_2q={item.config.n_2qb_gates} "
                f"success={bool(success)} counts={counts}"
            )

    recovered_points = len(recovered)
    expected_points = sum(int(request.get("requested_shots", 1)) for request in requests)
    complete = recovered_points == expected_points
    debug(f"recovered points={recovered_points}")
    debug(f"expected points={expected_points}")
    if not complete and not args.allow_partial:
        raise RuntimeError(
            f"Recovered {recovered_points}/{expected_points} pending point(s). "
            "Pass all execute job ids for this pending batch, or use "
            "--allow-partial if you deliberately want a partial recovery JSON."
        )

    added_points = append_recovered_records(previous_payload, recovered)
    recovered_items = [item for item, _ in recovered]
    update_experiment_metadata(
        previous_payload,
        pending=pending,
        previous_measurement=previous_measurement,
        recovered_items=recovered_items,
        job_ids=job_ids,
        matched_indices=matched_indices,
        project_name=project_name,
        device_name=device_name,
        complete=complete,
        recovered_points=recovered_points,
        expected_points=expected_points,
    )
    previous_payload["experiment"]["recovery"]["result_details"] = recovered_details

    output_path = args.output or output_measurement_path(
        run_folder,
        pending,
        partial=not complete,
    )
    debug(f"records after append={len(previous_payload.get('data', []))}")
    debug(f"added points={added_points}")
    if args.dry_run:
        debug(f"dry run: would write {output_path}")
        return output_path

    write_json(output_path, previous_payload)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-folder", type=Path, required=True)
    parser.add_argument(
        "--job-id",
        action="append",
        default=[],
        help="Execute job id. Repeat or comma-separate if the pending batch split.",
    )
    parser.add_argument("--project-name", default=None)
    parser.add_argument("--device-name", default=None)
    parser.add_argument("--previous-measurement", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--max-cost-per-run", type=float, default=None)
    parser.add_argument(
        "--submission-index",
        type=int,
        action="append",
        default=[],
        help="Reconstructed submission index for a job id, if auto-matching is ambiguous.",
    )
    parser.add_argument("--allow-circuit-mismatch", action="store_true")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--list-jobs", action="store_true")
    parser.add_argument("--n-recent-jobs", type=int, default=10)
    args = parser.parse_args()

    if args.list_jobs:
        pending = load_pending(args.run_folder)
        device_name = args.device_name or infer_device_name(pending, args.run_folder)
        project_name = args.project_name or infer_project_name(
            pending,
            args.run_folder,
            device_name,
        )
        list_recent_execute_jobs(
            project_name=project_name,
            device_name=device_name,
            n_recent_jobs=args.n_recent_jobs,
        )
        raise SystemExit(0)

    return args


def main() -> None:
    recover(parse_args())


if __name__ == "__main__":
    main()
