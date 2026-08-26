"""Recover one pending H-wrapper FLE batch from a completed QNexus job.

This script does not submit circuits. It reads the run folder's
``pending_backend_batch.json``, reconstructs the H-wrapped stitched programs,
downloads the completed QNexus execute job, destitches the results, appends the
recovered outcomes to the latest measurement checkpoint, and writes a new
``measurement_XXX_*.json`` in the same run folder.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
import qnexus as qnx
from pytket.backends.backendresult import BackendResult
from pytket.circuit import BitRegister, Circuit as PytketCircuit
from qnexus.client.auth import is_logged_in
from qnexus.models.filters import SortFilterEnum
from qnexus.models.job_status import JobStatusEnum
from qnexus.models.references import CircuitRef, ExecutionResultRef, JobType

from sympleq.applications.randomized_benchmarking.backends.quantinuum import (
    QuantinuumBackend,
)
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    measurement_rng,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.fantasy_levelset_settings import (
    FantasySettings,
    control_panel_settings_kwargs,
)
from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.h_wrapper_config import (
    HWrappedRMBConfig,
)
from sympleq.integrations.quantinuum.stitching import (
    MAX_QASM_PROGRAM_SIZE,
    circuit_stitching,
    destitch_results,
    estimate_qasm_program_size,
)
from sympleq.integrations.quantinuum.utils import (
    BASE_SIMULATION_COST,
    NATIVE_GATES_SET,
    pytket_bare_simulation_cost,
)


DEFAULT_RUN_FOLDER = Path(
    r"Personal\FLE\H_wrapper\H2_1\q56\seed_42\FLE_H_wrapper_20260819_134349"
)
DEFAULT_PROJECT_NAME = "FLE-benchmark-h2-1-Hwrap"
DEFAULT_DEVICE_NAME = "H2-1"


def debug(message: str) -> None:
    print(f"[hwrap-recovery] {message}", flush=True)


def settings_from_pending(pending_payload: dict) -> FantasySettings:
    kwargs = control_panel_settings_kwargs()
    snapshot = pending_payload.get("settings", {})
    if "rng_seed" in snapshot:
        kwargs["rng_seed"] = int(snapshot["rng_seed"])
    if "n_qubits" in snapshot:
        kwargs["n_qubits"] = int(snapshot["n_qubits"])
    if "n_gates_bounds" in snapshot:
        kwargs["n_gates_bounds"] = tuple(float(v) for v in snapshot["n_gates_bounds"])
    if "ratio_bounds" in snapshot:
        kwargs["ratio_bounds"] = tuple(float(v) for v in snapshot["ratio_bounds"])
    if "hqc_budget" in snapshot:
        kwargs["hqc_budget"] = float(snapshot["hqc_budget"])
    if "max_cost_per_run" in snapshot:
        kwargs["max_cost_per_run"] = float(snapshot["max_cost_per_run"])
    return FantasySettings(**kwargs)


def ensure_qnexus_login() -> None:
    if not is_logged_in():
        qnx.login()


def register_index(register: BitRegister) -> int:
    name = register.name
    if name.startswith("creg_"):
        return int(name.removeprefix("creg_"))
    return 10**9


def config_from_record(record: dict) -> HWrappedRMBConfig:
    return (
        HWrappedRMBConfig.default()
        .with_n_qubits(int(record["n_qubits"]))
        .with_n_1qb_gates(int(record["n_1qb_gates"]))
        .with_n_2qb_gates(int(record["n_2qb_gates"]))
        .with_random_elimination(float(record.get("random_elimination", 0.1)))
        .with_use_scrambler(bool(record.get("use_scrambler", True)))
        .with_gates_set(tuple(NATIVE_GATES_SET))
    )


def config_key(config: RMBConfig) -> tuple[int, int, int, float, bool]:
    return (
        int(config.n_1qb_gates),
        int(config.n_2qb_gates),
        int(config.n_qubits),
        float(config.random_elimination),
        bool(config.use_scrambler),
    )


def config_key_from_record(record: dict) -> tuple[int, int, int, float, bool]:
    return (
        int(record["n_1qb_gates"]),
        int(record["n_2qb_gates"]),
        int(record["n_qubits"]),
        float(record.get("random_elimination", 0.1)),
        bool(record.get("use_scrambler", True)),
    )


def pending_payload(run_folder: Path) -> dict:
    path = run_folder / "pending_backend_batch.json"
    if not path.exists():
        raise FileNotFoundError(f"No pending backend batch found: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
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


def latest_measurement_checkpoint(run_folder: Path) -> Path:
    checkpoints = sorted(
        path
        for path in run_folder.glob("measurement_*.json")
        if "_actual_gates" not in path.stem
    )
    if not checkpoints:
        raise FileNotFoundError(f"No measurement_*.json checkpoints in {run_folder}")
    return checkpoints[-1]


def output_measurement_path(run_folder: Path, pending: dict) -> Path:
    step = int(pending.get("step", 0))
    if step <= 0:
        step = 1
    phase = str(pending.get("phase", "globalsur"))
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    return run_folder / f"measurement_{step:03d}_{phase}_{timestamp}.json"


def append_recovered_records(
    payload: dict,
    recovered_by_config: dict[RMBConfig, dict[bool, int]],
) -> int:
    records = payload.setdefault("data", [])
    record_by_key = {config_key_from_record(record): record for record in records}
    added_points = 0

    for config, counts in recovered_by_config.items():
        key = config_key(config)
        record = record_by_key.get(key)
        if record is None:
            record = {
                "n_1qb_gates": int(config.n_1qb_gates),
                "n_2qb_gates": int(config.n_2qb_gates),
                "n_qubits": int(config.n_qubits),
                "random_elimination": float(config.random_elimination),
                "use_scrambler": bool(config.use_scrambler),
                "gates_set": [gate.name for gate in config.gates_set],
                "results": [],
            }
            records.append(record)
            record_by_key[key] = record

        existing_results = {
            bool(outcome): int(count)
            for outcome, count in record.get("results", [])
        }
        for outcome, count in counts.items():
            existing_results[bool(outcome)] = (
                existing_results.get(bool(outcome), 0) + int(count)
            )
            added_points += int(count)
        record["results"] = [[outcome, count] for outcome, count in existing_results.items()]

    return added_points


def build_pending_jobs(
    settings: FantasySettings,
    requests: list[dict],
) -> list[tuple[RMBConfig, PytketCircuit]]:
    if settings.rng_seed is None:
        raise ValueError("Cannot recover without a deterministic rng_seed.")

    jobs = []
    for request_index, request in enumerate(requests):
        config = config_from_record(request["config"])
        first_shot_index = int(request.get("first_shot_index", 0))
        shots = int(request.get("requested_shots", 1))
        debug(
            "pending request "
            f"{request_index:02d}: shots={shots} first_shot_index={first_shot_index} "
            f"n_1q={config.n_1qb_gates} n_2q={config.n_2qb_gates} "
            f"q={config.n_qubits} ratio={config.ratio_2_qb_gates:.6f}"
        )
        for shot_offset in range(shots):
            shot_index = first_shot_index + shot_offset
            circuit_rng = measurement_rng(int(settings.rng_seed), config, shot_index)
            circuit = QuantinuumBackend.compatible_circuit(config, circuit_rng)
            jobs.append((config, circuit))
            debug(
                "  built circuit "
                f"shot_index={shot_index} qubits={circuit.n_qubits} "
                f"bits={circuit.n_bits} gates={circuit.n_gates}"
            )
    return jobs


def pack_jobs(
    settings: FantasySettings,
    jobs: list[tuple[RMBConfig, PytketCircuit]],
) -> list[list[tuple[RMBConfig, PytketCircuit]]]:
    submissions = []
    current = []
    current_size = 0
    current_cost = float(BASE_SIMULATION_COST)

    for config, circuit in jobs:
        program_size = estimate_qasm_program_size(circuit)
        append_cost = pytket_bare_simulation_cost(circuit)
        if current:
            append_cost += config.n_qubits / 5000
            if (
                current_size + program_size > MAX_QASM_PROGRAM_SIZE
                or current_cost + append_cost > settings.max_cost_per_run
            ):
                submissions.append(current)
                current = []
                current_size = 0
                current_cost = float(BASE_SIMULATION_COST)
                append_cost = pytket_bare_simulation_cost(circuit)
        current.append((config, circuit))
        current_size += program_size
        current_cost += append_cost

    if current:
        submissions.append(current)

    debug(f"packed into {len(submissions)} stitched submission(s)")
    return submissions


def submission_cost(submission: list[tuple[RMBConfig, PytketCircuit]]) -> float:
    if not submission:
        return 0.0
    cost = float(BASE_SIMULATION_COST)
    for index, (config, circuit) in enumerate(submission):
        if index:
            cost += config.n_qubits / 5000
        cost += pytket_bare_simulation_cost(circuit)
    return float(cost)


def job_identifier_strings(job) -> list[str]:
    values = []
    for attr in ("id", "uuid", "item_id", "job_id"):
        value = getattr(job, attr, None)
        if value is not None:
            values.append(str(value))
    for attr in ("ref", "annotations"):
        value = getattr(job, attr, None)
        if value is not None:
            for nested in ("id", "uuid", "item_id", "job_id"):
                nested_value = getattr(value, nested, None)
                if nested_value is not None:
                    values.append(str(nested_value))
    values.append(str(job))
    return list(dict.fromkeys(values))


def result_identifier_strings(ref) -> list[str]:
    values = []
    for attr in ("id", "uuid", "item_id", "job_id", "job_item_id", "job_item_integer_id"):
        value = getattr(ref, attr, None)
        if value is not None:
            values.append(str(value))
    values.append(str(ref))
    return list(dict.fromkeys(values))


def find_execute_job_by_id(project_name: str, job_id: str):
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

    job_iter = qnx.jobs.get_all(
        project=project,
        job_type=[JobType.EXECUTE],
        sort_filters=[SortFilterEnum.CREATED_DESC],
    )
    for job in job_iter:
        if any(job_id in value for value in job_identifier_strings(job)):
            return job
        try:
            for ref in qnx.jobs.results(job, allow_incomplete=True):
                if any(job_id in value for value in result_identifier_strings(ref)):
                    return job
        except Exception:
            continue
    raise RuntimeError(f"No execute job matching id {job_id!r} in {project_name!r}.")


def download_program_results(job) -> list[tuple[PytketCircuit, BackendResult]]:
    fetched = []
    for ref in qnx.jobs.results(job):
        if not isinstance(ref, ExecutionResultRef):
            continue
        result = ref.download_result()
        if not isinstance(result, BackendResult):
            continue
        input_program = ref.get_input()
        if not isinstance(input_program, CircuitRef):
            continue
        fetched.append((input_program.download_circuit(), result))
    return fetched


def fetch_exact_job(
    *,
    project_name: str,
    job_id: str,
) -> list[tuple[PytketCircuit, BackendResult]]:
    ensure_qnexus_login()
    job = find_execute_job_by_id(project_name, job_id)
    system = None if getattr(job, "system", None) is None else job.system.name
    debug(
        f"using execute job id={job_id} "
        f"status={getattr(job, 'last_status', None)} system={system!r}"
    )
    if getattr(job, "last_status", None) != JobStatusEnum.COMPLETED:
        raise RuntimeError(f"Execute job {job_id!r} is not completed.")
    fetched = download_program_results(job)
    debug(f"downloaded {len(fetched)} stitched program result(s)")
    return fetched


def most_common_zero_outcome(result: BackendResult) -> bool | None:
    counts = result.get_empirical_distribution().as_counter()
    if not counts:
        return None
    outcome, _ = counts.most_common(1)[0]
    return all(bit == 0 for bit in outcome)


def recover_outcomes(
    submissions: list[list[tuple[RMBConfig, PytketCircuit]]],
    fetched: list[tuple[PytketCircuit, BackendResult]],
    *,
    allow_circuit_mismatch: bool,
) -> dict[RMBConfig, dict[bool, int]]:
    stitched_circuits = [
        circuit_stitching([circuit for _, circuit in submission])
        for submission in submissions
    ]

    debug("expected stitched program shapes:")
    for index, circuit in enumerate(stitched_circuits):
        debug(
            f"  expected {index:02d}: qubits={circuit.n_qubits} "
            f"bits={circuit.n_bits} gates={circuit.n_gates}"
        )
    debug("fetched stitched program shapes:")
    for index, (circuit, _) in enumerate(fetched):
        debug(
            f"  fetched  {index:02d}: qubits={circuit.n_qubits} "
            f"bits={circuit.n_bits} gates={circuit.n_gates}"
        )

    if len(fetched) != len(stitched_circuits):
        raise RuntimeError(
            f"Fetched {len(fetched)} program(s), expected {len(stitched_circuits)}."
        )

    for index, (expected, (actual, _)) in enumerate(zip(stitched_circuits, fetched)):
        if actual.n_qubits != expected.n_qubits or actual.n_gates != expected.n_gates:
            message = (
                f"program {index} mismatch: "
                f"actual=({actual.n_qubits}, {actual.n_gates}) "
                f"expected=({expected.n_qubits}, {expected.n_gates})"
            )
            if not allow_circuit_mismatch:
                raise RuntimeError(message)
            debug(f"WARNING: {message}")

    recovered: dict[RMBConfig, dict[bool, int]] = defaultdict(lambda: defaultdict(int))
    for submission_index, (submission, (_, result), stitched) in enumerate(
        zip(submissions, fetched, stitched_circuits)
    ):
        registers = sorted(stitched.c_registers, key=register_index)
        sorted_submission = sorted(
            submission,
            key=lambda item: item[1].n_qubits,
            reverse=True,
        )
        debug(
            f"destitching submission {submission_index:02d}: "
            f"registers={[register.name for register in registers]}"
        )
        for result_index, ((config, _), sub_result) in enumerate(
            zip(sorted_submission, destitch_results(result, registers))
        ):
            outcome = most_common_zero_outcome(sub_result)
            if outcome is None:
                debug(f"  result {result_index:02d}: empty counts, skipped")
                continue
            recovered[config][bool(outcome)] += 1
            debug(
                f"  result {result_index:02d}: q={config.n_qubits} "
                f"n_1q={config.n_1qb_gates} n_2q={config.n_2qb_gates} "
                f"success={bool(outcome)}"
            )
    return recovered


def recover_pending_batch(args: argparse.Namespace) -> Path:
    run_folder = args.run_folder
    pending = pending_payload(run_folder)
    requests = pending_requests(pending)
    settings = settings_from_pending(pending)
    previous_path = args.previous_measurement or latest_measurement_checkpoint(run_folder)
    previous_payload = json.loads(previous_path.read_text(encoding="utf-8"))

    debug(f"run folder={run_folder}")
    debug(f"pending batch={pending['_path']}")
    debug(f"previous measurement={previous_path}")
    debug(f"job id={args.job_id}")
    debug(f"project={args.project_name}")
    debug(f"device={args.device_name}")
    debug(f"rng_seed={settings.rng_seed}")
    debug(f"max_cost_per_run={settings.max_cost_per_run}")

    jobs = build_pending_jobs(settings, requests)
    submissions = pack_jobs(settings, jobs)
    batch_cost = sum(submission_cost(submission) for submission in submissions)
    debug(f"expected recovered batch cost={batch_cost:.6g} HQC")

    fetched = fetch_exact_job(project_name=args.project_name, job_id=args.job_id)
    recovered = recover_outcomes(
        submissions,
        fetched,
        allow_circuit_mismatch=args.allow_circuit_mismatch,
    )
    recovered_points = sum(sum(counts.values()) for counts in recovered.values())
    expected_points = sum(int(request.get("requested_shots", 1)) for request in requests)
    debug(f"recovered points={recovered_points}")
    debug(f"expected points={expected_points}")
    if recovered_points != expected_points:
        raise RuntimeError(
            f"Recovered {recovered_points} points, expected {expected_points}."
        )

    added_points = append_recovered_records(previous_payload, recovered)
    previous_experiment = previous_payload.get("experiment", {})
    previous_spent = float(previous_experiment.get("spent_hqc", 0.0))
    hqc_budget = float(previous_experiment.get("hqc_budget", settings.hqc_budget))
    spent_hqc = previous_spent + batch_cost

    previous_payload["experiment"] = {
        "phase": str(pending.get("phase", "globalsur")),
        "step": int(pending.get("step", previous_experiment.get("step", 0) + 1)),
        "spent_hqc": float(spent_hqc),
        "remaining_hqc": float(max(hqc_budget - spent_hqc, 0.0)),
        "hqc_budget": float(hqc_budget),
        "repair_stats": previous_experiment.get("repair_stats", {}),
        "circuit_variant": pending.get("settings", {}).get("circuit_variant"),
        "protected_h_wrapper": pending.get("settings", {}).get("protected_h_wrapper"),
        "h_wrapper_1q_gates_per_qubit": pending.get("settings", {}).get(
            "h_wrapper_1q_gates_per_qubit"
        ),
        "sent_configs": [request["config"] for request in requests],
        "proposal_configs": [request.get("proposal") for request in requests],
        "actual_after_elimination": [
            request.get("actual_after_elimination") for request in requests
        ],
        "recovery": {
            "source": "h_wrapper_pending_backend_batch_and_qnexus_job",
            "job_id": args.job_id,
            "project_name": args.project_name,
            "device_name": args.device_name,
            "pending_backend_batch": pending.get("_path"),
            "previous_measurement": str(previous_path),
            "recovered_batch_cost": float(batch_cost),
            "recovered_points": int(recovered_points),
        },
    }

    output_path = args.output or output_measurement_path(run_folder, pending)
    debug(f"records after append={len(previous_payload.get('data', []))}")
    debug(f"added points={added_points}")
    if args.dry_run:
        debug(f"dry run: would write {output_path}")
        return output_path

    output_path.write_text(json.dumps(previous_payload, indent=2), encoding="utf-8")
    debug(f"wrote recovered measurement JSON: {output_path}")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--run-folder", type=Path, default=DEFAULT_RUN_FOLDER)
    parser.add_argument("--project-name", default=DEFAULT_PROJECT_NAME)
    parser.add_argument("--device-name", default=DEFAULT_DEVICE_NAME)
    parser.add_argument("--previous-measurement", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--allow-circuit-mismatch", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    recover_pending_batch(args)


if __name__ == "__main__":
    main()
