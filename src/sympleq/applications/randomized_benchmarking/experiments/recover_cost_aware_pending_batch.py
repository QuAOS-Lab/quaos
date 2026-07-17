"""Recover a COST_AWARE batch that was submitted before a local disconnect.

The COST_AWARE loop writes a checkpoint with ``reason="before_backend"`` and a
``pending_batch`` before calling the remote Quantinuum backend.  If the Python
process loses its connection after the job was submitted, this script fetches a
completed QNexus execute job, destitches the results, records them against the
pending configs, appends the completed batch to checkpoint metadata, and clears
``pending_batch``.  After this, rerun ``common_run_models.py COST_AWARE`` and the
normal resume path will continue from the recovered data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import qnexus as qnx
from pytket.backends.backendresult import BackendResult
from qnexus.models.filters import SortFilterEnum
from qnexus.models.job_status import JobStatusEnum
from qnexus.models.references import CircuitRef, ExecutionResultRef, JobType

from sympleq.applications.randomized_benchmarking.RMB import RMB, resolve_data_path
from sympleq.applications.randomized_benchmarking.backends.base import MeasurementRequest
from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.applications.randomized_benchmarking.experiments.common import (
    BASE_SIMULATION_COST,
    measurement_rng,
    pytket_bare_simulation_cost,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_settings import (
    CostAwareSettings,
    backend_factory_for_name,
    control_panel_settings_kwargs,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import (
    _batch_checkpoint_record,
    _config_from_checkpoint_record,
    _measured_config_count,
    _write_json_atomic,
)
from sympleq.integrations.quantinuum.stitching import (
    MAX_QASM_PROGRAM_SIZE,
    circuit_stitching,
    destitch_results,
    estimate_qasm_program_size,
)
from sympleq.integrations.quantinuum.utils import fetch_recent_execute_jobs


GENERATED_PERSONAL_ROOT = (
    Path("scripts")
    / "personal"
    / "randomized_benchmarking_personal"
    / "Personal"
)


def _data_path_from_checkpoint(checkpoint_path: Path) -> Path:
    name = checkpoint_path.name
    if not name.endswith("_checkpoint.json"):
        raise ValueError(
            "Checkpoint path should end with '_checkpoint.json'; pass --data-path "
            "if this file uses a different naming convention."
        )
    return checkpoint_path.with_name(name.removesuffix("_checkpoint.json") + ".json")


def _requests_from_pending(meta: dict) -> tuple[list[MeasurementRequest], dict[RMBConfig, int]]:
    pending = meta.get("pending_batch")
    if not isinstance(pending, dict):
        raise ValueError("Checkpoint has no pending_batch to recover.")

    requests: list[MeasurementRequest] = []
    first_shot_indices: dict[RMBConfig, int] = {}
    for record in pending.get("requests", []):
        config_record = record.get("config")
        if not isinstance(config_record, dict):
            continue
        shots = int(record.get("requested_shots", 0))
        if shots <= 0:
            continue
        config = _config_from_checkpoint_record(config_record)
        requests.append(MeasurementRequest(config, shots))
        first_shot_indices[config] = int(record.get("first_shot_index", 0))
    if not requests:
        raise ValueError("pending_batch contains no positive-shot requests.")
    return requests, first_shot_indices


def _reconstruct_plain_quantinuum_submissions(
    backend,
    requests: list[MeasurementRequest],
    *,
    seed: int,
    first_shot_indices: dict[RMBConfig, int],
):
    """Rebuild the same submissions as QuantinuumBackend.fidelity_estimation.

    This supports the plain hardware path (H2-1/H2-2).  Gate-limited emulator
    targets have a separate packing rule and should not be recovered with this
    script unless that path is added explicitly.
    """
    jobs = []
    shot_counts: dict[RMBConfig, int] = {}
    for request in requests:
        for _ in range(max(0, int(request.shots))):
            index = shot_counts.get(request.config, 0)
            shot_counts[request.config] = index + 1
            circuit_rng = measurement_rng(
                seed,
                request.config,
                first_shot_indices[request.config] + index,
            )
            jobs.append(
                (
                    request.config,
                    backend.compatible_circuit(request.config, circuit_rng),
                )
            )
    if not jobs:
        return [], 0.0

    submissions = []
    current = []
    current_size = 0
    current_cost = float(BASE_SIMULATION_COST)
    total_cost = 0.0
    for config, circuit in jobs:
        program_size = estimate_qasm_program_size(circuit)
        append_cost = pytket_bare_simulation_cost(circuit)
        if current:
            append_cost += config.n_qubits / 5000
            if (
                current_size + program_size > MAX_QASM_PROGRAM_SIZE
                or current_cost + append_cost > backend.max_cost_per_run
            ):
                submissions.append(current)
                total_cost += current_cost
                current = []
                current_size = 0
                current_cost = float(BASE_SIMULATION_COST)
                append_cost = pytket_bare_simulation_cost(circuit)
        current.append((config, circuit))
        current_size += program_size
        current_cost += append_cost
    submissions.append(current)
    total_cost += current_cost
    return submissions, total_cost


def _job_identifier_strings(job) -> list[str]:
    """Best-effort extraction of stable ID strings from a QNexus job object."""
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


def _result_identifier_strings(ref) -> list[str]:
    values = []
    for attr in ("id", "uuid", "item_id", "job_id", "job_item_id", "job_item_integer_id"):
        value = getattr(ref, attr, None)
        if value is not None:
            values.append(str(value))
    values.append(str(ref))
    return list(dict.fromkeys(values))


def _find_execute_job_by_id(project_name: str, job_id: str):
    project = qnx.projects.get(name=project_name)
    try:
        job = qnx.jobs.get(
            id=job_id,
            project=project,
            job_type=[JobType.EXECUTE],
        )
        if job is not None:
            return job
    except Exception as exc:  # noqa: BLE001 - fall back to scanning.
        print(f"Direct qnx.jobs.get(id=...) lookup failed: {exc}")

    job_iter = qnx.jobs.get_all(
        project=project,
        job_type=[JobType.EXECUTE],
        sort_filters=[SortFilterEnum.CREATED_DESC],
    )
    for job in job_iter:
        if any(job_id in value for value in _job_identifier_strings(job)):
            return job
        try:
            for ref in qnx.jobs.results(job, allow_incomplete=True):
                if any(job_id in value for value in _result_identifier_strings(ref)):
                    return job
        except Exception:
            continue
    raise RuntimeError(f"No execute job matching id {job_id!r} found in project {project_name!r}.")


def _download_program_results_from_job(job):
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


def _recover_outcomes_from_fetched_programs(
    submissions,
    fetched,
    *,
    allow_circuit_mismatch: bool,
):
    stitched_circuits = [
        circuit_stitching([circuit for _, circuit in submission])
        for submission in submissions
    ]
    if len(fetched) < len(stitched_circuits):
        raise RuntimeError(
            f"Fetched {len(fetched)} result program(s), but expected "
            f"{len(stitched_circuits)} stitched program(s)."
        )
    fetched = fetched[: len(stitched_circuits)]

    for i, (expected, (actual, _)) in enumerate(zip(stitched_circuits, fetched)):
        if actual.n_qubits != expected.n_qubits or actual.n_gates != expected.n_gates:
            message = (
                f"Fetched program {i} does not match reconstructed pending program: "
                f"actual qubits/gates=({actual.n_qubits}, {actual.n_gates}), "
                f"expected=({expected.n_qubits}, {expected.n_gates})."
            )
            if not allow_circuit_mismatch:
                raise RuntimeError(message)
            print(f"WARNING: {message}")

    outcomes: dict[RMBConfig, list[bool]] = {}
    for submission, (_, result), stitched in zip(submissions, fetched, stitched_circuits):
        registers = sorted(
            stitched.c_registers,
            key=lambda register: int(register.name.removeprefix("creg_")),
        )
        sorted_submission = sorted(
            submission,
            key=lambda item: item[1].n_qubits,
            reverse=True,
        )
        for (config, _), sub_result in zip(
            sorted_submission,
            destitch_results(result, registers),
        ):
            counts = sub_result.get_empirical_distribution().as_counter()
            if not counts:
                continue
            top_outcome, _ = counts.most_common(1)[0]
            outcomes.setdefault(config, []).append(all(bit == 0 for bit in top_outcome))
    return outcomes


def _recover_outcomes_from_recent_or_exact_job(
    backend,
    submissions,
    *,
    job_id: str | None,
    n_recent_jobs: int,
    device_name: str | None,
    allow_circuit_mismatch: bool,
):
    if job_id is None:
        fetched = list(
            fetch_recent_execute_jobs(
                backend.project_name,
                n_recent_jobs,
                device_name=device_name,
            )
        )
    else:
        job = _find_execute_job_by_id(backend.project_name, job_id)
        system = None if getattr(job, "system", None) is None else job.system.name
        print(
            f"Using execute job id={job_id} status={getattr(job, 'last_status', None)} "
            f"system={system!r}"
        )
        if getattr(job, "last_status", None) != JobStatusEnum.COMPLETED:
            raise RuntimeError(f"Execute job {job_id!r} is not completed.")
        fetched = _download_program_results_from_job(job)

    return _recover_outcomes_from_fetched_programs(
        submissions,
        fetched,
        allow_circuit_mismatch=allow_circuit_mismatch,
    )


def list_recent_execute_jobs(
    *,
    project_name: str,
    n_recent_jobs: int,
) -> None:
    project = qnx.projects.get(name=project_name)
    job_iter = qnx.jobs.get_all(
        project=project,
        job_type=[JobType.EXECUTE],
        sort_filters=[SortFilterEnum.CREATED_DESC],
    )
    print(f"Recent execute jobs in project {project_name!r}:")
    for i, job in enumerate(job_iter):
        if i >= n_recent_jobs:
            break
        system = None if job.system is None else job.system.name
        identifiers = [
            value
            for value in _job_identifier_strings(job)
            if "object at 0x" not in value
        ]
        result_ids = []
        try:
            for ref in qnx.jobs.results(job, allow_incomplete=True):
                result_ids.extend(
                    value
                    for value in _result_identifier_strings(ref)
                    if "object at 0x" not in value
                )
        except Exception:
            pass
        print(
            f"{i:02d}: status={job.last_status} system={system!r} "
            f"name={getattr(job.annotations, 'name', None)!r} "
            f"created={getattr(job, 'created', None)!r} "
            f"ids={identifiers[:3]!r} result_ids={result_ids[:5]!r}"
        )


def recover_pending_batch(
    checkpoint_path: Path,
    data_path: Path | None,
    *,
    job_id: str | None,
    n_recent_jobs: int,
    project_name: str | None,
    device_name: str | None,
    no_device_filter: bool,
    list_jobs: bool,
    allow_circuit_mismatch: bool,
    dry_run: bool,
) -> None:
    checkpoint_path = checkpoint_path.resolve()
    data_path = (
        _data_path_from_checkpoint(checkpoint_path)
        if data_path is None
        else data_path.resolve()
    )
    meta = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if meta.get("reason") != "before_backend":
        raise ValueError(
            f"Checkpoint reason is {meta.get('reason')!r}, not 'before_backend'."
        )

    requests, first_shot_indices = _requests_from_pending(meta)
    backend_model = str(meta.get("settings", {}).get("backend_model", ""))
    seed = meta.get("settings", {}).get("rng_seed")
    if seed is None:
        raise ValueError("Cannot recover deterministic circuits because rng_seed is missing.")
    if backend_model.endswith("E") or backend_model == "emulator":
        raise NotImplementedError(
            "This recovery script currently supports plain H2-1/H2-2 packing. "
            "The gate-limited emulator packing path needs separate support."
        )

    kwargs = control_panel_settings_kwargs()
    kwargs["backend_model"] = backend_model
    kwargs["backend_factory"] = backend_factory_for_name(backend_model)
    kwargs["rng_seed"] = int(seed)
    kwargs["save_path"] = str(data_path)
    settings = CostAwareSettings(**kwargs)
    rng = np.random.default_rng(settings.rng_seed)
    backend = settings.backend_factory(settings, rng)
    if project_name is not None:
        backend.project_name = project_name
    if device_name is not None:
        backend.device_name = device_name

    if list_jobs:
        list_recent_execute_jobs(
            project_name=backend.project_name,
            n_recent_jobs=n_recent_jobs,
        )
        return

    device_filter = None if no_device_filter else backend.device_name

    submissions, cost = _reconstruct_plain_quantinuum_submissions(
        backend,
        requests,
        seed=int(seed),
        first_shot_indices=first_shot_indices,
    )
    print(
        f"Recovering iteration {meta['pending_batch'].get('iteration')} from "
        f"{len(submissions)} stitched submission(s), {sum(req.shots for req in requests)} "
        f"requested shot(s), expected cost {cost:.3f} HQC."
    )
    outcomes = _recover_outcomes_from_recent_or_exact_job(
        backend,
        submissions,
        job_id=job_id,
        n_recent_jobs=n_recent_jobs,
        device_name=device_filter,
        allow_circuit_mismatch=allow_circuit_mismatch,
    )
    recovered_shots = sum(len(values) for values in outcomes.values())
    print(f"Recovered {recovered_shots} shot outcome(s).")
    if dry_run:
        return

    rmb = RMB.load(data_path, rng=rng)
    for config, values in outcomes.items():
        estimator = rmb._data.setdefault(config, rmb.backend.default_estimator())
        for value in values:
            estimator.record(bool(value))
    rmb.save(data_path)

    batch_history = list(meta.get("batch_history", []))
    jobs = int(meta.get("jobs", 0)) + 1
    batch_history.append(
        _batch_checkpoint_record(
            batch_number=jobs,
            iteration=int(meta["pending_batch"].get("iteration", jobs)),
            cost_hqc=float(cost),
            requests=requests,
            outcomes=outcomes,
            first_shot_indices=first_shot_indices,
        )
    )

    spent = float(meta.get("spent_hqc", 0.0)) + float(cost)
    hqc_budget = float(meta.get("hqc_budget", settings.hqc_budget))
    meta["reason"] = "recovered_after_backend"
    meta["spent_hqc"] = spent
    meta["remaining_hqc"] = max(hqc_budget - spent, 0.0)
    meta["jobs"] = jobs
    meta["max_job_circuits"] = max(
        int(meta.get("max_job_circuits", 0)),
        int(sum(req.shots for req in requests)),
    )
    meta["submitted_configs"] = sum(
        len(batch.get("requests", [])) for batch in batch_history
    )
    meta["measured_configs"] = _measured_config_count(rmb._data)
    meta["batch_history"] = batch_history
    meta.pop("pending_batch", None)
    _write_json_atomic(checkpoint_path, meta)
    print(f"Updated {data_path}")
    print(f"Updated {checkpoint_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--checkpoint",
        type=Path,
        default=(
            GENERATED_PERSONAL_ROOT
            / "CostAware"
            / "seed_2026"
            / "CostAware_restartable_checkpoint.json"
        ),
        help="COST_AWARE checkpoint metadata with reason='before_backend'.",
    )
    parser.add_argument(
        "--data-path",
        type=Path,
        default=None,
        help="RMB data JSON. Defaults to checkpoint path without '_checkpoint'.",
    )
    parser.add_argument(
        "--job-id",
        default=None,
        help="Recover this exact QNexus execute job instead of newest completed jobs.",
    )
    parser.add_argument(
        "--n-recent-jobs",
        type=int,
        default=1,
        help="Number of most recent completed execute jobs to inspect.",
    )
    parser.add_argument(
        "--project-name",
        default=None,
        help="Override the project name inferred from cost_aware_settings.py.",
    )
    parser.add_argument(
        "--device-name",
        default=None,
        help="Override the device name inferred from the checkpoint backend_model.",
    )
    parser.add_argument(
        "--no-device-filter",
        action="store_true",
        help="Do not filter fetched execute jobs by device/system name.",
    )
    parser.add_argument(
        "--list-jobs",
        action="store_true",
        help="List recent execute jobs in the project and exit.",
    )
    parser.add_argument(
        "--allow-circuit-mismatch",
        action="store_true",
        help="Recover even if fetched stitched circuit dimensions differ.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Fetch and validate results without writing local files.",
    )
    args = parser.parse_args()
    recover_pending_batch(
        args.checkpoint,
        args.data_path,
        job_id=args.job_id,
        n_recent_jobs=args.n_recent_jobs,
        project_name=args.project_name,
        device_name=args.device_name,
        no_device_filter=args.no_device_filter,
        list_jobs=args.list_jobs,
        allow_circuit_mismatch=args.allow_circuit_mismatch,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()
