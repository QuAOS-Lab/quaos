from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import qnexus as qnx
from qnexus.client.auth import is_logged_in
from qnexus.models.filters import SortFilterEnum
from qnexus.models.job_status import JobStatusEnum
from qnexus.models.references import ExecutionResultRef, JobType


OUTPUT_ROOT = Path(r"Personal\FLE\Fancy_emulator\hardware_q26_nominal")


def ensure_qnexus_login() -> None:
    if not is_logged_in():
        qnx.login()


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
        print(f"[lookup] direct qnx.jobs.get failed: {exc}", flush=True)

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


def latest_manifest() -> Path:
    candidates = sorted(OUTPUT_ROOT.glob("nominal_q26_*/nominal_q26_job_manifest.json"))
    if not candidates:
        raise FileNotFoundError(f"No job manifest under {OUTPUT_ROOT}")
    return candidates[-1]


def summary(records: list[dict]) -> dict:
    by_target = defaultdict(int)
    by_source = defaultdict(int)
    agreement = 0
    compared = 0
    for record in records:
        by_target[record["target_device"]] += 1
        by_source[record["source_device"]] += 1
        if "emulator_success" in record:
            compared += 1
            if bool(record["hardware_success"]) == bool(record["emulator_success"]):
                agreement += 1
    return {
        "total_circuits": len(records),
        "completed_circuits": compared,
        "hardware_emulator_agreement": agreement,
        "by_source_device": dict(by_source),
        "by_target_device": dict(by_target),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Recover nominal q26 Fancy-emulator results from saved QNexus job ids."
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="Path to nominal_q26_job_manifest.json. Defaults to latest.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest_path = args.manifest or latest_manifest()
    print(f"[manifest] {manifest_path}", flush=True)

    with manifest_path.open("r", encoding="utf-8") as handle:
        manifest = json.load(handle)

    ensure_qnexus_login()

    circuits = {
        int(record["index"]): dict(record)
        for record in manifest["circuits"]
    }

    for batch in manifest["batches"]:
        project_name = batch["project_name"]
        job_id = batch["execute_job_id"]
        circuit_indices = [int(index) for index in batch["circuit_indices"]]
        print(
            "[fetch] "
            f"target={batch['target_device']} batch={batch['batch_index']:03d} "
            f"project={project_name} job_id={job_id}",
            flush=True,
        )

        job = find_execute_job(project_name, job_id)
        status = getattr(job, "last_status", None)
        print(f"[status] job_id={job_id} status={status}", flush=True)
        if status != JobStatusEnum.COMPLETED:
            raise RuntimeError(f"Execute job {job_id!r} is not completed.")

        refs = [
            ref
            for ref in qnx.jobs.results(job)
            if isinstance(ref, ExecutionResultRef)
        ]
        if len(refs) != len(circuit_indices):
            raise RuntimeError(
                f"Job {job_id!r}: got {len(refs)} result refs for "
                f"{len(circuit_indices)} circuits."
            )

        for circuit_index, ref in zip(circuit_indices, refs):
            result = ref.download_result()
            record = circuits[circuit_index]
            record["execute_job_id"] = job_id
            record["execute_job_project_name"] = project_name
            record["emulator_success"] = success_from_result(result)
            record["emulator_counts"] = outcome_counts(result)
            print(
                "[result] "
                f"{circuit_index:04d} target={record['target_device']} "
                f"hardware_success={record['hardware_success']} "
                f"emulator_success={record['emulator_success']} "
                f"counts={record['emulator_counts']}",
                flush=True,
            )

    records = [
        circuits[index]
        for index in sorted(circuits)
    ]
    output = dict(manifest)
    output["recovered_at"] = datetime.now().isoformat(timespec="seconds")
    output["summary"] = summary(records)
    output["circuits"] = records

    output_path = (
        manifest_path.parent
        / f"nominal_q26_emulator_results_from_jobs_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    )
    output_path.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(f"[write] {output_path}", flush=True)


if __name__ == "__main__":
    main()
