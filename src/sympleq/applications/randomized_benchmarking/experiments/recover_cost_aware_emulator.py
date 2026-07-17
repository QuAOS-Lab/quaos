"""Recover completed Quantinuum emulator results from a failed cost-aware batch.

Run this before rerunning ``common_run_models.py``.  It can recover only runs
whose checkpoint contains a ``pending_batch`` entry, which is written by current
cost-aware code just before calling the backend.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from numpy.random import default_rng
from pytket.qasm.qasm import circuit_to_qasm_str

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.base import (
    MeasurementRequest,
)
from sympleq.applications.randomized_benchmarking.backends.quantinuum import (
    QuantinuumBackend,
)
from sympleq.applications.randomized_benchmarking.experiments.common import (
    BASE_SIMULATION_COST,
    Budget,
    measurement_rng,
    pytket_bare_simulation_cost,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_settings import (
    CostAwareSettings,
    RNG_SEEDS,
    control_panel_settings_kwargs,
)
from sympleq.applications.randomized_benchmarking.experiments.cost_aware_surface_method import (
    _batch_checkpoint_record,
    _checkpoint_meta_path,
    _config_from_checkpoint_record,
    _save_measurement_checkpoint,
    _submitted_from_batch_history,
)
from sympleq.integrations.quantinuum.stitching import (
    MAX_QASM_PROGRAM_SIZE,
    circuit_stitching,
    destitch_results,
    estimate_qasm_program_size,
)
from sympleq.integrations.quantinuum.utils import fetch_recent_execute_jobs


RECENT_EXECUTE_JOBS = 30


def _restartable_cost_aware_save_path(seed: int | None = None) -> Path:
    seed_folder = "seed_unseeded" if seed is None else f"seed_{seed}"
    return Path("Personal") / "CostAware" / seed_folder / "CostAware_restartable.json"


def _settings_for_seed(seed: int | None) -> CostAwareSettings:
    kwargs = control_panel_settings_kwargs()
    kwargs["rng_seed"] = seed
    if kwargs.get("save_path") is None:
        kwargs["save_path"] = _restartable_cost_aware_save_path(seed=seed)
    return CostAwareSettings(**kwargs)


def _requests_from_pending(meta: dict) -> tuple[list[MeasurementRequest], dict]:
    pending = meta.get("pending_batch")
    if not isinstance(pending, dict):
        return [], {}
    requests = []
    first_shot_indices = {}
    for row in pending.get("requests", []):
        config = _config_from_checkpoint_record(row["config"])
        requests.append(MeasurementRequest(config, int(row["requested_shots"])))
        first_shot_indices[config] = int(row.get("first_shot_index", 0))
    return requests, first_shot_indices


def _qasm_key(circuit) -> str:
    return circuit_to_qasm_str(circuit, header="hqslib1")


def _submission_cost(submission) -> float:
    cost = float(BASE_SIMULATION_COST)
    for i, (_, circuit) in enumerate(submission):
        cost += pytket_bare_simulation_cost(circuit)
        if i:
            cost += circuit.n_qubits / 5000
    return cost


def _logical_cost(jobs) -> float:
    cost = float(BASE_SIMULATION_COST)
    for i, (config, circuit) in enumerate(jobs):
        cost += pytket_bare_simulation_cost(circuit)
        if i:
            cost += config.n_qubits / 5000
    return cost


def _prepare_submissions(backend, requests, rng, shot_rng):
    jobs = []
    shot_counts = {}
    for request in requests:
        for _ in range(max(0, request.shots)):
            index = shot_counts.get(request.config, 0)
            shot_counts[request.config] = index + 1
            circuit_rng = rng if shot_rng is None else shot_rng(request.config, index)
            jobs.append((request.config, backend.compatible_circuit(request.config, circuit_rng)))

    def stitched_gate_count(submission) -> int:
        return int(circuit_stitching([circuit for _, circuit in submission]).n_gates)

    def submission_size(submission) -> int:
        return int(sum(estimate_qasm_program_size(circuit) for _, circuit in submission))

    def fits(submission) -> bool:
        if not submission:
            return True
        max_cost = getattr(backend, "max_emulator_batch_cost", None)
        if max_cost is None:
            max_cost = getattr(backend, "max_cost_per_run", None)
        max_gates = getattr(backend, "max_stitched_gates", None)
        if submission_size(submission) > MAX_QASM_PROGRAM_SIZE:
            return False
        if max_cost is not None and _submission_cost(submission) > float(max_cost):
            return False
        return max_gates is None or stitched_gate_count(submission) <= int(max_gates)

    def balanced_submissions():
        if not jobs:
            return []
        if fits(jobs):
            return [jobs]
        max_gates = getattr(backend, "max_stitched_gates", None)
        if max_gates is None:
            chunks = []
            current = []
            for job in jobs:
                candidate = current + [job]
                if current and not fits(candidate):
                    chunks.append(current)
                    current = []
                current.append(job)
            if current:
                chunks.append(current)
            return chunks

        total_gates = stitched_gate_count(jobs)
        parts = max(1, -(-total_gates // int(max_gates)))
        while parts <= len(jobs):
            target = total_gates / parts
            chunks = []
            index = 0
            ok = True
            for remaining_parts in range(parts, 0, -1):
                current = []
                while index < len(jobs):
                    candidate = current + [jobs[index]]
                    if not fits(candidate):
                        if not current:
                            ok = False
                        break

                    must_leave_one_per_remaining = (
                        len(jobs) - (index + 1) < remaining_parts - 1
                    )
                    if current and not must_leave_one_per_remaining:
                        current_gates = stitched_gate_count(current)
                        candidate_gates = stitched_gate_count(candidate)
                        if abs(candidate_gates - target) > abs(current_gates - target):
                            break

                    current = candidate
                    index += 1
                    if len(jobs) - index == remaining_parts - 1:
                        break

                if not current:
                    ok = False
                    break
                chunks.append(current)

            if ok and index == len(jobs):
                return chunks
            parts += 1

        return [[job] for job in jobs]

    submissions = balanced_submissions()

    stitched_circuits = [
        circuit_stitching([circuit for _, circuit in submission])
        for submission in submissions
    ]
    return submissions, stitched_circuits, _logical_cost(jobs)


def recover_seed(seed: int | None) -> None:
    settings = _settings_for_seed(seed)
    if settings.save_path is None:
        raise ValueError("Cost-aware save_path is None; nothing to recover.")

    meta_path = _checkpoint_meta_path(settings.save_path)
    if not meta_path.exists():
        print(f"No checkpoint metadata found at {meta_path}")
        return
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    requests, first_shot_indices = _requests_from_pending(meta)
    if not requests:
        print(
            "No pending_batch found. Existing completed online jobs cannot be "
            "mapped back to configs without that metadata."
        )
        return

    rmb = RMB.load(settings.save_path)
    backend = settings.backend_factory(settings, rmb.rng)
    if not isinstance(backend, QuantinuumBackend):
        raise TypeError(f"Recovery expects QuantinuumBackend, got {type(backend).__name__}")

    def shot_rng(config, index):
        return measurement_rng(settings.rng_seed, config, first_shot_indices[config] + index)

    submissions, stitched_circuits, logical_cost = _prepare_submissions(
        backend, requests, default_rng(settings.rng_seed), shot_rng
    )
    pending_by_qasm = {
        _qasm_key(circuit): (submission, circuit, _submission_cost(submission))
        for submission, circuit in zip(submissions, stitched_circuits)
    }

    outcomes = {}
    matched = 0
    recovered_cost = 0.0
    for circuit, result in fetch_recent_execute_jobs(
        backend.project_name,
        RECENT_EXECUTE_JOBS,
        device_name=backend.device_name,
    ):
        match = pending_by_qasm.get(_qasm_key(circuit))
        if match is None:
            continue
        submission, stitched, cost = match
        registers = sorted(
            stitched.c_registers,
            key=lambda register: int(register.name.removeprefix("creg_")),
        )
        for (config, _), sub_result in zip(submission, destitch_results(result, registers)):
            counts = sub_result.get_empirical_distribution().as_counter()
            if counts:
                outcome, _ = counts.most_common()[0]
                outcomes.setdefault(config, []).append(all(bit == 0 for bit in outcome))
        matched += 1

    if not outcomes:
        print("No completed pending submissions were found in recent Qnexus jobs.")
        return
    recovered_cost = logical_cost

    for config, values in outcomes.items():
        estimator = rmb._data.setdefault(config, backend.default_estimator())
        for outcome in values:
            estimator.record(bool(outcome))

    recovered_requests = [
        MeasurementRequest(config, sum(counter.values()))
        for config, values in outcomes.items()
        for counter in [Counter(values)]
    ]
    batch_history = list(meta.get("batch_history", []))
    iteration = int(meta.get("pending_batch", {}).get("iteration") or meta.get("iteration") or 0)
    batch_history.append(
        _batch_checkpoint_record(
            batch_number=int(meta.get("jobs", 0)) + 1,
            iteration=iteration,
            cost_hqc=recovered_cost,
            requests=recovered_requests,
            outcomes=outcomes,
            first_shot_indices=first_shot_indices,
        )
        | {"source": "recovered_quantinuum_jobs", "matched_submissions": matched}
    )

    budget = Budget(
        remaining_hqc=max(
            float(meta.get("remaining_hqc", settings.hqc_budget)) - recovered_cost,
            0.0,
        ),
        spent_hqc=float(meta.get("spent_hqc", 0.0)) + recovered_cost,
        jobs=int(meta.get("jobs", 0)) + 1,
        max_job_circuits=max(
            int(meta.get("max_job_circuits", 0)),
            sum(request.shots for request in recovered_requests),
        ),
    )
    submitted = _submitted_from_batch_history(batch_history)
    _save_measurement_checkpoint(
        rmb,
        settings,
        budget,
        submitted,
        batch_history,
        iteration=iteration,
        reason="recovered_quantinuum_jobs",
    )
    print(
        f"Recovered {sum(len(v) for v in outcomes.values())} outcomes from "
        f"{matched} completed stitched submissions into {Path(settings.save_path)}"
    )


if __name__ == "__main__":
    for seed in RNG_SEEDS:
        recover_seed(seed)
