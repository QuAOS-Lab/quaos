"""Rewrite a cumulative 2D FLE checkpoint into actual-gate coordinates.

The ordinary 2D FLE checkpoints store cumulative measurement data under the
top-level ``data`` key, but ``experiment.actual_after_elimination`` only
describes the latest submitted batch. This script walks all measurement
checkpoints up to a target file, uses each checkpoint's batch metadata, and
writes a sibling ``*_actual_gates.json`` with cumulative outcomes recorded at
the actual post-elimination gate counts.

It does not submit jobs, recover from QNexus, or modify the source JSONs.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_RUN_FOLDER = Path(
    r"Personal\FLE\H2_2\q56\seed_42\FLE_20260804_175518"
)

DEFAULT_GATES_SET = [
    "S",
    "S_inv",
    "X",
    "Y",
    "Z",
    "V",
    "V_inv",
    "ZZP",
    "ZZP_inv",
]


def debug(message: str) -> None:
    print(f"[fle-2d-actual-gates] {message}")


def outcome_to_int(value) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return int(bool(value))
    text = str(value).strip().lower()
    if text in {"true", "1"}:
        return 1
    if text in {"false", "0"}:
        return 0
    raise ValueError(f"Cannot parse outcome value: {value!r}")


def measurement_step(path: Path, payload: dict | None = None) -> int:
    if payload is not None:
        step = payload.get("experiment", {}).get("step")
        if step is not None:
            return int(step)
    match = re.match(r"measurement_(\d+)_", path.name)
    return int(match.group(1)) if match else 0


def is_source_measurement(path: Path) -> bool:
    return (
        path.name.startswith("measurement_")
        and path.suffix == ".json"
        and "_actual_gates" not in path.stem
    )


def latest_measurement_json(folder: Path) -> Path:
    candidates = sorted(
        (path for path in folder.glob("measurement_*.json") if is_source_measurement(path)),
        key=lambda path: measurement_step(path),
    )
    if not candidates:
        raise FileNotFoundError(f"No source measurement JSONs found in {folder}")
    return candidates[-1]


def resolve_target_json(target: Path | None) -> Path:
    if target is None:
        return latest_measurement_json(DEFAULT_RUN_FOLDER)
    target = Path(target)
    if target.is_dir():
        return latest_measurement_json(target)
    return target


def nominal_key(record: dict) -> tuple[int, int, int, float, bool]:
    return (
        int(record["n_1qb_gates"]),
        int(record["n_2qb_gates"]),
        int(record["n_qubits"]),
        float(record.get("random_elimination", 0.0)),
        bool(record.get("use_scrambler", True)),
    )


def counts_by_nominal_key(payload: dict) -> dict[tuple, Counter]:
    counts: dict[tuple, Counter] = defaultdict(Counter)
    for record in payload.get("data", []):
        key = nominal_key(record)
        for outcome, count in record.get("results", []):
            counts[key][outcome_to_int(outcome)] += int(count)
    return counts


def gates_set_by_nominal_key(payload: dict) -> dict[tuple, list[str]]:
    gates_sets = {}
    for record in payload.get("data", []):
        gates_sets.setdefault(
            nominal_key(record),
            list(record.get("gates_set", DEFAULT_GATES_SET)),
        )
    return gates_sets


def delta_counts(
    current: dict[tuple, Counter],
    previous: dict[tuple, Counter],
) -> dict[tuple, Counter]:
    deltas: dict[tuple, Counter] = defaultdict(Counter)
    for key in set(current) | set(previous):
        for outcome in (0, 1):
            delta = int(current.get(key, Counter()).get(outcome, 0)) - int(
                previous.get(key, Counter()).get(outcome, 0)
            )
            if delta < 0:
                raise ValueError(f"Negative count delta for {key}, outcome={outcome}: {delta}")
            if delta:
                deltas[key][outcome] = delta
    return deltas


def pop_observed_outcome(counter: Counter, *, key: tuple, step: int) -> int:
    for outcome in (0, 1):
        if counter[outcome] > 0:
            counter[outcome] -= 1
            return outcome
    raise ValueError(f"No remaining observed outcome for key={key} at step={step}")


def actual_record_template(
    *,
    sent: dict,
    actual: dict,
    gates_set: list[str],
) -> dict:
    if not bool(actual.get("available", True)):
        debug("warning: actual metadata unavailable; falling back to nominal sent config")
        actual = sent

    return {
        "n_1qb_gates": int(actual["n_1qb_gates"]),
        "n_2qb_gates": int(actual["n_2qb_gates"]),
        "n_qubits": int(actual["n_qubits"]),
        "random_elimination": float(sent.get("random_elimination", 0.0)),
        "use_scrambler": bool(sent.get("use_scrambler", True)),
        "gates_set": list(gates_set),
    }


def actual_key(record: dict) -> tuple[int, int, int, float, bool, tuple[str, ...]]:
    return (
        int(record["n_1qb_gates"]),
        int(record["n_2qb_gates"]),
        int(record["n_qubits"]),
        float(record.get("random_elimination", 0.0)),
        bool(record.get("use_scrambler", True)),
        tuple(record.get("gates_set", DEFAULT_GATES_SET)),
    )


def checkpoint_paths(target_json: Path) -> list[Path]:
    target_payload = json.loads(target_json.read_text(encoding="utf-8"))
    target_step = measurement_step(target_json, target_payload)
    paths = []
    for path in target_json.parent.glob("measurement_*.json"):
        if not is_source_measurement(path):
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        step = measurement_step(path, payload)
        if 0 < step <= target_step:
            paths.append(path)
    return sorted(paths, key=lambda p: measurement_step(p))


def rewrite_to_actual_gates(
    target_json: Path,
    *,
    output_json: Path | None = None,
) -> Path:
    target_json = Path(target_json)
    output_json = output_json or target_json.with_name(f"{target_json.stem}_actual_gates.json")

    debug(f"target JSON: {target_json}")
    paths = checkpoint_paths(target_json)
    if not paths:
        raise FileNotFoundError(f"No measurement checkpoints found for {target_json}")
    debug(f"checkpoints used: {len(paths)}")
    debug(f"first checkpoint: {paths[0].name}")
    debug(f"last checkpoint:  {paths[-1].name}")

    previous_counts: dict[tuple, Counter] = defaultdict(Counter)
    actual_counts: dict[tuple, Counter] = defaultdict(Counter)
    actual_templates: dict[tuple, dict] = {}
    ambiguous_duplicate_batches = 0

    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        step = measurement_step(path, payload)
        experiment = payload.get("experiment", {})
        sent_configs = experiment.get("sent_configs") or []
        actual_configs = experiment.get("actual_after_elimination") or []
        if len(sent_configs) != len(actual_configs):
            raise ValueError(
                f"{path.name}: sent_configs ({len(sent_configs)}) and "
                f"actual_after_elimination ({len(actual_configs)}) differ."
            )

        current_counts = counts_by_nominal_key(payload)
        deltas = delta_counts(current_counts, previous_counts)
        gates_sets = gates_set_by_nominal_key(payload)

        batch_keys = [nominal_key(sent) for sent in sent_configs]
        duplicate_keys = [key for key, count in Counter(batch_keys).items() if count > 1]
        for key in duplicate_keys:
            actual_key_count = len({
                actual_key(
                    actual_record_template(
                        sent=sent,
                        actual=actual,
                        gates_set=gates_sets.get(nominal_key(sent), DEFAULT_GATES_SET),
                    )
                )
                for sent, actual in zip(sent_configs, actual_configs)
                if nominal_key(sent) == key
            })
            if actual_key_count > 1 and deltas[key][0] and deltas[key][1]:
                ambiguous_duplicate_batches += 1
                debug(
                    "warning: duplicate nominal config with mixed outcomes at "
                    f"step={step}; assigning outcomes in stored count order."
                )

        debug(
            f"step={step:03d} phase={experiment.get('phase')} "
            f"data_records={len(payload.get('data', []))} "
            f"batch_configs={len(sent_configs)} "
            f"new_observations={sum(sum(counter.values()) for counter in deltas.values())}"
        )

        for sent, actual in zip(sent_configs, actual_configs):
            key = nominal_key(sent)
            outcome = pop_observed_outcome(deltas[key], key=key, step=step)
            template = actual_record_template(
                sent=sent,
                actual=actual,
                gates_set=gates_sets.get(key, DEFAULT_GATES_SET),
            )
            key_actual = actual_key(template)
            actual_templates.setdefault(key_actual, template)
            actual_counts[key_actual][outcome] += 1

        leftover = sum(sum(counter.values()) for counter in deltas.values())
        if leftover:
            raise ValueError(f"{path.name}: {leftover} unassigned new observations remain.")

        previous_counts = current_counts

    target_payload = json.loads(target_json.read_text(encoding="utf-8"))
    output_records = []
    for key, template in actual_templates.items():
        counts = actual_counts[key]
        record = dict(template)
        record["results"] = [
            [bool(outcome), int(count)]
            for outcome, count in ((0, counts[0]), (1, counts[1]))
            if count
        ]
        output_records.append(record)

    output_payload = copy.deepcopy(target_payload)
    output_payload["data"] = output_records
    experiment = output_payload.setdefault("experiment", {})
    experiment["data_coordinate_system"] = "actual_after_elimination"
    experiment["actual_gates_rewrite"] = {
        "source_json": str(target_json),
        "output_json": str(output_json),
        "checkpoints_used": [path.name for path in paths],
        "input_data_records": len(target_payload.get("data", [])),
        "output_data_records": len(output_records),
        "input_total_observations": sum(
            sum(counter.values())
            for counter in counts_by_nominal_key(target_payload).values()
        ),
        "output_total_observations": sum(sum(counter.values()) for counter in actual_counts.values()),
        "ambiguous_duplicate_batches": int(ambiguous_duplicate_batches),
    }

    output_json.write_text(json.dumps(output_payload, indent=2), encoding="utf-8")
    debug(f"wrote actual-gates JSON: {output_json}")
    debug(f"actual-gates records: {len(output_records)}")
    return output_json


def config_from_record(record: dict):
    from sympleq.applications.randomized_benchmarking.config import RMBConfig

    return RMBConfig(
        n_1qb_gates=int(record["n_1qb_gates"]),
        n_2qb_gates=int(record["n_2qb_gates"]),
        n_qubits=int(record["n_qubits"]),
        random_elimination=float(record.get("random_elimination", 0.1)),
        use_scrambler=bool(record.get("use_scrambler", True)),
    )


def make_gp_grid(json_path: Path) -> Path:
    import torch

    from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.fantasy_levelset_estimation import (
        Observation,
        add_observation_to_strategy,
        build_strategy,
        choose_gp_device,
        point_from_config,
        refreshed_strategy_for_prediction,
        seed_fake_corners,
    )
    from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.fantasy_levelset_settings import (
        FantasySettings,
        control_panel_settings_kwargs,
    )
    from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.run_FLE import (
        save_gp_prediction_grid,
    )

    debug(f"building GP grid from actual-gates JSON: {json_path}")
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    settings = FantasySettings(**control_panel_settings_kwargs())
    device = choose_gp_device(settings)
    strategy = build_strategy(settings)
    observations: list[Observation] = []
    results_for_plot: list[tuple[float, float, int]] = []

    debug("seeding fake corners")
    seed_fake_corners(
        strategy,
        settings,
        observations,
        results_for_plot,
        device=device,
    )

    real_count = 0
    for record_index, record in enumerate(payload.get("data", []), start=1):
        config = config_from_record(record)
        for outcome, count in record.get("results", []):
            for _ in range(int(count)):
                obs = Observation(
                    x_cpu=point_from_config(config, device=torch.device("cpu")),
                    y=outcome_to_int(outcome),
                    source="actual_gates_json",
                )
                add_observation_to_strategy(strategy, obs, device=device)
                observations.append(obs)
                real_count += 1
        if record_index % 50 == 0:
            debug(f"replayed records: {record_index}/{len(payload.get('data', []))}")

    debug(f"real observations replayed: {real_count}")
    debug("refreshing strategy for prediction")
    prediction_strategy = refreshed_strategy_for_prediction(
        strategy,
        settings,
        observations,
        device=device,
    )
    grid_path = save_gp_prediction_grid(
        prediction_strategy,
        settings,
        device=device,
        json_path=json_path,
    )
    debug(f"wrote GP grid: {grid_path}")
    return grid_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "target_json",
        nargs="?",
        type=Path,
        default=None,
        help=(
            "Target cumulative measurement JSON, or a run folder. Default is the latest "
            "source measurement JSON in the current H2-2 q56 run folder."
        ),
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Optional output JSON path. Default is sibling *_actual_gates.json.",
    )
    parser.add_argument(
        "--make-grid",
        action="store_true",
        help="Also replay the actual-gates JSON into AEPsych and write *_gp_grid.npz.",
    )
    parser.add_argument(
        "--grid-only",
        action="store_true",
        help=(
            "Only replay target_json into AEPsych and write *_gp_grid.npz. "
            "Use this with an existing *_actual_gates.json."
        ),
    )
    args = parser.parse_args()

    target_json = resolve_target_json(args.target_json)
    if args.grid_only:
        make_gp_grid(target_json)
        return

    output_json = rewrite_to_actual_gates(
        target_json,
        output_json=args.output_json,
    )
    if args.make_grid:
        make_gp_grid(output_json)


if __name__ == "__main__":
    main()
