"""Combine an actual-coordinate recovery base with later mixed checkpoints.

Use this when a recovery run loaded ``*_actual_gates.json`` and then continued
saving ordinary measurement checkpoints. In that case the later checkpoints are
mixed: old records are already in actual-gate coordinates, while the latest
batch records are still nominal. This script keeps the actual-coordinate base,
converts only each later checkpoint delta using that checkpoint's
``actual_after_elimination`` metadata, and writes one combined actual-gates JSON.
"""

from __future__ import annotations

import argparse
import copy
import json
from collections import Counter, defaultdict
from pathlib import Path

from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.rewrite_fle_2d_actual_gates import (
    DEFAULT_GATES_SET,
    actual_key,
    actual_record_template,
    counts_by_nominal_key,
    delta_counts,
    gates_set_by_nominal_key,
    is_source_measurement,
    measurement_step,
    nominal_key,
    outcome_to_int,
    pop_observed_outcome,
)


DEFAULT_RUN_FOLDER = Path(
    r"Personal\FLE\H_wrapper\H2_1\q56\seed_42\FLE_H_wrapper_20260819_134349"
)
DEFAULT_BASE_ACTUAL_JSON = (
    DEFAULT_RUN_FOLDER
    / "measurement_011_globalsur_20260821_100536_518538_actual_gates.json"
)

# GP grid handles for combined actual-gates JSONs. These are intentionally
# explicit because the run search box can exclude the actual success region.
GRID_N_GATES_BOUNDS = (200.0, 1500.0)
GRID_RATIO_BOUNDS = (0.1, 0.75)


def debug(message: str) -> None:
    print(f"[hwrap-combine-actual] {message}", flush=True)


def latest_source_measurement(folder: Path) -> Path:
    candidates = sorted(
        (path for path in folder.glob("measurement_*.json") if is_source_measurement(path)),
        key=lambda path: measurement_step(path),
    )
    if not candidates:
        raise FileNotFoundError(f"No source measurement JSONs found in {folder}")
    return candidates[-1]


def source_checkpoints(base_actual_json: Path, target_json: Path) -> list[Path]:
    base_step = measurement_step(base_actual_json)
    target_payload = json.loads(target_json.read_text(encoding="utf-8"))
    target_step = measurement_step(target_json, target_payload)
    paths = []
    for path in target_json.parent.glob("measurement_*.json"):
        if not is_source_measurement(path):
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        step = measurement_step(path, payload)
        if base_step < step <= target_step:
            paths.append(path)
    return sorted(paths, key=lambda path: measurement_step(path))


def seed_actual_counts(base_payload: dict) -> tuple[dict[tuple, Counter], dict[tuple, dict]]:
    actual_counts: dict[tuple, Counter] = defaultdict(Counter)
    actual_templates: dict[tuple, dict] = {}
    for record in base_payload.get("data", []):
        template = {
            key: value
            for key, value in record.items()
            if key != "results"
        }
        template.setdefault("gates_set", list(DEFAULT_GATES_SET))
        key = actual_key(template)
        actual_templates.setdefault(key, template)
        for outcome, count in record.get("results", []):
            actual_counts[key][outcome_to_int(outcome)] += int(count)
    return actual_counts, actual_templates


def output_records(
    actual_counts: dict[tuple, Counter],
    actual_templates: dict[tuple, dict],
) -> list[dict]:
    records = []
    for key in sorted(actual_templates):
        counts = actual_counts[key]
        record = dict(actual_templates[key])
        record["results"] = [
            [bool(outcome), int(count)]
            for outcome, count in ((0, counts[0]), (1, counts[1]))
            if count
        ]
        records.append(record)
    return records


def config_from_actual_record(record: dict):
    from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.h_wrapper_config import (
        HWrappedRMBConfig,
    )
    from sympleq.integrations.quantinuum.utils import NATIVE_GATES_SET

    return (
        HWrappedRMBConfig.default()
        .with_n_qubits(int(record["n_qubits"]))
        .with_n_1qb_gates(int(record["n_1qb_gates"]))
        .with_n_2qb_gates(int(record["n_2qb_gates"]))
        .with_random_elimination(float(record.get("random_elimination", 0.1)))
        .with_use_scrambler(bool(record.get("use_scrambler", True)))
        .with_gates_set(tuple(NATIVE_GATES_SET))
    )


def make_hwrap_gp_grid(json_path: Path) -> Path:
    import torch

    from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.fantasy_levelset_estimation import (
        Observation,
        add_observation_to_strategy,
        build_strategy,
        choose_gp_device,
        point_from_config,
        refreshed_strategy_for_prediction,
        seed_fake_corners,
    )
    from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.fantasy_levelset_settings import (
        FantasySettings,
        RNG_SEEDS,
        control_panel_settings_kwargs,
    )
    from sympleq.applications.randomized_benchmarking.experiments.GP_Levelset_estimation.H_wrapper.run_FLE import (
        save_gp_prediction_grid,
    )

    debug(f"building H-wrapper GP grid from actual JSON: {json_path}")
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    settings_kwargs = control_panel_settings_kwargs()
    if RNG_SEEDS:
        settings_kwargs["rng_seed"] = int(RNG_SEEDS[0])
    settings_kwargs["n_gates_bounds"] = GRID_N_GATES_BOUNDS
    settings_kwargs["ratio_bounds"] = GRID_RATIO_BOUNDS
    settings = FantasySettings(**settings_kwargs)
    debug(f"grid n_gates_bounds={settings.n_gates_bounds}")
    debug(f"grid ratio_bounds={settings.ratio_bounds}")
    device = choose_gp_device(settings)
    strategy = build_strategy(settings)
    observations: list[Observation] = []
    results_for_plot: list[tuple[float, float, int]] = []

    seed_fake_corners(
        strategy,
        settings,
        observations,
        results_for_plot,
        device=device,
    )

    real_count = 0
    records = payload.get("data", [])
    for record_index, record in enumerate(records, start=1):
        config = config_from_actual_record(record)
        for outcome, count in record.get("results", []):
            for _ in range(int(count)):
                obs = Observation(
                    x_cpu=point_from_config(config, device=torch.device("cpu")),
                    y=outcome_to_int(outcome),
                    source="combined_actual_gates_json",
                )
                add_observation_to_strategy(strategy, obs, device=device)
                observations.append(obs)
                real_count += 1
        if record_index % 50 == 0:
            debug(f"replayed records: {record_index}/{len(records)}")

    debug(f"real observations replayed: {real_count}")
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
    debug(f"wrote H-wrapper GP grid: {grid_path}")
    return grid_path


def combine_mixed_actual_gates(
    target_json: Path,
    *,
    base_actual_json: Path,
    output_json: Path | None,
    dry_run: bool,
) -> Path:
    target_json = Path(target_json)
    base_actual_json = Path(base_actual_json)
    output_json = output_json or target_json.with_name(
        f"{target_json.stem}_combined_actual_gates.json"
    )

    base_payload = json.loads(base_actual_json.read_text(encoding="utf-8"))
    target_payload = json.loads(target_json.read_text(encoding="utf-8"))
    paths = source_checkpoints(base_actual_json, target_json)
    if not paths:
        raise FileNotFoundError(
            f"No source checkpoints after {base_actual_json.name} up to {target_json.name}"
        )

    actual_counts, actual_templates = seed_actual_counts(base_payload)
    previous_source_counts = counts_by_nominal_key(base_payload)

    debug(f"base actual JSON: {base_actual_json}")
    debug(f"target JSON:      {target_json}")
    debug(f"output JSON:      {output_json}")
    debug(f"later checkpoints used: {len(paths)}")
    debug(f"base actual observations: {sum(sum(c.values()) for c in actual_counts.values())}")

    converted_observations = 0
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

        current_source_counts = counts_by_nominal_key(payload)
        deltas = delta_counts(current_source_counts, previous_source_counts)
        gates_sets = gates_set_by_nominal_key(payload)
        new_observations = sum(sum(counter.values()) for counter in deltas.values())
        debug(
            f"step={step:03d} file={path.name} "
            f"batch_configs={len(sent_configs)} new_observations={new_observations}"
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
            converted_observations += 1

        leftover = sum(sum(counter.values()) for counter in deltas.values())
        if leftover:
            raise ValueError(f"{path.name}: {leftover} unassigned new observations remain.")
        previous_source_counts = current_source_counts

    records = output_records(actual_counts, actual_templates)
    output_payload = copy.deepcopy(target_payload)
    output_payload["data"] = records
    experiment = output_payload.setdefault("experiment", {})
    experiment["data_coordinate_system"] = "actual_after_elimination"
    experiment["actual_gates_rewrite"] = {
        "mode": "mixed_recovery_base_plus_later_deltas",
        "base_actual_json": str(base_actual_json),
        "target_json": str(target_json),
        "output_json": str(output_json),
        "checkpoints_used": [path.name for path in paths],
        "base_actual_records": len(base_payload.get("data", [])),
        "base_actual_observations": sum(
            sum(counter.values())
            for counter in counts_by_nominal_key(base_payload).values()
        ),
        "converted_later_observations": int(converted_observations),
        "output_data_records": len(records),
        "output_total_observations": sum(sum(counter.values()) for counter in actual_counts.values()),
    }

    debug(f"converted later observations: {converted_observations}")
    debug(f"output records: {len(records)}")
    debug(
        "output observations: "
        f"{sum(sum(counter.values()) for counter in actual_counts.values())}"
    )

    if dry_run:
        debug("dry run: not writing output")
        return output_json

    output_json.write_text(json.dumps(output_payload, indent=2), encoding="utf-8")
    debug(f"wrote combined actual-gates JSON: {output_json}")
    return output_json


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "target_json",
        nargs="?",
        type=Path,
        default=None,
        help="Latest mixed measurement JSON. Default: latest source measurement in H-wrapper run folder.",
    )
    parser.add_argument(
        "--base-actual-json",
        type=Path,
        default=DEFAULT_BASE_ACTUAL_JSON,
        help="Actual-gates JSON that recovery loaded as the base.",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=None,
        help="Output path. Default: sibling *_combined_actual_gates.json.",
    )
    parser.add_argument(
        "--make-grid",
        action="store_true",
        help="Also write the H-wrapper *_gp_grid.npz for the combined actual JSON.",
    )
    parser.add_argument(
        "--grid-only",
        action="store_true",
        help="Only write the H-wrapper *_gp_grid.npz for target_json.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    target_json = args.target_json or latest_source_measurement(DEFAULT_RUN_FOLDER)
    if args.grid_only:
        make_hwrap_gp_grid(target_json)
        return

    output_json = combine_mixed_actual_gates(
        target_json,
        base_actual_json=args.base_actual_json,
        output_json=args.output_json,
        dry_run=args.dry_run,
    )
    if args.make_grid and not args.dry_run:
        make_hwrap_gp_grid(output_json)


if __name__ == "__main__":
    main()
