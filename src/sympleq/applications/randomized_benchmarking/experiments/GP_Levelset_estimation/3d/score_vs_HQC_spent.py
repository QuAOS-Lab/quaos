"""Plot 3D success-side volume against recomputed HQC spent/Number of stitched submissions."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from score3d import score_grid
from sympleq.integrations.quantinuum.utils import BASE_SIMULATION_COST


FOLDER = Path(r"Personal\seed_2025")

LABELS = ['Emulator','Sympleq_Seed-2025', 'Sympleq_Seed-2026',
          'Sympleq_Seed-2027','Sympleq_Seed-2028', 'Sympleq_Seed-2029']
SAVE_FIG_PATH = Path(r"Personal\RMB_results_figs\Volume\score_vs_nss_20q.pdf")

# LABELS = ['Emulator (Upto 20q)']
# SAVE_FIG_PATH = Path(r"Personal\RMB_results_figs\Volume\score_vs_nss_emulator20.png")

def config_key(record: dict) -> tuple[int, int, int, bool]:
    return (
        int(record["n_qubits"]),
        int(record["n_1qb_gates"]),
        int(record["n_2qb_gates"]),
        bool(record["random_elimination"]),
    )


def counts_by_config(json_path: Path) -> dict[tuple[int, int, int, bool], int]:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    return {
        config_key(record): sum(int(count) for _, count in record["results"])
        for record in payload.get("data", [])
    }


def approximate_batch_hqc(requests: list[tuple[tuple[int, int, int, bool], int]]) -> float:
    if not requests:
        return 0.0
    bare = 0.0
    for (n_qubits, n_1q, n_2q, _), shots in requests:
        bare += shots * (0.0001 * n_1q + 0.001 * n_2q + n_qubits / 5000)
    return float(BASE_SIMULATION_COST + bare)


def sibling_grid(json_path: Path) -> Path:
    return json_path.parent / f"{json_path.stem}_gp_grid_3d.npz"


def plot_folder(folder: str | Path, use_x: str| None = None, label: str | None = None) -> bool:
    spent = 0.0
    previous: dict[tuple[int, int, int, bool], int] = {}
    xs: list[float] = []
    ys: list[float] = []
    sobol_done_x: float | None = None

    print(f"\n[{folder}]", flush=True)
    for json_path in sorted(Path(folder).rglob("*.json")):
        grid_path = sibling_grid(json_path)
        if not grid_path.exists():
            continue

        print(f"processing {json_path.name} ...", flush=True)
        current = counts_by_config(json_path)
        requests = [
            (config, count - previous.get(config, 0))
            for config, count in current.items()
            if count > previous.get(config, 0)
        ]
        spent += approximate_batch_hqc(requests)
        previous = current

        xs.append(spent)

        ys.append(score_grid(grid_path)["success_side_volume"])
        if "sobol" in json_path.stem.lower():
            sobol_done_x = xs[-1]
            sobol_done_mes = len(xs)
        print(f"{json_path.name}: x={xs[-1]:.6g}, y={ys[-1]:.6g}", flush=True)

    if not xs:
        print(f"No matching JSON + *_gp_grid_3d.npz pairs found in {folder}", flush=True)
        return False
    if use_x == "hqc":
        plt.plot(xs, ys, "o-", label=label or Path(folder).name)
        plt.xlabel("HQC spent")
        plt.ylabel("Success-side volume")
        if sobol_done_x is not None:
            plt.axvline(sobol_done_x, color="black", linestyle="--", linewidth=1.0, label='Initial Sobol')
    else:
        num_mes = [i+1 for i in range(len(ys))]
        if label is not None:
            plt.plot(num_mes, ys, "o-", label=label)
        else:
            plt.plot(num_mes, ys, "o-", label=Path(folder).name)
        plt.xlabel("Number of stitched submissions")
        plt.ylabel("Success-side volume")

        if sobol_done_x is not None:
            plt.axvline(sobol_done_mes, color="black",
                        linestyle="--", linewidth=1.0, label='Initial Sobol')
    return True


def main(folders: list[str | Path] | None = None) -> None:
    folders = [FOLDER] if not folders else folders
    plotted = False

    for index, folder in enumerate(folders):
        label = LABELS[index] if index < len(LABELS) else None
        plotted = plot_folder(folder, label=label) or plotted

    if not plotted:
        return

    plt.legend()
    plt.tight_layout()
    plt.savefig(SAVE_FIG_PATH, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":

    main([Path(arg) for arg in sys.argv[1:]])
