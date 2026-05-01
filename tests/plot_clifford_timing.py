import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
CSV_FILENAME = SCRIPT_DIR / "existing_clifford_timing_csv9.csv"
PLOT_FILENAME = SCRIPT_DIR / "clifford_timing_csv9.png"
MAX_QUBITS = 1000
# ONLY_MODELS = {"fermi_hubbard"}
ONLY_MODELS = None
# ONLY_METHODS = {"existing", "igraph"}
ONLY_METHODS = None


def load_results_csv(filename=CSV_FILENAME):
    def parse_int(value, default=None):
        return int(value) if value not in ("", "None", None) else default

    filename = Path(filename)
    if not filename.exists():
        raise FileNotFoundError(f"No timing CSV found: {filename}")

    results = []
    with open(filename, newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            row["nx"] = parse_int(row["nx"])
            row["ny"] = parse_int(row["ny"])
            row["qubits"] = parse_int(row["qubits"])
            row["paulis"] = parse_int(row["paulis"])
            row["found"] = row["found"] == "True"
            row["count"] = parse_int(row["count"], 0)
            row["checked"] = parse_int(row["checked"], 0)
            row["timeout"] = row["timeout"] == "True"
            row["seconds"] = float(row["seconds"])
            row["wall_seconds"] = float(row["wall_seconds"])
            results.append(row)
    return results


def plot_results(results, filename=PLOT_FILENAME):
    if not results:
        return
    results = [r for r in results if r["qubits"] is not None and r["qubits"] <= MAX_QUBITS]
    if ONLY_MODELS is not None:
        results = [r for r in results if r["model"] in ONLY_MODELS]
    if ONLY_METHODS is not None:
        results = [r for r in results if r["method"] in ONLY_METHODS]
    if not results:
        return

    model_colors = {
        "toric": "tab:blue",
        "ising_ladder": "tab:orange",
        "tv_chain": "tab:brown",
        "heisenberg_chain": "tab:green",
        "fermi_hubbard": "tab:red",
        "random_swap_symmetric": "tab:purple",
    }
    method_styles = {
        "existing": "-",
        "igraph": "--",
    }

    method_markers = {
        "existing": "o",
        "igraph": "s",
    }

    fig, ax = plt.subplots(figsize=(8, 5))
    for model_name in sorted({r["model"] for r in results}):
        for method in sorted({r["method"] for r in results if r["model"] == model_name}):
            method_results = sorted(
                [
                    r for r in results
                    if r["model"] == model_name and r["method"] == method and r["qubits"] is not None
                ],
                key=lambda r: r["qubits"],
            )
            if not method_results:
                continue

            qubits = [r["qubits"] for r in method_results]
            seconds = [r["seconds"] for r in method_results]
            ax.plot(
                qubits,
                seconds,
                marker=method_markers.get(method),
                linestyle=method_styles.get(method),
                color=model_colors.get(model_name),
                label=f"{model_name}: {method}",
            )
            timeout_results = [r for r in method_results if r["timeout"]]
            if timeout_results:
                ax.scatter(
                    [r["qubits"] for r in timeout_results],
                    [r["seconds"] for r in timeout_results],
                    marker="x",
                    s=70,
                    color="black",
                    zorder=5,
                    #label=f"{model_name}: {method} timeout",
                )

    ax.set_xlabel("number of qubits")
    ax.set_ylabel("seconds")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_ylim(0.5, 1e5)
    ax.set_xlim(1, 1100)
    ax.set_title("Clifford Symmetry Timing")
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')

    fig.tight_layout()
    fig.savefig(filename, dpi=160)
    plt.close(fig)
    print(f"Saved timing plot: {filename}")


def main():
    plot_results(load_results_csv())


if __name__ == "__main__":
    main()
