import multiprocessing as mp
import csv
from pathlib import Path
from queue import Empty
from time import perf_counter

import igraph as ig
import numpy as np

from sympleq.core.finite_field_solvers import get_linear_dependencies, _select_row_basis_indices
from sympleq.core.graphs.graph_automorphism import _LeafContext, _check_leaf_with_reason
from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.core.symmetries.clifford import find_clifford_symmetries
from sympleq.core.circuits import GATES
from sympleq.models.fermi_hubbard import disordered_tv_chain_model
from sympleq.models.heisenberg import all_to_all_heisenberg_hamiltonian
from sympleq.models.Ising import ising_2d_hamiltonian
# from sympleq.models.random_hamiltonian import random_gate_symmetric_hamiltonian
from sympleq.models.toric_code import ToricCode

TIMEOUT_SECONDS = 10**6
SCRIPT_DIR = Path(__file__).resolve().parent
CSV_FILENAME = SCRIPT_DIR / "clifford_timing_csv10.csv"
MAX_QUBITS = 1000


def build_model_pauli_sum(model_name, nx, ny, periodic):
    if model_name == "toric":
        return ToricCode(Nx=nx, Ny=ny, c_x=1.0, c_z=2.0, c_g=0.1, periodic=periodic).hamiltonian()
    if model_name == "ising_ladder":
        return ising_2d_hamiltonian(n_x=nx, n_y=ny, J_zz=1.0, h_x=0.5, periodic=periodic)
    if model_name == "tv_chain":

        return disordered_tv_chain_model(
            n_sites=ny,
            tunneling=1.0,
            coulomb=4.0,
            disorder_strength=0.5,
            periodic=periodic,
            seed=10_000 + ny,
        )
    if model_name == "heisenberg_chain":
        return all_to_all_heisenberg_hamiltonian(
            n=ny,
            J=1.0,
            delta_z=0.5,
        )
    if model_name == "fermi_hubbard":
        from sympleq.models.fermi_hubbard import fermi_hubbard_model
        return fermi_hubbard_model(
            x_dimension=ny,
            y_dimension=1,
            tunneling=1.0,
            coulomb=4.0,
            chemical_potential=0.0,
            periodic=periodic,
            spinless=False,
        )
    # if model_name == "random_swap_symmetric":
    #     return random_gate_symmetric_hamiltonian(
    #         GATES.SWAP,
    #         dimension=2,
    #         qudit_indices=(0, 1),
    #         n_qudits=ny,
    #         n_paulis=2 * ny,
    #         weight_mode="uniform",
    #         scrambled=True,
    #     )
    raise ValueError(f"Unknown model: {model_name}")


def existing_find_worker(model_name, nx, ny, periodic, queue):
    try:
        pauli_sum = build_model_pauli_sum(model_name, nx, ny, periodic)
        start = perf_counter()
        symmetries = find_clifford_symmetries(pauli_sum, num_symmetries=1)
        queue.put({
            "seconds": perf_counter() - start,
            "qubits": pauli_sum.n_qudits(),
            "paulis": pauli_sum.n_paulis(),
            "found": len(symmetries) > 0,
            "count": len(symmetries),
            "error": None,
        })
    except Exception as exc:
        queue.put({
            "seconds": None,
            "qubits": None,
            "paulis": None,
            "found": False,
            "count": 0,
            "error": f"{type(exc).__name__}: {exc}",
        })


def build_leaf_context(pauli_sum):
    pauli = pauli_sum.copy()
    pauli.weight_to_phase()

    independent_labels, dependencies = get_linear_dependencies(pauli.tableau, 2)
    labels = sorted(set(independent_labels) | set(dependencies.keys()))
    S_mod = pauli.symplectic_product_matrix()
    G, basis_order = pauli.matroid()

    pauli_standard = pauli.to_standard_form()
    pauli_standard.weight_to_phase()

    base_tableau = pauli.tableau.astype(int, copy=False)
    dims_array = np.asarray(pauli.dimensions, dtype=int)
    row_basis_cache = {}
    if dims_array.size and np.all(dims_array == dims_array[0]):
        p_uni = int(dims_array[0])
        cache = _select_row_basis_indices(base_tableau % p_uni, p_uni, base_tableau.shape[1])
        row_basis_cache["gf2" if p_uni == 2 else "gfp"] = cache

    return _LeafContext(
        p=int(pauli.lcm),
        two_lcm=2 * pauli_sum.lcm,
        n_qudits=pauli_sum.n_qudits(),
        identity_perm=np.arange(pauli.n_paulis(), dtype=np.int64),
        S_mod=S_mod,
        G=G,
        basis_order=basis_order,
        labels=labels,
        pauli_sum=pauli,
        ref_tableau=pauli_standard.tableau.astype(int, copy=False),
        ref_phases=np.asarray(pauli_standard.phases, dtype=int),
        ref_weights=np.asarray(pauli_standard.weights),
        base_tableau=base_tableau,
        base_weights=np.asarray(pauli.weights),
        base_phases=np.asarray(pauli.phases, dtype=int),
        basis_indices=np.asarray(independent_labels, dtype=int),
        basis_source_ps=pauli[np.asarray(independent_labels, dtype=int)],
        row_basis_cache=row_basis_cache,
    )


def build_subdivision_graph_from_s_mod(S_mod, vertex_colors):
    n = S_mod.shape[0]
    edges = []
    edge_colors = []
    for i in range(n):
        for j in range(i + 1, n):
            edges.append((i, j))
            edge_colors.append(int(S_mod[i, j]))

    unique_vertex_colors = sorted(set(vertex_colors), key=str)
    vertex_color_map = {c: i for i, c in enumerate(unique_vertex_colors)}
    h_colors = [vertex_color_map[int(c)] for c in vertex_colors]

    edge_color_offset = max(h_colors, default=-1) + 1
    unique_edge_colors = sorted(set(edge_colors), key=str)
    edge_color_map = {c: edge_color_offset + i for i, c in enumerate(unique_edge_colors)}

    h = ig.Graph(n=n + len(edges), directed=False)
    subdivided_edges = []
    for edge_idx, (u, v) in enumerate(edges):
        w = n + edge_idx
        subdivided_edges.append((u, w))
        subdivided_edges.append((w, v))
        h_colors.append(edge_color_map[edge_colors[edge_idx]])

    h.add_edges(subdivided_edges)
    h.vs["vertex_color"] = h_colors
    return h, h_colors


def permutation_to_tuple(permutation):
    if hasattr(permutation, "mapping"):
        return tuple(permutation.mapping)
    return tuple(permutation)


def compose(left, right):
    return tuple(left[i] for i in right)


def generated_group_permutations(generators, n_vertices):
    generators = [permutation_to_tuple(g) for g in generators]
    identity = tuple(range(n_vertices))
    seen = {identity}
    queue_items = [(identity, ())]

    while queue_items:
        current, word = queue_items.pop(0)
        for idx, generator in enumerate(generators):
            for side, candidate in (
                ("L", compose(generator, current)),
                ("R", compose(current, generator)),
            ):
                if candidate in seen:
                    continue
                seen.add(candidate)
                candidate_word = word + ((idx, side),)
                queue_items.append((candidate, candidate_word))
                yield candidate, candidate_word


def first_clifford_from_igraph_generators(pauli_sum, h, generators):
    ctx = build_leaf_context(pauli_sum)
    n_paulis = pauli_sum.n_paulis()
    identity_pauli_perm = tuple(range(n_paulis))
    checked = 0

    for idx, generator in enumerate(generators):
        pi = np.asarray(permutation_to_tuple(generator)[:n_paulis], dtype=np.int64)
        if tuple(pi) == identity_pauli_perm:
            continue
        checked += 1
        gate, reason = _check_leaf_with_reason(pi, ctx)
        if gate is not None:
            return True, 1, ((idx, "generator"),), checked, reason

    for full_perm, word in generated_group_permutations(generators, h.vcount()):
        pi_tuple = tuple(full_perm[:n_paulis])
        if pi_tuple == identity_pauli_perm:
            continue
        checked += 1
        gate, reason = _check_leaf_with_reason(np.asarray(pi_tuple, dtype=np.int64), ctx)
        if gate is not None:
            return True, 1, word, checked, reason

    return False, 0, None, checked, "no Clifford lift found"


def igraph_find_worker(model_name, nx, ny, periodic, queue):
    try:
        pauli_sum = build_model_pauli_sum(model_name, nx, ny, periodic)
        start = perf_counter()

        pauli_for_graph = pauli_sum.copy()
        pauli_for_graph.weight_to_phase()
        S_mod = pauli_for_graph.symplectic_product_matrix()
        colors, _ = _build_base_partition(
            S_mod,
            int(pauli_for_graph.lcm),
            coeffs=pauli_for_graph.weights,
            col_invariants=None,
            max_rounds=0,
            color_mode="wl",
        )
        h, h_colors = build_subdivision_graph_from_s_mod(S_mod, colors)
        automorphism_group = h.automorphism_group(color=h_colors)
        generators = automorphism_group.generators if hasattr(automorphism_group, "generators") else automorphism_group
        found, count, word, checked, reason = first_clifford_from_igraph_generators(pauli_sum, h, generators)

        queue.put({
            "seconds": perf_counter() - start,
            "qubits": pauli_sum.n_qudits(),
            "paulis": pauli_sum.n_paulis(),
            "found": found,
            "count": count,
            "checked": checked,
            "word": word,
            "error": None if found else reason,
        })
    except Exception as exc:
        queue.put({
            "seconds": None,
            "qubits": None,
            "paulis": None,
            "found": False,
            "count": 0,
            "checked": 0,
            "word": None,
            "error": f"{type(exc).__name__}: {exc}",
        })


def time_method(worker, method, model_name, nx, ny, periodic, timeout_seconds=TIMEOUT_SECONDS):
    def model_metadata():
        try:
            pauli_sum = build_model_pauli_sum(model_name, nx, ny, periodic)
            return pauli_sum.n_qudits(), pauli_sum.n_paulis(), None
        except Exception as exc:
            return None, None, f"{type(exc).__name__}: {exc}"

    queue = mp.Queue()
    process = mp.Process(
        target=worker,
        args=(model_name, nx, ny, periodic, queue),
    )

    wall_start = perf_counter()
    process.start()
    process.join(timeout_seconds)
    wall_seconds = perf_counter() - wall_start

    if process.is_alive():
        process.terminate()
        process.join()
        qubits, paulis, metadata_error = model_metadata()
        return {
            "method": method,
            "model": model_name,
            "nx": nx,
            "ny": ny,
            "periodic": periodic,
            "seconds": timeout_seconds,
            "wall_seconds": wall_seconds,
            "qubits": qubits,
            "paulis": paulis,
            "found": False,
            "count": 0,
            "checked": 0,
            "word": None,
            "timeout": True,
            "error": metadata_error,
        }

    try:
        result = queue.get_nowait()
    except Empty:
        result = None

    if result is not None:
        return {
            "method": method,
            "model": model_name,
            "nx": nx,
            "ny": ny,
            "periodic": periodic,
            "seconds": result["seconds"] if result["seconds"] is not None else wall_seconds,
            "wall_seconds": wall_seconds,
            "qubits": result["qubits"],
            "paulis": result["paulis"],
            "found": result["found"],
            "count": result["count"],
            "checked": result.get("checked", None),
            "word": result.get("word", None),
            "timeout": False,
            "error": result["error"],
        }

    qubits, paulis, metadata_error = model_metadata()
    error = f"worker exited with code {process.exitcode}"
    if metadata_error is not None:
        error = f"{error}; metadata error: {metadata_error}"
    return {
        "method": method,
        "model": model_name,
        "nx": nx,
        "ny": ny,
        "periodic": periodic,
        "seconds": wall_seconds,
        "wall_seconds": wall_seconds,
        "qubits": qubits,
        "paulis": paulis,
        "found": False,
        "count": 0,
        "checked": 0,
        "word": None,
        "timeout": False,
        "error": error,
    }


def time_existing_find_clifford(model_name, nx, ny, periodic, timeout_seconds=TIMEOUT_SECONDS):
    return time_method(
        existing_find_worker,
        "existing",
        model_name,
        nx,
        ny,
        periodic,
        timeout_seconds=timeout_seconds,
    )


def time_igraph_find_clifford(model_name, nx, ny, periodic, timeout_seconds=TIMEOUT_SECONDS):
    return time_method(
        igraph_find_worker,
        "igraph",
        model_name,
        nx,
        ny,
        periodic,
        timeout_seconds=timeout_seconds,
    )


def save_results_csv(results, filename="existing_clifford_timing.csv"):
    if not results:
        return
    filename = Path(filename)

    deduped = {}
    for result in results:
        key = (result["method"], result["model"], result["nx"], result["ny"])
        deduped[key] = result
    results = list(deduped.values())

    fieldnames = [
        "method",
        "model",
        "nx",
        "ny",
        "qubits",
        "paulis",
        "found",
        "count",
        "checked",
        "word",
        "timeout",
        "seconds",
        "wall_seconds",
        "error",
    ]
    with open(filename, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow({key: result[key] for key in fieldnames})

    print(f"Saved timing data: {filename}")


def load_results_csv(filename="existing_clifford_timing.csv"):
    filename = Path(filename)
    if not filename.exists():
        return []
    print(f"Data exists, loading: {filename}")

    def parse_int(value, default=None):
        return int(value) if value not in ("", "None", None) else default

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


def qubits_for_job(model_name, nx, ny, periodic):
    if model_name == "toric":
        return 2 * nx * ny if periodic else nx * (ny - 1) + ny * (nx - 1)
    if model_name == "ising_ladder":
        return nx * ny
    if model_name in ("tv_chain", "heisenberg_chain"):
        return ny
    # if model_name == "random_swap_symmetric":
    #     return ny
    if model_name == "fermi_hubbard":
        return 2 * ny
    raise ValueError(f"Unknown model: {model_name}")


def print_loaded_data_summary(results):
    if not results:
        return

    print("Existing data by method/model:")
    for method in sorted({r["method"] for r in results}):
        for model_name in sorted({r["model"] for r in results if r["method"] == method}):
            qubits = sorted(
                r["qubits"]
                for r in results
                if r["method"] == method and r["model"] == model_name and r["qubits"] is not None
            )
            print(f"{method} {model_name}: {qubits}")


def warmup_existing_find():
    print("Warm-up run...")
    pauli_sum = build_model_pauli_sum("ising_ladder", 2, 2, False)
    find_clifford_symmetries(pauli_sum, num_symmetries=1)
    print("Warm-up complete.")


def main():
    ny_range = [2, 4, 6, 8, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500]
    existing_toric_sizes = [(1, ny) for ny in range(2, 51)]
    existing_ising_ladder_sizes = [(2, ny) for ny in ny_range]
    existing_tv_chain_sizes = [(1, ny) for ny in ny_range]
    existing_heisenberg_sizes = [(1, ny) for ny in ny_range]
    existing_fermi_hubbard_sizes = [(1, ny) for ny in ny_range]
    # existing_random_swap_symmetric_sizes = [(1, ny) for ny in ny_range]
    igraph_toric_sizes = [(1, ny) for ny in range(2, 501,10)]
    igraph_ising_ladder_sizes = [(2, ny) for ny in ny_range]
    igraph_tv_chain_sizes = [(1, ny) for ny in ny_range]
    igraph_heisenberg_sizes = [(1, ny) for ny in ny_range]
    igraph_fermi_hubbard_sizes = [(1, ny) for ny in ny_range]
    # igraph_random_swap_symmetric_sizes = [(1, ny) for ny in ny_range]
    jobs = [
        # ("toric", True, existing_toric_sizes, igraph_toric_sizes),
        ("ising_ladder", False, existing_ising_ladder_sizes, igraph_ising_ladder_sizes),
        ("tv_chain", False, existing_tv_chain_sizes, igraph_tv_chain_sizes),
        ("heisenberg_chain", False, existing_heisenberg_sizes, igraph_heisenberg_sizes),
        ("fermi_hubbard", False, existing_fermi_hubbard_sizes, igraph_fermi_hubbard_sizes),
        # (
        #     "random_swap_symmetric",
        #     False,
        #     existing_random_swap_symmetric_sizes,
        #     igraph_random_swap_symmetric_sizes,
        # ),
    ]
    results = load_results_csv(CSV_FILENAME)
    print_loaded_data_summary(results)
    retry_timeout_rows = [r for r in results if r["timeout"]]
    for result in retry_timeout_rows:
        print(
            f"Timeout data exists for {result['method']} {result['model']} "
            f"at {result['qubits']} qubits; retrying this job."
        )
    results = [r for r in results if not r["timeout"]]
    completed = {(r["method"], r["model"], r["nx"], r["ny"]) for r in results}
    timeout_limits = {}
    for result in results:
        if result["timeout"] and result["qubits"] is not None:
            key = (result["method"], result["model"])
            timeout_limits[key] = min(result["qubits"], timeout_limits.get(key, result["qubits"]))

    has_pending = False
    for model_name, periodic, existing_sizes, igraph_sizes in jobs:
        for method, sizes in (("igraph", igraph_sizes),):
            timeout_limit = timeout_limits.get((method, model_name))
            for nx, ny in sizes:
                if (method, model_name, nx, ny) in completed:
                    continue
                qubits = qubits_for_job(model_name, nx, ny, periodic)
                if qubits > MAX_QUBITS:
                    break
                if timeout_limit is not None and qubits >= timeout_limit:
                    continue
                has_pending = True
                break
            if has_pending:
                break
        if has_pending:
            break

    if has_pending:
        # warmup_existing_find()
        pass

    print("method,model,nx,ny,qubits,paulis,found,count,checked,timeout,seconds,wall_seconds,word,error")
    for model_name, periodic, existing_sizes, igraph_sizes in jobs:
        method_jobs = [
            ("igraph", time_igraph_find_clifford, igraph_sizes),
            # ("existing", time_existing_find_clifford, existing_sizes),
        ]
        for method, timer, sizes in method_jobs:
            for nx, ny in sizes:
                if (method, model_name, nx, ny) in completed:
                    qubits = qubits_for_job(model_name, nx, ny, periodic)
                    print(f"Data exists for {method} {model_name} at {qubits} qubits; skipping.")
                    continue

                timeout_limit = timeout_limits.get((method, model_name))
                qubits = qubits_for_job(model_name, nx, ny, periodic)
                if qubits > MAX_QUBITS:
                    break
                if timeout_limit is not None and qubits >= timeout_limit:
                    break

                result = timer(model_name, nx, ny, periodic)
                results.append(result)
                completed.add((result["method"], result["model"], result["nx"], result["ny"]))
                print(
                    f"{result['method']},{result['model']},{result['nx']},{result['ny']},"
                    f"{result['qubits']},{result['paulis']},{result['found']},{result['count']},"
                    f"{result['checked']},{result['timeout']},{result['seconds']:.6f},"
                    f"{result['wall_seconds']:.6f},{result['word']},{result['error']}"
                )
                if result["timeout"]:
                    print(f"Timeout for {result['method']} {model_name} at nx={nx}, ny={ny}.")
                    if result["qubits"] is not None:
                        timeout_limits[(result["method"], result["model"])] = result["qubits"]
                    save_results_csv(results, CSV_FILENAME)
                    break

                save_results_csv(results, CSV_FILENAME)

    save_results_csv(results, CSV_FILENAME)
    print(f"Plot with: python tests\\plot_clifford_timing.py")


if __name__ == "__main__":
    main()
