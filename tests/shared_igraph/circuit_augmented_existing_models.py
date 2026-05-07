"""Run the augmented-generator leaf-check counters on existing timing models."""

from __future__ import annotations

import csv
import multiprocessing as mp
from collections import Counter
from pathlib import Path
from queue import Empty
from time import perf_counter

import numpy as np
import galois


from sympleq.core.circuits import Circuit, GATES, Gate
from sympleq.core.circuits.phase_correction import solve_phase_vector_h_from_residual
from sympleq.core.circuits.target import find_map_to_target_pauli_sum, get_phase_vector
from sympleq.core.paulis import PauliSum
from sympleq.core.graphs.graph_automorphism_code import check_code_automorphism
from sympleq.core.graphs.graph_automorphism_leaf import _solve_qubit_phase_family_correction
from sympleq.core.graphs.graph_automorphism_search import prepare_clifford_ga_search

from sympleq.models.fermi_hubbard import disordered_tv_chain_model
from sympleq.models.fermionic_chain import fermionic_chain_hamiltonian
from sympleq.models.heisenberg import (
    all_to_all_heisenberg_hamiltonian,
    heisenberg_2d_hamiltonian,
    modified_heisenberg_ladder_hamiltonian,
)
from sympleq.models.Ising import (
    ising_2d_hamiltonian,
    ising_chain_hamiltonian,
    ising_lower_triangular_hamiltonian,
    modified_ising_ladder_hamiltonian,
)
from sympleq.models.pxp import pxp_model
from sympleq.models.random_hamiltonian import random_gate_symmetric_hamiltonian
from sympleq.models.syk import syk4_majorana_hamiltonian
from sympleq.models.toric_code import ToricCode


from circuit_augmented_graph_builder import build_graphs


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
CSV_FILENAME = DATA_DIR / "circuit_augmented_existing_models.csv"
TIMEOUT_SECONDS = 100
MAX_QUBITS = 100

CHAIN_SIZES = [2, 4, 6, 8, 10, 20, 40, 60, 80, 100]
SMALL_CHAIN_SIZES = [2, 4, 6, 8, 10, 20, 40]
SQUARE_SIDES = [2, 3, 4, 5, 6, 8, 10]

JOBS = (
    [("toric_code", 1, ny, True) for ny in [2, 5, 10, 20, 30, 40, 50]]
    + [("ising_chain", 1, ny, False) for ny in CHAIN_SIZES]
    + [("open_ising_chain", 1, ny, False) for ny in CHAIN_SIZES]
    + [("square_ising", side, side, False) for side in SQUARE_SIDES]
    + [("ising_ladder", 2, ny, False) for ny in [2, 4, 6, 8, 10, 20, 30, 40, 50]]
    + [("triangle_ising", 1, side, False) for side in [2, 3, 4, 5, 8, 10, 12, 14]]
    + [("tv_chain", 1, ny, False) for ny in SMALL_CHAIN_SIZES]
    + [("heisenberg", 1, ny, False) for ny in SMALL_CHAIN_SIZES]
    + [("heisenberg_2d", 2, ny, False) for ny in [2, 3, 4, 5, 10, 20, 30, 40, 50]]
    + [("modified_heisenberg_ladder", ny, 2, False) for ny in [2, 3, 4, 5, 10, 20, 30, 40, 50]]
    + [("pxp", 1, ny, False) for ny in [3, 4, 6, 8, 10, 20, 40, 60, 80, 100]]
    + [("syk", 1, ny, False) for ny in SMALL_CHAIN_SIZES]
    + [("fermionic_chain", 1, ny, False) for ny in SMALL_CHAIN_SIZES]
    + [("random_fermionic_chain", 1, ny, False) for ny in SMALL_CHAIN_SIZES]
    + [("fermi_hubbard", 1, ny, False) for ny in [1, 2, 4, 8, 12, 16, 20, 25]]
    + [("random_swap_symmetric", 1, ny, False) for ny in CHAIN_SIZES]
)


def generator_leaf_outcome_counts(pauli_sum, generators, leaf_ctx) -> Counter:
    n_paulis = pauli_sum.n_paulis()
    identity_pauli_perm = tuple(range(n_paulis))
    outcome_counts = Counter()

    for generator in generators:
        pi_tuple = permutation_to_tuple(generator)[:n_paulis]
        if pi_tuple == identity_pauli_perm:
            outcome_counts["identity_skipped"] += 1
            continue

        _passed_checks, failed_check, gate = leaf_check_stages(
            np.asarray(pi_tuple, dtype=np.int64),
            leaf_ctx,
        )
        if gate is not None:
            outcome_counts["passed_all_checks"] += 1
        else:
            outcome_counts[f"failed_{failed_check}"] += 1

    return outcome_counts


def qubits_for_job(model_name, nx, ny, periodic) -> int:
    model_key = str(model_name).strip().lower()

    if model_key in ("toric", "toric_code"):
        return 2 * nx * ny if periodic else nx * (ny - 1) + ny * (nx - 1)
    if model_key in (
        "ising_ladder",
        "square_ising",
        "ising_2d",
        "heisenberg_2d",
        "square_heisenberg",
        "modified_heisenberg_ladder",
        "heisenberg_ladder_xxz",
    ):
        return nx * ny
    if model_key in (
        "ising_chain",
        "open_ising_chain",
        "ising_chain_open",
        "tv_chain",
        "heisenberg",
        "heisenberg_chain",
        "all_to_all_heisenberg",
        "pxp",
        "pxp_model",
        "syk",
        "syk4",
        "syk4_majorana",
        "fermionic_chain",
        "random_fermionic_chain",
        "random_swap_symmetric",
    ):
        return ny
    if model_key in ("triangle_ising", "triangle ising"):
        return ny * (ny + 1) // 2
    if model_key == "fermi_hubbard":
        return 4 * ny
    raise ValueError(f"Unknown model: {model_name}")

def build_model_pauli_sum(model_name, nx, ny, periodic):
    model_key = str(model_name).strip().lower()
    if model_key in ("toric", "toric_code"):
        return ToricCode(Nx=nx, Ny=ny, c_x=1.0, c_z=2.0, c_g=0.1, periodic=periodic).hamiltonian()
    if model_key == "ising_chain":
        return ising_chain_hamiltonian(n_spins=ny, J_zz=1.0, h_x=0.5, periodic=periodic)
    if model_key in ("open ising chain", "open_ising_chain", "ising_chain_open"):
        return ising_chain_hamiltonian(n_spins=ny, J_zz=1.0, h_x=0.5, periodic=False)
    if model_key in ("square ising", "square_ising", "ising_2d"):
        return ising_2d_hamiltonian(n_x=nx, n_y=ny, J_zz=1.0, h_x=0.5, periodic=periodic)
    if model_key in ("ising ladder", "ising_ladder"):
        return ising_2d_hamiltonian(n_x=nx, n_y=ny, J_zz=1.0, h_x=0.5, periodic=periodic)
    if model_key in ("modified_ising_ladder", "modified ising ladder"):
        return modified_ising_ladder_hamiltonian(n_x=nx, n_y=ny, J_zz=1.0, h_x=0.5)
    if model_key in ("triangle ising", "triangle_ising"):
        return ising_lower_triangular_hamiltonian(L=ny, J_zz=1.0, h_x=0.5)
    if model_key == "tv_chain":
        return disordered_tv_chain_model(
            n_sites=ny,
            tunneling=1.0,
            coulomb=4.0,
            disorder_strength=0.5,
            periodic=periodic,
            seed=10_000 + ny,
        )
    if model_key in ("heisenberg_chain", "heisenberg", "all_to_all_heisenberg"):
        return all_to_all_heisenberg_hamiltonian(
            n=ny,
            J=1.0,
            delta_z=0.5,
        )
    if model_key in ("heisenberg_2d", "square_heisenberg"):
        return heisenberg_2d_hamiltonian(n_x=nx, n_y=ny, J=1.0, h_z=0.5, periodic=periodic)
    if model_key in ("modified_heisenberg_ladder", "heisenberg_ladder_xxz", "modified heisenberg ladder"):
        return modified_heisenberg_ladder_hamiltonian(n_x=nx, n_y=ny, J=1.0, h_z=0.5)
    if model_key in ("pxp", "pxp_model"):
        return pxp_model(n_qubits=ny, coupling=1.0, periodic=periodic, normalized_projectors=False)
    if model_key in ("syk", "syk4", "syk4_majorana"):
        return syk4_majorana_hamiltonian(n_qubits=ny, J_scale=1.0, seed=1_234_567 + 1009 * ny)
    if model_key == "fermionic_chain":
        return fermionic_chain_hamiltonian(n=ny, J=1.0, V=1.0, D_vec=np.zeros(ny), periodic=periodic)
    if model_key == "random_fermionic_chain":
        rng = np.random.default_rng(20_000 + ny)
        return fermionic_chain_hamiltonian(
            n=ny,
            J=1.0,
            V=0.5,
            D_vec=rng.uniform(-0.5, 0.5, size=ny),
            periodic=periodic,
        )
    if model_key == "fermi_hubbard":
        from sympleq.models.fermi_hubbard import fermi_hubbard_model
        return fermi_hubbard_model(
            x_dimension=ny,
            y_dimension=2,
            tunneling=1.0,
            coulomb=4.0,
            chemical_potential=0.0,
            periodic=False,
            spinless=False,
        )
    if model_key == "random_swap_symmetric":
        p = 2
        n_qudits = int(ny)
        dims = [p] * n_qudits
        all_qudits = tuple(range(n_qudits))
        swap = Circuit.from_gates_and_qudits(dims, [GATES.SWAP], [(0, 1)]).composite_gate()
        return random_gate_symmetric_hamiltonian(
            swap,
            p,
            all_qudits,
            n_qudits,
            n_paulis=2 * n_qudits,
            weight_mode="uniform",
            scrambled=True,
        )
    raise ValueError(f"Unknown model: {model_name}")

def run_one(model_name, nx, ny, periodic, pauli_sum) -> dict:
    start = perf_counter()

    print(
        f"\nStarting model={model_name} nx={nx} ny={ny} periodic={periodic} "
        f"qubits={pauli_sum.n_qudits()} paulis={pauli_sum.n_paulis()}"
    )
    graph_result = build_graphs(pauli_sum)

    graph = graph_result["augmented_graph"]
    graph_colors = graph_result["augmented_graph_colors"]
    automorphism_group = graph.automorphism_group(color=graph_colors)
    generators = (
        automorphism_group.generators
        if hasattr(automorphism_group, "generators")
        else automorphism_group
    )

    prepared = prepare_clifford_ga_search(
        pauli_sum,
        extra_column_invariants="none",
        p2_bitset="auto",
        color_mode="none",
        max_wl_rounds=0,
    )
    outcome_counts = generator_leaf_outcome_counts(
        pauli_sum,
        generators,
        prepared.leaf_ctx,
    )
    clifford_symmetry_seconds = perf_counter() - start
    num_gen_passes = int(outcome_counts.get("passed_all_checks", 0))
    num_gen_fails_at_leaf = int(
        sum(
            count
            for outcome, count in outcome_counts.items()
            if outcome.startswith("failed_")
        )
    )

    return {
        "model": model_name,
        "nx": nx,
        "ny": ny,
        "periodic": periodic,
        "qubits": pauli_sum.n_qudits(),
        "paulis": pauli_sum.n_paulis(),
        "gf2_nullity": graph_result["gf2_nullity"],
        "dependency_sets": len(graph_result["dependency_sets"]),
        "num_circuits": len(graph_result["dependency_sets"]),
        "nullity": graph_result["gf2_nullity"],
        "unaugmented_vertices": graph_result["base_graph"].vcount(),
        "augmented_vertices": graph.vcount(),
        "num_generators": len(generators),
        "clifford_symmetry_seconds": clifford_symmetry_seconds,
        "num_gen_passes": num_gen_passes,
        "num_gen_fails_at_leaf": num_gen_fails_at_leaf,
        "outcome_counts": outcome_counts,
    }


def run_one_worker(model_name, nx, ny, periodic, pauli_sum, queue) -> None:
    try:
        queue.put(run_one(model_name, nx, ny, periodic, pauli_sum))
    except Exception as exc:
        queue.put(
            {
                "model": model_name,
                "nx": nx,
                "ny": ny,
                "periodic": periodic,
                "error": f"{type(exc).__name__}: {exc}",
            }
        )


def run_one_with_timeout(model_name, nx, ny, periodic, pauli_sum, timeout_seconds=TIMEOUT_SECONDS) -> dict:
    queue = mp.Queue()
    process = mp.Process(
        target=run_one_worker,
        args=(model_name, nx, ny, periodic, pauli_sum, queue),
    )

    wall_start = perf_counter()
    process.start()
    process.join(timeout_seconds)
    wall_seconds = perf_counter() - wall_start

    if process.is_alive():
        process.terminate()
        process.join()
        return {
            "model": model_name,
            "nx": nx,
            "ny": ny,
            "periodic": periodic,
            "qubits": qubits_for_job(model_name, nx, ny, periodic),
            "clifford_symmetry_seconds": timeout_seconds,
            "timeout": True,
            "error": f"timeout after {timeout_seconds} seconds",
        }

    try:
        result = queue.get_nowait()
    except Empty:
        result = {
            "model": model_name,
            "nx": nx,
            "ny": ny,
            "periodic": periodic,
            "qubits": qubits_for_job(model_name, nx, ny, periodic),
            "error": f"worker exited with code {process.exitcode}",
        }

    result["wall_seconds"] = wall_seconds
    result.setdefault("timeout", False)
    return result


def print_result(result: dict) -> None:
    if result.get("error") is not None:
        print(
            f"\nmodel={result['model']} nx={result['nx']} ny={result['ny']} "
            f"periodic={result['periodic']} qubits={result.get('qubits')} "
            f"error={result['error']}"
        )
        return

    print(
        f"\nmodel={result['model']} nx={result['nx']} ny={result['ny']} "
        f"periodic={result['periodic']}"
    )
    print(
        f"qubits={result['qubits']} paulis={result['paulis']} "
        f"gf2_nullity={result['gf2_nullity']} dependency_sets={result['dependency_sets']} "
        f"augmented_vertices={result['augmented_vertices']} "
        f"num_generators={result['num_generators']} "
        f"clifford_symmetry_seconds={result['clifford_symmetry_seconds']:.6f} "
        f"num_gen_passes={result['num_gen_passes']} "
        f"num_gen_fails_at_leaf={result['num_gen_fails_at_leaf']}"
    )
    if result["outcome_counts"]:
        for outcome, count in sorted(result["outcome_counts"].items()):
            print(f"{outcome}={count}")
    else:
        print("no_generators=0")


def save_results_csv(results, filename=CSV_FILENAME) -> None:
    fieldnames = [
        "model",
        "nx",
        "ny",
        "periodic",
        "qubits",
        "paulis",
        "gf2_nullity",
        "dependency_sets",
        "num_circuits",
        "nullity",
        "unaugmented_vertices",
        "augmented_vertices",
        "num_generators",
        "clifford_symmetry_seconds",
        "wall_seconds",
        "num_gen_passes",
        "num_gen_fails_at_leaf",
        "timeout",
        "error",
    ]

    filename.parent.mkdir(parents=True, exist_ok=True)
    with open(filename, "w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow({field: result.get(field) for field in fieldnames})

    print(f"\nSaved CSV: {filename}")


def load_results_csv(filename=CSV_FILENAME) -> list[dict]:
    filename = Path(filename)
    if not filename.exists():
        return []

    with open(filename, newline="") as csv_file:
        return list(csv.DictReader(csv_file))


def result_key(result: dict) -> tuple[str, int] | None:
    model = result.get("model")
    qubits = result.get("qubits")
    if model in ("", None) or qubits in ("", "None", None):
        return None
    return str(model), int(qubits)


def is_timeout_result(result: dict) -> bool:
    return str(result.get("timeout")).strip().lower() == "true"

def permutation_to_tuple(permutation) -> tuple[int, ...]:
    if hasattr(permutation, "mapping"):
        return tuple(permutation.mapping)
    return tuple(permutation)


def leaf_check_stages(pi: np.ndarray, leaf_ctx) -> tuple[list[str], str | None, object | None]:
    passed = []

    if np.array_equal(pi, leaf_ctx.identity_perm):
        return passed, "identity", None

    if not np.array_equal(leaf_ctx.S_mod[np.ix_(pi, pi)], leaf_ctx.S_mod):
        return passed, "commutation_matrix", None
    passed.append("commutation_matrix")

    if not check_code_automorphism(
        leaf_ctx.G,
        leaf_ctx.basis_order,
        leaf_ctx.labels,
        pi,
        leaf_ctx.G_mod2,
    ):
        return passed, "linear_code", None
    passed.append("linear_code")

    tgt_idx = pi[leaf_ctx.basis_indices]
    H_basis_tgt = PauliSum.from_tableau(
        leaf_ctx.base_tableau[tgt_idx],
        leaf_ctx.pauli_sum.dimensions,
        weights=leaf_ctx.base_weights[tgt_idx],
    )
    H_basis_tgt.set_phases(np.array(leaf_ctx.base_phases[tgt_idx], dtype=int, copy=True))
    H_basis_src = leaf_ctx.basis_source_ps

    if leaf_ctx.p == 2 and leaf_ctx.basis_src_inv_gf2 is not None:
        T = (leaf_ctx.base_tableau[tgt_idx] & 1).astype(np.uint8, copy=False)
        F = (leaf_ctx.basis_src_inv_gf2 @ T) & 1
        F = np.asarray(F, dtype=int)
    elif leaf_ctx.p != 2 and leaf_ctx.basis_src_inv_gfp is not None:
        GF = galois.GF(int(leaf_ctx.p))
        T = GF(leaf_ctx.base_tableau[tgt_idx] % leaf_ctx.p)
        F_gf = leaf_ctx.basis_src_inv_gfp @ T
        F = (np.asarray(F_gf, dtype=int) % leaf_ctx.p).astype(int, copy=False)
    else:
        try:
            F, _h0, _, _ = find_map_to_target_pauli_sum(H_basis_src, H_basis_tgt)
        except Exception:
            return passed, "symplectic_map", None
    passed.append("symplectic_map")

    h0 = get_phase_vector(F.T, int(leaf_ctx.pauli_sum.dimensions[0]))
    SG_F = Gate("Symmetry", F.T, np.asarray(h0, dtype=int))
    H_full_tg = leaf_ctx.pauli_sum.copy()[pi]
    H_full_F = SG_F.act(leaf_ctx.pauli_sum, tuple(range(leaf_ctx.n_qudits)))
    delta = (H_full_tg.phases - H_full_F.phases) % leaf_ctx.two_lcm

    if leaf_ctx.p == 2:
        h_lin = _solve_qubit_phase_family_correction(leaf_ctx.base_tableau, delta)
    else:
        h_lin = solve_phase_vector_h_from_residual(
            leaf_ctx.base_tableau,
            delta,
            leaf_ctx.pauli_sum.dimensions,
            debug=False,
            row_basis_cache=leaf_ctx.row_basis_cache,
        )

    if h_lin is None:
        return passed, "phase_correction", None
    passed.append("phase_correction")

    h0_mod = np.asarray(h0, dtype=int) % leaf_ctx.two_lcm
    h_lin_mod = np.asarray(h_lin, dtype=int) % leaf_ctx.two_lcm
    h_tot = (h0_mod + h_lin_mod) % leaf_ctx.two_lcm
    gate = Gate("Symmetry", F.T, h_tot)

    H_out_cf = gate.act(leaf_ctx.pauli_sum, tuple(range(leaf_ctx.n_qudits))).to_standard_form()
    H_out_cf.weight_to_phase()

    if not np.array_equal(H_out_cf.tableau, leaf_ctx.ref_tableau):
        return passed, "final_tableau_verification", None
    passed.append("final_tableau_verification")

    if not np.all((leaf_ctx.ref_phases - H_out_cf.phases) % leaf_ctx.two_lcm == 0):
        return passed, "final_phase_verification", None
    passed.append("final_phase_verification")

    if not np.all(np.isclose(H_out_cf.weights, leaf_ctx.ref_weights, atol=1e-8, rtol=0)):
        return passed, "final_weight_verification", None
    passed.append("final_weight_verification")

    return passed, None, gate


def main() -> None:
    results = load_results_csv()
    completed = {
        key
        for key in (result_key(result) for result in results)
        if key is not None
    }
    timed_out_models = {
        str(result["model"])
        for result in results
        if result.get("model") not in ("", None) and is_timeout_result(result)
    }
    if results:
        print(f"Loaded {len(results)} existing CSV rows from {CSV_FILENAME}")

    for model_name, nx, ny, periodic in JOBS:
        if model_name in timed_out_models:
            print(f"Skipping model={model_name}; earlier row timed out.")
            continue

        qubits = qubits_for_job(model_name, nx, ny, periodic)
        if qubits > MAX_QUBITS:
            continue
        if (model_name, qubits) in completed:
            print(f"Skipping existing row for model={model_name} qubits={qubits}")
            continue

        pauli_sum = build_model_pauli_sum(model_name, nx, ny, periodic)
        print(
            f"\nQueued model={model_name} nx={nx} ny={ny} periodic={periodic} "
            f"qubits={pauli_sum.n_qudits()} paulis={pauli_sum.n_paulis()} "
            f"cutoff_seconds={TIMEOUT_SECONDS}"
        )
        result = run_one_with_timeout(model_name, nx, ny, periodic, pauli_sum)
        result.setdefault("error", None)
        print_result(result)
        results.append(result)
        key = result_key(result)
        if key is not None:
            completed.add(key)
        save_results_csv(results)
        if result.get("timeout"):
            timed_out_models.add(model_name)


if __name__ == "__main__":
    main()
