"""Build a circuit-augmented graph from a Pauli Hamiltonian.

This is only a graph-construction sandbox. It does not run the Clifford graph
automorphism search. The flow is:

1. Build a small PauliSum model.
2. Build the ordinary colored graph from the symplectic-product matrix.
3. Extract GF(2) matroid dependency sets from the Pauli tableau.
4. Add one dependency node per dependency set and build the augmented graph.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from sympleq.core.graphs.graph_automorphism_circuits import (  # noqa: E402
    _gf2_nullspace_basis,
    augment_S_with_circuits,
    extract_circuits_from_nullspace_gf2,
)
from sympleq.core.graphs.graph_builder import (  # noqa: E402
    build_subdivision_graph_from_s_mod,
)
from sympleq.models import (  # noqa: E402
    ToricCode,
    all_to_all_heisenberg_hamiltonian,
    ising_2d_hamiltonian,
)


MODEL = "toric"
NX = 1
NY = 2
D = 2
PERIODIC = True
GAUGE = 0.1
FIELD = 0.7
DELTA_Z = 0.5
MAX_NULLITY = 1001
MAX_CIRCUITS = 5000
INCIDENCE_LABEL = 2


def _histogram(values: np.ndarray | list[int]) -> dict[int, int]:
    return dict(sorted(Counter(int(v) for v in values).items()))


def compress_vertex_colors(coeffs) -> list[int]:
    _, inverse = np.unique(np.asarray(coeffs), return_inverse=True)
    return inverse.astype(int).tolist()


def gf2_rank(matrix) -> int:
    matrix = (np.asarray(matrix, dtype=np.uint8) & 1).copy()
    rows, cols = matrix.shape
    rank = 0

    for col in range(cols):
        pivot = None
        for row in range(rank, rows):
            if matrix[row, col]:
                pivot = row
                break
        if pivot is None:
            continue

        if pivot != rank:
            matrix[[rank, pivot]] = matrix[[pivot, rank]]

        for row in range(rows):
            if row != rank and matrix[row, col]:
                matrix[row, :] ^= matrix[rank, :]

        rank += 1
        if rank == rows:
            break

    return rank


def build_model():
    if MODEL == "toric":
        return ToricCode(
            Nx=NX,
            Ny=NY,
            c_x=1.0,
            c_z=1.0,
            c_g=GAUGE,
            periodic=PERIODIC,
            d=D,
        ).hamiltonian()

    if MODEL == "ising":
        return ising_2d_hamiltonian(
            n_x=NX,
            n_y=NY,
            J_zz=1.0,
            h_x=FIELD,
            periodic=PERIODIC,
        )

    if MODEL == "heisenberg":
        return all_to_all_heisenberg_hamiltonian(
            n=NY,
            J=1.0,
            delta_z=DELTA_Z,
        )

    raise ValueError(f"Unknown model: {MODEL}")


def build_graphs(pauli_sum) -> dict:
    pauli = pauli_sum.copy()
    pauli.weight_to_phase()

    S_mod = pauli.symplectic_product_matrix()
    base_vertex_colors = compress_vertex_colors(pauli.weights)
    base_graph, base_graph_colors = build_subdivision_graph_from_s_mod(
        S_mod,
        base_vertex_colors,
    )

    dependency_sets: list[list[int]] = []
    nullspace_basis = []
    if int(pauli.lcm) == 2:
        nullity = pauli.tableau.shape[0] - gf2_rank(pauli.tableau)
        nullspace_basis = _gf2_nullspace_basis(pauli.tableau.T)
        if nullity > MAX_NULLITY:
            print(
                "WARNING: dependency-set extraction skipped because nullity="
                f"{nullity} exceeds MAX_NULLITY={MAX_NULLITY}."
            )

        dependency_sets = extract_circuits_from_nullspace_gf2(
            pauli.tableau,
            max_nullity=MAX_NULLITY,
            max_circuits=MAX_CIRCUITS,
        )
        if len(dependency_sets) >= MAX_CIRCUITS:
            print(
                "WARNING: dependency-set extraction reached MAX_CIRCUITS="
                f"{MAX_CIRCUITS}; the augmented graph may be truncated."
            )

    S_aug = S_mod
    augmented_vertex_colors = base_vertex_colors
    if dependency_sets:
        S_aug = augment_S_with_circuits(
            S_mod,
            dependency_sets,
            incidence_label=INCIDENCE_LABEL,
        )
        dependency_color_offset = max(base_vertex_colors, default=-1) + 1
        dependency_vertex_colors = [
            dependency_color_offset + len(dependency_set)
            for dependency_set in dependency_sets
        ]
        augmented_vertex_colors = base_vertex_colors + dependency_vertex_colors

    augmented_graph, augmented_graph_colors = build_subdivision_graph_from_s_mod(
        S_aug,
        augmented_vertex_colors,
    )

    return {
        "pauli": pauli,
        "S_mod": S_mod,
        "S_aug": S_aug,
        "gf2_nullity": pauli.tableau.shape[0] - gf2_rank(pauli.tableau)
        if int(pauli.lcm) == 2
        else None,
        "gf2_nullspace_basis": nullspace_basis,
        "dependency_sets": dependency_sets,
        "base_vertex_colors": base_vertex_colors,
        "base_graph": base_graph,
        "base_graph_colors": base_graph_colors,
        "augmented_vertex_colors": augmented_vertex_colors,
        "augmented_graph": augmented_graph,
        "augmented_graph_colors": augmented_graph_colors,
    }


def print_report(result: dict) -> None:
    pauli = result["pauli"]
    dependency_sets = result["dependency_sets"]
    base_graph = result["base_graph"]
    augmented_graph = result["augmented_graph"]

    print("=== Circuit augmented graph build ===")
    print(f"terms={pauli.n_paulis()} qudits={pauli.n_qudits()} lcm={pauli.lcm}")
    print(f"S_shape={tuple(result['S_mod'].shape)} S_aug_shape={tuple(result['S_aug'].shape)}")
    print(f"gf2_nullity={result['gf2_nullity']}")
    print("gf2_nullspace_basis=")
    for basis_vector in result["gf2_nullspace_basis"]:
        print(basis_vector.astype(int).tolist())
    print(f"dependency_sets_found={len(dependency_sets)}")

    if dependency_sets:
        sizes = [len(dependency_set) for dependency_set in dependency_sets]
        print(f"dependency_size_histogram={_histogram(sizes)}")
        print(f"first_dependency_sets={dependency_sets[:5]}")
    else:
        print("augmentation_status=no dependency nodes added")

    print(
        "base_graph="
        f"vertices:{base_graph.vcount()} edges:{base_graph.ecount()} "
        f"color_hist:{_histogram(result['base_graph_colors'])}"
    )
    print(
        "augmented_graph="
        f"vertices:{augmented_graph.vcount()} edges:{augmented_graph.ecount()} "
        f"color_hist:{_histogram(result['augmented_graph_colors'])}"
    )
    print(f"base_input_vertex_color_hist={_histogram(result['base_vertex_colors'])}")
    print(f"augmented_input_vertex_color_hist={_histogram(result['augmented_vertex_colors'])}")


def main() -> None:
    pauli_sum = build_model()
    result = build_graphs(pauli_sum)
    print_report(result)


if __name__ == "__main__":
    main()
