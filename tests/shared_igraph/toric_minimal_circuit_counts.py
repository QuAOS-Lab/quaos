"""Count minimal GF(2) dependency circuits for small toric-code Hamiltonians.

This intentionally does not build the base or augmented graph. It only builds
the toric-code PauliSum and counts the minimal dependent supports used by the
circuit-augmentation graph scripts.
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from sympleq.core.graphs.graph_automorphism_circuits import (
    _gf2_nullspace_basis,
    augment_S_with_circuits,
    extract_circuits_from_nullspace_gf2,
)
from sympleq.models.toric_code import ToricCode


CG = 0.1
CZ = 1.0
CX = 2.0
NX = 2
PERIODIC = True
MAX_NULLITY = 1001
MAX_CIRCUITS = 1_000_000
INCIDENCE_LABEL = 2


def gf2_rank(matrix: np.ndarray) -> int:
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


def converted_graph_vertex_count_from_s_mod(S_mod: np.ndarray) -> int:
    """Return the vertex count build_subdivision_graph_from_s_mod would create."""
    S_mod = np.asarray(S_mod)
    n = S_mod.shape[0]

    if np.all((S_mod == 0) | (S_mod == 1)):
        return n

    upper_rows, upper_cols = np.triu_indices(n, k=1)
    upper_values = S_mod[upper_rows, upper_cols].astype(int)
    values, counts = np.unique(upper_values, return_counts=True)
    no_edge_value = int(values[np.argmax(counts)])
    subdivision_vertices = int(np.count_nonzero(upper_values != no_edge_value))
    return n + subdivision_vertices


def toric_qudit_count(nx: int, ny: int, periodic: bool) -> int:
    if periodic:
        return 2 * nx * ny
    return nx * (ny - 1) + ny * (nx - 1)


def summarize_minimal_circuits(ny: int) -> dict:
    n_qubits = toric_qudit_count(NX, ny, PERIODIC)
    pauli_sum = ToricCode(
        Nx=NX,
        Ny=ny,
        c_x=CX,
        c_z=CZ,
        c_g=CG,
        periodic=PERIODIC,
        d=2,
    ).hamiltonian()
    pauli_sum.weight_to_phase()

    rank = gf2_rank(pauli_sum.tableau)
    nullity = pauli_sum.tableau.shape[0] - rank
    nullspace_basis = _gf2_nullspace_basis(pauli_sum.tableau.T)
    circuits = extract_circuits_from_nullspace_gf2(
        pauli_sum.tableau,
        max_nullity=MAX_NULLITY,
        max_circuits=MAX_CIRCUITS,
    )
    S_mod = pauli_sum.symplectic_product_matrix()
    S_aug = augment_S_with_circuits(S_mod, circuits, incidence_label=INCIDENCE_LABEL)
    circuit_sizes = [len(circuit) for circuit in circuits]

    return {
        "n_qubits": n_qubits,
        "nx": NX,
        "ny": ny,
        "periodic": PERIODIC,
        "cg": CG,
        "cz": CZ,
        "cx": CX,
        "n_paulis": pauli_sum.n_paulis(),
        "pauli_n_qudits": pauli_sum.n_qudits(),
        "tableau_shape": tuple(pauli_sum.tableau.shape),
        "lcm": int(pauli_sum.lcm),
        "gf2_rank": rank,
        "gf2_nullity": nullity,
        "nullspace_basis_count": len(nullspace_basis),
        "minimal_circuits": len(circuits),
        "circuit_size_histogram": dict(sorted(Counter(circuit_sizes).items())),
        "first_circuits": circuits[:10],
        "base_input_vertices": int(S_mod.shape[0]),
        "base_graph_vertices": converted_graph_vertex_count_from_s_mod(S_mod),
        "augmented_input_vertices": int(S_aug.shape[0]),
        "augmented_graph_vertices": converted_graph_vertex_count_from_s_mod(S_aug),
    }


def main() -> None:
    for ny in (7, 8):
        summary = summarize_minimal_circuits(ny)
        print("===")
        print(f"model=toric_code nx={summary['nx']} ny={summary['ny']} periodic={summary['periodic']}")
        print(f"coefficients cg={summary['cg']} cz={summary['cz']:g} cx={summary['cx']:g}")
        print(f"n_qubits={summary['n_qubits']} pauli_n_qudits={summary['pauli_n_qudits']}")
        print(f"n_paulis={summary['n_paulis']} tableau_shape={summary['tableau_shape']} lcm={summary['lcm']}")
        print(f"gf2_rank={summary['gf2_rank']} gf2_nullity={summary['gf2_nullity']}")
        print(f"nullspace_basis_count={summary['nullspace_basis_count']}")
        print(f"minimal_circuits={summary['minimal_circuits']}")
        print(f"circuit_size_histogram={summary['circuit_size_histogram']}")
        print(f"first_circuits={summary['first_circuits']}")
        print(f"base_input_vertices={summary['base_input_vertices']}")
        print(f"base_graph_vertices={summary['base_graph_vertices']}")
        print(f"augmented_input_vertices={summary['augmented_input_vertices']}")
        print(f"augmented_graph_vertices={summary['augmented_graph_vertices']}")


if __name__ == "__main__":
    main()
