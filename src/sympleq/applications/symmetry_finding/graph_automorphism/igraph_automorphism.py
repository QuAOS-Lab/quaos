from __future__ import annotations

from typing import Any

import numpy as np

from sympleq.core.circuits import Gate
from sympleq.core.circuits.helpers_solve_from_target import get_phase_vector
from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.core.paulis import PauliSum
from sympleq.core.phase_correction import solve_phase_vector_h_from_residual

from sympleq.core.paulis._typing import TableauType
from sympleq.core.circuits.helpers_solve_from_target import map_tableau_to_target_tableau

from .graph_builder import build_subdivision_graph_from_s_mod


def _permutation_to_tuple(permutation: Any) -> tuple[int, ...]:
    if hasattr(permutation, "mapping"):
        return tuple(int(x) for x in permutation.mapping)
    return tuple(int(x) for x in permutation)


def _compose_permutations(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[i] for i in right)


def _extract_pauli_permutation(full_perm: tuple[int, ...], n_paulis: int) -> tuple[int, ...] | None:
    pi = tuple(full_perm[:n_paulis])
    if any(x < 0 or x >= n_paulis for x in pi):
        return None
    if len(set(pi)) != n_paulis:
        return None
    return pi


def _generated_group_permutations(generators: list[Any], n_vertices: int):
    generators = [_permutation_to_tuple(g) for g in generators]
    identity = tuple(range(n_vertices))
    seen = {identity}
    queue_items = [identity]

    while queue_items:
        current = queue_items.pop(0)
        for generator in generators:
            for candidate in (
                _compose_permutations(generator, current),
                _compose_permutations(current, generator),
            ):
                if candidate in seen:
                    continue
                seen.add(candidate)
                queue_items.append(candidate)
                yield candidate


def _column_invariants(pauli_sum: PauliSum, mode: str) -> np.ndarray | None:
    toks = {t.strip().lower() for t in mode.replace(",", "+").split("+") if t.strip()}
    toks.discard("none")
    if not toks:
        return None
    if toks != {"hist"}:
        raise ValueError("igraph extra_invs currently supports only 'none' or 'hist'.")

    tableau = np.asarray(pauli_sum.tableau, dtype=int) % int(pauli_sum.lcm)
    p = int(pauli_sum.lcm)
    features = np.zeros((tableau.shape[0], min(p, 16)), dtype=np.int64)
    for i, row in enumerate(tableau):
        counts = np.bincount(row, minlength=p)
        features[i, : min(p, 16)] = counts[: min(p, 16)]
    return features


# This is a legacy function from find_symplectic.py,
# kept for compatibility with the igraph automorphism code.
# modified to use the helpers_solve_from_target.map_tableau_to_target_tableau
# function instead of the original map_paulisum_to_target_tableau.

def symplectic_from_pauli_permutation(
    paulisum_tableau: TableauType,
    permutation: TableauType,
    p: int = 2,
) -> TableauType:
    tableau = np.asarray(paulisum_tableau, dtype=int) % p
    pi = np.asarray(permutation, dtype=np.int64).reshape(-1)

    if tableau.ndim != 2:
        raise ValueError("paulisum_tableau must be a 2D tableau.")
    if pi.shape[0] != tableau.shape[0]:
        raise ValueError("permutation length must match the number of tableau rows.")
    if sorted(pi.tolist()) != list(range(tableau.shape[0])):
        raise ValueError("permutation must be a permutation of the tableau row indices.")

    F = map_tableau_to_target_tableau(tableau, tableau[pi], p=p)
    if F is None:
        raise ValueError("No symplectic map found for this Pauli-row permutation.")

    return F


def _gate_from_pauli_permutation(
    pauli_sum: PauliSum,
    pi: tuple[int, ...],
    *,
    lift_method: str = "inverse",
) -> Gate | None:
    p = int(pauli_sum.lcm)
    if not np.array_equal(pauli_sum.symplectic_product_matrix()[np.ix_(pi, pi)], pauli_sum.symplectic_product_matrix()):
        return None
    if not np.allclose(np.asarray(pauli_sum.weights)[list(pi)], np.asarray(pauli_sum.weights)):
        return None

    try:
        row_action = symplectic_from_pauli_permutation(
            pauli_sum.tableau,
            np.asarray(pi),
            p=p,
        )
    except Exception:
        return None

    gate_symplectic = np.asarray(row_action, dtype=int).T % p
    h0 = get_phase_vector(gate_symplectic, p)
    trial_gate = Gate("Symmetry", gate_symplectic, h0)

    target = pauli_sum[list(pi)]
    trial = trial_gate.act(pauli_sum, tuple(range(pauli_sum.n_qudits())))
    delta = (target.phases - trial.phases) % (2 * int(pauli_sum.lcm))

    h_lin = solve_phase_vector_h_from_residual(
        pauli_sum.tableau,
        delta,
        pauli_sum.dimensions,
    )
    if h_lin is None:
        return None

    gate = Gate("Symmetry", gate_symplectic, (h0 + h_lin) % (2 * int(pauli_sum.lcm)))
    out = gate.act(pauli_sum, tuple(range(pauli_sum.n_qudits()))).to_standard_form()
    ref = pauli_sum.to_standard_form()
    out.weight_to_phase()
    ref.weight_to_phase()
    if not np.array_equal(out.tableau, ref.tableau):
        return None
    if not np.array_equal(out.phases % (2 * int(pauli_sum.lcm)), ref.phases % (2 * int(pauli_sum.lcm))):
        return None
    if not np.allclose(out.weights, ref.weights, atol=1e-8, rtol=0):
        return None
    return gate


def find_igraph_clifford_symmetries(
    pauli_sum: PauliSum,
    *,
    num_symmetries: int | None = 1,
    lift_method: str = "inverse",
    extra_invs: str = "none",
    color_mode: str = "wl",
    max_wl_rounds: int = 10,
    circuit_augmented_graph: bool | str = False,
) -> tuple[list[Gate], int]:
    """Return up to ``num_symmetries`` Clifford symmetries using igraph automorphisms."""
    if num_symmetries is not None and int(num_symmetries) <= 0:
        return [], 0
    if circuit_augmented_graph:
        raise NotImplementedError(
            "circuit_augmented_graph is part of the custom graph automorphism branch, "
            "not the igraph-default PR."
        )

    pauli = pauli_sum.copy()
    pauli.weight_to_phase()
    if not np.all(pauli.dimensions == pauli.dimensions[0]):
        raise ValueError("igraph Clifford symmetry finding currently requires uniform qudit dimensions.")

    S_mod = pauli.symplectic_product_matrix()
    p = int(pauli.lcm)
    col_invariants = _column_invariants(pauli, str(extra_invs))
    base_colors, _ = _build_base_partition(
        S_mod,
        p,
        coeffs=np.asarray(pauli.weights),
        col_invariants=col_invariants if color_mode == "wl" else None,
        max_rounds=max_wl_rounds,
        color_mode=color_mode,
    )

    graph, graph_colors = build_subdivision_graph_from_s_mod(S_mod, base_colors)
    automorphism_group = graph.automorphism_group(color=graph_colors)
    generators = (
        automorphism_group.generators
        if hasattr(automorphism_group, "generators")
        else automorphism_group
    )

    n_paulis = pauli.n_paulis()
    identity_pauli_perm = tuple(range(n_paulis))
    checked = 0
    symmetries: list[Gate] = []
    seen_pauli_perms: set[tuple[int, ...]] = set()
    seen_gates: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()

    def wanted_enough() -> bool:
        return num_symmetries is not None and len(symmetries) >= int(num_symmetries)

    def add_if_clifford(pi_tuple: tuple[int, ...]) -> None:
        nonlocal checked
        if pi_tuple in seen_pauli_perms:
            return
        seen_pauli_perms.add(pi_tuple)
        checked += 1
        gate = _gate_from_pauli_permutation(pauli, pi_tuple, lift_method=lift_method)
        if gate is None:
            return
        key = (
            tuple(np.asarray(gate.symplectic, dtype=int).reshape(-1).tolist()),
            tuple(np.asarray(gate.phase_vector(), dtype=int).reshape(-1).tolist()),
        )
        if key in seen_gates:
            return
        seen_gates.add(key)
        symmetries.append(gate)

    for generator in generators:
        pi_tuple = _extract_pauli_permutation(_permutation_to_tuple(generator), n_paulis)
        if pi_tuple is None or pi_tuple == identity_pauli_perm:
            continue
        add_if_clifford(pi_tuple)
        if wanted_enough():
            return symmetries, checked

    for full_perm in _generated_group_permutations(generators, graph.vcount()):
        pi_tuple = _extract_pauli_permutation(full_perm, n_paulis)
        if pi_tuple is None or pi_tuple == identity_pauli_perm:
            continue
        add_if_clifford(pi_tuple)
        if wanted_enough():
            return symmetries, checked

    return symmetries, checked


__all__ = ["find_igraph_clifford_symmetries"]
