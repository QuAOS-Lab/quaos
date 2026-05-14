from __future__ import annotations

from typing import Any

import numpy as np

from sympleq.core.graphs.graph_automorphism_leaf import check_leaf
from sympleq.core.graphs.graph_automorphism_search import prepare_clifford_ga_search
from sympleq.core.graphs.graph_builder import build_subdivision_graph_from_s_mod


def _permutation_to_tuple(permutation: Any) -> tuple[int, ...]:
    if hasattr(permutation, "mapping"):
        return tuple(permutation.mapping)
    return tuple(permutation)


def _compose_permutations(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[i] for i in right)



def _extract_pauli_permutation(full_perm: tuple[int, ...], n_paulis: int) -> tuple[int, ...] | None:
    """Return the induced Pauli permutation, or None if Pauli vertices do not map to Pauli vertices.

    With correct vertex colours this should never fail, but the check protects
    against colour collisions or graph-builder changes.
    """
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
    queue_items = [(identity, ())]

    while queue_items:
        current, word = queue_items.pop(0)
        for idx, generator in enumerate(generators):
            for side, candidate in (
                ("L", _compose_permutations(generator, current)),
                ("R", _compose_permutations(current, generator)),
            ):
                if candidate in seen:
                    continue
                seen.add(candidate)
                candidate_word = word + ((idx, side),)
                queue_items.append((candidate, candidate_word))
                yield candidate, candidate_word


def find_igraph_clifford_symmetries(
    pauli_sum: Any,
    *,
    extra_invs: str = "none",
    circuit_augmented_graph: bool | str = False,
) -> tuple[list[Any], int]:
    prepared = prepare_clifford_ga_search(
        pauli_sum,
        extra_column_invariants=str(extra_invs),
        p2_bitset="auto",
        color_mode="wl",
        max_wl_rounds=0,
        circuit_augmented_graph=circuit_augmented_graph,
    )
    # Use the explicit augmented graph view when requested.  In "wl" mode
    # this remains the Pauli-only graph with WL-refined Pauli colours.
    graph_S_mod = getattr(prepared, "graph_S_mod", prepared.S_mod)
    graph_base_colors = getattr(prepared, "graph_base_colors", prepared.base_colors)
    graph, graph_colors = build_subdivision_graph_from_s_mod(graph_S_mod, graph_base_colors)
    automorphism_group = graph.automorphism_group(color=graph_colors)
    generators = (
        automorphism_group.generators
        if hasattr(automorphism_group, "generators")
        else automorphism_group
    )

    n_paulis = int(getattr(prepared, "n_pauli_vertices", pauli_sum.n_paulis()))
    identity_pauli_perm = tuple(range(n_paulis))
    checked = 0
    for generator in generators:
        full_perm = _permutation_to_tuple(generator)
        pi_tuple = _extract_pauli_permutation(full_perm, n_paulis)
        if pi_tuple is None or pi_tuple == identity_pauli_perm:
            continue
        checked += 1
        gate = check_leaf(np.asarray(pi_tuple, dtype=np.int64), prepared.leaf_ctx)
        if gate is not None:
            return [gate], checked

    for full_perm, _word in _generated_group_permutations(generators, graph.vcount()):
        pi_tuple = _extract_pauli_permutation(full_perm, n_paulis)
        if pi_tuple is None or pi_tuple == identity_pauli_perm:
            continue
        checked += 1
        gate = check_leaf(np.asarray(pi_tuple, dtype=np.int64), prepared.leaf_ctx)
        if gate is not None:
            return [gate], checked

    return [], checked


__all__ = ["find_igraph_clifford_symmetries"]
