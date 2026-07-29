"""
Helpers for mapping a Paulisum to a target Paulisum.

The helper functions completes the basis for the input and target tableaus
where they do not already contain 2n independent rows. The row-action Clifford map is obtained from
``input_basis^{-1} @ output_basis``.
"""
# TODO Mixed-dimension
# TODO Permutation

import numpy as np
from sympleq.core.paulis._typing import TableauType
from sympleq.core.finite_field_solvers import (
    get_linear_dependencies,
    gf_inv,
    solve_linear_system_over_gf,
)
from sympleq.core.circuits.utils import symplectic_form
from sympleq.core.symmetries.modular_helpers import nullspace_mod


def gf_rank(A: TableauType, p: int) -> int:
    """
    Return the rank of ``A`` over GF(p).
    """
    basis_indices, _ = get_linear_dependencies(
        np.asarray(A, dtype=int) % p,
        p,
        compute_dependencies=False,
    )
    return len(basis_indices)


def independent_solution(A: TableauType, b: TableauType, V: TableauType, p: int) -> TableauType:
    """
    Solve ``A @ y = b`` over GF(p), choosing ``y`` outside ``span(V)``.
    ``solve_linear_system_over_gf`` returns one particular solution. If that
    solution is dependent on the existing target rows, a nullspace direction of
    ``A`` is added without changing the equation ``A @ y = b``.
    """
    A = np.asarray(A, dtype=int) % p
    b = np.asarray(b, dtype=int) % p
    V = np.asarray(V, dtype=int) % p

    y0 = np.asarray(solve_linear_system_over_gf(A, b, p), dtype=int) % p

    if gf_rank(np.vstack([V, y0]), p) == gf_rank(V, p) + 1:
        return y0

    K = np.asarray(nullspace_mod(A, p), dtype=int) % p

    for k in K.T:
        y = (y0 + k) % p

        if gf_rank(np.vstack([V, y]), p) == gf_rank(V, p) + 1:
            return y

    raise ValueError("No independent compatible solution found.")


def complete_basis(input_tab: TableauType, target_tab: TableauType, p: int) -> tuple[TableauType, TableauType]:
    """
    Complete a compatible partial tableau map to full input/output bases.
    Starting from independent rows ``U -> V`` selected from the supplied tableaus,
    this adds standard input basis vectors ``e``. For each accepted ``e``, it
    chooses a target row ``y`` satisfying
    ``<V_i, y> = <U_i, e>`` for all existing rows, so the extended map preserves
    the symplectic product matrix.
    """
    input_tab = np.asarray(input_tab, dtype=int) % p
    target_tab = np.asarray(target_tab, dtype=int) % p

    d = input_tab.shape[1]
    Omega = symplectic_form(d // 2, p)

    basis_indices, _ = get_linear_dependencies(
        input_tab,
        p,
        compute_dependencies=False,
    )

    U = input_tab[basis_indices]
    V = target_tab[basis_indices]

    for e in np.eye(d, dtype=int):
        if gf_rank(np.vstack([U, e]), p) != gf_rank(U, p) + 1:
            continue

        A = (V @ Omega) % p
        b = (U @ Omega @ e) % p
        try:
            y = independent_solution(A, b, V, p)
        except ValueError:
            continue

        U = np.vstack([U, e % p])
        V = np.vstack([V, y])

        if gf_rank(U, p) == d and gf_rank(V, p) == d:
            return U, V

    raise ValueError("Could not complete basis.")


def map_paulisum_to_target_paulisum(
    paulisum_tableau: TableauType, target_paulisum_tableau: TableauType, p: int = 2
) -> TableauType | None:
    """
    Return a row-action symplectic map taking one tableau to another.
    If the input tableau already contains a complete independent basis, the map
    is computed directly from that basis. Otherwise the partial map is completed
    first via :func:`complete_basis`. The returned matrix ``F`` satisfies
    ``paulisum_tableau @ F == target_paulisum_tableau`` modulo ``p``.
    """
    input_tab = np.asarray(paulisum_tableau, dtype=int) % p
    output_tab = np.asarray(target_paulisum_tableau, dtype=int) % p

    if input_tab.ndim != 2 or output_tab.ndim != 2 or input_tab.shape != output_tab.shape:
        raise ValueError("Input and target tableaus must be 2-dimensional arrays with matching shape.")

    n_cols = input_tab.shape[1]
    if n_cols % 2 != 0:
        raise ValueError("Pauli tableau width must be even.")

    basis_indices, _ = get_linear_dependencies(input_tab, p, compute_dependencies=False)

    if len(basis_indices) < n_cols:
        input_basis, output_basis = complete_basis(input_tab, output_tab, p)

    else:
        input_basis = input_tab[basis_indices]
        output_basis = output_tab[basis_indices]

    try:
        F = (gf_inv(input_basis, p) @ output_basis) % p
    except ValueError:
        return None

    if not np.array_equal((input_tab @ F) % p, output_tab):
        raise ValueError("Complete-basis map failed to reproduce the target tableau.")

    return F
