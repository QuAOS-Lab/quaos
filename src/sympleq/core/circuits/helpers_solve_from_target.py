"""
Helpers for mapping a Paulisum to a target Paulisum.

The helper functions completes the basis for the input and target tableaus
where they do not already contain 2n independent rows.
The row-action Clifford map is obtained from
``input_basis^{-1} @ output_basis``.
"""
# TODO Mixed-dimension
# TODO Check why qubits are not working
# TODO Permutation

import numpy as np
from sympleq.core.paulis import PauliSum
from sympleq.core.paulis._typing import TableauType, PhasesType
from sympleq.core.finite_field_solvers import (
    get_linear_dependencies,
    gf_inv,
    solve_linear_system_over_gf,
)
from sympleq.core.circuits.utils import symplectic_form, symplectic_product_matrix
from sympleq.core.symmetries.modular_helpers import nullspace_mod
from sympleq.core.phase_correction.solve import solve_phase_vector_h_from_residual
from sympleq.core.circuits import Gate


def check_mappable_via_clifford(paulisum_tableau: TableauType,
                                target_paulisum_tableau: TableauType,
                                p: int = 2) -> bool:
    sym_check = np.all(
        symplectic_product_matrix(paulisum_tableau, p) == symplectic_product_matrix(target_paulisum_tableau, p)
    )
    if sym_check:
        return True

    return False


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
    V is the basis for the target tableau.``solve_linear_system_over_gf``
    returns one particular solution y0. If y0 is dependent on the existing
    target rows, k (nullspace direction of ``A``) is added to y0, and y=y0+k
    is checked to see if it is independent, and returned.
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
    Complete the partial basis of input and output tableaus to 2n independent entries.
    Selects independent rows ``U` -> V` from the supplied tableaus,
    Adds standard input basis vectors ``e``. For each accepted ``e``, it
    chooses a target row ``y`` satisfying ``<V_i, y> = <U_i, e>`` for all existing rows,
    so the extended map preserves the symplectic product matrix.
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


def map_tableau_to_target_tableau(
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


# To be removed??
# Could be useful for mixed??
def get_phase_vector(gate_symplectic: TableauType, dimension: int) -> PhasesType:
    """
    Calculate the phase vector for a gate given its symplectic matrix.

    See PRA 71, 042315 (2005) Eq. (10).
    Solves for h

    Args:
        gate_symplectic (np.ndarray): The symplectic matrix of the gate.
        dimension (int): The dimension of the qudit.

    Returns:
        np.ndarray: The phase vector of the gate.
    """
    n_qudits = gate_symplectic.shape[0] // 2

    U = np.zeros((2 * n_qudits, 2 * n_qudits), dtype=int)
    U[n_qudits:, :n_qudits] = np.eye(n_qudits, dtype=int)
    lhs = (dimension - 1) * np.diag(gate_symplectic.T @ U @ gate_symplectic) % 2  # Eq. (10) mod 2 is there for all d
    return lhs

# (- lhs) % (2 * dimension) TODO: d > 2 testing, do we need the minus?
# Shreya: For the version implemented here, we do not. Ideally we should solve for this equation:
# lhs + (dimension - 1) * np.diag(gate_symplectic.T @ U @ gate_symplectic) = 0 % 2,
# which can have more than one solution, implying more than one starting point for the phase-solver.
# Also, the implemented version returns a vector with all zero for qudits > 2, as (dimension -1) is even there.
# For qubits, it gives an initial parity vector, which matters. making it a qubit only function for now. However,
# could be useful for mixed??


"""
One more function is needed for a clean implementation of solve_from_target:
'from sympleq.core.phase_correction.solve import solve_phase_vector_h_from_residual'
However, this needs to be cleaned up, especially the 'solve linear equations parts' in it.
I am importing it for now before deciding if we are cleaning up
phase correction as a whole, or just for the 'solve from target'.

I am also defining the 'solve_from_target' function below; to be incorporated as the classmethod in Gates.py later.
Still not does do mixed qudits.
Odd primes work; qubits do not find a gate, fails at phase correction.
"""


def solve_from_target(#cls,
        input_pauli_sum: PauliSum,
        target_pauli_sum: PauliSum
) -> Gate:
    """
    Find a Clifford gate that maps an input PauliSum to the target PauliSum.

    Uses symplectic transvections to find a symplectic matrix F such that
    input_tableau @ F = target_tableau (mod p), with p=`dimension`.

    Parameters
    ----------
    input_pau : TableauLike
        Input Pauli tableau of shape (m, 2n) where m is the number of Paulis
        and n is the number of qudits.
    target_tableau : TableauLike
        Target Pauli tableau of the same shape.
    dimension : int
        Local Hilbert space dimension (e.g., 2 for qubits).

    Returns
    -------
    Gate
        A Clifford gate whose symplectic matrix performs the mapping.

    Raises
    ------
    ValueError
        If the tableaus have different shapes or are not mappable via Clifford,
        or does not have a phase correction.

    Notes
    -----
    The input and target must have matching symplectic
    product matrices for a Clifford mapping to exist.
    """

    input_tableau = input_pauli_sum.tableau
    target_tableau = target_pauli_sum.tableau

    if input_tableau.shape != target_tableau.shape:
        raise ValueError(
            f"Tableau shapes must match: {input_tableau.shape} vs {target_tableau.shape}"
        )

    if input_tableau.ndim == 1:
        input_tableau = input_tableau.reshape(1, -1)
        target_tableau = target_tableau.reshape(1, -1)

    n_qudits = input_tableau.shape[1] // 2
    p = int(input_pauli_sum.lcm)

    if check_mappable_via_clifford(input_tableau, target_tableau, p) == False:
        raise ValueError(
            f"Not mappable via Clifford: {input_tableau} ->  {target_tableau}."
        )
    else:
        F_total = map_tableau_to_target_tableau(input_tableau, target_tableau, p)

        if p == 2:
            h0 = get_phase_vector(F_total, p)
        else:
            h0 = np.asarray([0] * 2 * n_qudits, dtype=int)

        trial_gate = Gate("trial", F_total.T, h0)
        trial_pauli_sum = trial_gate.act(input_pauli_sum, tuple(range(n_qudits)))

        delta_2L = (target_pauli_sum.phases - trial_pauli_sum.phases) % (2 * p)
        delta_2L = np.asarray(delta_2L, dtype=int)

        if p == 2:
            h_lin = solve_phase_vector_h_from_residual(input_tableau, delta_2L, [p] * n_qudits, debug=True)
        else:
            h_lin = solve_phase_vector_h_from_residual(input_tableau, delta_2L, [p] * n_qudits)

        if h_lin is None:
            return None

        else:
            h_final = (h0 + h_lin) % (2 * p)

            final_gate = Gate("final", F_total.T, h_final)

        return final_gate

