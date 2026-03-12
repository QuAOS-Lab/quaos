"""
This module uses arXiv:1803.06987 to find symplectic solutions.

So far it only works for GF(2), as in the original paper. It could be extended to GF(p).
"""
from __future__ import annotations

import numpy as np
from sympleq.core.finite_field_solvers import solve_gf2
from sympleq.core.circuits.utils import transvection_matrix, symplectic_product_arrays, symplectic_product_matrix
from sympleq._typing import IntNDArray
from sympleq.core.paulis._typing import TableauType


def find_symplectic_solution(u: IntNDArray, v: IntNDArray) -> IntNDArray:
    """
    Find a binary vector w such that <u,w> = <v,w> = 1 in symplectic inner product.

    Args:
        u: Binary vector of length 2n.
        v: Binary vector of length 2n.

    Returns:
        Binary vector w of length 2n.

    Raises:
        ValueError: If u or v is the zero vector.
        Exception: If no solution exists in GF(2).
    """
    n = len(u) // 2

    # Check for impossible cases first
    if np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n)):
        # Cannot have symplectic inner product 1 with zero vector
        raise ValueError("Cannot find solution with zero vector input")

    # Check if u and v are symplectically independent
    symplectic_product_arrays_uv = symplectic_product_arrays(u, v)

    if symplectic_product_arrays_uv == 1:
        # u and v are symplectically independent - use direct construction
        return direct_construction(u, v)
    else:
        # u and v are symplectically orthogonal - use general linear system solver
        return solve_general_system(u, v)


def direct_construction(u: IntNDArray, v: IntNDArray) -> IntNDArray:
    """
    Direct construction when u and v are symplectically independent (<u,v> = 1).

    Uses a geometric construction based on symplectic orthogonality rather than
    solving a linear system.

    Args:
        u: Binary vector of length 2n.
        v: Binary vector of length 2n, symplectically independent from u.

    Returns:
        Binary vector w of length 2n such that <u,w> = <v,w> = 1.
    """
    n = len(u) // 2

    # For symplectically independent u and v, we can construct w geometrically
    # The key insight: if <u,v> = 1, then u and v span a 2D symplectic subspace

    # Strategy: construct w as a linear combination w = a*u + b*v + orthogonal_part
    # where orthogonal_part is symplectically orthogonal to both u and v

    # First, try the simplest approach: w = u + v
    w_candidate = (u + v) % 2
    if (symplectic_product_arrays(u, w_candidate) == 1 and symplectic_product_arrays(v, w_candidate) == 1):
        return w_candidate

    # If that doesn't work, try other simple combinations
    for a in [0, 1]:
        for b in [0, 1]:
            if a == 0 and b == 0:
                continue
            w_candidate = (a * u + b * v) % 2
            if (symplectic_product_arrays(u, w_candidate) == 1 and symplectic_product_arrays(v, w_candidate) == 1):
                return w_candidate

    # If simple combinations don't work, we need to add an orthogonal component
    # Find a vector orthogonal to both u and v, then add it to a base combination

    # Start with a base that might work
    w_base = u.copy()  # or v, or u+v

    # Try adding standard basis vectors to correct the inner products
    for i in range(2 * n):
        w_candidate = w_base.copy()
        w_candidate[i] = (w_candidate[i] + 1) % 2

        if (symplectic_product_arrays(u, w_candidate) == 1 and symplectic_product_arrays(v, w_candidate) == 1):
            return w_candidate

    # Fallback to the general linear solver if geometric construction fails
    A = np.zeros((2, 2 * n), dtype=int)
    A[0, :n] = u[n:]
    A[0, n:] = u[:n]
    A[1, :n] = v[n:]
    A[1, n:] = v[:n]
    b = np.array([1, 1])

    solution = solve_gf2(A, b)

    if solution is None:
        raise Exception("Could not find a solution in gf2.")

    return solution


def solve_general_system(u: IntNDArray, v: IntNDArray) -> IntNDArray:
    """
    Solve the general case when u and v are symplectically orthogonal (<u,v> = 0).

    Args:
        u: Binary vector of length 2n.
        v: Binary vector of length 2n, symplectically orthogonal to u.

    Returns:
        Binary vector w of length 2n such that <u,w> = <v,w> = 1.

    Raises:
        Exception: If no solution exists in GF(2).
    """
    n = len(u) // 2

    # Set up the linear system A @ w = b
    A = np.zeros((2, 2 * n), dtype=int)
    A[0, :n] = u[n:]  # X part of u multiplies Z part of w
    A[0, n:] = u[:n]  # Z part of u multiplies X part of w
    A[1, :n] = v[n:]  # X part of v multiplies Z part of w
    A[1, n:] = v[:n]  # Z part of v multiplies X part of w

    b = np.array([1, 1])

    solution = solve_gf2(A, b)

    if solution is None:
        raise Exception("Could not find a solution in gf2.")

    return solution


def find_symplectic_solution_extended(u: IntNDArray, v: IntNDArray,
                                      t_vectors: list | None = None) -> IntNDArray:
    """
    Find a binary vector w such that:
    - <u,w> = 1
    - <v,w> = 1
    - <t_i,w> = <t_i,v> for all t_i in t_vectors

    Args:
        u: Binary vector of length 2n.
        v: Binary vector of length 2n.
        t_vectors: List of binary vectors of length 2n for additional constraints.

    Returns:
        Binary vector w of length 2n.

    Raises:
        ValueError: If u or v is the zero vector.
        Exception: If no solution exists in GF(2).
    """
    if t_vectors is None or len(t_vectors) == 0:
        return find_symplectic_solution(u, v)

    n = len(u) // 2

    # Check for impossible cases first
    if np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n)):
        raise ValueError("Cannot find solution with zero vector input")

    # For extended system, we always use the general linear solver
    # since the additional constraints break the geometric structure
    return solve_extended_system(u, v, t_vectors)


def solve_extended_system(u: IntNDArray, v: IntNDArray, t_vectors: list) -> IntNDArray:
    """
    Solve the extended system with additional t_i constraints.

    Args:
        u: Binary vector of length 2n (primary constraint).
        v: Binary vector of length 2n (primary constraint).
        t_vectors: List of binary vectors of length 2n (additional constraints).

    Returns:
        Binary vector w of length 2n.

    Raises:
        Exception: If no solution exists in GF(2).
    """
    n = len(u) // 2
    k = len(t_vectors)

    # Set up the linear system A @ w = b
    # We have 2 + k constraints total
    A = np.zeros((2 + k, 2 * n), dtype=int)
    b = np.zeros(2 + k, dtype=int)

    # First constraint: <u, w> = 1
    A[0, :n] = u[n:]  # X part of u multiplies Z part of w
    A[0, n:] = u[:n]  # Z part of u multiplies X part of w
    b[0] = 1

    # Second constraint: <v, w> = 1
    A[1, :n] = v[n:]  # X part of v multiplies Z part of w
    A[1, n:] = v[:n]  # Z part of v multiplies X part of w
    b[1] = 1

    # Additional constraints: <t_i, w> = <t_i, v>
    for i, t in enumerate(t_vectors):
        row_idx = 2 + i
        A[row_idx, :n] = t[n:]  # X part of t_i multiplies Z part of w
        A[row_idx, n:] = t[:n]  # Z part of t_i multiplies X part of w
        b[row_idx] = symplectic_product_arrays(t, v)

    solution = solve_gf2(A, b)

    if solution is None:
        raise Exception("Could not find a solution in gf2.")

    return solution


def check_mappable_via_clifford(pauli_sum_tableau: TableauType,
                                target_pauli_sum_tableau: TableauType,
                                p: int = 2) -> bool:
    sym_check = np.all(
        symplectic_product_matrix(pauli_sum_tableau, p) == symplectic_product_matrix(target_pauli_sum_tableau, p)
    )
    if sym_check:
        return True

    return False


def map_single_pauli_string_to_target(pauli_string_tableau: TableauType, target_pauli_string_tableau: TableauType,
                                      constraint_paulis: list | None = None) -> TableauType:
    sp = symplectic_product_arrays(pauli_string_tableau, target_pauli_string_tableau)
    if sp == 1:
        h = pauli_string_tableau + target_pauli_string_tableau

        F_h = transvection_matrix(h)

        return F_h

    if sp == 0:
        w = find_symplectic_solution_extended(pauli_string_tableau, target_pauli_string_tableau, constraint_paulis)
        h_1 = target_pauli_string_tableau + w
        h_2 = pauli_string_tableau + w

        F_h_1 = transvection_matrix(h_1)
        F_h_2 = transvection_matrix(h_2)

        return (F_h_1 @ F_h_2) % 2

    else:
        raise Exception(f'sp = {sp}...This should never happen')


def map_pauli_sum_to_target_tableau(
        pauli_sum_tableau: TableauType, target_pauli_sum_tableau: TableauType) -> TableauType:
    """
    Map a Pauli sum to a target Pauli sum using symplectic transvections.
    """
    if not check_mappable_via_clifford(pauli_sum_tableau, target_pauli_sum_tableau):
        raise Exception(f'SPM not equal. Cannot map\n{pauli_sum_tableau} to\n{target_pauli_sum_tableau}')

    m = len(pauli_sum_tableau)
    n = len(pauli_sum_tableau[0]) // 2
    mapped_paulis = []
    F = np.eye(2 * n, dtype=int)
    for i in range(m):
        # update the starting point to whatever previous solutions mapped it to
        ps = (pauli_sum_tableau[i] @ F) % 2
        target_ps = target_pauli_sum_tableau[i]

        if np.array_equal(ps, target_ps):
            mapped_paulis.append(target_ps)  # these are now the constraints for the next iteration
            continue

        F_map = map_single_pauli_string_to_target(ps, target_ps, mapped_paulis)
        assert np.all((ps @ F_map) % 2 == target_ps), f"\n{F_map}\n{ps}\n{(ps @ F_map) % 2}\n{target_ps}"
        for mp in mapped_paulis:
            assert np.all((mp @ F_map) % 2 == mp), f"\n{F_map}\n{mp}\n{(mp @ F_map) % 2}"
        mapped_paulis.append(target_ps)  # these are now the constraints for the next iteration
        F = (F @ F_map) % 2

    return F
