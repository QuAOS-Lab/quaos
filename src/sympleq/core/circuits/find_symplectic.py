"""
This module uses arXiv:1803.06987 to find symplectic solutions.

So far it only works for GF(2), as in the original paper. It could be extended to GF(p).
"""
# TODO: move all functions that are required for the class method
#       `input_to_target` to a unique file, that works with mixed dimension.
# TODO: It seems that the logic can be simplified quite a bit here...
#       I would use S = input^-1 target (sorry if typos) plus some check before that ensure consistency...
#       This should be applicable to qudits as well and should be quite efficient?
#
#       Then it is only a matter to fix the phase, which should be always
#       (? please correct me if wrong) doable efficiently...
from __future__ import annotations

import numpy as np
import galois
from sympleq._typing import IntNDArray
from sympleq.core.circuits.utils import transvection_matrix, symplectic_product_arrays, symplectic_product_matrix
from sympleq.core.finite_field_solvers import get_linear_dependencies, gf2_inv, solve_gf2, solve_linear_system_over_gf
from sympleq.core.paulis._typing import TableauType


def _symplectic_constraint_row(vec: IntNDArray, p: int = 2) -> IntNDArray:
    """
    Build a row r such that r @ w = <vec, w> mod p for [x|z]-ordered vectors.
    """
    n = len(vec) // 2
    row = np.zeros(2 * n, dtype=int)
    row[:n] = (-vec[n:]) % p
    row[n:] = vec[:n] % p
    return row


def find_symplectic_solution(u: IntNDArray, v: IntNDArray, p: int = 2) -> IntNDArray:
    """
    Find w such that <u,w> = <v,w> = 1 over GF(p).

    Args:
        u: Vector of length 2n.
        v: Vector of length 2n.
        p: Prime dimension.

    Returns:
        Vector w of length 2n.

    Raises:
        ValueError: If u or v is the zero vector.
        Exception: If no solution exists in GF(p).
    """
    n = len(u) // 2

    # Check for impossible cases first
    if np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n)):
        # Cannot have symplectic inner product 1 with zero vector
        raise ValueError("Cannot find solution with zero vector input")

    # Check if u and v are symplectically independent
    symplectic_product_arrays_uv = symplectic_product_arrays(u, v, p)

    if p == 2 and symplectic_product_arrays_uv == 1:
        # Keep the fast geometric shortcut for the binary case.
        return direct_construction(u, v)

    # General case over GF(p): solve linear constraints directly.
    return solve_general_system(u, v, p=p)


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
    if symplectic_product_arrays(u, w_candidate) == 1 and symplectic_product_arrays(v, w_candidate) == 1:
        return w_candidate

    # If that doesn't work, try other simple combinations
    for a in [0, 1]:
        for b in [0, 1]:
            if a == 0 and b == 0:
                continue
            w_candidate = (a * u + b * v) % 2
            if symplectic_product_arrays(u, w_candidate) == 1 and symplectic_product_arrays(v, w_candidate) == 1:
                return w_candidate

    # If simple combinations don't work, we need to add an orthogonal component
    # Find a vector orthogonal to both u and v, then add it to a base combination

    # Start with a base that might work
    w_base = u.copy()  # or v, or u+v

    # Try adding standard basis vectors to correct the inner products
    for i in range(2 * n):
        w_candidate = w_base.copy()
        w_candidate[i] = (w_candidate[i] + 1) % 2

        if symplectic_product_arrays(u, w_candidate) == 1 and symplectic_product_arrays(v, w_candidate) == 1:
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

    return np.asarray(solution, dtype=int) % 2


def solve_general_system(u: IntNDArray, v: IntNDArray, p: int = 2) -> IntNDArray:
    """
    Solve linear constraints for <u,w> = 1 and <v,w> = 1 over GF(p).

    Args:
        u: Vector of length 2n.
        v: Vector of length 2n.
        p: Prime dimension.

    Returns:
        Vector w of length 2n such that <u,w> = <v,w> = 1.

    Raises:
        Exception: If no solution exists in GF(p).
    """
    n = len(u) // 2  # kept for shape clarity

    # Set up the linear system A @ w = b for:
    # <u,w> = 1, <v,w> = 1 (mod p)
    A = np.zeros((2, 2 * n), dtype=int)
    A[0] = _symplectic_constraint_row(u, p)
    A[1] = _symplectic_constraint_row(v, p)
    b = np.array([1, 1], dtype=int) % p

    if p == 2:
        solution = solve_gf2(A, b)
    else:
        try:
            solution = solve_linear_system_over_gf(A, b, p)
        except ValueError:
            solution = None

    if solution is None:
        raise Exception(f"Could not find a solution in GF({p}).")

    return np.asarray(solution, dtype=int) % p


def find_symplectic_solution_extended(
    u: IntNDArray, v: IntNDArray, t_vectors: list[IntNDArray] | None = None, p: int = 2
) -> IntNDArray:
    """
    Find w such that:
    - <u,w> = 1
    - <v,w> = 1
    - <t_i,w> = <t_i,v> for all t_i in t_vectors

    Args:
        u: Vector of length 2n.
        v: Vector of length 2n.
        t_vectors: Additional constraint vectors.
        p: Prime dimension.

    Returns:
        Vector w of length 2n.

    Raises:
        ValueError: If u or v is the zero vector.
        Exception: If no solution exists in GF(p).
    """
    if t_vectors is None or len(t_vectors) == 0:
        return find_symplectic_solution(u, v, p=p)

    n = len(u) // 2

    # Check for impossible cases first
    if np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n)):
        raise ValueError("Cannot find solution with zero vector input")

    # For extended system, we always use the general linear solver
    # since the additional constraints break the geometric structure
    return solve_extended_system(u, v, t_vectors, p=p)


def solve_extended_system(u: IntNDArray, v: IntNDArray, t_vectors: list[IntNDArray], p: int = 2) -> IntNDArray:
    """
    Solve the extended system with additional t_i constraints.

    Args:
        u: Vector of length 2n (primary constraint).
        v: Vector of length 2n (primary constraint).
        t_vectors: Additional vectors (additional constraints).
        p: Prime dimension.

    Returns:
        Vector w of length 2n.

    Raises:
        Exception: If no solution exists in GF(p).
    """
    n = len(u) // 2
    k = len(t_vectors)

    # Set up the linear system A @ w = b
    # We have 2 + k constraints total
    A = np.zeros((2 + k, 2 * n), dtype=int)
    b = np.zeros(2 + k, dtype=int)

    # First constraint: <u, w> = 1
    A[0] = _symplectic_constraint_row(u, p)
    b[0] = 1

    # Second constraint: <v, w> = 1
    A[1] = _symplectic_constraint_row(v, p)
    b[1] = 1

    # Additional constraints: <t_i, w> = <t_i, v>
    for i, t in enumerate(t_vectors):
        row_idx = 2 + i
        A[row_idx] = _symplectic_constraint_row(t, p)
        b[row_idx] = symplectic_product_arrays(t, v, p)

    A %= p
    b %= p
    if p == 2:
        solution = solve_gf2(A, b)
    else:
        try:
            solution = solve_linear_system_over_gf(A, b, p)
        except ValueError:
            solution = None

    if solution is None:
        raise Exception(f"Could not find a solution in GF({p}).")

    return np.asarray(solution, dtype=int) % p


def check_mappable_via_clifford(pauli_sum_tableau: TableauType,
                                target_pauli_sum_tableau: TableauType,
                                p: int = 2) -> bool:
    sym_check = np.all(
        symplectic_product_matrix(pauli_sum_tableau, p) == symplectic_product_matrix(target_pauli_sum_tableau, p)
    )
    if sym_check:
        return True

    return False


def map_single_pauli_string_to_target(
    pauli_string_tableau: TableauType,
    target_pauli_string_tableau: TableauType,
    constraint_paulis: list[TableauType] | None = None,
    p: int = 2,
) -> TableauType:
    if p != 2:
        from sympleq.core.circuits.find_symplectic_qudits import find_transvection_map_solve_extended

        constraints = [] if constraint_paulis is None else list(constraint_paulis)
        sps = [symplectic_product_arrays(t, target_pauli_string_tableau, p) for t in constraints]
        return find_transvection_map_solve_extended(
            pauli_string_tableau,
            target_pauli_string_tableau,
            constraints=constraints,
            sps=sps,
            p=p,
        )

    sp = symplectic_product_arrays(pauli_string_tableau, target_pauli_string_tableau, p)
    if sp == 1:
        h = pauli_string_tableau + target_pauli_string_tableau

        F_h = transvection_matrix(h, p)

        return F_h

    if sp == 0:
        w = find_symplectic_solution_extended(
            pauli_string_tableau,
            target_pauli_string_tableau,
            constraint_paulis,
            p=p,
        )
        h_1 = target_pauli_string_tableau + w
        h_2 = pauli_string_tableau + w

        F_h_1 = transvection_matrix(h_1, p)
        F_h_2 = transvection_matrix(h_2, p)

        return (F_h_1 @ F_h_2) % p

    else:
        raise Exception(f'sp = {sp}...This should never happen')


def map_pauli_sum_to_target_tableau(
    pauli_sum_tableau: TableauType,
    target_pauli_sum_tableau: TableauType,
    p: int = 2,
    method: str = "inverse",
) -> TableauType:
    """
    Map a Pauli sum tableau to a target tableau.

    ``method="inverse"`` is the default and computes the unique row-action map
    from a complete independent Pauli basis. ``method="transvection"`` uses the
    constructive transvection path. ``method="auto"`` keeps the old compatibility
    behavior: try the inverse method first, then fall back to transvections.
    """
    if not check_mappable_via_clifford(pauli_sum_tableau, target_pauli_sum_tableau, p=p):
        raise Exception(f'SPM not equal. Cannot map\n{pauli_sum_tableau} to\n{target_pauli_sum_tableau}')

    method_key = str(method).lower()
    if method_key not in {"inverse", "transvection", "auto"}:
        raise ValueError("method must be one of 'inverse', 'transvection', or 'auto'.")

    if method_key in {"inverse", "auto"}:
        complete_basis_map = _map_complete_basis_to_target(pauli_sum_tableau, target_pauli_sum_tableau, p=p)
        if complete_basis_map is not None:
            return complete_basis_map
        if method_key == "inverse":
            raise ValueError(
                "Inverse tableau mapping requires the input tableau to contain "
                "a complete independent Pauli basis."
            )

    return _map_pauli_sum_to_target_tableau_by_transvections(
        pauli_sum_tableau,
        target_pauli_sum_tableau,
        p=p,
    )


def _map_pauli_sum_to_target_tableau_by_transvections(
    pauli_sum_tableau: TableauType, target_pauli_sum_tableau: TableauType, p: int = 2
) -> TableauType:
    """Map a Pauli tableau to a target tableau using transvections."""
    if p != 2:
        from sympleq.core.circuits.find_symplectic_qudits import map_paulisum_to_target_paulisum

        return map_paulisum_to_target_paulisum(pauli_sum_tableau, target_pauli_sum_tableau, p)

    m = len(pauli_sum_tableau)
    n = len(pauli_sum_tableau[0]) // 2
    mapped_paulis: list[TableauType] = []
    F = np.eye(2 * n, dtype=int)
    for i in range(m):
        # update the starting point to whatever previous solutions mapped it to
        ps = (pauli_sum_tableau[i] @ F) % p
        target_ps = target_pauli_sum_tableau[i]

        if np.array_equal(ps, target_ps):
            mapped_paulis.append(target_ps)  # these are now the constraints for the next iteration
            continue

        F_map = map_single_pauli_string_to_target(ps, target_ps, mapped_paulis, p=p)
        assert np.all((ps @ F_map) % p == target_ps), f"\n{F_map}\n{ps}\n{(ps @ F_map) % p}\n{target_ps}"
        for mp in mapped_paulis:
            assert np.all((mp @ F_map) % p == mp), f"\n{F_map}\n{mp}\n{(mp @ F_map) % p}"
        mapped_paulis.append(target_ps)  # these are now the constraints for the next iteration
        F = (F @ F_map) % p

    return F


def symplectic_from_pauli_permutation(
    pauli_sum_tableau: TableauType,
    permutation: TableauType,
    p: int = 2,
    method: str = "inverse",
) -> TableauType:
    """
    Return the row-action symplectic map induced by a Pauli-row permutation.

    The returned matrix ``F`` satisfies ``tableau @ F == tableau[permutation]``
    modulo ``p`` when the permutation preserves the symplectic product matrix.
    This is graph-independent core logic used by higher-level symmetry finders.
    The default method is the inverse-basis lift.
    """
    tableau = np.asarray(pauli_sum_tableau, dtype=int) % p
    pi = np.asarray(permutation, dtype=np.int64).reshape(-1)
    if tableau.ndim != 2:
        raise ValueError("pauli_sum_tableau must be a 2D tableau.")
    if pi.shape[0] != tableau.shape[0]:
        raise ValueError("permutation length must match the number of tableau rows.")
    if sorted(pi.tolist()) != list(range(tableau.shape[0])):
        raise ValueError("permutation must be a permutation of the tableau row indices.")
    return map_pauli_sum_to_target_tableau(tableau, tableau[pi], p=p, method=method)


def _map_complete_basis_to_target(
    pauli_sum_tableau: TableauType, target_pauli_sum_tableau: TableauType, p: int = 2
) -> TableauType | None:
    """
    Fast full-rank basis map F = P_b^{-1} P'_b.

    The tableaus use row-vector action, so each row p maps as p @ F. When the
    input rows contain a complete independent basis, the image of that basis
    determines F uniquely.
    """
    input_tab = np.asarray(pauli_sum_tableau, dtype=int) % p
    output_tab = np.asarray(target_pauli_sum_tableau, dtype=int) % p

    if input_tab.ndim != 2 or output_tab.ndim != 2 or input_tab.shape != output_tab.shape:
        raise ValueError("Input and target tableaus must be 2-dimensional arrays with matching shape.")

    n_cols = input_tab.shape[1]
    if n_cols % 2 != 0:
        raise ValueError("Pauli tableau width must be even.")

    n_basis = n_cols
    basis_indices, _ = get_linear_dependencies(input_tab, p, compute_dependencies=False)
    if len(basis_indices) < n_basis:
        return None

    basis_indices = basis_indices[:n_basis]
    input_basis = input_tab[basis_indices]
    output_basis = output_tab[basis_indices]

    try:
        if p == 2:
            input_basis_inv = gf2_inv(input_basis)
            F = (input_basis_inv @ output_basis) & 1
            F = np.asarray(F, dtype=int)
        else:
            GF = galois.GF(int(p))
            F_gf = np.linalg.inv(GF(input_basis)) @ GF(output_basis)
            F = np.asarray(F_gf, dtype=int) % p
    except np.linalg.LinAlgError:
        return None

    if not np.array_equal((input_tab @ F) % p, output_tab):
        raise ValueError("Complete-basis map failed to reproduce the target tableau.")

    return F
