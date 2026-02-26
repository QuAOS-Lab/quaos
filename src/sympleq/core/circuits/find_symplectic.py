"""
This module uses arXiv:1803.06987 to find symplectic solutions.

So far it only works for GF(2), as in the original paper. It could be extended to GF(p).
"""

import numpy as np
from sympleq.core.finite_field_solvers import solve_gf2, solve_linear_system_over_gf
from sympleq.core.circuits.utils import transvection_matrix, symplectic_product_arrays, symplectic_product_matrix


def _symplectic_constraint_row(vec: np.ndarray, p: int = 2) -> np.ndarray:
    """
    Build a row r such that r @ w = <vec, w> mod p, with vec,w in [x|z] ordering.
    """
    n = len(vec) // 2
    row = np.zeros(2 * n, dtype=int)
    row[:n] = (-vec[n:]) % p   # coefficients for w_x
    row[n:] = vec[:n] % p      # coefficients for w_z
    return row


def find_symplectic_solution(u: np.ndarray, v: np.ndarray, p: int = 2) -> np.ndarray:
    """
    Find a binary vector w such that <u,w> = <v,w> = 1 in symplectic inner product.

    Args:
        u, v: Binary vectors of length 2n representing Pauli strings

    Returns:
        Binary vector w of length 2n, or None if no solution exists
    """
    n = len(u) // 2

    # Check for impossible cases first
    if np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n)):
        # Cannot have symplectic inner product 1 with zero vector
        raise ValueError("Cannot find solution with zero vector input")

    # Check if u and v are symplectically independent
    symplectic_product_arrays_uv = symplectic_product_arrays(u, v, p)

    if p == 2 and symplectic_product_arrays_uv == 1:
        # u and v are symplectically independent - use direct construction
        return direct_construction(u, v)

    # General case over GF(p): solve linear constraints directly.
    return solve_general_system(u, v, p=p)


def direct_construction(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    Direct construction when u and v are symplectically independent (<u,v> = 1).

    Uses a geometric construction based on symplectic orthogonality rather than
    solving a linear system.

    Args:
        u, v: Symplectically independent binary vectors

    Returns:
        Solution vector w
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


def solve_general_system(u: np.ndarray, v: np.ndarray, p: int = 2) -> np.ndarray:
    """
    Solve the general case when u and v are symplectically orthogonal (<u,v> = 0).

    Args:
        u, v: Binary vectors that are symplectically orthogonal

    Returns:
        Solution vector w or None if no solution exists
    """
    n = len(u) // 2

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


def find_symplectic_solution_extended(u: np.ndarray, v: np.ndarray,
                                      t_vectors: list | None = None,
                                      p: int = 2) -> np.ndarray:
    """
    Find a binary vector w such that:
    - <u,w> = 1
    - <v,w> = 1
    - <t_i,w> = <t_i,v> for all t_i in t_vectors

    Args:
        u, v: Binary vectors of length 2n representing primary Pauli strings
        t_vectors: List of binary vectors of length 2n for additional constraints

    Returns:
        Binary vector w of length 2n, or None if no solution exists
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


def solve_extended_system(u: np.ndarray, v: np.ndarray, t_vectors: list, p: int = 2) -> np.ndarray:
    """
    Solve the extended system with additional t_i constraints.

    Args:
        u, v: Primary constraint vectors
        t_vectors: Additional constraint vectors

    Returns:
        Solution vector w or None if no solution exists
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


def check_mappable_via_clifford(pauli_sum_tableau: np.ndarray,
                                target_pauli_sum_tableau: np.ndarray,
                                p: int = 2) -> bool:
    sym_check = np.all(
        symplectic_product_matrix(pauli_sum_tableau, p) == symplectic_product_matrix(target_pauli_sum_tableau, p)
    )
    if sym_check:
        return True

    return False


def map_single_pauli_string_to_target(pauli_string_tableau: np.ndarray, target_pauli_string_tableau: np.ndarray,
                                      constraint_paulis: list | None = None,
                                      p: int = 2):
    if p != 2:
        # Compatibility path: use the qudit-capable transvection solver while keeping
        # the legacy function name/signature shape.
        from sympleq.core.circuits.find_symplectic_qudits import Find_transvection_map_solve_extended

        constraints = [] if constraint_paulis is None else list(constraint_paulis)
        sps = [symplectic_product_arrays(t, target_pauli_string_tableau, p) for t in constraints]
        return Find_transvection_map_solve_extended(
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
    pauli_sum_tableau: np.ndarray,
    target_pauli_sum_tableau: np.ndarray,
    p: int = 2,
) -> np.ndarray:
    """
    Map a Pauli sum to a target Pauli sum using symplectic transvections.
    """
    if p != 2:
        from sympleq.core.circuits.find_symplectic_qudits import map_paulisum_to_target_paulisum

        return map_paulisum_to_target_paulisum(pauli_sum_tableau, target_pauli_sum_tableau, p)

    if not check_mappable_via_clifford(pauli_sum_tableau, target_pauli_sum_tableau, p=p):
        raise Exception(f'SPM not equal. Cannot map\n{pauli_sum_tableau} to\n{target_pauli_sum_tableau}')

    m = len(pauli_sum_tableau)
    n = len(pauli_sum_tableau[0]) // 2
    mapped_paulis = []
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
