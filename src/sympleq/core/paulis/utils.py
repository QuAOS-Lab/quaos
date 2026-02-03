from typing import Any
import numpy as np
import re
from .pauli_sum import PauliSum
import networkx as nx
from itertools import product
import sympy as sp


def string_to_symplectic(string: str
                         ) -> tuple[np.ndarray, int]:
    """
    Convert a string representation of a PauliString into its symplectic form.

    Parameters
    ----------
    string : str
        The string representation of the PauliString.

    Returns
    -------
    tuple[np.ndarray, int]
        The symplectic form of the PauliString and its phase.
    """
    substrings = string.split()
    local_symplectics = []
    phases = []
    for s in substrings:
        match = re.match(r'x(\d+)z(\d+)(?:p(\d+))?', s)
        if not match:
            raise ValueError(f"Invalid Pauli string: {s}")
        else:
            x = int(match.group(1))
            z = int(match.group(2))
            p = int(match.group(3)) if match.group(3) is not None else 0
            local_symplectics.append((x, z))
            phases.append(p)

    symplectic = np.array(local_symplectics).T
    return symplectic.flatten(), sum(phases)


def check_mappable_via_clifford(PauliSum: PauliSum,
                                target_PauliSum: PauliSum
                                ) -> bool:
    """
    Checks whether the given PauliSum can be mapped to the target PauliSum via Clifford operations.

    Parameters
    ----------
    PauliSum : PauliSum
        The PauliSum to check.
    target_PauliSum : PauliSum
        The target PauliSum to check against.

    Returns
    -------
    bool
        True if the PauliSum can be mapped to the target PauliSum, False otherwise.
    """
    return bool(
        np.all(
            PauliSum.symplectic_product_matrix() == target_PauliSum.symplectic_product_matrix()
        )
    )


def are_subsets_equal(PauliSum_1: PauliSum,
                      PauliSum_2: PauliSum,
                      subset_1: list[tuple[int, int]],
                      subset_2: list[tuple[int, int]] | None = None
                      ) -> bool:
    """
    Check if two subsets of two PauliSums objects are equal.
    I.e., the first `subset_1' of `PauliSum_1' must match the second `subset_2' of `PauliSum_2'.

    Parameters
    ----------
    PauliSum_1 : PauliSum
        The first PauliSum to compare.
    PauliSum_2 : PauliSum
        The second PauliSum to compare.
    subset_1 : list[tuple[int, int]]
        The indices of the first PauliSum to compare, given as a list of tuples
    subset_2 : list[tuple[int, int]] | None
        The indices of the second PauliSum to compare, given as a list of tuples

    Returns
    -------
    bool
        True if the subsets are equal, False otherwise.
    """
    if subset_2 is None:
        if subset_1 is None:
            # TODO: maybe set `subset_1' to all indices of `PauliSum_1'
            raise ValueError("At least one subset must be provided")
        else:
            subset_2 = subset_1
    else:
        if len(subset_1) != len(subset_2):
            raise ValueError("Subsets must be of the same length")
        if not all(isinstance(i, tuple) and len(i) == 2 for i in subset_1):
            raise ValueError("Subsets must be lists of tuples of length 2")
        if not all(isinstance(i, tuple) and len(i) == 2 for i in subset_2):
            raise ValueError("Subsets must be lists of tuples of length 2")

    for i in range(len(subset_1)):
        if PauliSum_1[subset_1[i]] != PauliSum_2[subset_2[i]]:
            return False
    return True


def commutation_graph(PauliSum: PauliSum,
                      axis: Any | None = None):
    """
    Plots the commutation graph of a PauliSum, based on the adjacency
    matrix given by the symplectic product matrix.

    Parameters
    ----------
    PauliSum : PauliSum
        The PauliSum to plot as a graph.
    labels : list[str] | None
        The labels for the nodes in the graph.
    axis : Any | None
        The axis to plot the graph on.

    Returns
    -------
    nx.Graph
        The commutation graph of the PauliSum.
    """
    adjacency_matrix = PauliSum.symplectic_product_matrix()
    rows, cols = np.where(adjacency_matrix == 1)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()
    all_rows = range(0, adjacency_matrix.shape[0])
    for n in all_rows:
        gr.add_node(n)
    gr.add_edges_from(edges)
    pos1 = nx.spring_layout(gr)

    nx.draw(gr, pos1, node_size=900, with_labels=True, ax=axis)
    return gr


def mod_inv(a: int,
            d: int
            ) -> int:
    """
    Compute the modular multiplicative inverse of an integer.

    Given integers `a` and `d`, this function finds an integer `i` such that
    `(a * i) % d == 1`. If no such integer exists, a `ValueError` is raised.

    Parameters
    ----------
    a : int
        The integer whose modular inverse is to be computed.
    d : int
        The modulus.

    Returns
    -------
    int
        The modular multiplicative inverse of `a` modulo `d`.

    Raises
    ------
    ValueError
        If the modular inverse does not exist (i.e., if `a` and `d` are not coprime).

    Examples
    --------
    >>> mod_inv(3, 11)
    4
    >>> mod_inv(10, 17)
    12
    """
    for i in range(1, d):
        if (a * i) % d == 1:
            return i
    raise ValueError(f"No inverse for {a} mod {d}")


def row_reduce_mod_d(A: np.ndarray,
                     d: int
                     ) -> tuple[np.ndarray, list[int], int]:
    """
    Performs row reduction (Gaussian elimination) of a matrix modulo a given integer.

    This function reduces the input matrix `A` to its row-echelon form over the integers modulo `d`.
    It also returns the list of pivot columns and the rank of the matrix modulo `d`.

    Parameters
    ----------
    A : numpy.ndarray
        The input matrix to be row reduced. Must be a 2D numpy array of integers (dtype=int).
    d : int
        The modulus for the arithmetic operations.

    Returns
    -------
    A_reduced : numpy.ndarray
        The row-reduced form of the input matrix modulo `d`.
    pivots : list of int
        List of column indices that are pivots in the reduced matrix.
    rank : int
        The rank of the matrix modulo `d`.

    Notes
    -----
    - The input matrix `A` must have integer dtype (e.g., np.int32, np.int64).
    - The function assumes that `mod_inv` is defined elsewhere and computes modular inverses.
    - The input matrix `A` is not modified; a copy is used internally.
    - All arithmetic is performed modulo `d`.

    TODO: This is a place where we should test if galois is faster and/or more stable.
    It has inbuilt row_reduce over GF(p)
    """
    A = A.copy() % d
    m, n = A.shape
    rank = 0
    pivots = []
    for col in range(n):
        for row in range(rank, m):
            if A[row, col] != 0:
                break
        else:
            continue
        if row != rank:
            A[[row, rank]] = A[[rank, row]]
        inv = mod_inv(A[rank, col], d)
        A[rank] = (A[rank] * inv) % d
        for r in range(m):
            if r != rank and A[r, col] != 0:
                A[r] = (A[r] - A[r, col] * A[rank]) % d
        pivots.append(col)
        rank += 1
    return A, pivots, rank


def solve_mod_d(A: np.ndarray,
                b: np.ndarray,
                d: int,
                max_solutions: int = 1000
                ) -> list[np.ndarray]:
    """
    Solve a system of linear equations modulo d.
    Given a matrix equation A x = b (mod d), this function finds all possible solutions x
    over the integers modulo d, up to a maximum number of solutions.

    Parameters
    ----------
    A : np.ndarray
        The coefficient matrix of shape (m, n), where m is the number of equations and n is the number of variables.
    b : np.ndarray
        The right-hand side vector of shape (m,).
    d : int
        The modulus for the system of equations.
    max_solutions : int, optional
        The maximum number of solutions to return. Default is 1000.

    Returns
    -------
    list[np.ndarray]
        A list of solutions, where each solution is a numpy array of shape (n,) representing a solution vector x
        such that A @ x % d == b % d.

    Notes
    -----
    - The function uses sympy for symbolic computation and Gaussian elimination modulo d.
    - If the system has infinitely many solutions, only up to `max_solutions` are returned.
    - Free variables are enumerated exhaustively, which may be slow for large systems or large d.

    Examples
    --------
    >>> import numpy as np
    >>> A = np.array([[1, 2], [3, 4]])
    >>> b = np.array([1, 0])
    >>> solve_mod_d(A, b, 5)
    [array([3, 4]), array([0, 2]), array([2, 0]), array([4, 3]), array([1, 1])]
    """

    A_sym = sp.Matrix(A.tolist())
    b_sym = sp.Matrix(b.tolist())
    A_aug = A_sym.row_join(b_sym)
    A_mod = A_aug.applyfunc(lambda x: x % d)

    Ab_rref, pivot_cols = A_mod.rref(iszerofunc=lambda x: x % d == 0, simplify=True)
    n_vars = A.shape[1]

    pivot_cols = [p for p in pivot_cols if p < n_vars]
    free_vars = [i for i in range(n_vars) if i not in pivot_cols]

    solutions = []
    for free_vals in product(range(d), repeat=len(free_vars)):
        sol = [0] * n_vars
        for i, val in zip(free_vars, free_vals):
            sol[i] = val

        for i, pivot in enumerate(pivot_cols):
            rhs = Ab_rref[i, -1]
            lhs = sum(Ab_rref[i, j] * sol[j] for j in range(n_vars)) % d
            sol[pivot] = int((rhs - lhs) % d)

        solutions.append(np.array(sol, dtype=int))
        if len(solutions) >= max_solutions:
            break

    return solutions
