from __future__ import annotations
import numpy as np
import galois
from typing import TYPE_CHECKING

'''
if TYPE_CHECKING:
    from sympleq.core.paulis import PauliSum
    from sympleq.core.paulis import PauliObject
'''


def bases_to_int(base, dimensions) -> int:
    """
    Converts a list of integers (base) given the dimensions to an integer. Base can be thought of as a number
    in basis of the dimensions which is converted to a number in base 10.

    The base is a list of integers, where the i-th element is the value
    in the i-th dimension. The dimensions parameter is a list of integers,
    where the i-th element is the size of the i-th dimension.

    The function returns the integer that corresponds to the input base in
    the given dimensions.

    Parameters
    ----------
    base : list of int
        The base to be converted.
    dimensions : list of int
        The dimensions of the base.

    Returns
    -------
    int
        The integer that corresponds to the input base in the given
        dimensions.
    """
    # FIXME: maybe there is a way to avoid flipping twice?
    dimensions = np.flip(dimensions)
    base = np.flip(base)
    number = base[0] + sum([base[qudit] * np.prod(dimensions[:qudit])
                           for qudit in range(1, len(dimensions))])
    dimensions = np.flip(dimensions)
    base = np.flip(base)
    return number


def int_to_bases(number: int, dimensions: int | list[int] | np.ndarray) -> np.ndarray:
    """
    Converts an integer to a list of integers given the dimensions. The returned list of integers can be thought of
    as a number in basis of the dimensions which is converted from a number in base 10.

    The function takes two parameters, an integer and a list of integers, where the i-th element of the list is the
    size of the i-th dimension.

    The function returns the list of integers that corresponds to the input number in the given dimensions.

    Parameters
    ----------
    number : int
        The number to be converted.
    dimensions : list of int
        The dimensions of the base.

    Returns
    -------
    np.ndarray
        The list of integers that corresponds to the input number in the given dimensions.
    """
    if isinstance(dimensions, int):
        dimensions = [dimensions]

    # FIXME: maybe there is a way to avoid flipping twice?
    dims = np.flip(dimensions)
    base = [number % dims[0]]
    for i in range(1, len(dimensions)):
        s0 = base[0] + sum([base[i1] * dims[i1 - 1]
                           for i1 in range(1, i)])
        s1 = np.prod(dims[:i])
        base.append(((number - s0) // s1) % dims[i])
    return np.flip(np.array(base, dtype=int))


def get_linearly_independent_rows(A: np.ndarray, d: int) -> list[int]:
    """
    Returns the pivot column indices for the row-reduced form of matrix A over a Galois field.

    Args:
        A (galois.FieldArray): Input matrix over GF(p).

    Returns:
        List[int]: List of pivot column indices.
    """

    GF = galois.GF(d)
    A = GF(A)
    R = A.row_reduce()
    pivots = []
    for row in R:
        nz_indices = np.nonzero(row)[0]
        if nz_indices.size > 0:
            pivots.append(nz_indices[0])
    return pivots


def complex_phase_value(phase: int, dimension: int, full: bool = True) -> complex:
    """
    Computes the a-th eigenvalue of a pauli with dimension d.

    Args:
        phase (int): The integer to compute the eigenvalue for.
        dimension (int): The dimension of the pauli to use.

    Returns:
        complex: The computed eigenvalue.
    """
    if full:
        return np.exp(2 * np.pi * 1j * phase / (2 * dimension))
    else:
        return np.exp(2 * np.pi * 1j * phase / dimension)


'''
def rand_state(dimension):
    """
    Generate a random quantum state vector for a system of dimension d.

    Args:
        d (int): Dimension of the quantum system.

    Returns:
        np.ndarray: A normalized random state vector in the complex space of size d.
    """
    gamma_sample = np.random.gamma(1, 1, int(dimension))
    phases = np.random.uniform(0, 2 * np.pi, int(dimension))
    normalized_state = np.sqrt(
        gamma_sample / np.sum(gamma_sample)) * np.exp(1j * phases)
    return normalized_state
'''
