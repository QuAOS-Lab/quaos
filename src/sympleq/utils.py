from __future__ import annotations
from typing import Sequence
import numpy as np
import galois

from sympleq._typing import ComplexNDArray, IntArrayLike, IntNDArray


def bases_to_int(base: IntArrayLike, dimensions: IntArrayLike) -> int:
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


def int_to_bases(number: int, dimensions: IntArrayLike) -> IntNDArray:
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


# TODO: This function is not used anywhere in the codebase. Remove it in a following PR.
def get_linearly_independent_rows(matrix: IntNDArray, dimension: int) -> list[int]:
    """
    Returns the pivot column indices for the row-reduced form of `matrix` over a Galois field.

    Parameters
    ----------
    matrix : IntNDArray
        Input matrix over GF(dimension).
    dimension : int
        The prime (or prime power) defining the Galois field GF(dimension).

    Returns
    -------
    list of int
        List of pivot column indices.
    """

    field = galois.GF(dimension)
    matrix = field(matrix)
    reduced_matrix = matrix.row_reduce()
    pivots = []
    for row in reduced_matrix:
        nz_indices = np.nonzero(row)[0]
        if nz_indices.size > 0:
            pivots.append(nz_indices[0])
    return pivots


def complex_phase_value(phase: int, dimension: int) -> complex:
    """
    Roots of unity (varying `phase`) with respect to (twice a) chosen dimension `dimension`.

    The "twice" is for taking into account the qubit case (`dimension = 2`), where X*Z = i Y.
    For details, see: `IEEE International Symposium on Information Theory (ISIT), pp. 791-795.
    IEEE (2018) <https://doi.org/10.1109/ISIT.2018.8437652>`_

    Parameters
    ----------
    phase : int
        The integer to compute the eigenvalue for.
    dimension : int
        The dimension of the pauli to use.

    Returns
    -------
    complex
        The computed eigenvalue.
    """
    phase = phase % (2 * dimension)

    # Avoid roundoff for the common quadrant roots 1, i, -1, -i.
    if (2 * phase) % dimension == 0:
        quadrant = ((2 * phase) // dimension) % 4
        return 1j ** quadrant
    else:
        return np.exp(2 * np.pi * 1j * phase / (2 * dimension))


def multi_kron(matrices: Sequence[ComplexNDArray]) -> ComplexNDArray:
    """
    Compute the Kronecker product of multiple matrices.

    Parameters
    ----------
    matrices : Sequence[ComplexNDArray]
        A sequence of complex square matrices to compute the Kronecker product of.

    Returns
    -------
    ComplexNDArray
        The Kronecker product of the input matrices.

    Raises
    ------
    ValueError
        If `matrices` is empty, or if any matrix is not a 2-dimensional square
        array with a complex dtype.
    """
    if not matrices:
        raise ValueError("At least one matrix must be provided.")

    for m in matrices:
        if not isinstance(m, np.ndarray) or m.ndim != 2 or m.shape[0] != m.shape[1]:
            shape = getattr(m, "shape", None)
            raise ValueError(f"Each matrix must be a square 2-dimensional array, got shape {shape}.")
        if not np.issubdtype(m.dtype, np.complexfloating):
            raise ValueError(f"Each matrix must have a complex dtype, got {m.dtype}.")

    if len(matrices) == 1:
        return matrices[0]
    M = np.kron(matrices[0], matrices[1])
    for i in range(2, len(matrices)):
        M = np.kron(M, matrices[i])
    return M
