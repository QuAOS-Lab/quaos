import numpy as np


def bases_to_int(aa: list[int] | np.ndarray,
                 dims: list[int] | np.ndarray) -> int:
    """
    Convert a basis vector in a given dimension to its integer representation.

    Parameters
    ----------
    aa : array-like
        The basis vector to convert. Should be a list or numpy array of integers.
    dims : array-like
        The dimensions of each basis. Should be a list or numpy array of integers.

    Returns
    -------
    int
        The integer representation of the basis vector.

    Examples
    --------
    >>> bases_to_int([1, 0, 2], [2, 2, 3])
    8 # Explanation: (2 * 3^0) + (0 * 2^0*3) + (1 * 2^0*2*3)

    Notes
    -----
    The function interprets `aa` as a vector of coefficients, each referring to an integer basis with dimension
    specified by `dims`: sum(aa[i] * prod(dims[:i]) for i in range(len(aa)))
    """

    # FIXME: Avoid doing extra flips.
    dims = np.flip(dims)
    aa = np.flip(aa)
    a = aa[0] + sum([aa[i1] * np.prod(dims[:i1]) for i1 in range(1, len(dims))])
    # TODO: If deprecated, remove the following lines
    # Following lines commented because they do nothing
    # dims = np.flip(dims)
    # aa = np.flip(aa)
    return a


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
