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
