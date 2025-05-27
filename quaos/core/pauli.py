"""
Module: pauli.py

This module provides functionalities for managing and manipulating sets of Pauli operators using symplectic matrices.
"""

import numpy as np
import itertools
import re
import math 
import scipy


class pauli:
    """A class for storing sets of Pauli operators as pairs of symplectic matrices"""	

    def __init__(self, 
                 X: np.ndarray, 
                 Z: np.ndarray, 
                 dims: int | list[int] = 2, 
                 phases: int | list[int] = 0):
        """
        Initialize the Pauli object with symplectic matrices.

        Args:
            X (numpy.array): X-part of Pauli in symplectic form with shape (p, q).
            Z (numpy.array): Z-part of Pauli in symplectic form with shape (p, q).
            dims (int or list): Dimensions of the qudits, default is 2.
            phases (int or list): Phases for the Pauli operators, default is 0.

        Raises:
            Exception: If X and Z do not have the same shape.
        """
        if X.shape != Z.shape:
            raise Exception("X- and Z-parts must have same shape")
        if type(dims) == int:
            dims = np.array([dims] * X.shape[1])
        elif type(dims) == list:
            dims = np.array(dims)
        if type(phases) == int:
            phases = np.array([phases] * X.shape[0])
        elif type(phases) == list:
            phases = np.array(phases)
        self.X = X % dims
        self.Z = Z % dims
        self.lcm = np.lcm.reduce(dims)
        self.dims = dims
        self.phases = phases % self.lcm

    def is_IX(self) -> bool:
        """Check whether self has only X component.

        Returns:
            bool: True if self has only X component, False otherwise.
        """
        return not np.any(self.Z)

    def is_IZ(self) -> bool:
        """Check whether self has only Z component.

        Returns:
            bool: True if self has only Z component, False otherwise.
        """
        return not np.any(self.X)

    def is_commuting(self) -> bool:
        """Check whether the set of Paulis are pairwise commuting.

        Returns:
            bool: True if self is pairwise commuting set of Paulis, False otherwise.
        """
        p = self.paulis()
        PP = [self.a_pauli(i) for i in range(p)]
        return not any(symplectic_inner_product(PP[i0], PP[i1]) for i0, i1 in itertools.combinations(range(p), 2))

    def is_quditwise_commuting(self) -> bool:
        """
        Check whether the set of Paulis are pairwise commuting on every qudit.

        Returns:
            bool: True if self is a pairwise quditwise commuting set of Paulis.
        """
        p = self.paulis()
        PP = [self.a_pauli(i) for i in range(p)]
        return not any(quditwise_inner_product(PP[i0], PP[i1]) for i0, i1 in itertools.combinations(range(p), 2))

    def a_pauli(self, a: int) -> 'pauli':
        """
        Pull out the ath Pauli from self.

        Args:
            a (int): Index of Pauli to be returned.

        Returns:
            pauli: The ath Pauli in self.
        """
        return pauli(np.array([self.X[a, :]]), 
                     np.array([self.Z[a, :]]),
                     self.dims,
                     np.array([self.phases[a]]))

    def paulis(self) -> int:
        """Count the number of Paulis in self.

        Returns:
            int: The number of Paulis in self.
        """
        return self.X.shape[0]

    def qudits(self) -> int:
        """Count the number of qudits in self.

        Returns:
            int: The number of qudits in self.
        """
        return self.X.shape[1]

    def delete_paulis_(self, aa: int | list[int]) -> 'pauli':
        """
        Delete Paulis indexed by aa.

        Args:
            aa (list of int): Indices of the Paulis to be deleted.

        Returns:
            pauli: The updated pauli object after deletion.
        """
        if type(aa) is int:
            self.X = np.delete(self.X, aa, axis=0)
            self.Z = np.delete(self.Z, aa, axis=0)
            self.phases = np.delete(self.phases, aa)
        else:
            for a in sorted(aa, reverse=True):
                self.X = np.delete(self.X, a, axis=0)
                self.Z = np.delete(self.Z, a, axis=0)
                self.phases = np.delete(self.phases, a)
        return self

    def delete_qudits_(self, aa: int | list[int]) -> None:
        """
        Delete qudits indexed by aa.

        Args:
            aa (list of int): Indices of the qudits to be deleted.

        Returns:
            pauli: The updated pauli object after deletion.
        """
        if type(aa) is int:
            self.X = np.delete(self.X, aa, axis=1)
            self.Z = np.delete(self.Z, aa, axis=1)
            self.dims = np.delete(self.dims, aa)
        else:
            for a in sorted(aa, reverse=True):
                self.X = np.delete(self.X, a, axis=1)
                self.Z = np.delete(self.Z, a, axis=1)
                self.dims = np.delete(self.dims, a)
        self.lcm = np.lcm.reduce(self.dims)

    def copy(self) -> 'pauli':
        """
        Return deep copy of self.

        Returns:
            pauli: A deep copy of self.
        """
        X = np.array([[self.X[i0, i1] for i1 in range(self.qudits())] for i0 in range(self.paulis())])
        Z = np.array([[self.Z[i0, i1] for i1 in range(self.qudits())] for i0 in range(self.paulis())])
        dims = np.array([self.dims[i] for i in range(self.qudits())])
        phases = np.array([self.phases[i] for i in range(self.paulis())])
        return pauli(X, Z, dims, phases)

    def print(self) -> None:
        """Print string representation of self."""
        sss, dims, phases = pauli_to_string(self)
        for ss in sss:
            print(ss)

    def print_symplectic(self) -> None:
        """Print the symplectic representation of the Pauli object."""
        print(''.join(str(int(i)) for i in self.dims), 
              ''.join(str(int(i)) for i in self.dims))

        print('-' * self.qudits(), '-' * self.qudits())

        for i in range(self.paulis()):
            print(''.join(str(int(i1)) for i1 in self.X[i, :]) + ' ' + ''.join(str(int(i1)) for i1 in self.Z[i, :]) + ' ' + str(self.phases[i]) + '/' + str(self.lcm))


def pauli_to_matrix(P: 'pauli') -> scipy.sparse.csr_matrix:
    """Convert a Pauli object to its matrix representation.

    Parameters
    ----------
    P : pauli
        Pauli object to be converted. Must have shape (1,q).

    Returns
    -------
    scipy.sparse.csr_matrix
        Matrix representation of input Pauli.
    """
    if P.paulis() != 1:
        raise Exception("Matrix can only be constructed for a single Pauli")
    X, Z, dims, phase = P.X[0], P.Z[0], P.dims, P.phases[0]
    return math.e**(phase * 2 * math.pi * 1j / P.lcm) * tensor([XZ_mat(dims[i], X[i], Z[i]) for i in range(P.qudits())])


def string_to_pauli(sss: str | list[str], 
                    dims: int | list[int] = 2, 
                    phases: int | list[int] = 0) -> 'pauli':
    """
    Convert a collection of strings (or single string) to a pauli object.

    Parameters
    ----------
    sss : str or list[str]
        String representation of Pauli. If a single string, it is a single
        Pauli. If a list of strings, each string is a Pauli.

    Returns
    -------
    pauli
        Pauli corresponding to input string(s).
    """
    if type(sss) is str:
        X = np.array([[re.split("x|z", s)[1] for s in sss.split()]], dtype=int)
        Z = np.array([[re.split("x|z", s)[2] for s in sss.split()]], dtype=int)
        return pauli(X, Z, dims, phases)
    else:
        X = np.array([[re.split('x|z', s)[1] for s in ss.split()] for ss in sss], dtype=int)
        Z = np.array([[re.split('x|z', s)[2] for s in ss.split()] for ss in sss], dtype=int)
        return pauli(X, Z, dims, phases)


def pauli_to_string(P: 'pauli') -> str | list[str]:
    """Convert a pauli object to a collection of strings (or single string).

    Parameters
    ----------
    P : pauli
        Pauli to be stringified

    Returns
    -------
    str or list[str]
        String representation of Pauli
    """
    X, Z = P.X, P.Z
    return [''.join('x' + str(X[i0, i1]) + 'z' + str(Z[i0, i1]) + ' ' for i1 in range(P.qudits()))[:-1] for i0 in range(P.paulis())], P.dims, P.phases


def symplectic_inner_product(P0: 'pauli', P1: 'pauli') -> bool:
    """Compute the symplectic inner product of two pauli objects (each with a single Pauli).
    For two Pauli operators :math:`P_0` and :math:`P_1`, each of which is represented in terms of its X and Z components:
    
    .. math::

        P = (X, Z)

    the symplectic inner product is defined as:

    .. math::

        \\langle P_0, P_1 \\rangle = X_0 \\cdot Z_1 - Z_0 \\cdot X_1

    where :math:`X` and :math:`Z` are binary vectors representing the presence of Pauli :math:`X` and :math:`Z` components.
    This function checks whether the result is **odd (True) or even (False)**, which determines whether two Pauli operators **commute or anti-commute**:

    - **0 (False)** → Operators commute.
    - **1 (True)** → Operators anti-commute.

    Args:
        P0 (pauli): A Pauli object with shape (1, q).
        P1 (pauli): A Pauli object with shape (1, q).

    Returns:
        bool: True if the symplectic inner product is odd (anti-commuting), False otherwise (commuting).

    References:
        Bandyopadhyay, et al. *"A new proof for the existence of mutually unbiased bases."* 
        Available at: `arXiv:quant-ph/0103162 <https://arxiv.org/abs/quant-ph/0103162>`_
    """
    if (P0.paulis() != 1) or (P1.paulis() != 1):
        raise Exception("Symplectic inner product only works with pair of single Paulis")
    if any(P0.dims - P1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    tmp = np.sum((P0.X * P1.Z - P0.Z * P1.X) * np.array([P0.lcm] * len(P0.dims)) // P0.dims)
    return bool(tmp % P0.lcm)


def quditwise_inner_product(P0: 'pauli', P1: 'pauli') -> bool:
    """The quditwise inner product of two pauli objects (each with a single Pauli).

    Parameters
    ----------
    P0 : pauli
        Pauli object with shape (1,q)
    P1 : pauli
        Pauli object with shape (1,q)

    Returns
    -------
    bool
        Quditwise inner product of Paulis
    """
    if (P0.paulis() != 1) or (P1.paulis() != 1):
        raise Exception("Quditwise inner product only works with pair of single Paulis")
    if any(P0.dims - P1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    return any(np.sum(P0.X[0, i] * P1.Z[0, i] - P0.Z[0, i] * P1.X[0, i]) % P0.dims[i] for i in range(P0.qudits()))


def pauli_product(P0: 'pauli', P1: 'pauli') -> 'pauli':
    """
    The product of two pauli objects.

    Parameters
    ----------
    P0 : pauli
        Pauli object with shape (1,q)
    P1 : pauli
        Pauli object with shape (1,q)

    Returns
    -------
    pauli
        Product of Paulis
    """
    if P0.paulis() != 1 or P1.paulis() != 1:
        raise Exception("Product can only be calculated for single Paulis")
    if any(P0.dims - P1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    return pauli(P0.X + P1.X, P0.Z + P1.Z, P0.dims, P0.phases + P1.phases)


def tensor(mm: list[scipy.sparse.csr_matrix]) -> scipy.sparse.csr_matrix:
    """Function for computing the tensor product of a list of matrices.

    Parameters
    ----------
    mm : list{scipy.sparse.csr_matrix}
        List of matrices to tensor

    Returns
    -------
    scipy.sparse.csr_matrix
        Tensor product of matrices
    """
    if len(mm) == 0:
        return matrix([])  # ???
    elif len(mm) == 1:
        return mm[0]
    else:
        return scipy.sparse.kron(mm[0], tensor(mm[1:]), format="csr")


def XZ_mat(d: int, aX: int, aZ: int) -> scipy.sparse.csr_matrix:
    """Function for creating generalized Pauli matrix.

    Parameters
    ----------
    d : int
        Dimension of the qudit
    aX : int
        X-part of the Pauli matrix
    aZ : int
        Z-part of the Pauli matrix

    Returns
    -------
    scipy.sparse.csr_matrix
        Generalized Pauli matrix
    """
    omega = math.e ** (2 * math.pi * 1j / d)
    aa0 = np.array([1 for i in range(d)])
    aa1 = np.array([i for i in range(d)])
    aa2 = np.array([(i - aX) % d for i in range(d)])
    X = scipy.sparse.csr_matrix((aa0, (aa1, aa2)))
    aa0 = np.array([omega**(i * aZ) for i in range(d)])
    aa1 = np.array([i for i in range(d)])
    aa2 = np.array([i for i in range(d)])
    Z = scipy.sparse.csr_matrix((aa0, (aa1, aa2)))
    if (d == 2) and (aX % 2 == 1) and (aZ % 2 == 1):
        return 1j * (X @ Z)
    return X @ Z
