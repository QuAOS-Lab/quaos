import numpy as np
from numpy.typing import NDArray
from .pauli_object import PauliObject
from .pauli_sum import PauliSum


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
    source_symplectic = PauliSum.symplectic_product_matrix()
    target_symplectic = target_PauliSum.symplectic_product_matrix()

    if source_symplectic.shape != target_symplectic.shape:
        return False

    return bool(np.all(source_symplectic == target_symplectic))


def mod_inv(a: int | np.int64 | NDArray | np.integer,
            d: int | np.int64 | NDArray | np.integer
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
    if not isinstance(a, int):
        a = int(a)
    if not isinstance(d, int):
        d = int(d)
    inv = pow(a, -1, d)

    return inv


# PHYSICS FUNCTIONS
def hamiltonian_mean(P: PauliObject, psi: np.ndarray) -> float:
    """Returns the mean of a Hamiltonian with a given state.

    Args:
        P: pauli, Paulis of Hamiltonian
        psi: numpy.array, state for mean

    Returns:
        numpy.float64, mean sum(c*<psi|P|psi>)
    """
    mu = np.real(np.transpose(np.conjugate(psi)) @ P.to_hilbert_space() @ psi)
    # FIXME: better modify the input, saying psi is complex array?
    # FIXME: should not be necessary to specify float, should be fixed with new formatting PR.
    return float(mu)


def covariance_matrix(P: PauliObject, psi: np.ndarray) -> np.ndarray:
    """
    Computes the covariance matrix for a given set of Pauli operators and a quantum state.

    Args:
        P (PauliSum): The set of Pauli operators, represented as a PauliSum object, with associated weights.
        psi (np.ndarray): The state vector for which the covariance matrix is computed.

    Returns:
        np.ndarray: A 2D numpy array representing the covariance matrix of the Pauli operators with respect to
                    the given state. Each element [i, j] corresponds to the covariance between the i-th and j-th
                    Pauli operators.
    """
    n_paulis = P.n_paulis()
    # Recall that weights and phases are already included into the elements of pauli_strings!
    pauli_strings = [P.to_hilbert_space(i) for i in range(n_paulis)]
    psi_dag = psi.conj().T
    covariance_matrix = np.array(
        [
            [
                (psi_dag @ pauli_strings[i0].conj().T @ pauli_strings[i1] @ psi) -
                (psi_dag @ pauli_strings[i0].conj().T @ psi) * (psi_dag @ pauli_strings[i1] @ psi)
                for i1 in range(n_paulis)]
            for i0 in range(n_paulis)]
    )
    return covariance_matrix
