import numpy as np
from sympleq.core.paulis import PauliSum


def true_mean(H: PauliSum, psi) -> float:
    """
    Calculate the true mean of a Hamiltonian given a PauliSum object and a state vector.

    Parameters:
        H (PauliSum): The PauliSum object representing the Hamiltonian.
        psi (numpy.ndarray): The state vector.

    Returns:
        float: The true mean of the Hamiltonian.
    """
    mu = np.real(np.transpose(np.conjugate(psi)) @ H.to_hilbert_space() @ psi)
    return mu


def true_covariance_graph(H: PauliSum, psi) -> np.ndarray:
    """
    Calculate the true covariance matrix of a Hamiltonian given a PauliSum object and a state vector.

    Parameters:
        H (PauliSum): The PauliSum object representing the Hamiltonian.
        psi (numpy.ndarray): The state vector.

    Returns:
        numpy.ndarray: The true covariance matrix of the Hamiltonian.
    """
    p = H.n_paulis()
    mm = [H.to_hilbert_space(i) for i in range(p)]
    psi_dag = psi.conj().T
    cc1 = [psi_dag @ mm[i] @ psi for i in range(p)]
    cc2 = [psi_dag @ mm[i].conj().T @ psi for i in range(p)]
    cm = np.zeros((p, p), dtype=complex)
    for i0 in range(p):
        for i1 in range(p):
            cov = (psi_dag @ mm[i0].conj().T @ mm[i1] @ psi) - cc2[i0] * cc1[i1]
            cm[i0, i1] = cov
    return cm
