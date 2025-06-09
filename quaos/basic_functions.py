import numpy as np
from quaos.paulis import PauliSum


def read_luca_test_2(path: str, dims: list[int] | int = 2, spaces: bool = True):
    """Reads a Hamiltonian file and parses the Pauli strings and coefficients.

    Args:
        path (str): Path to the Hamiltonian file.
        dims (Union[int, List[int]]): Dimension(s) of the qudits, default is 2.
        spaces (bool): Whether to expect spaces in the Pauli string format, default is True.

    Returns:
        Tuple: Parsed Pauli operators and corresponding coefficients.
    """
    with open(path, "r") as f:
        lines = f.readlines()

    pauli_strings = []
    coefficients = []

    for line in lines:
        pauli_list = line.split(', {') if spaces else line.split(',{')
        coeff = pauli_list[0][1:].replace(" ", "").replace('*I', 'j')
        coefficients.append(complex(coeff))

        pauli_str = ' '.join(f"x{item.count('X')}z{item.count('Z')}" for item in pauli_list[1:])
        pauli_strings.append(pauli_str.strip())

    return PauliSum(pauli_strings, weights=coefficients, dimensions=dims)


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
    dimensions = np.flip(dimensions)
    base = np.flip(base)
    number = base[0] + sum([base[qudit] * np.prod(dimensions[:qudit]) for qudit in range(1, len(dimensions))])
    dimensions = np.flip(dimensions)
    base = np.flip(base)
    return number


def int_to_bases(number, dimensions) -> np.ndarray:
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
    dimensions = np.flip(dimensions)
    base = [number % dimensions[0]]
    for i in range(1, len(dimensions)):
        s0 = base[0] + sum([base[i1] * dimensions[i1 - 1] for i1 in range(1, i)])
        s1 = np.prod(dimensions[:i])
        base.append(((number - s0) // s1) % dimensions[i])
    dimensions = np.flip(dimensions)
    return np.flip(np.array(base))


# PHYSICS FUNCTIONS


def Hamiltonian_Mean(P: PauliSum, psi: np.ndarray) -> float:
    """Returns the mean of a Hamiltonian with a given state.

    Args:
        P: pauli, Paulis of Hamiltonian
        psi: numpy.array, state for mean

    Returns:
        numpy.float64, mean sum(c*<psi|P|psi>)
    """
    p = P.n_paulis()
    psi_dag = psi.conj().T
    return sum(P.weights[i] * psi_dag @ P.matrix_form(i) @ psi for i in range(p))


def covariance_matrix(P: PauliSum, psi: np.ndarray) -> np.ndarray:
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
    p = P.n_paulis()
    cc = P.weights
    mm = [P.matrix_form(i) for i in range(p)]
    psi_dag = psi.conj().T
    cc1 = [psi_dag @ mm[i] @ psi for i in range(p)]
    cc2 = [psi_dag @ mm[i].conj().T @ psi for i in range(p)]
    return np.array([[np.conj(cc[i0]) * cc[i1] * ((psi_dag @ mm[i0].conj().T @ mm[i1] @ psi) - cc2[i0] * cc1[i1]) for i1 in range(p)] for i0 in range(p)])


def commutation_graph(P: PauliSum) -> np.ndarray:
    """
    Computes the commutation graph for a given set of Pauli operators.

    Args:
        P (PauliSum): A set of Pauli operators represented as a PauliSum object.

    Returns:
        np.ndarray: A 2D numpy array representing the commutation graph. Each element [i, j] is 1 if the i-th and 
                    j-th Pauli operators commute, otherwise 0.
    """
    p = P.n_paulis()
    return np.array([[int(P[i0].commute(P[i1])) for i1 in range(p)] for i0 in range(p)])


def complex_phase_value(phase, dimension):
    """
    Computes the a-th eigenvalue of a pauli with dimension d.

    Args:
        a (int): The integer to compute the eigenvalue for.
        d (int): The dimension of the pauli to use.

    Returns:
        complex: The computed eigenvalue.
    """
    return np.exp(2 * np.pi * 1j * phase / dimension)


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
    normalized_state = np.sqrt(gamma_sample / np.sum(gamma_sample)) * np.exp(1j * phases)
    return normalized_state


