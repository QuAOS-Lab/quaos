import random
import numpy as np
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits import Gate
from sympleq.utils import int_to_bases, bases_to_int


def random_pauli_hamiltonian(num_paulis, qudit_dims, mode='rand'):
    """
    Generates a random Pauli Hamiltonian with the given number of Pauli operators plus their Hermitian conjugate pairs.

    Parameters:
        num_paulis (int): Number of Pauli operators to generate.
        qudit_dims (list): List of dimensions for each qudit.
        mode (str): 'rand' or 'uniform' - dictates form of weights in the PauliSum

    Returns:
        tuple: A random PauliSum
    """

    '''
    n_qudits = len(qudit_dims)
    pauli_strings = []
    coefficients = []

    for p in range(num_paulis):
        # np.random.randint(qudit_dims, size=n_qudits)
        x_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]
        # np.random.randint(qudit_dims, size=n_qudits)
        z_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]
        while np.all(np.array(x_exp) == 0) and np.all(np.array(z_exp) == 0):
            x_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]
            z_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]

        x_exp_H = np.zeros_like(x_exp)
        z_exp_H = np.zeros_like(z_exp)
        phase_factor = 1
        pauli_str = ''
        pauli_str_H = ''
    '''
    q2 = np.repeat(qudit_dims, 2)
    available_paulis = list(np.arange(int(np.prod(q2))))

    pauli_strings = []
    coefficients = []

    for _ in range(num_paulis):
        pauli_index = random.choice(available_paulis)
        available_paulis.remove(pauli_index)

        exponents = int_to_bases(pauli_index, q2)
        exponents_H = np.zeros_like(exponents)
        phase_factor = 1
        pauli_str = ' '
        pauli_str_H = ' '

        for j in range(len(qudit_dims)):
            r, s = int(exponents[2 * j]), int(exponents[2 * j + 1])
            pauli_str += f"x{r}z{s} "
            exponents_H[2 * j] = (-r) % qudit_dims[j]
            exponents_H[2 * j + 1] = (-s) % qudit_dims[j]
            pauli_str_H += f"x{exponents_H[2 * j]}z{exponents_H[2 * j + 1]} "

            omega = np.exp(2 * np.pi * 1j / qudit_dims[j])
            phase_factor *= omega**(r * s)

        pauli_strings.append(PauliString.from_string(pauli_str.strip(), dimensions=qudit_dims))
        if mode == 'rand' or mode == 'random':
            coeff = np.random.normal(0, 1) + 1j * np.random.normal(0, 1)
        elif mode == 'uniform' or mode == 'one':
            coeff = 1 + 0 * 1j
        elif mode[0:7] == 'randint':
            # mode is 'randint2', 'randint3', etc.
            d = int(mode[7:])
            coeff = np.random.randint(1, d + 1)

        if not np.array_equal(exponents, exponents_H):
            # random string not Hermitian, add conjugate pair
            conjugate_index = bases_to_int(exponents_H, q2)
            coefficients.append(coeff)
            coefficients.append(np.conj(coeff) * phase_factor)
            available_paulis.remove(conjugate_index)
            pauli_strings.append(PauliString.from_string(pauli_str_H.strip(), dimensions=qudit_dims))
        else:
            coefficients.append(coeff.real)

    rand_ham = PauliSum.from_pauli_strings(pauli_strings, weights=coefficients)
    # rand_ham.combine_equivalent_paulis()
    return rand_ham


def random_pauli_symmetry_hamiltonian(n_qudits: int, n_paulis: int, n_redundant=0,
                                      n_conditional=0, weight_mode='uniform', phase_mode='zero'):
    # 0: I, 1: X, 2: Z, 3: Y
    """
    Generate a random Pauli Hamiltonian with n_qudits qudits and n_paulis Pauli strings,
    with n_redundant redundant qubits and n_conditional conditional qubits.

    The Pauli strings are chosen randomly from the set of strings with the given
    structure. The weights are chosen randomly from the set of real numbers.

    Parameters
    ----------
    n_qudits : int
        The number of qudits in the Hamiltonian.
    n_paulis : int
        The number of Pauli strings in the Hamiltonian.
    n_redundant : int, optional
        The number of redundant qubits. Default is 0.
    n_conditional : int, optional
        The number of conditional qubits. Default is 0.
    weight_mode : str, optional
        The mode of the weights. Can be 'uniform' (default) or 'random'.
    phase_mode : str, optional
        The mode of the phases. Can be 'zero' (default) or 'random'.
    Returns
    -------
    P : PauliSum
        The random Pauli Hamiltonian.

    Examples
    --------
    >>> from sympleq.models.random_hamiltonian import random_pauli_symmetry_hamiltonian
    >>> random_pauli_symmetry_hamiltonian(2, 4)
    PauliSum of size 4x2 with 4 terms and 0 redundant or conditional qubits.
    """
    # TODO: Implementation for Qudits
    # TODO: Make sure that remaining paulis are always unique

    n_rest = n_qudits - n_redundant - n_conditional
    if n_paulis < 2 * n_rest:
        raise ValueError('Too few paulis for full basis with this number of independent qubits')
    # create general structure of the Pauli Hamiltonian before scrambling it with clifford gates
    # redundant qubits
    P = np.zeros((n_paulis, n_qudits), dtype=int)

    # conditional qubits
    for i in range(n_conditional - n_redundant):
        q = np.arange(n_redundant, n_conditional)[i]
        P[2 * n_rest + i, q] = 2  # Z
        for j in range(1, n_paulis - (2 * n_rest + i)):
            P[2 * n_rest + i + j, q] = np.random.choice([0, 2])  # I, Z

    # remaining qubits
    for i in range(n_qudits - (n_redundant + n_conditional)):
        q = np.arange(n_redundant + n_conditional, n_qudits)[i]
        P[i, q] = 1
        P[n_rest + i, q] = 2
        for j in range(1, n_rest - i):
            P[i + j, q] = np.random.choice([0, 1, 2, 3])
            P[n_rest + i + j, q] = np.random.choice([0, 1, 2, 3])

        for j in range(n_paulis - 2 * n_rest):
            P[2 * n_rest + j, q] = np.random.choice([0, 1, 2, 3])

    # Turn P-Matrix into PauliSum
    pauli_strings = ['' for _ in range(n_paulis)]
    for i in range(n_paulis):
        for j in range(n_qudits):
            if P[i, j] == 0:
                pauli_strings[i] += 'x0z0'
            elif P[i, j] == 1:
                pauli_strings[i] += 'x1z0'
            elif P[i, j] == 2:
                pauli_strings[i] += 'x0z1'
            elif P[i, j] == 3:
                pauli_strings[i] += 'x1z1'

            pauli_strings[i] += ' '
        pauli_strings[i] = pauli_strings[i].strip()

    if weight_mode == 'uniform':
        weights = np.ones(n_paulis, dtype=float)
    elif weight_mode == 'random':
        weights = np.random.rand(n_paulis)

    if phase_mode == 'zero':
        phases = np.zeros(n_paulis, dtype=int)
    elif phase_mode == 'random':
        phases = np.random.randint(0, 2, size=n_paulis, dtype=int)

    P = PauliSum.from_string(pauli_strings, dimensions=[2] * n_qudits, weights=weights, phases=phases)

    g = Gate.from_random(n_qudits, 2)
    P = g.act(P, 2)

    return P


def random_gate_symmetric_hamiltonian(G: Gate,
                                      n_qudits: int,
                                      n_paulis: int | None = None,
                                      weight_mode: str = 'uniform',
                                      scrambled: bool = False):
    """
    Generate a random symmetric Hamiltonian from a gate G.

    Parameters
    ----------
    G : Gate
        The gate for which to generate the symmetric Hamiltonian.
    n_qudits : int
        The number of qudits in the resulting Hamiltonian. If None, it is set to G.dimension + 1.
    n_paulis : int
        The number of Pauli strings in the resulting Hamiltonian. If None, it is set to 2 * n_qudits.
    weight_mode : str
        Whether to use 'uniform' or 'random' weights in the Hamiltonian.

    Returns
    -------
    P_sym : PauliSum
        The symmetric Hamiltonian as a PauliSum.

    Notes
    -----
    This function first generates a random Pauli string Hamiltonian, then applies the gate and its inverse to it.
    The sum of the two is the symmetric Hamiltonian. The weights are rounded to 10 decimal places and Pauli strings
    with zero weight are removed.
    """

    if n_qudits <= 0:
        raise ValueError("n_qudits must be greater than zero.")

    if n_paulis is None:
        n_paulis = 2 * n_qudits
    P = random_pauli_symmetry_hamiltonian(n_qudits, n_paulis, 0, 0, weight_mode=weight_mode)
    G_inv = G.inverse()
    qudit_index = np.random.randint(0, n_qudits - 1)
    P_prime = G_inv.act(P, qudit_index)
    P_sym = P + P_prime
    P_sym.phase_to_weight()
    P_sym.combine_equivalent_paulis()
    P_sym.set_weights(np.around(P_sym.weights, decimals=10))
    P_sym.remove_zero_weight_paulis()
    if scrambled is True:
        g = Gate.from_random(n_qudits, 2)
        P_sym = g.act(P_sym, qudit_index)

    return P_sym
