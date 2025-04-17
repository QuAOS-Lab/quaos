import sys
import numpy as np
sys.path.append("./")

from quaos.symplectic import Pauli, PauliString, PauliSum
from quaos.hamiltonian import symplectic_pauli_reduction, random_pauli_hamiltonian


def near_hadamard_symmetries(hamiltonian, threshold=0.1):
    """
    Input Hamiltonian should be of reduced form

    """
    
    # check if there are equal coefficients up to a threshold
    coeffs = hamiltonian.coeffs
    coeffs = np.array(coeffs)
    coeffs = np.abs(coeffs)
    list_of_near_equal_coeffs = [(i, j) for i in range(len(coeffs) - 1) for j in range(i + 1, len(coeffs))
                                 if np.abs(coeffs[i] - coeffs[j]) < threshold]

    # check if these have a Hadamard symmetry

    pass


def find_leading_x(hamiltonian):
    leading_xs = []
    if np.any(hamiltonian.x_exp):
        for i, p in enumerate(hamiltonian.pauli_strings):
            xs = p.x_exp
            if np.any(xs):
                leading_x = len(xs) - 1 - np.argmax(xs[::-1] != 0)  # furthest x to right hand side
            else:
                leading_x = None
            leading_xs.append(leading_x)
    lx_tuple = [(i, v) for i, v in enumerate(leading_xs) if v is not None]
    max_val = max(v for _, v in lx_tuple)
    indices = [i for i, v in lx_tuple if v == max_val]
    values = [v for i, v in lx_tuple if v == max_val]
    return indices, values


def near_pauli_symmetries(hamiltonian):
    """
    Input Hamiltonian should be of reduced form

    """

    # Check for leading X and remove them
    leading_x_string_indexes, leading_x_qudit_indexes = find_leading_x(hamiltonian)
    if max(leading_x_qudit_indexes) == 0:
        # no leading x
        print('No leading x')
        return hamiltonian
    # remove leading x
    hamiltonian_no_lx = hamiltonian.copy()
    leading_x_strings = [hamiltonian.select_pauli_string(lx) for lx in leading_x_string_indexes]
    for lx in leading_x_string_indexes:
        hamiltonian_no_lx.delete_paulis_(lx)
    reduction_circuit = symplectic_pauli_reduction(hamiltonian_no_lx)
    leading_x_strings = [reduction_circuit.act(pauli) for pauli in leading_x_strings]
    hamiltonian_reduced = reduction_circuit.act(hamiltonian_no_lx)
    approximate_hamiltonian = hamiltonian_reduced + np.sum(leading_x_strings)
    return approximate_hamiltonian


if __name__ == "__main__":

    n_qubits = 4
    n_paulis = 4
    dimension_list = [2] * n_qubits
    h = random_pauli_hamiltonian(n_paulis, dimension_list, mode='uniform')
    C = symplectic_pauli_reduction(h)
    
    print('Hamiltonian: \n')
    print(h)

    h_r = C.act(h)
    print('Reduced Hamiltonian: \n')
    print(h_r)

    h_approx = near_pauli_symmetries(h_r)
    print('Approximate Hamiltonian: \n')
    print(h_approx)
