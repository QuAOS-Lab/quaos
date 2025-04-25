import sys
import numpy as np

from symplectic import Pauli, PauliString, PauliSum
from hamiltonian import symplectic_pauli_reduction, random_pauli_hamiltonian
from circuit_utils import find_circuit


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
    return indices, 


def near_pauli_symmetries_naive(hamiltonian):
    """
    Input Hamiltonian should be of reduced form

    """

    # Check for leading X and remove them
    leading_x_string_indexes, leading_x_qudit_index = find_leading_x(hamiltonian)
    if leading_x_qudit_index == 0:
        # no leading x
        print('No leading x')
        return hamiltonian
    # remove leading x
    hamiltonian_no_lx = hamiltonian.copy()
    leading_x_strings = [hamiltonian.select_pauli_string(lx) for lx in leading_x_string_indexes]
    for lx in leading_x_string_indexes:
        hamiltonian_no_lx._delete_paulis(lx)
    reduction_circuit = symplectic_pauli_reduction(hamiltonian_no_lx)
    leading_x_strings = [reduction_circuit.act(pauli) for pauli in leading_x_strings]
    hamiltonian_reduced = reduction_circuit.act(hamiltonian_no_lx)
    approximate_hamiltonian = hamiltonian_reduced
    for pauli in leading_x_strings:
        approximate_hamiltonian += pauli
    return approximate_hamiltonian


def near_pauli_symmetries(hamiltonian):
        
    leading_x_string_indexes, leading_x_qudit = find_leading_x(hamiltonian)

    norms = np.abs(hamiltonian.coeffs)**2
    norm_indexes = np.argsort(norms)
    
    # the goal is now to move the leading x to the minimal norm pauli
    if len(leading_x_string_indexes) == 0:
        # no leading x
        print('No leading x')
        return hamiltonian
    elif len(leading_x_string_indexes) == 1:
        # only one leading x
        min_index = norm_indexes[0]
        leading_x_string = leading_x_string_indexes[0]
        for qudit in range(0, leading_x_qudit):  # from the first qudit to the qudit before the leading x
            local_pauli_sum = hamiltonian[leading_x_string, qudit] @ hamiltonian[leading_x_string, leading_x_qudit] + hamiltonian[min_index, qudit] @ hamiltonian[min_index, leading_x_qudit]
            target_pauli_sum = 
    else:
        # multiple leading x
        pass




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
