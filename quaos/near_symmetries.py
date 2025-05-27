import sys
import numpy as np
from gates import Circuit
from symplectic import Pauli, PauliString, PauliSum
from hamiltonian import symplectic_pauli_reduction, random_pauli_hamiltonian
from circuit_utils import find_circuit, find_agnostic_circuit
from pauli_utils import concatenate_pauli_sums

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
        return indices, int(max_val)
    else:
        print('no leading x')
        return [], None


def find_leading_z(hamiltonian):
    leading_zs = []
    if np.any(hamiltonian.z_exp):
        for i, p in enumerate(hamiltonian.pauli_strings):
            zs = p.z_exp
            if np.any(zs):
                leading_x = len(zs) - 1 - np.argmax(zs[::-1] != 0)  # furthest x to right hand side
            else:
                leading_x = None
            leading_zs.append(leading_x)
    lx_tuple = [(i, v) for i, v in enumerate(leading_zs) if v is not None]
    max_val = max(v for _, v in lx_tuple)
    indices = [i for i, v in lx_tuple if v == max_val]
    return indices, int(max_val)


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


def near_pauli_symmetries(hamiltonian: PauliSum) -> Circuit | None:
            
    leading_x_string_indexes, leading_x_qudit = find_leading_x(hamiltonian)
    if leading_x_qudit is not None:
        if leading_x_qudit == 0:
            # no leading x
            print('No leading x')
            return None

        dim_leading_x = hamiltonian.dimensions[leading_x_qudit]
            
        norms = np.abs(hamiltonian.weights)**2
        norm_indexes = np.argsort(norms)
        
        # the goal is now to move the leading x to the minimal norm pauli
        nearly_symmetrising_circuit = Circuit(dimensions=hamiltonian.dimensions, gates=[])  # empty circuit of the correct dimensions
        if len(leading_x_string_indexes) == 1:
            # only one leading x
            min_index = int(norm_indexes[0])
            if min_index == leading_x_string_indexes[0]:
                print('Already in the right place')
                return nearly_symmetrising_circuit
            else:
                leading_x_string = leading_x_string_indexes[0]
                for qudit in range(0, leading_x_qudit):  # from the first qudit to the qudit before the leading x

                    local_pauli_sum = concatenate_pauli_sums([hamiltonian[0:hamiltonian.n_paulis(), qudit],
                                                              hamiltonian[0:hamiltonian.n_paulis(), leading_x_qudit]])
                    local_pauli_sum.weights = hamiltonian.weights

                    print('local pauli \n', local_pauli_sum)
                    target_pauli_list = [Pauli(0, 0, dim_leading_x)] * hamiltonian.n_paulis()
                    target_pauli_list[min_index] = Pauli(1, 0, dim_leading_x)
                    target_pauli_list[leading_x_string] = Pauli(0, 1, dim_leading_x)
                    target_indices = [(i, 1) for i in range(hamiltonian.n_paulis())]
                    C = find_agnostic_circuit(local_pauli_sum.copy(), target_pauli_list, target_indices, 100)

                    if C is not None:
                        print('found circuit \n')
                        print('local pauli \n', local_pauli_sum)
                        print(C[0])
                        print(C[0].act(local_pauli_sum))
                        nearly_symmetrising_circuit.embed_circuit(C[0], [qudit, leading_x_qudit])
                        return nearly_symmetrising_circuit
        else:
            raise NotImplementedError('More than one leading x not implemented yet')
        print('No circuit found')
    return nearly_symmetrising_circuit


if __name__ == "__main__":

    n_qubits = 4
    n_paulis = 4
    dimension_list = [2] * n_qubits
    h = random_pauli_hamiltonian(n_paulis, dimension_list, mode='random')
    C = symplectic_pauli_reduction(h)
    
    print('Hamiltonian: \n')
    print(h)

    h_r = C.act(h)
    print('Reduced Hamiltonian: \n')
    print(h_r)

    C_ns = near_pauli_symmetries(h_r)
    print(C_ns)
    print('Standard form \n')
    print(h_r)
    print('Approximate Hamiltonian: \n')

    print(C_ns.act(h_r) if C_ns is not None else 'No circuit found')
