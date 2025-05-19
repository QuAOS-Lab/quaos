from quspin.operators import hamiltonian, quantum_operator
from quspin.basis import spin_basis_1d
import numpy as np



def pauli_sum_to_quspin(hamiltonian_ps):
    if not np.all(hamiltonian_ps.dimensions == 2):
        raise ValueError("PauliString must be a 2-dimensional PauliString")
    
    L = hamiltonian_ps.n_qudits()
    basis = spin_basis_1d(L, kblock=0, kshift=0, pbc=False)
    for ps in hamiltonian_ps.pauli_strings:
        h_eff_statis = pauli_string_to_quspin(ps)




def pauli_string_to_quspin(pauli_string):
    """
    Convert a Pauli string to a list of lists of the form to input to the quspin hamiltonian constructor.

    this needs to be a list of lists of the form J_list = [[coeff, i, j, ...], ...]
    where coeff is the coefficient of the Pauli string, and i, j, ... are the indices of the Pauli operators.

    then need the operator list
    [['op', J_list], ...]

    where 'op' is 'xx', 'yy', 'zz', 'x', 'y', 'z' or 'I' or even larger n-body terms like 'xyzz...zzy...'
    """
    L = pauli_string.n_qudits()
    