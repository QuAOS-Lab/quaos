import numpy as np
import sys
sys.path.append("./")
from quaos.symplectic import PauliSum, Pauli, PauliString
from quaos.gates import GateOperation, Circuit, Hadamard as H, SUM as CX, PHASE as S
import random


def ground_state(P):
    """Returns the ground state of a given Hamiltonian

    Args:
        P: pauli, Paulis for Hamiltonian
        cc: list[int], coefficients for Hamiltonian
    
    Returns:
        numpy.array: eigenvector corresponding to lowest eigenvalue of Hamiltonian
    """

    m = P.matrix_form()

    m = m.toarray()
    val, vec = np.linalg.eig(m)
    val = np.real(val)
    vec = np.transpose(vec)

    tmp_index = val.argmin(axis=0)

    gs = vec[tmp_index]
    gs = np.transpose(gs)
    gs = gs / np.linalg.norm(gs)

    if abs(min(val) - np.transpose(np.conjugate(gs)) @ m @ gs) > 10**-10:
        print("ERROR with the GS!!!")

    return gs


def random_pauli_hamiltonian(num_paulis, qudit_dims, mode='rand', seed=None):
    """
    Generates a random Pauli Hamiltonian with the given number of Pauli operators plus their Hermitian conjugate pairs.
    
    Parameters:
        num_paulis (int): Number of Pauli operators to generate.
        qudit_dims (list): List of dimensions for each qudit.
        mode (str): 'rand' or 'uniform' - dictates form of weights in the PauliSum
    
    Returns:
        tuple: A random PauliSum
    """
    n_qudits = len(qudit_dims)
    pauli_strings = []
    coefficients = []
    
    for p in range(num_paulis):
        x_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]  # np.random.randint(qudit_dims, size=n_qudits)
        z_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]  # np.random.randint(qudit_dims, size=n_qudits)
        x_exp_H = np.zeros_like(x_exp)
        z_exp_H = np.zeros_like(z_exp)
        phase_factor = 1
        pauli_str = ''
        pauli_str_H = ''

        for j in range(len(qudit_dims)):
            r, s = x_exp[j], z_exp[j]
            pauli_str += f"x{r}z{s} "
            x_exp_H[j] = (-r) % qudit_dims[j]
            z_exp_H[j] = (-s) % qudit_dims[j]
            pauli_str_H += f"x{x_exp_H[j]}z{z_exp_H[j]} "
            
            omega = np.exp(2 * np.pi * 1j / qudit_dims[j])
            phase_factor *= omega**(r * s)
        
        pauli_strings.append(PauliString(pauli_str.strip(), dimensions=qudit_dims))
        if mode == 'rand' or mode == 'random':
            coeff = np.random.normal(0, 1) + 1j * np.random.normal(0, 1)
        elif mode == 'uniform' or mode == 'one':
            coeff = 1 + 0 * 1j
        
        if (not np.array_equal(x_exp, x_exp_H)) and (not np.array_equal(z_exp, z_exp_H)):
            # random string not Hermitian, add conjugate pair
            coefficients.append(coeff)
            coefficients.append(np.conj(coeff) * phase_factor)
            pauli_strings.append(PauliString(pauli_str_H.strip(), dimensions=qudit_dims))
        else:
            coefficients.append(coeff.real)

    rand_ham = PauliSum(pauli_strings, weights=coefficients)
    return rand_ham


def pauli_eigenvalue(index, dimension):
    """
    Computes the a-th eigenvalue of a pauli with dimension d
    """
    return np.exp(2 * np.pi * 1j * index / dimension)


def symplectic_reduction_qudit(P):
    d = P.dimensions
    q = P.n_qudits()
    P1 = P.copy()
    C = Circuit(d)
    pivots = []

    for i in range(P.n_qudits()):
        C, pivots = symplectic_reduction_iter_qudit_(P1.copy(), C, pivots, i)
    P1 = C.act(P1)

    removable_qubits = set(range(q)) - set([pivot[1] for pivot in pivots])
    conditional_qubits = sorted(set(range(q)) - removable_qubits - set([pivot[1] for pivot in pivots if pivot[2] == 'Z']))
    if any(conditional_qubits):

        for cq in conditional_qubits:
            g = H(cq, d[cq])
            C.add_gate(g)
        P1 = g.act(P1)

    return C, sorted(pivots, key=lambda x: x[1])


def number_of_SUM_X(r_control, r_target, d):
    N = 1
    while (r_target + N * r_control) % d != 0:
        if N > d:
            raise Exception('Error in Exponents r_control = ' + str(r_control) + ' r_target = ' + str(r_target))
        N += 1
        
    return N


def number_of_SUM_Z(s_control, s_target, d):
    N = 1
    while (s_control - N * s_target) % d != 0:
        if N > d:
            raise Exception('Error in Exponents s_control = ' + str(s_control) + ' s_target = ' + str(s_target))
        N += 1
        
    return N


def number_of_S(x_exp, z_exp, d):
    N = 1
    while (x_exp * N + z_exp) % d != 0:
        if N > d:
            raise Exception('Error in Exponents x_exp = ' + str(x_exp) + ' z_exp = ' + str(z_exp))
        N += 1
        
    return N


def cancel_X(pauli_sum, qudit, pauli_index, C, q_max):
    list_of_gates = []
    for i in range(qudit + 1, q_max):
        # print(pauli_sum.x_exp[pauli_index, :])
        if pauli_sum.x_exp[pauli_index, i]:
            # print('num', i, number_of_SUM_X(pauli_sum.x_exp[pauli_index, qudit],
            #                                                                                pauli_sum.x_exp[pauli_index, i],
            #                                                                                pauli_sum.dimensions[i]))
            list_of_gates += [CX(qudit, i, pauli_sum.dimensions[qudit])] * number_of_SUM_X(pauli_sum.x_exp[pauli_index, qudit],
                                                                                           pauli_sum.x_exp[pauli_index, i],
                                                                                           pauli_sum.dimensions[i])
    C.add_gate(list_of_gates)
    for g in list_of_gates:
        pauli_sum = g.act(pauli_sum)
    return pauli_sum, C


def cancel_Z(pauli_sum, qudit, pauli_index, C, q_max):

    list_of_gates = []
    list_of_gates += [H(qudit, pauli_sum.dimensions[qudit])]
    for i in range(qudit + 1, q_max):
        if pauli_sum.z_exp[pauli_index, i]:
            list_of_gates += [CX(i, qudit, pauli_sum.dimensions[qudit])] * number_of_SUM_Z(pauli_sum.z_exp[pauli_index, i],
                                                                                           pauli_sum.x_exp[pauli_index, qudit],
                                                                                           pauli_sum.dimensions[i])
    list_of_gates += [H(qudit, pauli_sum.dimensions[qudit])]
    C.add_gate(list_of_gates)
    for g in list_of_gates:
        pauli_sum = g.act(pauli_sum)
    return pauli_sum, C


def cancel_Y(pauli_sum, qudit, pauli_index, C):
    list_of_gates = [S(qudit, pauli_sum.dimensions[qudit])] * number_of_S(pauli_sum.x_exp[pauli_index, qudit],
                                                                          pauli_sum.z_exp[pauli_index, qudit],
                                                                          pauli_sum.dimensions[qudit])
    C.add_gate(list_of_gates)
    for g in list_of_gates:
        pauli_sum = g.act(pauli_sum)
    return pauli_sum, C


def cancel_pauli(P, current_qudit, pauli_index, circuit, n_q_max):
    # add CX gates to cancel out all non-zero X-parts on Pauli pauli_index, i > qudit
    if any(P.x_exp[pauli_index, i] for i in range(current_qudit + 1, n_q_max)):
        P, circuit = cancel_X(P, current_qudit, pauli_index, circuit, n_q_max)

    # add CZ gates to cancel out all non-zero Z-parts on Pauli pauli_index, i > qudit
    if any(P.z_exp[pauli_index, i] for i in range(current_qudit + 1, n_q_max)):
        P, circuit = cancel_Z(P, current_qudit, pauli_index, circuit, n_q_max)

    # if indexed Pauli, qudit is Y, add S gate to make it X
    if P.z_exp[pauli_index, current_qudit] and P.x_exp[pauli_index, current_qudit]:
        P, circuit = cancel_Y(P, current_qudit, pauli_index, circuit)
    return P, circuit


def symplectic_reduction_iter_qudit_(P, C, pivots, current_qudit):
    n_p, n_q = P.n_paulis(), P.n_qudits()
    P = C.act(P)
    n_q_max = n_q
    for i in range(n_q - current_qudit):
        if P.dimensions[current_qudit + i] != P.dimensions[current_qudit]:
            n_q_max = current_qudit + i - 1
            break

    if any(P.x_exp[:, current_qudit]) or any(P.z_exp[:, current_qudit]):
        if not any(P.x_exp[:, current_qudit]):
            g = H(current_qudit, P.dimensions[current_qudit])
            C.add_gate(g)
            P = g.act(P)

        current_pauli = min(i for i in range(n_p) if P.x_exp[i, current_qudit])  # first Pauli that has an x-component
        pivots.append((current_pauli, current_qudit, 'X'))

        P, C = cancel_pauli(P, current_qudit, current_pauli, C, n_q_max)

    if any(P.z_exp[:, current_qudit]):
        current_pauli = min(i for i in range(n_p) if P.z_exp[i, current_qudit])  # first Pauli that has a z-component
        pivots.append((current_pauli, current_qudit, 'Z'))

        g = H(current_qudit, P.dimensions[current_qudit])
        C.add_gate(g)
        P = g.act(P)

        P, C = cancel_pauli(P, current_qudit, current_pauli, C, n_q_max)  ## ### ## 

        g = H(current_qudit, P.dimensions[current_qudit])
        C.add_gate(g)
        P = g.act(P)
    
    return C, pivots


def symplectic_pauli_reduction(hamiltonian: PauliSum) -> PauliSum:
    C, pivots = symplectic_reduction_qudit(hamiltonian)
    return C, pivots


if __name__ == "__main__":
    n_paulis = 5
    n_qudits = 5
    dimension = 2

    for i in range(10):
        print(i)

        ham = random_pauli_hamiltonian(n_paulis, [dimension] * n_qudits, mode='random')
        # print(ham, '/n')

        circuit, pivots = symplectic_pauli_reduction(ham)
        reduced_hamiltonian = circuit.act(ham)
        # print(circuit)
        print(reduced_hamiltonian)
        # circuit.show()
