import numpy as np
import itertools
import random
import sys
sys.path.append("./")
from quaos.paulis import PauliSum, PauliString
from quaos.gates import Circuit, Hadamard as H, SUM as CX, PHASE as S


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


def number_of_SUM_X(r_control, r_target, d):
    """
    Find the smallest positive integer N such that:
        (r_target + N * r_control) % d == 0

    This counts the number N of SUM gates needed to cancel out the X part of a Pauli operator.
    """
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
        if pauli_sum.x_exp[pauli_index, i]:
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
    """
    Needs an x component on current_qudit

    P -> p_1 ... p_current_qudit  I I ... I p_n_q_max p.... p_n_paulis
    
    """
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
    if len(conditional_qubits) > 0:
        for cq in conditional_qubits:
            g = H(cq, d[cq])
            C.add_gate(g)
        P1 = g.act(P1)
    return C, sorted(pivots, key=lambda x: x[1])


def symplectic_reduction_iter_qudit_(P, C, pivots, current_qudit):
    n_p, n_q = P.n_paulis(), P.n_qudits()
    P = C.act(P)
    n_q_max = n_q
    # find n_q_max, the last qudit of the same dimension as current_qudit
    for i in range(n_q - current_qudit):
        if P.dimensions[current_qudit + i] != P.dimensions[current_qudit]:
            n_q_max = current_qudit + i - 1
            break

    # does the current qudit have any X or Z components?
    if any(P.x_exp[:, current_qudit]) or any(P.z_exp[:, current_qudit]):
        if not any(P.x_exp[:, current_qudit]):  # If it is z we need to add a Hadamard gate to make it an X
            g = H(current_qudit, P.dimensions[current_qudit])
            C.add_gate(g)
            P = g.act(P)

        current_pauli = min(i for i in range(n_p) if P.x_exp[i, current_qudit])  # first Pauli that has an x-component
        pivots.append((current_pauli, current_qudit, 'X'))

        P, C = cancel_pauli(P, current_qudit, current_pauli, C, n_q_max)

    # If there was previously a y we need to cancel the left over z parts
    if any(P.z_exp[:, current_qudit]):
        current_pauli = min(i for i in range(n_p) if P.z_exp[i, current_qudit])  # first Pauli that has a z-component
        pivots.append((current_pauli, current_qudit, 'Z'))

        g = H(current_qudit, P.dimensions[current_qudit])
        C.add_gate(g)
        P = g.act(P)

        P, C = cancel_pauli(P, current_qudit, current_pauli, C, n_q_max)

        g = H(current_qudit, P.dimensions[current_qudit])
        C.add_gate(g)
        P = g.act(P)
    return C, pivots


def symplectic_pauli_reduction(hamiltonian: PauliSum) -> Circuit:
    C, pivots = symplectic_reduction_qudit(hamiltonian)
    return C


def pauli_reduce(hamiltonian: PauliSum) -> tuple[PauliSum, list[PauliSum], Circuit, list]:
    """
    Reduces the Hamiltonian to a smaller number of qudits by removing leading X and Z operators.
    
    This returns a list of reduced Hamiltonians, each corresponding to a different symmetry sector of the Z symmetries.

    """
    # hamiltonian.remove_trivial_qudits()
    C = symplectic_pauli_reduction(hamiltonian)

    h_red = C.act(hamiltonian)
    # first we remove any qudits with only identities
    h_red.remove_trivial_qudits()

    # build list of z symmetries as those qubits with only z
    list_of_z_symmetries = []
    list_of_phases = []
    n_sectors = 1
    z_symmetric_qudits = set()
    for i in range(h_red.n_qudits()):
        if not any(h_red.x_exp[:, i]):  # z only
            list_of_z_symmetries.append((i, np.where(h_red.z_exp[:, i] != 0)[0]))
            z_symmetric_qudits.add(i)
            list_of_phases += np.arange(h_red.dimensions[i]).tolist()
            n_sectors *= h_red.dimensions[i]

    n_symmetries = len(list_of_z_symmetries)
    all_phases = [list(bits) for bits in itertools.product(list_of_phases)]
    # z symmetries can simply alter the phase of the Paulis
    conditioned_hamiltonians = []

    for sector in range(min(n_sectors, len(all_phases))):

        conditioned_hamiltonian = h_red.copy()
        phase_factor = np.zeros(h_red.n_paulis(), dtype=int)
        for i, z_symmetry in enumerate(list_of_z_symmetries):

            phase_factor[list_of_z_symmetries[i][1]] += all_phases[sector]
        conditioned_hamiltonian.phases += phase_factor
        conditioned_hamiltonian.phases = conditioned_hamiltonian.phases % conditioned_hamiltonian.lcm
        conditioned_hamiltonian.delete_qudits(list(z_symmetric_qudits))
        conditioned_hamiltonian.combine_equivalent_paulis()
        conditioned_hamiltonians.append(conditioned_hamiltonian)
    return h_red, conditioned_hamiltonians, C, all_phases


if __name__ == "__main__":

    # Example from the paper

    # ham = ['x1z0 x1z0 x0z0', 'x0z0 x1z0 x1z0', 'x0z1 x0z0 x0z1']
    # ham = PauliSum(ham, weights=[1, 1, 1], dimensions=[2, 2, 2])
    # print(ham, '/n')
    # circuit = symplectic_pauli_reduction(ham)
    # h_reduced, conditioned_hams, reducing_circuit, eigenvalues = pauli_reduce(ham)
    # print(h_reduced)
    # for h in conditioned_hams:
    #     print(h)

    # random hamiltonian example

    # n_qudits = 10
    # n_paulis = 10
    # dimension = 2
    # ham = random_pauli_hamiltonian(n_paulis, [dimension] * n_qudits, mode='uniform')
    # circuit = symplectic_pauli_reduction(ham)
    # print(ham)
    # h_reduced, conditioned_hams, reducing_circuit, eigenvalues = pauli_reduce(ham)
    # print(h_reduced)
    # for h in conditioned_hams:
    #     print(h)


    ps = ['x1z0 x1z0',
          'x1z0 x0z1',
          'x1z0 x0z0',
          'x1z0 x1z1'
          ]
    
    ps = PauliSum(ps, dimensions=[2, 2], standardise=True)
    print(ps)
    circuit = symplectic_pauli_reduction(ps)
    h_reduced, conditioned_hams, reducing_circuit, eigenvalues = pauli_reduce(ps)

    print(circuit.act(ps))
    for h in conditioned_hams:
        print(h)
