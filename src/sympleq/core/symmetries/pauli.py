import numpy as np
import itertools
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Circuit, Hadamard as H, SUM as CX, PHASE as S


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
            number_of_sum_x = number_of_SUM_X(pauli_sum.x_exp[pauli_index, qudit],
                                              pauli_sum.x_exp[pauli_index, i],
                                              pauli_sum.dimensions[i])
            list_of_gates += [CX(qudit, i, pauli_sum.dimensions[qudit])] * number_of_sum_x
    C.add_gate(list_of_gates)
    for g in list_of_gates:
        pauli_sum = g.act(pauli_sum)
    return pauli_sum, C


def cancel_Z(pauli_sum, qudit, pauli_index, C, q_max):

    list_of_gates = []
    list_of_gates += [H(qudit, pauli_sum.dimensions[qudit])]
    for i in range(qudit + 1, q_max):
        if pauli_sum.z_exp[pauli_index, i]:
            number_of_sum_z = number_of_SUM_Z(pauli_sum.z_exp[pauli_index, i],
                                              pauli_sum.x_exp[pauli_index, qudit],
                                              pauli_sum.dimensions[i])
            list_of_gates += [CX(i, qudit, pauli_sum.dimensions[qudit])] * number_of_sum_z
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
    pivot_qudits = set([pivot[1] for pivot in pivots if pivot[2] == 'Z'])
    conditional_qubits = sorted(set(range(q)) - removable_qubits - pivot_qudits)
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

    _ = len(list_of_z_symmetries)
    all_phases = [list(bits) for bits in itertools.product(list_of_phases)]
    # z symmetries can simply alter the phase of the Paulis
    conditioned_hamiltonians = []

    for sector in range(min(n_sectors, len(all_phases))):

        conditioned_hamiltonian = h_red.copy()
        phase_factor = np.zeros(h_red.n_paulis(), dtype=int)
        for i, z_symmetry in enumerate(list_of_z_symmetries):

            phase_factor[list_of_z_symmetries[i][1]] += all_phases[sector]
        conditioned_hamiltonian.set_phases((conditioned_hamiltonian.phases +
                                           phase_factor) % conditioned_hamiltonian.lcm)
        # TODO: evaluate if it is correct to set _delete_qudits as internal methods

        conditioned_hamiltonian._delete_qudits(list(z_symmetric_qudits))
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
    from sympleq.models.random_hamiltonian import random_pauli_hamiltonian, random_pauli_symmetry_hamiltonian

    n_qudits = 12
    n_paulis = 24
    dimension = 2
    ham = random_pauli_symmetry_hamiltonian(n_qudits, n_paulis, n_redundant=0, n_conditional=2)
    circuit = symplectic_pauli_reduction(ham)
    print(ham)
    h_reduced, conditioned_hams, reducing_circuit, eigenvalues = pauli_reduce(ham)
    print(h_reduced)
    print(len(conditioned_hams))

    # for h in conditioned_hams:
    #     print(h)

    # ps = ['x1z0 x1z0',
    #       'x1z0 x0z1',
    #       'x1z0 x0z0',
    #       'x1z0 x1z1'
    #       ]

    # ps = PauliSum(ps, dimensions=[2, 2], standardise=True)
    # print(ps)
    # circuit = symplectic_pauli_reduction(ps)
    # h_reduced, conditioned_hams, reducing_circuit, eigenvalues = pauli_reduce(ps)

    # print(circuit.act(ps))
    # for h in conditioned_hams:
    #     print(h)
