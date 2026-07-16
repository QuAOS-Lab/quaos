import numpy as np
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits import Circuit, GATES, PauliGate
from sympleq.core.symmetries.conditional_hamiltonian import ConditionalHamiltonian2


def number_of_SUM_X(r_control: int, r_target: int, d: int) -> int:
    """
    Return the number of SUM gates needed to cancel out the X part of a Pauli operator.

    Parameters
        ----------
        r_control : int
            The exponent of the control qudit.
        r_target : int
            The exponent of the target qudit.
        d : int
            The dimension of the qudits.

    Returns
    -------
    int
        The number of SUM gates needed to cancel out the X part of a Pauli operator.
    """
    N = 1
    while (r_target + N * r_control) % d != 0:
        if N > d:
            raise Exception('Error in Exponents r_control = ' + str(r_control) + ' r_target = ' + str(r_target))
        N += 1

    return N


def number_of_SUM_Z(s_control: int, s_target: int, d: int) -> int:
    """
    Return the number of SUM gates needed to cancel out the Z part of a Pauli operator.

    Parameters
        ----------
        s_control : int
            The exponent of the control qudit.
        s_target : int
            The exponent of the target qudit.
        d : int
            The dimension of the qudits.

    Returns
    -------
    int
        The number of SUM gates needed to cancel out the Z part of a Pauli operator.
    """
    N = 1
    while (s_control - N * s_target) % d != 0:
        if N > d:
            raise Exception('Error in Exponents s_control = ' + str(s_control) + ' s_target = ' + str(s_target))
        N += 1

    return N


def number_of_S(x_exp: int, z_exp: int, d: int) -> int:
    """
    Return the number of PHASE gates needed to cancel out the Z part of a Pauli operator.

    Parameters
        ----------
        x_exp : int
            The exponent of the X part of the Pauli operator.
        z_exp : int
            The exponent of the Z part of the Pauli operator.
        d : int
            The dimension of the qudits.

    Returns
    -------
    int
        The number of PHASE gates needed to cancel out the Z part of a Pauli operator.
    """
    N = 1
    while (x_exp * N + z_exp) % d != 0:
        if N > d:
            raise Exception('Error in Exponents x_exp = ' + str(x_exp) + ' z_exp = ' + str(z_exp))
        N += 1

    return N


def cancel_X(pauli_sum: PauliSum, qudit: int, pauli_index: int, C: Circuit, q_max: int) -> tuple[PauliSum, Circuit]:
    """
    Cancel out the X part of a Pauli operator.

    Parameters
    ----------
    pauli_sum : PauliSum
        The PauliSum that we are cancelling out the X part of.
    qudit : int
        The qudit that we are currently operating on.
    pauli_index : int
        The index of the Pauli operator that we are cancelling out the X part of.
    C : Circuit
        The circuit that we are adding the gates to.
    q_max : int
        The maximum number of qudits.

    Returns
    -------
    pauli_sum : PauliSum
        The PauliSum after the X part of the Pauli operator has been cancelled out.
    C : Circuit
        The circuit after the gates have been added to cancel out the X part of the Pauli operator.
    """
    for i in range(qudit + 1, q_max):
        if pauli_sum.x_exp[pauli_index, i]:
            number_of_sum_x = number_of_SUM_X(pauli_sum.x_exp[pauli_index, qudit],
                                              pauli_sum.x_exp[pauli_index, i],
                                              pauli_sum.dimensions[i])
            for _ in range(number_of_sum_x):
                C.add_gate(GATES.CX, qudit, i)
                pauli_sum = GATES.CX.act(pauli_sum, (qudit, i))
    return pauli_sum, C


def cancel_Z(pauli_sum: PauliSum, qudit: int, pauli_index: int, C: Circuit, q_max: int) -> tuple[PauliSum, Circuit]:
    """
    Cancel out the Z part of a Pauli operator.

    Parameters
    ----------
    pauli_sum : PauliSum
        The PauliSum that we are cancelling out the Z part of.
    qudit : int
        The qudit that we are currently operating on.
    pauli_index : int
        The index of the Pauli operator that we are cancelling out the Z part of.
    C : Circuit
        The circuit that we are adding the gates to.
    q_max : int
        The maximum number of qudits.

    Returns
    -------
    pauli_sum : PauliSum
        The PauliSum after the Z part of the Pauli operator has been cancelled out.
    C : Circuit
        The circuit after the gates have been added to cancel out the Z part of the Pauli operator.
    """
    C.add_gate(GATES.H, qudit)
    pauli_sum = GATES.H.act(pauli_sum, qudit)
    for i in range(qudit + 1, q_max):
        if pauli_sum.z_exp[pauli_index, i]:
            number_of_sum_z = number_of_SUM_Z(pauli_sum.z_exp[pauli_index, i],
                                              pauli_sum.z_exp[pauli_index, qudit],
                                              pauli_sum.dimensions[i])
            for _ in range(number_of_sum_z):
                C.add_gate(GATES.CX, i, qudit)
                pauli_sum = GATES.CX.act(pauli_sum, (i, qudit))
    C.add_gate(GATES.H, qudit)
    pauli_sum = GATES.H.act(pauli_sum, qudit)
    return pauli_sum, C


def cancel_Y(pauli_sum: PauliSum, qudit: int, pauli_index: int, C: Circuit) -> tuple[PauliSum, Circuit]:
    """
    Cancel out the Y part of a Pauli operator.

    Parameters
    ----------
    pauli_sum : PauliSum
        The PauliSum that we are cancelling out the Y part of.
    qudit : int
        The qudit that we are currently operating on.
    pauli_index : int
        The index of the Pauli operator that we are cancelling out the Y part of.
    C : Circuit
        The circuit that we are adding the gates to.

    Returns
    -------
    pauli_sum : PauliSum
        The PauliSum after the Y part of the Pauli operator has been cancelled out.
    C : Circuit
        The circuit after the gates have been added to cancel out the Y part of the Pauli operator.
    """
    number_of_phase = number_of_S(pauli_sum.x_exp[pauli_index, qudit], pauli_sum.z_exp[pauli_index, qudit],
                                  pauli_sum.dimensions[qudit])
    for _ in range(number_of_phase):
        C.add_gate(GATES.S, qudit)
        pauli_sum = GATES.S.act(pauli_sum, qudit)
    return pauli_sum, C


def cancel_pauli(P: PauliSum, current_qudit: int,
                 pauli_index: int, circuit: Circuit, n_q_max: int) -> tuple[PauliSum, Circuit]:
    """
    Cancel out all non-zero X and Z parts of a Pauli operator.

    Parameters
    ----------
    P : PauliSum
        The PauliSum that we are cancelling out the Pauli operator from.
    current_qudit : int
        The qudit that we are currently operating on.
    pauli_index : int
        The index of the Pauli operator that we are cancelling out.
    circuit : Circuit
        The circuit that we are adding the gates to.
    n_q_max : int
        The maximum number of qudits.

    Returns
    -------
    P : PauliSum
        The PauliSum after the Pauli operator has been cancelled out.
    circuit : Circuit
        The circuit after the gates have been added to cancel out the Pauli operator.
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


def symplectic_reduction_qudit(P) -> tuple[Circuit, list]:
    """
    Applies the symplectic reduction algorithm to a PauliSum.

    This algorithm will reduce the number of qudits in the PauliSum by
    cancelling out any non-zero X and Z parts of the Pauli operators
    that are not part of the symplectic group.

    Parameters
    ----------
    P : PauliSum
        The PauliSum that we are applying the symplectic reduction algorithm to.

    Returns
    -------
    C : Circuit
        The circuit that implements the symplectic reduction algorithm.
    pivots : list
        A list of the pivots of the symplectic reduction algorithm.
    """

    d = P.dimensions
    q = P.n_qudits()
    P1 = P.copy()
    C = Circuit.empty(d)
    pivots = []

    for i in range(P.n_qudits()):
        C, pivots = symplectic_reduction_iter_qudit_(P1.copy(), C, pivots, i)
    P1 = C.act(P1)

    removable_qubits = set(range(q)) - set([pivot[1] for pivot in pivots])
    pivot_qudits = set([pivot[1] for pivot in pivots if pivot[2] == 'Z'])
    conditional_qubits = sorted(set(range(q)) - removable_qubits - pivot_qudits)
    if len(conditional_qubits) > 0:
        for cq in conditional_qubits:
            C.add_gate(GATES.H, cq)
            P1 = GATES.H.act(P1, (cq,))
    return C, sorted(pivots, key=lambda x: x[1])


def symplectic_pauli_reduction(P) -> tuple[Circuit, list]:
    """
    Backward-compatible alias for the qudit symplectic reduction routine.
    """
    return symplectic_reduction_qudit(P)


def symplectic_reduction_iter_qudit_(P, C, pivots, current_qudit) -> tuple[Circuit, list]:
    """
    Applies one iteration of the symplectic reduction algorithm to a PauliSum.

    Parameters
    ----------
    P : PauliSum
        The PauliSum that we are applying the symplectic reduction algorithm to.
    C : Circuit
        The circuit that implements the symplectic reduction algorithm.
    pivots : list
        A list of the pivots of the symplectic reduction algorithm.
    current_qudit : int
        The qudit that we are currently operating on.

    Returns
    -------
    C : Circuit
        The circuit after one iteration of the symplectic reduction algorithm.
    pivots : list
        The list of pivots after one iteration of the symplectic reduction algorithm.
    """
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
            C.add_gate(GATES.H, current_qudit)
            P = GATES.H.act(P, (current_qudit,))

        current_pauli = min(i for i in range(n_p) if P.x_exp[i, current_qudit])  # first Pauli that has an x-component
        pivots.append((current_pauli, current_qudit, 'X'))

        P, C = cancel_pauli(P, current_qudit, current_pauli, C, n_q_max)

    # If there was previously a y we need to cancel the left over z parts
    if any(P.z_exp[:, current_qudit]):
        current_pauli = min(i for i in range(n_p) if P.z_exp[i, current_qudit])  # first Pauli that has a z-component
        pivots.append((current_pauli, current_qudit, 'Z'))

        C.add_gate(GATES.H, current_qudit)
        P = GATES.H.act(P, (current_qudit,))

        P, C = cancel_pauli(P, current_qudit, current_pauli, C, n_q_max)

        C.add_gate(GATES.H, current_qudit)
        P = GATES.H.act(P, (current_qudit,))
    return C, pivots


def pauli_reduce(hamiltonian: PauliSum) -> tuple[PauliSum, list[PauliSum], Circuit, list]:
    """
    Applies the symplectic reduction algorithm to a PauliSum and returns the reduced hamiltonian and the
    conditioned hamiltonians.

    Parameters
    ----------
    hamiltonian : PauliSum
        The PauliSum that we are applying the symplectic reduction algorithm to.

    Returns
    -------
    h_red : PauliSum
        The PauliSum after applying the symplectic reduction algorithm.
    conditioned_hamiltonians : list[PauliSum]
        The list of PauliSum's that are the conditioned hamiltonians.
    C : Circuit
        The Circuit that implements the symplectic reduction algorithm.
    all_phases : list
        The list of all possible phases of the conditioned hamiltonians.
    """
    C, _ = symplectic_reduction_qudit(hamiltonian)

    h_red = C.act(hamiltonian)
    # first we remove any qudits with only identities
    h_red.remove_trivial_qudits()

    # build list of z symmetries as those qubits with only z
    z_symmetric_qudits = set()
    for i in range(h_red.n_qudits()):
        if not any(h_red.x_exp[:, i]):  # z only
            z_symmetric_qudits.add(i)

    # build the conditional hamiltonian for the z symmetries
    z_exp = np.array([1 if i in z_symmetric_qudits else 0 for i in range(h_red.n_qudits())])
    Sym = PauliString.from_exponents(np.zeros(h_red.n_qudits()), z_exp, h_red.dimensions)
    S = PauliGate(Sym)
    T = Circuit.empty(h_red.dimensions).composite_gate()
    conditional_hamiltonian = ConditionalHamiltonian2(h_red, S, T)

    return conditional_hamiltonian
