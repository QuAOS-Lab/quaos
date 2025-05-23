from quaos.gates import SUM as CX, PHASE as S, Hadamard as H
from quaos.circuits import Circuit
from quaos.paulis import PauliString, PauliSum
from quaos.circuits.utils import solve_modular_linear


def add_phase(xz_pauli_sum: PauliSum, qudit_index: int, qudit_index_2: int, phase_key: str) -> Circuit:
    """
    Uses phase key to alter phase of a pauli sum based on the qudit index 2 pauli.

    Acts like identity on qudit index 1.
    Assumes pauli_sum has the form:
    qudit_index_1  | qudit_index_2
          X        |      *
          Z        |      *

    key = S -> same phase D -> different phase
    order IXZY

    e.g 'SDSD' keeps the same phase for all but X and Y paulis on qudit index 2.

    """

    if xz_pauli_sum.dimensions[qudit_index] != 2 or xz_pauli_sum.dimensions[qudit_index] != 2:
        raise ValueError("Pauli dimensions must be equal to 2")

    if phase_key == 'SSSS':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())])
        return C
    elif phase_key == 'DDSS':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())],
                    gates=[CX(qudit_index, qudit_index_2, 2), S(qudit_index_2, 2), H(qudit_index_2, 2),
                           S(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2)])
        return C
    elif phase_key == 'DSDS':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())],
                    gates=[CX(qudit_index, qudit_index_2, 2), S(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2),
                           H(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2), S(qudit_index, 2)])
        return C
    elif phase_key == 'DSSD':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())],
                    gates=[S(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2), S(qudit_index_2, 2),
                           H(qudit_index_2, 2), S(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2)])
        return C
    elif phase_key == 'SDDS':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())],
                    gates=[S(qudit_index_2, 2), H(qudit_index_2, 2), CX(qudit_index_2, qudit_index, 2),
                           S(qudit_index, 2), CX(qudit_index_2, qudit_index, 2), H(qudit_index_2, 2),
                           CX(qudit_index, qudit_index_2, 2), S(qudit_index, 2)])
        return C
    elif phase_key == 'SDSD':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())],
                    gates=[CX(qudit_index_2, qudit_index, 2), S(qudit_index, 2), CX(qudit_index_2, qudit_index, 2),
                           H(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2), S(qudit_index, 2)])
        return C
    elif phase_key == 'SSDD':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())],
                    gates=[H(qudit_index_2, 2), CX(qudit_index_2, qudit_index, 2), S(qudit_index, 2),
                           CX(qudit_index_2, qudit_index, 2), H(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2),
                           S(qudit_index, 2)])
        return C
    elif phase_key == 'DDDD':
        C = Circuit(dimensions=[2 for i in range(xz_pauli_sum.n_qudits())],
                    gates=[CX(qudit_index_2, qudit_index, 2), S(qudit_index, 2), CX(qudit_index_2, qudit_index, 2),
                           H(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2), S(qudit_index, 2),
                           CX(qudit_index, qudit_index_2, 2), S(qudit_index_2, 2), H(qudit_index_2, 2),
                           S(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2)])
        return C
    else:
        raise ValueError("Invalid phase key. Must be one of 'SSSS', 'DDSS', 'DSDS', 'DSSD', 'SDDS', 'SDSD', 'SSDD', 'DDDD'")


def add_s2(pauli_sum: PauliSum, qudit_index_1: int, qudit_index_2: int) -> Circuit:
    """
    xr1zs1 xr2zs2 -> xr1+s2 zs1+s2  *
    """
    C = Circuit(dimensions=[2 for i in range(pauli_sum.n_qudits())],
                gates=[CX(qudit_index_1, qudit_index_2, 2), H(qudit_index_2, 2), CX(qudit_index_2, qudit_index_1, 2),
                       S(qudit_index_2, 2), H(qudit_index_2, 2)])
    return C


def add_r2(pauli_sum: PauliSum, qudit_index_1: int, qudit_index_2: int) -> Circuit:
    """
    xr1zs1 xr2zs2 -> xr1+r2 zs1+r2  *
    """
    C = Circuit(dimensions=[2 for i in range(pauli_sum.n_qudits())],
                gates=[S(qudit_index_1, 2), CX(qudit_index_2, qudit_index_1, 2), S(qudit_index_1, 2)])
    return C


def add_r2s2(pauli_sum: PauliSum, qudit_index_1: int, qudit_index_2: int) -> Circuit:
    """
    xr1zs1 xr2zs2 -> xr1+r2+s2 zs1+r2+s2  *
    """
    C = Circuit(dimensions=[2 for i in range(pauli_sum.n_qudits())],
                gates=[S(qudit_index_2, 2), CX(qudit_index_1, qudit_index_2, 2), H(qudit_index_2, 2),
                       CX(qudit_index_2, qudit_index_1, 2), S(qudit_index_2, 2), H(qudit_index_2, 2)])
    return C


def ensure_zx_components(pauli_sum: PauliSum, pauli_index_x: int, pauli_index_z: int, target_qubit: int) -> tuple[Circuit, PauliSum]:
    """
    Assumes anti-commutation between pauli_index_x and pauli_index_z.
    brings pauli_sum to the form:

                    target_qubit
    pauli_index_x |  xr1zs1
    pauli_index_z |  xr2zs2
    where r1 and s2 are always non-zero

    """

    if not pauli_sum[pauli_index_x, target_qubit:].commute(pauli_sum[pauli_index_z, target_qubit:]):
        raise ValueError("ensure_zx_components requires anti-commutation between pauli_index_x and pauli_index_z beyond target_qubit")

    C = Circuit(dimensions=pauli_sum.dimensions)
    # prepare anti-commuting pauli strings with the same absolute coefficients for test of hadamard Symmetry
    # prime pauli pi and pj for cancel_pauli
    if pauli_sum.x_exp[pauli_index_x, target_qubit] == 1 and pauli_sum.z_exp[pauli_index_z, target_qubit] == 1:  # x,z
        px = pauli_index_x
    elif pauli_sum.z_exp[pauli_index_x, target_qubit] == 1 and pauli_sum.x_exp[pauli_index_z, target_qubit] == 1:  # z,x
        px = pauli_index_z
    elif pauli_sum.x_exp[pauli_index_x, target_qubit] == 1 and pauli_sum.z_exp[pauli_index_z, target_qubit] == 0:  # x,id or x,x
        if any(pauli_sum.z_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_z, i]]), 2)
        elif any(pauli_sum.x_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = H(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]), 2)
            pauli_sum = g.act(pauli_sum)
            C.add_gate(g)
            g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_z, i]]), 2)
        C.add_gate(g)
        pauli_sum = g.act(pauli_sum)
        px = pauli_index_x
    elif pauli_sum.z_exp[pauli_index_x, target_qubit] == 1 and pauli_sum.x_exp[pauli_index_z, target_qubit] == 0:  # z,id or z,z
        if any(pauli_sum.x_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]), target_qubit, 2)
        elif any(pauli_sum.z_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = H(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_z, i]]), 2)
            pauli_sum = g.act(pauli_sum)
            C.add_gate(g)
            g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]), target_qubit, 2)
        C.add_gate(g)
        pauli_sum = g.act(pauli_sum)
        px = pauli_index_z
    elif pauli_sum.x_exp[pauli_index_x, target_qubit] == 0 and pauli_sum.z_exp[pauli_index_z, target_qubit] == 1:  # id,z
        if any(pauli_sum.x_exp[pauli_index_x, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_x, i]]), target_qubit, 2)
        elif any(pauli_sum.z_exp[pauli_index_x, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = H(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_x, i]]), 2)
            pauli_sum = g.act(pauli_sum)
            C.add_gate(g)
            g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_x, i]]), target_qubit, 2)
        C.add_gate(g)
        pauli_sum = g.act(pauli_sum)
        px = pauli_index_x
    elif pauli_sum.x_exp[pauli_index_x, target_qubit] == 0 and pauli_sum.x_exp[pauli_index_z, target_qubit] == 1:   # id,x
        if any(pauli_sum.z_exp[pauli_index_x, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_x, i]]), 2)
        elif any(pauli_sum.x_exp[pauli_index_x, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = H(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_x, i]]), 2)
            pauli_sum = g.act(pauli_sum)
            C.add_gate(g)
            g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_x, i]]), 2)
        C.add_gate(g)
        pauli_sum = g.act(pauli_sum)
        px = pauli_index_z
    else:  # id,id
        if any(pauli_sum.x_exp[pauli_index_x, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_x, i]]), target_qubit, 2)
            pauli_sum = g.act(pauli_sum)
            C.add_gate(g)
            if any(pauli_sum.z_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
                g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_z, i]]), 2)
            elif any(pauli_sum.x_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
                g = H(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]), 2)
                pauli_sum = g.act(pauli_sum)
                C.add_gate(g)
                g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_z, i]]), 2)
            C.add_gate(g)
            pauli_sum = g.act(pauli_sum)
            px = pauli_index_x
        elif any(pauli_sum.z_exp[pauli_index_x, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_x, i]]), 2)
            pauli_sum = g.act(pauli_sum)
            C.add_gate(g)
            if any(pauli_sum.x_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
                g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]), target_qubit, 2)
            elif any(pauli_sum.z_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
                g = H(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_z, i]]), 2)
                pauli_sum = g.act(pauli_sum)
                C.add_gate(g)
                g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]), target_qubit, 2)
            C.add_gate(g)
            pauli_sum = g.act(pauli_sum)
            px = pauli_index_z
    if px == pauli_index_z:
        C.add_gate(H(target_qubit, 2))
        pauli_sum = H(target_qubit, 2).act(pauli_sum)

    return C, pauli_sum


def to_ix(pauli_string: PauliString, target_index: int, ignore: int | list[int] | None = None) -> Circuit:
    """Finds a circuit to turn a PauliString to III...IXI...II where the X is at target_index
    
    for qudits X = xrz0 for any r != 0, I = x0z0

    ignore is a list of qudits to not try to turn into an I
    """

    if ignore is None:
        ignore = []
    if isinstance(ignore, int):
        ignore = [ignore]
    if target_index in ignore:
        raise Exception("target_index must not be in ignore")

    p_string_in = pauli_string.copy()
    # First we make the target qudit an X
    circuit = to_x(pauli_string, target_index, ignore)
    pauli_string = circuit.act(p_string_in)
    n_q = pauli_string.n_qudits()
    for q in range(n_q):
        if q != target_index and q not in ignore:
            if pauli_string[q].x_exp == 0 and pauli_string[q].z_exp != 0:
                circuit.add_gate(H(q, pauli_string.dimensions[q]))
                pauli_string = circuit.act(p_string_in)
            if pauli_string[q].x_exp != 0:
                if pauli_string[q].z_exp != 0:
                    # use the x to cancel the z of q with S gates
                    n_s = solve_modular_linear(pauli_string[q].x_exp, pauli_string[q].z_exp, pauli_string.dimensions[q])
                    for i in range(n_s):
                        circuit.add_gate(S(q, pauli_string.dimensions[q]))
                    pauli_string = circuit.act(p_string_in)

                # use cnot to cancel the x of q with the x of target. n_cnot = n where x_q + x+_target) % d= 0
                n_cnot = solve_modular_linear(pauli_string[q].x_exp, pauli_string[target_index].x_exp, pauli_string.dimensions[target_index])
                if n_cnot is None:
                    raise Exception("Weird")
                for i in range(n_cnot):
                    circuit.add_gate(CX(target_index, q, pauli_string.dimensions[target_index]))
                pauli_string = circuit.act(p_string_in)

    else:
        return circuit


def to_x(pauli_string: PauliString, target_index: int, ignore: int | list[int] | None = None) -> Circuit:
    """Finds a circuit to turn a PauliString to ****X*** where the X is at target_index"""
    if ignore is None:
        ignore = []
    if isinstance(ignore, int):
        ignore = [ignore]
    if target_index in ignore:
        raise Exception("target_index must not be in ignore")

    if target_index < 0:
        target_index += pauli_string.n_qudits()
    if target_index > pauli_string.n_qudits():
        raise Exception(f"target_index {target_index} out of range {pauli_string.n_qudits()}")

    n_q = pauli_string.n_qudits()
    if PauliString.n_identities == n_q:
        raise Exception("PauliString is identity - cannot be converted to X")

    circuit = Circuit(dimensions=pauli_string.dimensions)  # Empty circuit of correct dimensions
    
    # First we check if we can get there from a single qudit gate
    # target is X
    dim_target = pauli_string.dimensions[target_index]
    if pauli_string[target_index].x_exp != 0 and pauli_string[target_index].z_exp == 0:  # already x
        return circuit
    elif pauli_string[target_index].x_exp != 0 and pauli_string[target_index].z_exp != 0:
        # We remove the Z by performing S gates
        n_s_gates = solve_modular_linear(pauli_string[target_index].z_exp, pauli_string[target_index].x_exp, dim_target)
        for i in range(n_s_gates):
            circuit.add_gate(S(target_index, dim_target))  # xrzs -> xrz(s+r) until s+r = 0
        return circuit
    elif pauli_string[target_index].x_exp == 0 and pauli_string[target_index].z_exp != 0:
        circuit.add_gate(H(target_index, dim_target))  # swap x and z
        return circuit
    elif pauli_string[target_index].x_exp == 0 and pauli_string[target_index].z_exp == 0:  # must be x0z0
        # We have to use a multi-qubit gate - we loop through the qudits until we find a resource
        for q in range(n_q - 1, -1, -1):
            if q != target_index and q not in ignore:
                if pauli_string[q].x_exp != 0:
                    # move x of q to target
                    circuit.add_gate(CX(q, target_index, dim_target))
                    return circuit
                if pauli_string[q].x_exp == 0 and pauli_string[q].z_exp != 0:
                    # move z of q to target_index then Hadamard
                    circuit.add_gate(CX(target_index, q, dim_target))
                    circuit.add_gate(H(target_index, dim_target))
                    return circuit
    else:
        raise Exception("A case has been missed from the above loop")
    raise Exception(f"No circuit found to convert {pauli_string} pauli to X")
    
