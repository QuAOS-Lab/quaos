from gates import Circuit, SUM as CX, PHASE as S, Hadamard as H, GateOperation
from symplectic import Pauli, PauliString, PauliSum, symplectic_product
import numpy as np


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
                    gates=[CX(qudit_index,qudit_index_2, 2), S(qudit_index_2, 2), CX(qudit_index, qudit_index_2, 2),
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
        pz = pauli_index_z
    elif pauli_sum.z_exp[pauli_index_x, target_qubit] == 1 and pauli_sum.x_exp[pauli_index_z, target_qubit] == 1:  # z,x
        px = pauli_index_z
        pz = pauli_index_x
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
        pz = pauli_index_z
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
        pz = pauli_index_x
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
        pz = pauli_index_z
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
        pz = pauli_index_x
    else: # id,id
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
            pz = pauli_index_z
        elif any(pauli_sum.z_exp[pauli_index_x, i] for i in range(target_qubit, pauli_sum.n_qudits())):
            g = CX(target_qubit, min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_x, i]]),2)
            pauli_sum = g.act(pauli_sum)
            C.add_gate(g)
            if any(pauli_sum.x_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
                g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]), target_qubit, 2)
            elif any(pauli_sum.z_exp[pauli_index_z, i] for i in range(target_qubit, pauli_sum.n_qudits())):
                g = H(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.z_exp[pauli_index_z, i]]), 2)
                pauli_sum = g.act(pauli_sum)
                C.add_gate(g)
                g = CX(min([i for i in range(target_qubit, pauli_sum.n_qudits()) if pauli_sum.x_exp[pauli_index_z, i]]),target_qubit, 2)
            C.add_gate(g)
            pauli_sum = g.act(pauli_sum)
            px = pauli_index_z
    if px == pauli_index_z:
        C.add_gate(H(target_qubit, 2))
        pauli_sum = H(target_qubit, 2).act(pauli_sum)

    return C, pauli_sum
