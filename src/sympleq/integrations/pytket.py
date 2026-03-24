from __future__ import annotations

import numpy as np
from pytket.circuit import Circuit as PytketCircuit, OpType

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES, Gate


# Map SympleQ gate singletons to pytket OpTypes (qubit-only).
_GATE_MAP: dict[Gate, OpType] = {
    GATES.H: OpType.H,
    GATES.H_inv: OpType.H,       # H is self-inverse up to global phase; pytket has no H_inv
    GATES.S: OpType.S,
    GATES.S_inv: OpType.Sdg,
    GATES.CX: OpType.CX,
    GATES.CX_inv: OpType.CX,     # CX is self-inverse
    GATES.SWAP: OpType.SWAP,
    GATES.CZ: OpType.CZ,
}

_REVERSE_MAP: dict[OpType, Gate] = {
    OpType.H: GATES.H,
    OpType.S: GATES.S,
    OpType.Sdg: GATES.S_inv,
    OpType.CX: GATES.CX,
    OpType.SWAP: GATES.SWAP,
    OpType.CZ: GATES.CZ,
}


def to_pytket_circuit(circuit: Circuit) -> PytketCircuit:
    """
    Convert a SympleQ Circuit to a pytket Circuit.

    Only qubit circuits (all dimensions equal to 2) are supported.

    Parameters
    ----------
    circuit : Circuit
        The SympleQ circuit to convert.

    Returns
    -------
    PytketCircuit
        The equivalent pytket circuit.

    Raises
    ------
    ValueError
        If any qudit dimension is not 2 (non-qubit).
    ValueError
        If a gate has no known pytket mapping.
    """
    if not np.all(circuit.dimensions == 2):
        raise ValueError(
            f"Only qubit circuits (dimension=2) are supported, "
            f"got dimensions={list(circuit.dimensions)}."
        )

    tk_circuit = PytketCircuit(circuit.n_qudits())

    for gate, qudits in zip(circuit.gates, circuit.qudit_indices):
        op_type = _GATE_MAP.get(gate)
        if op_type is None:
            raise ValueError(f"No pytket mapping for gate '{gate.name}'.")
        tk_circuit.add_gate(op_type, list(qudits))

    return tk_circuit


def from_pytket_circuit(tk_circuit: PytketCircuit) -> Circuit:
    """
    Convert a pytket Circuit to a SympleQ Circuit.

    Only a subset of pytket gates is supported (H, S, Sdg, CX, SWAP, CZ).

    Parameters
    ----------
    tk_circuit : PytketCircuit
        The pytket circuit to convert.

    Returns
    -------
    Circuit
        The equivalent SympleQ circuit.

    Raises
    ------
    ValueError
        If the pytket circuit contains an unsupported gate.
    """

    n_qubits = tk_circuit.n_qubits
    dimensions = np.array([2] * n_qubits, dtype=int)

    gates: list[Gate] = []
    qudit_indices: list[tuple[int, ...]] = []

    for command in tk_circuit.get_commands():
        op_type = command.op.type
        if op_type == OpType.Barrier:
            continue

        sympleq_gate = _REVERSE_MAP.get(op_type)
        if sympleq_gate is None:
            raise ValueError(f"No SympleQ mapping for pytket gate '{op_type}'.")

        qubits = tuple(qubit.index[0] for qubit in command.qubits)
        gates.append(sympleq_gate)
        qudit_indices.append(qubits)

    return Circuit(dimensions, gates, qudit_indices)
