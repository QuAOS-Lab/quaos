from __future__ import annotations

import numpy as np
from pytket.circuit import Circuit as PytketCircuit, OpType

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES, Gate


# Map SympleQ gate singletons to (pytket OpType, params) pairs (qubit-only).
# `params` are angle parameters in half-turns for parameterized gates.
_GATE_MAP: dict[Gate, tuple[OpType, list[float]]] = {
    GATES.Id: (OpType.noop, []),
    GATES.H: (OpType.H, []),
    GATES.H_inv: (OpType.H, []),       # H is self-inverse for qubits
    GATES.S: (OpType.S, []),
    GATES.S_inv: (OpType.Sdg, []),
    GATES.X: (OpType.X, []),
    GATES.X_inv: (OpType.X, []),       # X is self-inverse for qubits
    GATES.Y: (OpType.Y, []),
    GATES.Y_inv: (OpType.Y, []),       # Y is self-inverse for qubits
    GATES.Z: (OpType.Z, []),
    GATES.Z_inv: (OpType.Z, []),       # Z is self-inverse for qubits
    GATES.CX: (OpType.CX, []),
    GATES.CX_inv: (OpType.CX, []),     # CX is self-inverse
    GATES.SWAP: (OpType.SWAP, []),
    GATES.CZ: (OpType.CZ, []),
    # ZZPhase(α) = exp(-i α π/2 Z⊗Z); α=±0.5 corresponds to our Clifford ZZPhase(±π/4).
    GATES.ZZPhase: (OpType.ZZPhase, [0.5]),
    GATES.ZZPhase_inv: (OpType.ZZPhase, [-0.5]),
}

_REVERSE_MAP: dict[OpType, Gate] = {
    OpType.noop: GATES.Id,
    OpType.H: GATES.H,
    OpType.S: GATES.S,
    OpType.Sdg: GATES.S_inv,
    OpType.X: GATES.X,
    OpType.Y: GATES.Y,
    OpType.Z: GATES.Z,
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
        mapping = _GATE_MAP.get(gate)
        if mapping is None:
            raise ValueError(f"No pytket mapping for gate '{gate.name}'.")
        op_type, params = mapping
        if params:
            tk_circuit.add_gate(op_type, params, list(qudits))
        else:
            tk_circuit.add_gate(op_type, list(qudits))

    return tk_circuit


def from_pytket_circuit(tk_circuit: PytketCircuit) -> Circuit:
    """
    Convert a pytket Circuit to a SympleQ Circuit.

    Only a subset of pytket gates is supported (H, S, Sdg, X, Y, Z, CX, SWAP, CZ, noop).

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

        if op_type == OpType.ZZPhase:
            # ZZPhase(α) has period 4 in half-turns; only α ≡ ±0.5 (mod 4) is the
            # Clifford ZZPhase(±π/4) we support.
            angle = float(command.op.params[0]) % 4.0
            if np.isclose(angle, 0.5):
                sympleq_gate: Gate = GATES.ZZPhase
            elif np.isclose(angle, 3.5):
                sympleq_gate = GATES.ZZPhase_inv
            else:
                raise ValueError(
                    f"Only Clifford ZZPhase(±π/4) is supported, got angle={angle} half-turns."
                )
        else:
            mapped = _REVERSE_MAP.get(op_type)
            if mapped is None:
                raise ValueError(f"No SympleQ mapping for pytket gate '{op_type}'.")
            sympleq_gate = mapped

        qubits = tuple(qubit.index[0] for qubit in command.qubits)
        gates.append(sympleq_gate)
        qudit_indices.append(qubits)

    return Circuit(dimensions, gates, qudit_indices)
