from __future__ import annotations

import numpy as np
import qnexus as qnx
from pytket.backends.backendresult import BackendResult
from pytket.circuit import Circuit as PytketCircuit, OpType
from qnexus.models.filters import SortFilterEnum
from qnexus.models.job_status import JobStatusEnum
from qnexus.models.references import (
    CircuitRef,
    ExecutionResultRef,
    JobType
)

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES, Gate


# https://docs.quantinuum.com/systems/user_guide/hardware_user_guide/h2.html#native-gate-set

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
    # SympleQ Y = X·Z = -i·σ_y, so Y² = -I and Y_inv = -Y. Both differ from
    # OpType.Y only by a global phase, which is invisible to measurement
    # (and stays global through any subsequent unconditional gates).
    GATES.Y: (OpType.Y, []),
    GATES.Y_inv: (OpType.Y, []),
    GATES.Z: (OpType.Z, []),
    GATES.Z_inv: (OpType.Z, []),       # Z is self-inverse for qubits
    GATES.CX: (OpType.CX, []),
    GATES.CX_inv: (OpType.CX, []),     # CX is self-inverse
    GATES.SWAP: (OpType.SWAP, []),
    GATES.CZ: (OpType.CZ, []),
    GATES.ZZMax: (OpType.ZZPhase, [0.5]),
    # ZZMax is fixed-angle (no params) so the inverse uses ZZPhase(-0.5).
    GATES.ZZMax_inv: (OpType.ZZPhase, [-0.5]),
    # V = √X = PhasedX(0.5, 0); V^{-1} = √X^{-1} = PhasedX(-0.5, 0).
    GATES.V: (OpType.PhasedX, [0.5, 0]),
    GATES.V_inv: (OpType.PhasedX, [-0.5, 0]),
}

# SympleQ gates that map 1:1 to a single H2 native gate (no rebase needed).
# Use this as a ``gates_set`` to keep the compiled circuit native.
# {GATES.S, GATES.V, GATES.ZZMax} alone already generates the Clifford group.
NATIVE_GATES_SET: list[Gate] = [
    GATES.S, GATES.S_inv,
    GATES.X, GATES.Y, GATES.Z,
    GATES.V, GATES.V_inv,
    GATES.ZZMax, GATES.ZZMax_inv,
]


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
    OpType.ZZMax: GATES.ZZMax,
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
    from pytket.passes import AutoRebase
    AutoRebase({OpType.PhasedX, OpType.Rz, OpType.ZZPhase}).apply(tk_circuit)
    # from pytket.passes import (
    #     SequencePass, SquashRzPhasedX, RemoveRedundancies, CommuteThroughMultis,
    # )
    # SequencePass([
    #     CommuteThroughMultis(),
    #     SquashRzPhasedX(),
    #     RemoveRedundancies(),
    # ]).apply(tk_circuit)
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
            # Only the Clifford angle θ=-0.5 maps back (the inverse of ZZMax).
            theta = float(command.op.params[0]) % 4.0
            if np.isclose(theta, 3.5):
                sympleq_gate: Gate = GATES.ZZMax_inv
            elif np.isclose(theta, 0.5):
                sympleq_gate = GATES.ZZMax
            else:
                raise ValueError(
                    f"Only Clifford ZZPhase(±0.5) is supported, got θ={theta} half-turns."
                )
        elif op_type == OpType.PhasedX:
            # Only the Clifford √X axis (φ=0) maps back to a sympleq gate.
            theta = float(command.op.params[0]) % 4.0
            phi = float(command.op.params[1]) % 2.0
            if not np.isclose(phi, 0.0):
                raise ValueError(
                    f"Only PhasedX with φ=0 is supported, got φ={phi} half-turns."
                )
            if np.isclose(theta, 0.5):
                sympleq_gate = GATES.V
            elif np.isclose(theta, 3.5):
                sympleq_gate = GATES.V_inv
            else:
                raise ValueError(
                    f"Only Clifford PhasedX(±0.5, 0) is supported, got θ={theta} half-turns."
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


def fetch_recent_execute_jobs(
    project_name: str,
    n: int,
    device_name: str | None = None,
) -> list[tuple[PytketCircuit, BackendResult]]:
    """Return up to ``n`` ``(circuit, result)`` pairs from the most recent execute jobs.

    Parameters
    ----------
    project_name : str
        Name of the Nexus project to query.
    n : int
        Maximum number of execute jobs to fetch.
    device_name : str | None
        If given, only jobs whose ``system.name`` equals ``device_name``
        are kept (filtered client-side, newest-first).

    Returns
    -------
    list[tuple[PytketCircuit, BackendResult]]
        Flat list of ``(circuit, result)`` pairs across the last ``n`` jobs,
        ordered newest-first.
    """
    project = qnx.projects.get(name=project_name)
    job_iter = qnx.jobs.get_all(
        project=project,
        job_type=[JobType.EXECUTE],
        sort_filters=[SortFilterEnum.CREATED_DESC],
    )

    jobs = []
    for job in job_iter:
        if device_name is not None and (job.system is None or job.system.name != device_name):
            continue
        if job.last_status != JobStatusEnum.COMPLETED:
            continue
        jobs.append(job)
        if len(jobs) >= n:
            break

    print(f"Fetched {len(jobs)} jobs.")
    pairs: list[tuple[PytketCircuit, BackendResult]] = []
    for job in jobs:
        for ref in qnx.jobs.results(job):
            if not isinstance(ref, ExecutionResultRef):
                continue
            result = ref.download_result()
            if not isinstance(result, BackendResult):
                continue
            input_program = ref.get_input()
            if not isinstance(input_program, CircuitRef):
                continue
            pairs.append((input_program.download_circuit(), result))
        print(len(pairs))
    print(f"Fetched {len(pairs)} results.")
    return pairs
