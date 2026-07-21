from __future__ import annotations

from functools import lru_cache
from typing import Iterator

import numpy as np
import qnexus as qnx
from pytket.backends.backendresult import BackendResult
from pytket.circuit import Circuit as PytketCircuit, OpType
from pytket.extensions.quantinuum.backends.quantinuum import QuantinuumBackend
from pytket.extensions.quantinuum.backends.api_wrappers import QuantinuumAPI
from pytket.passes import BasePass
from qnexus.models.filters import SortFilterEnum
from qnexus.models.job_status import JobStatusEnum
from qnexus.models.references import (
    CircuitRef,
    CompilationResultRef,
    ExecutionResultRef,
    JobType
)

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES, Gate

BASE_SIMULATION_COST: int = 5

# https://docs.quantinuum.com/systems/user_guide/hardware_user_guide/h2.html#native-gate-set

# Map SympleQ gate singletons to (pytket OpType, params) pairs (qubit-only).
# `params` are angle parameters in half-turns for parameterized gates.
_GATE_MAP: dict[Gate, tuple[OpType, list[float]]] = {
    # Id is encoded as Rz(0), a native H-series op the Nexus rebase keeps,
    # so 1-qubit gate counts survive compilation (OpType.noop gets removed).
    GATES.Id: (OpType.Rz, [0.0]),
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


# Clifford X-axis rotations: index k gives Rx(k/2 half-turns), modulo 8 (i.e., 4 half-turns).
# Up-to-global-phase, since Rx(2 half-turns) = -I.
_RX_CLIFFORD: dict[int, list[Gate]] = {
    0: [],                  # Rx(0)
    1: [GATES.V],           # Rx(0.5) = √X
    2: [GATES.X],           # Rx(1)   = X (up to phase)
    3: [GATES.V_inv],       # Rx(1.5)
    4: [],                  # Rx(2)   = -I
    5: [GATES.V],           # Rx(2.5)
    6: [GATES.X],           # Rx(3)
    7: [GATES.V_inv],       # Rx(3.5)
}

# Clifford Z-axis rotations.
_RZ_CLIFFORD: dict[int, list[Gate]] = {
    0: [],
    1: [GATES.S],
    2: [GATES.Z],
    3: [GATES.S_inv],
    4: [],
    5: [GATES.S],
    6: [GATES.Z],
    7: [GATES.S_inv],
}


def _clifford_rotation(half_turns: float, table: dict[int, list[Gate]], axis: str) -> list[Gate]:
    k = round(half_turns * 2)
    if not np.isclose(half_turns, k / 2.0):
        raise ValueError(f"Non-Clifford R{axis} angle: {half_turns} half-turns.")
    return table[k % 8]


@lru_cache(maxsize=None)
def _nexus_compilation_pass(optimisation_level: int = 0) -> BasePass:
    """The compilation pass Nexus runs server-side for H-series devices.

    Built from offline H2-1 device data, whose gate set is identical to the
    one H2-Emulator reports.
    """
    backend = QuantinuumBackend(device_name="H2-1", api_handler=QuantinuumAPI())
    return backend.default_compilation_pass(optimisation_level=optimisation_level)


def to_pytket_circuit(circuit: Circuit) -> PytketCircuit:
    """
    Convert a SympleQ Circuit to a pytket Circuit.

    All qubits are measured and the result is compiled with the same
    default pass Nexus applies server-side (optimisation level 0).

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

    tk_circuit.measure_all()

    _nexus_compilation_pass().apply(tk_circuit)
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
        if op_type == OpType.Measure:
            continue

        decomposition: list[Gate]
        if op_type == OpType.ZZPhase:
            # Only the Clifford angle θ=-0.5 maps back (the inverse of ZZMax).
            theta = float(command.op.params[0]) % 4.0
            if np.isclose(theta, 3.5):
                decomposition = [GATES.ZZMax_inv]
            elif np.isclose(theta, 0.5):
                decomposition = [GATES.ZZMax]
            else:
                raise ValueError(
                    f"Only Clifford ZZPhase(±0.5) is supported, got θ={theta} half-turns."
                )
        elif op_type == OpType.Rz:
            half_turns = float(command.op.params[0])
            k = round(half_turns * 2)
            if np.isclose(half_turns, k / 2.0) and k % 8 == 0:
                # Rz(0) encodes GATES.Id (see _GATE_MAP); keep it explicit so
                # gate counts survive the roundtrip.
                decomposition = [GATES.Id]
            else:
                decomposition = _clifford_rotation(half_turns, _RZ_CLIFFORD, "z")
        elif op_type == OpType.Rx:
            decomposition = _clifford_rotation(float(command.op.params[0]), _RX_CLIFFORD, "x")
        elif op_type == OpType.PhasedX:
            # PhasedX(θ, φ) = Rz(φ)·Rx(θ)·Rz(-φ). For Clifford (θ, φ) — both
            # multiples of 0.5 half-turns — decompose into a sequence of
            # {V, V_inv, X} (for Rx) and {S, S_inv, Z} (for Rz).
            theta = float(command.op.params[0])
            phi = float(command.op.params[1])
            k = round(theta * 2)
            if np.isclose(theta, k / 2.0) and k % 4 == 2:
                # Odd θ is a π rotation about an axis in the XY plane, equal to
                # Rz(2φ)·X up to global phase. This covers the diagonal-axis
                # Cliffords (φ an odd multiple of 0.25) the squash emits.
                rz_gates = _clifford_rotation(2 * phi, _RZ_CLIFFORD, "z")
                if rz_gates == [GATES.Z]:
                    # Z·X ~ Y up to global phase; keep it a single gate so
                    # counts are preserved (Y compiles to PhasedX(1, 0.5)).
                    decomposition = [GATES.Y]
                else:
                    decomposition = [GATES.X] + rz_gates
            else:
                decomposition = (
                    _clifford_rotation(-phi, _RZ_CLIFFORD, "z") +
                    _clifford_rotation(theta, _RX_CLIFFORD, "x") +
                    _clifford_rotation(phi, _RZ_CLIFFORD, "z")
                )
        else:
            mapped = _REVERSE_MAP.get(op_type)
            if mapped is None:
                raise ValueError(f"No SympleQ mapping for pytket gate '{op_type}'.")
            decomposition = [mapped]

        qubits = tuple(qubit.index[0] for qubit in command.qubits)
        gates.extend(decomposition)
        qudit_indices.extend([qubits] * len(decomposition))

    # Rebasing with implicit swaps enabled absorbs SWAPs into a virtual wire
    # permutation (initial -> final) instead of emitting gates; restore it as
    # explicit SWAP gates so the reconstructed circuit is equivalent.
    permutation = {
        qubit.index[0]: target.index[0]
        for qubit, target in tk_circuit.implicit_qubit_permutation().items()
    }
    for pair in _permutation_to_swaps(permutation):
        gates.append(GATES.SWAP)
        qudit_indices.append(pair)

    return Circuit(dimensions, gates, qudit_indices)


def _permutation_to_swaps(permutation: dict[int, int]) -> list[tuple[int, int]]:
    """Decompose a wire permutation (initial -> final) into a SWAP sequence."""
    swaps: list[tuple[int, int]] = []
    visited: set[int] = set()
    for start in permutation:
        if start in visited:
            continue
        cycle = [start]
        nxt = permutation[start]
        while nxt != start:
            cycle.append(nxt)
            nxt = permutation[nxt]
        visited.update(cycle)
        for i in range(len(cycle) - 1, 0, -1):
            swaps.append((cycle[i - 1], cycle[i]))
    return swaps


def fetch_recent_compile_jobs(
    project_name: str,
    n: int,
    device_name: str | None = None,
) -> Iterator[PytketCircuit]:
    """Yield ``(circuit, result)`` pairs from the most recent compilation jobs.

    Parameters
    ----------
    project_name : str
        Name of the Nexus project to query.
    n : int
        Maximum number of compile jobs to fetch.
    device_name : str | None
        If given, only jobs whose ``system.name`` equals ``device_name``
        are kept (filtered client-side, newest-first).

    Yields
    ------
    PytketCircuit
        ``circuit`` across the last ``n`` jobs, ordered
        newest-first.
    """
    project = qnx.projects.get(name=project_name)
    job_iter = qnx.jobs.get_all(
        project=project,
        job_type=[JobType.COMPILE],
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

    print(f"Fetched {len(jobs)} compile jobs.")

    for job in jobs:
        for ref in qnx.jobs.results(job):
            if not isinstance(ref, CompilationResultRef):
                continue
            result = ref.get_output()
            print(ref.get_passes())
            input()
            yield result.download_circuit()


def fetch_recent_execute_jobs(
    project_name: str,
    n: int,
    device_name: str | None = None,
) -> Iterator[tuple[PytketCircuit, BackendResult]]:
    """Yield ``(circuit, result)`` pairs from the most recent execute jobs.

    Parameters
    ----------
    project_name : str
        Name of the Nexus project to query.
    n : int
        Maximum number of execute jobs to fetch.
    device_name : str | None
        If given, only jobs whose ``system.name`` equals ``device_name``
        are kept (filtered client-side, newest-first).

    Yields
    ------
    tuple[PytketCircuit, BackendResult]
        ``(circuit, result)`` pairs across the last ``n`` jobs, ordered
        newest-first.
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

    print(f"Fetched {len(jobs)} execution jobs.")
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
            yield input_program.download_circuit(), result


def get_remaining_quotas() -> dict[str, float | None]:
    """Return remaining Quantinuum quota for each metered resource.

    Returns
    -------
    dict[str, float | None]
        Mapping of quota name to ``quota - usage``, or ``None`` for
        uncapped resources (``quota == "NO_QUOTA_SET"``).
    """
    remaining: dict[str, float | None] = {}
    for q in qnx.quotas.get_all():
        if isinstance(q.quota, str):
            remaining[q.name] = None
        else:
            remaining[q.name] = q.quota - q.usage
    return remaining


def _cost(shot_count: int, one_qubit: int, two_qubits: int, measurements: int) -> float:
    return round(BASE_SIMULATION_COST + shot_count * (one_qubit + 10 * two_qubits + 5 * measurements) / 5000, 2)


def pytket_simulation_cost(circuit: PytketCircuit, shot_count: int = 1) -> float:
    # Assume that num of measurements is circuit.n_qubits. If we have intermediate measurements this is not true.
    return _cost(shot_count, circuit.n_1qd_gates(), circuit.n_2qd_gates(), circuit.n_qubits)


def pytket_bare_simulation_cost(circuit: PytketCircuit, shot_count: int = 1) -> float:
    # Assume that num of measurements is circuit.n_qubits. If we have intermediate measurements this is not true.
    return _cost(shot_count, circuit.n_1qd_gates(), circuit.n_2qd_gates(), circuit.n_qubits) - BASE_SIMULATION_COST
