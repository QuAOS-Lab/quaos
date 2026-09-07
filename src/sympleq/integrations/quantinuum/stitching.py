"""
Modified from https://docs.quantinuum.com/systems/trainings/knowledge_articles/circuit_stitching.html
"""

import sys

from pytket.circuit import BitRegister, Circuit, CircBox
from pytket.passes import DecomposeBoxes
from pytket.qasm.qasm import circuit_to_qasm_str
from pytket.backends.backendresult import BackendResult
from pytket.utils.outcomearray import OutcomeArray


MAX_QASM_PROGRAM_SIZE: int = 6 * 1024 * 1024


def reset_operations(
    n_qubits: int
) -> CircBox:
    r"""Generate a n-qubit CircBox instance containing OpType.Reset operations.

    :param n_qubits: Number of qubits used in the CircBox instance
    :param_type int:
    :returns: CircBox
    """
    circuit = Circuit(n_qubits)
    circuit.name = "Reset"
    for q in circuit.qubits:
        circuit.Reset(q)
    return CircBox(circuit)


def circuit_stitching(
    input_circuits: list[Circuit],
) -> Circuit:
    r"""Generate a stitched circuit based on a list of input circuits.
    The circuit is depth-wise.

    :param input_circuits: Circuit instances to stitch.
    :param backend_info: BackendInfo instance containing information on number
        of device qubits, number of allowed classical register and maximum width
        for each classical register.
    :returns: Circuit
    """

    n_qubits = max([c.n_qubits for c in input_circuits])

    # Put circuits in descending order of number of qubits.
    input_circuits = sorted(input_circuits, key=lambda c: c.n_qubits, reverse=True)

    sum_circuit = Circuit(n_qubits)
    # reset_box = reset_operations(s_circuit.n_qubits)

    creg_index = 0
    for idx in range(len(input_circuits)):
        s_circuit = input_circuits[idx]
        # An input may carry several classical registers (e.g. when it is
        # itself an already-stitched circuit), so wire every one of its bit
        # registers, in the lexicographic order add_circbox_regwise expects.
        cregs = []
        for src_creg in sorted(s_circuit.c_registers, key=lambda r: r.name):
            cregs.append(sum_circuit.add_c_register(f"creg_{creg_index}", src_creg.size))
            creg_index += 1

        qreg = s_circuit.q_registers
        sum_circuit.add_circbox_regwise(CircBox(s_circuit), qreg, cregs)
        if idx == len(input_circuits) - 1:
            continue
        local_reset_box = reset_operations(s_circuit.n_qubits)
        sum_circuit.add_circbox(local_reset_box, s_circuit.qubits)

    # We can also reset the whole circuit after stitching, as the circuits are stitched in
    # descending order of number of qubits. However, the local reset option is more general.
    # To reset the whole circuit:
    # uncomment line 52 (reset_box = reset_operations(s_circuit.n_qubits)),
    # remove line 69 i.e., local_reset_box = reset_operations(s_circuit.n_qubits)), and
    # replace line 70 with: sum_circuit.add_circbox(reset_box, sum_circuit.qubits)

    # Flatten the CircBoxes into native gates so the stitched circuit is a
    # single genuine circuit. This is what lets gate-count-based cost and
    # QASM-size estimates see the stitched contents - gates inside an
    # undecomposed CircBox are invisible to n_1qd_gates()/n_2qd_gates(), so
    # without this the running cost never grows as more circuits are stitched.
    DecomposeBoxes().apply(sum_circuit)

    return sum_circuit


def destitch_results(
    stitched_result: BackendResult,
    registers: list[BitRegister],
) -> list[BackendResult]:
    r"""Split a stitched result into one result per classical register.

    Each stitched sub-circuit writes its measurements to its own classical
    register, so de-stitching reads back the shots restricted to each
    register's bits.

    :param stitched_result: Result of running the stitched circuit.
    :param registers: Classical registers, one per stitched sub-circuit.
    :returns: One :class:`BackendResult` per register, in the given order.
    """
    destitched_results = []
    for register in registers:
        bits = [register[i] for i in range(register.size)]
        outcome_array = OutcomeArray.from_readouts(stitched_result.get_shots(cbits=bits))
        destitched_results.append(BackendResult(shots=outcome_array))
    return destitched_results


def estimate_qasm_program_size(
    circuit: Circuit
) -> int:
    qasm_str = circuit_to_qasm_str(circuit, header="hqslib1")
    return sys.getsizeof(qasm_str) // 1024**2
