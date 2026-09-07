"""Local tests for stitching circuits with different qubit counts.

This file checks the local stitching order: circuits should be stitched
from largest to smallest qubit count, so the resulting classical registers
should have descending sizes.

An ordering issue can appear when stitched registers are later
de-stitched with more than 10 circuits. Lexicographic sorting would place
``creg_10`` before ``creg_2``.  The``destitch_results`` helper trusts the register
list it is passed and does not sort internally. So the desticthcing must obey the 'stitch order'
(Line 114 of quantinuum.py, in the "application/rmb" branch).

The RMB algorithms (Cost-Aware and FLE) guard against this by sorting stitched registers
numerically before calling ``destitch_results``.  Same is implemented here in line 58.
"""
from pytket.circuit import Circuit
from sympleq.integrations.quantinuum.stitching import circuit_stitching


def _dummy_circuit(n_qubits: int, seed: int) -> Circuit:
    """Build a circuit n qubits."""

    circuit = Circuit(n_qubits, n_qubits)

    for qubit in range(n_qubits):
        circuit.Measure(qubit, qubit)

    return circuit


class TestCircuitStitching:
    """
    Local stitching checks that circuits with different qubit counts are
    stitched in descending qubit order.
    """

    def test_stitches_circuits_in_descending_qubit_order(self):
        circuits = [_dummy_circuit(n_qubits=1 + index % 10, seed=index)
                    for index in range(9)
                    ]
        stitched = circuit_stitching([circuit for circuit in circuits])

        registers = sorted(
            stitched.c_registers,
            key=lambda register: int(register.name.removeprefix("creg_")),
        )

        expected_register_sizes = [
            circuit.n_qubits
            for circuit in sorted(circuits, key=lambda item: item.n_qubits, reverse=True)
        ]

        assert [register.size for register in registers] == expected_register_sizes
