import sys
import random
from pytket.circuit import BitRegister, Circuit, CircBox
from pytket.passes import DecomposeBoxes
from pytket.qasm.qasm import circuit_to_qasm_str
from pytket.backends.backendresult import BackendResult
from pytket.utils.outcomearray import OutcomeArray


from sympleq.integrations.quantinuum.stitching import (reset_operations, circuit_stitching, destitch_results)
from sympleq.applications.randomized_benchmarking.experiments.common import (quantinuum_emulator_backend_factory)
from sympleq.integrations.quantinuum.workflow import run_circuits_on_device


def random_x_circuit(nq: int, n_x: int, seed: int | None = None) -> Circuit:
    rng = random.Random(seed)

    circ = Circuit(nq)
    nx_list = rng.sample(range(nq), n_x)
    for q in nx_list:
        circ.X(q)
    circ.measure_all()

    return circ

def random_x_circuit1(nq: int, n_x: int, seed: int | None = None) -> Circuit:
    rng = random.Random(seed)

    circ = Circuit(nq)
    nx_list = rng.sample(range(nq), n_x)
    for q in nx_list:
        circ.X(q)
    circ.ZZPhase(0.5, 0, 1)
    circ.measure_all()

    return circ


if __name__ == "__main__":
    circ1 = random_x_circuit1(nq=2, n_x=1, seed=123)
    circ2 = random_x_circuit(nq=6, n_x=2, seed=456)

    print('Circuit-1')
    for command in circ1.get_commands():
        print(command)
    print()
    print('Circuit-2')
    for command in circ2.get_commands():
        print(command)
    print()

    stitched = circuit_stitching([circ1, circ2])
    print('Stitched')
    for cmd in stitched:
        print(cmd)

    results = run_circuits_on_device(
        [stitched],
        n_shots=1,
        device_name="H2-Emulator",
        project_name="level-benchmark",
        verbose=True,
    )
    print(results)

    stitched_result = results[0]
    registers = sorted(stitched.c_registers, key=lambda r: r.name)
    destitched = destitch_results(stitched_result, registers)

    for i, result in enumerate(destitched):
        print(f"circuit {i}")
        print(result.get_counts())
