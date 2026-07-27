import sys
import random
import numpy as np
from pytket.circuit import BitRegister, Circuit, CircBox
from pytket.passes import DecomposeBoxes
from pytket.qasm.qasm import circuit_to_qasm_str
from pytket.backends.backendresult import BackendResult
from pytket.utils.outcomearray import OutcomeArray


from sympleq.integrations.quantinuum.stitching import (reset_operations, circuit_stitching, destitch_results)
from sympleq.applications.randomized_benchmarking.experiments.common import (quantinuum_emulator_backend_factory)
from sympleq.integrations.quantinuum.workflow import run_circuits_on_device


def random_native_deterministic_circuit(
    nq: int,
    seed: int,
    max_x: int | None = None,
) -> tuple[Circuit, tuple[int, ...]]:
    """
    Make a deterministic circuit using X, ZZPhase and measurements.

    The expected output is known exactly:
    a qubit gives 1 iff we applied X to it an odd number of times.

    ZZPhase is included because it is diagonal, so it should not change
    computational-basis measurement results.
    """
    rng = random.Random(seed)

    c = Circuit(nq, nq)

    expected = [0] * nq

    if max_x is None:
        max_x = max(1, nq)

    n_x = rng.randint(0, max_x)

    for _ in range(n_x):
        q = rng.randrange(nq)
        c.X(q)
        expected[q] ^= 1

    # Add some diagonal two-qubit gates. These should not change bit outcomes.
    if nq >= 2:
        n_zz = rng.randint(0, nq - 1)
        for _ in range(n_zz):
            q0, q1 = rng.sample(range(nq), 2)
            c.ZZPhase(0.5, q0, q1)

    for q in range(nq):
        c.Measure(q, q)

    return c, tuple(expected)


def one_shot_tuple(result) -> tuple[int, ...]:
    """
    Convert a one-shot BackendResult into a plain tuple of bits.
    """
    shots = np.asarray(result.get_shots())

    assert shots.shape[0] == 1, f"Expected one shot, got shape {shots.shape}"

    return tuple(int(x) for x in shots[0])


def print_first_commands(circ: Circuit, n: int = 20, title: str = "Circuit"):
    print(f"\n{title}: first {n} commands")
    for i, cmd in enumerate(circ.get_commands()[:n]):
        print(f"{i:03d}: {cmd}")


if __name__ == "__main__":
    ' Include more than 10 circuits to test lexicographic ordering'
    n_circuits = 11

    circuits_sent = []

    for i in range(n_circuits):
        # Different qubit numbers: 1, 2, ..., 10, repeated.
        nq = 1 + (i % 10)
        circ, expected = random_native_deterministic_circuit(
            nq=nq,
            seed=10_000 + i,
        )

        circuits_sent.append((i, circ, expected))

    # print(f"qubits in sent_circuits: {[c[1].n_qubits for c in circuits_sent]}")

    stitched = circuit_stitching([c[1] for c in circuits_sent])

    raw_registers = list(stitched.c_registers)

    raw_names = [r.name for r in raw_registers]
    lex_names = [r.name for r in sorted(raw_registers, key=lambda r: r.name)]
    num_names = [r.name for r in sorted(raw_registers, key=lambda r: int(r.name.replace("creg_", "")))]

    # print("\nRaw stitched register order:")
    # print(raw_names[:20])

    # print("\nLexicographic sorted register order:")
    # print(lex_names[:20])

    # print("\nNumeric sorted register order:")
    # print(num_names[:20])

    # Expected_results:
    # print(f"outcomes: {[c[2] for c in circuits_sent]}")

    circuits_with_expected_sorted = sorted(
        circuits_sent,
        key=lambda item: item[1].n_qubits,
        reverse=True,
    )

    # print_first_commands(stitched, n=30, title="Stitched circuit") # Check the descending order of qubits
    # without sending to the emulator

    answer = input(
        "\nThis will submit a job to the Quantinuum emulator and may use about "
        "10 minutes of emulator time. Continue? [y/N]: "
    ).strip().lower()

    if answer not in {"y", "yes"}:
        print("Cancelled before submitting to Quantinuum emulator.")
        sys.exit(0)

    results = run_circuits_on_device(
        [stitched],
        n_shots=1,
        device_name="H2-Emulator",
        project_name="Stitching Test",
        verbose=True,
    )

    # print('lexicographic_check')

    stitched_result = results[0]

    # print(stitched_result)

    registers = sorted(stitched.c_registers,
                       key=lambda register: int(register.name.removeprefix("creg_")))  # Line 117 quantinuum.py

    destitched = destitch_results(stitched_result, registers)

    assert len(destitched) == n_circuits, (
        f"Expected {n_circuits} destitched results, got {len(destitched)}"
    )

    failures = []

    for destitched_index, result in enumerate(destitched):
        original_index, circ, expected = circuits_with_expected_sorted[destitched_index]

        actual = one_shot_tuple(result)
        print(actual)
        print(one_shot_tuple(result))

        if actual != expected:
            failures.append(
                {
                    "destitched_index": destitched_index,
                    "original_index": original_index,
                    "n_qubits": circ.n_qubits,
                    "expected": expected,
                    "actual": actual,
                    "counts": result.get_counts(),
                }
            )

    if failures:
        print(f"FAILED: {len(failures)} circuits did not match.")

        for failure in failures[:10]:
            print(failure)

        raise AssertionError("Destitched results do not match expected circuit outputs.")

    print(f"PASSED: all {n_circuits} stitched/destitched circuit results matched.")

