import random
import numpy as np
from pytket.circuit import Circuit


from sympleq.applications.randomized_benchmarking.experiments.common import (quantinuum_emulator_backend_factory)


from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.integrations.quantinuum.utils import NATIVE_GATES_SET
from numpy.random import default_rng
from sympleq.applications.randomized_benchmarking.backends.base import MeasurementRequest
from sympleq.applications.randomized_benchmarking.experiments.common import (
    quantinuum_emulator_backend_factory
)


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

    if nq >= 2:
        n_zz = rng.randint(0, nq - 1)
        for _ in range(n_zz):
            q0, q1 = rng.sample(range(nq), 2)
            c.ZZPhase(0.5, q0, q1)

    for q in range(nq):
        c.Measure(q, q)

    return c, tuple(expected)


class DummySettings:
    max_cost_per_run: float = 0.0


def main():
    rng = default_rng(1234)
    settings = DummySettings()

    backend = quantinuum_emulator_backend_factory(settings, rng)

    configs = [
        (
            RMBConfig.default()
            .with_n_qubits(nq)
            .with_n_1qb_gates(2 * nq)
            .with_n_2qb_gates(2)
            .with_random_elimination(0.0)
            .with_use_scrambler(True)
            .with_gates_set(tuple(NATIVE_GATES_SET))
        )
        for nq in [2, 5, 3, 8]
    ]
    print("\n[request order]")
    for i, config in enumerate(configs):
        print(
            f"request_index={i:02d} "
            f"n_qubits={config.n_qubits} "
            f"n_gates={config.n_gates}"
        )

    expected_stitch_order = sorted(
        configs,
        key=lambda config: config.n_qubits,
        reverse=True,
    )

    print("\n[expected stitch order]")
    for i, config in enumerate(expected_stitch_order):
        print(
            f"stitched_index={i:02d} "
            f"n_qubits={config.n_qubits} "
            f"n_gates={config.n_gates}"
        )

    requests = [
        MeasurementRequest(config=config, shots=1)
        for config in configs
    ]

    outcomes = backend.fidelity_estimation(requests, rng)

    print("\n[outcome dict order]")
    for i, (config, values) in enumerate(outcomes.outcomes.items()):
        print(
            f"outcome_index={i:02d} "
            f"n_qubits={config.n_qubits} "
            f"n_gates={config.n_gates} "
            f"outcomes={values}"
        )


if __name__== "__main__":
    main()
