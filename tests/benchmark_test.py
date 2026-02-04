import pytest
import numpy as np

from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES
from sympleq.core.paulis.pauli_string import PauliString
from sympleq.core.paulis.pauli_sum import PauliSum


@pytest.mark.benchmark(group="PauliSum")
def test_paulisum_sum(benchmark):
    n_paulis = 20
    n_qudits = 100
    dimensions = np.random.randint(2, 8, size=n_qudits)
    ps1 = PauliSum.from_random(n_paulis, dimensions, rand_weights=True)
    ps2 = PauliSum.from_random(n_paulis, dimensions, rand_weights=True)

    def sum():
        _ = ps1 + ps2

    benchmark(sum)


@pytest.mark.benchmark(group="PauliSum")
def test_paulisum_multiplication(benchmark):
    n_paulis = 20
    n_qudits = 100
    dimensions = np.random.randint(2, 8, size=n_qudits)
    ps1 = PauliSum.from_random(n_paulis, dimensions, rand_weights=True)
    ps2 = PauliSum.from_random(n_paulis, dimensions, rand_weights=True)

    def multiply():
        _ = ps1 * ps2

    benchmark(multiply)


@pytest.mark.benchmark(group="PauliSum")
def test_paulistring_amend(benchmark):
    n_qudits = 100
    dimensions = np.random.randint(2, 8, size=n_qudits)
    ps = PauliString.from_random(dimensions)

    def amend_all():
        for i in range(n_qudits):
            new_x = np.random.randint(0, dimensions[i])
            new_z = np.random.randint(0, dimensions[i])
            ps.amend(i, new_x, new_z)

    benchmark(amend_all)


@pytest.mark.benchmark(group="PauliSum")
def test_paulisum_amend(benchmark):
    n_paulis = 20
    n_qudits = 100
    dimensions = np.random.randint(2, 8, size=n_qudits)
    ps = PauliSum.from_random(n_paulis, dimensions, rand_weights=True)

    def amend_all():
        for i in range(n_paulis):
            idx = np.random.randint(0, n_qudits)
            new_x = np.random.randint(0, dimensions[idx])
            new_z = np.random.randint(0, dimensions[idx])
            ps.select_pauli_string(i).amend(idx, new_x, new_z)

    benchmark(amend_all)


@pytest.mark.benchmark(group="PauliSum")
def test_paulisum_delete_qudits(benchmark):
    n_paulis = 20
    n_qudits = 100
    dimensions = np.random.randint(2, 8, size=n_qudits)
    ps = PauliSum.from_random(n_paulis, dimensions, rand_weights=True)

    def delete_random():
        qudit_indices = np.random.randint(0, ps.n_qudits(), size=ps.n_qudits() // 10).tolist()
        ps._delete_qudits(qudit_indices)

    benchmark(delete_random)


@pytest.mark.benchmark(group="Gate")
def test_hadamard_paulisum_benchmark(benchmark):
    n_paulis = 20
    n_qudits = 100
    dimensions = np.random.randint(2, 8, size=n_qudits)

    gate = GATES.H

    ps = PauliSum.from_random(n_paulis, dimensions, rand_weights=True)

    def apply_gate():
        _ = gate.act(ps, 0)

    benchmark(apply_gate)


@pytest.mark.benchmark(group="Circuit")
@pytest.mark.parametrize("n_qudits", [10, 100])
@pytest.mark.parametrize("n_gates", [50, 200])
def test_circuit_composite_gate(benchmark, n_qudits: int, n_gates: int):
    dimensions = [2] * n_qudits
    circuit = Circuit.from_random(n_gates, dimensions)

    def comp_gate():
        _ = circuit.composite_gate()

    benchmark(comp_gate)
