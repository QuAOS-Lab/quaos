import sys
import qnexus as qnx
from qnexus.client.auth import is_logged_in
from qnexus.models.references import IncompleteJobItemRef, CircuitRef
from pytket.backends.backendresult import BackendResult
from pytket.utils.distribution import EmpiricalDistribution
from numpy.random import default_rng

from sympleq.core.circuits import Circuit
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.integrations.quantinuum.utils import to_pytket_circuit


def _info(msg: str):
    print(f"\033[0;37m{msg}\033[0m")


def _ok(msg: str):
    print(f"\033[1;32m{msg}\033[0m")


def _err(msg: str):
    print(f"\033[1;31m  [ERROR]\033[0m {msg}", file=sys.stderr)


def setup(project_name: str):
    if not is_logged_in():
        qnx.login()
        _ok("Logged in")
    active = qnx.context.get_active_project()
    if active is not None and active.annotations.name == project_name:
        return
    project = qnx.projects.get_or_create(name=project_name)
    qnx.context.set_active_project(project)
    _ok(f"Project: {project.annotations.name}")


def build_and_compile_circuit(circuit: Circuit, backend_config: qnx.QuantinuumConfig) -> CircuitRef:
    pytket_circuit = to_pytket_circuit(circuit)
    pytket_circuit.measure_all()
    ref = qnx.circuits.upload(circuit=pytket_circuit, name="Test-Circuit")
    _ok("Circuit built and uploaded")

    ref_compile_job = qnx.start_compile_job(
        programs=[ref],
        backend_config=backend_config,
        optimisation_level=0,
        name="compilation-job",
    )
    qnx.jobs.wait_for(ref_compile_job)

    ref_result = qnx.jobs.results(ref_compile_job)[0]
    if isinstance(ref_result, IncompleteJobItemRef):
        _err("Compilation failed")
        sys.exit(1)
    ref_compiled_circuit = ref_result.get_output()
    compiled_circuit = ref_compiled_circuit.download_circuit()
    _info(f"Sympleq  circuit: qubits={circuit.n_qudits()} gates={circuit.n_gates()} ")
    _info(f"Pytket   circuit: qubits={pytket_circuit.n_qubits} gates={pytket_circuit.n_gates} ")
    _info(f"Compiled circuit: qubits={compiled_circuit.n_qubits} gates={circuit.n_gates()} ")
    _ok("Compilation complete")

    return ref_compiled_circuit


def run_compiled_circuit(ref_compiled_circuit: CircuitRef,
                         n_shots: int,
                         backend_config: qnx.QuantinuumConfig,
                         syntax_checker: str | None = None) -> BackendResult:

    # Syntax checker not available on emulator
    if not backend_config.device_name.endswith("E") and syntax_checker is not None:
        execution_cost = qnx.circuits.cost(
            circuit_ref=ref_compiled_circuit,
            n_shots=n_shots,
            backend_config=backend_config,
            syntax_checker=syntax_checker,
        )

        _info(f"Execution cost: {execution_cost}")
        input("Do you want to proceed?")

    ref_execute_job = qnx.start_execute_job(
        programs=[ref_compiled_circuit],
        n_shots=n_shots,
        backend_config=backend_config,
        name="execution-job",
    )
    _info(f"Status: {qnx.jobs.status(ref_execute_job)}")
    qnx.jobs.wait_for(ref_execute_job)

    ref_result = qnx.jobs.results(ref_execute_job)[0]
    if isinstance(ref_result, IncompleteJobItemRef):
        _err("Execution failed")
        raise Exception(f"Qnexus execution failed: {ref_result}")
    _ok("Execution complete")

    backend_result = ref_result.download_result()
    assert isinstance(backend_result, BackendResult)
    return backend_result


def run_circuit_on_device(circuit: Circuit,
                          n_shots: int,
                          device_name: str,
                          project_name: str,
                          verbose: bool = False) -> EmpiricalDistribution:
    setup(project_name)

    backend_config = qnx.QuantinuumConfig(device_name=device_name, no_opt=True, allow_implicit_swaps=False)
    ref_circuit = build_and_compile_circuit(circuit, backend_config)
    backend_result = run_compiled_circuit(ref_circuit, n_shots, backend_config)

    distribution = backend_result.get_empirical_distribution()

    if verbose:
        counts = distribution.as_counter()
        total = distribution.total
        for state, count in counts.most_common():
            _info(f"{state}: {count / total:.4f} ({count}/{total})")

    return distribution


if __name__ == "__main__":
    device_name = "H2-2E"
    n_shots = 100
    n_gates = 24
    n_qubits = 6
    dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qubits
    rng = default_rng(10)
    circuit = Circuit.from_random(n_gates, dimensions, two_qudit_gate_ratio=0.25, rng=rng)
    circuit = circuit + circuit.inverse()

    run_circuit_on_device(circuit, n_shots, device_name, "Benchmark", verbose=True)
