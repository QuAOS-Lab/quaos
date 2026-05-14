from pytket import Circuit
import qnexus as qnx
from qnexus.client.auth import is_logged_in
from qnexus.models.references import IncompleteJobItemRef, CircuitRef
from pytket.backends.backendresult import BackendResult
from pytket.utils.distribution import EmpiricalDistribution
from numpy.random import default_rng

from sympleq.core.circuits import Circuit as SympleqCircuit
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.integrations.quantinuum.utils import to_pytket_circuit


def _info(msg: str):
    print(f"\033[0;37m{msg}\033[0m")


def _ok(msg: str):
    print(f"\033[1;32m{msg}\033[0m")


def default_backend_config(device_name: str) -> qnx.QuantinuumConfig:
    return qnx.QuantinuumConfig(
        device_name=device_name,
        no_opt=True,
        allow_implicit_swaps=False,
        leakage_detection=False,
        attempt_batching=False)


def stabilizer_backend_config(device_name: str) -> qnx.QuantinuumConfig:
    """Backend config selecting Quantinuum's stabilizer emulator (Clifford-only)."""
    return qnx.QuantinuumConfig(
        device_name=device_name,
        simulator="stabilizer",
        no_opt=True,
        allow_implicit_swaps=False,
        leakage_detection=False,
        attempt_batching=False)


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


def build_and_compile_circuits(circuits: list[Circuit],
                               backend_config: qnx.QuantinuumConfig | None = None,
                               name: str | None = None) -> list[CircuitRef]:

    programs = []
    name = "" if name is None else f"{name}-"
    for idx, circuit in enumerate(circuits):
        ref = qnx.circuits.upload(circuit=circuit, name=f"{name}Circuit_{idx}")
        programs.append(ref)

    _ok("Circuits built and uploaded")

    if backend_config is None:
        backend_config = default_backend_config(device_name)

    ref_compile_job = qnx.start_compile_job(
        programs=programs,
        backend_config=backend_config,
        optimisation_level=0,
        name="compilation-job",
    )
    qnx.jobs.wait_for(ref_compile_job)

    ref_results = qnx.jobs.results(ref_compile_job)
    if isinstance(ref_results, IncompleteJobItemRef):
        raise Exception(f"Compilation failed: {ref_results}")

    outputs = []
    for ref_result in ref_results:
        if isinstance(ref_result, IncompleteJobItemRef):
            continue

        ref_compiled_circuit = ref_result.get_output()
        compiled_circuit = ref_compiled_circuit.download_circuit()
        # NOTE: The discrepancy between Pytket and Sympleq is due to the measure_all() operation from above,
        # which adds extra measurement gates.
        _info(f"Pytket   circuit: qubits={circuit.n_qubits} gates={circuit.n_gates}")
        _info(f"Compiled circuit: qubits={compiled_circuit.n_qubits} gates={compiled_circuit.n_gates}")

        outputs.append(ref_compiled_circuit)
    _ok("Compilation complete")

    return outputs


def run_compiled_circuits(ref_circuits: list[CircuitRef],
                          n_shots: int,
                          backend_config: qnx.QuantinuumConfig,
                          syntax_checker: str | None = None) -> list[BackendResult]:

    # Syntax checker not available on emulator
    if syntax_checker is not None:
        try:
            execution_cost = qnx.circuits.cost(
                circuit_ref=ref_circuits,
                n_shots=n_shots,
                backend_config=backend_config,
                syntax_checker=syntax_checker,
            )

            _info(f"Execution cost: {execution_cost}")
        except Exception as err:
            _info(f"Syntax checker not available, cannot estimate cost: {err}.")

        input("Do you want to proceed?")

    ref_execute_job = qnx.start_execute_job(
        programs=ref_circuits,  # type: ignore
        n_shots=n_shots,
        backend_config=backend_config,
        name="execution-job",
    )
    qnx.jobs.wait_for(ref_execute_job)

    ref_results = qnx.jobs.results(ref_execute_job)
    if isinstance(ref_results, IncompleteJobItemRef):
        raise Exception(f"Qnexus execution failed: {ref_results}")
    _ok("Execution complete")

    backend_results = []
    for ref_result in ref_results:
        if isinstance(ref_result, IncompleteJobItemRef):
            continue
        backend_results.append(ref_result.download_result())
    return backend_results


def run_circuits_on_device(circuits: list[Circuit],
                           n_shots: int,
                           device_name: str,
                           project_name: str,
                           verbose: bool = False) -> list[EmpiricalDistribution]:
    setup(project_name)

    backend_config = default_backend_config(device_name)
    ref_circuits = build_and_compile_circuits(circuits, backend_config, project_name)
    backend_results = run_compiled_circuits(ref_circuits, n_shots, backend_config)

    distributions = [backend_result.get_empirical_distribution() for backend_result in backend_results]

    if verbose:
        for distribution in distributions:
            counts = distribution.as_counter()
            total = distribution.total
            for state, count in counts.most_common():
                _info(f"{state}: {count / total:.4f} ({count}/{total})")

    return distributions


if __name__ == "__main__":
    device_name = "H2-2E"
    n_shots = 100
    n_gates = 24
    n_qubits = 6
    dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qubits
    rng = default_rng(10)
    circuit = SympleqCircuit.from_random(n_gates, dimensions, two_qudit_gate_ratio=0.25, rng=rng)
    circuit = circuit + circuit.inverse()
    p_circuit = to_pytket_circuit(circuit)

    run_circuits_on_device([p_circuit], n_shots, device_name, "Benchmark", verbose=True)
