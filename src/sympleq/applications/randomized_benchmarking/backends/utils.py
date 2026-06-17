import datetime
import os
from pathlib import Path
from typing import Callable, Iterator
from pytket.circuit import Circuit as PytketCircuit, OpType

from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.integrations.quantinuum.utils import fetch_recent_execute_jobs, to_pytket_circuit
from sympleq.integrations.quantinuum.workflow import build_and_compile_circuits, default_backend_config, setup


def config_from_pytket_circuit(circuit: PytketCircuit) -> RMBConfig:
    n_total = circuit.n_gates - circuit.n_gates_of_type(OpType.Measure)
    if n_total == 0:
        raise ValueError("Invalid input circuit")
    return RMBConfig(
        n_1qb_gates=circuit.n_1qb_gates(),
        n_2qb_gates=circuit.n_2qb_gates(),
        n_qubits=circuit.n_qubits,
    )


_RMB_DATA_DIR = Path(__file__).resolve().parent.parent / "rmb_data"


def fetch_and_store_circuits(
        project_name: str,
        n: int,
        folder_name: str,
        device_name: str | None = None) -> Path:
    """Fetch the last ``n`` execute jobs and store the converted SympleQ circuits.

    Parameters
    ----------
    project_name : str
        Name of the Nexus project to query.
    n : int
        Maximum number of execute jobs to fetch.
    folder_name : str
        Subdirectory name created (or reused) under ``rmb_data/``.
    device_name : str | None
        If given, only jobs whose ``system.name`` equals ``device_name``
        are kept.

    Returns
    -------
    Path
        The directory the circuits were written to.
    """
    folder = _RMB_DATA_DIR / folder_name
    folder.mkdir(parents=True, exist_ok=True)

    for i, (tk_circuit, _) in enumerate(
            fetch_recent_execute_jobs(project_name, n, device_name=device_name)):
        (folder / f"circuit_{i}.json").write_text(tk_circuit.to_json())

    return folder


def load_pytket_circuits(folder_name: str) -> list[PytketCircuit]:
    circuits = []
    folder = _RMB_DATA_DIR / folder_name
    for filename in os.listdir(folder):
        with open(folder / filename, "r") as f:
            p_circuit_json = f.read()
        p_circuit = PytketCircuit.from_json(p_circuit_json)
        circuits.append(p_circuit)

    return circuits


def data_from_pytket_circuit_results(
        jobs: Iterator[tuple[PytketCircuit, bool]] | list[tuple[PytketCircuit, bool]],
        default_estimator: Callable[[], BayesianEstimator]) -> RMBData:
    """Build :type:`RMBData` from the last ``n`` execute jobs in ``self.project_name``.

    Each fetched circuit's ``n_qubits``, non-measurement gate count, and
    two-qudit gate ratio are used to synthesise an :class:`RMBConfig`;
    remaining fields fall back to defaults. The outcome recorded is
    ``True`` iff the most-common measured bitstring is all-zeros.

    Parameters
    ----------
    n : int
        Maximum number of execute jobs to fetch.

    Returns
    -------
    RMBData
        Mapping of inferred configs to their estimators.
    """
    data: RMBData = {}
    for circuit, result in jobs:
        config = config_from_pytket_circuit(circuit)
        if config is None:
            continue

        estimator = data.setdefault(config, default_estimator())
        estimator.record(result)
    return data


def generate_random_pytket_circuits(
        config: RMBConfig, n_circuits: int,
        folder_name: str,
        project_name: str = "Benchmark",
        device_name: str = "H2-1LE") -> list[PytketCircuit]:
    folder = _RMB_DATA_DIR / folder_name
    folder.mkdir(parents=True, exist_ok=True)
    circuits = [to_pytket_circuit(config.random_circuit()) for _ in range(n_circuits)]
    setup(project_name)
    backend_config = default_backend_config(device_name)
    p_circuits_ref = build_and_compile_circuits(circuits, backend_config)

    timestamp = f"{datetime.datetime.now():%Y-%m-%dT%H-%M-%S}"
    p_circuits = []
    for i, circuit_ref in enumerate(p_circuits_ref):
        p_circuit = circuit_ref.download_circuit()
        p_circuits.append(p_circuit)
        (folder / f"{timestamp}-circuit_{i}.json").write_text(p_circuit.to_json())

    return p_circuits
