from __future__ import annotations
from numpy.random import Generator as RNGGenerator
from pytket.circuit import OpType

from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.integrations.quantinuum.utils import (
    NATIVE_GATES_SET,
    fetch_recent_execute_jobs,
    to_pytket_circuit,
)


class QuantinuumBackend(RMBBackend):
    """
    Quantinuum hardware/emulator backend.

    Builds a random circuit from the config, submits it to the
    Quantinuum device identified by ``device_name`` for ``n_shots``
    shots, and returns whether the most-common measurement outcome is
    the all-zeros bitstring.

    Parameters
    ----------
    device_name : str
        Quantinuum device identifier (e.g. ``"H2-1LE"``).
    n_shots : int
        Number of shots to request from the backend.
    """

    type = "quantinuum"

    def __init__(self, device_name: str = "H2-2E", n_shots: int = 1, project_name: str = "Benchmark") -> None:
        self.device_name = device_name
        self.n_shots = n_shots
        self.project_name = project_name
        self.batch_size = 11

    def fidelity_estimation(self, config: RMBConfig, rng: RNGGenerator) -> list[bool]:
        from sympleq.integrations.quantinuum.workflow import run_circuits_on_device
        circuits = []
        while len(circuits) < self.batch_size:
            circuit = to_pytket_circuit(config.random_circuit(rng=rng))
            if circuit.n_gates == 0:
                continue
            two_qudit_ratio = circuit.n_2qb_gates() / circuit.n_gates
            if not (config.min_two_qubit_gate_ratio
                    <= two_qudit_ratio
                    <= config.max_two_qubit_gate_ratio):
                continue
            circuit.measure_all()
            circuits.append(circuit)

        distributions = run_circuits_on_device(
            circuits, self.n_shots, self.device_name, self.project_name, verbose=True, attempt_batching=False)

        results = []
        for distribution in distributions:
            counts = distribution.as_counter()
            if not counts:
                continue

            outcome, count = counts.most_common()[0]
            assert count <= self.n_shots
            # FIXME: compare to initial state
            results.append(all(bit == 0 for bit in outcome))

        return results

    @classmethod
    def default_config(cls) -> RMBConfig:
        return RMBConfig.default()\
            .with_depth(20)\
            .with_random_elimination(0.1)\
            .with_n_qubits(6)\
            .with_two_qubit_gate_ratio(0.4, 0.6)\
            .with_scrambling_probability(0.5)\
            .with_gates_set(tuple(NATIVE_GATES_SET))

    def default_estimator(self) -> BayesianEstimator:
        """Return a fresh :class:`BayesianEstimator` with ``threshold=0.1`` and ``min_runs=self.batch_size``."""
        return BayesianEstimator(threshold=10**(-1), min_runs=self.batch_size)

    def populate_from_recent_jobs(self, n: int) -> RMBData:
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
        for circuit, result in fetch_recent_execute_jobs(
                self.project_name, n, device_name=self.device_name):
            n_total = circuit.n_gates - circuit.n_gates_of_type(OpType.Measure)
            if n_total == 0:
                continue
            n_2q = circuit.n_2qb_gates()

            counts = result.get_empirical_distribution().as_counter()
            if not counts:
                continue
            outcome, _ = counts.most_common()[0]

            ratio = round(n_2q / n_total, 2)
            config = RMBConfig(
                depth=n_total,
                n_qubits=circuit.n_qubits,
                min_two_qubit_gate_ratio=ratio,
                max_two_qubit_gate_ratio=ratio,
            )
            estimator = data.setdefault(config, self.default_estimator())
            estimator.record(all(bit == 0 for bit in outcome))
        return data

    def to_dict(self) -> dict:
        return {
            "type": self.type,
            "device_name": self.device_name,
            "n_shots": self.n_shots,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> QuantinuumBackend:
        return cls(device_name=payload["device_name"], n_shots=payload["n_shots"])
