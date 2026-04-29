from __future__ import annotations
from numpy.random import Generator as RNGGenerator

from sympleq.applications.randomized_benchmarking.backends.base import RMBBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig


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
        self.batch_size = 10

    def fidelity_estimation(self, config: RMBConfig, rng: RNGGenerator) -> list[bool]:
        from sympleq.integrations.quantinuum.workflow import run_circuits_on_device
        circuits = [config.random_circuit(rng=rng) for _ in range(self.batch_size)]
        distributions = run_circuits_on_device(
            circuits, self.n_shots, self.device_name, self.project_name, verbose=True)

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

    def to_dict(self) -> dict:
        return {
            "type": self.type,
            "device_name": self.device_name,
            "n_shots": self.n_shots,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> QuantinuumBackend:
        return cls(device_name=payload["device_name"], n_shots=payload["n_shots"])
