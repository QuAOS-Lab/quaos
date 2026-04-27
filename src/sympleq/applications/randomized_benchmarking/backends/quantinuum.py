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

    def __init__(self, device_name: str = "H2-1LE", n_shots: int = 1, project_name: str = "Benchmark") -> None:
        self.device_name = device_name
        self.n_shots = n_shots
        self.project_name = project_name

    def fidelity_estimation(self, config: RMBConfig, rng: RNGGenerator) -> bool:
        from sympleq.integrations.quantinuum.workflow import run_circuit_on_device
        circuit = config.random_circuit(rng=rng)
        distribution = run_circuit_on_device(circuit, self.n_shots, self.device_name, self.project_name, verbose=True)
        counts = distribution.as_counter()
        if not counts:
            return False

        # FIXME: update to get self.n_shots most common
        outcomes = counts.most_common()
        outcome, count = outcomes[0]
        assert count <= self.n_shots
        return all(bit == 0 for bit in outcome)

    def to_dict(self) -> dict:
        return {
            "type": self.type,
            "device_name": self.device_name,
            "n_shots": self.n_shots,
        }

    @classmethod
    def from_dict(cls, payload: dict) -> QuantinuumBackend:
        return cls(device_name=payload["device_name"], n_shots=payload["n_shots"])
