from __future__ import annotations

from dataclasses import dataclass

from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.applications.randomized_benchmarking.config import RMBConfig
from sympleq.core.circuits import Circuit
from sympleq.core.circuits.gates import GATES


H_WRAPPER_VARIANT = "native_v_wrapper"
H_WRAPPER_1Q_GATES_PER_QUBIT = 2


def native_v_layer(dimensions) -> Circuit:
    """Native-V wrapper layer, applied to every qubit."""

    layer = Circuit.empty(dimensions)
    for q_idx in range(len(dimensions)):
        layer.add_gate(GATES.V, q_idx)
    return layer


def h_wrapper_1q_gates(n_qubits: int) -> int:
    return H_WRAPPER_1Q_GATES_PER_QUBIT * int(n_qubits)


def _config_key(config: RMBConfig) -> tuple:
    return (
        int(config.n_1qb_gates),
        int(config.n_2qb_gates),
        tuple(gate.name for gate in config.gates_set),
        int(config.n_qubits),
        float(config.random_elimination),
        bool(config.use_scrambler),
    )


@dataclass(frozen=True, eq=False)
class HWrappedRMBConfig(RMBConfig):
    """RMBConfig variant whose generated circuit has protected native-H wrappers."""

    def __eq__(self, other) -> bool:
        if isinstance(other, RMBConfig):
            return _config_key(self) == _config_key(other)
        return NotImplemented

    def __hash__(self) -> int:
        return hash(_config_key(self))

    def circuit_metadata(self) -> dict:
        return {
            "circuit_variant": H_WRAPPER_VARIANT,
            "protected_h_wrapper": True,
            "native_h_decomposition": "V,V_inv",
            "h_wrapper_1q_gates": h_wrapper_1q_gates(self.n_qubits),
            "h_wrapper_1q_gates_per_qubit": H_WRAPPER_1Q_GATES_PER_QUBIT,
        }

    def random_circuit(self, rng: RNGGenerator | None = None) -> Circuit:
        if rng is None:
            rng = default_rng()

        target_n_2qb_gates = self.n_2qb_gates // 2

        single_1q_layer = Circuit.empty(self.dimensions)
        scrambler_gates = self.n_qubits if self.use_scrambler else 0
        if self.use_scrambler:
            scrambling_gates = [GATES.X, GATES.Y, GATES.Z]
            for q_idx in range(self.n_qubits):
                gate = scrambling_gates[rng.integers(0, len(scrambling_gates))]
                single_1q_layer.add_gate(gate, q_idx)

        core_1qb_total = self.n_1qb_gates - h_wrapper_1q_gates(self.n_qubits)
        target_n_1qb_gates = core_1qb_total // 2 - scrambler_gates
        if target_n_1qb_gates < 0:
            raise ValueError(
                "n_1qb_gates is too small for the protected H wrapper and "
                f"single-qubit scrambler layer: n_1qb_gates={self.n_1qb_gates}, "
                f"required_min={h_wrapper_1q_gates(self.n_qubits) + 2 * scrambler_gates}."
            )
        body = Circuit.from_number_of_gates(
            target_n_1qb_gates,
            target_n_2qb_gates,
            self.dimensions,
            gates_set=self.gates_set,
            rng=rng,
        )

        left_wrapper = native_v_layer(self.dimensions)
        right_wrapper = left_wrapper.inverse()
        if self.use_scrambler:
            core = single_1q_layer + body + body.inverse() + single_1q_layer.inverse()
        else:
            core = body + body.inverse()

        circuit = left_wrapper + core + right_wrapper

        left_h_len = len(left_wrapper.gates)
        core_len = len(core.gates)
        right_h_len = len(right_wrapper.gates)
        protected_indices = set(range(left_h_len))
        protected_indices.update(
            range(left_h_len + core_len, left_h_len + core_len + right_h_len)
        )

        if self.random_elimination > 0.0:
            pauli = self.initial_state()
            for idx, (gate, q_idxs) in enumerate(zip(circuit.gates, circuit.qudit_indices)):
                if gate.n_qudits > 1:
                    continue
                intermediate = gate.act(pauli, q_idxs)
                if idx not in protected_indices and pauli == intermediate:
                    circuit.gates[idx] = GATES.Id
                pauli = intermediate

        return circuit
