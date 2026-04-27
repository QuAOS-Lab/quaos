from __future__ import annotations
from dataclasses import dataclass, field, replace
from typing import Callable
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.core.bayesian_estimation import BayesianEstimator
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import GATES, Gate
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.core.paulis.pauli_sum import PauliSum

type RMBData = dict[RMBConfig, BayesianEstimator]
type UpdateStrategy = Callable[[RMBData, RMBConfig], RMBConfig]


@dataclass(frozen=True)  # Frozen to avoid possible mistakes in using references
class RMBConfig:
    """
    Configuration for a randomized benchmarking experiment.

    Frozen dataclass: instances are immutable. Use the ``with_*`` methods
    to obtain a new configuration with one field changed.

    Parameters
    ----------
    depth : int
        Number of gate layers in the random circuit (must be >= 1).
    scrambling_probability : float
        Probability of inserting an X gate (vs Id) on each qudit in the
        scrambler layer wrapping the random circuit. In ``[0, 1]``.
    two_qudit_gate_ratio : float
        Fraction of gates in the random circuit drawn from the two-qudit
        subset of ``gates_set``. In ``[0, 1]``.
    gates_set : tuple[Gate, ...]
        Gates available to sample from when building the random circuit.
    n_qubits : int
        Number of qudits in the system (must be >= 1).
    random_elimination : float
        Probability used to randomly replace single-qudit gates with
        identities in order to make the circuit asymmetric. In ``[0, 1]``.
    """
    depth: int = 1
    scrambling_probability: float = 0.0
    two_qudit_gate_ratio: float = 0.0
    gates_set: tuple[Gate, ...] = (GATES.H, GATES.S, GATES.ZZPhase)
    n_qubits: int = 1
    random_elimination: float = 0.0
    dimensions: np.ndarray = field(init=False, compare=False, hash=False, repr=False)
    _initial_state: PauliSum = field(init=False, compare=False, hash=False, repr=False)

    def __post_init__(self) -> None:
        """Validate fields and initialize the derived ``dimensions`` array and ``_initial_state``."""
        if self.depth < 1:
            raise ValueError(
                f"Invalid depth, it should be larger than 0 (got {self.depth}).")
        if not 0.0 <= self.scrambling_probability <= 1.0:
            raise ValueError(
                f"Invalid scrambling_probability, it should be between 0 and 1 (got {self.scrambling_probability}).")
        if not 0.0 <= self.two_qudit_gate_ratio <= 1.0:
            raise ValueError(
                f"Invalid two_qudit_gate_ratio, it should be between 0 and 1 (got {self.two_qudit_gate_ratio}).")
        if not 0.0 <= self.random_elimination <= 1.0:
            raise ValueError(
                f"Invalid random_elimination, it should be between 0 and 1 (got {self.random_elimination}).")
        if self.n_qubits < 1:
            raise ValueError(
                f"Invalid n_qubits, it should be larger than 0 (got {self.n_qubits}).")
        object.__setattr__(self, "dimensions",
                           np.asarray([DEFAULT_QUDIT_DIMENSION] * self.n_qubits, dtype=int))

        pauli_strings = []
        for p_idx in range(self.n_qubits):
            pauli_string = ""
            for q_idx in range(self.n_qubits):
                if p_idx == q_idx:
                    pauli_string += "x0z1"
                else:
                    pauli_string += "x0z0"
            pauli_strings.append(pauli_string)
        object.__setattr__(self, "_initial_state",
                           PauliSum.from_string(pauli_strings, self.dimensions))

    @classmethod
    def default(cls) -> RMBConfig:
        """Return a sensible default configuration for an RMB sweep."""
        return cls(depth=10, random_elimination=0.1, n_qubits=2, two_qudit_gate_ratio=0.3)

    def with_depth(self, depth: int) -> RMBConfig:
        """Return a copy of this config with ``depth`` replaced."""
        return replace(self, depth=depth)

    def with_scrambling_probability(self, scrambling_probability: float) -> RMBConfig:
        """Return a copy of this config with ``scrambling_probability`` replaced."""
        return replace(self, scrambling_probability=scrambling_probability)

    def with_two_qudit_gate_ratio(self, two_qudit_gate_ratio: float) -> RMBConfig:
        """Return a copy of this config with ``two_qudit_gate_ratio`` replaced."""
        return replace(self, two_qudit_gate_ratio=two_qudit_gate_ratio)

    def with_n_qubits(self, n_qubits: int) -> RMBConfig:
        """Return a copy of this config with ``n_qubits`` replaced."""
        return replace(self, n_qubits=n_qubits)

    def with_random_elimination(self, random_elimination: float) -> RMBConfig:
        """Return a copy of this config with ``random_elimination`` replaced."""
        return replace(self, random_elimination=random_elimination)

    def initial_state(self) -> PauliSum:
        """
        Return the initial state for the benchmark.

        The state is a Pauli sum whose stabilizers are ``Z`` on each
        individual qudit (and identity elsewhere), i.e. the all-zeros
        computational basis state.

        Returns
        -------
        PauliSum
            Initial state encoded as a sum of Pauli strings.
        """
        return self._initial_state

    def random_circuit(self, rng: RNGGenerator | None = None) -> Circuit:
        """
        Generate a random benchmarking circuit.

        Builds a random circuit of depth ``self.depth`` from
        ``self.gates_set`` with the configured two-qudit gate ratio,
        wraps it with a scrambler layer (X or Id per qudit, sampled
        with ``self.scrambling_probability``) and its inverse plus
        the inverse of the random circuit, applies the configured
        noise models, and finally optionally turns matching single-qudit
        gates into identities according to ``self.random_elimination``.

        Parameters
        ----------
        rng : numpy.random.Generator | None
            Random number generator. If ``None``, a fresh ``default_rng()``
            is used.

        Returns
        -------
        Circuit
            The constructed randomized benchmarking circuit.
        """
        if rng is None:
            rng = default_rng()

        _circuit = Circuit.from_depth(self.depth,
                                      self.dimensions,
                                      gates_set=self.gates_set,
                                      two_qudit_gate_ratio=self.two_qudit_gate_ratio,
                                      rng=rng)
        _scrambler = Circuit.empty(self.dimensions)
        for q_idx in range(self.n_qubits):
            if rng.random() <= self.scrambling_probability:
                _scrambler.add_gate(GATES.X, q_idx)
            else:
                _scrambler.add_gate(GATES.Id, q_idx)

        circuit = _scrambler + _circuit + _circuit.inverse() + _scrambler.inverse()

        _initial_state = self.initial_state()

        # Eliminate and insert identity gates from and to the circuit to make it asymmetric.
        # This step is performed without applying errors.
        if self.random_elimination > 0.0:
            pauli = _initial_state
            for idx, (gate, q_idxs) in enumerate(zip(circuit.gates, circuit.qudit_indices)):
                if gate.n_qudits > 1:
                    continue
                intermediate = gate.act(pauli, q_idxs)
                if pauli == intermediate:
                    circuit.gates[idx] = GATES.Id

                pauli = intermediate

        return circuit

    @classmethod
    def default_update_strategy(cls, data: RMBData, current_config: RMBConfig) -> RMBConfig:
        """
        Default rule for advancing an RMB sweep.

        Parameters
        ----------
        data : RMBData
            Mapping from previously seen configurations to their
            Bayesian estimators.
        current_config : RMBConfig
            Configuration that was just run (or about to be run).

        Returns
        -------
        RMBConfig
            New configuration to run next.
        """
        if current_config in data:
            estimator = data[current_config]
        else:
            return RMBConfig.default()
        # probability(True) is the fidelity
        if estimator.probability(True) >= 0.8:
            new_config = current_config.with_depth(current_config.depth + 4)
        elif estimator.probability(True) >= 0.6:
            new_config = current_config.with_depth(current_config.depth + 2)
        elif estimator.probability(True) >= 0.5:
            new_config = current_config.with_depth(current_config.depth + 1)
        else:
            new_value = round(current_config.two_qudit_gate_ratio * 0.9, 2)
            new_config = current_config.with_two_qudit_gate_ratio(new_value)

        return new_config
