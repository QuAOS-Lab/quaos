from __future__ import annotations
from typing import Generator, overload
import json
import numpy as np
import scipy.sparse as sp
from numpy.random import Generator as RNGGenerator, default_rng
from pathlib import Path
from collections import defaultdict
import warnings

from sympleq.core.noise.noise_model import NoiseModel
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from sympleq.core.paulis._typing import TableauType, DimensionsLike, DimensionsType, PhasesType, HilbertOperator

from .utils import embed_unitary
from .gates import Gate, GATES, _GenericGate
from .utils import embed_symplectic
from sympleq.core.paulis import PauliSum, PauliString, PauliObject


# Type alias for from_tuples input: (Gate, qudit_idx1, qudit_idx2, ...)
GateSpec = tuple[Gate, *tuple[int, ...]]


class Circuit:
    """
    A quantum circuit consisting of gates applied to specific qudits.

    The circuit stores:
    - dimensions: the dimension of each qudit (e.g., [2, 2, 2] for 3 qubits)
    - gates: list of Gate objects (dimension-independent)
    - qudit_indices: list of tuples indicating which qudits each gate acts on

    Gates and qudit indices are stored separately because gates are now dimension-independent
    singletons that don't store their target qudits.
    """

    def __init__(self, dimensions: DimensionsType,
                 gates: list[Gate],
                 qudit_indices: list[tuple[int, ...]],
                 use_unitary_cache: bool = True):
        """
        Initialize the Circuit.

        Parameters
        ----------
        dimensions : DimensionsType
            The dimension of each qudit in the circuit.
        gates : list[Gate]
            List of Gate objects.
        qudit_indices : list[tuple[int, ...]]
            List of tuples indicating which qudits each gate acts on.
            Must have the same length as gates.
        """
        self.dimensions = np.asarray(dimensions, dtype=int)
        self.dimensions.setflags(write=False)

        self._gates = gates
        self._qudit_indices = list(qudit_indices)

        self._noise_model_per_gate: list[NoiseModel | None] = [None] * self.n_gates()
        self._use_unitary_cache = use_unitary_cache
        self._unitary_cache: dict[tuple[Gate, tuple[int, ...]], tuple[HilbertOperator, HilbertOperator]] = {}

    @property
    def gates(self) -> list[Gate]:
        """List of gates in the circuit."""
        return self._gates

    @property
    def qudit_indices(self) -> list[tuple[int, ...]]:
        """List of qudit index tuples for each gate."""
        return self._qudit_indices

    @property
    def noise_model_per_gate(self) -> list[NoiseModel | None]:
        """List of qudit index tuples for each gate."""
        return self._noise_model_per_gate

    @classmethod
    def empty(cls, dimensions: DimensionsLike) -> Circuit:
        """
        Create an empty circuit with no gates.

        Parameters
        ----------
        dimensions : DimensionsLike
            The dimension of each qudit.

        Returns
        -------
        Circuit
            An empty circuit.
        """
        dimensions = np.asarray(dimensions, dtype=int)
        C = cls(dimensions, [], [])
        C._sanity_check()

        return C

    @classmethod
    def from_random(cls, n_gates: int,
                    dimensions: DimensionsLike,
                    gates_set: tuple[Gate] | list[Gate] | set[Gate] | None = None,
                    two_qudit_gate_ratio: float = 0.3,
                    rng: RNGGenerator | None = None) -> Circuit:
        """
        Creates a random circuit with the given number of gates.

        Parameters
        ----------
        n_gates : int
            Number of gates in the circuit.
        dimensions : DimensionsLike
            The dimension of each qudit.
        gates_set : tuple[Gate] | list[Gate] | set[Gate]
            The set of gates from which to draw.
        two_qudit_gate_ratio : float
            Probability of choosing a two-qudit gate vs single-qudit gate.
        rng : numpy.random.Generator or None, optional
            Random number generator. If ``None``, a default generator is used.

        Returns
        -------
        Circuit
            A new random Circuit.
        """
        def index_lists(lst):
            groups = defaultdict(list)
            for i, val in enumerate(lst):
                groups[val].append(i)
            return list(groups.values())

        if rng is None:
            rng = default_rng()

        index_sets = index_lists(dimensions)  # list of lists of indexes for each dimension
        n_dims = len(index_sets)  # number of different dimensions

        dimensions = np.asarray(dimensions, dtype=int)
        index_sets = index_lists(dimensions)  # list of lists of indexes for each dimension
        n_dims = len(index_sets)

        if gates_set is None:
            single_qudit_gates: list[Gate] = [GATES.H, GATES.S]
            two_qudit_gates: list[Gate] = [GATES.CX, GATES.SWAP]
        elif isinstance(gates_set, (tuple, list, set)):
            single_qudit_gates = [gate for gate in gates_set if gate.n_qudits == 1]
            two_qudit_gates = [gate for gate in gates_set if gate.n_qudits == 2]
            if any([gate.n_qudits > 2 for gate in gates_set]):
                warnings.warn("Only single qudit and 2-qudits gates are used to generate the circuit.")
        else:
            raise ValueError("Invalid gates_set type.")

        gates = []
        qudit_indices = []

        for _ in range(n_gates):
            set_idx = rng.integers(0, n_dims)
            if rng.random() < two_qudit_gate_ratio and len(index_sets[set_idx]) > 1:
                indices = tuple(int(idx) for idx in rng.choice(index_sets[set_idx], 2, replace=False))
                gate = two_qudit_gates[rng.integers(0, len(two_qudit_gates))]
                gates.append(gate)
                qudit_indices.append(indices)
            else:
                index = int(rng.choice(index_sets[set_idx]))
                gate = single_qudit_gates[rng.integers(0, len(single_qudit_gates))]
                gates.append(gate)
                qudit_indices.append((index,))

        C = cls(dimensions, gates, qudit_indices)
        C._sanity_check()

        return C

    @classmethod
    def from_tuples(cls, dimensions: DimensionsLike,
                    data: list[GateSpec] | GateSpec) -> Circuit:
        """
        Creates a circuit from a list of (gate, qudit_indices...) tuples.

        Parameters
        ----------
        dimensions : DimensionsLike
            The dimension of each qudit.
        data : list of tuples
            Each tuple contains (Gate, qudit_idx1, qudit_idx2, ...).

        Returns
        -------
        Circuit
            A new Circuit.

        Example
        -------
        >>> Circuit.from_tuples([2, 2], [(GATES.H, 0), (GATES.CX, 0, 1)])
        """
        if isinstance(data, tuple) and isinstance(data[0], Gate):
            data = [data]

        assert isinstance(data, list)

        dimensions = np.asarray(dimensions, dtype=int)
        gates = [d[0] for d in data]
        qudit_indices = [d[1:] for d in data]

        C = cls(dimensions, gates, qudit_indices)
        C._sanity_check()

        return C

    @classmethod
    def from_gates_and_qudits(cls, dimensions: DimensionsLike,
                              gates: list[Gate], qudit_indices: list[tuple[int, ...]]) -> Circuit:
        """
        Creates a circuit from a list of (gate, qudit_indices...) tuples.

        Parameters
        ----------
        dimensions : DimensionsLike
            The dimension of each qudit.
        gates : list[Gate]
            List of Gate objects.
        qudit_indices : list[tuple[int, ...]]
            List of tuples indicating which qudits each gate acts on.
            Must have the same length as gates.

        Returns
        -------
        Circuit
            A new Circuit.

        Example
        -------
        >>> Circuit.from_gates_and_qudits([2, 2], [GATES.H], [(0,)])
        """

        if len(gates) != len(qudit_indices):
            raise ValueError("Gates and qudit indices must have the same length.")

        dimensions = np.asarray(dimensions, dtype=int)
        C = cls(dimensions, gates, qudit_indices)
        C._sanity_check()

        return C

    @classmethod
    def from_depth(cls,
                   depth: int,
                   dimensions: DimensionsLike,
                   gates_set: tuple[Gate] | list[Gate] | set[Gate] | None = None,
                   two_qudit_gate_ratio: float = 0.3,
                   rng: RNGGenerator | None = None) -> Circuit:
        """
        Creates a random circuit with the given depth. Similar to from_random, but the circuit
        is created with layers of gates, acting on each qudit.

        Parameters
        ----------
        depth : int
            Number of gates per qudit in the circuit.
        dimensions : DimensionsLike
            The dimension of each qudit.
        gates_set : tuple[Gate] | list[Gate] | set[Gate] | None, default None
            The set of gates from which to draw. If ``None`` a default gate set is used.
        two_qudit_gate_ratio : float, default 0.3
            Probability of choosing a two-qudit gate vs single-qudit gate.
        rng : numpy.random.Generator or None, default None
            Random number generator. If ``None``, a default generator is used.

        Returns
        -------
        Circuit
            A new random Circuit.
        """

        if gates_set is None:
            single_qudit_gates: list[Gate] = [GATES.H, GATES.S]
            two_qudit_gates: list[Gate] = [GATES.CX, GATES.CZ, GATES.SWAP]
        elif isinstance(gates_set, (tuple, list, set)):
            single_qudit_gates = [gate for gate in gates_set if gate.n_qudits == 1]
            two_qudit_gates = [gate for gate in gates_set if gate.n_qudits == 2]
            if any([gate.n_qudits > 2 for gate in gates_set]):
                warnings.warn("Only single qudit and 2-qudits gates are used to generate the circuit.")
        else:
            raise ValueError("Invalid gates_set type.")

        if two_qudit_gate_ratio < 0.0 or two_qudit_gate_ratio > 1.0:
            raise ValueError(
                f"Invalid two_qudit_gate_ratio, it should be between 0 and 1 (got {two_qudit_gate_ratio}).")

        if rng is None:
            rng = default_rng()

        dimensions = np.asarray(dimensions, dtype=int)
        # generate list of lists of indexes for each dimension
        groups: dict[int, list[int]] = defaultdict(list)
        for i, val in enumerate(dimensions):
            groups[val].append(i)
        two_qudits_gates_index_sets = [v for v in groups.values() if len(v) >= 2]
        num_index_sets = len(two_qudits_gates_index_sets)

        n_qudits = len(dimensions)
        num_gates = n_qudits * depth
        # Divide by two since each gate applies to 2 qudits
        num_two_qudits_gates = int(two_qudit_gate_ratio * num_gates) // 2

        # First assign only 1-qudit gates
        layers: list[dict[tuple[int, ...], Gate]] = []
        for _ in range(depth):
            layer = {}

            for q in range(n_qudits):
                gate = single_qudit_gates[rng.integers(0, len(single_qudit_gates))]
                layer[(q,)] = gate
            layers.append(layer)

        # Distribute num_two_qudits_gates over depth layers randomly,
        # with max max_two_qudits_gates_per_layer per layer.
        if num_two_qudits_gates > 0 and two_qudit_gates:
            for l_idx in rng.choice(range(depth), num_two_qudits_gates):
                layer = layers[int(l_idx)]
                set_idx = rng.integers(0, num_index_sets)
                first_pick = int(set_idx)
                while True:
                    # Restrict to qudits who had at most `depth` 2-qudit gates
                    available_qudits = [q for q in two_qudits_gates_index_sets[set_idx] if (q,) in layer]
                    if len(available_qudits) >= 2:
                        break
                    else:
                        set_idx = (set_idx + 1) % num_index_sets
                        # Ensure this terminates eventually
                        if set_idx == first_pick:
                            # FIXME: instead pick a different layer and raise only if none is available
                            raise ValueError("Could not find a gate layout to satisfy the 2-qudits requirements.")
                indices = tuple(int(idx) for idx in rng.choice(available_qudits, 2, replace=False))
                for idx in indices:
                    layer.pop((idx,))
                gate = two_qudit_gates[rng.integers(0, len(two_qudit_gates))]
                layer[indices] = gate

        gates = [g for layer in layers for g in layer.values()]
        qudit_indices = [idxs for layer in layers for idxs in layer.keys()]

        C = cls(dimensions, gates, qudit_indices)
        C._sanity_check()

        return C

    @classmethod
    def from_string(cls, s: str) -> Circuit:
        """
        Create a Circuit from a JSON string.

        The string should be a JSON object with:
        - "data": list of gate operations, each as [gate_name, [qudit_indices], noise]

        Parameters
        ----------
        s : str
            JSON string representing the circuit.

        Returns
        -------
        Circuit
            The deserialized circuit.

        Example
        -------
        >>> s = '{"data": [["H", [0]], ["CX", [0, 1]]]}'
        >>> circuit = Circuit.from_string(s)
        """
        data = json.loads(s)
        dimensions = data["dimensions"]
        gate_data = data["data"]

        # Map gate names to gate singletons
        gate_map = {
            "H": GATES.H,
            "H_inv": GATES.H_inv,
            "S": GATES.S,
            "S_inv": GATES.S_inv,
            "CX": GATES.CX,
            "CX_inv": GATES.CX_inv,
            "SWAP": GATES.SWAP,
            "CZ": GATES.CZ,
        }

        gates = []
        qudit_indices = []
        noise_model_per_gate = []

        for gate_spec in gate_data:
            gate_name = gate_spec[0]
            indices = tuple(gate_spec[1])

            if gate_name not in gate_map:
                raise ValueError(f"Unknown gate name: {gate_name}")

            gates.append(gate_map[gate_name])
            qudit_indices.append(indices)

            noise_str = gate_spec[2] if len(gate_spec) > 2 else "None"
            noise = NoiseModel.from_string(noise_str)
            noise_model_per_gate.append(noise)

        dimensions = np.asarray(dimensions, dtype=int)
        C = cls(dimensions, gates, qudit_indices).with_noise(noise_model_per_gate)
        C._sanity_check()

        return C

    @classmethod
    def from_file(cls, file_path: str | Path) -> Circuit:
        """
        Create a Circuit from a JSON file.

        Parameters
        ----------
        file_path : str | Path
            Path to the JSON file.

        Returns
        -------
        Circuit
            The deserialized circuit.
        """
        file_path = Path(file_path)
        with open(file_path, 'r') as f:
            return cls.from_string(f.read())

    def with_noise(self, noise_model: NoiseModel | list[NoiseModel | None]) -> Circuit:
        """
        Attach a noise model and return self for chaining.

        Parameters
        ----------
        noise_model : NoiseModel | list[NoiseModel | None]
            A single noise model applied after every gate, or a list of
            per-gate noise models (with ``None`` for noiseless gates).
            A list must have length equal to the number of gates.

        Returns
        -------
        Circuit
            This circuit instance (for method chaining).
        """

        if isinstance(noise_model, NoiseModel):
            self._noise_model_per_gate = [noise_model] * self.n_gates()
        elif isinstance(noise_model, list):
            if len(noise_model) != self.n_gates():
                raise ValueError("Invalid noise input. A list of noise models should have exactly \
                                 one element per gate.")
            self._noise_model_per_gate = noise_model
        else:
            raise ValueError("Invalid noise input: must be either a NoiseModel (set global noise) or \
                             a list of noise models with length equal to the number of gates.")

        return self

    def set_noise(self, noise_model: NoiseModel | None | list[NoiseModel | None]):
        """
        Attach a noise model to this circuit.

        Parameters
        ----------
        noise_model : NoiseModel | None | list[NoiseModel | None]
            A single noise model applied after every gate, ``None`` to
            remove noise, or a list of per-gate noise models.
        """
        if noise_model is None:
            self._noise_model_per_gate = [None] * self.n_gates()
        else:
            self = self.with_noise(noise_model)

    def _sanity_check(self):
        """
        Validate internal consistency of the Circuit.

        Raises
        ------
        ValueError
            If gates and qudit indices are not consistent.
        """
        if len(self._gates) != len(self._qudit_indices):
            raise ValueError(f"gates and qudit_indices must have the same length, "
                             f"got {len(self._gates)} gates and {len(self._qudit_indices)} qudit tuples.")

        if np.any(self.dimensions < DEFAULT_QUDIT_DIMENSION):
            bad_dims = self.dimensions[self.dimensions < DEFAULT_QUDIT_DIMENSION]
            raise ValueError(f"Dimensions {bad_dims} are less than {DEFAULT_QUDIT_DIMENSION}")

        for gate, idxs in zip(self._gates, self._qudit_indices):
            if len(idxs) == 0:
                raise ValueError("Gate cannot act on no qudit.")

            if len(idxs) != gate.n_qudits:
                raise ValueError(f"Gate and qudit indices do not match. Gate acts on {gate.n_qudits} qudits, "
                                 f"but {len(idxs)} qudit indices were provided.")

            if len(idxs) != len(set(idxs)):
                raise ValueError(f"Qudit indices must all differ, got {idxs}.")

            relevant_dimensions = self.dimensions[list(idxs)]
            if not np.all(relevant_dimensions == relevant_dimensions[0]):
                raise ValueError("Gate cannot act on qudits with different dimensions.")

    def add_gate(self, gate: Gate, *qudit_indices: int, noise_model: NoiseModel | None = None):
        """
        Appends a gate acting on the specified qudits.

        Parameters
        ----------
        gate : Gate
            The gate to add.
        qudit_indices : int
            The indices of the qudits the gate acts on.
        """

        if len(qudit_indices) != gate.n_qudits:
            raise ValueError(f"Gate {gate.name} acts on {gate.n_qudits} qudits, "
                             f"but {len(qudit_indices)} indices provided.")

        for idx in qudit_indices:
            if idx < 0 or idx >= len(self.dimensions):
                raise IndexError(f"Qudit index {idx} out of range for circuit with {len(self.dimensions)} qudits.")

        affected_dimensions = [self.dimensions[i] for i in qudit_indices]
        if len(set(affected_dimensions)) != 1:
            raise ValueError(f"Gate must act on qudits with the same dimensions (found {set(affected_dimensions)}).")

        self._gates.append(gate)
        self._qudit_indices.append(tuple(qudit_indices))
        self._noise_model_per_gate.append(noise_model)

    def remove_gate(self, index: int):
        """Removes the gate at the specified index."""
        self._gates.pop(index)
        self._qudit_indices.pop(index)
        self._noise_model_per_gate.pop(index)

    def n_qudits(self) -> int:
        """Returns the number of qudits in the circuit."""
        return len(self.dimensions)

    def n_gates(self) -> int:
        """
        Returns the number of gates in the circuit.
        """
        return len(self.gates)

    def set_use_unitary_cache(self, value: bool):
        """
        Enable or disable caching of the circuit's unitary matrix.

        Parameters
        ----------
        value : bool
            If ``True``, the unitary is cached after the first computation.
        """
        self._use_unitary_cache = value

    @property
    def lcm(self) -> int:
        """Returns the LCM of all qudit dimensions."""
        return int(np.lcm.reduce(self.dimensions))

    def __add__(self, other: Circuit) -> Circuit:
        """Concatenates two circuits."""
        if not isinstance(other, Circuit):
            raise TypeError("Can only add another Circuit object.")

        if not np.array_equal(self.dimensions, other.dimensions):
            raise ValueError("Cannot concatenate circuits with different dimensions.")

        new_gates = self._gates + other._gates
        new_qudits = self._qudit_indices + other._qudit_indices
        new_noise_model_per_gate = self._noise_model_per_gate + other._noise_model_per_gate
        C = Circuit(self.dimensions, new_gates, new_qudits)
        C._noise_model_per_gate = new_noise_model_per_gate

        return C

    def __eq__(self, other: Circuit) -> bool:
        if not isinstance(other, Circuit):
            return False
        if not np.array_equal(self.dimensions, other.dimensions):
            return False
        if len(self._gates) != len(other._gates):
            return False
        for i in range(len(self._gates)):
            if self._gates[i] is not other._gates[i]:  # Compare by identity for singletons
                return False
            if self._qudit_indices[i] != other._qudit_indices[i]:
                return False
            if self._noise_model_per_gate[i] != other._noise_model_per_gate[i]:
                return False
        return True

    def __len__(self) -> int:
        return len(self._gates)

    def __str__(self) -> str:
        lines = [f"Circuit on {self.n_qudits()} qudits (dims={list(self.dimensions)}):"]
        for gate, qudits, noise in zip(self._gates, self._qudit_indices, self._noise_model_per_gate):
            noise_str = " " + noise.__str__() if noise else " "
            lines.append(f"  {gate.name} {qudits}{noise_str}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"Circuit(dimensions={list(self.dimensions)}, n_gates={len(self._gates)})"

    @overload
    def act(self, pauli: PauliString) -> PauliString:
        ...

    @overload
    def act(self, pauli: PauliSum) -> PauliSum:
        ...

    def act(self, pauli: PauliObject) -> PauliObject:
        """
        Apply all gates in the circuit to a Pauli object.

        If a noise model is attached, it is applied after each gate.

        Parameters
        ----------
        pauli : PauliObject
            The Pauli object to transform.

        Returns
        -------
        PauliObject
            The transformed Pauli object after all gates (and noise) are applied.
        """
        for gate, qudits, noise in zip(self._gates, self._qudit_indices, self._noise_model_per_gate):
            pauli = gate.act(pauli, qudits)
            if noise is not None:
                pauli = noise.apply_quantum_trajectory(pauli, qudits)

        return pauli

    @overload
    def act_iter(self, pauli: PauliString) -> Generator[PauliString, None, None]:
        ...

    @overload
    def act_iter(self, pauli: PauliSum) -> Generator[PauliSum, None, None]:
        ...

    def act_iter(self, pauli: PauliObject) -> Generator[PauliObject, None, None]:
        """
        Yield the Pauli object after each gate application.

        If a noise model is attached, it is applied after each gate.

        Parameters
        ----------
        pauli : PauliObject
            The Pauli object to transform.

        Yields
        ------
        PauliObject
            The Pauli object after each successive gate (and noise) application.
        """
        for gate, qudits, noise in zip(self._gates, self._qudit_indices, self._noise_model_per_gate):
            pauli = gate.act(pauli, qudits)
            if noise is not None:
                pauli = noise.apply_quantum_trajectory(pauli, qudits)
            yield pauli

    def act_in_hilbert_space(self, rho: HilbertOperator) -> HilbertOperator:
        """
        Apply all gates in the circuit to a density matrix in Hilbert space.

        Gate unitaries are cached for reuse. If a noise model is attached,
        it is applied after each gate.

        Parameters
        ----------
        rho : HilbertOperator
            The input density matrix.

        Returns
        -------
        HilbertOperator
            The density matrix after all gates (and noise) are applied.
        """
        for gate, qudits, noise in zip(self._gates, self._qudit_indices, self._noise_model_per_gate):
            key = (gate, qudits)
            if key not in self._unitary_cache:
                # FIXME: we should delegate to gate.act_in_hilbert_space.
                # Problem is that it is unclear how to use the cache AND delegate to Gate method.
                U = embed_unitary(gate.local_unitary(self.dimensions[qudits[0]]), qudits, self.dimensions)
                if self._use_unitary_cache:
                    self._unitary_cache[key] = (U, U.conj().T)
            U, U_dag = self._unitary_cache[key]

            rho = U @ rho @ U_dag
            if noise is not None:
                # FIXME: Add cache also for noise gates.
                # Not trivial since noise modle does not know about circuit dimensions.
                # A more consistent way would be to standardize the cache, thus using also the dimensions
                # in the key.
                rho = noise.act_in_hilbert_space(rho, qudits, self.dimensions)

        return rho

    def act_in_hilbert_space_iter(self, rho: HilbertOperator) -> Generator[HilbertOperator, None, None]:
        """
        Yield the density matrix after each gate application in Hilbert space.

        Gate unitaries are cached for reuse. If a noise model is attached,
        it is applied after each gate.

        Parameters
        ----------
        rho : HilbertOperator
            The input density matrix.

        Yields
        ------
        HilbertOperator
            The density matrix after each successive gate (and noise) application.
        """
        for gate, qudits, noise in zip(self._gates, self._qudit_indices, self._noise_model_per_gate):
            key = (gate, qudits)
            if key not in self._unitary_cache:
                U = embed_unitary(gate.local_unitary(self.dimensions[qudits[0]]), qudits, self.dimensions)
                self._unitary_cache[key] = (U, U.conj().T)
            U, U_dag = self._unitary_cache[key]

            rho = U @ rho @ U_dag
            if noise is not None:
                rho = noise.act_in_hilbert_space(rho, qudits, self.dimensions)

            yield rho

    def copy(self) -> Circuit:
        """Returns a shallow copy of the circuit."""
        return Circuit(self.dimensions, self._gates.copy(), self._qudit_indices.copy())

    def _composite_phase_vector(self, F_1: TableauType, F_2: TableauType, h_2: PhasesType) -> PhasesType:
        """
        Faster equivalent of:
            U = [[0,0],[I,0]]
            Uc = F2.T @ U @ F2
            p2 = diag(F1 @ ( (2*triu(Uc)-diag(diag(Uc))) @ F1.T ))
        Uses:
            p2 = sum( (F1 @ Q) * F1, axis=1 )   (one matmul)
        and block construction of Uc.
        """
        F1 = np.ascontiguousarray(F_1, dtype=int)
        F2 = np.ascontiguousarray(F_2, dtype=int)
        h2 = np.asarray(h_2, dtype=int)

        n2 = F1.shape[0]
        assert n2 % 2 == 0 and F1.shape == (n2, n2) and F2.shape == (n2, n2)
        n = n2 // 2

        # ---- Build U_conjugated = F2^T U F2 using blocks (A,B,C,D) ----
        A = F2[:n, :n]
        B = F2[:n, n:]
        C = F2[n:, :n]
        D = F2[n:, n:]

        # Uc = [[C^T A, C^T B],
        #       [D^T A, D^T B]]
        Uc = np.empty((n2, n2), dtype=int)
        Uc[:n, :n] = C.T @ A
        Uc[:n, n:] = C.T @ B
        Uc[n:, :n] = D.T @ A
        Uc[n:, n:] = D.T @ B

        # diag(Uc) without forming np.diag(Uc)
        diag_uc = np.empty(n2, dtype=int)
        diag_uc[:n] = np.sum(C * A, axis=0)
        diag_uc[n:] = np.sum(D * B, axis=0)

        # ---- p1 and p3 ----
        p1 = F1 @ h2
        p3 = F1 @ diag_uc

        # ---- Build Q = 2*triu(Uc) - diag(diag(Uc)) ----
        # After this, Uc == Q
        tril_i, tril_j = np.tril_indices(n2, -1)
        Uc[tril_i, tril_j] = 0          # keep only upper triangle + diag
        Uc *= 2                         # doubles diag too (we will restore)
        np.fill_diagonal(Uc, diag_uc)   # restore diagonal to original (not doubled)

        # ---- p2 with ONE matmul ----
        F1Q = F1 @ Uc                   # only one matmul
        p2 = np.sum(F1Q * F1, axis=1)   # diag(F1 Q F1^T)

        return p1 + p2 - p3

    def composite_gate(self) -> Gate:
        """
        Composes all gates into a single equivalent gate.

        Returns a generic Gate (not a singleton) representing the full circuit transformation.
        """
        n_qudits = self.n_qudits()
        total_symplectic = np.eye(2 * n_qudits, dtype=int)
        lcm = self.lcm
        total_phase_vector = np.zeros(2 * n_qudits, dtype=int)

        for i, (gate, qudits) in enumerate(zip(self._gates, self._qudit_indices)):
            # Get the phase vector for the relevant dimension
            relevant_dim = int(np.lcm.reduce(self.dimensions[list(qudits)]))
            phase_vec = gate.phase_vector(relevant_dim)

            # Embed the local symplectic into the full space
            F, h = embed_symplectic(gate.symplectic, phase_vec, qudits, n_qudits)

            if i == 0:
                total_phase_vector = h
            else:
                total_phase_vector = np.mod(
                    total_phase_vector + self._composite_phase_vector(total_symplectic, F, h),
                    2 * lcm
                )

            total_symplectic = np.mod(total_symplectic @ F.T, lcm)

        total_symplectic = total_symplectic.T
        return _GenericGate('CompositeGate', total_symplectic, total_phase_vector)

    def inverse(self) -> Circuit:
        """Returns the inverse circuit (gates in reverse order, each inverted)."""
        inv_gates = [g.inverse() for g in reversed(self._gates)]
        inv_qudits = list(reversed(self._qudit_indices))
        return Circuit(self.dimensions, inv_gates, inv_qudits)

    def full_symplectic(self) -> TableauType:
        """Returns the full symplectic matrix of the composite gate."""
        return self.composite_gate().symplectic

    def unitary(self) -> HilbertOperator:
        """
        Compute the unitary matrix of the full circuit.

        Returns a sparse matrix of shape (D, D) where D = prod(dimensions).
        Gates are applied in sequence, with each gate's unitary embedded
        into the full Hilbert space.

        For single-qudit gates, the gate's dimension is taken from the target qudit.
        For multi-qudit gates, all target qudits must have the same dimension.

        Returns
        -------
        HilbertOperator
            The unitary matrix of the circuit.

        Raises
        ------
        ValueError
            If a multi-qudit gate acts on qudits with different dimensions.
        """

        D = int(np.prod(self.dimensions))
        U_total = sp.eye(D, format='csr')

        for gate, qudits in zip(self._gates, self._qudit_indices):
            # Get the dimension(s) for this gate
            gate_dims = self.dimensions[list(qudits)]

            if gate.n_qudits > 1:
                # Multi-qudit gate: all qudits must have the same dimension
                if not np.all(gate_dims == gate_dims[0]):
                    raise ValueError(
                        f"Gate {gate.name} acts on qudits with different dimensions {gate_dims}. "
                        "Multi-qudit gates require equal dimensions."
                    )

            d = int(gate_dims[0])
            U_local = gate.local_unitary(d)

            # Embed into the full Hilbert space
            U_embedded = embed_unitary(U_local, list(qudits), self.dimensions)

            # Compose: circuit is applied left-to-right, so U_total = U_embedded @ U_total
            U_total = U_embedded @ U_total

        return U_total

    to_hilbert_space = unitary

    def to_string(self) -> str:
        """
        Serialize the circuit to a JSON string.

        Returns
        -------
        str
            JSON string representation of the circuit.

        Example
        -------
        >>> circuit = Circuit.from_tuples([(GATES.H, 0), (GATES.CX, 0, 1)])
        >>> circuit.to_string()
        '{"dimensions": [2, 3], "data": [["H", [0]], ["CX", [0, 1]]]}'
        """
        gate_data = []
        for gate, qudits, noise in zip(self._gates, self._qudit_indices, self._noise_model_per_gate):
            gate_data.append([gate.name, [int(q) for q in qudits], noise.__str__()])
        return json.dumps({"dimensions": [int(d) for d in self.dimensions], "data": gate_data})

    def save_to_file(self, file_path: str | Path) -> None:
        """
        Save the circuit to a JSON file.

        Parameters
        ----------
        file_path : str | Path
            Path to the output file.

        Example
        -------
        >>> circuit = Circuit.from_tuples([(GATES.H, 0)])
        >>> circuit.save_to_file("my_circuit.json")
        """
        file_path = Path(file_path)
        with open(file_path, 'w') as f:
            f.write(self.to_string())

    def cleanup(self):
        # TODO If two gates are the inverse of each other and next to each other, remove them both. This happens
        # in a few algorithms
        raise NotImplementedError

    def gates_layout(self,
                     with_qudit_indices: bool = False,
                     with_input: PauliSum | None = None,
                     with_output: PauliSum | None = None,
                     wires: str | list[str] | None = None,
                     wrap: bool = True) -> str:
        """
        Returns a visual circuit diagram of the RMB.

        Renders the circuit as ASCII art with gates displayed as boxes
        connected by wires.

        Parameters
        ----------
        with_qudit_indices : bool, default False
            If True, display qudit indices on the left of each wire.
        with_input : PauliSum | None, default None
            If provided, display input phases on the left.
        with_output : PauliSum | None, default None
            If provided, display output phases on the right.
        wires : str | list[str] | None, default None
            If provided, overrides default wires strin. If a list is provided, its length must match circuit n_qudits.
        wrap : bool, default True
            If True, wrap the output to fit the terminal width by splitting
            at gate boundaries.

        Returns
        -------
        str
            A string representation of the circuit diagram.
        """

        def gate_name(gate: Gate) -> str:
            return gate.name.replace("_inv", "*")[:gate_name_len].center(gate_name_len)

        n_qudits = self.n_qudits()
        lines: list[str] = ["" for _ in range(3 * n_qudits)]
        gate_num: list[int] = [0 for _ in range(n_qudits)]

        if wires is None:
            wires = ["="] * n_qudits
        elif isinstance(wires, str):
            wires = wires * n_qudits
        elif len(wires) != n_qudits:
            raise ValueError("Wires list length must match the circuit number of qudits.")

        gate_name_len = 5
        gate_len = gate_name_len + 4

        if with_qudit_indices:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] = " " * 4
                lines[3 * l_idx + 1] = f"{l_idx:>2}: "
                lines[3 * l_idx + 2] = " " * 4
        else:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] = ""
                lines[3 * l_idx + 1] = ""
                lines[3 * l_idx + 2] = ""

        if with_input is None:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += " " * 2 + wires[l_idx] * 2
                lines[3 * l_idx + 2] += " " * 4
        else:
            # Put initial state phase on the left
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += f"{with_input.phases[l_idx]:<2}" + wires[l_idx] * 2
                lines[3 * l_idx + 2] += " " * 4

        for gate, qudit_indices in zip(self.gates, self.qudit_indices):
            if gate.n_qudits == 1:
                l_idx = qudit_indices[0]
                gate_num[l_idx] += 1

                lines[3 * l_idx + 0] += " ┌" + "─" * gate_name_len + "┐ "
                lines[3 * l_idx + 1] += wires[l_idx] + "│" + f"{gate_name(gate)}" + "│" + wires[l_idx]
                lines[3 * l_idx + 2] += " └" + "─" * gate_name_len + "┘ "
            # 2-qudit gate
            else:
                # Get max line length of affected qudits
                max_num_gate_affected_qudits = max(
                    [gate_num[idx] for idx in qudit_indices])

                for l_idx in qudit_indices:
                    while gate_num[l_idx] < max_num_gate_affected_qudits:
                        gate_num[l_idx] += 1
                        lines[3 * l_idx + 0] += " " * gate_len
                        lines[3 * l_idx + 1] += wires[l_idx] * gate_len
                        lines[3 * l_idx + 2] += " " * gate_len

                    gate_num[l_idx] += 1

                    is_top_qudit = l_idx == min(qudit_indices)
                    is_btm_qudit = l_idx == max(qudit_indices)

                    if is_top_qudit:
                        lines[3 * l_idx + 0] += " ┌" + "─" * gate_name_len + "┐ "
                    else:
                        lines[3 * l_idx + 0] += " ┌" + "─" * (gate_name_len // 2) + \
                            "┴" + "─" * (gate_name_len // 2) + "┐ "

                    lines[3 * l_idx + 1] += wires[l_idx] + "│" + f"{gate_name(gate)}" + "│" + wires[l_idx]
                    if is_btm_qudit:
                        lines[3 * l_idx + 2] += " └" + "─" * gate_name_len + "┘ "
                    else:
                        lines[3 * l_idx + 2] += " └" + "─" * (gate_name_len // 2) + \
                            "┬" + "─" * (gate_name_len // 2) + "┘ "

        max_num_gate = max(gate_num)
        for l_idx in range(n_qudits):
            while gate_num[l_idx] < max_num_gate:
                gate_num[l_idx] += 1
                lines[3 * l_idx + 0] += " " * gate_len
                lines[3 * l_idx + 1] += wires[l_idx] * gate_len
                lines[3 * l_idx + 2] += " " * gate_len

        if with_output is None:
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += wires[l_idx] * 2 + " " * 2
                lines[3 * l_idx + 2] += " " * 4
        else:
            # Put final state phase on the right
            for l_idx in range(n_qudits):
                lines[3 * l_idx + 0] += " " * 4
                lines[3 * l_idx + 1] += wires[l_idx] * 2 + " " + f"{with_output.phases[l_idx]}"
                lines[3 * l_idx + 2] += " " * 4

        if not wrap:
            return "\n".join(lines)

        # Wrap output to fit terminal width
        import shutil
        import re
        term_width = shutil.get_terminal_size().columns
        line_len = len(lines[0])  # Top border has no ANSI codes

        if line_len <= term_width:
            return "\n".join(lines)

        # Strip ANSI codes from wire lines for correct slicing
        ansi_pattern = re.compile(r'\033\[[0-9;]*m')
        plain_wire_lines = {
            l_idx: ansi_pattern.sub('', lines[3 * l_idx + 1])
            for l_idx in range(n_qudits)
        }

        sections = []
        pos = 0

        while pos < line_len:
            if pos == 0:
                # First section includes prefix from input state
                prefix_len = 8 if with_qudit_indices else 4
                gates_per_section = max(1, (term_width - prefix_len - 2) // gate_len)
                end = min(prefix_len + gates_per_section * gate_len, line_len)
                section_lines = []
                for l_idx in range(n_qudits):
                    section_lines.append(lines[3 * l_idx + 0][pos:end])
                    wire_slice = plain_wire_lines[l_idx][pos:end]
                    wire_slice = wire_slice.replace("=", wires[l_idx])
                    section_lines.append(wire_slice)
                    section_lines.append(lines[3 * l_idx + 2][pos:end])
            else:
                # Continuation sections: add wire prefix for visual continuity
                extra_prefix_side = " " * 4 if with_qudit_indices else ""
                prefix_len = len(extra_prefix_side)
                gates_per_section = max(1, (term_width - prefix_len - 2) // gate_len)
                end = min(pos + gates_per_section * gate_len, line_len)
                section_lines = []

                for l_idx in range(n_qudits):
                    extra_prefix_central = f"{l_idx:>2}: " if with_qudit_indices else ""
                    section_lines.append(extra_prefix_side + "  " + lines[3 * l_idx + 0][pos:end])
                    wire_slice = plain_wire_lines[l_idx][pos:end]
                    wire_slice = wire_slice.replace("=", wires[l_idx])
                    section_lines.append(extra_prefix_central + wires[l_idx] * 2 + wire_slice)
                    section_lines.append(extra_prefix_side + "  " + lines[3 * l_idx + 2][pos:end])

            sections.append("\n".join(section_lines))
            pos = end

        return "\n\n\n".join(sections)
