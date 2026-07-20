"""
Noise models for quantum channels, based on the Kraus operator formalism.

Reference: "Enhancing Quantum Computation via Superposition of Quantum Gates",
Phys. Rev. A 108, 062604 (2023). DOI: 10.1103/PhysRevA.108.062604

Overview
--------
A quantum noise channel transforms a density matrix PS as:

    PS_out = Σ_i K_i PS K†_i

where each Kraus operator K_i is a linear combination of Paulis:

    K_i = Σ_j α_{ij} σ_j.

The probability of each Kraus operator (used for quantum trajectory sampling) is:

    p_i = Σ_j |α_{ij}|² = ||K_i||²

We take two simplifying assumptions:
**Uncorrelated noise**
If the noise model is defined on a single qudit, when acting on a multiple-qudits gate it acts
on each qudit independently. It is still possible to define correlated noise models ab initio.

**Clifford noise**
For now, we assume that Kraus operators are Clifford. This simplifies greatly how they act on PauliSums
as they act as Gate, one for each trajectory (see apply_quantum_trajectory).
"""

from __future__ import annotations
from abc import ABC, abstractmethod
import itertools
import re
import warnings
import numpy as np
from numpy.random import Generator as RNGGenerator, default_rng

from sympleq.core.circuits.gates import GATES, Gate
from sympleq.core.circuits.utils import embed_unitary
from sympleq.core.paulis._typing import DimensionsType, HilbertOperator
from sympleq.core.paulis.pauli_object import PauliObject
from sympleq.core.paulis.pauli_sum import PauliSum


class NoiseModel(ABC):
    """
    Abstract base class for quantum noise models.

    Subclasses must implement:
    - n_qudits(): number of qudits the noise model acts on
    - kraus_gates(): list of Gates corresponding to Kraus operators
    - kraus_probabilities(): probabilities p_i = ||K_i||²
    """
    def __init__(self, rng: RNGGenerator) -> None:
        self.rng = rng

    @abstractmethod
    def n_qudits(self) -> int:
        """
        Return the number of qudits this noise model acts on.

        Returns
        -------
        int
            Number of qudits.
        """
        pass

    @abstractmethod
    def kraus_gates(self) -> list[Gate]:
        """
        Return the list of Gates corresponding to each Kraus operator.

        Returns
        -------
        list[Gate]
            Gates G_i such that the Kraus action is ps_out = G_i.act(ps, qudit_indices).
        """
        pass

    @abstractmethod
    def kraus_probabilities(self) -> list[float]:
        """
        Return the Kraus probabilities for the noise channel.

        p_i = sum_j |alpha_{ij}|^2, i.e. the squared norm of each Kraus operator's
        coefficients. The probabilities sum to 1.

        Returns
        -------
        list[float]
            Probabilities p_i for each Kraus operator.
        """
        pass

    def __str__(self) -> str:
        return self.__class__.__name__

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, NoiseModel):
            return NotImplemented
        return str(self) == str(other)

    def __hash__(self) -> int:
        return hash(str(self))

    def to_dict(self) -> dict:
        """
        Serialize this noise model to a JSON-friendly dict.

        Returns
        -------
        dict
            A dict with a ``"type"`` key naming the concrete subclass and
            additional fields needed to rebuild the instance via
            :meth:`from_dict`.
        """
        raise NotImplementedError(
            f"to_dict is not implemented for {type(self).__name__}.")

    @classmethod
    def from_dict(cls, payload: dict, rng: RNGGenerator | None = None) -> NoiseModel:
        """
        Rebuild a noise model from a dict produced by :meth:`to_dict`.

        Parameters
        ----------
        payload : dict
            Dict with a ``"type"`` key identifying the concrete subclass.
        rng : numpy.random.Generator | None
            RNG to attach to the rebuilt model. ``None`` uses ``default_rng()``.
        """
        if rng is None:
            rng = default_rng()
        type_name = payload["type"]
        if type_name == "Noiseless":
            return Noiseless()
        if type_name == "DephasingNoise":
            return DephasingNoise(error_rate=payload["error_rate"], rng=rng)
        if type_name == "DepolarizingNoise":
            return DepolarizingNoise(error_rate=payload["error_rate"], rng=rng)
        if type_name == "CompositeNoise":
            inner = [cls.from_dict(p, rng=rng) for p in payload["noise_models"]]
            return CompositeNoise(inner, rng=rng)
        if type_name == "GenericNoise":
            gates = [getattr(GATES, name) for name in payload["gates"]]
            probabilities = [float(p) for p in payload["probabilities"]]
            return GenericNoise(probabilities, gates, rng=rng)
        raise ValueError(f"Unknown noise model type: {type_name}")

    @classmethod
    def from_string(cls, s: str) -> NoiseModel | None:
        if s == "None":
            return None
        if s == "Noiseless":
            return Noiseless()

        m = re.match(r"CompositeNoise\(\[(.+)\]\)$", s)
        if m is not None:
            inner = m.group(1)
            # Split on ", " but only at the top level (not inside nested parens)
            models = []
            depth = 0
            current = ""
            for ch in inner:
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                if ch == ',' and depth == 0:
                    models.append(cls.from_string(current.strip()))
                    current = ""
                else:
                    current += ch
            if current.strip():
                models.append(cls.from_string(current.strip()))
            return CompositeNoise.from_noise_models(models)

        m = re.match(r"(\w+)\(error_rate=([\d.]+)\)", s)
        if m is None:
            raise ValueError(f"Cannot parse noise model from string: {s}")

        name, error_rate = m.group(1), float(m.group(2))
        if name == "DephasingNoise":
            return DephasingNoise(error_rate)
        elif name == "DepolarizingNoise":
            return DepolarizingNoise(error_rate)
        else:
            raise ValueError(f"Unknown noise model: {name}")

    def apply_quantum_trajectory(self, pauli_sum: PauliObject, qudits: tuple[int, ...]) -> PauliObject:
        """
        Sample a single Kraus operator and apply it to the Pauli sum.

        Uses the quantum trajectory method: a single Kraus operator is sampled
        according to the Kraus probabilities and applied to the Pauli sum.

        Parameters
        ----------
        pauli_sum : PauliObject
            The Pauli object to apply noise to.
        qudits : tuple[int, ...]
            Qudit indices on which the noise acts.

        Returns
        -------
        PauliObject
            The transformed Pauli object after applying the sampled Kraus operator.
        """
        n_qudits = len(qudits)
        if n_qudits > 2:
            warnings.warn(f"Warning: Noise not implemented for gates acting on more than 2 qudits (got {n_qudits}).")
            return pauli_sum

        kraus_n_qudits = self.n_qudits()
        if n_qudits < kraus_n_qudits:
            return pauli_sum

        # Get probabilities to select one possible quantum trajectory
        # Given probabilities [p0, p1, p2, p3], cumsum gives [p0, p0+p1, p0+p1+p2, 1.0].
        # This allows O(log n) sampling via searchsorted with a uniform random number.
        probs = np.cumsum(self.kraus_probabilities())

        if kraus_n_qudits == n_qudits:
            idx = int(np.searchsorted(probs, self.rng.random()))
            kraus_gate = self.kraus_gates()[idx]
            pauli_sum = kraus_gate.act(pauli_sum, qudits)
        elif kraus_n_qudits == 1:
            # Act independently on each qudit
            for qudit in qudits:
                idx = int(np.searchsorted(probs, self.rng.random()))
                kraus_gate = self.kraus_gates()[idx]
                pauli_sum = kraus_gate.act(pauli_sum, qudit)

        # FIXME: handle generic case, else:...
        # this is for kraus_n_qudits > 1, e.g. 2-qudits noise on 3-qudits gate.

        return pauli_sum

    def act_in_hilbert_space(self, rho: HilbertOperator,
                             qudits: tuple[int, ...], dimensions: DimensionsType) -> HilbertOperator:
        """
        Apply the full noise channel to a density matrix.

        Computes rho_out = sum_i K_i rho K_i^dagger by iterating over all
        Kraus operators (or their tensor products for multi-qudit gates).

        Parameters
        ----------
        rho : HilbertOperator
            The input density matrix.
        qudits : tuple[int, ...]
            Qudit indices on which the noise acts.
        dimensions : DimensionsType
            Local Hilbert space dimensions for each qudit.

        Returns
        -------
        HilbertOperator
            The density matrix after applying the noise channel.
        """

        def _apply_unitary(output_rho: HilbertOperator | None, unitary: HilbertOperator,
                           probability: float) -> HilbertOperator:
            if output_rho is None:
                output_rho = probability * (unitary @ rho @ unitary.conjugate().transpose())
            else:
                output_rho += probability * (unitary @ rho @ unitary.conjugate().transpose())
            assert output_rho is not None
            return output_rho

        n_qudits = len(qudits)
        if n_qudits == 0:
            raise ValueError("Noise model must act on at least one qudit.")

        dimension = dimensions[qudits[0]]
        if not np.all(dimensions[qudits] == dimension):
            raise ValueError(f"Noise model must act on qudits with equal dimensions ( got {dimensions[qudits]}).")

        kraus_n_qudits = self.n_qudits()
        if n_qudits < kraus_n_qudits:
            return rho

        output_rho: HilbertOperator | None = None

        # E.g.: single qudit noise on single-qudit gate
        cum_prob = 0
        if n_qudits == kraus_n_qudits:
            unitaries = [embed_unitary(kraus_gate.local_unitary(dimension), qudits, dimensions)
                         for kraus_gate in self.kraus_gates()]
            return np.sum([probability * (unitary @ rho @ unitary.conjugate().transpose())
                           for unitary, probability in zip(unitaries, self.kraus_probabilities())])

        # E.g.: single qudit noise on multiple-qudits gate.
        # Act independently on each qudit.
        if kraus_n_qudits == 1:
            # Get all gates combinations
            all_gates_combinations = list(itertools.product(self.kraus_gates(), repeat=n_qudits))
            # Pray the order is correct, like, do it
            all_probabilities_combinations = list(itertools.product(self.kraus_probabilities(), repeat=n_qudits))
            for kraus_gates, probabilities in zip(all_gates_combinations, all_probabilities_combinations):
                unitary = embed_unitary(kraus_gates[0].local_unitary(dimension), (qudits[0],), dimensions)
                for kraus_gate, qudit in zip(kraus_gates[1:], qudits[1:]):
                    unitary = embed_unitary(kraus_gate.local_unitary(dimension), (qudit,), dimensions) @ unitary

                probability = float(np.prod(probabilities))
                output_rho = _apply_unitary(output_rho, unitary, probability)
                cum_prob += probability

            assert np.abs(cum_prob - 1.0) < 10**(-10), f"cum prob {cum_prob}"
            assert output_rho is not None
            return output_rho

        # FIXME: extend this for self.n_qudits > 1, e.g. 2-qudits noise on 3-qudits gate.

        return rho


class Noiseless(NoiseModel):
    def __init__(self) -> None:
        super().__init__(rng=default_rng())
        self._probabilities = [1.0]

    def n_qudits(self) -> int:
        return 1

    def kraus_gates(self) -> list[Gate]:
        return [GATES.Id]

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities

    def to_dict(self) -> dict:
        return {"type": "Noiseless"}

    def apply_quantum_trajectory(self, pauli_sum: PauliSum, qudits: tuple[int, ...]) -> PauliSum:
        return pauli_sum

    def act_in_hilbert_space(self, rho: HilbertOperator,
                             qudits: tuple[int, ...], dimensions: DimensionsType) -> HilbertOperator:
        return rho


class DephasingNoise(NoiseModel):
    def __str__(self) -> str:
        return f"DephasingNoise(error_rate={1.0 - self.p0:.4f})"

    def __init__(self, error_rate: float, rng: RNGGenerator | None = None) -> None:
        # The dephasing channel has tow Klaus operators: the identity and Z.
        # A single parameter p0 models the noise probability,
        # representing the probability of having no error.
        # Kraus operators:
        # K0 = sqrt(p0) 1
        # K1 = sqrt(1 − p0) Z
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        self._probabilities = [self.p0, 1.0 - self.p0]

        if rng is None:
            rng = default_rng()
        super().__init__(rng)

    def n_qudits(self) -> int:
        return 1

    def kraus_gates(self) -> list[Gate]:
        return [GATES.Id, GATES.Z]

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities

    def to_dict(self) -> dict:
        return {"type": "DephasingNoise", "error_rate": 1.0 - self.p0}


class DepolarizingNoise(NoiseModel):
    def __str__(self) -> str:
        return f"DepolarizingNoise(error_rate={(1.0 - self.p0):.4f})"

    def __init__(self, error_rate: float, rng: RNGGenerator | None = None) -> None:
        # The depolarizing channel has four Klaus operators: the 4 paulis.
        # A single parameter p0 models the noise probability,
        # representing the probability of having no error.
        # Kraus operators:
        # K0 = sqrt(p0) 1
        # Ki =  sqrt(1 − p0/3) σi
        if error_rate > 1.0 or error_rate < 0.0:
            raise ValueError(f"Error rate should be between 0.0 and 1.0 (got {error_rate}).")
        self.p0 = 1.0 - error_rate
        self._probabilities = [self.p0, (1.0 - self.p0) / 3, (1.0 - self.p0) / 3, (1.0 - self.p0) / 3]

        if rng is None:
            rng = default_rng()
        super().__init__(rng)

    def n_qudits(self) -> int:
        return 1

    def kraus_gates(self) -> list[Gate]:
        return [GATES.Id, GATES.X, GATES.Y, GATES.Z]

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities

    def to_dict(self) -> dict:
        return {"type": "DepolarizingNoise", "error_rate": 1.0 - self.p0}


class CompositeNoise(NoiseModel):
    """
    A noise model that combines multiple noise models by merging their Kraus operators.

    Operators representing the same Pauli are combined by summing their weights
    and normalizing by the total number of models.
    """

    def __str__(self) -> str:
        models_str = ", ".join(str(m) for m in self.noise_models)
        return f"CompositeNoise([{models_str}])"

    def __init__(self, noise_models: list[NoiseModel], rng: RNGGenerator) -> None:
        self.noise_models = noise_models
        self._n_models = len(noise_models)
        combined_gates: dict[Gate, float] = self._combine_gates()
        self._probabilities = list(combined_gates.values())
        self._gates = list(combined_gates.keys())

        super().__init__(rng)

    @classmethod
    def from_noise_models(cls, noise_models: list[NoiseModel],
                          rng: RNGGenerator | None = None) -> CompositeNoise:
        if not noise_models:
            raise ValueError("noise_models list cannot be empty")
        if rng is None:
            rng = default_rng()
        return cls(noise_models, rng)

    def _combine_gates(self, compound_identity: bool = True) -> dict[Gate, float]:
        """
        Collect and combine Kraus operators from all underlying models.

        Returns a dict mapping gates to combined probabilities:
        Operators with the same gates have their probability summed and normalized.
        """
        gates_probabilities: dict[Gate, float] = {}

        for model in self.noise_models:
            for gate, probability in zip(model.kraus_gates(), model.kraus_probabilities()):
                if not compound_identity and gate == GATES.Id:
                    continue

                # If operator is unity, don't add to the weight.
                # The noiseless channel will have a weight equal to 1 - Sum other_channels
                if gate in gates_probabilities:
                    gates_probabilities[gate] += probability
                else:
                    gates_probabilities[gate] = probability

        # Normalize weights so probabilities (|w|^2) sum to 1
        norm = sum(gates_probabilities.values())
        for gate in gates_probabilities:
            gates_probabilities[gate] /= norm

        return gates_probabilities

    def n_qudits(self) -> int:
        return self.noise_models[0].n_qudits()

    def kraus_gates(self) -> list[Gate]:
        return self._gates

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities

    def to_dict(self) -> dict:
        return {
            "type": "CompositeNoise",
            "noise_models": [m.to_dict() for m in self.noise_models],
        }


class GenericNoise(NoiseModel):
    """
    A noise model defined by explicit Kraus gate–probability pairs.

    Unlike :class:`CompositeNoise`, which derives its operators by merging
    other noise models, ``GenericNoise`` accepts the gates and their
    probabilities directly.
    """

    def __str__(self) -> str:
        pairs = ", ".join(f"({p:.4f}, {g.name})" for p, g in zip(self._probabilities, self._gates))
        return f"GenericNoise([{pairs}])"

    def __init__(self, probabilities: list[float], gates: list[Gate], rng: RNGGenerator) -> None:
        """
        Initialize a generic noise model.

        Parameters
        ----------
        probabilities : list[float]
            Kraus probabilities for each gate. Must sum to 1.
        gates : list[Gate]
            Kraus gates corresponding to each probability.
        rng : numpy.random.Generator
            Random number generator for trajectory sampling.

        Raises
        ------
        ValueError
            If ``probabilities`` and ``gates`` have different lengths.
        """
        if len(probabilities) != len(gates):
            raise ValueError(
                f"Probabilities and gates lists must have the same length (got {len(probabilities)} and {len(gates)}).")
        self._probabilities = probabilities
        self._gates = gates
        self._n_qudits = max(g.n_qudits for g in self._gates)

        super().__init__(rng)

    @classmethod
    def from_probabilities_and_gates(cls, probabilities: list[float], gates: list[Gate],
                                     rng: RNGGenerator | None = None) -> GenericNoise:
        """
        Create a GenericNoise from probabilities and gates.

        Parameters
        ----------
        probabilities : list[float]
            Kraus probabilities for each gate.
        gates : list[Gate]
            Kraus gates corresponding to each probability.
        rng : numpy.random.Generator or None, optional
            Random number generator. If ``None``, a default is used.

        Returns
        -------
        GenericNoise
            A new GenericNoise instance.
        """
        if rng is None:
            rng = default_rng()
        return cls(probabilities, gates, rng)

    @classmethod
    def from_paulis(cls, probabilities: list[float], rng: RNGGenerator | None = None) -> GenericNoise:
        """
        Create a GenericNoise from probabilities using the Pauli gates as default gates set.
        probabilities[i] refer to Pauli X, Y, and Z for i=1, 2, and 3, respectively

        Parameters
        ----------
        probabilities : list[float]
            Kraus probabilities for each gate.
        rng : numpy.random.Generator or None, optional
            Random number generator. If ``None``, a default is used.

        Returns
        -------
        GenericNoise
            A new GenericNoise instance.
        """
        if rng is None:
            rng = default_rng()

        p0 = 1.0 - sum(probabilities)
        probabilities = [p0, probabilities[0], probabilities[1], probabilities[2]]
        gates = [GATES.Id, GATES.X, GATES.Y, GATES.Z]
        return cls(probabilities, gates, rng)

    def n_qudits(self) -> int:
        return self._n_qudits

    def kraus_gates(self) -> list[Gate]:
        return self._gates

    def kraus_probabilities(self) -> list[float]:
        return self._probabilities

    def to_dict(self) -> dict:
        return {
            "type": "GenericNoise",
            "probabilities": list(self._probabilities),
            "gates": [g.name for g in self._gates],
        }
