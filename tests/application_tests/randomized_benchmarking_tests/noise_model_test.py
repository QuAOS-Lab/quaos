import numpy as np
import pytest
from numpy.random import default_rng

from sympleq.applications.randomized_benchmarking.noise_model import (
    NoiseModel,
    Noiseless,
    DephasingNoise,
    DepolarizingNoise,
)
from sympleq.applications.randomized_benchmarking.RMB import RMB, pauli_to_rho
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import PHASE, Hadamard
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION


# All noise model classes to test
NOISE_MODELS = [
    Noiseless(),
    DephasingNoise(error_rate=0.1),
    DepolarizingNoise(error_rate=0.1),
]


class TestNoiseModelDimensions:
    """Test that output dimensions are consistent with n_kraus_operators and n_qudits."""

    @pytest.mark.parametrize("model", NOISE_MODELS)
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_probabilities_and_operators_length(self, model: NoiseModel, n_qudits: int):
        """Verify that probabilities and Kraus operators have length n_kraus^n_qudits."""
        n_kraus = model.n_kraus_operators()
        expected_length = n_kraus ** n_qudits

        # Test probabilities length
        probs = model.kraus_probabilities(n_qudits)
        assert len(probs) == expected_length, \
            f"Expected {expected_length} probabilities, got {len(probs)}"

        # Test Kraus operators length
        dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits
        qudit_indices = list(range(n_qudits))
        operators = model.kraus_operators(dimensions, qudit_indices)
        assert len(operators) == expected_length, \
            f"Expected {expected_length} Kraus operators, got {len(operators)}"

    @pytest.mark.parametrize("model", NOISE_MODELS)
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_process_matrix_shape(self, model: NoiseModel, n_qudits: int):
        """Verify process matrix has shape (n^M, n^M)."""
        n_kraus = model.n_kraus_operators()
        expected_size = n_kraus ** n_qudits

        matrix = model.process_matrix(n_qudits)

        assert matrix.shape == (expected_size, expected_size), \
            f"Expected shape ({expected_size}, {expected_size}), got {matrix.shape}"


class TestNoiseModelProbabilities:
    """Test probability properties: sum to 1, non-negative."""

    @pytest.mark.parametrize("model", NOISE_MODELS)
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_probabilities_sum_to_one(self, model: NoiseModel, n_qudits: int):
        probs = model.kraus_probabilities(n_qudits)
        assert np.isclose(probs.sum(), 1.0), \
            f"Probabilities sum to {probs.sum()}, expected 1.0"

    @pytest.mark.parametrize("model", NOISE_MODELS)
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_probabilities_non_negative(self, model: NoiseModel, n_qudits: int):
        probs = model.kraus_probabilities(n_qudits)
        assert np.all(probs >= 0), "Found negative probabilities"

    @pytest.mark.parametrize("error_rate", [0.0, 0.1, 0.5, 1.0])
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_dephasing_probabilities_sum_to_one_various_rates(self, error_rate: float, n_qudits: int):
        model = DephasingNoise(error_rate=error_rate)
        probs = model.kraus_probabilities(n_qudits)
        assert np.isclose(probs.sum(), 1.0)

    @pytest.mark.parametrize("error_rate", [0.0, 0.1, 0.5, 1.0])
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_depolarizing_probabilities_sum_to_one_various_rates(self, error_rate: float, n_qudits: int):
        model = DepolarizingNoise(error_rate=error_rate)
        probs = model.kraus_probabilities(n_qudits)
        assert np.isclose(probs.sum(), 1.0)


class TestNoiseModelProcessMatrix:
    """Test process matrix properties."""

    @pytest.mark.parametrize("model", [DephasingNoise(0.1), DepolarizingNoise(0.1)])
    @pytest.mark.parametrize("n_qudits", [1, 2])
    def test_process_matrix_diagonal(self, model: NoiseModel, n_qudits: int):
        """These noise models should have diagonal process matrices."""
        matrix = model.process_matrix(n_qudits)

        off_diagonal = matrix - np.diag(np.diag(matrix))
        assert np.allclose(off_diagonal, 0), "Process matrix has non-zero off-diagonal elements"

    @pytest.mark.parametrize("model", NOISE_MODELS)
    @pytest.mark.parametrize("n_qudits", [1, 2])
    def test_process_matrix_trace_equals_one(self, model: NoiseModel, n_qudits: int):
        """For diagonal noise models, trace of process matrix should be 1."""
        matrix = model.process_matrix(n_qudits)
        assert np.isclose(np.trace(matrix), 1.0), \
            f"Trace = {np.trace(matrix)}, expected 1.0"


class TestNoiseModelEdgeCases:
    """Test edge cases and error handling."""

    def test_dephasing_invalid_error_rate_negative(self):
        with pytest.raises(ValueError):
            DephasingNoise(error_rate=-0.1)

    def test_dephasing_invalid_error_rate_above_one(self):
        with pytest.raises(ValueError):
            DephasingNoise(error_rate=1.1)

    def test_depolarizing_invalid_error_rate_negative(self):
        with pytest.raises(ValueError):
            DepolarizingNoise(error_rate=-0.1)

    def test_depolarizing_invalid_error_rate_above_one(self):
        with pytest.raises(ValueError):
            DepolarizingNoise(error_rate=1.1)

    def test_dephasing_zero_error_rate(self):
        """Zero error rate should give probability 1 for identity Kraus operator."""
        model = DephasingNoise(error_rate=0.0)
        probs = model.kraus_probabilities(1)

        assert np.isclose(probs[0], 1.0)
        assert np.isclose(probs[1], 0.0)

    def test_depolarizing_zero_error_rate(self):
        """Zero error rate should give probability 1 for identity Kraus operator."""
        model = DepolarizingNoise(error_rate=0.0)
        probs = model.kraus_probabilities(1)

        assert np.isclose(probs[0], 1.0)
        assert np.allclose(probs[1:], 0.0)


class TestRMB:
    """Tests for RMB class, moved from RMB.py."""

    def test_simple_circuit_with_dephasing(self):
        """Test RMB with a simple circuit and dephasing noise."""
        n_qudits = 2
        dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits
        dimension = dimensions[0]

        circuit = Circuit(dimensions, [Hadamard(0, dimension), PHASE(0, dimension)])
        rmb = RMB(circuit, False, 0.0, 0.0, [DephasingNoise(0.1)], default_rng(42))

        # Basic sanity checks
        assert rmb.n_qudits == n_qudits
        assert rmb.initial_state is not None

        # Check that rho has trace 1
        rho = pauli_to_rho(rmb.initial_state)
        assert np.isclose(np.trace(rho), 1.0)

    def test_random_rmb_with_multiple_noise_models(self):
        """Test RMB.from_random with multiple noise models."""
        n_qudits = 3
        gate_density = 2.5
        dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits

        rmb = RMB.from_random(
            dimensions,
            gate_density,
            noise_models=[DephasingNoise(0.25), DepolarizingNoise(0.15)],
            with_random_elimination=0.0,
            rng=default_rng(42)
        )

        assert rmb.n_qudits == n_qudits

        # Test rho_exact returns valid density matrix
        rho_exact = rmb.rho_exact()
        assert np.isclose(np.trace(rho_exact), 1.0), \
            f"rho_exact trace = {np.trace(rho_exact)}, expected 1.0"

    def test_rho_exact_vs_rho_average_convergence(self):
        """Test that rho_average converges to rho_exact with enough samples."""
        n_qudits = 2
        dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits

        rmb = RMB.from_random(
            dimensions,
            gate_density=1.0,
            noise_models=[DephasingNoise(0.1)],
            rng=default_rng(42)
        )

        rho_exact = rmb.rho_exact()
        rho_average = rmb.rho_average(n_runs=1000)

        # Check traces are both 1
        assert np.isclose(np.trace(rho_exact), 1.0)
        assert np.isclose(np.trace(rho_average), 1.0)

        # Check convergence (with some tolerance for stochastic sampling)
        error = np.sum(np.abs(rho_exact - rho_average))
        assert error < 0.5, f"Error between rho_exact and rho_average: {error}"
