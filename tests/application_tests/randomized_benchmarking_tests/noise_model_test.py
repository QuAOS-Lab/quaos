import numpy as np
import pytest

from sympleq.applications.randomized_benchmarking.noise_model import (
    NoiseModel,
    Noiseless,
    DephasingNoise,
    DepolarizingNoise,
    CompositeNoise,
)
from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION


# All noise model classes to test
NOISE_MODELS = [
    Noiseless(),
    DephasingNoise(error_rate=0.1),
    DepolarizingNoise(error_rate=0.1),
    CompositeNoise.from_noise_models([DephasingNoise(0.1), DepolarizingNoise(0.1)]),
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

    @pytest.mark.parametrize("n_qudits", [2, 4, 6])
    def test_rho_exact_vs_rho_average_convergence(self, n_qudits: int):
        """Test that rho_average converges to rho_exact with enough samples."""
        dimensions = [DEFAULT_QUDIT_DIMENSION] * n_qudits
        noise_model = CompositeNoise.from_noise_models([DephasingNoise(0.2), DepolarizingNoise(0.0)])

        rmb = RMB.from_random(
            dimensions,
            gate_density=1.5,
            noise_model=noise_model
        )

        n_runs = 1000
        rho_average = rmb.rho_average(n_runs=n_runs)
        rho_exact = rmb.rho_exact()

        # Check traces are both 1
        assert np.isclose(np.trace(rho_exact), 1.0)
        assert np.isclose(np.trace(rho_average), 1.0)

        # Check convergence: error scales with state space size
        max_error = np.max(np.abs(rho_exact - rho_average))
        tolerance = 30 * n_qudits / n_runs
        assert max_error < tolerance, f"Max error {max_error} exceeds tolerance {tolerance}"


class TestCompositeNoise:
    """Tests specific to CompositeNoise class."""

    def test_from_noise_models_empty_list_raises(self):
        """Empty noise_models list should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            CompositeNoise.from_noise_models([])

    def test_from_noise_models_single_model(self):
        """CompositeNoise with single model should work."""
        model = CompositeNoise.from_noise_models([DephasingNoise(0.1)])
        assert model.n_kraus_operators() == 2  # I and Z

    def test_n_kraus_operators_combines_unique(self):
        """n_kraus_operators should return count of unique Paulis."""
        # Dephasing has I, Z (2 operators)
        # Depolarizing has I, X, Y, Z (4 operators)
        # Combined should have I, X, Y, Z (4 unique operators)
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.1),
            DepolarizingNoise(0.1)
        ])
        assert composite.n_kraus_operators() == 4

    def test_kraus_operators_length_matches_n_kraus(self):
        """kraus_operators should return n_kraus_operators() operators."""
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.1),
            DepolarizingNoise(0.1)
        ])
        operators = composite.kraus_operators([2], [0])
        assert len(operators) == composite.n_kraus_operators()

    def test_probabilities_sum_to_one(self):
        """Combined probabilities should be normalized to sum to 1."""
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.1),
            DepolarizingNoise(0.1)
        ])
        probs = composite.kraus_probabilities(1)
        assert np.isclose(probs.sum(), 1.0)

    def test_probabilities_non_negative(self):
        """All probabilities should be non-negative."""
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.5),
            DepolarizingNoise(0.5)
        ])
        probs = composite.kraus_probabilities(1)
        assert np.all(probs >= 0)

    def test_weights_are_averaged(self):
        """Weights for same Pauli should be averaged across models."""
        # Two identical models should give same result as one
        single = DephasingNoise(0.1)
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.1),
            DephasingNoise(0.1)
        ])

        single_ops = single.kraus_operators([2], [0], weighted=True)
        composite_ops = composite.kraus_operators([2], [0], weighted=True)

        # Same operators, same weights (since averaging identical models)
        assert len(single_ops) == len(composite_ops)
        for s_op, c_op in zip(single_ops, composite_ops):
            assert np.isclose(s_op.weights[0], c_op.weights[0])
