import numpy as np
import pytest

from sympleq.core.noise.noise_model import (
    NoiseModel,
    Noiseless,
    DephasingNoise,
    DepolarizingNoise,
    CompositeNoise,
)

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
        """Verify that probabilities and Kraus gates have the same length."""
        assert len(model.kraus_probabilities()) == len(model.kraus_gates())


class TestNoiseModelProbabilities:
    """Test probability properties: sum to 1, non-negative."""

    @pytest.mark.parametrize("model", NOISE_MODELS)
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_probabilities_consistency(self, model: NoiseModel, n_qudits: int):
        probs = model.kraus_probabilities()
        assert all([0 <= p <= 1 for p in probs])
        assert np.isclose(sum(probs), 1.0)

    @pytest.mark.parametrize("error_rate", [0.0, 0.1, 0.5, 1.0])
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_dephasing_probabilities_consistency_various_rates(self, error_rate: float, n_qudits: int):
        model = DephasingNoise(error_rate=error_rate)
        probs = model.kraus_probabilities()
        assert all([0 <= p <= 1 for p in probs])
        assert np.isclose(sum(probs), 1.0)

    @pytest.mark.parametrize("error_rate", [0.0, 0.1, 0.5, 1.0])
    @pytest.mark.parametrize("n_qudits", [1, 2, 3])
    def test_depolarizing_probabilities_consistency_various_rates(self, error_rate: float, n_qudits: int):
        model = DepolarizingNoise(error_rate=error_rate)
        probs = model.kraus_probabilities()
        assert all([0 <= p <= 1 for p in probs])
        assert np.isclose(sum(probs), 1.0)


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
        probs = model.kraus_probabilities()

        assert np.isclose(probs[0], 1.0)
        assert np.isclose(probs[1], 0.0)

    def test_depolarizing_zero_error_rate(self):
        """Zero error rate should give probability 1 for identity Kraus operator."""
        model = DepolarizingNoise(error_rate=0.0)
        probs = model.kraus_probabilities()

        assert np.isclose(probs[0], 1.0)
        assert np.allclose(probs[1:], 0.0)


class TestCompositeNoise:
    """Tests specific to CompositeNoise class."""

    def test_from_noise_models_empty_list_raises(self):
        """Empty noise_models list should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            CompositeNoise.from_noise_models([])

    def test_from_noise_models_single_model(self):
        """CompositeNoise with single model should work."""
        model = CompositeNoise.from_noise_models([DephasingNoise(0.1)])
        assert len(model.kraus_gates()) == 2  # I and Z

    def test_n_kraus_operators_combines_unique(self):
        """n_kraus_operators should return count of unique Paulis."""
        # Dephasing has I, Z (2 operators)
        # Depolarizing has I, X, Y, Z (4 operators)
        # Combined should have I, X, Y, Z (4 unique operators)
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.1),
            DepolarizingNoise(0.1)
        ])
        assert len(composite.kraus_gates()) == 4

    def test_probabilities_consistency(self):
        """All probabilities should be non-negative."""
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.5),
            DepolarizingNoise(0.5)
        ])
        probs = composite.kraus_probabilities()
        assert all([0 <= p <= 1 for p in probs])
        assert np.isclose(sum(probs), 1.0)

    def test_weights_are_averaged(self):
        """Weights for same Pauli should be averaged across models."""
        # Two identical models should give same result as one
        single = DephasingNoise(0.1)
        composite = CompositeNoise.from_noise_models([
            DephasingNoise(0.1),
            DephasingNoise(0.1)
        ])

        single_ops = single.kraus_probabilities()
        composite_ops = composite.kraus_probabilities()

        # Same operators, same weights (since averaging identical models)
        assert len(single_ops) == len(composite_ops)
        for s_op, c_op in zip(single_ops, composite_ops):
            assert np.isclose(s_op, c_op)
