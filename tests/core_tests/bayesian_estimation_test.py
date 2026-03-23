import numpy as np
import pytest

from sympleq.core.bayesian_estimation import BayesianEstimator


class TestBayesianEstimator:

    def test_fair_coin(self):
        """Estimator recovers fair coin probabilities."""
        rng = np.random.default_rng(42)
        estimator = BayesianEstimator(threshold=1e-3, min_runs=500)
        estimator.run(lambda: rng.choice(["H", "T"]))

        assert abs(estimator.probability("H") - 0.5) < 0.05
        assert abs(estimator.probability("T") - 0.5) < 0.05

    def test_biased_coin(self):
        """Estimator recovers biased coin probabilities."""
        rng = np.random.default_rng(123)
        estimator = BayesianEstimator(threshold=1e-3, min_runs=1000)
        estimator.run(lambda: "H" if rng.random() < 0.8 else "T")

        assert abs(estimator.probability("H") - 0.8) < 0.05
        assert abs(estimator.probability("T") - 0.2) < 0.05

    def test_probabilities_sum_to_one(self):
        """All estimated probabilities sum to approximately 1."""
        rng = np.random.default_rng(7)
        estimator = BayesianEstimator(threshold=1e-3, min_runs=500)
        estimator.run(lambda: rng.choice([1, 2, 3]))

        total = sum(estimator.probability(k) for k in estimator.results())
        assert abs(total - 1.0) < 1e-6

    def test_variances_below_threshold(self):
        """All variances are at or below the threshold after convergence."""
        rng = np.random.default_rng(0)
        threshold = 1e-3
        estimator = BayesianEstimator(threshold=threshold, min_runs=100)
        estimator.run(lambda: rng.choice(["a", "b"]))

        for key in estimator.results():
            assert estimator.variance(key) <= threshold

    def test_min_runs_respected(self):
        """Estimator runs at least min_runs times."""
        rng = np.random.default_rng(0)
        min_runs = 200
        estimator = BayesianEstimator(threshold=1.0, min_runs=min_runs)
        estimator.run(lambda: rng.choice([0, 1]))

        assert estimator.num_runs() >= min_runs

    def test_max_runs_respected(self):
        """Estimator stops at max_runs even if not converged."""
        rng = np.random.default_rng(0)
        estimator = BayesianEstimator(threshold=1e-10, min_runs=10, max_runs=50)
        estimator.run(lambda: rng.choice([0, 1]))

        assert estimator.num_runs() <= 51  # may overshoot by 1 due to check order

    def test_single_outcome(self):
        """Estimator handles a deterministic callable."""
        estimator = BayesianEstimator(threshold=1e-3, min_runs=100)
        estimator.run(lambda: "only")

        assert estimator.probability("only") == pytest.approx(1.0, abs=0.02)
        assert estimator.results() == ["only"]
