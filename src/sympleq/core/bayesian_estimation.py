import math
from typing import Callable, Hashable


class BayesianEstimator:
    def __init__(self, threshold: float = 10**(-2), min_runs: int = 100) -> None:
        self._results: dict[Hashable, int] = {}
        self._probabilities: dict[Hashable, float] = {}
        self._probabilities_squared: dict[Hashable, float] = {}
        self._variances: dict[Hashable, float] = {}

        self._counts_tot = 0
        self._a = []
        self._a_tot = 0

        self.threshold = threshold
        self.min_runs = min_runs

    def probability(self, key: Hashable) -> float:
        return self._probabilities[key]

    def variance(self, key: Hashable) -> float:
        return self._variances[key]

    def run(self, callable: Callable[[], Hashable]):
        base = 1

        for _ in range(self.min_runs):
            result = callable()
            self._results[result] = self._results.get(result, 0) + 1

        counts_tot = self.min_runs
        a_tot = base * len(self._results)

        self._update_estimates(base, counts_tot, a_tot)

        while not all(v <= self.threshold for v in self._variances.values()):
            result = callable()
            is_new = result not in self._results
            self._results[result] = self._results.get(result, 0) + 1
            counts_tot += 1
            if is_new:
                a_tot += base

            self._update_estimates(base, counts_tot, a_tot)

    def __str__(self) -> str:
        return (f"BayesianEstimator(threshold={self.threshold}, "
                f"min_runs={self.min_runs})")

    def report(self) -> str:
        total = sum(self._results.values())
        lines = [f"BayesianEstimator report ({total} samples, {len(self._results)} outcomes)"]
        for key in sorted(self._results, key=lambda k: -self._probabilities.get(k, 0)):
            p = self._probabilities[key]
            var = self._variances[key]
            count = self._results[key]
            lines.append(f"  {key}: p={p:.4f} ± {var:.4f} (n={count})")
        return "\n".join(lines)

    def _update_estimates(self, base: float, counts_tot: int, a_tot: float):
        denom = counts_tot + a_tot
        for key, count in self._results.items():
            p = (count + base) / denom
            p2 = (count + base) * (count + base + 1) / (denom * (denom + 1))
            self._probabilities[key] = p
            self._probabilities_squared[key] = p2
            self._variances[key] = math.sqrt(p2 - p**2)
