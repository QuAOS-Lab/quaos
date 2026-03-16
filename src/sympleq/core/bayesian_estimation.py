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
        for _ in range(self.min_runs):
            result = callable()
            if result not in self._results:
                self._results[result] = 1
            else:
                self._results[result] += 1

        base = 1

        tot = len(self._results)
        counts_tot = self.min_runs
        a_tot = base * tot

        for key, value in self._results.items():
            self._probabilities[key] = (value + base) / (counts_tot + a_tot)
            self._probabilities_squared[key] = (value + base) / (counts_tot + a_tot) * \
                (value + base + 1) / (counts_tot + a_tot + 1)
            self._variances[key] = math.sqrt(self._probabilities_squared[key] - self._probabilities[key]**2)

        while not all(self._variances.values()) <= self.threshold:
            result = callable()
            if result not in self._results:
                self._results[result] = 1
            else:
                self._results[result] += 1

            self._probabilities[result] = (value + base) * (counts_tot + a_tot) / (counts_tot + 1 + a_tot + base)
            self._probabilities_squared[result] = \
                (value + base) * (counts_tot + a_tot) / (counts_tot + 1 + a_tot + base) * \
                (value + base + 1) * (counts_tot + a_tot + 1) / (counts_tot + 1 + a_tot + 1 + base)
            self._variances[key] = math.sqrt(self._probabilities_squared[key] - self._probabilities[key]**2)

            tot += 1
            counts_tot += 1
            a_tot += base
