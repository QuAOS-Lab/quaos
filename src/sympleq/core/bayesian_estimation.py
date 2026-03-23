from typing import Any, Callable, Generator, Hashable
import numpy as np

from sympleq.core.paulis._typing import HilbertOperator


def fidelity(rho: np.ndarray | HilbertOperator, sigma: np.ndarray | HilbertOperator) -> float:
    def _sparse_fidelity(rho: HilbertOperator, sigma: HilbertOperator) -> float:
        from scipy.sparse.linalg import eigsh
        lam, V = eigsh(rho)
        # sqrt(rho) = V @ diag(sqrt(lam)) @ V^dag
        sqrt_lam = np.sqrt(np.maximum(lam, 0))
        # M = sqrt(rho) @ sigma @ sqrt(rho)
        # = V @ diag(sqrt_lam) @ V^dag @ sigma @ V @ diag(sqrt_lam) @ V^dag
        inner = (V * sqrt_lam) .T.conj() @ sigma @ (V * sqrt_lam)
        # Eigendecompose M to get sqrt(M)
        mu = np.linalg.eigvalsh(inner)
        return float(np.real(np.sum(np.sqrt(np.maximum(mu, 0)))) ** 2)

    if isinstance(rho, HilbertOperator) and isinstance(sigma, HilbertOperator):
        return _sparse_fidelity(rho, sigma)

    from scipy.linalg import sqrtm
    sq_rho = sqrtm(rho)
    tmp_matrix = sq_rho @ sigma @ sq_rho
    return float(np.real(np.trace(sqrtm(tmp_matrix)) ** 2))


class BayesianEstimator:
    def __init__(self, threshold: float = 10**(-2), min_runs: int = 100) -> None:
        self._results: dict[Hashable, int] = {}
        self._probabilities: dict[Hashable, float] = {}
        self._probabilities_squared: dict[Hashable, float] = {}
        self._variances: dict[Hashable, float] = {}

        self._base = 1
        self._counts_tot = 0
        self._a_tot = 0

        self.threshold = threshold
        self.min_runs = min_runs

    def __str__(self) -> str:
        return (f"BayesianEstimator(threshold={self.threshold}, "
                f"min_runs={self.min_runs})")

    def results(self) -> list[Any]:
        return list(self._results.keys())

    def num_runs(self) -> int:
        return self._counts_tot

    def probability(self, key: Hashable) -> float:
        return self._probabilities[key]

    def variance(self, key: Hashable) -> float:
        return self._variances[key]

    def _update_estimates(self, base: float, counts_tot: int, a_tot: float):
        denom = counts_tot + a_tot
        for key, count in self._results.items():
            p = (count + base) / denom
            p2 = (count + base) * (count + base + 1) / (denom * (denom + 1))
            self._probabilities[key] = p
            self._probabilities_squared[key] = p2
            self._variances[key] = (p2 - p**2)

    def _record_result(self, result: Hashable):
        is_new = result not in self._results
        self._results[result] = self._results.get(result, 0) + 1
        self._counts_tot += 1
        if is_new:
            self._a_tot += self._base
        self._update_estimates(self._base, self._counts_tot, self._a_tot)

    def _converged(self) -> bool:
        return (self._counts_tot >= self.min_runs and all(v <= self.threshold for v in self._variances.values()))

    def run(self, callable: Callable[[], Hashable]):
        for _ in self.run_iter(callable):
            pass

    def run_iter(self, callable: Callable[[], Hashable]) -> Generator[None, None, None]:
        while not self._converged():
            result = callable()
            self._record_result(result)
            yield

    def report(self) -> str:
        total = sum(self._results.values())
        lines = [f"BayesianEstimator report ({total} samples, {len(self._results)} outcomes)"]
        for key in sorted(self._results, key=lambda k: -self._probabilities.get(k, 0)):
            p = self._probabilities[key]
            var = self._variances[key]
            count = self._results[key]
            lines.append(f"  {key}: p={p:.4f} ± {var:.4f} (n={count})")
        return "\n".join(lines)
