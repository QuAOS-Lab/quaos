from typing import Any, Callable, Generator, Hashable


class BayesianEstimator:
    def __init__(self, threshold: float = 10**(-2), min_runs: int = 100, max_runs: int | None = None) -> None:
        """
        Initialize the Bayesian estimator.

        Parameters
        ----------
        threshold : float, optional
            Variance threshold for convergence. The estimator runs until
            all outcome variances are below this value.
        min_runs : int, optional
            Minimum number of samples before checking convergence.
        max_runs : int or None, optional
            Maximum number of samples. If ``None``, no upper limit is imposed.
        """
        self._results: dict[Hashable, int] = {}
        self._probabilities: dict[Hashable, float] = {}
        self._probabilities_squared: dict[Hashable, float] = {}
        self._variances: dict[Hashable, float] = {}

        self._base = 1
        self._counts_tot = 0
        self._a_tot = 0

        self.threshold = threshold
        self.min_runs = min_runs
        self.max_runs = max_runs

    def __str__(self) -> str:
        return (f"BayesianEstimator(threshold={self.threshold}, "
                f"min_runs={self.min_runs}, max_runs={self.max_runs})")

    def results(self) -> list[Any]:
        """
        Return the list of distinct observed outcomes.

        Returns
        -------
        list[Any]
            Observed outcomes in insertion order.
        """
        return list(self._results.keys())

    def num_runs(self) -> int:
        """
        Return the total number of samples collected so far.

        Returns
        -------
        int
            Total sample count.
        """
        return self._counts_tot

    def probability(self, key: Hashable) -> float:
        """
        Return the estimated probability of an outcome.

        Parameters
        ----------
        key : Hashable
            The outcome to query.

        Returns
        -------
        float
            Bayesian posterior mean probability.
        """
        if key not in self._probabilities:
            return 0.0
        return self._probabilities[key]

    def variance(self, key: Hashable) -> float:
        """
        Return the estimated variance of an outcome's probability.

        Parameters
        ----------
        key : Hashable
            The outcome to query.

        Returns
        -------
        float
            Bayesian posterior variance.
        """
        if key not in self._probabilities:
            return 0.0
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

    def run(self, callable: Callable[[], Hashable], verbose: bool = False):
        """
        Run the estimator to convergence.

        Repeatedly calls ``callable`` and records results until the
        variance threshold is met or ``max_runs`` is reached.

        Parameters
        ----------
        callable : Callable[[], Hashable]
            A zero-argument function that returns a hashable outcome.
        """
        for _ in self.run_iter(callable, verbose):
            pass

    def run_iter(self, callable: Callable[[], Hashable], verbose: bool = False) -> Generator[None, None, None]:
        """
        Run the estimator, yielding after each sample.

        Same convergence logic as :meth:`run`, but yields control after
        each sample so the caller can inspect intermediate results.

        Parameters
        ----------
        callable : Callable[[], Hashable]
            A zero-argument function that returns a hashable outcome.

        Yields
        ------
        None
            Yields after each sample is recorded.
        """
        if verbose:
            import time
            now = time.time()
            n_printed = 0

        while not self._converged():
            if self.max_runs and self.num_runs() > self.max_runs:
                break
            result = callable()
            self._record_result(result)
            if verbose:
                import numpy as np
                if n_printed > 0:
                    print(f"\033[{n_printed}A", end="")

                n_printed = 1
                print(f"Threshold={self.threshold} - {time.time() - now:.2f}s")

                results: list = self.results()
                for res in results:
                    p = self.probability(res)
                    std = np.sqrt(self.variance(res))
                    print(f"\033[K{res}: p={p:.5f} ± {std:.5f}")
                n_printed += len(results)

                n_runs = self.num_runs()
                print(f"n_runs={n_runs}\n")
                n_printed += 2
            yield

    def report(self) -> str:
        """
        Return a human-readable summary of the estimation results.

        Returns
        -------
        str
            Multi-line string listing each outcome with its estimated
            probability, variance, and sample count.
        """
        total = sum(self._results.values())
        lines = [f"BayesianEstimator report ({total} samples, {len(self._results)} outcomes)"]
        for key in sorted(self._results, key=lambda k: -self._probabilities.get(k, 0)):
            p = self._probabilities[key]
            var = self._variances[key]
            count = self._results[key]
            lines.append(f"  {key}: p={p:.4f} ± {var:.4f} (n={count})")
        return "\n".join(lines)
