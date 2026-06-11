from dataclasses import dataclass

from sympleq.applications.randomized_benchmarking.RMB import RMB
from sympleq.applications.randomized_benchmarking.backends.sympleq import SympleqBackend
from sympleq.applications.randomized_benchmarking.config import RMBConfig, RMBData
from sympleq.core.noise.noise_model import GenericNoise


import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

from viarregio2 import (
    config_depth,
    config_from_parameters,
    config_two_qubit_gate_ratio,
)


def default_update_strategy(data: RMBData, current_config: RMBConfig) -> RMBConfig:
    """
    Default rule for advancing an RMB sweep (formerly
    ``randomized_benchmarking.update_strategy``; this script is its only user).
    """
    if current_config in data:
        estimator = data[current_config]
    else:
        return RMBConfig.default()
    # probability(True) is the fidelity
    if estimator.probability(True) >= 0.8:
        new_config = current_config.with_n_1qb_gates(current_config.n_1qb_gates + 4 * current_config.n_qubits)
    elif estimator.probability(True) >= 0.6:
        new_config = current_config.with_n_1qb_gates(current_config.n_1qb_gates + 2 * current_config.n_qubits)
    elif estimator.probability(True) >= 0.5:
        new_config = current_config.with_n_1qb_gates(current_config.n_1qb_gates + 1 * current_config.n_qubits)
    else:
        new_config = current_config.with_n_2qb_gates(max(0, current_config.n_2qb_gates - 1 * current_config.n_qubits))

    return new_config


@dataclass(frozen=True)
class FidelityBoundaryFit:
    """
    Monotone logistic fit for the RMB fidelity landscape.

    The model is

        Pr(fidelity=True | depth, ratio) = sigmoid(alpha - gamma_depth * depth
                                                  - gamma_ratio * ratio)

    so the estimated fidelity=0.5 line is alpha - gamma @ theta = 0.
    """
    alpha: float
    gamma: np.ndarray

    def probability(self, theta: np.ndarray) -> np.ndarray:
        theta = np.asarray(theta, dtype=float)
        return expit(self.alpha - theta @ self.gamma)

    def boundary_ratio(self, depth: np.ndarray | float) -> np.ndarray | float:
        """
        Return the two-qubit ratio on the fitted fidelity=0.5 boundary.
        """
        if len(self.gamma) != 2:
            raise ValueError("boundary_ratio is only defined for [depth, ratio] fits.")
        if abs(self.gamma[1]) < 1e-12:
            raise ValueError("Cannot solve boundary for ratio: ratio coefficient is zero.")
        return (self.alpha - self.gamma[0] * np.asarray(depth)) / self.gamma[1]

    def boundary_depth(self, ratio: np.ndarray | float) -> np.ndarray | float:
        """
        Return the depth on the fitted fidelity=0.5 boundary.
        """
        if len(self.gamma) != 2:
            raise ValueError("boundary_depth is only defined for [depth, ratio] fits.")
        if abs(self.gamma[0]) < 1e-12:
            raise ValueError("Cannot solve boundary for depth: depth coefficient is zero.")
        return (self.alpha - self.gamma[1] * np.asarray(ratio)) / self.gamma[0]

    def report(self) -> str:
        if len(self.gamma) == 2:
            lines = [
                "Fitted monotone logistic model:",
                f"  p(True) = sigmoid({self.alpha:.6g} "
                f"- {self.gamma[0]:.6g} * depth "
                f"- {self.gamma[1]:.6g} * two_qubit_ratio)",
                "Estimated fidelity=0.5 boundary:",
                f"  {self.gamma[0]:.6g} * depth "
                f"+ {self.gamma[1]:.6g} * two_qubit_ratio = {self.alpha:.6g}",
            ]
            if abs(self.gamma[1]) >= 1e-12:
                slope = -self.gamma[0] / self.gamma[1]
                intercept = self.alpha / self.gamma[1]
                lines.append(
                    f"  two_qubit_ratio = {intercept:.6g} {slope:+.6g} * depth"
                )
            if abs(self.gamma[0]) >= 1e-12:
                slope = -self.gamma[1] / self.gamma[0]
                intercept = self.alpha / self.gamma[0]
                lines.append(
                    f"  depth = {intercept:.6g} {slope:+.6g} * two_qubit_ratio"
                )
            return "\n".join(lines)

        return (
            "Fitted monotone logistic model:\n"
            f"  alpha={self.alpha:.6g}\n"
            f"  gamma={self.gamma}"
        )


def config_to_theta(config: RMBConfig) -> np.ndarray:
    """
    Convert an RMBConfig into numerical parameters.

    Adjust this if you want more parameters in the model.
    """
    return np.array(
        [
            float(config_depth(config)),
            float(config_two_qubit_gate_ratio(config)),
        ]
    )


def theta_to_config(theta: np.ndarray, template: RMBConfig) -> RMBConfig:
    """
    Convert numerical parameters back into an RMBConfig.

    Here we keep the two-qubit gate ratio as a narrow interval.
    """
    depth = max(1, int(round(theta[0])))

    ratio = float(np.clip(theta[1], 0.0, 1.0))
    return config_from_parameters(template=template, depth=depth, ratio=ratio)


def fit_monotone_logistic(
    X: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray | None = None,
    l2: float = 1e-3,
) -> FidelityBoundaryFit:
    """
    Fit

        Pr(True | theta) = sigmoid(alpha - gamma @ theta)

    with gamma_i >= 0.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)

    _, d = X.shape
    if weights is None:
        weights = np.ones_like(y)
    else:
        weights = np.asarray(weights, dtype=float)
        weights = weights / np.mean(weights)

    def loss_and_grad(w):
        alpha = w[0]
        gamma = w[1:]

        z = alpha - X @ gamma
        p = expit(z)

        eps = 1e-12
        nll = -np.sum(weights * (y * np.log(p + eps) + (1 - y) * np.log(1 - p + eps)))
        nll += 0.5 * l2 * np.sum(gamma**2)

        r = weights * (p - y)
        grad_alpha = np.sum(r)
        grad_gamma = -X.T @ r + l2 * gamma

        return nll, np.concatenate([[grad_alpha], grad_gamma])

    w0 = np.zeros(d + 1)

    bounds = [(None, None)] + [(0.0, None)] * d

    result = minimize(
        fun=lambda w: loss_and_grad(w)[0],
        x0=w0,
        jac=lambda w: loss_and_grad(w)[1],
        bounds=bounds,
        method="L-BFGS-B",
    )

    if not result.success:
        raise RuntimeError(result.message)

    return FidelityBoundaryFit(alpha=float(result.x[0]), gamma=result.x[1:])


def fit_fidelity_boundary(
    data: RMBData,
    *,
    min_points: int = 4,
    skip_incomplete: bool = False,
    l2: float = 1e-3,
) -> FidelityBoundaryFit:
    """
    Fit and return the estimated fidelity=0.5 boundary for a set of RMB data.

    Each estimator contributes its posterior mean `probability(True)` as the
    fidelity estimate and is weighted by the number of recorded samples.
    """
    X = []
    y = []
    weights = []

    for config, estimator in data.items():
        if skip_incomplete and not estimator.is_converged():
            continue
        if estimator.num_runs() == 0:
            continue
        X.append(config_to_theta(config))
        y.append(estimator.probability(True))
        weights.append(max(1, estimator.num_runs()))

    if len(X) < min_points:
        raise ValueError(
            f"Need at least {min_points} data points to fit the boundary, got {len(X)}."
        )

    return fit_monotone_logistic(
        np.asarray(X, dtype=float),
        np.asarray(y, dtype=float),
        weights=np.asarray(weights, dtype=float),
        l2=l2,
    )


def plot_fidelity_boundary(
    data: RMBData,
    *,
    skip_incomplete: bool = False,
    show: bool = True,
) -> list:
    """
    Plot the RMB data and overlay a fitted fidelity=0.5 boundary per n_qubits.
    """
    import matplotlib.pyplot as plt
    from sympleq.applications.randomized_benchmarking.experiments.plots import plot_data

    axes = plot_data(data, show=False, skip_incomplete=skip_incomplete)
    if not axes:
        return []

    groups: dict[int, RMBData] = {}
    for config, estimator in data.items():
        if skip_incomplete and not estimator.is_converged():
            continue
        if estimator.num_runs() == 0:
            continue
        groups.setdefault(config.n_qubits, {})[config] = estimator

    for ax, (n_qubits, group) in zip(axes, sorted(groups.items())):
        if len(group) < 4:
            continue
        try:
            fit = fit_fidelity_boundary(group, skip_incomplete=skip_incomplete)
        except (RuntimeError, ValueError):
            continue
        x_min, x_max = ax.get_xlim()
        depths = np.linspace(max(1.0, x_min), x_max, 200)
        ratios = fit.boundary_ratio(depths)
        in_view = (0.0 <= ratios) & (ratios <= 1.0)
        if np.any(in_view):
            ax.plot(
                depths[in_view],
                ratios[in_view],
                color="black",
                linewidth=2,
                label="fidelity = 0.5 fit",
                zorder=4,
            )
            ax.legend(loc="best")
        ax.set_title(f"# Qubits = {n_qubits}")

    if show:
        plt.show()

    return axes


def logistic_boundary_update_strategy(
    data: RMBData,
    current_config: RMBConfig,
) -> RMBConfig:
    """
    Active-learning update strategy.

    Uses all previous RMB data to fit

        Pr(True | theta) = sigmoid(alpha - gamma @ theta),

    then proposes a new config close to the estimated p=True = 0.5
    boundary.
    """

    # Need at least a few points before fitting a useful model.
    if len(data) < 4:
        return default_update_strategy(data, current_config)

    try:
        fit = fit_fidelity_boundary(data)
    except (RuntimeError, ValueError):
        return default_update_strategy(data, current_config)
    alpha = fit.alpha
    gamma = fit.gamma

    # Work around current point, but move it onto the fitted p=0.5 boundary:
    #
    #     alpha - gamma @ theta = 0.
    #
    theta = config_to_theta(current_config)

    if np.linalg.norm(gamma) < 1e-12:
        return default_update_strategy(data, current_config)

    signed_distance = (alpha - gamma @ theta) / (np.dot(gamma, gamma))
    theta_boundary = theta + signed_distance * gamma

    # Add a small exploration move along the boundary.
    # For 2 parameters, a vector perpendicular to gamma is easy.
    if len(gamma) == 2:
        tangent = np.array([gamma[1], -gamma[0]])
        tangent_norm = np.linalg.norm(tangent)

        if tangent_norm > 1e-12:
            tangent = tangent / tangent_norm

            # Alternate directions to avoid always walking one way.
            direction = 1.0 if len(data) % 2 == 0 else -1.0
            theta_boundary = theta_boundary + direction * tangent

    new_config = theta_to_config(theta_boundary, current_config)

    # Avoid returning an already-run config if possible.
    if new_config in data:
        return default_update_strategy(data, current_config)

    return new_config


if __name__ == "__main__":
    backend = SympleqBackend()

    noise_model = GenericNoise.from_paulis([0.00075, 0.00075, 0.00075])   # 000025
    two_qubit_noise_model = GenericNoise.from_paulis([0.005, 0.005, 0.005])  # 00079

    backend = SympleqBackend(noise_model=noise_model, two_qubit_noise_model=two_qubit_noise_model)

    rmb = RMB.default().with_backend(backend).with_update_strategy(logistic_boundary_update_strategy)  # takes a config and maps to bayesian eztimator

    template = (
        RMBConfig.default()
        .with_n_qubits(3)
        .with_scrambling_probability(0.85)
    )
    config = config_from_parameters(template=template, depth=10, ratio=0.3)

    rmb.run(config) 

    # rmb.save('viarregio1')

    fit = fit_fidelity_boundary(rmb._data)
    print(fit.report())

    plot_fidelity_boundary(rmb._data)
