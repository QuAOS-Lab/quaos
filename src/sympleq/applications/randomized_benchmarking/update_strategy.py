from typing import Protocol

from .config import RMBConfig, RMBData


class UpdateStrategy(Protocol):
    def __call__(self, data: RMBData, current_config: RMBConfig) -> RMBConfig:
        ...


def default_update_strategy(data: RMBData, current_config: RMBConfig) -> RMBConfig:
    """
    Default rule for advancing an RMB sweep.

    Parameters
    ----------
    data : RMBData
        Mapping from previously seen configurations to their Bayesian estimators.
    current_config : RMBConfig
        Configuration that was just run (or about to be run).

    Returns
    -------
    RMBConfig
        New configuration to run next.
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
