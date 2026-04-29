from typing import Callable

from .config import RMBConfig, RMBData

type UpdateStrategy = Callable[[RMBData, RMBConfig], RMBConfig]


def default_update_strategy(data: RMBData, current_config: RMBConfig) -> RMBConfig:
    """
    Default rule for advancing an RMB sweep.

    Parameters
    ----------
    data : RMBData
        Mapping from previously seen configurations to their
        Bayesian estimators.
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
        new_config = current_config.with_depth(current_config.depth + 4)
    elif estimator.probability(True) >= 0.6:
        new_config = current_config.with_depth(current_config.depth + 2)
    elif estimator.probability(True) >= 0.5:
        new_config = current_config.with_depth(current_config.depth + 1)
    else:
        new_value = round(current_config.two_qubit_gate_ratio * 0.9, 2)
        new_config = current_config.with_two_qubit_gate_ratio(new_value)

    return new_config
