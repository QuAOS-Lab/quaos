from typing import Callable

from .config import RMBConfig, RMBData

type UpdateStrategy = Callable[[RMBData, RMBConfig], RMBConfig]


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
        new_config = current_config.with_depth(current_config.depth + 4)
    elif estimator.probability(True) >= 0.6:
        new_config = current_config.with_depth(current_config.depth + 2)
    elif estimator.probability(True) >= 0.5:
        new_config = current_config.with_depth(current_config.depth + 1)
    else:
        delta = current_config.max_two_qubit_gate_ratio - current_config.min_two_qubit_gate_ratio
        new_max = current_config.min_two_qubit_gate_ratio
        new_min = round(max(0.0, new_max - delta), 2)
        new_config = current_config.with_two_qubit_gate_ratio(new_min, new_max)

    return new_config
