import numpy as np


def example_results_calibration(ss, circuit_list_total, algorithm_variables, mode='Null', p_noise=0.01):
    P, cc, pauli_block_sizes, X, S, D, CG, aaa, k_phases, general_commutation, allocation_mode, N_chain, N_mcmc, N_mcmc_max, mcmc_shot_scale, shots_total, circuit_list_total, xxx_total = algorithm_variables
    p = P.paulis()
    dims = P.dims
    q = P.qudits()
    rr = []
    if mode == 'Null':
        for i in range(len(ss)):
            if np.random.uniform() <= p_noise:
                rr.append(np.array([np.random.randint(dims[j])
                          for j in range(q)]))
            else:
                rr.append(np.zeros(q))

    return (rr)
