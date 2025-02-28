import pytest
import numpy as np

import sys
from pathlib import Path
root = Path(__file__).parent.parent
sys.path.append(str(root))

from quaos.prime_Functions_quditV2 import *


@pytest.mark.benchmark
def test_main_function(benchmark):
    P,cc = random_pauli_hamiltonian(20,[3,3,3])
    p,q = P.paulis(),P.qudits()

    # sort Hamiltonian for calculation
    P,cc,pauli_block_sizes = sort_hamiltonian(P,cc)

    # state
    psi = ground_state(P,cc)

    # simulation settings
    shots = 12801
    D = {}
    part_func = weighted_vertex_covering_maximal_cliques
    full_simulation = False
    general_commutation = True
    intermediate_results_list = [6,12,25,50,100,200,400,800,1600,3200,6400,12800]
    update_steps = np.array([6,12,25,50,100,200,400,800,1600,3200,6400,12800,25600,51200,102400])
    allocation_mode = 'set'
    mcmc_shot_scale = 0
    N = 500
    N_max = 4001
    p_noise = 0

    # simulation
    benchmark(bucket_filling_qudit, 
              P, cc, psi, shots, part_func, pauli_block_sizes,
              full_simulation=full_simulation,
              update_steps=update_steps,
              general_commutation=general_commutation,
              D=D,
              M_list=intermediate_results_list,
              allocation_mode=allocation_mode,
              mcmc_shot_scale=mcmc_shot_scale,
              N_mcmc=N, 
              N_mcmc_max=N_max,
              p_noise=p_noise,
              Q_progress_bar=False)
