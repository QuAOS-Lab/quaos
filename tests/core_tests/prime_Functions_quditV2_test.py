import pytest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)

from quaos.core.prime_Functions_quditV2 import (
    read_luca_test_2, random_pauli_hamiltonian, pauli_hermitian, 
    sort_hamiltonian, xi, rand_state, truncated_exponential_sample, 
    get_p_matrix, get_psi, get_p, mcmc_starting_point, psi_sample, 
    log_posterior_ratio, mcmc_covariance_estimate, geweke_test, 
    gelman_rubin_test, update_chain, mcmc_integration, get_alpha, 
    bayes_covariance_estimation, bayes_Var_estimate, variance_graph, 
    bayes_covariance_graph, noise_adder_sim, clique_sampling, levi_civita, 
    quditwise_inner_product, get_phase_matrix, perform_measurements, 
    bucket_filling_qudit, construct_circuit_list, 
    construct_diagonalization_circuit, update_X, choose_measurement, 
    bfq_experiment_initial, bfq_experiment, bfq_estimation, example_results, 
    noise_adder, diagnosis_states, example_results_calibration, 
    error_calibration, bfq_error_correction, error_correction_estimation, 
    example_calibration,
)


class TestPrimeFunctionQuditV2:

    def test_read_luca_test_2(self) -> None:
        # TODO: test_read_luca_test_2
        assert True

    def test_random_pauli_hamiltonian(self) -> None:
        # TODO: test_random_pauli_hamiltonian
        assert True

    def test_pauli_hermitian(self) -> None:
        # TODO: test_pauli_hermitian
        assert True

    def test_sort_hamiltonian(self) -> None:
        # TODO: test_sort_hamiltonian
        assert True

    def test_xi(self) -> None:
        # TODO: test_f
        assert True

    def test_rand_state(self) -> None:
        # TODO: test_rand_state
        assert True

    def test_truncated_exponential_sample(self) -> None:
        # TODO: test_truncated_exponential_sample
        assert True

    def test_get_p_matrix(self) -> None:
        # TODO: test_get_p_matrix
        assert True

    def test_get_psi(self) -> None:
        # TODO: test_get_psi
        assert True

    def test_get_p(self) -> None:
        # TODO: test_get_p
        assert True

    def test_mcmc_starting_point(self) -> None:
        # TODO: test_mcmc_starting_point
        assert True

    def test_psi_sample(self) -> None:
        # TODO: test_psi_sample
        assert True

    def test_log_posterior_ratio(self) -> None:
        # TODO: test_log_posterior_ratio
        assert True

    def test_mcmc_covariance_estimate(self) -> None:
        # TODO: test_mcmc_covariance_estimate
        assert True

    def test_geweke_test(self) -> None:
        # TODO: test_geweke_test
        assert True

    def test_gelman_rubin_test(self) -> None:
        # TODO: test_gelman_rubin_test
        assert True

    def test_update_chain(self) -> None:
        # TODO: test_update_chain
        assert True

    def test_mcmc_integration(self) -> None:
        # TODO: test_mcmc_integration
        assert True

    def test_get_alpha(self) -> None:
        # TODO: test_get_alpha
        assert True

    def test_bayes_covariance_estimation(self) -> None:
        # TODO: test_bayes_covariance_estimation
        assert True

    def test_bayes_Var_estimate(self) -> None:
        # TODO: test_bayes_Var_estimate
        assert True

    def test_variance_graph(self) -> None:
        # TODO: test_variance_graph
        assert True

    def test_bayes_covariance_graph(self) -> None:
        # TODO: test_bayes_covariance_graph
        assert True

    def test_noise_adder_sim(self) -> None:
        # TODO: test_noise_adder_sim
        assert True

    def test_clique_sampling(self) -> None:
        # TODO: test_clique_sampling
        assert True

    def test_levi_civita(self) -> None:
        # TODO: test_levi_civita
        assert True

    def test_quditwise_inner_product(self) -> None:
        # TODO: test_quditwise_inner_product
        assert True

    def test_get_phase_matrix(self) -> None:
        # TODO: test_get_phase_matrix
        assert True

    def test_perform_measurements(self) -> None:
        # TODO: test_perform_measurements
        assert True

    def test_bucket_filling_qudit(self) -> None:
        # TODO: test_bucket_filling_qudit
        assert True

    def test_construct_circuit_list(self) -> None:
        # TODO: test_construct_circuit_list
        assert True

    def test_construct_diagonalization_circuit(self) -> None:
        # TODO: test_construct_diagonalization_circuit
        assert True

    def test_update_X(self) -> None:
        # TODO: test_update_X
        assert True

    def test_choose_measurement(self) -> None:
        # TODO: test_choose_measurement
        assert True

    def test_bfq_experiment_initial(self) -> None:
        # TODO: test_bfq_experiment_initial
        assert True

    def test_bfq_experiment(self) -> None:
        # TODO: test_bfq_experiment
        assert True

    def test_bfq_estimation(self) -> None:
        # TODO: test_bfq_estimation
        assert True

    def test_example_results(self) -> None:
        # TODO: test_example_results
        assert True

    def test_noise_adder(self) -> None:
        # TODO: test_noise_adder
        assert True

    def test_diagnosis_states(self) -> None:
        # TODO: test_diagnosis_states
        assert True

    def test_example_results_calibration(self) -> None:
        # TODO: test_example_results_calibration
        assert True

    def test_error_calibration(self) -> None:
        # TODO: test_error_calibration
        assert True

    def test_bfq_error_correction(self) -> None:
        # TODO: test_bfq_error_correction
        assert True

    def test_error_correction_estimation(self) -> None:
        # TODO: test_error_correction_estimation
        assert True

    def test_example_calibration(self) -> None:
        # TODO: test_example_calibration
        assert True
