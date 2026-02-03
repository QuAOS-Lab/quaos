import numpy as np
import json
from unittest.mock import patch, MagicMock
import warnings
import pytest

from sympleq.applications.measurement.aquire import Aquire, AquireConfig, simulate_measurement
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.applications.measurement.covariance_graph import commutation_graph, all_maximal_cliques
from sympleq.applications.measurement.allocation import construct_circuit_list
from sympleq.core.statistic_utils import true_covariance_graph


class TestAquire:

    def compare_loaded_file(self, P, psi, correct_commutation_graph, correct_xxx, correct_circuit_list,
                            correct_pdf_list, simulated_measurement_results, correct_results, true_variance_graph):
        assert P.is_hermitian(), "P not hermitian"
        assert np.abs(np.linalg.norm(psi) - 1) < 10**(-6), "psi not normalized"
        com_graph = commutation_graph(P)
        assert np.all(com_graph.adj == correct_commutation_graph), "commutation graphs do not match"
        clique_covering_unsorted = list(all_maximal_cliques(com_graph))
        clique_covering = [sorted(clique) for clique in clique_covering_unsorted]
        xxx = sorted(clique_covering)
        assert xxx == correct_xxx, "Clique coverings do not match"
        circuit_list, circuit_dictionary = construct_circuit_list(P, xxx, {})
        for i, c in enumerate(circuit_list):
            comp_c = correct_circuit_list[i]
            for j, g in enumerate(c.gates):
                assert g.name == comp_c[j][0] and np.all(g.qudit_indices == comp_c[j][1]), f"Circuit {i} does not match"

        true_cov_graph = true_covariance_graph(P, psi) * com_graph.adj
        assert np.allclose(true_cov_graph, true_variance_graph, atol=10**(-6)), "true covariance graphs do not match"

        for i, aa in enumerate(xxx):
            P1, C, _ = circuit_dictionary[str(aa)]
            psi_diag = C.unitary() @ psi
            pdf = np.abs(psi_diag * psi_diag.conj())
            assert np.allclose(pdf, correct_pdf_list[i]), f"pdf does not match for clique {i}, {aa}"

        d = P.lcm

        for i, aa in enumerate(xxx):
            (P1, C, k_dict) = circuit_dictionary[str(aa)]
            p1, q1, phases1 = P1.n_paulis(), P1.n_qudits(), P1.phases

            bases_a1 = simulated_measurement_results[i]

            ss = [(sum((bases_a1[i1] * P1.z_exp[i0, i1] * P1.lcm) // P1.dimensions[i1]
                       for i1 in range(q1))) % P1.lcm for i0 in range(p1)]

            comp_correct_results = correct_results[i]
            j = 0
            for j0, s0 in enumerate(ss):
                for a0, a1, s1 in k_dict[str(j0)]:
                    assert (phases1[j0] + s1) % 2 != 1, f"Eigenvalue for {i} not even (not able to sort in correctly)"
                    assert (a0 == comp_correct_results[j][0] and a1 == comp_correct_results[j][1])
                    comp = (int(s0 + (phases1[j0] + s1) / 2) % d == comp_correct_results[j][2])
                    assert comp, f"Eigenvalue for clique {i} does not match"
                    j += 1
        pass

    def AEQuO_comparison(self, filename):
        # comparison is with respect to code written by Andrew J. Jena Plinsky
        with open(f'./tests/application_tests/measurement_tests/comparison_json/{filename}' + ".json", "r") as f:
            comparison_data = json.load(f)

        P = PauliSum.from_string(comparison_data["strings"],
                                 weights=np.array(comparison_data["weights_real"]) +
                                 1j * np.array(comparison_data["weights_imag"]),
                                 dimensions=comparison_data["dimensions"],
                                 phases=comparison_data["phases"])
        P = P.XZ_to_Y()
        psi = np.array(comparison_data["psi_real"]) + 1j * np.array(comparison_data["psi_imag"])
        correct_commutation_graph = np.array(comparison_data["commutation_graph"])
        correct_xxx = comparison_data["clique_covering"]
        correct_circuit_list = comparison_data["circuit_list"]
        correct_pdf_list = comparison_data["pdf_list"]
        simulated_measurement_results = comparison_data["simulated_measurement_results"]
        correct_results = comparison_data["comparison_results"]
        true_variance_graph = np.array(comparison_data["true_variance_graph_real"]) + \
            1j * np.array(comparison_data["true_variance_graph_imag"])

        self.compare_loaded_file(P, psi, correct_commutation_graph, correct_xxx, correct_circuit_list,
                                 correct_pdf_list, simulated_measurement_results, correct_results, true_variance_graph)

    def random_comparison_hamiltonian(self, num_paulis, dimensions, mode='rand'):
        paulistrings = []
        weights = []
        phases = []
        selected_paulistring_indexes = []
        while len(paulistrings) < num_paulis:
            ps = PauliString.from_random(dimensions)
            ps_index = ps._to_int()
            if ps_index not in selected_paulistring_indexes:
                # paulistring
                paulistrings.append(ps)

                # weights
                if mode == 'rand':
                    if ps.H().has_equal_tableau(ps):
                        weights.append(np.random.normal(0, 1))
                    else:
                        weights.append(np.random.normal(0, 1) + 1j * np.random.normal(0, 1))
                elif mode == 'uniform' or mode == 'one':
                    weights.append(1)
                else:
                    raise ValueError("mode must be 'rand', 'uniform' or 'one'")

                # phases
                if 2 in dimensions:
                    num_Ys = sum([ps.x_exp[i] * ps.z_exp[i] for i in range(len(ps.x_exp)) if dimensions[i] == 2])
                    phases.append(num_Ys * ps.lcm / 2)
                else:
                    phases.append(0)

                selected_paulistring_indexes.append(ps_index)
                if not ps.H().has_equal_tableau(ps):
                    ps_conj = PauliSum.from_pauli_strings([ps], weights=[weights[-1]], phases=[phases[-1]]).H()
                    ps_conj.set_phases([(ps_conj.phases[0] - phases[-1]) % (2 * ps.lcm)])
                    ps_conj.phase_to_weight()

                    weights.append(ps_conj.weights[0])
                    phases.append(phases[-1])
                    ps_conj_str = ps_conj[0].copy()
                    paulistrings.append(ps_conj_str)
                    selected_paulistring_indexes.append(ps_conj_str._to_int())

        return PauliSum.from_pauli_strings(paulistrings, weights=weights, phases=phases)

    def test_aquire_for_AEQuO_mono_qubit_d222_p10(self):
        self.AEQuO_comparison("mono_qubit_d222_p10")
        pass

    def test_aquire_for_AEQuO_mono_qudit_d333_p20(self):
        self.AEQuO_comparison("mono_qudit_d333_p20")
        pass

    def test_aquire_for_AEQuO_mixed_d2233_p20(self):
        self.AEQuO_comparison("mixed_d2233_p20")
        pass

    def test_aquire_for_AEQuO_mixed_d235_p20(self):
        self.AEQuO_comparison("mixed_d235_p20")
        pass

    def test_aquire_for_AEQuO_mono_qubit_d2222_p40(self):
        self.AEQuO_comparison("mono_qubit_d2222_p40")
        pass

    def test_aquire_for_AEQuO_mono_qudit_d555_p16(self):
        self.AEQuO_comparison("mono_qudit_d555_p16")
        pass

    def test_aquire_run_no_errors(self):
        update_steps = [25, 50]
        dim_list = [[2, 2, 2], [3, 3, 3], [2, 2, 3]]
        for dims in dim_list:
            # general commutation mode
            P = self.random_comparison_hamiltonian(6, dims, mode='rand')
            psi = np.random.rand(int(np.prod(dims))) + 1j * np.random.rand(int(np.prod(dims)))
            psi = psi / np.linalg.norm(psi)

            model = Aquire(H=P, psi=psi)

            model.config.set_params(commutation_mode='general',
                                    calculate_true_values=True,
                                    save_covariance_graph_checkpoints=False,
                                    auto_update_covariance_graph=True,
                                    auto_update_settings=True,
                                    verbose=True)

            model.simulate_observable(update_steps=update_steps)

            # bitwise commutation mode
            P = self.random_comparison_hamiltonian(6, dims, mode='rand')
            psi = np.random.rand(int(np.prod(dims))) + 1j * np.random.rand(int(np.prod(dims)))
            psi = psi / np.linalg.norm(psi)

            model = Aquire(H=P, psi=psi)

            model.config.set_params(commutation_mode='bitwise',
                                    calculate_true_values=False,
                                    save_covariance_graph_checkpoints=False,
                                    auto_update_covariance_graph=True,
                                    auto_update_settings=True,
                                    verbose=True)

            model.simulate_observable(update_steps=update_steps)

            # with noise
            P = self.random_comparison_hamiltonian(6, dims, mode='rand')
            psi = np.random.rand(int(np.prod(dims))) + 1j * np.random.rand(int(np.prod(dims)))
            psi = psi / np.linalg.norm(psi)

            model = Aquire(H=P, psi=psi)

            model.config.set_params(commutation_mode='general',
                                    calculate_true_values=False,
                                    enable_simulated_hardware_noise=True,
                                    enable_diagnostics=True,
                                    save_covariance_graph_checkpoints=False,
                                    auto_update_covariance_graph=True,
                                    auto_update_settings=True,
                                    verbose=True)

            model.simulate_observable(update_steps=update_steps)

    @pytest.mark.system
    def test_aquire_mean_distance(self):
        update_steps = [6, 12, 25, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800]
        dim_list = [[2, 2, 2], [3, 3, 3], [5, 5, 5], [2, 3, 5], [2, 2, 3, 3]]
        for dims in dim_list:
            P = self.random_comparison_hamiltonian(20, dims, mode='rand')
            D = int(np.prod(dims))
            psi = np.random.rand(D) + 1j * np.random.rand(D)
            psi = psi / np.linalg.norm(psi)

            model = Aquire(H=P, psi=psi)

            model.config.set_params(commutation_mode='general',
                                    calculate_true_values=True,
                                    save_covariance_graph_checkpoints=False,
                                    auto_update_covariance_graph=True,
                                    auto_update_settings=True,
                                    verbose=True)

            model.simulate_observable(update_steps=update_steps)
            mean_distance = np.abs(model.true_mean_value - model.estimated_mean[-1])
            distance_in_sigma = mean_distance / np.sqrt(model.statistical_variance[-1])
            assert distance_in_sigma < 5, f"Mean estimate too far from true value for dims {dims}"

    @pytest.mark.system
    def test_aquire_mean_distance_with_noise(self):
        update_steps = [6, 12, 25, 50, 100, 200, 400, 800, 1600]
        dim_list = [[2, 2, 2], [2, 2, 3], [3, 3, 3]]
        for dims in dim_list:
            P = self.random_comparison_hamiltonian(10, dims, mode='rand')
            D = int(np.prod(dims))
            psi = np.random.rand(D) + 1j * np.random.rand(D)
            psi = psi / np.linalg.norm(psi)

            model = Aquire(H=P, psi=psi)

            model.config.set_params(commutation_mode='general',
                                    calculate_true_values=True,
                                    enable_simulated_hardware_noise=True,
                                    enable_diagnostics=True,
                                    save_covariance_graph_checkpoints=False,
                                    auto_update_covariance_graph=True,
                                    auto_update_settings=True,
                                    verbose=True)

            model.simulate_observable(update_steps=update_steps)
            mean_distance = np.abs(model.true_mean_value - model.estimated_mean[-1])
            error = np.sqrt(model.statistical_variance[-1] + model.systematic_variance[-1])
            distance_in_sigma = mean_distance / error
            assert distance_in_sigma - 1 < 5, f"Mean estimate too far from true value for dims {dims}"

    # AQUIRE CONFIG TESTS
    def test_aquire_state_hamiltonian_mismatch_validation(self):
        # check that the function that compares state and hamiltonian dimension raises an error correctly
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        psi = np.random.rand(10) + 1j * np.random.rand(10)
        psi = psi / np.linalg.norm(psi)
        with pytest.raises(ValueError):
            Aquire(H=P, psi=psi)

    def test_aquire_state_normalization_validation(self):
        # check that the function that compares state and hamiltonian dimension raises an error correctly
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        psi = np.random.rand(int(np.prod([3, 3, 3]))) + 1j * np.random.rand(int(np.prod([3, 3, 3])))  # not normalized
        with pytest.raises(ValueError):
            Aquire(H=P, psi=psi)

    def test_aquire_commutation_validation(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        commutation_rule = 'test'
        with pytest.raises(ValueError):
            Aquire(H=P, commutation_mode=commutation_rule)

    def test_aquire_allocation_mode_validation(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        allocation_mode = 'test'
        with pytest.raises(ValueError):
            Aquire(H=P, allocation_mode=allocation_mode)

    def test_aquire_diagnostic_mode_validation(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        diagnostic_mode = 'test'
        with pytest.raises(ValueError):
            Aquire(H=P, diagnostic_mode=diagnostic_mode)
        diagnostic_mode = 'informed'
        with pytest.raises(NotImplementedError):
            Aquire(H=P, diagnostic_mode=diagnostic_mode)

    def test_aquire_true_value_mode_validation(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        calculate_true_values = True
        with pytest.warns(UserWarning):
            model = Aquire(H=P, calculate_true_values=calculate_true_values)
        assert not model.config.calculate_true_values

    def test_aquire_config_update_all(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        model = Aquire(H=P, calculate_true_values=False, verbose=False, auto_update_settings=True)
        model.config.set_params(commutation_mode='local')
        assert model.config.commutation_mode == 'local'
        with pytest.raises(ValueError):
            model.config.set_params(commutation_mode='test')
        model.config.set_params(commutation_mode='local', auto_update_settings=False)
        with pytest.warns(UserWarning):
            model.config.set_params(commutation_mode='general')

    def test_aquire_config_update_set_params(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        model = Aquire(H=P, calculate_true_values=False, verbose=False, auto_update_settings=True)
        with pytest.raises(AttributeError):
            model.config.set_params(test='test')

    def test_aquire_config_verbose(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        model = Aquire(H=P, calculate_true_values=False, verbose=True, auto_update_settings=True)
        with pytest.warns(UserWarning):
            model.config.update_all()

    def test_aquire_config_string_representation(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        psi = np.random.rand(int(np.prod([3, 3, 3]))) + 1j * np.random.rand(int(np.prod([3, 3, 3])))
        psi = psi / np.linalg.norm(psi)
        model = Aquire(H=P, calculate_true_values=False, verbose=False, auto_update_settings=True)
        assert type(str(model.config)) == str
        model = Aquire(H=P, psi=psi, calculate_true_values=True, mcmc_initial_samples_per_chain=500,
                       mcmc_max_samples_per_chain=2000, verbose=False, auto_update_settings=True)
        assert type(str(model.config)) == str

    def test_aquire_config_mcmc_validation(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        model = Aquire(H=P, calculate_true_values=False, verbose=False, auto_update_settings=True)

        def n_test1(shots):
            return shots + 10

        def n_max_test1(shots):
            return int(shots * 0.9) + 100

        with pytest.warns(UserWarning):
            model.config.set_params(mcmc_initial_samples_per_chain=n_test1,
                                    mcmc_max_samples_per_chain=n_max_test1)
            model.config.test_mcmc_settings()

        with pytest.raises(ValueError):
            model.config.set_params(mcmc_initial_samples_per_chain=101.1,
                                    mcmc_max_samples_per_chain=2000)
            model.config.test_mcmc_settings()

        with pytest.raises(ValueError):
            model.config.set_params(mcmc_initial_samples_per_chain=100,
                                    mcmc_max_samples_per_chain=2000.2)
            model.config.test_mcmc_settings()

        def n_test2(shots):
            return shots + 10 + 0.1

        def n_max_test2(shots):
            return shots + 100 + 0.2

        with pytest.warns(UserWarning):
            model.config.set_params(mcmc_initial_samples_per_chain=n_test2,
                                    mcmc_max_samples_per_chain=n_max_test2)
            model.config.test_mcmc_settings()

        with pytest.raises(ValueError):
            model.config.set_params(mcmc_initial_samples_per_chain=-10,
                                    mcmc_max_samples_per_chain=2000)
            model.config.test_mcmc_settings()

        with pytest.raises(ValueError):
            with pytest.warns(UserWarning):
                model.config.set_params(mcmc_initial_samples_per_chain=10,
                                        mcmc_max_samples_per_chain=-2000)
            model.config.test_mcmc_settings()

    def test_aquire_noise_and_error_function_validation(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        model = Aquire(H=P, calculate_true_values=False, verbose=False, auto_update_settings=True)

        def faulty_noise_probability_function(_, *args, **kwargs):
            return -0.1

        def faulty_error_function(results, *args, **kwargs):
            return np.array([1.2 for _ in results])

        with pytest.raises(ValueError):
            model.config.set_params(noise_probability_function=faulty_noise_probability_function)
            model.config.test_noise_and_error_function()
        with pytest.raises(ValueError):
            model.config.set_params(error_function=faulty_error_function)
            model.config.test_noise_and_error_function()

        def faulty_error_function2(results, *args, **kwargs):
            return np.array([10 for _ in results])

        with pytest.raises(ValueError):
            model.config.set_params(error_function=faulty_error_function2)
            model.config.test_noise_and_error_function()

    def test_aquire_config_info(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        model = Aquire(H=P, calculate_true_values=False, verbose=False, auto_update_settings=True)
        model.config.info()
        model.config.info(name="Hamiltonian")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            model.config.info(name='Invalid info')
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)

    # AQUIRE CONFIG TESTS END

    # AQUIRE FUNCTIONALITY TESTS

    def test_aquire_identity_warning(self):
        P = PauliSum.from_string(['x0z0 x0z0', 'x1z1 x1z1'], dimensions=[2, 2], weights=[1, 1], phases=[0, 0])
        with pytest.warns(UserWarning):
            Aquire(H=P, calculate_true_values=False)

    @patch("builtins.input", return_value="y")
    def test_aquire_input_with_wrong_phase_y(self, mock_input):
        P = PauliSum.from_string(['x1z0 x0z0', 'x1z1 x1z0'], dimensions=[2, 2], weights=[1, 1], phases=[0, 0])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Aquire(H=P, calculate_true_values=False)
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)

    @patch("builtins.input", return_value="n")
    def test_aquire_input_with_wrong_phase_n(self, mock_input):
        P = PauliSum.from_string(['x1z0 x0z0', 'x1z1 x1z0'], dimensions=[2, 2], weights=[1, 1], phases=[0, 0])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            Aquire(H=P, calculate_true_values=False)
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)

    @patch("builtins.input", return_value="invalid input")
    def test_aquire_input_with_wrong_phase_invalid(self, mock_input):
        P = PauliSum.from_string(['x1z0 x0z0', 'x1z1 x1z0'], dimensions=[2, 2], weights=[1, 1], phases=[0, 0])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with pytest.raises(Exception):
                Aquire(H=P, calculate_true_values=False)
            assert len(w) == 1
            assert issubclass(w[0].category, UserWarning)

    @patch("builtins.input", side_effect=["y", "y"])
    def test_aquire_non_hermitian_input_yy(self, mock_input):
        P = PauliSum.from_string(['x1z0 x0z0', 'x1z1 x1z0'], dimensions=[3, 3], weights=[1, 1], phases=[0, 0])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            model = Aquire(H=P, calculate_true_values=False)
            assert len(w) == 2
            assert issubclass(w[0].category, UserWarning)
            assert issubclass(w[1].category, UserWarning)
            assert model.H.is_hermitian()

    @patch("builtins.input", side_effect=["y", "n"])
    def test_aquire_non_hermitian_input_yn(self, mock_input):
        P = PauliSum.from_string(['x1z0 x0z0', 'x1z1 x1z0'], dimensions=[3, 3], weights=[1, 1], phases=[0, 0])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with pytest.raises(Exception):
                _ = Aquire(H=P, calculate_true_values=False)
            assert len(w) == 2
            assert issubclass(w[0].category, UserWarning)
            assert issubclass(w[1].category, UserWarning)

    def test_aquire_with_input_config_file(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        AConfig = AquireConfig(P, calculate_true_values=False)
        _ = Aquire(H=P, config=AConfig)

    @patch("builtins.input", side_effect=["y", "y"])   # answer to input(), input2()
    @patch("sympleq.applications.measurement.allocation.sort_hamiltonian")
    def test_make_hermitian_fails(self, mock_sort, mock_input):
        # --- Create a fake P with .is_hermitian() behavior ---
        P_before = MagicMock()               # P returned by sort_hamiltonian
        P_before.is_hermitian.return_value = False

        P_after = MagicMock()                # P returned by make_hermitian
        P_after.is_hermitian.return_value = False   # <-- forces the error branch!

        mock_sort.return_value = (P_before, None)

        H = MagicMock()                      # mock Hamiltonian with trivial attributes
        H.n_paulis.return_value = 0
        H.weights = []
        H.__getitem__ = lambda *args: MagicMock(is_identity=lambda: False)

        # --- Capture the exception ---
        with pytest.raises(Exception):
            Aquire(H)

    def test_aquire_manual_inputs(self):
        P = self.random_comparison_hamiltonian(6, [3, 3, 3], mode='rand')
        psi = np.random.rand(int(np.prod([3, 3, 3]))) + 1j * np.random.rand(int(np.prod([3, 3, 3])))
        psi = psi / np.linalg.norm(psi)

        model = Aquire(H=P, commutation_mode='general',
                       calculate_true_values=False,
                       enable_diagnostics=True, auto_update_settings=False,
                       diagnostic_mode='zero', enable_debug_checks=True, save_covariance_graph_checkpoints=True)
        mock_circuit_list = model.allocate_measurements(shots=100)
        mock_diagnostic_circuit_list, mock_dsp_circuit_list = model.construct_diagnostic_circuits()

        assert model.shots_since_last_update() == 100
        assert model.total_shots() == 100
        assert len(model.cliques_since_last_update()) == 100
        assert len(model.measurement_circuits_since_last_update()) == 100
        assert len(model.diagnostic_circuits_since_last_update()) == 100
        assert len(model.diagnostic_states_since_last_update()) == 100
        assert len(model.diagnostic_state_preparation_circuits_since_last_update()) == 100

        mock_results = []
        for c in mock_circuit_list:
            mock_results.append(simulate_measurement(P, psi, circuit=c))
        mock_diagnostic_results = []
        for i, c in enumerate(mock_diagnostic_circuit_list):
            mock_diagnostic_results.append(simulate_measurement(P, model.diagnostic_states[i], circuit=c))

        model.input_measurement_data(mock_results)
        model.input_diagnostic_data(mock_diagnostic_results)

        model.data_at_shot(100)
        model.data_at_shot([25, 50, 75, 100])

        model.scaling_matrix_at_shot(100)
        model.scaling_matrix_at_shot([25, 50, 75, 100])

        model.covariance_graph_at_shot(10)
        model.covariance_graph_at_shot([10, 20])

        model.diagnostic_data_at_shot(10)
        model.diagnostic_data_at_shot([10, 20])

        str(model)
        model.info()
        model.info(name="circuits")
        fig, ax = model.plot()
