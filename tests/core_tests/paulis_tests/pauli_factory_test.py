import numpy as np
import pytest
import math
import random
from pathlib import Path
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from tests import PRIME_LIST


class TestPauliSumFactories:

    def test_pauli_sum_from_pauli_strings_single(self):
        ps = PauliString.from_exponents([0, 1], [1, 0], [3, 2])
        assert ps.shape() == (1, 2)

        P = PauliSum.from_pauli_strings(ps, weights=2.0, phases=[1])
        assert P.tableau.shape == (1, 4)
        assert P.weights[0] == 2.0
        assert P.phases[0] == 1

    def test_pauli_sum_from_pauli_strings_multiple(self):
        ps1 = PauliString.from_exponents([0, 1], [1, 2], [5, 3])
        ps2 = PauliString.from_exponents([1, 1], [0, 1], [5, 3])
        with pytest.warns(UserWarning, match="Phases are disregarded"):
            P = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0, 2.0], inherit_phases=True, phases=[0, 1])

        assert P.tableau.shape == (2, ps1.tableau.shape[1])
        assert np.allclose(P.weights, [1.0, 2.0])

    def test_pauli_sum_from_pauli_strings_dimension_mismatch(self):
        ps1 = PauliString.from_exponents([0], [1], [3])
        ps2 = PauliString.from_exponents([0, 1], [1, 0], [2, 2])
        with pytest.raises(ValueError):
            _ = PauliSum.from_pauli_strings([ps1, ps2])

    def test_pauli_sum_from_pauli_strings_invalid_input(self):
        dimensions = [3, 3]
        ps1 = PauliSum.from_random(2, dimensions)
        ps2 = PauliString.from_exponents([0, 1], [1, 0], dimensions)
        with pytest.raises(ValueError):
            _ = PauliSum.from_pauli_strings([ps1, ps2])  # type: ignore

        with pytest.raises(ValueError):
            _ = PauliSum.from_pauli_strings([])

    def test_pauli_sum_from_pauli_objects(self):
        ps = PauliSum.from_random(3, 3, rand_phases=True)
        pauli_objects = [
            PauliString.from_string('x1z0', dimensions=3),
            PauliString.from_string('x1z2', dimensions=3),
            ps
        ]
        phases = np.concatenate([np.array([0, 0], dtype=int), ps.phases])
        P = PauliSum.from_pauli_objects(pauli_objects, inherit_phases=True)
        assert P.shape() == (5, 1)
        assert np.array_equal(P.phases, phases)

        ps1 = PauliSum.from_string(['x0z1 x1z2 x4z4', 'x0z0 x2z2 x3z0'], [2, 3, 5], weights=[2, 0.5], phases=[1, 0])
        ps2 = PauliSum.from_string(['x0z1 x1z1 x0z4', 'x0z0 x2z1 x2z1', 'x1z0 x2z1 x1z1'],
                                   [2, 3, 5], weights=[2, 3.25, 1j], phases=[0, 11, 2])

        P = PauliSum.from_pauli_objects(ps1, inherit_weights=True, inherit_phases=True)
        P2 = PauliSum.from_string(
            ['x0z1 x1z2 x4z4', 'x0z0 x2z2 x3z0'],
            [2, 3, 5], weights=[2, 0.5], phases=[1, 0])
        assert P == P2

        P = PauliSum.from_pauli_objects([ps1, ps2], inherit_weights=True, inherit_phases=True)
        P2 = PauliSum.from_string(
            ['x0z1 x1z2 x4z4', 'x0z0 x2z2 x3z0', 'x0z1 x1z1 x0z4', 'x0z0 x2z1 x2z1', 'x1z0 x2z1 x1z1'],
            [2, 3, 5], weights=[2, 0.5, 2, 3.25, 1j], phases=[1, 0, 0, 11, 2])
        assert P == P2

        with pytest.warns(UserWarning, match="Phases are disregarded"):
            _ = PauliSum.from_pauli_objects(ps1, phases=[0, 0], inherit_phases=True)

        with pytest.warns(UserWarning, match="Weights are disregarded"):
            _ = PauliSum.from_pauli_objects(ps1, weights=[0, 0], inherit_weights=True)

    def test_pauli_sum_from_pauli_objects_exceptions(self):
        with pytest.raises(ValueError):
            _ = PauliSum.from_pauli_objects([])

        ps1 = PauliSum.from_string(['x0z1 x1z2 x4z4', 'x0z0 x2z2 x3z0'], [2, 3, 5])
        ps2 = PauliSum.from_string(['x0z1 x1z1 x0z4', 'x0z0 x2z1 x2z1', 'x1z0 x2z1 x1z1'],
                                   [2, 3, 7])
        with pytest.raises(ValueError):
            _ = PauliSum.from_pauli_objects([ps1, ps2])

    def test_pauli_sum_from_string(self):
        ps = PauliString.from_exponents([1, 2], [3, 1], [5, 7])
        s = str(ps)  # Use the same formatting
        P = PauliSum.from_string(s, [5, 7])

        assert P.tableau.shape[0] == 1

    def test_pauli_sum_from_string_list(self):
        s1 = "x1z0 x0z2"
        P = PauliSum.from_string([s1], 3)
        assert P.tableau.shape[0] == 1
        assert P.dimensions.tolist() == [3, 3]

        s1 = "x1z0 x0z1"
        P = PauliSum.from_string([s1])
        assert P.tableau.shape[0] == 1
        assert P.dimensions.tolist() == [DEFAULT_QUDIT_DIMENSION, DEFAULT_QUDIT_DIMENSION]

    def test_pauli_sum_from_tableau_roundtrip(self):
        ps = PauliString.from_exponents([1, 0], [2, 0], [5, 7])
        P = PauliSum.from_pauli_strings(ps)
        P2 = PauliSum.from_tableau(P.tableau, P.dimensions, P.weights, P.phases)

        assert np.array_equal(P2.tableau, P.tableau)
        assert np.array_equal(P2.weights, P.weights)
        assert np.array_equal(P2.phases, P.phases)

    def test_pauli_sum_from_random_small_population_no_duplicates(self):
        dims = [2, 3]  # max_n_paulis = 4 * 9 = 36 < 1e6
        N = math.prod([d**2 for d in dims])
        P = PauliSum.from_random(N, dims, seed=123)
        P.combine_equivalent_paulis()

        assert P.n_paulis() == N

    def test_pauli_sum_from_random_large_population_allows_duplicates(self):
        dims = [2, 3, 5, 7, 11, 13, 17]
        N = 300
        P = PauliSum.from_random(N, dims, seed=123)
        P.combine_equivalent_paulis()

        assert P.n_paulis() <= N

    def test_pauli_sum_from_random_reproducibility(self):
        dims = [3, 5]
        P1 = PauliSum.from_random(10, dims, seed=555)
        P2 = PauliSum.from_random(10, dims, seed=555)

        # Same seeds should give identical tableau/weights/phases
        assert np.array_equal(P1.tableau, P2.tableau)
        assert np.allclose(P1.weights, P2.weights)
        assert np.array_equal(P1.phases, P2.phases)

    def test_pauli_sum_from_random_too_many(self):
        dims = [2, 2]  # max = 16
        with pytest.raises(ValueError):
            PauliSum.from_random(17, dims)

    def test_pauli_sum_from_random_decoding_correctness_small_population(self):
        dims = [2, 3]
        n = 5
        P = PauliSum.from_random(n, dims, seed=999)

        pauli_strings = [P.select_pauli_string(i) for i in range(P.n_paulis())]
        for ps in pauli_strings:
            assert len(ps.x_exp) == len(dims)
            assert len(ps.z_exp) == len(dims)
            for xe, ze, d in zip(ps.x_exp, ps.z_exp, dims):
                assert 0 <= xe < d
                assert 0 <= ze < d

    def test_pauli_sum_to_and_from_file(self, tmp_path: Path):
        path = tmp_path / "paulisum.data"
        dimensions = [2, 3, 5, 7, 11]

        for _ in range(50):
            P = PauliSum.from_random(10, dimensions)
            P.to_file(path)
            P_from_file = PauliSum.from_file(path, dimensions)
            assert P == P_from_file

    def test_pauli_sum_remove_trivial_qudits(self):
        P = PauliSum.from_string(["x0z2 x1z1 x0z0 x2z3", "x0z1 x0z0 x0z0 x0z1"], dimensions=7)
        P.remove_trivial_qudits()
        assert P == PauliSum.from_string(["x0z2 x1z1 x2z3", "x0z1 x0z0 x0z1"], dimensions=7)

    def test_pauli_sum_remove_trivial_paulis(self):
        P = PauliSum.from_string(["x0z0 x0z0 x0z0 x0z0", "x0z1 x0z0 x0z0 x0z1"], dimensions=7)
        P.remove_trivial_paulis()
        assert P == PauliSum.from_string(["x0z1 x0z0 x0z0 x0z1"], dimensions=7)

    def test_pauli_sum_remove_zero_weight_paulis(self):
        P = PauliSum.from_string(["x0z2 x1z1 x0z0 x2z3", "x0z1 x0z0 x0z0 x0z1"], dimensions=7, weights=[0, 1e-13])
        P.remove_zero_weight_paulis()
        assert P == PauliSum.from_string(["x0z1 x0z0 x0z0 x0z1"], dimensions=7) / 1e13

        P = PauliSum.from_string(["x0z2 x1z1 x0z0 x2z3", "x0z1 x0z0 x0z0 x0z1"], dimensions=7, weights=[1, 1e-14])
        P.remove_zero_weight_paulis()
        assert P == PauliSum.from_string(["x0z2 x1z1 x0z0 x2z3"], dimensions=7)

    def test_pauli_sum_weight_to_phase(self):
        strings = ["x0z2 x1z2 x0z0 x2z2 x3z0", "x0z0 x3z2 x1z1 x0z0 x0z1", "x0z1 x0z0 x0z0 x0z1 x2z4"]
        dimensions = [2, 3, 2, 5, 3]
        P = PauliSum.from_string(strings, dimensions, weights=[1j, 0, 1 + 2j], phases=[0, 1, 5])
        P.weight_to_phase()
        Pres = PauliSum.from_string(strings, dimensions, weights=[1, 0, 2.233 - 0.1j], phases=[15, 1, 16])
        assert P.is_close(Pres, threshold=3)

    def test_pauli_string_from_int(self):
        ps = PauliString.from_exponents(2, 1, 3)
        assert ps.shape() == (1, 1)

        P = PauliSum.from_pauli_strings(ps, weights=2.0, phases=[1])
        assert P.tableau.shape == (1, 2)
        assert P.weights[0] == 2.0
        assert P.phases[0] == 1

    def test_pauli_string_from_int_mismatch(self):
        with pytest.raises(ValueError):
            _ = PauliString.from_exponents(2, 1, [3, 2])

        with pytest.raises(ValueError):
            _ = PauliString.from_exponents(2, [1, 2], 3)

    def test_pauli_string_from_tableau(self):
        tableau = np.asarray([0, 0, 1, 1], dtype=int)
        _ = PauliString.from_tableau(tableau)

        tableau = np.asarray([[0, 0, 1, 1]], dtype=int)
        _ = PauliString.from_tableau(tableau)

        tableau = [0, 0, 1, 1]
        _ = PauliString.from_tableau(tableau)

        with pytest.raises(ValueError):
            tableau = np.asarray([[[0, 0, 1, 1]]], dtype=int)
            _ = PauliString.from_tableau(tableau)

        with pytest.raises(ValueError):
            tableau = np.asarray([[0, 0, 1, 1], [1, 1, 0, 1]], dtype=int)
            _ = PauliString.from_tableau(tableau)

        for _ in range(100):
            dimension = random.choices(PRIME_LIST)[0]
            x = random.randint(0, dimension - 1)
            z = random.randint(0, dimension - 1)
            tableau = [x, z]
            p = PauliString.from_tableau(tableau, dimension)
            assert f"{p}" == f"x{x}z{z} "

    def test_pauli_string_from_random(self):
        ps1 = PauliString.from_random(3)
        ps2 = PauliString.from_random(3, 42)
        assert ps1.shape() == ps2.shape()

        ps3 = PauliString.from_random([2, 3, 5])
        ps4 = PauliString.from_random([2, 3, 5], 12345)
        assert ps3.shape() == ps4.shape()
