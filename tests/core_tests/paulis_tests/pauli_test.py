import numpy as np
import random
import pytest
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.paulis import PauliSum, PauliString, Pauli
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION
from tests import PRIME_LIST, choose_random_dimensions


N_tests = 30


class TestPaulis:

    def test_pauli_multiplication(self):
        for dim in PRIME_LIST:
            x1 = Pauli.Xnd(1, dim)
            y1 = Pauli.Ynd(1, dim)
            z1 = Pauli.Znd(1, dim)
            id = Pauli.Idnd(dim)

            # REMARK: phases do not matter, since these are Pauli objects
            assert x1 * z1 == y1, 'Error in Pauli multiplication (x * z = y) ' + (x1 * z1).__str__()
            assert x1**dim == id, 'Error in Pauli exponentiation (x**dim = id) ' + (x1**dim).__str__()
            assert y1**dim == id, 'Error in Pauli exponentiation (y**dim = id) ' + (y1**dim).__str__()
            assert z1**dim == id, 'Error in Pauli exponentiation (z**dim = id)  ' + (z1**dim).__str__()
            assert x1 * y1 == x1**2 * z1, 'Error in Pauli multiplication (x * y = x**2 * z) ' + (x1 * y1).__str__()
            assert y1 * z1 == x1 * z1**2, 'Error in Pauli multiplication (y * z = x**2 * z) ' + (y1 * z1).__str__()
            assert z1 * x1 == y1, 'Error in Pauli multiplication (z * x = y) ' + (z1 * x1).__str__()
            assert x1 * id == x1, 'Error in Pauli multiplication (x * id = x) ' + (x1 * id).__str__()
            assert y1 * id == y1, 'Error in Pauli multiplication (y * id = y) ' + (y1 * id).__str__()
            assert z1 * id == z1, 'Error in Pauli multiplication (z * id = z) ' + (z1 * id).__str__()

        for dim in PRIME_LIST:
            for _ in range(N_tests):
                s1 = np.random.randint(0, dim)
                r1 = np.random.randint(0, dim)
                s2 = np.random.randint(0, dim)
                r2 = np.random.randint(0, dim)
                p1 = Pauli.from_exponents(r1, s1, dim)
                p2 = Pauli.from_exponents(r2, s2, dim)
                p3 = p1 * p2
                assert p3.x_exp == (p1.x_exp + p2.x_exp) % dim, 'Error in Pauli multiplication (x_exp)'
                assert p3.z_exp == (p1.z_exp + p2.z_exp) % dim, 'Error in Pauli multiplication (z_exp)'
                assert p3.dimension == dim, 'Error in Pauli multiplication (dimension)'

    def test_pauli_string_multiplication(self):
        for dim in PRIME_LIST:
            for _ in range(N_tests):
                r11 = np.random.randint(0, dim)
                r12 = np.random.randint(0, dim)
                s11 = np.random.randint(0, dim)
                s12 = np.random.randint(0, dim)

                r21 = np.random.randint(0, dim)
                r22 = np.random.randint(0, dim)
                s21 = np.random.randint(0, dim)
                s22 = np.random.randint(0, dim)

                input_str1 = f"x{r11}z{s11} x{r12}z{s12}"
                input_str2 = f"x{r21}z{s21} x{r22}z{s22}"

                left = f"x{(r11 + r21) % dim}z{(s11 + s21) % dim}"
                right = f"x{(r12 + r22) % dim}z{(s12 + s22) % dim}"
                output_str_correct = f"{left} {right}"

                input_ps1 = PauliString.from_string(input_str1, dimensions=[dim, dim])
                input_ps2 = PauliString.from_string(input_str2, dimensions=[dim, dim])
                output_ps = input_ps1 * input_ps2

                assert output_ps.has_equal_tableau(PauliString.from_string(
                    output_str_correct, dimensions=[dim, dim]
                )), 'Error in PauliString multiplication'

    def test_pauli_string_tensor_product(self):
        for dimensions in PRIME_LIST:
            for _ in range(N_tests):
                r1 = np.random.randint(0, dimensions)
                r2 = np.random.randint(0, dimensions)
                s1 = np.random.randint(0, dimensions)
                s2 = np.random.randint(0, dimensions)

                r12 = np.random.randint(0, dimensions)
                r22 = np.random.randint(0, dimensions)
                s12 = np.random.randint(0, dimensions)
                s22 = np.random.randint(0, dimensions)

                input_str1 = f"x{r1}z{s1} x{r2}z{s2}"
                input_str2 = f"x{r12}z{s12} x{r22}z{s22}"
                output_str_correct = f"x{r1}z{s1} x{r2}z{s2} x{r12}z{s12} x{r22}z{s22}"
                input_ps1 = PauliString.from_string(input_str1, dimensions)
                input_ps2 = PauliString.from_string(input_str2, dimensions)
                output_ps = input_ps1 @ input_ps2

                assert output_ps == PauliString.from_string(
                    output_str_correct, dimensions), 'Error in PauliString tensor product'

    def test_pauli_string_indexing(self):
        for dim in PRIME_LIST:
            for _ in range(N_tests):

                r1 = np.random.randint(0, dim)
                r2 = np.random.randint(0, dim)
                s1 = np.random.randint(0, dim)
                s2 = np.random.randint(0, dim)

                p_string1 = PauliString.from_string(f"x{r1}z{s1} x{r2}z{s2}", dimensions=[dim, dim])
                ps0 = Pauli.from_string(f"x{r1}z{s1}", dimension=dim)
                ps1 = Pauli.from_string(f"x{r2}z{s2}", dimension=dim)

                assert p_string1[0] == ps0, 'Error in PauliString indexing (first PauliString)'
                assert p_string1[1] == ps1, 'Error in PauliString indexing'

    def test_pauli_sum_multiplication(self):
        for dim in PRIME_LIST:
            for _ in range(N_tests):
                # Notice I generate numbers beyond dim to check modulo operations :)
                r11 = np.random.randint(0, 2 * dim)
                r21 = np.random.randint(0, 2 * dim)
                s11 = np.random.randint(0, 2 * dim)
                s21 = np.random.randint(0, 2 * dim)
                p_string1 = PauliString.from_string(f"x{r11}z{s11} x{r21}z{s21}", dimensions=[dim, dim])

                r12 = np.random.randint(0, 2 * dim)
                r22 = np.random.randint(0, 2 * dim)
                s12 = np.random.randint(0, 2 * dim)
                s22 = np.random.randint(0, dim)
                p_string2 = PauliString.from_string(f"x{r12}z{s12} x{r22}z{s22}", dimensions=[dim, dim])

                random_pauli_sum = PauliSum.from_pauli_strings([p_string1, p_string2])

                # Test multiplication of PauliSum with PauliString
                r13 = np.random.randint(0, 2 * dim)
                r23 = np.random.randint(0, 2 * dim)
                s13 = np.random.randint(0, 2 * dim)
                s23 = np.random.randint(0, 2 * dim)
                input_ps1 = PauliString.from_string(
                    f"x{r13}z{s13} x{r23}z{s23}",
                    dimensions=[dim, dim]
                )

                output_str_correct = (
                    f"x{(r11 + r13) % dim}z{(s11 + s13) % dim} "
                    f"x{(r21 + r23) % dim}z{(s21 + s23) % (2 * dim)}"
                )
                output_str2_correct = (
                    f"x{(r12 + r13) % dim}z{(s12 + s13) % dim} "
                    f"x{(r22 + r23) % dim}z{(s22 + s23) % (2 * dim)}"
                )
                output_ps = random_pauli_sum * input_ps1

                phase1 = 2 * (s11 * r13 + s21 * r23) % (2 * dim)
                phase2 = 2 * (s12 * r13 + s22 * r23) % (2 * dim)
                output_phases = [phase1, phase2]
                output_ps_correct = PauliSum.from_pauli_strings(
                    [PauliString.from_string(output_str_correct, dimensions=[dim, dim]),
                     PauliString.from_string(output_str2_correct, dimensions=[dim, dim])],
                    phases=output_phases)

                assert output_ps == output_ps_correct, (
                    'Error in PauliSum multiplication with PauliString\n' +
                    output_ps.__str__() +
                    '\n' +
                    output_ps_correct.__str__()
                )

                # Test multiplication of PauliSum with PauliSum
                random_pauli_sum2 = PauliSum.from_pauli_strings([input_ps1])
                output_ps2 = random_pauli_sum * random_pauli_sum2

                assert output_ps2 == output_ps_correct, 'Error in PauliSum multiplication with PauliSum'

    def test_pauli_sum_tensor_product(self):
        for dim in PRIME_LIST:
            for _ in range(N_tests):
                r11 = np.random.randint(0, 2 * dim)
                s11 = np.random.randint(0, 2 * dim)
                r12 = np.random.randint(0, 2 * dim)
                s12 = np.random.randint(0, 2 * dim)
                p_string1 = PauliString.from_string(f"x{r11}z{s11} x{r12}z{s12}", dimensions=[dim, dim])

                r21 = np.random.randint(0, 2 * dim)
                s21 = np.random.randint(0, 2 * dim)
                r22 = np.random.randint(0, 2 * dim)
                s22 = np.random.randint(0, 2 * dim)
                p_string2 = PauliString.from_string(f"x{r21}z{s21} x{r22}z{s22}", dimensions=[dim, dim])

                r31 = np.random.randint(0, 2 * dim)
                s31 = np.random.randint(0, 2 * dim)
                r32 = np.random.randint(0, 2 * dim)
                s32 = np.random.randint(0, 2 * dim)
                p_string3 = PauliString.from_string(f"x{r31}z{s31} x{r32}z{s32}", dimensions=[dim, dim])

                random_pauli_sum = PauliSum.from_pauli_strings([p_string1, p_string2])
                random_pauli_sum2 = PauliSum.from_pauli_strings([p_string3])
                ps_out = random_pauli_sum @ random_pauli_sum2
                ps_out_correct = PauliSum.from_pauli_strings(
                    [p_string1 @ p_string3, p_string2 @ p_string3])
                assert ps_out == ps_out_correct, 'Error in PauliSum tensor product'

    def test_combine_equivalent_paulis_merges_duplicates(self):
        """
        Duplicate tableau rows should be merged even when their weights differ, with
        the resulting weight equal to the sum of the duplicates.
        """
        rng = np.random.default_rng()
        dims = [DEFAULT_QUDIT_DIMENSION, DEFAULT_QUDIT_DIMENSION]

        for _ in range(25):
            x_exp = rng.integers(0, dims[0], size=len(dims))
            z_exp = rng.integers(0, dims[0], size=len(dims))

            ps = PauliString.from_exponents(x_exp, z_exp, dims)

            # Two identical PauliStrings with different weights
            weights = np.array([rng.normal(), rng.normal()])
            P = PauliSum.from_pauli_strings([ps, ps], weights=weights, phases=[0, 0])

            P.combine_equivalent_paulis()

            assert P.n_paulis() == 1, "Duplicates were not collapsed"
            assert np.isclose(P.weights[0], weights.sum()), "Weights did not sum when collapsing duplicates"
            assert np.array_equal(P.tableau[0], ps.tableau[0]), "Tableau changed when collapsing duplicates"

    def test_symplectic_matrix_single_pauli(self):
        for _ in range(N_tests):
            dim = random.choice(PRIME_LIST)
            expected_tableau = np.array([[random.randint(0, dim - 1), random.randint(0, dim - 1)]])
            pauli_list = [f'x{expected_tableau[0, 0]}z{expected_tableau[0, 1]}']
            weights = np.array([1.])
            ps = PauliSum.from_string(pauli_list, dimensions=[dim], weights=weights)

            np.testing.assert_array_equal(ps.tableau, expected_tableau)

    def test_symplectic_matrix_multiple_paulis(self):
        for _ in range(N_tests):
            dimensions = choose_random_dimensions(max_product=30)
            n_paulis = 10
            # choose tableau
            expected_tableau = np.array([
                (lambda x, z: x + z)(
                    [random.randint(0, dim - 1) for dim in dimensions],
                    [random.randint(0, dim - 1) for dim in dimensions]
                )
                for _ in range(n_paulis)
            ], dtype=int)
            # create pauli_list from tableau
            pauli_list = [
                ' '.join([f'x{expected_tableau[i, j]}z{expected_tableau[i, len(dimensions) + j]}'
                          for j in range(len(dimensions))])
                for i in range(n_paulis)
            ]
            # randomly generate weights
            weights = np.array([random.normalvariate(0, 1) for _ in range(n_paulis)])
            # PauliSum
            ps = PauliSum.from_string(pauli_list, dimensions=dimensions, weights=weights)

            np.testing.assert_array_equal(ps.tableau, expected_tableau)

    def test_basic_pauli_relations(self):
        for d in PRIME_LIST:
            x_exp = random.randint(1, d - 1)
            z_exp = random.randint(1, d - 1)
            x1 = Pauli.from_string(f'x{x_exp}z0', dimension=d)
            z1 = Pauli.from_string(f'x0z{z_exp}', dimension=d)
            y1 = Pauli.from_string(f'x{x_exp}z{z_exp}', dimension=d)
            id = Pauli.Idnd(dimension=d)

            assert x1 * z1 == y1, f'Error in Pauli multiplication for d={d}'
            assert x1**d == id, f'Error in Pauli exponentiation (x**{d} = id) for d={d}'
            assert y1**d == id, f'Error in Pauli exponentiation (y**{d} = id) for d={d}'
            assert z1**d == id, f'Error in Pauli exponentiation (z**{d} = id) for d={d}'
            assert x1 * id == x1, f'Error in Pauli multiplication (x * id = x) for d={d}'
            assert id * z1 == z1, f'Error in Pauli multiplication (id * z = z) for d={d}'

    def test_pauli_string_construction(self):
        for _ in range(N_tests):
            dimensions = choose_random_dimensions(max_product=60)

            # test 1
            exps = [[random.randint(0, dim - 1) for dim in dimensions],
                    [random.randint(0, dim - 1) for dim in dimensions]]
            string = ' '.join(f'x{xe}z{ze}' for xe, ze in zip(exps[0], exps[1]))
            s_1 = PauliString.from_string(string, dimensions)
            s_2 = PauliString.from_exponents(exps[0], exps[1], dimensions)

            assert s_1 == s_2

        # test 2
        dims = [3, 3]
        x1y0 = PauliString.from_exponents([1, 1], [1, 0])
        x1y0_2 = PauliString.from_string('x1z1 x1z0', dimensions=DEFAULT_QUDIT_DIMENSION)
        x1y0_3 = PauliString.from_string('x1z1 x1z0', dimensions=dims)

        assert x1y0 == x1y0_2
        assert x1y0 != x1y0_3

        # Test 3 with minimum allowed qudit dimension.
        with pytest.raises(ValueError):
            _ = PauliString.from_string('x0z0 x0z0', dimensions=DEFAULT_QUDIT_DIMENSION - 1)

    def test_pauli_sum_addition(self):

        for dim in PRIME_LIST:
            for _ in range(N_tests):
                r11 = np.random.randint(0, 2 * dim)
                s11 = np.random.randint(0, 2 * dim)
                r12 = np.random.randint(0, 2 * dim)
                s12 = np.random.randint(0, 2 * dim)
                p_string1 = PauliString.from_string(f"x{r11}z{s11} x{r12}z{s12}", dimensions=[dim, dim])

                r21 = np.random.randint(0, 2 * dim)
                s21 = np.random.randint(0, 2 * dim)
                r22 = np.random.randint(0, 2 * dim)
                s22 = np.random.randint(0, 2 * dim)
                p_string2 = PauliString.from_string(f"x{r21}z{s21} x{r22}z{s22}", dimensions=[dim, dim])

                r31 = np.random.randint(0, 2 * dim)
                s31 = np.random.randint(0, 2 * dim)
                r32 = np.random.randint(0, 2 * dim)
                s32 = np.random.randint(0, 2 * dim)
                p_string3 = PauliString.from_string(f"x{r31}z{s31} x{r32}z{s32}", dimensions=[dim, dim])

                random_pauli_sum_1 = PauliSum.from_pauli_strings([p_string1, p_string2])
                random_pauli_sum_2 = PauliSum.from_pauli_strings([p_string3])
                ps_out = random_pauli_sum_1 + random_pauli_sum_2
                ps_out_correct = PauliSum.from_pauli_strings([p_string1, p_string2, p_string3])
                assert ps_out == ps_out_correct, 'Error in PauliSum addition'

        for _ in range(N_tests):
            dimensions = choose_random_dimensions(max_product=30)
            n_paulis = 10

            # Generate random PauliStrings
            pauli_strings = []
            for _ in range(n_paulis):
                # exponents intentionally out of bounds to test modulus operation
                x_exps = [np.random.randint(0, 2 * dim) for dim in dimensions]
                z_exps = [np.random.randint(0, 2 * dim) for dim in dimensions]
                ps = PauliString.from_exponents(x_exps, z_exps, dimensions=dimensions)
                pauli_strings.append(ps)

            # Create PauliSum by adding individual PauliSums
            psum = PauliSum.from_pauli_strings([pauli_strings[0]])
            for ps in pauli_strings[1:]:
                psum = psum + PauliSum.from_pauli_strings([ps])
            psum.standardise()

            # Expected result: all PauliStrings with default weights and phases
            expected = PauliSum.from_pauli_strings(pauli_strings, weights=[1 + 0j] * n_paulis, phases=[0] * n_paulis)
            expected.standardise()

        assert isinstance(psum, PauliSum) and isinstance(expected, PauliSum), (
            f"PauliSum used are not of type PauliSum, being them \n{psum}\n and \n{expected}\n,"
            f"with dimensions {dimensions}"
        )
        assert psum == expected, (
            f"PauliSum addition failed, \n obtained \n{psum}\n expected \n{expected}\n,"
            f"with dimensions {dimensions}"
        )

    def test_phase_and_dot_product(self):

        for _ in range(N_tests):
            d = random.choice(PRIME_LIST)
            x = PauliString.from_string('x1z0', dimensions=[d])
            z = PauliString.from_string('x0z1', dimensions=[d])

            assert z.acquired_phase(x) == 2.0, 'Expected phase to be 2.0, got {}'.format(z.acquired_phase(x))
            assert x.acquired_phase(z) == 0.0, 'Expected phase to be 0.0, got {}'.format(x.acquired_phase(z))

        for _ in range(N_tests):
            dims = [random.choice(PRIME_LIST) for _ in range(2)]
            x1x1 = PauliSum.from_pauli_strings(PauliString.from_string('x1z0 x1z0', dimensions=dims))
            x1y1 = PauliSum.from_pauli_strings(PauliString.from_string('x1z0 x1z1', dimensions=dims))

            # define random weight
            random_weight = random.normalvariate(0, 1)
            # phase correction for this specific case
            phases_correction = 2 * x1x1.lcm // dims[-1]

            s1 = x1x1 + x1y1 * random_weight
            s2 = x1x1 + x1x1

            s3 = PauliSum.from_string(['x2z0 x2z0', 'x2z0 x2z0', 'x2z0 x2z1', 'x2z0 x2z1'],
                                      weights=[1, 1, random_weight, random_weight],
                                      phases=[0, 0, phases_correction, phases_correction],
                                      dimensions=dims)

            assert (
                s1 * s2).has_equal_tableau(s3), f'Expected s1 * s2 to equal s3, got \n{s1 * s2}\n instead of \n{s3}\n'

    def test_tensor_product_distributivity(self):
        for _ in range(N_tests):
            dimensions = choose_random_dimensions(max_product=30)

            tmp1 = PauliString.from_exponents(
                x_exp=[random.randint(0, d - 1) for d in dimensions],
                z_exp=[random.randint(0, d - 1) for d in dimensions],
                dimensions=dimensions
            )
            tmp2 = PauliString.from_exponents(
                x_exp=[random.randint(0, d - 1) for d in dimensions],
                z_exp=[random.randint(0, d - 1) for d in dimensions],
                dimensions=dimensions
            )

            p1 = PauliSum.from_pauli_strings(tmp1)
            p2 = PauliSum.from_pauli_strings(tmp2)

            s1 = p1 * random.normalvariate(0, 1) + p2 * random.normalvariate(0, 1)
            s2 = p2 * random.normalvariate(0, 1) + p2 * random.normalvariate(0, 1)

            left = (s1 + s2) @ s2
            right = s1 @ s2 + s2 @ s2

            assert left == right

    def test_pauli_sum_indexing(self):
        for _ in range(N_tests):
            dimensions = choose_random_dimensions(max_product=60)
            n_paulis = 10
            n_qudits = len(dimensions)

            # Generate random PauliStrings
            pauli_strings: list[PauliString] = []
            for _ in range(n_paulis):
                ps = PauliString.from_exponents(
                    x_exp=[random.randint(0, d - 1) for d in dimensions],
                    z_exp=[random.randint(0, d - 1) for d in dimensions],
                    dimensions=dimensions
                )
                pauli_strings.append(ps)

            # Generate random weights and phases
            weights = [random.normalvariate(0, 1) + 1j * random.normalvariate(0, 1) for _ in range(n_paulis)]
            phases = [random.randint(0, 2 * np.lcm.reduce(dimensions) - 1) for _ in range(n_paulis)]

            # Create PauliSum
            ps = PauliSum.from_pauli_strings(pauli_strings, weights=weights, phases=phases)

            # Test single indexing
            for i in range(n_paulis):
                assert ps.select_pauli_string(i).has_equal_tableau(pauli_strings[i]), \
                    f'Error in PauliSum single indexing at position {i}'

            # Test slice indexing
            if n_paulis >= 2:
                slice_idx = slice(0, min(2, n_paulis))
                expected_slice = PauliSum.from_pauli_strings(
                    pauli_strings[slice_idx],
                    weights=weights[slice_idx],
                    phases=phases[slice_idx]
                )
                assert ps[slice_idx].has_equal_tableau(expected_slice), 'Error in PauliSum slice indexing'

            # Test list indexing (select multiple Pauli strings)
            if n_paulis >= 2:
                idx_list = [0, min(n_paulis - 1, random.randint(1, n_paulis - 1))]
                expected_list = PauliSum.from_pauli_strings(
                    [pauli_strings[i] for i in idx_list],
                    weights=[weights[i] for i in idx_list],
                    phases=[phases[i] for i in idx_list]
                )
                assert ps[idx_list].has_equal_tableau(expected_list), 'Error in PauliSum list indexing'

            # Test combined indexing (select Pauli strings and qudits)
            if n_paulis >= 2 and n_qudits >= 1:
                idx_list = [0, min(n_paulis - 1, random.randint(1, n_paulis - 1))]
                qudit_idx = random.randint(0, n_qudits - 1)

                expected_combined = PauliSum.from_pauli_objects(
                    [pauli_strings[i][qudit_idx] for i in idx_list],
                    weights=[weights[i] for i in idx_list],
                    phases=[phases[i] for i in idx_list]
                )
                assert ps[idx_list, qudit_idx].has_equal_tableau(expected_combined), \
                    'Error in PauliSum combined indexing (list, single qudit)'

            # Test combined indexing with multiple qudits
            if n_paulis >= 2 and n_qudits >= 2:
                idx_list = [0, min(n_paulis - 1, random.randint(1, n_paulis - 1))]
                qudit_list = [0, random.randint(1, n_qudits - 1)]

                expected_multi = PauliSum.from_pauli_strings(
                    [PauliString.from_exponents(
                        x_exp=[pauli_strings[i].x_exp[j] for j in qudit_list],
                        z_exp=[pauli_strings[i].z_exp[j] for j in qudit_list],
                        dimensions=[dimensions[j] for j in qudit_list]
                    ) for i in idx_list],
                    weights=[weights[i] for i in idx_list],
                    phases=[phases[i] for i in idx_list]
                )
                assert ps[idx_list, qudit_list].has_equal_tableau(expected_multi), \
                    'Error in PauliSum combined indexing (list, list)'

    def test_pauli_sum_amend(self):
        for _ in range(N_tests):
            dims = random.choices(PRIME_LIST, k=2)
            while dims[-1] == 2:
                dims = random.choices(PRIME_LIST, k=2)
            # p1 = X on qubit, p2 = Z on qutrit, p3 = XZ on qutrit
            p1 = PauliString.from_exponents(x_exp=[1, 0], z_exp=[0, 0], dimensions=dims)
            p2 = PauliString.from_exponents(x_exp=[0, 0], z_exp=[0, 1], dimensions=dims)
            p3 = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 1], dimensions=dims)
            ps = PauliSum.from_pauli_strings([p1, p2, p3])

            ps[0] = ps.select_pauli_string(0).amend(0, 0, 1)
            ps[1] = ps.select_pauli_string(1).amend(1, 1, 1)
            ps[2] = ps.select_pauli_string(2).amend(1, 1, 2)  # Try assignment via __setitem__

            # p1 = Z on qubit, p2 = XZ on qutrit, p3 = XZ^2 on qutrit
            new_p1 = PauliString.from_exponents(x_exp=[0, 0], z_exp=[1, 0], dimensions=dims)
            new_p2 = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 1], dimensions=dims)
            new_p3 = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 2], dimensions=dims)
            new_ps = PauliSum.from_pauli_strings([new_p1, new_p2, new_p3])

            assert ps.select_pauli_string(0) == new_p1
            assert ps.select_pauli_string(1) == new_p2
            assert ps.select_pauli_string(2) == new_p3

            assert ps == new_ps

    def test_ordering(self):
        # check that the symplectic basis gives the identity when ordered
        for _ in range(N_tests):
            n_qudits = random.randint(2, 50)
            d = random.choices(PRIME_LIST, k=n_qudits)
            n_paulis = 2 * n_qudits
            symplectic_basis = np.eye(2 * n_qudits, dtype=int)
            basis_list = [symplectic_basis[i, :] for i in range(n_paulis)]
            shuffled_basis = basis_list.copy()
            np.random.shuffle(shuffled_basis)
            ps = PauliSum.from_tableau(np.array(shuffled_basis), d)

            assert np.all(ps.to_standard_form().tableau == symplectic_basis)

    def test_pauli_sum_product_mixed_species(self):
        # Test multiplication
        P1 = PauliSum.from_string(['x1z1 x0z0'],
                                  dimensions=[3, 2],
                                  weights=[1], phases=[0])

        P2 = PauliSum.from_string(['x2z2 x0z0'],
                                  dimensions=[3, 2],
                                  weights=[1], phases=[0])

        product = P2 * P2
        assert product.phases == [4]

        product = P1.H() * P2
        assert product.phases == [8]

        dimensions = choose_random_dimensions(max_product=30)
        max_n_paulis = int(np.prod(dimensions)**2)
        for _ in range(int(np.ceil(np.sqrt(N_tests)))):
            for _ in range(int(np.ceil(np.sqrt(N_tests)))):
                P1 = PauliSum.from_random(random.randint(1, max_n_paulis), dimensions)
                P2 = PauliSum.from_random(random.randint(1, max_n_paulis), dimensions)
                P_res = P1 * P2
                phase_symplectic = P_res.phases[0]

                phase_computed = 0
                for j in range(len(dimensions)):
                    s1 = P1.z_exp[0, j]
                    r2 = P2.x_exp[0, j]
                    phase_computed += ((s1 * r2) % dimensions[j]) * P1.lcm / dimensions[j]
                phase_computed = phase_computed * 2 % (2 * P1.lcm)
                assert phase_symplectic == phase_computed

    def test_pauli_sum_delete_qudits(self):
        dims_pool = PRIME_LIST
        for _ in range(N_tests):
            n_qudits = random.randint(2, 8)
            dims = [random.choice(dims_pool) for _ in range(n_qudits)]
            n_terms = random.randint(2, 6)

            # Build original PauliStrings
            original_ps = []
            for _ in range(n_terms):
                x_exp = [random.randint(0, d - 1) for d in dims]
                z_exp = [random.randint(0, d - 1) for d in dims]
                original_ps.append(PauliString.from_exponents(x_exp=x_exp, z_exp=z_exp, dimensions=dims))

            psum = PauliSum.from_pauli_strings(original_ps)

            # Choose qudits to delete (at least one, at most n_qudits - 1)
            k_delete = random.randint(1, n_qudits - 1)
            delete_indices = sorted(random.sample(range(n_qudits), k_delete))

            # Perform deletion
            psum_deleted = psum.copy()
            psum_deleted._delete_qudits(delete_indices)

            # Build expected PauliSum
            keep_indices = [i for i in range(n_qudits) if i not in delete_indices]
            expected_dims = [dims[i] for i in keep_indices]
            expected_ps_list = []
            for ps in original_ps:
                x_keep = [ps.x_exp[i] for i in keep_indices]
                z_keep = [ps.z_exp[i] for i in keep_indices]
                expected_ps_list.append(PauliString.from_exponents(x_keep, z_keep, dimensions=expected_dims))
            expected_psum = PauliSum.from_pauli_strings(expected_ps_list)

            assert psum_deleted == expected_psum, (
                f"Deletion mismatch.\nDims={dims}\nDeleted={delete_indices}\n"
                f"Expected={expected_psum}\nGot={psum_deleted}"
            )
            assert np.all(expected_psum.dimensions == expected_dims), (
                f"Deletion mismatch.\nDims={dims}\nDeleted={delete_indices}\n"
                f"Expected dimensions={expected_dims}\nGot={psum_deleted.dimensions}"
            )
            assert all(len(ps.x_exp) == len(expected_dims) for ps in expected_psum), (
                f"Deletion mismatch.\nDims={dims}\nDeleted={delete_indices}\n"
                f"Expected number of qudits={len(expected_dims)}\nGot={len(ps.x_exp)}"
            )

            # Check sequential single-index deletion equals batch deletion
            psum_seq = PauliSum.from_pauli_strings(original_ps)
            for idx in sorted(delete_indices, reverse=True):
                psum_seq._delete_qudits([idx])
            assert psum_seq == expected_psum, "Sequential deletions differ from batch deletion"

    def test_symplectic_product(self):
        for _ in range(N_tests):
            dim = random.choices(PRIME_LIST, k=2)
            while dim[-1] < 3:
                dim = random.choices(PRIME_LIST, k=2)
            P1 = PauliString.from_string('x1z0', dimensions=[dim[0]])
            P2 = PauliString.from_string('x0z1', dimensions=[dim[0]])
            assert P1.symplectic_product(P2) == 1

            P1 = PauliString.from_string('x1z0', dimensions=[dim[0]])
            P2 = PauliString.from_string('x1z0', dimensions=[dim[0]])
            assert P1.symplectic_product(P2) == 0

            P1 = PauliString.from_string('x0z1', dimensions=[dim[0]])
            P2 = PauliString.from_string('x0z1', dimensions=[dim[0]])
            assert P1.symplectic_product(P2) == 0

            P1 = PauliString.from_string('x1z0 x1z0', dimensions=dim)
            P2 = PauliString.from_string('x0z1 x0z1', dimensions=dim)
            expected = 2 if dim[0] == dim[1] else np.sum(dim) % P1.lcm
            assert P1.symplectic_product(P2) == expected

            P1 = PauliString.from_string('x1z0 x0z1', dimensions=dim)
            P2 = PauliString.from_string('x1z0 x1z0', dimensions=dim)
            assert P1.symplectic_product(P2) == (dim[1] - 1) * P1.lcm / dim[1]

            P1 = PauliString.from_string('x1z0', dimensions=[dim[-1]])
            P2 = PauliString.from_string('x2z0', dimensions=[dim[-1]])
            assert P1.symplectic_product(P2) == 0

            P1 = PauliString.from_string('x1z2', dimensions=[dim[-1]])
            P2 = PauliString.from_string('x2z1', dimensions=[dim[-1]])
            assert P1.symplectic_product(P2) == (dim[-1] - 3) % dim[-1]

            P1 = PauliString.from_string('x1z1 x1z2', dimensions=dim)
            P2 = PauliString.from_string('x1z1 x2z1', dimensions=dim)
            assert P1.symplectic_product(P2) == (dim[-1] - 3) * P1.lcm / dim[1]

    def test_hermitian_generation(self):
        for _ in range(N_tests):
            dim = random.choices(PRIME_LIST, k=1)
            P1 = PauliString.from_string('x1z0', dimensions=dim)
            P2 = PauliString.from_string(f'x{dim[0] - 1}z0', dimensions=dim)
            assert P1.H().has_equal_tableau(P2)
        P1 = PauliString.from_string('x1z1', dimensions=dim)
        P2 = PauliString.from_string(f'x{dim[0] - 1}z{dim[0] - 1}', dimensions=dim)
        assert P1.H().has_equal_tableau(P2)

        P1 = PauliString.from_string(f'x{dim[0] - 1}z1', dimensions=dim)
        P2 = PauliString.from_string(f'x1z{dim[0] - 1}', dimensions=dim)
        assert P1.H().has_equal_tableau(P2)

        P1 = PauliString.from_string('x0z0', dimensions=dim)
        P2 = PauliString.from_string('x0z0', dimensions=dim)
        assert P1.H().has_equal_tableau(P2)
        assert P1.is_hermitian()

        dim = random.choices(PRIME_LIST, k=2)
        P1 = PauliString.from_string('x0z1 x1z1', dimensions=dim)
        P2 = PauliString.from_string(f'x0z{dim[0] - 1} x{dim[1] - 1}z{dim[1] - 1}', dimensions=dim)
        assert P1.H().has_equal_tableau(P2)

        P1 = PauliString.from_string('x0z1 x1z0', dimensions=dim)
        P2 = PauliString.from_string(f'x0z{dim[0] - 1} x{dim[1] - 1}z0', dimensions=dim)
        assert P1.H().has_equal_tableau(P2)

        dim = random.choices(PRIME_LIST, k=2)
        while dim[0] == 2:
            dim = random.choices(PRIME_LIST, k=2)
        P1 = PauliString.from_string('x2z1 x0z1', dimensions=dim)
        P2 = PauliString.from_string(f'x{dim[0] - 2}z{dim[0] - 1} x0z{dim[1] - 1}', dimensions=dim)
        assert P1.H().has_equal_tableau(P2)

    def test_pauli_sum_is_hermitian(self):
        for run in range(int(np.ceil(np.sqrt(N_tests)))):
            dimensions = choose_random_dimensions(max_product=30)
            D = int(np.prod(dimensions))

            # Generate random matrix
            H_e = (-1 + 2 * np.random.rand(D, D)) + 1j * (-1 + 2 * np.random.rand(D, D))
            # Make it Hermitian
            H_e = H_e + H_e.transpose().conjugate()
            assert np.array_equal(H_e, H_e.conjugate().transpose())

            # Sometimes add a non-Hermitian term
            if run % 3 == 0:
                H_e = H_e + 0.15 * (np.random.rand(D, D) + 1j * np.random.rand(D, D))
                assert not np.array_equal(H_e, H_e.conjugate().transpose())

            pauli_sum = PauliSum.from_hilbert_space(H_e, dimensions)
            tolerance = 1e-12
            assert np.max(np.abs((pauli_sum.to_hilbert_space().toarray() - H_e))) < tolerance
            assert pauli_sum.is_hermitian() == np.array_equal(H_e, H_e.conjugate().transpose())

    def test_qubit_XZ_phase_is_minus_one(self):
        # Single qubit (dimension 2): X * Z = (-1) Z * X  => scalar exponent r = 1 mod 2
        for _ in range(N_tests):
            dims = random.choices(PRIME_LIST, k=1)
            psX = PauliString.from_exponents(x_exp=[1], z_exp=[0], dimensions=dims)  # X
            psZ = PauliString.from_exponents(x_exp=[0], z_exp=[1], dimensions=dims)  # Z

            # per-qudit residue
            rj = psX.symplectic_residues(psZ)
            assert rj.tolist() == [1]

            # scalar phase mod LCM
            L = int(np.lcm.reduce(np.array(dims, int)))
            r_scalar = psX.symplectic_product(psZ, as_scalar=True)
            assert r_scalar == 1 % L

    def test_mixed_dims_XZ_phase(self):
        for _ in range(N_tests):
            # dims = [2, 3]; use site 1 (qutrit) to get omega_3
            dims = random.choices(PRIME_LIST, k=2)
            P = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 0], dimensions=dims)  # X on qutrit
            Q = PauliString.from_exponents(x_exp=[0, 0], z_exp=[0, 1], dimensions=dims)  # Z on qutrit

            # residues: [0 mod dim[0], 1 mod dim[1]]
            assert P.symplectic_residues(Q).tolist() == [0, 1]

            # scalar: r = 0*lcm/dim[0] + 1*lcm/dim[1] = lcm/dim[1] mod lcm
            assert P.symplectic_product(Q, as_scalar=True) == P.lcm // dims[1]

    def test_mixed_dims_all_sites_X_vs_Z_product(self):
        for _ in range(N_tests):
            dims = random.choices(PRIME_LIST, k=3)
            P = PauliString.from_exponents(x_exp=[1, 1, 1], z_exp=[0, 0, 0], dimensions=dims)  # X on all
            Q = PauliString.from_exponents(x_exp=[0, 0, 0], z_exp=[1, 1, 1], dimensions=dims)  # Z on all

            # residues per site are all 1
            assert P.symplectic_residues(Q).tolist() == [1, 1, 1]

            # scalar: r = lcm/dims[0] + lcm/dims[1] + lcm/dims[2]) mod lcm
            assert P.symplectic_product(Q, as_scalar=True) == (
                P.lcm // dims[0] + P.lcm // dims[1] + P.lcm // dims[2]) % P.lcm

    def test_bilinearity_scalar_mode(self):
        # Check <v1+v2, w> = <v1, w> + <v2, w>  (phase-preserving scalar)
        for _ in range(N_tests):
            dims = random.choices(PRIME_LIST, k=2)
            v1 = PauliString.from_exponents(x_exp=[1, 0], z_exp=[0, 0], dimensions=dims)  # X on qubit
            v2 = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 0], dimensions=dims)  # X on qutrit
            w = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 1], dimensions=dims)  # XZ on qutrit

            # Build v12 by adding tableaux modulo dims
            x12 = (np.array([1, 0]) + np.array([0, 1])) % np.array(dims)
            z12 = (np.array([0, 0]) + np.array([0, 0])) % np.array(dims)
            v12 = PauliString.from_exponents(x_exp=x12.tolist(), z_exp=z12.tolist(), dimensions=dims)

            L = int(np.lcm.reduce(np.array(dims, int)))
            lhs = v12.symplectic_product(w, as_scalar=True)
            rhs = (v1.symplectic_product(w, as_scalar=True) +
                   v2.symplectic_product(w, as_scalar=True)) % L
            assert lhs == rhs

    def test_antisymmetry_residues_and_scalar(self):
        for _ in range(N_tests):
            dims = random.choices(PRIME_LIST, k=2)
            P = PauliString.from_exponents(x_exp=[1, 0], z_exp=[0, 1], dimensions=dims)    # X0 Z1
            Q = PauliString.from_exponents(x_exp=[1, 1], z_exp=[0, 0], dimensions=dims)    # X0 X1

            rP_Q = P.symplectic_residues(Q)
            rQ_P = Q.symplectic_residues(P)
            dims_arr = np.array(dims, int)

            # r(P,Q) == - r(Q,P) (mod d_j) per qudit
            assert np.all((rP_Q + rQ_P) % dims_arr == 0)

            # scalar version mod L
            L = int(np.lcm.reduce(np.array(dims, int)))
            sPQ = P.symplectic_product(Q, as_scalar=True)
            sQP = Q.symplectic_product(P, as_scalar=True)
            assert (sPQ + sQP) % L == 0

    def test_symplectic_product_matrix_matches_pairwise_mixed_dims(self):
        for _ in range(N_tests):
            dims = random.choices(PRIME_LIST, k=2)
            # P1 = X on first qudit, P2 = Z on second qudit, P3 = XZ on second qudit
            P1 = PauliString.from_exponents(x_exp=[1, 0], z_exp=[0, 0], dimensions=dims)
            P2 = PauliString.from_exponents(x_exp=[0, 0], z_exp=[0, 1], dimensions=dims)
            P3 = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 1], dimensions=dims)
            S = PauliSum.from_pauli_strings([P1, P2, P3], weights=None, phases=None)

            SPM = S.symplectic_product_matrix()
            L = int(np.lcm.reduce(np.array(dims, int)))

            # Compute expected lower triangle via pairwise scalar products
            expect = np.zeros_like(SPM)
            ps = [P1, P2, P3]
            for i in range(len(ps)):
                for j in range(i):
                    expect[i, j] = ps[i].symplectic_product(ps[j], as_scalar=True)
            expect = (expect + expect.T) % L
            np.fill_diagonal(expect, 0)

            assert np.array_equal(SPM % L, expect % L)

    def test_symplectic_product_matrix_properties(self):
        for _ in range(N_tests):
            # Symmetry and zero diagonal
            dims = random.choices(PRIME_LIST, k=3)
            P1 = PauliString.from_exponents(x_exp=[1, 0, 0], z_exp=[0, 0, 0], dimensions=dims)  # X on first qudit
            P2 = PauliString.from_exponents(x_exp=[0, 1, 0], z_exp=[0, 1, 0], dimensions=dims)  # XZ on second qudit
            P3 = PauliString.from_exponents(x_exp=[0, 0, 1], z_exp=[0, 0, 0], dimensions=dims)  # X on third qudit
            P4 = PauliString.from_exponents(x_exp=[0, 0, 0], z_exp=[1, 0, 1],
                                            dimensions=dims)  # Z on first & third qudit
            S = PauliSum.from_pauli_strings([P1, P2, P3, P4])

            SPM = S.symplectic_product_matrix()
            L = S.lcm

            # Symmetric and zero diagonal modulo L
            assert np.array_equal(SPM % L, SPM.T % L)
            assert np.all((np.diag(SPM) % L) == 0)

            # Spot-check against pairwise method:

            ps1 = S.select_pauli_string(1)
            ps3 = S.select_pauli_string(3)
            assert SPM[1, 3] % L == ps1.symplectic_product(ps3, as_scalar=True) % L

    def test_pauli_sum_commutation_with_matrix(self):
        """
        Generate a random PauliSum with mixed dimensions and verify that
        pairwise commutation (via matrices) matches the symplectic scalar product.
        """
        for _ in range(int(np.ceil(np.sqrt(N_tests)))):
            dimensions = choose_random_dimensions(max_product=30)
            n_paulis = 5

            # generate random PauliSum
            P = PauliSum.from_random(n_paulis, dimensions, rand_phases=True)
            L = P.lcm
            # check commutation relations pairwise
            for i in range(n_paulis):
                for j in range(i + 1, n_paulis):
                    psi, psj = P[[i, j]]
                    # scalar symplectic product
                    s = psi.symplectic_product(psj, as_scalar=True)

                    # matrix commutator
                    Mi = psi.to_hilbert_space().toarray()
                    Mj = psj.to_hilbert_space().toarray()
                    comm = Mi @ Mj - Mj @ Mi
                    is_commuting = np.allclose(comm, np.zeros_like(comm), atol=1e-12)
                    # they commute iff scalar symplectic product is 0 (mod L)
                    assert is_commuting == (s == 0), (
                        f"Mismatch commutation for pair ({i},{j}): scalar={s} mod {L}, "
                        f"is_commuting={is_commuting}"
                        f"PauliStrings:\n{psi}\n{psj}"
                    )

                    # define product as a PauliString then get its matrix
                    prod_ps = psi * psj
                    M_prod_from_ps = prod_ps.to_hilbert_space().toarray()

                    # directly from matrices
                    M_prod_direct = Mi @ Mj

                    # check for phases being handled appropriately
                    assert np.allclose(M_prod_from_ps, M_prod_direct, atol=1e-12), (
                        f"Product matrix mismatch for pair ({i},{j}): via PauliString vs direct multiplication\n"
                        f"element i = {str(psi)}\n "
                        f"element j = {str(psj)}\n"
                        f"product i*j = {str(prod_ps)}\n"
                    )

    def test_pauli_sum_phase_invariance_matrix_equivalence(self):
        """
        Create a random PauliSum, then modify both the integer phases by adding
        multiples of 2*L and the complex weights by integer multiples of 2π
        (which are unity) so the resulting PauliSum should represent the same
        matrix. Verify matrices are equal.
        """
        for _ in range(int(np.ceil(np.sqrt(N_tests)))):
            dimensions = choose_random_dimensions(max_product=30)
            n_paulis = 5

            P = PauliSum.from_random(n_paulis, dimensions, rand_phases=False)
            n_terms = P.n_paulis()

            L = P.lcm
            # random integer shifts: add up to 2*L - 1
            phase_shifts = np.random.randint(0, 2 * L - 1, size=n_terms)
            new_phases = (P.phases + phase_shifts).tolist()

            # fix coefficients accordingly
            new_weights = (P.weights * np.exp(-2j * np.pi * phase_shifts / (2 * L))).tolist()

            # construct new PauliSum
            P_dephased = PauliSum.from_tableau(P.tableau, weights=new_weights,
                                               phases=new_phases, dimensions=P.dimensions)

            # remove phases with "phase_to_weight"
            P_rephased = P_dephased.copy()
            P_rephased.phase_to_weight()

            # check that the weights of P_rephased match those of P
            assert np.allclose(P_rephased.weights, P.weights, atol=1e-12), "Weights differ after phase_to_weight"

            # check that all PauliSums are the same:
            assert P.is_close(P_dephased, literal=False), "PauliSums differ after applying trivial phase/weight shifts"
            assert P.is_close(P_rephased, literal=False), "PauliSums differ after phase_to_weight"

            # compare matrices
            M1 = P.to_hilbert_space().toarray()
            M2 = P_dephased.to_hilbert_space().toarray()

            assert np.allclose(M1, M2, atol=1e-12), "Matrices differ after applying trivial phase/weight shifts"

    def test_is_close(self):
        # literal = True
        dims = [2, 3]
        ps1 = PauliString.from_exponents(x_exp=[1, 0], z_exp=[0, 0], dimensions=dims)  # X on qubit
        ps2 = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 1], dimensions=dims)  # XZ on qutrit
        ps3 = PauliString.from_exponents(x_exp=[1, 1], z_exp=[0, 1], dimensions=dims)  # X on qubit, XZ on qutrit

        psum1 = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0, 0.5], phases=[0, 1])
        psum2 = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0 + 1e-10, 0.5 - 1e-10], phases=[0, 1])
        psum3 = PauliSum.from_pauli_strings([ps1, ps3], weights=[1.0, 0.5], phases=[0, 1])

        assert psum1.is_close(psum2, threshold=9)
        assert not psum1.is_close(psum3, threshold=9)

        # literal = False
        psum1 = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0, 0.5], phases=[0, 1])
        psum1.phase_to_weight()
        psum2 = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0, 0.5], phases=[0, 1])
        assert not psum1.is_close(psum2, literal=True)
        assert psum1.is_close(psum2, literal=False)

        psum1 = PauliSum.from_pauli_strings([ps2, ps1], weights=[0.5, 1.0], phases=[1, 0])
        psum1.phase_to_weight()
        psum2 = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0 + 1e-12, 0.5 - 1e-12], phases=[0, 1])
        assert psum1.is_close(psum2, literal=False)

        psum1 = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0, 0.5], phases=[0, 2])
        psum1.phase_to_weight()
        psum2 = PauliSum.from_pauli_strings([ps1, ps2], weights=[1.0 + 1e-12, 0.5 - 1e-12], phases=[0, 1])
        assert not psum1.is_close(psum2, literal=False)

    def test_pauli_object_invalid_setters(self):
        p = Pauli.Xnd(1, 2)
        with pytest.raises(Exception):
            p.lcm = 2
        with pytest.raises(Exception):
            p.dimensions = np.array([2], dtype=int)

        p = PauliString.from_random([2, 3, 5])
        with pytest.raises(Exception):
            p.lcm = 12
        with pytest.raises(Exception):
            p.dimensions = np.array([2], dtype=int)

        p = PauliSum.from_random(4, [2, 3, 5])

        # Cannot set lcm
        with pytest.raises(Exception):
            p.lcm = 24

        # Cannot set dimensions
        with pytest.raises(Exception):
            p.dimensions = np.array([2], dtype=int)

        # Dimensions is read-only
        with pytest.raises(Exception):
            p.dimensions[0] = 2

        # Invalid phases length
        with pytest.raises(ValueError):
            p.phases[:] = np.array([2.4, 1, 0])

        # Cannot convert complex to int
        with pytest.raises(TypeError):
            p.phases[0] = 2.4j

        # Invalid weights length
        with pytest.raises(ValueError):
            p.weights[:] = np.array([2.4, 1, 0])

    def test_pauli_object_sum(self):
        dimension = 4
        pauli_objects = [
            Pauli.Xnd(1, dimension),
            PauliString.from_string('x2z3', dimension),
            PauliSum.from_random(3, dimensions=dimension)
        ]
        P = sum(pauli_objects, start=PauliSum.from_random(1, dimensions=dimension))
        assert isinstance(P, PauliSum)
        assert P.shape() == (6, 1)

        with pytest.raises(ValueError):
            _ = P + PauliString.from_exponents([2, 3], [0, 0], dimensions=dimension)

        ps1 = PauliString.from_exponents([0, 0], [1, 1])
        ps2 = PauliString.from_exponents([1, 0], [1, 1])
        P = ps1 + ps2
        assert P.shape() == (2, 2)
        assert P.x_exp.tolist() == [[0, 0], [1, 0]]

    def test_pauli_object_sub(self):
        ps1 = PauliString.from_exponents([0, 0, 1], [1, 1, 1])
        ps2 = PauliString.from_exponents([1, 0, 0], [1, 1, 0])
        P = ps1 - ps2
        assert P.shape() == (2, 3)
        assert P.x_exp.tolist() == [[0, 0, 1], [1, 0, 0]]
        assert P.weights.tolist() == [1, -1]

        with pytest.raises(ValueError):
            _ = P - PauliString.from_exponents([2, 3], [0, 0], dimensions=[4, 5])

    def test_pauli_ordering(self):
        p1 = Pauli.from_string("x0z1")
        p2 = Pauli.from_string("x1z0")
        assert p1 > p2

        ps1 = PauliString.from_string("x1z0 x0z1")
        ps2 = PauliString.from_string("x0z1 x1z0")
        assert ps1 < ps2
        ps3 = ps2.copy()
        assert not ps3 > ps2 and not ps3 < ps2

        psum1 = PauliSum.from_string(["x1z0 x0z1"])
        psum2 = PauliSum.from_string(["x0z1 x1z0"])
        assert psum1 < psum2

        psum1 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])
        psum2 = PauliSum.from_string(["x1z0 x0z1", "x1z0 x1z2"], [2, 3])
        with pytest.raises(Exception):
            assert psum1 < psum2
        with pytest.raises(Exception):
            assert psum1 > psum2

    def test_pauli_ordering_dimensional_mismatch(self):
        ps1 = PauliString.from_string("x1z0 x0z1")
        ps2 = PauliString.from_string("x0z1 x1z0", dimensions=[2, 3])
        with pytest.raises(Exception):
            assert ps1 < ps2

    def test_pauli_phase_setters(self):
        p1 = Pauli.from_string("x0z1")
        p2 = Pauli.from_string("x0z1")
        p3 = Pauli.from_string("x0z1")
        p4 = Pauli.from_string("x0z1")

        p1.phases[0] = 1.9
        p2.set_phases([1])
        p4.phases = [1]
        assert p1 == p2
        assert p1 == p4

        p1.reset_phases()
        assert p1 == p3
        assert p2 != p3

        ps1 = PauliString.from_string("x1z0 x0z1")
        ps2 = PauliString.from_string("x1z0 x0z1")
        ps3 = PauliString.from_string("x1z0 x0z1")
        ps4 = PauliString.from_string("x1z0 x0z1")

        ps1.phases[0] = 3
        ps2.set_phases([3])
        ps4.phases = np.asarray([3.4])
        assert ps1 == ps2
        assert ps1 == ps4

        ps1.reset_phases()
        assert ps1 == ps3
        assert ps2 != ps3

        psum1 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])
        psum2 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])
        psum3 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])
        psum4 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])

        psum1.phases[0:2] = (1, 4)
        psum2.set_phases([1, 4])
        psum4.phases = np.asarray([1.9, 4.15])
        assert psum1 == psum2
        assert psum1 == psum4

        psum1.reset_phases()
        assert psum1 == psum3
        assert psum2 != psum3

    def test_pauli_weight_setters(self):
        p1 = Pauli.from_string("x0z1")
        p2 = Pauli.from_string("x0z1")
        p3 = Pauli.from_string("x0z1")
        p4 = Pauli.from_string("x0z1")

        p1.weights[0] = 1.9 + 2j
        p2.set_weights([1.9 + 2j])
        p4.weights = [1.9 + 2j]
        assert p1 == p2
        assert p1 == p4

        p1.reset_weights()
        assert p1 == p3
        assert p2 != p3

        ps1 = PauliString.from_string("x1z0 x0z1")
        ps2 = PauliString.from_string("x1z0 x0z1")
        ps3 = PauliString.from_string("x1z0 x0z1")
        ps4 = PauliString.from_string("x1z0 x0z1")

        ps1.weights[0] = -1
        ps2.set_weights([-1])
        ps4.weights = np.asarray([-1])
        assert ps1 == ps2
        assert ps1 == ps4

        ps1.reset_weights()
        assert ps1 == ps3
        assert ps2 != ps3

        psum1 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])
        psum2 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])
        psum3 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])
        psum4 = PauliSum.from_string(["x1z2 x2z1", "x0z0 x1z2"], [2, 3])

        psum1.weights[0:2] = (1.9, 4.15)
        psum2.set_weights([1.9, 4.15])
        psum4.weights = np.asarray([1.9, 4.15])
        assert psum1 == psum2
        assert psum1 == psum4

        psum1.reset_weights()
        assert psum1 == psum3
        assert psum2 != psum3

    def test_make_hermitian(self):
        available_dimensions = [2, 3, 5, 7]
        # with random phases
        for i in range(50):
            n_qubits = np.random.randint(1, 10)
            dims = [np.random.choice(available_dimensions) for _ in range(n_qubits)]
            n_paulis = np.random.randint(1, np.min([4**n_qubits - 1, 10]))
            P = PauliSum.from_random(n_paulis=n_paulis,
                                     dimensions=dims,
                                     rand_weights=True)
            P.set_phases(np.random.randint(0, 2 * P.lcm, n_paulis))
            P1 = P.make_hermitian()
            assert P1.is_hermitian()

        # without random phases
        for i in range(50):
            n_qubits = np.random.randint(1, 10)
            dims = [np.random.choice(available_dimensions) for _ in range(n_qubits)]
            n_paulis = np.random.randint(1, np.min([4**n_qubits - 1, 10]))
            P = PauliSum.from_random(n_paulis=n_paulis,
                                     dimensions=dims,
                                     rand_weights=True)
            P1 = P.make_hermitian()
            assert P1.is_hermitian()

    def test_ordered_eigenspectrum(self):
        for _ in range(N_tests):
            dimensions = [2, 3, 5, 7]
            n_paulis = len(dimensions)

            P = PauliSum.from_random(n_paulis, dimensions, rand_weights=False).make_hermitian()
            assert P.is_hermitian()

            m = np.around(P.to_hilbert_space().toarray(), 10)
            energies, states = P.ordered_eigenspectrum()

            assert len(energies) == len(states)
            assert len(states) == np.prod(dimensions)

            # Check: the eigenvectors give raise to the correct eigenvalues
            for energy, state in zip(energies, states):
                check_energy = np.around(state.conjugate().transpose() @ m @ state, 10)
                assert np.isclose(
                    check_energy, energy), f"eigenvalue mismatch for state {state}: {energy} vs {check_energy}."

    def test_ordered_eigenspectrum_raise_non_hermitian(self):
        for _ in range(N_tests):
            dimensions = [2, 3, 5, 7]
            n_paulis = len(dimensions)

            P = PauliSum.from_random(n_paulis, dimensions)
            if P.is_hermitian():
                continue

            with pytest.raises(ValueError):
                _, _ = P.ordered_eigenspectrum()

    @pytest.mark.skip()
    def test_stabilizer_to_hilbert_space(self):
        for _ in range(N_tests):
            dimensions = [2] * 5
            n_qubits = len(dimensions)

            # P = PauliSum.from_random(n_paulis, dimensions, rand_weights=False).make_hermitian()
            # P.reset_weights()
            # assert P.is_hermitian()

            # Test both less and equal number of paulis than qudits
            n_paulis = np.random.randint(1, n_qubits)

            # Create a stabilizer IIIZ, IIZI, IZII, leaving identities at the beginning if n_paulis < n_qubits
            tableau = np.zeros((n_paulis, 2 * n_qubits), dtype=int)
            tableau[:, n_qubits:] = np.eye(n_paulis, n_qubits, dtype=int)

            # Create stabilizer in PauliSum form
            stabilizer = PauliSum.from_tableau(tableau,
                                               weights=np.ones(n_paulis),
                                               dimensions=dimensions)

            # Initialize stabilizer to random computational state
            phases = [2 * np.random.randint(0, stabilizer.lcm - 1) for _ in range(n_paulis)]
            stabilizer.set_phases(phases)

            # Random Clifford circuit
            n_gates = 10 * n_paulis**2
            C = Circuit.from_random(n_gates, dimensions)

            # Act with Clifford on stabilizer
            stabilizer_shuffled = C.act(stabilizer)
            # Get hilbert space representations
            stabilizer_shuffled_hilbert = stabilizer_shuffled.stabilizer_to_hilbert_space()

            # Ensure the state stabilizes all PauliStrings in the shuffled stabilizer
            for idx in range(stabilizer_shuffled.n_paulis()):

                phases = stabilizer_shuffled.phases.copy()
                phases[idx] = 0
                phase_to_test: int = stabilizer_shuffled.phases[idx]

                ps = stabilizer_shuffled[[idx]]
                ps_hilbert = ps.to_hilbert_space()

                lhs = np.exp(
                    (1j * 2 * np.pi * phase_to_test) / (2 * stabilizer.lcm)
                )

                rhs = (
                    stabilizer_shuffled_hilbert @ ps_hilbert @ stabilizer_shuffled_hilbert.T
                ).trace()

                print(np.around(lhs - rhs, 10))

                assert np.around(lhs - rhs, 10) == 0
