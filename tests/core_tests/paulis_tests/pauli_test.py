import numpy as np
import pytest
from sympleq.core.paulis import PauliSum, PauliString, Pauli
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION


class TestPaulis:

    def random_pauli_string(self, dim: int) -> tuple[PauliString, int, int, int, int]:
        r1 = np.random.randint(0, dim)
        r2 = np.random.randint(0, dim)
        s1 = np.random.randint(0, dim)
        s2 = np.random.randint(0, dim)
        return PauliString.from_string(f"x{r1}z{s1} x{r2}z{s2}", dimensions=[dim, dim]), r1, r2, s1, s2

    def test_pauli_multiplication(self):
        for dim in [2]:
            x1 = Pauli.Xnd(1, dim)
            y1 = Pauli.Ynd(1, dim)
            z1 = Pauli.Znd(1, dim)
            id = Pauli.Idnd(dim)

            assert x1 * z1 == y1, 'Error in Pauli multiplication (x * z = y) ' + (x1 * z1).__str__()
            assert x1**dim == id, 'Error in Pauli exponentiation (x**dim = id) ' + (x1**dim).__str__()
            assert y1**dim == id, 'Error in Pauli exponentiation (y**dim = id) ' + (y1**dim).__str__()
            assert z1**dim == id, 'Error in Pauli exponentiation (z**dim = id)  ' + (z1**dim).__str__()
            assert x1 * y1 == z1, 'Error in Pauli multiplication (x * y = z) ' + (x1 * y1).__str__()
            assert y1 * z1 == x1, 'Error in Pauli multiplication (y * z = x) ' + (y1 * z1).__str__()
            assert z1 * x1 == y1, 'Error in Pauli multiplication (z * x = y) ' + (z1 * x1).__str__()
            assert x1 * id == x1, 'Error in Pauli multiplication (x * id = x) ' + (x1 * id).__str__()
            assert y1 * id == y1, 'Error in Pauli multiplication (y * id = y) ' + (y1 * id).__str__()
            assert z1 * id == z1, 'Error in Pauli multiplication (z * id = z) ' + (z1 * id).__str__()

        for dim in [3, 5, 7, 11, 13, 17]:
            for i in range(100):
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
        for dim in [2, 3, 5, 7, 11, 13, 15, 17]:
            for i in range(100):
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
                output_str_correct = f"x{(r11 + r21) % dim}z{(s11 + s21) % dim} x{(r12 + r22) % dim}z{(s12 + s22) % dim}"
                input_ps1 = PauliString.from_string(input_str1, dimensions=[dim, dim])
                input_ps2 = PauliString.from_string(input_str2, dimensions=[dim, dim])
                output_ps = input_ps1 * input_ps2
                assert output_ps == PauliString.from_string(output_str_correct, dimensions=[
                                                            dim, dim]), 'Error in PauliString multiplication'

    def test_pauli_string_tensor_product(self):
        for dimensions in [2, 3, 5, 7, 11]:
            for i in range(100):
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
        for dim in [2, 3, 5, 7]:
            for _ in range(100):
                p_string1, r1, r2, s1, s2 = self.random_pauli_string(dim)
                ps0 = Pauli.from_string(f"x{r1}z{s1}", dimension=dim)
                ps1 = Pauli.from_string(f"x{r2}z{s2}", dimension=dim)

                assert p_string1[0] == ps0, 'Error in PauliString indexing (first PauliString)'
                assert p_string1[1] == ps1, 'Error in PauliString indexing'

    def test_pauli_sum_multiplication(self):
        for dim in [2, 3, 5]:

            for i in range(100):
                p_string1, r1, r2, s1, s2 = self.random_pauli_string(dim)
                p_string2, r12, r22, s12, s22 = self.random_pauli_string(dim)

                random_pauli_sum = PauliSum.from_pauli_strings([p_string1, p_string2])

                # Test multiplication of PauliSum with PauliString
                input_ps1, r13, r23, s13, s23 = self.random_pauli_string(dim)

                output_str_correct = (
                    f"x{(r1 + r13) % dim}z{(s1 + s13) % dim} "
                    f"x{(r2 + r23) % dim}z{(s2 + s23) % (2 * dim)}"
                )
                output_str2_correct = (
                    f"x{(r12 + r13) % dim}z{(s12 + s13) % dim} "
                    f"x{(r22 + r23) % dim}z{(s22 + s23) % (2 * dim)}"
                )
                output_ps = random_pauli_sum * input_ps1

                phase1 = 2 * (s1 * r13 + s2 * r23) % (2 * dim)
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

        for dim in [2, 3, 5, 17]:
            for _ in range(50):
                p_string1, r1, r2, s1, s2 = self.random_pauli_string(dim)
                p_string2, r12, r22, s12, s22 = self.random_pauli_string(dim)
                p_string3, r13, r23, s13, s23 = self.random_pauli_string(dim)

                random_pauli_sum = PauliSum.from_pauli_strings([p_string1, p_string2])
                random_pauli_sum2 = PauliSum.from_pauli_strings([p_string3])
                ps_out = random_pauli_sum @ random_pauli_sum2
                ps_out_correct = PauliSum.from_pauli_strings(
                    [p_string1 @ p_string3, p_string2 @ p_string3])
                assert ps_out == ps_out_correct, 'Error in PauliSum tensor product'

    def test_symplectic_matrix_single_pauli(self):
        pauli_list = ['x1z0']
        weights = np.array([1.])
        sp = PauliSum.from_string(pauli_list, dimensions=[2], weights=weights)
        expected_matrix = np.array([[1, 0]])

        np.testing.assert_array_equal(sp.tableau, expected_matrix)

    def test_symplectic_matrix_multiple_paulis(self):
        pauli_list = ['x1z0', 'x0z1', 'x1z1']
        weights = np.array([1, 2, 3])
        sp = PauliSum.from_string(pauli_list, dimensions=[2], weights=weights)
        expected_matrix = np.array([
            [1, 0],
            [0, 1],
            [1, 1]
        ])

        np.testing.assert_array_equal(sp.tableau, expected_matrix)

    def test_basic_pauli_relations(self):
        dims = 3
        x1 = Pauli.from_string('x1z0', dimension=dims)
        y1 = Pauli.from_string('x1z1', dimension=dims)
        z1 = Pauli.from_string('x0z1', dimension=dims)
        id = Pauli.from_string('x0z0', dimension=dims)

        assert x1 * z1 == y1
        assert x1 * x1 * x1 == id

    def test_pauli_string_construction(self):
        dims = [3, 3]
        x1x1 = PauliString.from_string('x1z0 x1z0', dimensions=dims)
        x1x1_2 = PauliString.from_exponents([1, 1], [0, 0], dimensions=dims)

        assert x1x1 == x1x1_2

        x1y1 = PauliString.from_exponents([1, 1], [0, 1], dimensions=dims)
        x1y1_2 = PauliString.from_string('x1z0 x1z1', dimensions=dims)

        assert x1y1 == x1y1_2

        x1y0 = PauliString.from_exponents([1, 1], [1, 0])
        x1y0_2 = PauliString.from_string('x1z1 x1z0', dimensions=DEFAULT_QUDIT_DIMENSION)
        x1y0_3 = PauliString.from_string('x1z1 x1z0', dimensions=dims)

        assert x1y0 == x1y0_2
        assert x1y0 != x1y0_3

        # Test minimum allowed qudit dimension.
        with pytest.raises(ValueError):
            _ = PauliString.from_string('x0z0 x0z0', dimensions=DEFAULT_QUDIT_DIMENSION - 1)

    def test_pauli_sum_addition(self):

        for dimensions in [2, 3, 5]:
            for _ in range(100):
                p_string1, r1, r2, s1, s2 = self.random_pauli_string(dimensions)
                p_string2, r12, r22, s12, s22 = self.random_pauli_string(dimensions)
                p_string3, r13, r23, s13, s23 = self.random_pauli_string(dimensions)

                random_pauli_sum = PauliSum.from_pauli_strings([p_string1, p_string2])
                random_pauli_sum2 = PauliSum.from_pauli_strings([p_string3])
                ps_out = random_pauli_sum + random_pauli_sum2
                ps_out_correct = PauliSum.from_pauli_strings([p_string1, p_string2, p_string3])
                assert ps_out == ps_out_correct, 'Error in PauliSum addition'

        dimensions = [3, 3]
        x1x1 = PauliSum.from_pauli_strings(PauliString.from_string('x1z0 x1z0', dimensions))
        x1y1 = PauliSum.from_pauli_strings(PauliString.from_string('x1z0 x1z1', dimensions))

        psum = x1x1 + x1y1
        psum.standardise()
        expected = PauliSum.from_pauli_strings(
            [PauliString.from_string('x1z0 x1z0', dimensions), PauliString.from_string('x1z0 x1z1', dimensions)],
            weights=[1, 1], phases=[0, 0])

        assert psum == expected

    def test_phase_and_dot_product(self):
        d = 7
        x = PauliString.from_string('x1z0', dimensions=[d])
        z = PauliString.from_string('x0z1', dimensions=[d])

        assert z.acquired_phase(x) == 2.0, 'Expected phase to be 2.0, got {}'.format(z.acquired_phase(x))
        assert x.acquired_phase(z) == 0.0, 'Expected phase to be 0.0, got {}'.format(x.acquired_phase(z))

        dims = [3, 3]
        x1x1 = PauliSum.from_pauli_strings(PauliString.from_string('x1z0 x1z0', dimensions=dims))
        x1y1 = PauliSum.from_pauli_strings(PauliString.from_string('x1z0 x1z1', dimensions=dims))

        s1 = x1x1 + x1y1 * 0.5
        s2 = x1x1 + x1x1
        print("1", s1.tableau)
        print("2", s2.tableau)

        s3 = PauliSum.from_string(['x2z0 x2z0', 'x2z0 x2z0', 'x2z0 x2z1', 'x2z0 x2z1'],
                                  weights=[1, 1, 0.5, 0.5],
                                  phases=[0, 0, 2, 2],
                                  dimensions=dims)

        print("3", s3.tableau)

        assert s1 * s2 == s3, 'Expected s1 * s2 to equal s3, got {}'.format(s1 * s2) + '\n' + s3.__str__()

    def test_tensor_product_distributivity(self):
        dimensions = [3, 3]
        x1x1 = PauliSum.from_string('x1z0 x1z0', dimensions)
        x1y1 = PauliSum.from_string('x1z0 x1z1', dimensions)

        s1 = x1x1 + x1y1 * 0.5
        s2 = x1x1 + x1x1

        left = (s1 + s2) @ s2
        right = s1 @ s2 + s2 @ s2

        assert left == right

    def test_pauli_sum_indexing(self):
        dims = [3, 3, 3]
        ps = PauliSum.from_string(['x2z0 x2z0 x1z1', 'x2z0 x2z0 x0z0', 'x2z0 x2z1 x2z0', 'x2z0 x2z1 x1z1'],
                                  weights=[1, 1, 0.5, 0.5],
                                  phases=[0, 0, 1, 1],
                                  dimensions=dims)

        ps0 = PauliString.from_string('x2z0 x2z0 x1z1', dimensions=dims)
        ps1 = PauliString.from_string('x2z0 x2z0 x0z0', dimensions=dims)
        ps2 = PauliString.from_string('x2z0 x2z1 x2z0', dimensions=dims)
        ps3 = PauliString.from_string('x2z0 x2z1 x1z1', dimensions=dims)
        assert ps.select_pauli_string(0) == ps0, f'{ps[0].__str__()}\n{ps0.__str__()}'
        assert ps.select_pauli_string(1) == ps1, f'{ps[1].__str__()}\n{ps1.__str__()}'
        assert ps.select_pauli_string(2) == ps2, f'{ps[2].__str__()}\n{ps2.__str__()}'
        assert ps.select_pauli_string(3) == ps3, f'{ps[3].__str__()}\n{ps3.__str__()}'
        assert ps[0:2] == PauliSum.from_string(['x2z0 x2z0 x1z1', 'x2z0 x2z0 x0z0'], dimensions=dims)
        assert ps[[0, 3]] == PauliSum.from_string(['x2z0 x2z0 x1z1', 'x2z0 x2z1 x1z1'], weights=[1, 0.5], phases=[0, 1],
                                                  dimensions=dims)
        assert ps[[0, 2], 1] == PauliSum.from_string(['x2z0', 'x2z1'], weights=[1, 0.5], phases=[0, 1], dimensions=[3])
        assert ps[[0, 2], [0, 2]] == PauliSum.from_string(['x2z0 x1z1', 'x2z0 x2z0'], weights=[1, 0.5], phases=[0, 1],
                                                          dimensions=[3, 3])

    def test_pauli_sum_amend(self):
        dims = [2, 3]
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
        n_qudits = 10
        d = 2
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

        N = 250
        dimensions = [2, 3, 5]
        for n_paulis in range(1, 4):
            for _ in range(N):
                P1 = PauliSum.from_random(n_paulis, dimensions)
                P2 = PauliSum.from_random(n_paulis, dimensions)
                P_res = P1 * P2
                phase_symplectic = P_res.phases[0]

                phase_computed = 0
                for j in range(3):
                    s1 = P1.z_exp[0, j]
                    r2 = P2.x_exp[0, j]
                    phase_computed += ((s1 * r2) % dimensions[j]) * P1.lcm / dimensions[j]
                phase_computed = phase_computed * 2 % (2 * P1.lcm)
                assert phase_symplectic == phase_computed

    def test_pauli_sum_delete_qudits(self):
        dims = [2, 3, 5, 6, 7]

        ps1 = PauliString.from_exponents(x_exp=[1, 2, 0, 0, 3], z_exp=[0, 1, 4, 4, 5], dimensions=dims)
        ps2 = PauliString.from_exponents(x_exp=[0, 1, 3, 4, 3], z_exp=[1, 0, 2, 4, 5], dimensions=dims)
        psum = PauliSum.from_pauli_strings([ps1, ps2])
        psum._delete_qudits([1, 3])

        expected_ps1 = PauliString.from_exponents(x_exp=[1, 0, 3], z_exp=[0, 4, 5], dimensions=[2, 5, 7])
        expected_ps2 = PauliString.from_exponents(x_exp=[0, 3, 3], z_exp=[1, 2, 5], dimensions=[2, 5, 7])
        expected_psum = PauliSum.from_pauli_strings([expected_ps1, expected_ps2])

        assert psum == expected_psum, f"Expected {expected_psum}, got {psum}"

    def test_symplectic_product(self):
        P1 = PauliString.from_string('x1z0', dimensions=[2])
        P2 = PauliString.from_string('x0z1', dimensions=[2])
        assert P1.symplectic_product(P2) == 1

        P1 = PauliString.from_string('x1z0', dimensions=[2])
        P2 = PauliString.from_string('x1z0', dimensions=[2])
        assert P1.symplectic_product(P2) == 0

        P1 = PauliString.from_string('x0z1', dimensions=[2])
        P2 = PauliString.from_string('x0z1', dimensions=[2])
        assert P1.symplectic_product(P2) == 0

        P1 = PauliString.from_string('x1z0 x1z0', dimensions=[2, 2])
        P2 = PauliString.from_string('x0z1 x0z1', dimensions=[2, 2])
        assert P1.symplectic_product(P2) == 0

        P1 = PauliString.from_string('x1z0 x0z1', dimensions=[2, 2])
        P2 = PauliString.from_string('x1z0 x1z0', dimensions=[2, 2])
        assert P1.symplectic_product(P2) == 1

        P1 = PauliString.from_string('x1z0', dimensions=[3])
        P2 = PauliString.from_string('x2z0', dimensions=[3])
        assert P1.symplectic_product(P2) == 0

        P1 = PauliString.from_string('x1z2', dimensions=[3])
        P2 = PauliString.from_string('x2z1', dimensions=[3])
        assert P1.symplectic_product(P2) == 0

        P1 = PauliString.from_string('x1z2 x1z1', dimensions=[3, 2])
        P2 = PauliString.from_string('x2z1 x1z1', dimensions=[3, 2])
        assert P1.symplectic_product(P2) == 0

    def test_hermitian_generation(self):
        P1 = PauliString.from_string('x1z0', dimensions=[3])
        P2 = PauliString.from_string('x2z0', dimensions=[3])
        assert P1.H() == P2

        P1 = PauliString.from_string('x1z1', dimensions=[3])
        P2 = PauliString.from_string('x2z2', dimensions=[3])
        assert P1.H() == P2

        P1 = PauliString.from_string('x2z1', dimensions=[3])
        P2 = PauliString.from_string('x1z2', dimensions=[3])
        assert P1.H() == P2

        P1 = PauliString.from_string('x0z0', dimensions=[3])
        P2 = PauliString.from_string('x0z0', dimensions=[3])
        assert P1.H() == P2
        assert P1.is_hermitian()

        P1 = PauliString.from_string('x0z1 x1z1', dimensions=[3, 2])
        P2 = PauliString.from_string('x0z2 x1z1', dimensions=[3, 2])
        P3 = PauliString.from_exponents([0, 1], [2, 1], dimensions=[3, 2])
        assert P1.H() == P3

        P1 = PauliString.from_string('x0z1 x1z0', dimensions=[5, 2])
        P2 = PauliString.from_string('x0z4 x1z0', dimensions=[5, 2])
        assert P1.H() == P2

        P1 = PauliString.from_string('x2z1 x0z1', dimensions=[5, 2])
        P2 = PauliString.from_string('x3z4 x0z1', dimensions=[5, 2])
        assert P1.H() == P2

    def test_pauli_sum_is_hermitian(self):
        for run in range(40):
            # Generate random dimensions array {d_1, d_2, d_3,...} such that \prod_i d_i <= 10
            dimensions = []
            while True:
                import random
                d = random.sample([2, 3, 5], 1)[0]
                if d * int(np.prod(dimensions)) > 10:
                    break
                dimensions.append(d)

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
            tolerance = 10**(-12)
            assert np.max(np.abs((pauli_sum.to_hilbert_space().toarray() - H_e))) < tolerance
            assert pauli_sum.is_hermitian() == np.array_equal(H_e, H_e.conjugate().transpose())

    def test_qubit_XZ_phase_is_minus_one(self):
        # Single qubit (dimension 2): X * Z = (-1) Z * X  => scalar exponent r = 1 mod 2
        dims = [2]
        psX = PauliString.from_exponents(x_exp=[1], z_exp=[0], dimensions=dims)  # X
        psZ = PauliString.from_exponents(x_exp=[0], z_exp=[1], dimensions=dims)  # Z

        # per-qudit residue
        rj = psX.symplectic_residues(psZ)
        assert rj.tolist() == [1]

        # scalar phase mod LCM
        L = int(np.lcm.reduce(np.array(dims, int)))
        r_scalar = psX.symplectic_product(psZ, as_scalar=True)
        assert r_scalar == 1 % L

    def test_mixed_dims_qutrit_XZ_phase(self):
        # dims = [2, 3]; use site 1 (qutrit) to get omega_3
        dims = [2, 3]
        P = PauliString.from_exponents(x_exp=[0, 1], z_exp=[0, 0], dimensions=dims)  # X on qutrit
        Q = PauliString.from_exponents(x_exp=[0, 0], z_exp=[0, 1], dimensions=dims)  # Z on qutrit

        # residues: [0 mod 2, 1 mod 3]
        assert P.symplectic_residues(Q).tolist() == [0, 1]

        # scalar: L = 6, weights = [3, 2], r = 0*3 + 1*2 = 2 mod 6
        assert P.symplectic_product(Q, as_scalar=True) == 2

    def test_mixed_dims_all_sites_X_vs_Z_product(self):
        dims = [2, 3, 5]
        P = PauliString.from_exponents(x_exp=[1, 1, 1], z_exp=[0, 0, 0], dimensions=dims)  # X on all
        Q = PauliString.from_exponents(x_exp=[0, 0, 0], z_exp=[1, 1, 1], dimensions=dims)  # Z on all

        # residues per site are all 1
        assert P.symplectic_residues(Q).tolist() == [1, 1, 1]

        # scalar: L=30, weights=[15,10,6], r = 15+10+6 = 31 ≡ 1 (mod 30)
        assert P.symplectic_product(Q, as_scalar=True) == 1

    def test_bilinearity_scalar_mode(self):
        # Check <v1+v2, w> = <v1, w> + <v2, w>  (phase-preserving scalar)
        dims = [2, 3]
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
        dims = [2, 5]
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
        dims = [2, 3]
        # P1 = X on qubit, P2 = Z on qutrit, P3 = XZ on qutrit
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
        # Symmetry and zero diagonal
        dims = [2, 3, 5]
        P1 = PauliString.from_exponents(x_exp=[1, 0, 0], z_exp=[0, 0, 0], dimensions=dims)  # X on qubit
        P2 = PauliString.from_exponents(x_exp=[0, 1, 0], z_exp=[0, 1, 0], dimensions=dims)  # XZ on qutrit
        P3 = PauliString.from_exponents(x_exp=[0, 0, 1], z_exp=[0, 0, 0], dimensions=dims)  # X on ququint
        P4 = PauliString.from_exponents(x_exp=[0, 0, 0], z_exp=[1, 0, 1], dimensions=dims)  # Z on qubit & ququint
        S = PauliSum.from_pauli_strings([P1, P2, P3, P4])
        print("S", S.tableau.shape)

        SPM = S.symplectic_product_matrix()
        L = S.lcm

        # Symmetric and zero diagonal modulo L
        assert np.array_equal(SPM % L, SPM.T % L)
        assert np.all((np.diag(SPM) % L) == 0)

        # Spot-check against pairwise method:
        ps1 = S.select_pauli_string(1)
        ps3 = S.select_pauli_string(3)
        assert SPM[1, 3] % L == ps1.symplectic_product(ps3, as_scalar=True) % L

    # Comprehensive tests
    # commutation relations for mixed dimensions
    # phases for mixed dimensions
