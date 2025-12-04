import numpy as np
import random
from sympleq.core.circuits import SUM, SWAP, Hadamard, PHASE, Gate, Circuit, PauliGate
from sympleq.core.circuits.utils import is_symplectic
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits.random_symplectic import (symplectic_gf2, symplectic_group_size,
                                                     symplectic_random_transvection)


class TestGates():

    @staticmethod
    def states(dim):
        ps1 = PauliString.from_string('x1z0 x0z0', dimensions=[dim, dim])
        ps2 = PauliString.from_string('x0z0 x0z1', dimensions=[dim, dim])
        ps3 = PauliString.from_string('x1z0 x1z0', dimensions=[dim, dim])
        ps4 = PauliString.from_string('x0z1 x0z1', dimensions=[dim, dim])
        ps5 = PauliString.from_string('x1z1 x0z0', dimensions=[dim, dim])
        p_sum = PauliSum.from_string(['x1z0 x0z0 x1z1', 'x0z0 x0z1 x1z0',
                                      'x1z1 x1z0 x0z0'], dimensions=[dim, dim, dim])
        return ps1, ps2, ps3, ps4, ps5, p_sum

    def random_pauli_sum(self, dim, n_paulis=10):
        # Generates a random PauliSum with n_paulis random PauliStrings of dimension dim
        # TODO: remove this and use from_random, make sure that the same PauliStrings are not generated in from_random
        ps_list = []
        element_list = [(0, 0, 0, 0)]  # to keep track of already generated PauliStrings. Avoids identity and duplicates
        for _ in range(n_paulis):
            ps, r1, r2, s1, s2 = self.random_pauli_string(dim)
            element_list.append((r1, r2, s1, s2))
            while (r1, r2, s1, s2) in element_list:
                ps, r1, r2, s1, s2 = self.random_pauli_string(dim)
            ps_list.append(ps)
        return PauliSum.from_pauli_strings(ps_list)

    def random_pauli_string(self, dim):
        # Generates a random PauliString of dimension dim
        r1 = np.random.randint(0, dim)
        r2 = np.random.randint(0, dim)
        s1 = np.random.randint(0, dim)
        s2 = np.random.randint(0, dim)
        return PauliString.from_string(f'x{r1}z{s1} x{r2}z{s2}', dimensions=[dim, dim]), r1, r2, s1, s2

    def test_SUM(self):
        # acts the SUM gate on a bunch of random pauli strings and pauli sums
        for d in [2, 5, 11]:
            # test pauli_strings
            for i in range(1000):
                r1 = np.random.randint(0, d)
                r2 = np.random.randint(0, d)
                s1 = np.random.randint(0, d)
                s2 = np.random.randint(0, d)

                input_str = f"x{r1}z{s1} x{r2}z{s2}"
                output_str_correct = f"x{r1}z{(s1 - s2) % d} x{(r2 + r1) % d}z{s2}"

                input_ps = PauliString.from_string(input_str, dimensions=[d, d])
                output_ps = SUM(0, 1, d).act(input_ps)
                assert output_ps == PauliString.from_string(output_str_correct, dimensions=[d, d]), 'Error in SUM gate'

            # test pauli_sums
            for i in range(10):
                ps_list_in = []
                ps_list_out_correct = []
                ps_phase_out_correct = []
                for j in range(10):
                    r1 = np.random.randint(0, d)
                    r2 = np.random.randint(0, d)
                    s1 = np.random.randint(0, d)
                    s2 = np.random.randint(0, d)

                    input_str = f"x{r1}z{s1} x{r2}z{s2}"
                    output_str_correct = f"x{r1}z{(s1 - s2) % d} x{(r2 + r1) % d}z{s2}"

                    input_ps = PauliString.from_string(input_str, dimensions=[d, d])
                    output_ps_correct = PauliString.from_string(output_str_correct, dimensions=[d, d])
                    ps_list_in.append(input_ps)
                    ps_list_out_correct.append(output_ps_correct)
                    ps_phase_out_correct.append(0)

                input_psum = PauliSum.from_pauli_strings(ps_list_in)
                output_psum_correct = PauliSum.from_pauli_strings(ps_list_out_correct, phases=ps_phase_out_correct)

                output_psum = SUM(0, 1, d).act(input_psum)
                assert output_psum == output_psum_correct, (
                    'Error in SUM gate: \n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_SWAP(self):
        # acts the SUM gate on a bunch of random pauli strings and pauli sums
        # C = np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
        # U = np.zeros((4, 4), dtype=int)
        # U[2:, :2] = np.eye(2, dtype=int)
        for d in [2, 5, 11]:
            # test pauli_strings
            for i in range(100):
                r1 = np.random.randint(0, d)
                r2 = np.random.randint(0, d)
                s1 = np.random.randint(0, d)
                s2 = np.random.randint(0, d)

                input_str = f"x{r1}z{s1} x{r2}z{s2}"
                output_str_correct = f"x{r2}z{(s2) % d} x{r1}z{s1}"

                input_ps = PauliString.from_string(input_str, dimensions=[d, d])
                output_ps = SWAP(0, 1, d).act(input_ps)
                assert output_ps == PauliString.from_string(output_str_correct, dimensions=[d, d]), 'Error in SUM gate'

            # test pauli_sums
            for i in range(100):
                ps_list_in = []
                ps_list_out_correct = []
                ps_phase_out_correct = []
                for j in range(10):
                    r1 = np.random.randint(0, d)
                    r2 = np.random.randint(0, d)
                    s1 = np.random.randint(0, d)
                    s2 = np.random.randint(0, d)

                    input_str = f"x{r1}z{s1} x{r2}z{s2}"
                    output_str_correct = f"x{r2}z{s2} x{r1}z{s1}"

                    input_ps = PauliString.from_string(input_str, dimensions=[d, d])
                    output_ps_correct = PauliString.from_string(output_str_correct, dimensions=[d, d])
                    ps_list_in.append(input_ps)
                    ps_list_out_correct.append(output_ps_correct)
                    # a = np.array([r1, r2, s1, s2])

                    # phase_out = np.diag(C.T @ U @ C) @ a
                    # + 2 * a @ np.triu(C.T @ U @ C) @ a - a @ np.diag(np.diag(C.T @ U @ C)) @ a
                    ps_phase_out_correct.append(0)

                input_psum = PauliSum.from_pauli_strings(ps_list_in)
                output_psum_correct = PauliSum.from_pauli_strings(ps_list_out_correct,
                                                                  phases=ps_phase_out_correct)

                output_psum = SWAP(0, 1, d).act(input_psum)
                assert output_psum == output_psum_correct, (
                    'Error in SUM gate: \n' +
                    input_psum.__str__() + '\n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_Hadamard(self):
        # TODO: Be certain of inverse convention - ultimately arbitrary but should match prevalent literature
        for d in [2, 5, 11]:
            # test pauli_strings
            gate = Hadamard(0, d, inverse=False)  # Hadamard on qubit 0

            for i in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{(-s1) % d}z{r1} x{r2}z{s2}"

                output_ps = gate.act(input_ps)
                assert output_ps.has_equal_tableau(PauliString.from_string(
                    output_str_correct, dimensions=[d, d]
                )), 'Error in Hadamard gate 0'

            gate = Hadamard(1, d, inverse=False)  # Hadamard on qudit 1
            for i in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{r1}z{s1} x{(-s2) % d}z{r2}"

                output_ps = gate.act(input_ps)
                assert output_ps.has_equal_tableau(PauliString.from_string(
                    output_str_correct, dimensions=[d, d]
                )), 'Error in Hadamard gate 1'
            # test pauli_sums
            gate = Hadamard(0, d, inverse=False)  # Hadamard on qubit 0

            for i in range(100):
                ps_list_in = []
                ps_list_out_correct = []
                ps_phase_out_correct = []
                for j in range(10):
                    input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                    output_str_correct = f"x{(-s1) % d}z{r1} x{r2}z{s2}"

                    ps_list_in.append(input_ps)
                    ps_list_out_correct.append(PauliString.from_string(output_str_correct, dimensions=[d, d]))
                    ps_phase_out_correct.append(-2 * r1 * s1)

                input_psum = PauliSum.from_pauli_strings(ps_list_in)
                output_psum = gate.act(input_psum)
                output_psum_correct = PauliSum.from_pauli_strings(
                    ps_list_out_correct, phases=ps_phase_out_correct)
                assert output_psum.has_equal_tableau(output_psum_correct), (
                    'Error in Hadamard gate: \n' +
                    input_psum.__str__() + '\n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_PHASE(self):
        # TODO: Better approach - is there another approach that we can test the first one with like Hadamard?
        for d in [2, 5, 11]:
            gate = PHASE(0, d)  # PHASE on qubit 0

            for i in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{r1}z{(r1 + s1) % d} x{r2}z{s2}"

                output_ps = gate.act(input_ps)
                assert output_ps.has_equal_tableau(PauliString.from_string(
                    output_str_correct, dimensions=[d, d])), 'Error in H gate 0'

            gate = PHASE(1, d)  # PHASE on qudit 1
            for i in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{r1}z{s1} x{r2}z{(r2 + s2) % d}"

                output_ps = gate.act(input_ps)
                assert output_ps.has_equal_tableau(PauliString.from_string(
                    output_str_correct, dimensions=[d, d])), 'Error in H gate 1'
            # test pauli_sums
            gate = PHASE(0, d)  # PHASE on qubit 0
            h = gate.phase_vector
            C = gate.symplectic
            for i in range(100):
                ps_list_in = []
                ps_list_out_correct = []
                ps_phase_out_correct = []
                for j in range(10):
                    input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                    output_str_correct = f"x{r1}z{(r1 + s1) % d} x{r2}z{s2}"

                    ps_list_in.append(input_ps)
                    ps_list_out_correct.append(PauliString.from_string(output_str_correct, dimensions=[d, d]))
                    a1 = np.array([r1, s1])
                    U = np.zeros((input_ps.n_qudits(), input_ps.n_qudits()), dtype=int)
                    U[1, 0] = 1

                    p = (
                        np.dot(h, a1) -
                        np.diag(C.T @ U @ C) @ a1 +
                        2 * a1 @ np.triu(C.T @ U @ C) @ a1 -
                        a1 @ np.diag(np.diag(C.T @ U @ C)) @ a1
                    )
                    ps_phase_out_correct.append(p)

                input_psum = PauliSum.from_pauli_strings(ps_list_in)
                output_psum = gate.act(input_psum)
                output_psum_correct = PauliSum.from_pauli_strings(
                    ps_list_out_correct, phases=ps_phase_out_correct)
                assert output_psum.has_equal_tableau(output_psum_correct), (
                    'Error in Hadamard gate: \n' +
                    input_psum.__str__() + '\n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_group_homomorphism(self):
        # Tests all gates and all combinations of two pauli strings for the group homomorphism property:
        # gate.act(p1) * gate.act(p2) == gate.act(p1 * p2)

        gates = [SUM(0, 1, 2), SUM(1, 0, 2), SWAP(0, 1, 2), Hadamard(0, 2), Hadamard(1, 2), PHASE(0, 2), PHASE(1, 2)]  #
        i = 0

        # test all for dim = 2
        for gate in gates:
            for x0 in range(2):
                for z0 in range(2):
                    for x1 in range(2):
                        for z1 in range(2):
                            for x0p in range(2):
                                for z0p in range(2):
                                    for x1p in range(2):
                                        for z1p in range(2):
                                            i += 1
                                            p1 = PauliSum.from_string([f'x{x0}z{z0} x{x1}z{z1}'], dimensions=[2, 2])
                                            p2 = PauliSum.from_string([f'x{x0p}z{z0p} x{x1p}z{z1p}'], dimensions=[2, 2])
                                            err0 = 'In: \n' + p1.__str__() + '\n' + p2.__str__()
                                            err = (
                                                'Out: \n' +
                                                (gate.act(p1) * gate.act(p2)).__str__() +
                                                '\n' +
                                                gate.act(p1 * p2).__str__()
                                            )
                                            assert gate.act(p1) * gate.act(p2) == gate.act(p1 * p2), err0 + err + '\n'

        for dim in [3, 5, 7, 15]:
            gates = [SUM(0, 1, dim), SUM(1, 0, dim), SWAP(0, 1, dim), Hadamard(0, dim),
                     Hadamard(1, dim), PHASE(0, dim), PHASE(1, dim)]
            for gate in gates:

                for _ in range(100):
                    p1 = self.random_pauli_sum(dim, n_paulis=1)
                    p2 = self.random_pauli_sum(dim, n_paulis=1)
                    err0 = 'In: \n' + p1.__str__() + '\n' + p2.__str__()
                    err = 'Out: \n' + (gate.act(p1) * gate.act(p2)).__str__() + '\n' + gate.act(p1 * p2).__str__()
                    assert gate.act(p1) * gate.act(p2) == gate.act(p1 * p2), err0 + err + '\n'

    def test_is_symplectic(self):

        # Tests if the symplectic matrix of the gate is symplectic
        gates = [SUM(0, 1, 2), SWAP(0, 1, 2), Hadamard(0, 2), Hadamard(1, 2), PHASE(0, 2), PHASE(1, 2)]
        for gate in gates:
            assert is_symplectic(gate.symplectic, 2), (
                f"Gate {gate.name} is not symplectic. \n" +
                gate.symplectic.__str__()
            )

    def test_random_symplectic(self, num_tests=20, max_n=10, primes=[2, 3, 5, 7]):
        """
        Test random_symplectic() across several n, p values using assertions only.
        """
        for n in range(2, max_n + 1):
            for d in primes:
                for i in range(num_tests):
                    if n < 6 and d == 2:
                        index = random.randint(0, symplectic_group_size(n) - 1)
                        F = symplectic_gf2(index, n)
                    else:
                        F = symplectic_random_transvection(n, dimension=d)
                    assert is_symplectic(F, d), f"Failed symplectic check: n={n}, test {i}"

    # def test_find_symplectic_map(self):
    #     # this just tests the underlying solver, not the Gate or Pauli... implementation
    #     # see test below for implementation

    #     for n in [5, 10, 20, 50]:  # n_qudits
    #         for n_p in range(1, 10):
    #             for _ in range(100):
    #                 index = random.randint(0, symplectic_group_size(n))
    #                 F_true = symplectic_gf2(index, n)
    #                 X = np.random.randint(0, 2, size=(n_p, 2 * n), dtype=np.uint8)
    #                 Y = (X @ F_true) % 2
    #                 F_found = map_pauli_sum_to_target_tableau(X, Y)
    #                 assert np.array_equal((X @ F_found) % 2, Y)

    def test_gate_from_target(self):
        n_qudits = 4
        n_paulis = 2
        dim = 2
        dimensions = [dim] * n_qudits

        input_ps = PauliSum.from_random(n_paulis, dimensions=dimensions)
        input_ps.remove_trivial_paulis()
        input_ps.combine_equivalent_paulis()
        print(input_ps)
        circuit = Circuit.from_random(100, dimensions)
        target_ps = circuit.act(input_ps)
        # FIXME: create API to change phases and weights
        target_ps.set_phases(np.zeros(n_paulis))

        gate_from_solver = Gate.solve_from_target('ArbGate', input_ps, target_ps, dimensions)
        output_ps = gate_from_solver.act(input_ps)
        output_ps.set_phases(np.zeros(n_paulis))
        assert output_ps == target_ps, (
            'Error in gate_from_solver\n input:\n' + input_ps.__str__() +
            '\n target:\n' + target_ps.__str__() + '\n output:\n' + output_ps.__str__())

        for d in [2]:  # only solves on GF(2) for now...
            for i in range(100):
                input_ps = self.random_pauli_sum(d, n_paulis=6)
                target_ps = self.random_pauli_sum(d, n_paulis=6)

                if np.all(input_ps.symplectic_product_matrix() == target_ps.symplectic_product_matrix()):
                    print(i)
                    print('input')
                    print(input_ps.tableau)
                    print('target')
                    print(target_ps.tableau)
                    gate = Gate.solve_from_target('ArbGate', input_ps, target_ps, dimensions)
                    output_ps = gate.act(input_ps)
                    output_ps.set_phases(input_ps.phases)  # So far it does not solve for phases as well
                    print('output')
                    print(output_ps.tableau)
                    print('gate')
                    print(gate.symplectic)
                    print('check')
                    print(input_ps.tableau @ gate.symplectic % 2)

                    assert output_ps == target_ps, (
                        f'Error test {i} \n In: \n' +
                        input_ps.__str__() +
                        '\n Out: \n' +
                        output_ps.__str__() +
                        '\n Target: \n' +
                        target_ps.__str__()
                    )

    def test_gate_transvection_symplecticity(self):
        for _ in range(100):
            g = Gate.from_random(5, 2)
            gt = g.transvection(np.random.randint(0, 1, size=10))
            assert is_symplectic(gt.symplectic, 2), 'Error in transvection'

    def test_gate_inverse(self):
        # TODO: test mixed dimensions
        n_qudits = 4
        n_paulis = 4
        for dimension in [2, 3, 5, 7, 11]:
            for i in range(100):
                g = Gate.from_random(n_qudits, dimension, seed=1)
                gt = g.inv()
                rps = PauliSum.from_random(n_paulis, [dimension] * n_qudits, False, seed=1)  #
                assert rps == gt.act(
                    g.act(rps)), f'Inversion Error run {i}:\n' + rps.__str__() + '\n' + g.act(gt.act(rps)).__str__()

    def phase_table_local(self, G: Gate) -> tuple[np.ndarray, np.ndarray]:
        d = G.dimensions[0]
        phase_table_unitary = np.zeros((d, d))
        phase_table_symplectic = np.zeros((d, d))
        for i in range(d):
            for j in range(d):
                pauli_string = 'x' + str(i) + 'z' + str(j)
                ps = PauliSum.from_string([pauli_string],
                                          dimensions=[d],
                                          weights=[1], phases=[0])
                ps_m = ps.to_hilbert_space()

                ps_res = G.act(ps)
                phase_table_symplectic[i, j] = ps_res.phases[0] % (2 * d)
                # FIXME: create new PSum or open API to change phases
                ps_res.set_phases([0])
                ps_res_m = ps_res.to_hilbert_space().toarray()

                ps_m_res = (G.unitary() @ ps_m @ G.unitary().conj().T).toarray()

                mask = (ps_res_m != 0)
                factors = np.around(ps_m_res[mask] / ps_res_m[mask], 14)
                factor = factors[0]
                phase_table_unitary[i, j] = np.around(d * np.angle(factor) / (np.pi), 6) % (2 * d)

        return phase_table_unitary, phase_table_symplectic

    def phase_table_entangling(self, G: Gate) -> tuple[np.ndarray, np.ndarray]:
        d = G.dimensions[0]
        phase_table_unitary = np.zeros((d**2, d**2))
        phase_table_symplectic = np.zeros((d**2, d**2))
        U = G.unitary()
        for i1 in range(d):
            for j1 in range(d):
                for i2 in range(d):
                    for j2 in range(d):
                        pauli_string = 'x' + str(i1) + 'z' + str(j1) + ' ' + 'x' + str(i2) + 'z' + str(j2)
                        ps = PauliSum.from_string([pauli_string], dimensions=[d, d])
                        ps_m = ps.to_hilbert_space()

                        ps_res = G.act(ps)
                        phase_table_symplectic[i1 * d + i2, j1 * d + j2] = ps_res.phases[0] % (2 * d)
                        ps_res.set_phases([0])
                        ps_res_m = ps_res.to_hilbert_space().toarray()
                        ps_m_res = (U @ ps_m @ U.conj().T).toarray()
                        mask = (ps_res_m != 0)
                        factors = np.around(ps_m_res[mask] / ps_res_m[mask], 14)
                        factor = factors[0]
                        phase_table_unitary[i1 * d + i2,
                                            j1 * d + j2] = np.around(d * np.angle(factor) / (np.pi), 6) % (2 * d)

        return phase_table_unitary, phase_table_symplectic

    def test_phase_table(self):
        for d in [2, 3, 5, 7]:  # 11, 13 These work as well, but it's very slow to run
            G = Hadamard(0, d)
            phase_table_unitary, phase_table_symplectic = self.phase_table_local(G)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(diff_m), 'Symplectic phase table does not match unitary phase table for Hadamard gate'

            G = PHASE(0, d)
            phase_table_unitary, phase_table_symplectic = self.phase_table_local(G)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(diff_m), 'Symplectic phase table does not match unitary phase table for Phase gate'

            G = SUM(0, 1, d)
            phase_table_unitary, phase_table_symplectic = self.phase_table_entangling(G)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(diff_m), 'Symplectic phase table does not match unitary phase table for SUM[0,1] gate'

            G = SUM(1, 0, d)
            phase_table_unitary, phase_table_symplectic = self.phase_table_entangling(G)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(diff_m), 'Symplectic phase table does not match unitary phase table for SUM[1,0] gate'

    def test_SUM_unitary(self):
        n_tests = 100
        n_qudits = 2
        n_paulis = 2

        for d in [2, 3, 5, 7]:
            for _ in range(n_tests):
                ps = PauliSum.from_random(n_paulis, [d] * n_qudits, False, seed=1)
                ps_m = ps.to_hilbert_space()

                G = SUM(0, 1, d)
                ps_res = G.act(ps)
                ps_res_m = ps_res.to_hilbert_space()

                U_G = G.unitary()
                ps_m_res = U_G @ ps_m @ U_G.conj().T

                diff_m = np.around(ps_res_m.toarray() - ps_m_res.toarray(), 10)
                assert not np.any(diff_m), 'Symplectic phase table does not match unitary phase table for SUM gate'

    def test_two_qudit_unitary(self):
        n_tests = 10
        n_qudits = 4
        n_paulis = 2
        dims = [2, 3, 5, 7]
        gates = [SUM, SWAP]  #
        for gate in gates:
            for d in dims:
                for _ in range(n_tests):
                    dimensions = [d] * n_qudits
                    ps = PauliSum.from_random(n_paulis, dimensions, False, seed=1)
                    i, j = random.sample(range(n_qudits), 2)
                    ps_m = ps.to_hilbert_space()

                    G = gate(i, j, d)
                    ps_res = G.act(ps)
                    ps_res_m = ps_res.to_hilbert_space()

                    U_G = G.unitary(dimensions)
                    ps_m_res = U_G @ ps_m @ U_G.conj().T

                    diff_m = np.around(ps_res_m.toarray() - ps_m_res.toarray(), 10)
                    assert not np.any(diff_m), 'failed for ' + gate.__name__ + ' dimension ' + str(d)

    def test_one_qudit_unitary(self):
        n_tests = 100
        n_qudits = 1
        n_paulis = 2
        dims = [2, 3, 5, 7, 11]
        gates = [Hadamard, PHASE]
        for gate in gates:
            for d in dims:
                for _ in range(n_tests):
                    ps = PauliSum.from_random(n_paulis, [d] * n_qudits, False, seed=1)
                    ps_m = ps.to_hilbert_space()

                    G = gate(0, d)
                    ps_res = G.act(ps)
                    ps_res_m = ps_res.to_hilbert_space()

                    U_G = G.unitary()
                    ps_m_res = U_G @ ps_m @ U_G.conj().T

                    diff_m = np.around(ps_res_m.toarray() - ps_m_res.toarray(), 10)
                    assert not np.any(diff_m), 'failed for ' + gate.__name__ + ' dimension ' + str(d)

    def test_pauli_gate_unitary(self):
        n_tests = 100
        n_qudits = 1
        n_paulis = 2
        dims = [2, 3, 5, 7, 11]
        for d in dims:
            for _ in range(n_tests):
                ps = PauliSum.from_random(n_paulis, [d] * n_qudits, False, seed=1)
                ps_m = ps.to_hilbert_space()

                G = PauliGate(PauliString.from_random([d]))
                ps_res = G.act(ps)
                ps_res_m = ps_res.to_hilbert_space().toarray()

                U_G = G.unitary()
                ps_m_res = (U_G @ ps_m @ U_G.conj().T).toarray()

                diff_m = np.around(ps_res_m - ps_m_res, 10)
                print(diff_m)
                assert not np.any(diff_m), 'failed for PauliGate dimension ' + str(d)
