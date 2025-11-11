import numpy as np
import random

from sympleq.core.circuits import Gate, Circuit
from sympleq.core.circuits.gates import GATES
from sympleq.core.circuits.utils import is_symplectic
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits.random_symplectic import symplectic_gf2, \
    symplectic_group_size, symplectic_random_transvection
from sympleq.core.paulis.constants import DEFAULT_QUDIT_DIMENSION


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

    def random_pauli_string(self, dim):
        # Generates a random PauliString of dimension dim
        r1 = np.random.randint(0, dim)
        r2 = np.random.randint(0, dim)
        s1 = np.random.randint(0, dim)
        s2 = np.random.randint(0, dim)
        return PauliString.from_string(f'x{r1}z{s1} x{r2}z{s2}', dimensions=[dim, dim]), r1, r2, s1, s2

    def test_SUM(self):
        # acts the GATES.sum gate on a bunch of random pauli strings and pauli sums
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
                output_ps = GATES.sum.act(input_ps, (0, 1))
                assert output_ps == PauliString.from_string(output_str_correct, dimensions=[
                                                            d, d]), 'Error in GATES.sum gate'

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

                output_psum = GATES.sum.act(input_psum, (0, 1))
                assert output_psum == output_psum_correct, (
                    'Error in GATES.sum gate: \n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_SWAP(self):
        # acts the GATES.sum gate on a bunch of random pauli strings and pauli sums
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
                output_ps = GATES.swap.act(input_ps, (0, 1))
                assert output_ps == PauliString.from_string(output_str_correct, dimensions=[
                                                            d, d]), 'Error in GATES.sum gate'

            # test pauli_sums
            for _ in range(100):
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

                output_psum = GATES.swap.act(input_psum, (0, 1))
                assert output_psum == output_psum_correct, (
                    'Error in GATES.sum gate: \n' +
                    input_psum.__str__() + '\n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_Hadamard(self):
        # TODO: Be certain of inverse convention - ultimately arbitrary but should match prevalent literature
        for d in [2, 5, 11]:
            for _ in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{(-s1) % d}z{r1} x{r2}z{s2}"

                output_ps = GATES.H.act(input_ps, 0)
                assert output_ps == PauliString.from_string(
                    output_str_correct, dimensions=[d, d]
                ), 'Error in GATES.H gate 0'

            for _ in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{r1}z{s1} x{(-s2) % d}z{r2}"

                output_ps = GATES.H.act(input_ps, 1)
                assert output_ps == PauliString.from_string(
                    output_str_correct, dimensions=[d, d]
                ), 'Error in GATES.H gate 1'

            # test pauli_sums
            for _ in range(100):
                ps_list_in = []
                ps_list_out_correct = []
                ps_phase_out_correct = []
                for _ in range(10):
                    input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                    output_str_correct = f"x{(-s1) % d}z{r1} x{r2}z{s2}"

                    ps_list_in.append(input_ps)
                    ps_list_out_correct.append(PauliString.from_string(output_str_correct, dimensions=[d, d]))
                    ps_phase_out_correct.append(-2 * r1 * s1)

                input_psum = PauliSum.from_pauli_strings(ps_list_in)
                output_psum = GATES.H.act(input_psum, 0)
                output_psum_correct = PauliSum.from_pauli_strings(
                    ps_list_out_correct, phases=ps_phase_out_correct)
                assert output_psum == output_psum_correct, (
                    'Error in GATES.H gate: \n' +
                    input_psum.__str__() + '\n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_PHASE(self):
        # TODO: Better approach - is there another approach that we can test the first one with like GATES.H?
        for d in [2, 5, 11]:

            for i in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{r1}z{(r1 + s1) % d} x{r2}z{s2}"

                output_ps = GATES.phase.act(input_ps, 0)
                assert output_ps == PauliString.from_string(output_str_correct, dimensions=[d, d]), 'Error in H gate 0'

            for i in range(100):
                input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                output_str_correct = f"x{r1}z{s1} x{r2}z{(r2 + s2) % d}"

                output_ps = GATES.phase.act(input_ps, 1)
                assert output_ps == PauliString.from_string(output_str_correct, dimensions=[d, d]), 'Error in H gate 1'

            # test pauli_sums
            gate = GATES.phase
            h = gate.phase_vector(d)
            C = gate.symplectic()
            for i in range(100):
                ps_list_in = []
                ps_list_out_correct = []
                ps_phase_out_correct = []
                for _ in range(10):
                    input_ps, r1, r2, s1, s2 = self.random_pauli_string(d)
                    output_str_correct = f"x{r1}z{(r1 + s1) % d} x{r2}z{s2}"

                    ps_list_in.append(input_ps)
                    ps_list_out_correct.append(PauliString.from_string(output_str_correct, dimensions=[d, d]))
                    a1 = np.array([r1, s1])
                    U = np.zeros((input_ps.n_qudits(), input_ps.n_qudits()), dtype=int)
                    U[1, 0] = 1

                    p = int(
                        np.dot(h, a1) -
                        np.diag(C.T @ U @ C) @ a1 +
                        2 * a1 @ np.triu(C.T @ U @ C) @ a1 -
                        a1 @ np.diag(np.diag(C.T @ U @ C)) @ a1
                    )
                    ps_phase_out_correct.append(p)

                input_psum = PauliSum.from_pauli_strings(ps_list_in)
                output_psum = gate.act(input_psum, 0)
                output_psum_correct = PauliSum.from_pauli_strings(
                    ps_list_out_correct, phases=ps_phase_out_correct)

                assert output_psum == output_psum_correct, (
                    'Error in GATES.H gate: \n' +
                    input_psum.__str__() + '\n' +
                    output_psum.__str__() + '\n' +
                    output_psum_correct.__str__()
                )

    def test_group_homomorphism(self):
        # Tests all gates and all combinations of two pauli strings for the group homomorphism property:
        # gate.act(p1) * gate.act(p2) == gate.act(p1 * p2)

        gates_data = [(GATES.sum, 0, 1), (GATES.sum, 1, 0), (GATES.swap, 0, 1),
                      (GATES.H, 0), (GATES.H, 1), (GATES.phase, 0), (GATES.phase, 1)]
        i = 0

        # test all for dim = 2
        for (gate, *qudits) in gates_data:
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
                                                (gate.act(p1, qudits) * gate.act(p2, qudits)).__str__() +
                                                '\n' +
                                                gate.act(p1 * p2, qudits).__str__()
                                            )
                                            assert gate.act(p1, qudits) * \
                                                gate.act(p2, qudits) == gate.act(p1 * p2, qudits), err0 + err + '\n'

        for dimensions in [3, 5, 7, 15]:
            for (gate, *qudits) in gates_data:
                for _ in range(100):
                    p1 = PauliSum.from_random(1, [dimensions] * 2)
                    p2 = PauliSum.from_random(1, [dimensions] * 2)
                    err0 = 'In: \n' + p1.__str__() + '\n' + p2.__str__()
                    err = 'Out: \n' + (gate.act(p1, qudits) * gate.act(p2, qudits)
                                       ).__str__() + '\n' + gate.act(p1 * p2, qudits).__str__()
                    assert gate.act(p1, qudits) * gate.act(p2, qudits) == gate.act(p1 * p2, qudits), err0 + err + '\n'

    def test_is_symplectic(self):
        # Tests if the symplectic matrix of the gate is symplectic
        gates_data = [(GATES.sum, 0, 1), (GATES.swap, 0, 1), (GATES.H, 0),
                      (GATES.H, 1), (GATES.phase, 0), (GATES.phase, 1)]
        for (gate, *qudits) in gates_data:
            assert is_symplectic(gate.symplectic(), 2), (
                f"Gate {gate.name()} is not symplectic. \n" +
                gate.symplectic().__str__()
            )

    def test_random_symplectic(self, num_tests=20, max_n=10, primes=[2, 3, 5, 7]):
        """
        Test random_symplectic() across several n, p values using assertions only.
        """
        for n in range(2, max_n + 1):
            for d in primes:
                for i in range(num_tests):
                    if n < 6 and d == 2:
                        index = random.randint(0, symplectic_group_size(n))
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
        print(input_ps.dimensions())
        circuit = Circuit.from_random(n_qudits, 100)
        target_ps = circuit.act(input_ps)
        print(target_ps.dimensions())
        target_ps.set_phases(np.zeros(n_paulis))

        gate_from_solver = Gate.solve_from_target('ArbGate', input_ps, target_ps)
        qudits = tuple(range(n_qudits))
        output_ps = gate_from_solver.act(input_ps, qudits)
        output_ps.set_phases(np.zeros(n_paulis))
        assert output_ps == target_ps, (
            'Error in gate_from_solver\n input:\n' + input_ps.__str__() +
            '\n target:\n' + target_ps.__str__() + '\n output:\n' + output_ps.__str__())

        for dimensions in [2]:  # only solves on GF(2) for now...
            for i in range(10):
                input_ps = PauliSum.from_random(6, [dimensions] * n_qudits)
                target_ps = PauliSum.from_random(6, [dimensions] * n_qudits)

                if np.all(input_ps.symplectic_product_matrix() == target_ps.symplectic_product_matrix()):
                    print(i)
                    print('input')
                    print(input_ps.tableau())
                    print('target')
                    print(target_ps.tableau())
                    gate = Gate.solve_from_target('ArbGate', input_ps, target_ps)
                    qudits = tuple(range(n_qudits))
                    output_ps = gate.act(input_ps, qudits)
                    output_ps.set_phases(input_ps.phases())  # So far it does not solve for phases as well
                    print('output')
                    print(output_ps.tableau())
                    print('gate')
                    print(gate.symplectic)
                    print('check')
                    print(input_ps.tableau() @ gate.symplectic % 2)

                    assert output_ps == target_ps, (
                        f'Error test {i} \n In: \n' +
                        input_ps.__str__() +
                        '\n Out: \n' +
                        output_ps.__str__() +
                        '\n Target: \n' +
                        target_ps.__str__()
                    )

    def test_gate_transvection(self):

        for _ in range(100):
            g = Gate.from_random(5, 2)
            gt = g.transvection(np.random.randint(0, 1, size=10))
            assert is_symplectic(gt.symplectic(), 2), 'Error in transvection'

    # def test_gate_inverse(self):
    #     n_qudits = 2
    #     n_paulis = 3
    #     dimension = 3
    #     for _ in range(1):
    #         g = Gate.from_random(n_qudits, dimension)
    #         gt = g.inv()
    #         rps = PauliSum.from_random(n_paulis, n_qudits, [dimension] * n_qudits, False, seed=1)
    #         print(rps)
    #         assert rps == g.act(gt.act(rps)), 'Inversion Error:\n' + \
    #               rps.__str__() + '\n' + g.act(gt.act(rps)).__str__()

    def phase_table_local(self, G: Gate, dimensions: int, *qudits: int) -> tuple[np.ndarray, np.ndarray]:
        # dimensions = random.randint(DEFAULT_QUDIT_DIMENSION, 7)
        lcm = np.lcm.reduce(dimensions)

        phase_table_unitary = np.zeros((dimensions, dimensions))
        phase_table_symplectic = np.zeros((dimensions, dimensions))
        for i in range(dimensions):
            for j in range(dimensions):
                pauli_string = 'x' + str(i) + 'z' + str(j)
                ps = PauliSum.from_string(pauli_string, dimensions, weights=[1], phases=[0])
                ps_m = ps.to_hilbert_space()

                ps_res = G.act(ps, qudits)
                phase_table_symplectic[i, j] = ps_res.phases()[0] % (2 * lcm)
                ps_res.set_phases([0])
                ps_res_m = ps_res.to_hilbert_space().toarray()

                # FIXME: accept non-list inputs in unitary
                gate_unitary = G.unitary(qudits, [dimensions] * G.n_qudits())
                ps_m_res = (gate_unitary @ ps_m @ gate_unitary.conj().T).toarray()

                mask = (ps_res_m != 0)
                factors = np.around(ps_m_res[mask] / ps_res_m[mask], 14)
                factor = factors[0]
                phase_table_unitary[i, j] = np.around(dimensions * np.angle(factor) / (np.pi), 6) % (2 * lcm)

        return phase_table_unitary, phase_table_symplectic

    def phase_table_entangling(self, G: Gate, *qudits: int) -> tuple[np.ndarray, np.ndarray]:
        dimensions = random.randint(DEFAULT_QUDIT_DIMENSION, 7)
        lcm = np.lcm.reduce(dimensions)

        phase_table_unitary = np.zeros((dimensions**2, dimensions**2))
        phase_table_symplectic = np.zeros((dimensions**2, dimensions**2))
        U = G.unitary(qudits, [dimensions] * G.n_qudits())
        for i1 in range(dimensions):
            for j1 in range(dimensions):
                for i2 in range(dimensions):
                    for j2 in range(dimensions):
                        pauli_string = 'x' + str(i1) + 'z' + str(j1) + ' ' + 'x' + str(i2) + 'z' + str(j2)
                        ps = PauliSum.from_string([pauli_string], dimensions=[dimensions, dimensions])
                        ps_m = ps.to_hilbert_space()

                        ps_res = G.act(ps, qudits)
                        phase_table_symplectic[i1 * dimensions + i2, j1 *
                                               dimensions + j2] = ps_res.phases()[0] % (2 * lcm)
                        ps_res.set_phases([0])
                        ps_res_m = ps_res.to_hilbert_space().toarray()
                        ps_m_res = (U @ ps_m @ U.conj().T).toarray()
                        mask = (ps_res_m != 0)
                        factors = np.around(ps_m_res[mask] / ps_res_m[mask], 14)
                        factor = factors[0]
                        phase_table_unitary[i1 * dimensions + i2, j1 * dimensions + j2] = \
                            np.around(dimensions * np.angle(factor) / (np.pi), 6) % (2 * lcm)

        return phase_table_unitary, phase_table_symplectic

    def test_phase_table(self):
        for d in [2, 3, 5, 7]:  # 11, 13 These work as well, but it's very slow to run
            phase_table_unitary, phase_table_symplectic = self.phase_table_local(GATES.H, d, 0)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(diff_m), 'Symplectic phase table does not match unitary phase table for GATES.H gate'

            phase_table_unitary, phase_table_symplectic = self.phase_table_local(GATES.phase, d, 0)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(diff_m), 'Symplectic phase table does not match unitary phase table for Phase gate'

            phase_table_unitary, phase_table_symplectic = self.phase_table_entangling(GATES.sum, 0, 1)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(
                diff_m), 'Symplectic phase table does not match unitary phase table for GATES.sum[0,1] gate'

            phase_table_unitary, phase_table_symplectic = self.phase_table_entangling(GATES.sum, 1, 0)
            diff_m = np.around(phase_table_unitary - phase_table_symplectic, 10)
            assert not np.any(
                diff_m), 'Symplectic phase table does not match unitary phase table for GATES.sum[1,0] gate'
