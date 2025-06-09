import sys
sys.path.append("./")
from quaos.circuits import Gate, SUM, SWAP, Hadamard, PHASE, CNOT
from quaos.paulis import PauliSum, PauliString
import numpy as np


class TestGates():

    @staticmethod
    def states(dim):
        ps1 = PauliString('x1z0 x0z0', dimensions=[dim, dim])
        ps2 = PauliString('x0z0 x0z1', dimensions=[dim, dim])
        ps3 = PauliString('x1z0 x1z0', dimensions=[dim, dim])
        ps4 = PauliString('x0z1 x0z1', dimensions=[dim, dim])
        ps5 = PauliString('x1z1 x0z0', dimensions=[dim, dim])
        p_sum = PauliSum(['x1z0 x0z0 x1z1', 'x0z0 x0z1 x1z0', 'x1z1 x1z0 x0z0'], dimensions=[dim, dim, dim])
        return ps1, ps2, ps3, ps4, ps5, p_sum

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

                input_ps = PauliString(input_str, dimensions=[d, d])
                output_ps = SUM(0, 1, d).act(input_ps)
                assert output_ps == PauliString(output_str_correct, dimensions=[d, d]), 'Error in SUM gate'

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

                    input_ps = PauliString(input_str, dimensions=[d, d])
                    output_ps_correct = PauliString(output_str_correct, dimensions=[d, d])
                    ps_list_in.append(input_ps)
                    ps_list_out_correct.append(output_ps_correct)
                    ps_phase_out_correct.append((r1 * s2) % d)

                input_psum = PauliSum(ps_list_in, dimensions=[d, d], standardise=False)
                output_psum_correct = PauliSum(ps_list_out_correct, phases=ps_phase_out_correct, dimensions=[d, d], standardise=False)
                
                output_psum = SUM(0, 1, d).act(input_psum)
                assert output_psum == output_psum_correct, 'Error in SUM gate: \n' + output_psum.__str__() + '\n' + output_psum_correct.__str__()

    def test_SWAP(self):
        pass


    def test_Hadamard(self):
        pass


    def test_PHASE(self):
        pass

    def test_arbitrary_gate(self):
        pass

    def test_group_homomorphism(self):
        gate = SUM(0, 1, 2)
        for x0 in range(2):
            for z0 in range(2):
                for x1 in range(2):
                    for z1 in range(2):
                        for x0p in range(2):
                            for z0p in range(2):
                                for x1p in range(2):
                                    for z1p in range(2):
                                        p1 = PauliSum([f'x{x0}z{z0} x{x1}z{z1}'], dimensions=[2, 2])
                                        p2 = PauliSum([f'x{x0p}z{z0p} x{x1p}z{z1p}'], dimensions=[2, 2])
                                        err0 = 'In: \n' + p1.__str__() + '\n' + p2.__str__()
                                        err = 'Out: \n' + (gate.act(p1) * gate.act(p2)).__str__() + '\n' + gate.act(p1 * p2).__str__()
                                        assert gate.act(p1) * gate.act(p2) == gate.act(p1 * p2), err0 + err


if __name__ == "__main__":
    # import pytest

    # # Run the tests in this file
    # pytest.main(["-v", __file__])
    # dim = 2

    # sum_gate = SUM(0, 1, dimension=dim)
    # print("Gate name:", sum_gate.name)
    # print("Symplectic matrix:\n", sum_gate.symplectic)

    # swap_gate = SWAP(0, 1, dimension=dim)
    # print("Gate name:", swap_gate.name)
    # print("Symplectic matrix:\n", swap_gate.symplectic)

    # # Example PauliString
    # # ps = PauliString("x0z1 x0z2", dimensions=[dim, dim])
    # # print("PauliString:", ps)
    # # print("Sum ps:", sum_gate.act(ps)[0])

    # # # testing gates

    # ps1 = PauliString('x1z0 x0z0', dimensions=[dim, dim])
    # ps2 = PauliString('x0z0 x0z1', dimensions=[dim, dim])
    # ps3 = PauliString('x1z0 x1z0', dimensions=[dim, dim])
    # ps4 = PauliString('x0z1 x0z1', dimensions=[dim, dim])
    # ps5 = PauliString('x1z1 x0z0', dimensions=[dim, dim])

    # print(ps1, '->', sum_gate.act(ps1))
    # print(ps2, '->', sum_gate.act(ps2))
    # print(ps3, '->', sum_gate.act(ps3))
    # print(ps4, '->', sum_gate.act(ps4))

    # print('Hadamard')
    # H = Hadamard(0, 2)
    # print(ps1, '->', H.act(ps1))
    # print(ps2, '->', H.act(ps2))
    # print(ps3, '->', H.act(ps3))
    # print(ps4, '->', H.act(ps4))
    # print(ps5, '->', H.act(ps5))

    # print('SWAP')
    # print(ps1, '->', swap_gate.act(ps1))
    # print(ps2, '->', swap_gate.act(ps2))
    # print(ps3, '->', swap_gate.act(ps3))
    # print(ps4, '->', swap_gate.act(ps4))
                # test pauli_sums
    d = 2
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
            print(input_str, '->', output_str_correct)

            input_ps = PauliString(input_str, dimensions=[d, d])
            output_ps_correct = PauliString(output_str_correct, dimensions=[d, d])
            ps_list_in.append(input_ps)
            ps_list_out_correct.append(output_ps_correct)
            ps_phase_out_correct.append((r1 * s2) % d)

        input_psum = PauliSum(ps_list_in, dimensions=[d, d], standardise=False)
        output_psum_correct = PauliSum(ps_list_out_correct, phases=ps_phase_out_correct, dimensions=[d, d], standardise=False)
        
        # print(input_psum)
        print('correct: \n')
        print(output_psum_correct)
        output_psum = SUM(0, 1, d).act(input_psum)
        print('output: \n')
        print(output_psum)
        assert output_psum == output_psum_correct  # , ('\n' + output_psum.__str__() + '\n' + PauliSum(ps_list_out_correct, phases=ps_phase_out_correct, dimensions=[d, d]).__str__())
