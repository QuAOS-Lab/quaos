import numpy as np
import pytest
from numpy.random import default_rng
from sympleq import bases_to_int, int_to_bases
from sympleq.utils import complex_phase_value, get_linearly_independent_rows, tensor
from tests import PRIME_LIST, choose_random_dimensions


N_tests = 30
rng = default_rng()


class TestUtils:
    def test_bases_to_int(self):
        for _ in range(N_tests):
            dimensions = choose_random_dimensions(2500)
            bases = [rng.integers(0, d - 1) for d in dimensions]
            digit = bases_to_int(bases, dimensions)

            assert np.all(int_to_bases(digit, dimensions) == bases), (
                f"Expected int_to_bases to be the inverse of bases_to_int.\n"
                f"Dimensions: {dimensions}\n"
                f"Bases: {bases}\n"
                f"Digit: {digit}"
            )

    def test_int_to_bases(self):
        for _ in range(N_tests):
            dimensions = choose_random_dimensions(2500)
            max_val = np.prod(dimensions)
            digit = rng.integers(0, int(max_val - 1))
            bases = int_to_bases(digit, dimensions)

            assert np.all(bases_to_int(bases, dimensions) == digit), (
                f"Expected bases_to_int to be the inverse of int_to_bases.\n"
                f"Dimensions: {dimensions}\n"
                f"Bases: {bases}\n"
                f"Digit: {digit}"
            )

    def test_get_linearly_independent_rows(self):
        n = 10
        d = 5
        for n in [10, 15, 20]:
            for d in PRIME_LIST:
                id = np.eye(n, dtype=int)

                assert get_linearly_independent_rows(id, d) == np.arange(n).tolist()

                for i in range(100):
                    # add dependent rows to id, check the original independent rows are the only ones obtained
                    id = np.vstack([id, rng.integers(0, d, n)])
                    assert get_linearly_independent_rows(id, d) == np.arange(n).tolist()

    def test_complex_phase_value_matches_reference_formula(self):
        for dimension in PRIME_LIST:
            for phase in range(2 * dimension):
                expected = np.exp(2 * np.pi * 1j * phase / (2 * dimension))
                assert np.isclose(complex_phase_value(phase, dimension), expected)

    def test_complex_phase_value_periodicity(self):
        for dimension in PRIME_LIST:
            for phase in range(-2 * dimension, 2 * dimension):
                assert np.isclose(
                    complex_phase_value(phase, dimension),
                    complex_phase_value(phase + 2 * dimension, dimension),
                )

    def test_complex_phase_value_exact_quadrant_roots(self):
        # dimension=2 is the qubit case: (2 * phase) % dimension == 0 always holds,
        # so every phase must land exactly (no floating-point roundoff) on a quadrant root.
        exact_roots = {1 + 0j, 1j, -1 + 0j, -1j}
        for phase in range(8):
            value = complex_phase_value(phase, 2)
            assert value in exact_roots

        expected_by_phase = {0: 1 + 0j, 1: 1j, 2: -1 + 0j, 3: -1j}
        for phase, expected in expected_by_phase.items():
            assert complex_phase_value(phase, 2) == expected

    def test_tensor_single_matrix_passthrough(self):
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        assert np.array_equal(tensor([X]), X)

    def test_tensor_matches_sequential_np_kron(self):
        rng_local = default_rng(0)
        matrices = [rng_local.normal(size=(2, 2)) + 1j * rng_local.normal(size=(2, 2)) for _ in range(4)]

        expected = matrices[0]
        for m in matrices[1:]:
            expected = np.kron(expected, m)

        assert np.allclose(tensor(matrices).toarray(), expected)

    def test_tensor_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="At least one matrix"):
            tensor([])

    def test_tensor_raises_on_non_square_matrix(self):
        with pytest.raises(ValueError, match="square"):
            tensor([np.zeros((2, 3), dtype=complex)])

    def test_tensor_raises_on_non_2d_array(self):
        with pytest.raises(ValueError, match="square"):
            tensor([np.zeros(4, dtype=complex)])

    # Do we need these tests? Why do we need only complex datatypes?

    # def test_multi_kron_raises_on_non_complex_dtype(self):
    #     with pytest.raises(ValueError, match="complex dtype"):
    #         multi_kron([np.eye(2, dtype=float)])

    # def test_multi_kron_raises_on_mixed_valid_and_invalid_matrices(self):
    #     valid = np.eye(2, dtype=complex)
    #     invalid = np.eye(2, dtype=float)
    #     with pytest.raises(ValueError, match="complex dtype"):
    #         multi_kron([valid, invalid])
