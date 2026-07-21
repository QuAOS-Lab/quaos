import numpy as np
from numpy.random import default_rng
from sympleq import bases_to_int, int_to_bases
from sympleq.utils import get_linearly_independent_rows
from tests import choose_random_dimensions


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
            for d in [2, 5, 11, 17]:
                id = np.eye(n, dtype=int)

                assert get_linearly_independent_rows(id, d) == np.arange(n).tolist()

                for i in range(100):
                    # add dependent rows to id, check the original independent rows are the only ones obtained
                    id = np.vstack([id, rng.integers(0, d, n)])
                    assert get_linearly_independent_rows(id, d) == np.arange(n).tolist()
