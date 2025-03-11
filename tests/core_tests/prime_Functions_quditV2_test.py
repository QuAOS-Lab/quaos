import unittest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)

from quaos.core.prime_Functions_Andrew import (
    ground_state
)


class TestPrimeFunctionAndrew(unittest.TestCase):

    def __get_matrix_x(self, dimension: int = 2) -> np.ndarray:
        pass
