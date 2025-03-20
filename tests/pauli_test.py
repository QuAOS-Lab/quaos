import unittest
import sys
import os

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)
from quaos.pauli import string_to_pauli


class TestPauli(unittest.TestCase):

    def test_string_to_pauli(self):
        pauli_from_string = string_to_pauli('x1z0')
        self.assertEqual(
            pauli_from_string.X.shape[0], 1, msg="Wrong Pauli X from string"
        )
        self.assertTrue(
            pauli_from_string.is_IX(), msg="Pauli should be X!"
        )
