import unittest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)
from quaos.core.pauli import (
    pauli, 
    pauli_to_matrix,
    string_to_pauli,
)


class TestPauli(unittest.TestCase):

    def __get_matrix_x(self, dimension: int = 2) -> np.ndarray:
        pauli_x = np.zeros((dimension, dimension), dtype=complex)
        for i in range(dimension):
            pauli_x[i][(dimension - 1 + i) % dimension] = 1
        return pauli_x

    def __get_matrix_z(self, dimension: int = 2) -> np.ndarray:
        pauli_z = np.zeros((dimension, dimension), dtype=complex)
        for i in range(dimension):
            pauli_z[i][i] = np.exp(2 * np.pi * 1j * i / dimension)
        return pauli_z

    def __get_pauli_x(self, dimension: int = 2) -> pauli:
        return string_to_pauli("x1z0", dimension)

    def __get_pauli_z(self, dimension: int = 2) -> pauli:
        return string_to_pauli("x0z1", dimension)

    def test_pauli_to_matrix(self):
        pauli_x = self.__get_pauli_x(4)
        actual_matrix_x = pauli_to_matrix(pauli_x)

        pauli_z = self.__get_pauli_z(4)
        actual_matrix_z = pauli_to_matrix(pauli_z)

        expected_matrix_x = self.__get_matrix_x(4)
        self.assertTrue(np.all(expected_matrix_x == actual_matrix_x), 
                        msg="Wrong Pauli X matrix")
        expected_matrix_z = self.__get_matrix_z(4)
        self.assertTrue(np.all(expected_matrix_z == actual_matrix_z), 
                        msg="Wrong Pauli Z matrix")

    def test_pauli_z_from_string(self):
        test_pauli = string_to_pauli('x0z1')

        self.assertEqual(
            test_pauli.Z.shape[0], 1, msg="Wrong Pauli Z from string")

        self.assertTrue(
            test_pauli.is_IZ(), msg="Pauli should be Z!")

    def test_pauli_x_from_string(self):
        test_pauli = string_to_pauli('x1z0')

        self.assertEqual(
            test_pauli.X.shape[0], 1, msg="Wrong Pauli X from string")

        self.assertTrue(
            test_pauli.is_IX(), msg="Pauli should be X!")
