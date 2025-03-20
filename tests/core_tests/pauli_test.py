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
    XZ_mat
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

    def __get_complex_pauli(self) -> np.ndarray:
        return string_to_pauli(
            ["x1z0", "x1z0"], 
            [2, 3, 4]
        ) 

    def __get_heterogeneous_pauli_matrix(self) -> np.ndarray:
        return None

    def __get_pauli_x(self, dimension: int = 2) -> pauli:
        return string_to_pauli("x1z0", dimension)

    def __get_pauli_z(self, dimension: int = 2) -> pauli:
        return string_to_pauli("x0z1", dimension)

    def test_XZ_mat(self):
        # x_pauli_matrix = XZ_mat(4, 1, 0)
        # z_pauli_matrix = XZ_mat(4, 0, 1)
        # xz_pauli_matrix = XZ_mat(4, 1, 1)
        self.assertTrue(True)

    def test_pauli_from_string(self):
        # pauli_z = string_to_pauli('x0z1')
        # pauli_x = string_to_pauli('x1z0')
        # complex_pauli = string_to_pauli(
        #     ["x1z1", "x1z0", "x0z1"], 
        #     [2, 3, 4]
        # )
# 
        # self.assertEqual(pauli_z.Z.shape[0], 1, msg="Wrong Pauli Z from string")
        # self.assertTrue(pauli_z.is_IZ(), msg="Pauli should be Z!")
# 
        # self.assertEqual(pauli_x.X.shape[0], 1, msg="Wrong Pauli X from string")
        # self.assertTrue(pauli_x.is_IX(), msg="Pauli should be X!")
# 
        # expected_complex_pauli_matrix = self.__get_complex_pauli_matrix()
        # self.assertEqual(complex_pauli.X.shape[0], 2, 
        #                  msg="Wrong Pauli from string")
        # complex_pauli_matrix = self.__get_complex_pauli_matrix()
        # self.assertTrue(
        #     np.allclose(complex_pauli_matrix, expected_complex_pauli_matrix)
        # )
        self.assertTrue(True)

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
        
    def test_is_IX(self):
        pauli_x = self.__get_pauli_x(4)
        self.assertTrue(pauli_x.is_IX(), msg="Pauli should be X!")
        
    def test_is_IZ(self):
        pauli_z = self.__get_pauli_z(4)
        self.assertTrue(pauli_z.is_IZ(), msg="Pauli should be Z!")
        
    def test_is_commuting(self):
        # TODO: Implement test_is_commuting
        self.assertTrue(True)
        
    def test_is_quditwise_commuting(self):
        # TODO: Implement test_is_quditwise_commuting
        self.assertTrue(True)
        
    def test_a_pauli(self):
        # TODO: Implement test_a_pauli
        self.assertTrue(True)
        
    def test_paulis(self):
        # TODO: Implement test_paulis
        self.assertTrue(True)
        
    def test_qudits(self):
        # TODO: Implement test_qudits
        self.assertTrue(True)
        
    def test_delete_paulis_(self):
        # TODO: Implement test_delete_paulis_
        self.assertTrue(True)
        
    def test_delete_qudits_(self):
        # TODO: Implement test_delete_qudits_
        self.assertTrue(True)
        
    def test_copy(self):
        # TODO: Implement test_copy
        self.assertTrue(True)
        
    def test_print(self):
        # TODO: Implement test_print
        self.assertTrue(True)
        
    def test_print_symplectic(self):
        # TODO: Implement test_print_symplectic
        self.assertTrue(True)
        
    def test_string_to_pauli(self):
        # TODO: Implement test_string_to_pauli
        self.assertTrue(True)
        
    def test_pauli_to_string(self):
        # TODO: Implement test_pauli_to_string
        self.assertTrue(True)
        
    def test_symplectic_inner_product(self):
        # TODO: Implement test_symplectic_inner_product
        self.assertTrue(True)
        
    def test_quditwise_inner_product(self):
        # TODO: Implement test_quditwise_inner_product
        self.assertTrue(True)
        
    def test_pauli_product(self):
        # TODO: Implement test_pauli_product
        self.assertTrue(True)
        
    def test_tensor(self):
        # TODO: Implement test_tensor
        self.assertTrue(True)
