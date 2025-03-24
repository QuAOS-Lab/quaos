import pytest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)
from quaos.core.pauli import (
    pauli, pauli_to_matrix, string_to_pauli, XZ_mat, tensor, pauli_to_string,
    symplectic_inner_product, quditwise_inner_product, pauli_product,
)


class TestPauli:

    def __get_matrix_x(self, dimension: int = 2, exponent: int = 1) -> np.ndarray:
        pauli_x = np.zeros((dimension, dimension), dtype=complex)
        for col in range(dimension): 
            pauli_x[col][(col - exponent) % dimension] = 1
        return pauli_x

    def __get_matrix_z(self, dimension: int = 2, exponent: int = 1) -> np.ndarray:
        pauli_z = np.zeros((dimension, dimension), dtype=complex)
        for i in range(dimension):
            pauli_z[i][i] = np.exp(2 * np.pi * 1j * i * exponent / dimension)
        return pauli_z
    
    def __get_matrix_xz(self, dimension: int = 2,
                        exponent_x: int = 1,
                        exponent_z: int = 1) -> np.ndarray:
        x_matrix = self.__get_matrix_x(dimension, exponent_x)
        z_matrix = self.__get_matrix_z(dimension, exponent_z)
        return x_matrix @ z_matrix

    def __get_heterogeneous_pauli(self) -> pauli:
        return pauli(
            np.array([[1, 0, 3], [1, 1, 1]]), 
            np.array([[0, 2, 4], [1, 1, 1]]), 
            [2, 3, 4]
        ) 
    
    def __get_commuting_pauli(self) -> pauli:
        return pauli(
            np.array([[0, 2, 0], [1, 0, 1]]), 
            np.array([[1, 0, 3], [0, 1, 0]]), 
            [2, 3, 4]
        ) 

    def __get_pauli_x(self, dimension: int = 2) -> pauli:
        return pauli(np.array([[1]]), np.array([[0]]), [dimension]) 

    def __get_pauli_z(self, dimension: int = 2) -> pauli:
        return pauli(np.array([[0]]), np.array([[1]]), [dimension]) 
    
    def __assertPauliEqual(self, 
                           expected_pauli: pauli, 
                           actual_pauli: pauli, 
                           msg: str = "Paulis are not equal"):
        assert np.allclose(expected_pauli.X, actual_pauli.X), msg
        assert np.allclose(expected_pauli.Z, actual_pauli.Z), msg
        assert np.allclose(expected_pauli.phases, actual_pauli.phases), msg
        assert np.allclose(expected_pauli.dims, actual_pauli.dims), msg
        assert expected_pauli.lcm == actual_pauli.lcm, msg  

    def test_XZ_mat(self):
        expected_x_matrix = self.__get_matrix_x(4, 1)
        expected_z_matrix = self.__get_matrix_z(4, 2)
        expected_xz_matrix = self.__get_matrix_xz(4, 3, 4)
        
        x_pauli_matrix = XZ_mat(4, 1, 0).toarray()
        z_pauli_matrix = XZ_mat(4, 0, 2).toarray()
        xz_pauli_matrix = XZ_mat(4, 3, 4).toarray()
        
        assert np.allclose(x_pauli_matrix, expected_x_matrix), "Wrong X matrix"
        assert np.allclose(z_pauli_matrix, expected_z_matrix),  "Wrong Z matrix"
        assert np.allclose(xz_pauli_matrix, expected_xz_matrix), "Wrong XZ matrix"

    def test_tensor(self):
        x_pauli_csr_matrix = XZ_mat(3, 1, 0)
        z_pauli_csr_matrix = XZ_mat(2, 0, 2)
        xz_pauli_csr_matrix = XZ_mat(3, 1, 1)
        x_pauli_matrix = x_pauli_csr_matrix.toarray()
        z_pauli_matrix = z_pauli_csr_matrix.toarray()
        xz_pauli_matrix = xz_pauli_csr_matrix.toarray()

        expected_result = np.kron(
            x_pauli_matrix, np.kron(z_pauli_matrix, xz_pauli_matrix)
        )
        actual_pauli_product = tensor(
            [x_pauli_csr_matrix, z_pauli_csr_matrix, xz_pauli_csr_matrix]
        ).toarray()

        assert np.allclose(actual_pauli_product, expected_result), "Wrong tensor product"

    def test_pauli_to_matrix(self):
        pauli_x = self.__get_pauli_x(4)
        pauli_z = self.__get_pauli_z(4)
        heterogeneous_pauli = self.__get_heterogeneous_pauli()
        expected_matrix_x = self.__get_matrix_x(4)
        expected_matrix_z = self.__get_matrix_z(4)
        
        actual_matrix_x = pauli_to_matrix(pauli_x).toarray()
        actual_matrix_z = pauli_to_matrix(pauli_z).toarray()

        assert np.allclose(expected_matrix_x, actual_matrix_x), "Wrong Pauli X matrix"
        assert np.allclose(expected_matrix_z, actual_matrix_z), "Wrong Pauli Z matrix"
        with pytest.raises(Exception):
            pauli_to_matrix(heterogeneous_pauli), "Passing multiple pauli shall raise an exception"

    def test_is_IX(self):
        pauli_x = self.__get_pauli_x(4)
        assert pauli_x.is_IX(), "Pauli should be X"
        
    def test_is_IZ(self):
        pauli_z = self.__get_pauli_z(4)
        assert pauli_z.is_IZ(), "Pauli should be Z"
        
    def test_pauli_from_string(self):
        expected_heterogenous_pauli = self.__get_heterogeneous_pauli()
        
        pauli_z = string_to_pauli('x0z1')
        pauli_x = string_to_pauli('x1z0')
        actual_heterogeneous_pauli = string_to_pauli(
            ["x1z0 x0z2 x3z4", "x1z1 x1z1 x1z1"], 
            [2, 3, 4]
        )

        assert pauli_z.Z.shape[0] == 1, "Wrong Pauli Z from string"
        assert pauli_z.is_IZ(), "Pauli should be Z"

        assert pauli_x.X.shape[0] == 1, "Wrong Pauli X from string"
        assert pauli_x.is_IX(), "Pauli should be X"

        assert actual_heterogeneous_pauli.X.shape[0] == 2, "Wrong Pauli from string"
        self.__assertPauliEqual(actual_heterogeneous_pauli,
                                expected_heterogenous_pauli,
                                msg="Wrong Pauli from string")

    def test_symplectic_inner_product(self):
        not_single_pauli = self.__get_heterogeneous_pauli()
        heterogeneous_pauli = self.__get_heterogeneous_pauli()

        with pytest.raises(Exception):
            symplectic_inner_product(not_single_pauli, not_single_pauli), 
            "Passing multiple pauli shall raise an exception"
        with pytest.raises(Exception):
            symplectic_inner_product(heterogeneous_pauli, heterogeneous_pauli),
            "Passing paulis with different dimensions shall raise an exception"
        # assert np.allclose(expected_matrix_x, actual_matrix_x), 
        #                 msg="Wrong Pauli X matrix")

    def test_is_commuting(self):
        commuting_pauli = self.__get_commuting_pauli()
        not_commuting_pauli = self.__get_heterogeneous_pauli()

        assert commuting_pauli.is_commuting(), "Pauli should be commuting"
        assert not not_commuting_pauli.is_commuting(), "Pauli should not be commuting"
        
    def test_is_quditwise_commuting(self):
        # TODO: Implement test_is_quditwise_commuting
        assert True
        
    def test_a_pauli(self):
        # TODO: Implement test_a_pauli
        assert True
        
    def test_paulis(self):
        # TODO: Implement test_paulis
        assert True
        
    def test_qudits(self):
        # TODO: Implement test_qudits
        assert True
        
    def test_delete_paulis_(self):
        # TODO: Implement test_delete_paulis_
        assert True
        
    def test_delete_qudits_(self):
        # TODO: Implement test_delete_qudits_
        assert True
        
    def test_copy(self):
        # TODO: Implement test_copy
        assert True
        
    def test_print(self):
        # TODO: Implement test_print
        assert True
        
    def test_print_symplectic(self):
        # TODO: Implement test_print_symplectic
        assert True
        
    def test_string_to_pauli(self):
        # TODO: Implement test_string_to_pauli
        assert True
        
    def test_pauli_to_string(self):
        # TODO: Implement test_pauli_to_string
        assert True
        
    def test_quditwise_inner_product(self):
        # TODO: Implement test_quditwise_inner_product
        assert True
        
    def test_pauli_product(self):
        # TODO: Implement test_pauli_product
        assert True
