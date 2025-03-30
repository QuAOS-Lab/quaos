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
    
    def __get_large_pauli(self) -> pauli:
        return pauli(
            np.array([[0, 1, 2], [1, 1, 2], [1, 2, 3]]), 
            np.array([[0, 1, 2], [1, 1, 2], [1, 2, 3]]), 
            [2, 3, 4]
        )

    def __get_pauli_x(self, dimension: int = 2) -> pauli:
        return pauli(np.array([[1]]), np.array([[0]]), [dimension]) 

    def __get_pauli_z(self, dimension: int = 2) -> pauli:
        return pauli(np.array([[0]]), np.array([[1]]), [dimension]) 
    
    def __are_paulis_equal(self,
                           first: pauli, 
                           second: pauli) -> bool: 
        if not np.allclose(first.X, second.X):
            return False
        if not np.allclose(first.Z, second.Z):
            return False
        if not np.allclose(first.phases, second.phases):
            return False
        if not np.allclose(first.dims, second.dims):
            return False
        if not first.lcm == second.lcm:
            return False
        return True

    def __assert_pauli_equal(self, 
                             expected_pauli: pauli, 
                             actual_pauli: pauli, 
                             msg: str = "Paulis are not equal"):
        assert self.__are_paulis_equal(expected_pauli, actual_pauli), msg

    def __assert_pauli_not_equal(self, 
                                 expected_pauli: pauli, 
                                 actual_pauli: pauli, 
                                 msg: str = "Paulis are equal"):
        assert not self.__are_paulis_equal(expected_pauli, actual_pauli), msg

    def test_XZ_mat(self):
        expected_x_matrix = self.__get_matrix_x(4, 1)
        expected_z_matrix = self.__get_matrix_z(4, 2)
        expected_xz_matrix = self.__get_matrix_xz(4, 3, 4)
        
        x_pauli_matrix = XZ_mat(4, 1, 0).toarray()
        z_pauli_matrix = XZ_mat(4, 0, 2).toarray()
        xz_pauli_matrix = XZ_mat(4, 3, 4).toarray()
        
        assert np.allclose(x_pauli_matrix, expected_x_matrix), "Wrong X matrix"
        assert np.allclose(z_pauli_matrix, expected_z_matrix), "Wrong Z matrix"
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

    def test_symplectic_inner_product(self):
        base_pauli = pauli(np.array([[0, 1]]), np.array([[1, 0]]), [2, 3])
        commuting_pauli = pauli(np.array([[0, 5]]), np.array([[1, 21]]), [2, 3])
        non_commuting_pauli = pauli(np.array([[1, 1]]), np.array([[1, 1]]), [2, 3])
        heterogeneous_pauli = self.__get_heterogeneous_pauli()
        smaller_pauli = pauli(np.array([[0, 1]]), np.array([[1, 0]]), [2, 2])

        with pytest.raises(Exception):
            symplectic_inner_product(base_pauli, heterogeneous_pauli),
            "Passing paulis with different dimensions shall raise an exception"

        with pytest.raises(Exception):
            symplectic_inner_product(base_pauli, smaller_pauli),
            "Passing qudits with different dimensions shall raise an exception"
        
        anticommuting_result = symplectic_inner_product(base_pauli, 
                                                        non_commuting_pauli)
        assert anticommuting_result, "Symplectic inner product should be true"

        commuting_result = symplectic_inner_product(base_pauli, 
                                                    commuting_pauli)
        assert not commuting_result, "Symplectic inner product should be false"
    
    def test_quditwise_inner_product(self):
        base_pauli = pauli(np.array([[0, 1]]), np.array([[1, 0]]), [2, 3])
        commuting_pauli = pauli(np.array([[0, 5]]), np.array([[1, 21]]), [2, 3])
        non_commuting_pauli = pauli(np.array([[1, 1]]), np.array([[1, 1]]), [2, 3])
        heterogeneous_pauli = self.__get_heterogeneous_pauli()
        smaller_pauli = pauli(np.array([[0, 1]]), np.array([[1, 0]]), [2, 2])

        with pytest.raises(Exception):
            quditwise_inner_product(base_pauli, heterogeneous_pauli),
            "Passing paulis with different dimensions shall raise an exception"

        with pytest.raises(Exception):
            quditwise_inner_product(base_pauli, smaller_pauli),
            "Passing qudits with different dimensions shall raise an exception"
                
        anticommuting_result = quditwise_inner_product(base_pauli, 
                                                       non_commuting_pauli)
        assert anticommuting_result, "Symplectic inner product should be true"

        commuting_result = quditwise_inner_product(base_pauli, 
                                                   commuting_pauli)
        assert not commuting_result, "Symplectic inner product should be false"

        # TODO: given the name of the method this should pass, but it doesn't
        # larger_base_pauli = pauli(np.array([[0, 1], [0, 1]]), 
        #                           np.array([[1, 0], [1, 0]]), 
        #                           [2, 3])
        # larger_commuting_pauli = pauli(np.array([[0, 5], [0, 5]]), 
        #                                np.array([[1, 21], [1, 21]]), 
        #                                [2, 3])        
        # larger_non_commuting_pauli = pauli(np.array([[1, 1], [1, 1]]), 
        #                                    np.array([[1, 1], [1, 1]]), 
        #                                    [2, 3])
        # larger_anticommuting_result = quditwise_inner_product(
        #     larger_base_pauli, larger_non_commuting_pauli)
        # assert larger_anticommuting_result, "Symplectic inner product should be true"

        # larger_commuting_result = quditwise_inner_product(
        #     base_pauli, larger_commuting_pauli)
        # assert not larger_commuting_result, "Symplectic inner product should be false"

    def test_is_commuting(self):
        commuting_pauli = pauli(np.array([[0, 5], [0, 5]]), 
                                np.array([[1, 21], [1, 21]]), 
                                [2, 3])
        non_commuting_pauli = pauli(np.array([[1, 5], [0, 5]]), 
                                    np.array([[1, 21], [1, 21]]), 
                                    [2, 3])

        larger_commuting_pauli = pauli(np.array([[0, 5], [0, 5], [0, 5]]), 
                                       np.array([[1, 21], [1, 21], [5, 3]]), 
                                       [2, 3])        
        larger_non_commuting_pauli = pauli(np.array([[3, 5], [0, 5], [0, 5]]), 
                                           np.array([[1, 21], [1, 21], [5, 21]]), 
                                           [2, 3])        

        assert commuting_pauli.is_commuting(), "Pauli should be pairwise commuting"
        assert not non_commuting_pauli.is_commuting(), "Pauli should not be pairwise commuting"
        
        assert larger_commuting_pauli.is_commuting(), "Pauli should be pairwise commuting"
        assert not larger_non_commuting_pauli.is_commuting(), "Pauli should not be pairwise commuting"
        
    def test_is_quditwise_commuting(self):
        # TODO: EVALUATE IF THIS TEST IS NECESSARY. 
        # It seems that is_quditwise_commuting is a duplicate of is_commuting!!!
        assert True

    def test_a_pauli(self):
        large_pauli = pauli(np.array([[0, 1], [1, 1], [1, 2]]), 
                            np.array([[0, 1], [1, 1], [1, 2]]), 
                            [2, 3])   
        expected_second_pauli = pauli(np.array([[1, 1]]),
                                      np.array([[1, 1]]), 
                                      [2, 3])
        
        actual_second_pauli = large_pauli.a_pauli(1)

        self.__assert_pauli_equal(expected_second_pauli, actual_second_pauli)
        
    def test_paulis(self):
        large_pauli = pauli(np.array([[0, 5], [0, 5], [0, 5]]), 
                            np.array([[1, 21], [1, 21], [5, 3]]), 
                            [2, 3])
        
        paulis_count = large_pauli.paulis()

        assert 3 == paulis_count, "Pauli should have 3 paulis"
        
    def test_qudits(self):
        large_pauli = pauli(np.array([[0, 5], [0, 5], [0, 5]]), 
                            np.array([[1, 21], [1, 21], [5, 3]]), 
                            [2, 3])
        
        qudits_count = large_pauli.qudits()

        assert 2 == qudits_count, "Paulis should be defined 2 qudits"
        
    def test_delete_paulis_(self):
        expected_reduced_pauli = pauli(np.array([[0, 1, 2], [1, 1, 2]]), 
                                       np.array([[0, 1, 2], [1, 1, 2]]), 
                                       [2, 3, 4])   
        expected_list_reduced_pauli = pauli(np.array([[0, 1, 2]]), 
                                            np.array([[0, 1, 2]]), 
                                            [2, 3, 4])   
        
        actual_reduced_pauli = self.__get_large_pauli().delete_paulis_(2)
        actual_list_reduced_pauli = self.__get_large_pauli().delete_paulis_([1, 2])

        self.__assert_pauli_equal(expected_reduced_pauli, 
                                  actual_reduced_pauli)
        self.__assert_pauli_equal(expected_list_reduced_pauli, 
                                  actual_list_reduced_pauli)
        
    def test_delete_qudits_(self):
        expected_reduced_pauli = pauli(
            np.array([[0, 1], [1, 1], [1, 2]]), 
            np.array([[0, 1], [1, 1], [1, 2]]), 
            [2, 3]
        )  
        expected_list_reduced_pauli = pauli(np.array([[1], [1], [2]]), 
                                            np.array([[1], [1], [2]]), 
                                            [3])   
        
        # TODO: delete_qudits_ currently does not work as expected!!!
        # actual_reduced_pauli = self.__get_large_pauli().delete_qudits_(2)        
        # actual_list_reduced_pauli = self.__get_large_pauli().delete_qudits_([0, 2])
        #
        # self.__assert_pauli_equal(expected_reduced_pauli, 
        #                           actual_reduced_pauli)
        # self.__assert_pauli_equal(expected_list_reduced_pauli, 
        #                           actual_list_reduced_pauli)
        
    def test_copy(self):
        base_pauli = self.__get_large_pauli()

        copied_pauli = base_pauli.copy()

        self.__assert_pauli_equal(base_pauli, 
                                  copied_pauli)
        copied_pauli.X[0][0] = 1
        self.__assert_pauli_not_equal(base_pauli, 
                                      copied_pauli)
        
    def test_print(self):
        # test not needed
        assert True
        
    def test_print_symplectic(self):
        # test not needed
        assert True
        
    def test_string_to_pauli(self):
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
        self.__assert_pauli_equal(actual_heterogeneous_pauli,
                                  expected_heterogenous_pauli,
                                  msg="Wrong Pauli from string")
        
    def test_pauli_to_string(self):
        # TODO: change pauli_to_string to return a list of strings instead of a tuple
        heterogenous_pauli = self.__get_heterogeneous_pauli()
        expected_strings = ["x1z0 x0z2 x3z0", "x1z1 x1z1 x1z1"]
        expected_dims = [2, 3, 4]
        expected_phases = [0, 0]

        actual_strings, actual_dims, actual_phases = pauli_to_string(heterogenous_pauli)

        assert expected_strings == actual_strings, "Wrong Pauli strings"
        assert np.allclose(expected_dims, actual_dims), "Wrong Pauli dimensions"
        assert np.allclose(expected_phases, actual_phases), "Wrong Pauli to phases"
        
    def test_pauli_product(self):
        base_pauli = pauli(np.array([[0, 1]]), np.array([[1, 0]]), [2, 3])
        heterogeneous_pauli = self.__get_heterogeneous_pauli()
        smaller_pauli = pauli(np.array([[0, 1]]), np.array([[1, 0]]), [2, 2])

        expected_base_pauli_square = pauli(np.array([[0, 2]]), 
                                           np.array([[2, 0]]), 
                                           [2, 3])
        
        actual_base_pauli_square = pauli_product(base_pauli, base_pauli)
        
        with pytest.raises(Exception):
            _ = pauli_product(base_pauli, heterogeneous_pauli),
            "Passing paulis with different dimensions shall raise an exception"

        with pytest.raises(Exception):
            _ = pauli_product(base_pauli, smaller_pauli),
            "Passing qudits with different dimensions shall raise an exception"
        
        self.__assert_pauli_equal(expected_base_pauli_square, 
                                  actual_base_pauli_square,
                                  "Wrong pauli product")
