import pytest
import numpy as np
import sys
sys.path.append("./")
from quaos.symplectic import PauliSum

def test_symplectic_matrix_single_pauli():
    pauli_list = ['X']
    weights = [1]
    sp = PauliSum(pauli_list, weights)
    expected_matrix = np.array([[1, 1, 0, 0]])
    np.testing.assert_array_equal(sp.symplectic_matrix(), expected_matrix)

def test_symplectic_matrix_multiple_paulis():
    pauli_list = ['X', 'Z', 'Y']
    weights = [1, 2, 3]
    sp = PauliSum(pauli_list, weights)
    expected_matrix = np.array([
        [1, 1, 0, 0],
        [2, 0, 1, 0],
        [3, 1, 1, 1]
    ])
    np.testing.assert_array_equal(sp.symplectic_matrix(), expected_matrix)

def test_symplectic_matrix_with_zero_weights():
    pauli_list = ['X', 'Z', 'Y']
    weights = [0, 0, 0]
    sp = PauliSum(pauli_list, weights)
    expected_matrix = np.array([
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 1, 1, 1]
    ])
    np.testing.assert_array_equal(sp.symplectic_matrix(), expected_matrix)

def test_symplectic_matrix_with_different_weights():
    pauli_list = ['X', 'Z', 'Y']
    weights = [1, 2, 3]
    sp = PauliSum(pauli_list, weights)
    expected_matrix = np.array([
        [1, 1, 0, 0],
        [2, 0, 1, 0],
        [3, 1, 1, 1]
    ])
    np.testing.assert_array_equal(sp.symplectic_matrix(), expected_matrix)

if __name__ == "__main__":
    pytest.main()