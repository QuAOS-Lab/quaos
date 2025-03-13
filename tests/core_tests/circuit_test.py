import unittest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)
from quaos.core.prime_Functions_Andrew import (
    circuit
)


class TestCircuit(unittest.TestCase):

    def test_circuit(self):
        # TODO: Implement test_circuit
        self.assertTrue(True)

    def test_length(self):
        # TODO: Implement test_length
        self.assertTrue(True)

    def test_unitary(self):
        # TODO: Implement test_unitary
        self.assertTrue(True)

    def test_add_gates_(self):
        # TODO: Implement test_add_gates_
        self.assertTrue(True)

    def test_insert_gates_(self):
        # TODO: Implement test_insert_gates_
        self.assertTrue(True)

    def test_delete_gates_(self):
        # TODO: Implement test_delete_gates_
        self.assertTrue(True)

    def test_copy(self):
        # TODO: Implement test_copy
        self.assertTrue(True)

    def test_print(self):
        # TODO: Implement test_print
        self.assertTrue(True)
