import pytest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)
from quaos.core.prime_Functions_Andrew import (
    gate
)


class TestGate:

    def test_gate(self):
        # TODO: Implement test_gate  
        assert True

    def test_name_string(self):
        # TODO: Implement test_name_string  
        assert True

    def test_copy(self):
        # TODO: Implement test_copy  
        assert True

    def test_print(self):
        # TODO: Implement test_print  
        assert True
