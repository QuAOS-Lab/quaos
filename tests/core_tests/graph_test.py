import pytest
import sys
import os
import numpy as np

test_dir = os.path.dirname(os.path.abspath(__file__))
module_dir = os.path.join(test_dir, "..")
sys.path.append(module_dir)
sys.path.append(test_dir)
from quaos.core.prime_Functions_Andrew import (
    graph
)


class TestGraph:

    def test_graph(self):
        # TODO: Implement test_gate  
        assert True

    def test_add_vertex_(self):
        # TODO: Implement test_add_vertex_
        assert True

    def test_lade_vertex_(self):
        # TODO: Implement test_lade_vertex_
        assert True

    def test_lade_edge_(self):
        # TODO: Implement test_lade_edge_
        assert True

    def test_neighbors(self):
        # TODO: Implement test_neighbors
        assert True

    def test_edges(self):
        # TODO: Implement test_edges
        assert True

    def test_clique(self):
        # TODO: Implement test_clique
        assert True

    def test_degree(self):
        # TODO: Implement test_degree
        assert True

    def test_ord(self):
        # TODO: Implement test_ord
        assert True

    def test_print(self):
        # TODO: Implement test_print
        assert True

    def test_print_neighbors(self):
        # TODO: Implement test_print_neighbors
        assert True

    def test_copy(self):
        # TODO: Implement test_copy
        assert True
