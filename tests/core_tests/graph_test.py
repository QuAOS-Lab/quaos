import unittest
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


class TestGraph(unittest.TestCase):

    def test_graph(self):
        # TODO: Implement test_gate  
        self.assertTrue(True)

    def test_add_vertex_(self):
        # TODO: Implement test_add_vertex_
        self.assertTrue(True)

    def test_lade_vertex_(self):
        # TODO: Implement test_lade_vertex_
        self.assertTrue(True)

    def test_lade_edge_(self):
        # TODO: Implement test_lade_edge_
        self.assertTrue(True)

    def test_neighbors(self):
        # TODO: Implement test_neighbors
        self.assertTrue(True)

    def test_edges(self):
        # TODO: Implement test_edges
        self.assertTrue(True)

    def test_clique(self):
        # TODO: Implement test_clique
        self.assertTrue(True)

    def test_degree(self):
        # TODO: Implement test_degree
        self.assertTrue(True)

    def test_ord(self):
        # TODO: Implement test_ord
        self.assertTrue(True)

    def test_print(self):
        # TODO: Implement test_print
        self.assertTrue(True)

    def test_print_neighbors(self):
        # TODO: Implement test_print_neighbors
        self.assertTrue(True)

    def test_copy(self):
        # TODO: Implement test_copy
        self.assertTrue(True)
