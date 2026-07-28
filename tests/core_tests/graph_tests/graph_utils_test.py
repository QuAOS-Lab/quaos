import numpy as np
import pytest

from sympleq.core.graphs.utils import qudit_coupling_graph


class TestGraphUtils:
    def test_qudit_coupling_graph_identity_has_no_edges(self):
        graph = qudit_coupling_graph(np.eye(4, dtype=int))

        assert graph == {0: set(), 1: set()}

    def test_qudit_coupling_graph_detects_two_qudit_coupling(self):
        F = np.eye(4, dtype=int)
        F[1, 0] = 1

        graph = qudit_coupling_graph(F)

        assert graph == {0: {1}, 1: {0}}

    def test_qudit_coupling_graph_rejects_invalid_shapes(self):
        with pytest.raises(ValueError, match="square"):
            qudit_coupling_graph(np.ones((2, 3), dtype=int))

        with pytest.raises(ValueError, match="even dimension"):
            qudit_coupling_graph(np.eye(3, dtype=int))
