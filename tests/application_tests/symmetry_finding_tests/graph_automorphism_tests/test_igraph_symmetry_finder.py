import numpy as np
import pytest

igraph = pytest.importorskip("igraph")

from sympleq.applications.symmetry_finding.graph_automorphism import (
    find_igraph_clifford_symmetries,
)
from sympleq.core.circuits.gates import GATES
from sympleq.models.random_hamiltonian import random_gate_symmetric_hamiltonian


def test_igraph_symmetry_finder_returns_valid_clifford_symmetry():
    pauli_sum = random_gate_symmetric_hamiltonian(
        GATES.SWAP,
        dimension=2,
        qudit_indices=(0, 1),
        n_qudits=2,
        n_paulis=4,
        scrambled=False,
    )
    pauli_sum.weight_to_phase()

    symmetries, checked = find_igraph_clifford_symmetries(pauli_sum, num_symmetries=1)

    assert checked > 0
    assert len(symmetries) == 1

    out = symmetries[0].act(pauli_sum, (0, 1)).to_standard_form()
    ref = pauli_sum.to_standard_form()
    out.weight_to_phase()
    ref.weight_to_phase()

    assert np.array_equal(out.tableau, ref.tableau)
    assert np.array_equal(out.phases, ref.phases)
    assert np.allclose(out.weights, ref.weights)
