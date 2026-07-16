import pytest
from sympleq.core.circuits.gates import GATES
from sympleq.core.symmetries.conditional_hamiltonian import ConditionalHamiltonian
from sympleq.models.random_hamiltonian import random_gate_symmetric_hamiltonian
from sympleq.core.symmetries.pauli import pauli_reduce
from sympleq.core.symmetries.clifford import min_qudit_clifford_symmetry
from sympleq.core.circuits.gate_decomposition_to_circuit import gate_to_circuit
import numpy as np
from numpy.random import default_rng

rng = default_rng()


class TestConditionalHamiltonianFinder:

    def generate_symmetry(self, n_qudits, n_paulis):
        P_sym = random_gate_symmetric_hamiltonian(GATES.H, dimension=2, qudit_indices=tuple([0]),
                                                  n_paulis=n_paulis, n_qudits=n_qudits)
        conditional_hamiltonian = pauli_reduce(P_sym)
        h_red = conditional_hamiltonian.original_hamiltonian
        F, Sy, T = min_qudit_clifford_symmetry(P_sym)
        C_F = gate_to_circuit(F, dimensions=[2 for i in range(n_qudits)])
        symmetrised = T.inverse().act(h_red, tuple(np.arange(n_qudits)))
        assert C_F.act(h_red).is_close(h_red, literal=False)
        assert Sy.act(symmetrised, tuple(np.arange(n_qudits))).is_close(symmetrised, literal=False)
        return h_red, Sy, T

    def test_conditional_hamiltonians(self):
        n_test = 10
        n_qudits = 3
        n_paulis = 16
        for _ in range(n_test):
            Ham, sym, tra = self.generate_symmetry(n_qudits=n_qudits, n_paulis=n_paulis)
            cond_hamiltonian = ConditionalHamiltonian(Ham, sym, tra)
            assert cond_hamiltonian.test_conditional_hamiltonian()
