import sys
sys.path.append("./")

from quaos.hamiltonian import symplectic_pauli_reduction, pauli_reduce
from quaos.paulis import PauliSum


ham = ['x1z0 x1z0 x0z0', 'x0z0 x1z0 x1z0', 'x0z1 x0z0 x0z1']
ham = PauliSum(ham, weights=[1, 1, 1], dimensions=[2, 2, 2])
print(ham, '/n')
circuit = symplectic_pauli_reduction(ham)
h_reduced, conditioned_hams, reducing_circuit, eigenvalues = pauli_reduce(ham)
print(h_reduced)
for h in conditioned_hams:
    print(h)

# random hamiltonian example
# n_qudits = 4
# n_paulis = 4
# dimension = 2
# ham = random_pauli_hamiltonian(n_paulis, [dimension] * n_qudits, mode='uniform')
# circuit = symplectic_pauli_reduction(ham)
# print(ham)
# h_reduced, conditioned_hams, reducing_circuit, eigenvalues = pauli_reduce(ham)
# print(h_reduced)
# for h in conditioned_hams:
#     print(h)
