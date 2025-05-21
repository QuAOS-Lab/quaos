import sys
sys.path.append("./")
from quaos.paulis import PauliString, PauliSum
from quaos.circuits import Circuit, Gates
from quaos.known_circuits import to_x, to_ix
from quaos.hamiltonian import random_pauli_hamiltonian


def find_anticommuting_pairs(pauli_sum: PauliSum) -> list[tuple[int, int]]:
    """
    Pairs the pauli strings that anticommute
    """
    ps = pauli_sum.copy()
    spm = ps.symplectic_product_matrix()
    anticommuting_pairs = []
    used_paulis = []
    for i in range(pauli_sum.n_paulis()):
        if i not in used_paulis:
            for j in range(pauli_sum.n_paulis()):
                if i != j and j not in used_paulis and i not in used_paulis:
                    if spm[i, j] == 1:
                        anticommuting_pairs.append((i, j))
                        used_paulis.append(i)
                        used_paulis.append(j)
                        

    return anticommuting_pairs


def to_basis(pauli_sum: PauliSum, anticommuting_pairs: list[tuple[int, int]]) -> Circuit:
    
    # now we loop through the pairs making them XIIII, ZIIII, IXIIII, IZIIII, IIIXII, IIIZII, ....

    if len(anticommuting_pairs) > 2 * pauli_sum.n_qudits():
        anticommuting_pairs = anticommuting_pairs[0:2 * pauli_sum.n_qudits()]

    ps = pauli_sum.copy()
    c = Circuit(dimensions = pauli_sum.dimensions)
    for qudit_number, pair in enumerate(anticommuting_pairs):
        pass



def pauli_sum_to_basis(pauli_sum: PauliSum) -> Circuit:
    """
    Puts the PauliSum as close to the form

    XIIIII...
    ZIIIII...
    IXIIII...
    IZIIII...
    IIIXII...
    IIIZII...
    .
    . 
    . 

    as possible
    """
    pairs = find_anticommuting_pairs(pauli_sum)
    if len(pairs) >= 2 * pauli_sum.n_qudits():
        commuting_basis_pairs = 0
    else:
        commuting_basis_pairs = 2 * pauli_sum.n_qudits() - len(pairs)

    if commuting_basis_pairs == 0:
        return to_complete_basis(pauli_sum, pairs)



if __name__ == "__main__":

    dims = [2, 2, 2, 2]
    ham = random_pauli_hamiltonian(4, dims, mode='uniform')

    print(ham)
    pairs = find_anticommuting_pairs(ham)
    order = [x for tup in pairs for x in tup]
    print(ham.weights)
    print(order)
    print(pairs)
    ham.reorder(order)
    print(ham.weights)
    print(ham)
    print(order)

    print(pairs)
