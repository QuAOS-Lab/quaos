import sys
sys.path.append("./")
from quaos.paulis import PauliString, PauliSum
from quaos.circuits import Circuit
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
    c = Circuit(dimensions=pauli_sum.dimensions)
    ignore = []
    for qudit_number, pair in enumerate(anticommuting_pairs):
        print(qudit_number, pair)

        c_temp = to_ix(ps[pair[0]], qudit_number, ignore)
        ps = c_temp.act(ps)
        c += c_temp
        c_temp = to_ix(ps[pair[1]], qudit_number, ignore)
        ps = c_temp.act(ps)
        c += c_temp
        print(ps)
        ignore.append(qudit_number)

    anticommuting = [x for tup in pairs for x in tup]
    remaining = [p for p in range(pauli_sum.n_paulis()) if p not in anticommuting]

    for qudit_number in range(pauli_sum.n_paulis()):
        if len(remaining) > 0:
            if qudit_number not in ignore and len(ignore) < pauli_sum.n_qudits() - 2:
                c_temp = to_ix(ps[remaining[0]], qudit_number, ignore)
                ps = c_temp.act(ps)
                c += c_temp
                remaining.pop(0)
    
    return c

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

def basis_anticommuting_iteration(pauli_sum: PauliSum, anticommuting_pair: tuple[int, int]) -> Circuit:
    # First we get the form
    # Z IIII
    # X IIII
    # whatever
    ps = pauli_sum.copy()
    c = to_ix(ps[anticommuting_pair[0]], 0)
    ps = c.act(ps)
    c = to_ix(ps[anticommuting_pair[1]], 0)
    ps = c.act(ps)
    # Now we want
    # X IIII
    # Z IIII
    # I ...
    # I ...

    return c


if __name__ == "__main__":

    dims = [2, 2, 2, 2]
    ham = random_pauli_hamiltonian(4, dims, mode='uniform')
    # ham = PauliSum(['x0z0 x0z1 x1z0 x1z1', 'x0z1 x0z1 x1z0 x0z1',
    #                 'x1z0 x0z1 x1z0 x1z0', 'x0z1 x0z0 x0z0 x1z1'], dimensions=dims)
    print(ham)

    print(ham)
    pairs = find_anticommuting_pairs(ham)
    order = [x for tup in pairs for x in tup]

    # ham.reorder(order)
    c = to_basis(ham, pairs)
    print(ham)
    print(c.act(ham))
    print(order)
    print(pairs)
