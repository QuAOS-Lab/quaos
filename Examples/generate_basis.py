import sys
import numpy as np
sys.path.append("./")
from quaos.paulis import PauliString, PauliSum
from quaos.circuits import Circuit
from quaos.circuits.utils import solve_modular_linear
from quaos.known_circuits import to_x, to_ix
from quaos.hamiltonian import random_pauli_hamiltonian, pauli_reduce


def find_anticommuting_pairs(pauli_sum: PauliSum) -> list[tuple[int, int]]:
    """
    Pairs the pauli strings that anticommute
    """
    ps = pauli_sum.copy()
    spm = ps.symplectic_product_matrix()
    anticommuting_pairs = []
    used_paulis = []
    for i in range(pauli_sum.n_paulis()):
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

    n_pairs = len(anticommuting_pairs)
    n_q = pauli_sum.n_qudits()
    n_p = pauli_sum.n_paulis()
    n_unpaired_paulis = n_p - 2 * n_pairs
    paired_paulis = [x for tup in anticommuting_pairs for x in tup]
    remaining = [x for x in range(n_p) if x not in paired_paulis]

    for qudit_number, pair in enumerate(anticommuting_pairs):
        if qudit_number <= n_q:
            c_temp = to_ix(ps[pair[0]], qudit_number)
            ps = c_temp.act(ps)
            c += c_temp
            c_temp = to_ix(ps[pair[1]], qudit_number)
            ps = c_temp.act(ps)
            c += c_temp

    i = 0
    print(remaining)
    for qudit_number in range(n_pairs, n_pairs + n_unpaired_paulis):
        if qudit_number <= n_q:
            c_temp = to_ix(ps[remaining[i]], qudit_number)
            i += 1
            ps = c_temp.act(ps)
            c += c_temp
    
    return c


def is_ix(pauli_string: PauliString) -> bool:
    if np.any(pauli_string.z_exp != 0):
        return False
    if np.count_nonzero(pauli_string.x_exp) == 1:
        return True
    else:
        return False


def is_iz(pauli_string: PauliString) -> bool:
    if np.any(pauli_string.x_exp != 0):
        return False
    if np.count_nonzero(pauli_string.z_exp) == 1:
        return True
    else:
        return False
    

def find_ix_iz(pauli_sum: PauliSum) -> tuple(list[int], list[int]):
    ixs = []
    izs = []
    for i in range(pauli_sum.n_paulis()):
        if is_ix(pauli_sum.pauli_strings[i]):
            ixs.append(i)
        elif is_iz(pauli_sum.pauli_strings[i]):
            izs.append(i)
    return ixs, izs


def use_ix_remove_x(pauli_sum: PauliSum, ixs: list[int]):
    new_ps = pauli_sum.copy()
    multiplied_paulis = []
    for ix in ixs:  # the ixth string is the ix - multiply others by this to remove their x terms on x_qudit
        x_qudit = np.where(new_ps[ix].x_exp != 0)[0][0]
        x_exp = new_ps[ix].x_exp[x_qudit]
        for i in range(pauli_sum.n_paulis()):
            if i != ix:
                if pauli_sum[i].x_exp[x_qudit] != 0:
                    n = solve_modular_linear(pauli_sum[i].x_exp[x_qudit], x_exp, pauli_sum.dimensions[x_qudit])
                    new_ps[i] = pauli_sum[i] * pauli_sum[ix]**n
                multiplied_paulis.append(())



def basis_from_standard_form(pauli_sum: PauliSum):
    ixs, izs = find_ix_iz(pauli_sum)

    if len(ixs) > 0:
        use_ix_remove_x(pauli_sum, ixs)



if __name__ == "__main__":
    from quaos.paulis import commutation_graph
    import matplotlib.pyplot as plt
    n_qudits = 7
    dims = [2] * n_qudits
    n_paulis = 10
    ham = random_pauli_hamiltonian(n_paulis, dims, mode='uniform')

    h_red, conditioned_hamiltonians, C, all_phases = pauli_reduce(ham)

    # choose a symmetry subsector and simplify
    h = conditioned_hamiltonians[0]
    # anticommuting_pairs = find_anticommuting_pairs(h)
    # print(anticommuting_pairs)
    # print(h)
    # c = to_basis(h, anticommuting_pairs)
    # print(c.act(h))
    h_red2, conditioned_hamiltonians2, C2, all_phases2 = pauli_reduce(h)

    print(ham)
    print(h)
    print(conditioned_hamiltonians2[0])

    # print(h_red == h_red2)

    # c = to_basis(h_red, find_anticommuting_pairs(h_red))


    plt.show()
