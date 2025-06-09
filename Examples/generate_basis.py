import sys
import numpy as np
sys.path.append("./")
from quaos.paulis import PauliString, PauliSum
from quaos.circuits import Circuit
from quaos.circuits.utils import solve_modular_linear
from quaos.circuits.known_circuits import to_x, to_ix
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
    
    linear_operations = []  # list of tuples of (pauli1, pauli2) of the the paulis multiplied together


    return c


def is_ix(pauli_string: PauliString, qudit: int | None = None) -> bool:
    if qudit is None:
        if np.any(pauli_string.z_exp != 0):
            return False
        if np.count_nonzero(pauli_string.x_exp) == 1:
            return True
        else:
            return False
    else:
        if pauli_string.z_exp[qudit] != 0:
            return False
        if pauli_string.x_exp[qudit] >= 1 and np.count_nonzero(pauli_string.x_exp) == 1:
            return True
        else:
            return False


def is_iz(pauli_string: PauliString, qudit: int | None = None) -> bool:
    if qudit is None:
        if np.any(pauli_string.x_exp != 0):
            return False
        if np.count_nonzero(pauli_string.z_exp) == 1:
            return True
        else:
            return False
    else:
        if pauli_string.x_exp[qudit] != 0:
            return False
        if pauli_string.z_exp[qudit] >= 1 and np.count_nonzero(pauli_string.z_exp) == 1:
            return True
        else:
            return False


def find_ix_iz(pauli_sum: PauliSum, qudit: int | None = None) -> tuple[list[int], list[int]]:
    ixs = []
    izs = []
    for i in range(pauli_sum.n_paulis()):
        if is_ix(pauli_sum.pauli_strings[i], qudit):
            ixs.append(i)
        elif is_iz(pauli_sum.pauli_strings[i], qudit):
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


def standard_form_to_basis(pauli_sum: PauliSum) -> tuple[PauliSum, list[tuple[int, int, int]]]:
    new_ps = pauli_sum.copy()
    multiplied_paulis = []

    for q in range(pauli_sum.n_qudits()):
    
        ixs, izs = find_ix_iz(new_ps, q)
        if len(ixs) == 0 and len(izs) == 0:
            print('No ix or iz found for qudit ', q)
            continue
        ixs = ixs[0] if len(ixs) > 0 else None
        izs = izs[0] if len(izs) > 0 else None
        if ixs is None:
            print('No ix found for qudit ', q)
        if izs is None:
            print('No iz found for qudit ', q)
        for p in range(new_ps.n_paulis()):  # min( , 2 * new_ps.n_qudits())
            if new_ps[p].x_exp[q] != 0 and ixs is not None and ixs != p:
                n = solve_modular_linear(new_ps[p].x_exp[q], new_ps[ixs].x_exp[q], new_ps.dimensions[q])
                new_ps[p] = new_ps[p] * new_ps[ixs]**n
                multiplied_paulis.append((p, ixs, n))
            if new_ps[p].z_exp[q] != 0 and izs is not None and izs != p:
                n = solve_modular_linear(new_ps[p].z_exp[q], new_ps[izs].z_exp[q], new_ps.dimensions[q])
                new_ps[p] = new_ps[p] * new_ps[izs]**n
                multiplied_paulis.append((p, izs, n))
    return new_ps, multiplied_paulis


def multiply_paulis(pauli_sum: PauliSum, multiplier_list: list[tuple[int, int, int]]) -> PauliSum:
    """
    Multiply the pauli strings in the pauli sum by the given multipliers.

    :param pauli_sum: The PauliSum to multiply.
    :param multiplier_list: A list of tuples (pauli_index, multiplier_index, multiplier_value) where
    the pauli at pauli_index is multiplied by the pauli at multiplier_index raised to the power of multiplier_value.

    :return: The PauliSum after multiplication.
    """
    new_pauli_sum = pauli_sum.copy()
    for p, m, n in multiplier_list:
        new_pauli_sum[p] = new_pauli_sum[p] * new_pauli_sum[m]**n
    return new_pauli_sum


def is_basis(pauli_sum: PauliSum) -> tuple[bool, list[int]]:
    ixs, izs = find_ix_iz(pauli_sum)
    if len(ixs) + len(izs) == 2 * pauli_sum.n_qudits():
        return True, ixs + izs
    elif len(ixs) + len(izs) > 2 * pauli_sum.n_qudits():
        # over complete - pick only the first 2 * n_qudits independent pauli strings
        ix_qudits = []
        ixs_ = []
        iz_qudits = []
        izs_ = []
        for ix in ixs:
            ixq = np.where(pauli_sum[ix].x_exp != 0)[0][0]
            if ixq not in ix_qudits:
                ix_qudits.append(ixq)
                ixs_.append(ix)
        for iz in izs:
            izq = np.where(pauli_sum[iz].z_exp != 0)[0][0]
            if izq not in iz_qudits:
                iz_qudits.append(izq)
                izs_.append(iz)
        return True, ixs_ + izs_
    else:
        # Check if the remaining pauli strings can be made up of the ixs and izs - if so it is an incomplete basis
        remaining = [i for i in range(pauli_sum.n_paulis()) if i not in ixs and i not in izs]
        missing_ix = []
        missing_iz = []
        for q in range(pauli_sum.n_qudits()):
            ixs_q, izs_q = find_ix_iz(pauli_sum, q)
            if len(ixs_q) == 0:
                missing_ix.append(q)
            if len(izs_q) == 0:
                missing_iz.append(q)
        for r in remaining:
            r_x_exp = pauli_sum[r].x_exp
            r_z_exp = pauli_sum[r].z_exp
            for q in missing_ix:
                if r_x_exp[q] != 0:
                    # If the pauli string has an x term on a qudit where we are missing an ix, it cannot be a basis
                    return False, []
            for q in missing_iz:
                if r_z_exp[q] != 0:
                    # If the pauli string has a z term on a qudit where we are missing an iz, it cannot be a basis
                    return False, []
        return True, ixs + izs


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

    print(h)
    output_pauli, multipliers = standard_form_to_basis(ham)
    output_pauli2 = multiply_paulis(h, multipliers)
    print(output_pauli)
    print(output_pauli2)
    print(output_pauli == output_pauli2)

    basis_check, basis_indices = is_basis(output_pauli)
    print("Is basis:", basis_check)
    print("Basis indices:", basis_indices)
    # c = to_basis(h_red, find_anticommuting_pairs(h_red))
