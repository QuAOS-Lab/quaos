import numpy as np
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Circuit
from collections import defaultdict


def int_to_bases(a, dims):
    dims = np.flip(dims)
    aa = [a % dims[0]]
    for i in range(1, len(dims)):
        s0 = aa[0] + sum([aa[i1] * dims[i1 - 1] for i1 in range(1, i)])
        s1 = np.prod(dims[:i])
        aa.append(((a - s0) // s1) % dims[i])
    dims = np.flip(dims)
    return np.flip(np.array(aa))


def group_indices(lst):
    """
    Groups indices of the same value in a list into sub-lists.
    For example, if the input list is [1, 2, 1, 3, 2], the output will be [[0, 2], [1, 4], [3]].
    """
    index_dict = defaultdict(list)
    for idx, value in enumerate(lst):
        index_dict[value].append(idx)
    return [indices for indices in index_dict.values()]


def Hadamard_Symmetric_PauliSum(n_paulis, n_qubits, n_sym_q, q_print=False, random_coefficients=False, seed=None):
    # create coefficients
    if seed is not None:
        np.random.seed(seed)
    c_int_bands = np.sort(np.random.randint(n_paulis, size=n_paulis))
    c_bands = group_indices(c_int_bands)

    if random_coefficients:
        coefficients = np.zeros(n_paulis)
        sym_bands = []
        for i, b in enumerate(c_bands):
            coefficients[b] = np.random.normal(0, 1)
            if len(b) != 1:
                sym_bands.append(b)
    else:
        coefficients = np.ones(n_paulis)
        sym_bands = []
        for i, b in enumerate(c_bands):
            if len(b) != 1:
                sym_bands.append(b)

    n_extra_bands = np.sum([np.floor(len(b) / 2) for b in sym_bands]) - n_sym_q
    pauli_strings = ['' for i in range(n_paulis)]

    all_x = []
    all_z = []
    for i in range(n_sym_q):
        x_pauli = []
        z_pauli = []
        if len(sym_bands) >= 1:
            b_ind = np.random.randint(len(sym_bands))
            b = sym_bands[b_ind]
            x_ind = np.random.randint(len(b))
            x_pauli.append(b[x_ind])
            b.pop(x_ind)
            z_ind = np.random.randint(len(b))
            z_pauli.append(b[z_ind])
            b.pop(z_ind)
            if len(b) < 2:
                sym_bands.pop(b_ind)
            else:
                sym_bands[b_ind] = b

            # randomly adding extra x and zs if possible
            if n_extra_bands > 0 and len(sym_bands) >= 1:
                extras = np.random.randint(n_extra_bands)
                n_extra_bands -= extras

                for j in range(extras):
                    b_ind = np.random.randint(len(sym_bands))
                    b = sym_bands[b_ind]
                    x_ind = np.random.randint(len(b))
                    x_pauli.append(b[x_ind])
                    b.pop(x_ind)
                    z_ind = np.random.randint(len(b))
                    z_pauli.append(b[z_ind])
                    b.pop(z_ind)
                    if len(b) < 2:
                        sym_bands.pop(b_ind)
                    else:
                        sym_bands[b_ind] = b
        if q_print:
            print(coefficients[x_pauli])
        for j in range(n_paulis):
            if j in x_pauli:
                pauli_strings[j] += 'x1z0 '
            elif j in z_pauli:
                pauli_strings[j] += 'x0z1 '
            else:
                pauli_strings[j] += 'x0z0 '
        all_x += x_pauli
        all_z += z_pauli
    if q_print:
        print(all_x, all_z)
    non_sym_paulis = [i for i in range(n_paulis) if i not in all_x and i not in all_z]
    q_dims = [2 for i in range(2 * (n_qubits - n_sym_q))]
    available_paulis = list(np.arange(int(np.prod(q_dims))))
    for i, x in enumerate(all_x):
        pauli_index = np.random.choice(available_paulis)
        available_paulis.remove(pauli_index)
        exponents = int_to_bases(pauli_index, q_dims)
        for j in range(n_qubits - n_sym_q):
            r, s = int(exponents[2 * j]), int(exponents[2 * j + 1])
            pauli_strings[x] += f"x{r}z{s} "
            pauli_strings[all_z[i]] += f"x{r}z{s} "

        pauli_strings[x].strip()
        pauli_strings[all_z[i]].strip()

    for i in non_sym_paulis:
        pauli_index = np.random.choice(available_paulis)
        available_paulis.remove(pauli_index)
        exponents = int_to_bases(pauli_index, q_dims)
        for j in range(n_qubits - n_sym_q):
            r, s = int(exponents[2 * j]), int(exponents[2 * j + 1])
            pauli_strings[i] += f"x{r}z{s} "
        pauli_strings[i].strip()

    P = PauliSum.from_string(pauli_strings, weights=coefficients, dimensions=[2] * n_qubits)

    # construct random Clifford circuit
    n_qubits = P.n_qudits()
    if n_qubits < 2:
        return P, None

    C = Circuit.from_random(n_gates=100, dimensions=P.dimensions)

    phases = P.phases
    weights = P.weights
    weights *= np.array([-1] * n_paulis) ** phases
    P = PauliSum.from_string(pauli_strings, weights=weights, dimensions=P.dimensions)

    return P, C


def SWAP_symmetric_PauliSum(n_paulis, n_qubits):
    ps = PauliSum.from_random(n_paulis, dimensions=[2] * n_qubits, rand_weights=False)
    ps[:, 1] = ps[:, 0]  # make qudit 2 equal to qudit 1

    C = Circuit.from_random(n_gates=100, dimensions=ps.dimensions)

    return C.act(ps)


# def random_clifford(n_qubits, depth=100) -> Circuit:

#     C = Circuit(dimensions=[2 for i in range(n_qubits)])
#     gate_list = [H, S, CX]
#     gg = []
#     for i in range(depth):
#         g_i = np.random.randint(3)
#         if g_i == 2:
#             aa = list(random.sample(range(n_qubits), 2))
#             gg += [gate_list[g_i](aa[0], aa[1], 2)]
#         else:
#             aa = list(random.sample(range(n_qubits), 1))
#             gg += [gate_list[g_i](aa[0], 2)]

#     C.add_gate(gg)
#     return C


if __name__ == "__main__":
    # Test the Hadamard_Symmetric_PauliSum function
    P, C = Hadamard_Symmetric_PauliSum(9, 4, 2)
    print(P)
