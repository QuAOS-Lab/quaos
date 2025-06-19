import numpy as np
import scipy
import scipy.sparse
import scipy.sparse.linalg
import itertools
import random
import math
import networkx as nx
from functools import reduce
from operator import mul
from .pauli import (
    pauli, string_to_pauli, pauli_to_matrix, pauli_to_string,
    symplectic_inner_product, quditwise_inner_product
)

# NAMING CONVENTIONS

# a - int
# aa - list{int}
# aaa - list{list{int}}
# b - bool
# c - float (constant)
# r - random
# i - indexing
# ss - str
# m - matrix (numpy.array or scipy.sparse.csr_matrix)
# P - Pauli
# C - circuit
# G - gate
# A - graph (adjacency matrix)
# p - number of paulis
# q - number of qudits
# d - dimensions
# functions or methods which end with underscores modify the inputs


def I_mat(d):
    """Return the d-dimensional Identity matrix."""
    return scipy.sparse.csr_matrix(np.diag([1] * d))


def H_mat(d):
    """Return the d-dimensional Hadarmard matrix."""
    omega = math.e**(2 * math.pi * 1j / d)
    H_mat = np.array([[omega**(i0 * i1) for i0 in range(d)] for i1 in range(d)])
    H_mat = 1 / np.sqrt(d) * H_mat
    return scipy.sparse.csr_matrix(H_mat)


def S_mat(d):
    """Return the d-dimensional S Clifford matrix."""
    if d == 2:
        return scipy.sparse.csr_matrix(np.diag([1, 1j]))

    omega = math.e**(2 * math.pi * 1j / d)
    S_mat = np.diag([omega ** (i * (i - 1) / 2) for i in range(d)])
    return scipy.sparse.csr_matrix(S_mat)


# READING, WRITING, AND MISCELLANEOUS
def read_Luca_Test(path, dims=2):
    """Reads a Hamiltonian file and parses the Pauli strings and coefficients.

    Args:
        path (str): Path to the Hamiltonian file.

    Returns:
        tuple: Number of Paulis in Hamiltonian, number of qudits in Hamiltonian,
            set of Paulis in Hamiltonian, and coefficients in Hamiltonian.
    """
    f = open(path, "r")
    sss0 = f.readlines()
    f.close()
    sss1 = []
    cc = []
    for i0 in range(len(sss0)):
        pauli_list = sss0[i0].split(', {')
        coefficient = pauli_list[0][1:].replace(" ", "").replace('*I', 'j')
        # print(coefficient)
        # print(coefficient[0])
        cc.append(np.complex(coefficient))
        ss = ''
        for i1 in range(1, len(pauli_list)):
            ss += 'x' + str(pauli_list[i1].count('X'))
            ss += 'z' + str(pauli_list[i1].count('Z'))
            ss += ' '
        sss1.append(ss[:-1])
    # if type(dims) == int:
    #     dims = np.array([dims]*(len(sss0[0].split(', {'))-1))
    # elif type(dims) == list:
    #     dims = np.array(dims)

    # cc1 = []
    # for i0 in range(len(sss1)):
    #     pauli_list = sss1[i0].split(' ')
    #     pauli_list_conj = ''
    #     for i1 in range(len(dims)):
    #         ss = pauli_list[i1]
    #         ss_conj = 'x'+str((-int(ss[1]))%int(dims[i1]))+'z'+str((-int(ss[3]))%int(dims[i1]))+' '
    #         pauli_list_conj += ss_conj
    #     cc1.append(1/2*cc[i0]+1/2*np.conj(cc[sss1.index(pauli_list_conj[:-1])]))

    # # cc1 = cc.copy()

    return string_to_pauli(sss1, dims), cc


def loading_bar(runs, length=50, scalings=[]):
    """
    Call within a loop to produce a loading bar which can be scaled for various timings.

    Args:
        runs (list[tuple[int, int]]): List of pairs: current iteration, total iterations.
        length (int): Length of bar.
        scaling (list[function]): Scale speed of progress.
    """
    a0 = len(runs)
    scalings += [lambda x:x] * (a0 - len(scalings))
    a1 = scalings[0](runs[0][0]) + sum(scalings[i](runs[i][0]) * scalings[i - 1](runs[i - 1][1]) for i in range(1, a0))
    a2 = reduce(mul, [scalings[i](runs[i][1]) for i in range(a0)])
    # for scale in scalings[::-1]:
    #     a1,a2 = scale(a1),scale(a2)
    ss0 = ("{0:.1f}").format(int(1000 * a1 / a2) / 10)
    ss1 = ' '.join(str(runs[i][0]) + '/' + str(runs[i][1]) for i in range(a0))
    bar = '█' * int(length * a1 / a2) + '-' * (length - int(length * a1 / a2))
    print(f' |{bar}| {ss0}% {ss1}', end="\r")
    if runs[0][0] >= runs[0][1] - 1:
        print(" " * (length + 6 + len(ss0) + len(ss1)), end="\r")


# GATES & CIRCUITS
class gate:
    """A class for storing quantum gates as a name and a list of qudits.

    Attributes:
        name (function): Function representing the action of the gate.
        aa (list[int]): List of qudits acted upon.
    """
    def __init__(self, name, aa):
        """
        Args:
            name (function): Function representing the action of the gate.
            aa (list[int]): List of qudits acted upon.
        """
        self.name = name
        self.aa = aa

    def name_string(self):
        """Returns the name of the gate as a string.

        Returns:
            str: Name of function as a string.
        """
        return self.name.__name__

    def copy(self):
        """
        Deep copy of self.

        Returns:
            gate: A deep copy of gate.
        """
        return gate(self.name, [a for a in self.aa])

    def print(self):
        """Print self as a name and a tuple of qudits.

        Prints the name of the gate and a tuple of the qudits it acts on.
        """
        print("%s(%s)" % (self.name_string(), str(self.aa).replace(' ', '')[1:-1]))


class circuit:
    """
    A class for storing quantum circuits as a collection of gates along with the dimension.

    Attributes:
        dims (list[int]): Dimensions of the qudits in the circuit.
        gg (list[gate]): List of gates in the circuit.

    Methods:
        length(): Returns the number of gates in the circuit.
        unitary(): Converts the circuit to its corresponding unitary matrix.
        add_gates_(C): Appends gates to the end of the circuit.
    """
    def __init__(self, dims):
        """A class for storing quantum circuits as a collection of gates along with the dimension.

        Args:
            qudits (int): Number of qudits in circuit.
        """
        self.dims = dims
        self.gg = []

    def length(self):
        """Returns the number of gates in the circuit.

        Returns:
            int: The number of gates in the circuit.
        """
        return len(self.gg)

    def unitary(self):
        """
        Convert self to its corresponding unitary matrix.

        Returns:
            scipy.sparse.csr_matrix: Unitary matrix representation of self.
        """
        q = len(self.dims)
        m = scipy.sparse.csr_matrix(([1] * (np.prod(self.dims)), (range(np.prod(self.dims)), range(np.prod(self.dims)))))
        for g in self.gg:
            m = globals()[g.name_string() + '_unitary'](g.aa, self.dims) @ m
        return m

    def add_gates_(self, C):
        """Append gates to the end of self.

        Args:
            C (circuit or list[gate] or gate): Gates to be appended to self.
        """
        if type(C) is circuit:
            self.gg += C.gg
        elif type(C) is gate:
            self.gg.append(C)
        else:
            self.gg += C

    def insert_gates_(self, C, a):
        """
        Insert gates at a given timestep in self.

        Args:
            C (circuit or list[gate] or gate): Gates to be inserted into self.
            a (int): Index for insertion.
        """
        if type(C) is gate:
            self.gg.insert(a, C)
        elif type(C) is circuit:
            self.gg[a:a] = C.gg
        else:
            self.gg[a:a] = C

    def delete_gates_(self, aa):
        """
        Delete gates at specific timesteps.

        Args:
            aa (int or list[int]): Indices where gates should be deleted.
        """
        if type(aa) is int:
            del self.gg[aa]
        else:
            self.gg = [self.gg[i] for i in range(length(self.dims)) if i not in aa]

    def copy(self):
        """Deep copy of self.

        Returns:
            circuit: A deep copy of self.
        """
        return circuit(self.dims, [g.copy() for g in self.gg])

    def print(self):
        """Print self as gates on consecutive lines."""
        for g in self.gg:
            g.print()


def act(P, C):
    """Returns the outcome of a given circuit (or gate or list of gates) acting on a given pauli.

    Args:
        P (pauli): Pauli to be acted upon
        C (gate or circuit or list[gate]): gates to act on Pauli

    Returns:
        pauli: result of C acting on P by conjugation
    """
    if P is None:
        return P
    elif type(C) is gate:
        return C.name(P, C.aa)
    elif type(C) is circuit:
        return act(P, C.gg)
    elif len(C) == 0:
        return P
    elif len(C) == 1:
        return act(P, C[0])
    else:
        return act(act(P, C[0]), C[1:])


def H(P, aa):
    """
    Function for the gate representation of Hadamard gate. Transformation:
    X: Z, Z: -X

    Args:
        P (pauli): Pauli to be acted upon.
        aa (list[int]): Qudits to be acted upon.

    Returns:
        pauli: Result of H(aa) acting on P.
    """
    Q = P.copy()
    X, Z = Q.X, Q.Z
    for a in aa:
        for i in range(Q.paulis()):
            Q.phases[i] += -1 * X[i, a] * Z[i, a] * P.lcm // P.dims[a]
        X[:, a], Z[:, a] = -Z[:, a].copy(), X[:, a].copy()
    return pauli(X, Z, Q.dims, Q.phases)


def S(P, aa):
    """
    Function for the gate representation of phase gate. Phase gate transformation: X: XZ, Z: Z

    Args:
        P (pauli): Pauli to be acted upon.
        aa (list[int]): Qudits to be acted upon.

    Returns:
        pauli: Result of S(aa) acting on P.
    """
    Q = P.copy()
    X, Z = Q.X, Q.Z
    for a in aa:
        for i in range(Q.paulis()):
            if P.dims[a] == 2:
                Q.phases[i] += X[i, a] * Z[i, a] * P.lcm // P.dims[a]
            else:
                Q.phases[i] += math.comb(X[i, a], 2) * P.lcm // P.dims[a]
        Z[:, a] += X[:, a]
    return pauli(X, Z, Q.dims, Q.phases)


def CX(P, aa):
    """
    Function for the gate representation of CNOT gate. Transformation: X*I: X*X I*X: I*X Z*I: Z*I I*Z: -Z*Z

    Args:
        P (pauli): Pauli to be acted upon.
        aa (list[int]): Control aa[0] and target aa[1].

    Returns:
        pauli: Result of CNOT(aa[0],aa[1]) acting on P.
    """
    a0, a1 = aa[0], aa[1]
    if P.dims[a0] != P.dims[a1]:
        raise Exception("Entangling gates must be between two qudits of equal dimensions")
    Q = P.copy()
    X, Z = Q.X, Q.Z
    for i in range(Q.paulis()):
        if P.dims[a0] == 2:
            Q.phases[i] += X[i, a0] * Z[i, a1] * (X[i, a1] + Z[i, a0] + 1) * P.lcm // P.dims[a0]
    X[:, a1] += X[:, a0]
    Z[:, a0] -= Z[:, a1]
    return pauli(X, Z, Q.dims, Q.phases)


def SWAP(P, aa):
    """Function for the gate representation of SWAP gate. Gate transformation:
    X*I -> I*X, I*X -> X*I, Z*I -> I*Z, I*Z -> Z*I

    Args:
        P (pauli): Pauli to be acted upon.
        aa (list[int]): Targets aa[0] and aa[1].

    Returns:
        pauli: result of SWAP(aa[0],aa[1]) acting on P.

    """
    X, Z = P.X, P.Z
    a0, a1 = aa[0], aa[1]
    X[:, a0], X[:, a1] = X[:, a1].copy(), X[:, a0].copy()
    Z[:, a0], Z[:, a1] = Z[:, a1].copy(), Z[:, a0].copy()
    return pauli(X, Z)


def H_unitary(aa, dims):
    """Function for mapping Hadamard gate to corresponding unitary matrix

    Args:
        aa (list{int}): indices for Hadamard tensors
        dims (list{int}): dimensions of the qudits

    Returns:
        scipy.sparse.csr_matrix: matrix representation of q-dimensional H(aa)
    """
    return tensor([H_mat(dims[i]) if i in aa else I_mat(dims[i]) for i in range(len(dims))])


def S_unitary(aa, dims):
    """
    Function for mapping phase gate to corresponding unitary matrix.

    Args:
        aa (list[int]): indices for phase gate tensors
        dims (list[int]): dimensions of the qudits

    Returns:
        scipy.sparse.csr_matrix: matrix representation of q-dimensional S(aa)
    """
    return tensor([S_mat(dims[i]) if i in aa else I_mat(dims[i]) for i in range(len(dims))])


def bases_to_int(aa, dims):
    dims = np.flip(dims)
    aa = np.flip(aa)
    a = aa[0] + sum([aa[i1] * np.prod(dims[:i1]) for i1 in range(1, len(dims))])
    dims = np.flip(dims)
    aa = np.flip(aa)
    return a


def int_to_bases(a, dims):
    dims = np.flip(dims)
    aa = [a % dims[0]]
    for i in range(1, len(dims)):
        s0 = aa[0] + sum([aa[i1] * dims[i1 - 1] for i1 in range(1, i)])
        s1 = np.prod(dims[:i])
        aa.append(((a - s0) // s1) % dims[i])
    dims = np.flip(dims)
    return np.flip(np.array(aa))


def CX_func(i, a0, a1, dims):
    aa = int_to_bases(i, dims)
    aa[a1] = (aa[a1] + aa[a0]) % dims[a1]
    return bases_to_int(aa, dims)


def CX_unitary(aa, dims):
    """
    Function for mapping CNOT gate to corresponding unitary matrix.

    Args:
        aa (list[int]): Control aa[0] and target aa[1].
        dims (list[int]): Dimensions of the qudits.

    Returns:
        scipy.sparse.csr_matrix: Matrix representation of q-dimensional CNOT(aa[0], aa[1]).
    """
    _ = len(dims)
    D = np.prod(dims)
    a0 = aa[0]
    a1 = aa[1]
    aa2 = np.array([1 for i in range(D)])
    aa3 = np.array([CX_func(i, a0, a1, dims) for i in range(D)])
    aa4 = np.array([i for i in range(D)])
    return scipy.sparse.csr_matrix((aa2, (aa3, aa4)))


def SWAP_func(i, a0, a1, dims):
    aa = int_to_bases(i, dims)
    aa[a0], aa[a1] = aa[a1], aa[a0]
    return sum([aa[i] * int(np.prod(dims[:i])) for i in range(len(aa))])


def SWAP_unitary(aa, dims):
    """
    Function for mapping SWAP gate to corresponding unitary matrix.

    Args:
        aa (list[int]): Targets aa[0] and aa[1].
        dims (list[int]): Dimensions of the qudits.

    Returns:
        scipy.sparse.csr_matrix: Matrix representation of q-dimensional SWAP(aa[0],aa[1]).
    """
    D = np.prod(dims)
    a0 = q - 1 - aa[0]
    a1 = q - 1 - aa[1]
    aa2 = np.array([1 for i in range(D)])
    aa3 = np.array([i for i in range(D)])
    aa4 = np.array([SWAP_func(i, a0, a1, dims) for i in range(D)])
    return scipy.sparse.csr_matrix((aa2, (aa3, aa4)))


def diagonalize(P):
    """
    Returns the circuit which diagonalizes a pairwise commuting pauli object.

    Args:
        P: pauli
            Pauli to be diagonalized

    Returns:
        circuit
            Circuit which diagonalizes P
    """
    dims = P.dims
    q = P.qudits()

    if not P.is_commuting():
        raise Exception("Paulis must be pairwise commuting to be diagonalized")
    P1 = P.copy()
    C = circuit(dims)

    if P.is_quditwise_commuting():
        # for each dimension, call diagonalize_iter_quditwise_ on the qudits of the same dimension
        for d in sorted(set(dims)):
            aa = [i for i in range(q) if dims[i] == d]
            while aa:
                C = diagonalize_iter_quditwise_(P1, C, aa)
        P1 = act(P1, C)
    else:
        # for each dimension, call diagonalize_iter_ on the qudits of same dimension
        for d in sorted(set(dims)):
            aa = [i for i in range(q) if dims[i] == d]
            while aa:
                C = diagonalize_iter_(P1, C, aa)
        P1 = act(P1, C)

    # if any qudits are X rather than Z, apply H to make them Z
    if [i for i in range(q) if any(P1.X[:, i])]:
        g = gate(H, [i for i in range(q) if any(P1.X[:, i])])
        C.add_gates_(g)
        P1 = act(P1, g)
    return C


def diagonalize_iter_(P, C, aa):
    """An iterative function called within diagonalize()

    Args:
        P (pauli): Pauli to be diagonalized
        C (circuit): circuit which diagonalizes first a-1 qudits
        aa (list[int]): remaining qudits of same dimension

    Returns:
        circuit: circuit which diagonalizes first a qudits
    """
    p, q = P.paulis(), P.qudits()
    P = act(P, C)
    a = aa.pop(0)

    # if all Paulis have no X-part on qudit a, return C
    if not any(P.X[:, a]):
        return C

    # set a1 to be the index of the minimum Pauli with non-zero X-part of qudit a
    a1 = min(i for i in range(p) if P.X[i, a])

    # add CNOT gates to cancel out all non-zero X-parts on Pauli a1, qudits in aa
    while any(P.X[a1, i] for i in aa):
        gg = [gate(CX, [a, i]) for i in aa if P.X[a1, i]]
        C.add_gates_(gg)
        P = act(P, gg)

    # check whether there are any non-zero Z-parts on Pauli a1, qudits in aa
    while any(P.Z[a1, i] for i in aa):

        # if Pauli a1, qudit a is X, apply S gate to make it Y
        if not P.Z[a1, a]:
            g = gate(S, [a])
            C.add_gates_(g)
            P = act(P, g)

        # add backwards CNOT gates to cancel out all non-zero Z-parts on Pauli a1, qudits in aa
        gg = [gate(CX, [i, a]) for i in aa if P.Z[a1, i]]
        C.add_gates_(gg)
        P = act(P, gg)

    # if Pauli a1, qudit a is Y, add S gate to make it X
    while P.Z[a1, a]:
        g = gate(S, [a])
        C.add_gates_(g)
        P = act(P, g)
    return C


def diagonalize_iter_quditwise_(P, C, aa):
    """An iterative function called within diagonalize()

    Args:
        P (pauli): Pauli to be diagonalized
        C (circuit): circuit which diagonalizes first a-1 qudits
        aa (list[int]): remaining qudits of same dimension

    Returns:
        circuit: circuit which diagonalizes first a qudits
    """
    p, q = P.paulis(), P.qudits()
    P = act(P, C)
    a = aa.pop(0)

    # if all Paulis have no X-part on qudit a, return C
    if not any(P.X[:, a]):
        return C

    # set a1 to be the index of the minimum Pauli with non-zero X-part of qudit a
    a1 = min(i for i in range(p) if P.X[i, a])

    # if Pauli a1, qudit a is Y, add S gate to make it X
    while P.Z[a1, a]:
        g = gate(S, [a])
        C.add_gates_(g)
        P = act(P, g)
    return C


def is_diagonalizing_circuit(P, C, aa):
    """
    Checks whether the circuit properly diagonalizes a subset of Paulis.

    Args:
        P (pauli): The Pauli object to be diagonalized.
        C (circuit): The circuit that should diagonalize the Pauli object.
        aa (list[int]): Indices of the subset of Paulis to be considered.

    Returns:
        bool: True if the circuit properly diagonalizes the subset of Paulis, False otherwise.
    """
    P1 = P.copy()
    P1.delete_paulis_([i for i in range(P.paulis()) if i not in aa])
    P1 = act(P1, C)
    return P1.is_IZ()


class graph:
    """A class for storing graphs as adjacency matrices
    Since we are dealing with covariance matrices with both vertex and edge weights,
    this is a suitable format to capture that complexity
    """
    def __init__(self, adj_mat=np.array([]), dtype=complex):
        """A class for storing graphs as adjacency matrices
        Since we are dealing with covariance matrices with both vertex and edge weights,
        this is a suitable format to capture that complexity

        Args:
            adj_mat (numpy.array): (weighted) adjacency matrix of graph
            dtype   (numpy.dtype): data type of graph weights
        """
        self.adj = adj_mat.astype(dtype)

    def add_vertex_(self, c=1):
        """Adds a vertex to self.

        Args:
            c (float): vertex weight
        """
        if len(self.adj) == 0:
            self.adj = np.array([c])
        else:
            m0 = np.zeros((len(self.adj), 1))
            m1 = np.zeros((1, len(self.adj)))
            m2 = np.array([[c]])
            self.adj = np.block([[self.adj, m0], [m1, m2]])

    def lade_vertex_(self, a, c):
        """Weight a vertex
        Args:
            a - (int)   - vertex to be weighted
            c - (float) - vertex weight
        """
        self.adj[a, a] = c

    def lade_edge_(self, a0, a1, c):
        """Weight an edge.

        Args:
            a0 (int): First vertex.
            a1 (int): Second vertex.
            c (float): Vertex weight.
        """
        self.adj[a0, a1] = c
        self.adj[a1, a0] = c

    def neighbors(self, a):
        """Returns a set of the neighbors of a given vertex

        Args:
            a (int): vertex for which neighbors should be returned

        Returns:
            set[int]: set of neighbors of vertex a
        """
        aa1 = set([])
        for i in range(self.ord()):
            if (a != i) and (self.adj[a, i] != 0):
                aa1.add(i)
        return aa1

    def edges(self):
        """Returns list of all edges in self.

        Returns:
            list[list[int]]: list of edges in self
        """
        aaa = []
        for i0, i1 in itertools.combinations(range(self.ord()), 2):
            if i1 in self.neighbors(i0):
                aaa.append([i0, i1])
        return aaa

    def clique(self, aa):
        """Check whether a collection of vertices is a clique in self

        Args:
            aa (list[int]): list of vertices to be checked for clique

        Returns:
            bool: True if aa is a clique in self; False otherwise
        """
        for i0, i1 in itertools.combinations(aa, 2):
            if self.adj[i0, i1] == 0:
                return False
        return True

    def degree(self, a):
        """
        Returns the degree of a given vertex.

        Args:
            a (int): Vertex for which degree should be returned.

        Returns:
            int: Degree of vertex a.
        """
        return np.count_nonzero(self.adj[a, :])

    def ord(self):
        """Returns the number of vertices in self.

        Returns:
            int: number of vertices in self
        """
        return self.adj.shape[0]

    def print(self):
        """Print adjacency matrix representation of self"""
        for i0 in range(self.ord()):
            print('[', end=' ')
            for i1 in range(self.ord()):
                s = self.adj[i0, i1]
                if str(s)[0] == '-':
                    print(f'{self.adj[i0,i1]:.2f}', end=" ")
                else:
                    print(' ' + f'{self.adj[i0,i1]:.2f}', end=" ")
            print(']')

    def print_neighbors(self):
        """Print self as a list of vertices together with their neighbors"""
        for i0 in range(self.ord()):
            print(i0, end=": ")
            for i1 in self.neighbors(i0):
                print(i1, end=" ")
            print()

    def copy(self):
        """Return a deep copy of self.

        Returns:
            graph: A deep copy of self.
        """
        return graph(np.array([[self.adj[i0, i1] for i1 in range(self.ord())] for i0 in range(self.ord)]))


def nonempty_cliques(A):
    """
    Returns all non-empty cliques in a graph.

    Args:
        A (graph): Graph for which all cliques should be found.

    Returns:
        list[list[int]]: A list containing all non-empty cliques in A.
    """
    p = A.ord()
    aaa = set([frozenset([])])
    for i in range(p):
        iset = set([i])
        inter = A.neighbors(i)
        aaa |= set([frozenset(iset | (inter & aa)) for aa in aaa])
    aaa.remove(frozenset([]))
    return list([list(aa) for aa in aaa])


def all_maximal_cliques(A):
    """
    Returns a generator over all maximal cliques in a graph.

    Args:
        A (graph): graph for which all cliques should be found

    Returns:
        list[int]: a maximal clique in A
    """
    p = A.ord()
    N = {}
    for i in range(p):
        N[i] = A.neighbors(i)
    nxG = nx.Graph()
    nxG.add_nodes_from([i for i in range(p)])
    nxG.add_edges_from([(i0, i1) for i0 in range(p) for i1 in N[i0]])
    return nx.algorithms.clique.find_cliques(nxG)


def weighted_vertex_covering_maximal_cliques(A, A1=None, cc=None, k=1):
    """Returns a clique covering of a graph which hits every vertex at least a certain number of times.

    Args:
        A (graph): commutation graph for which covering should be found
        A1 (graph, optional): variance graph for which covering should be found
        cc (list[int], optional): coefficients of the Hamiltonian
        k (int, optional): number of times each vertex should be covered

    Returns:
        list[list[int]]: a list containing cliques which cover A
    """
    p = A.ord()
    if A1 is None and cc is None:
        return vertex_covering_maximal_cliques(A, k=k)
    elif A1 is None:
        cc2 = [np.abs(cc[i])**2 for i in range(p)]
        N = {}
        for i in range(p):
            N[i] = A.neighbors(i)
        aaa = []
        for i0 in range(p):
            for i1 in range(k):
                aa0 = [i0]
                aa1 = list(N[i0])
                while aa1:
                    c1 = sum(cc2[a0] for a0 in aa0)
                    cc1 = [c1 + sum(cc2[a2] for a2 in N[a1].intersection(aa1)) for a1 in aa1]
                    if sum(cc1) == 0:
                        cc1 = [1 for a in aa1]
                    r = random.choices(aa1, cc1)[0]
                    aa0.append(r)
                    aa1 = list(N[r].intersection(aa1))
                aaa.append(aa0)
        return [sorted(list(aa1)) for aa1 in set([frozenset(aa) for aa in aaa])]
    else:
        V1 = A1.adj
        N = {}
        for i in range(p):
            N[i] = A.neighbors(i)
        N2 = {}
        for i in range(p):
            N2[i] = A.neighbors(i) | set([i])
        aaa = []
        for i0 in range(p):
            for i1 in range(k):
                aa0 = [i0]
                aa1 = list(N[i0])
                aa2 = aa0 + aa1
                while aa1:
                    cc1 = [V1[list(N2[a1].intersection(aa2))][:, list(N2[a1].intersection(aa2))].sum() for a1 in aa1]
                    if sum(cc1) == 0:
                        cc1 = [1 for a in aa1]
                    r = random.choices(aa1, cc1)[0]
                    aa0.append(r)
                    aa1 = list(N[r].intersection(aa1))
                    aa2 = aa0 + aa1
                aaa.append(aa0)
        return [sorted(list(aa1)) for aa1 in set([frozenset(aa) for aa in aaa])]


def vertex_covering_maximal_cliques(A, k=1):
    """Returns a clique covering of a graph which hits every vertex at least a certain number of times.

    Args:
        A - (graph) - commutation graph for which covering should be found
        k - (int)   - number of times each vertex must be covered

    Returns:
        (list{list{int}}) - a list containing cliques which cover A
    """
    p = A.ord()
    N = {}
    for i in range(p):
        N[i] = A.neighbors(i)
    aaa = []
    for i0 in range(p):
        for i1 in range(k):
            aa0 = [i0]
            aa1 = list(N[i0])
            while aa1:
                cc = [len(N[a1].intersection(aa1)) for a1 in aa1]
                if sum(cc) == 0:
                    cc = [1 for a in aa1]
                r = random.choices(aa1, cc)[0]
                aa0.append(r)
                aa1 = list(N[r].intersection(aa1))
            aaa.append(aa0)
    return [sorted(list(aa1)) for aa1 in set([frozenset(aa) for aa in aaa])]


def post_process_cliques(A, aaa, k=1):
    """Reduces a clique covering of a graph by removing cliques with the lowest weight.

    Args:
        A (graph): Variance graph from which weights of cliques can be obtained.
        aaa (list[list[int]]): A clique covering of the Hamiltonian.
        k (int): Number of times each vertex must be covered.

    Returns:
        list[list[int]]: A list containing cliques which cover A.
    """
    p = A.ord()
    V = A.adj
    s = np.array([sum([i in aa for aa in aaa]) for i in range(p)])
    D = {}
    for aa in aaa:
        D[str(aa)] = V[aa][:, aa].sum()
    aaa1 = aaa.copy()
    aaa1 = list(filter(lambda x: all(a >= (k + 1) for a in s[aa]), aaa1))
    while aaa1:
        aa = min(aaa1, key=lambda x: D[str(x)])
        aaa.remove(aa)
        aaa1.remove(aa)
        s -= np.array([int(i in aa) for i in range(p)])
        aaa1 = list(filter(lambda x: all(a >= (k + 1) for a in s[aa]), aaa1))
    return aaa


def LDF(A):
    """Returns a largest-degree-first clique partition of a graph.

    Args:
        A (graph): graph for which partition should be found

    Returns:
        list[list[int]]: a list containing cliques which partition A
    """
    p = A.ord()
    remaining = set(range(p))
    N = {}
    for i in range(p):
        N[i] = A.neighbors(i)
    aaa = []
    while remaining:
        a = max(remaining, key=lambda x: len(N[x] & remaining))
        aa0 = set([a])
        aa1 = N[a] & remaining
        while aa1:
            a2 = max(aa1, key=lambda x: len(N[x] & aa1))
            aa0.add(a2)
            aa1 &= N[a2]
        aaa.append(aa0)
        remaining -= aa0
    return [sorted(list(aa)) for aa in aaa]


# PHYSICS FUNCTIONS
def tensor(mm):
    """Returns the tensor product of a list of matrices.

    Args:
        mm: list of scipy.sparse.csr_matrix, matrices to tensor

    Returns:
        scipy.sparse.csr_matrix, tensor product of matrices
    """
    if len(mm) == 0:
        return matrix([])
    elif len(mm) == 1:
        return mm[0]
    else:
        return scipy.sparse.kron(mm[0], tensor(mm[1:]), format="csr")


def Mean(P, psi):
    """Returns the mean of a single Pauli with a given state.

    Args:
        P: pauli, Pauli for mean
        psi: numpy.array, state for mean

    Returns:
        numpy.float64, mean <psi|P|psi>
    """
    m = pauli_to_matrix(P)
    psi_dag = psi.conj().T
    mean = psi_dag @ m @ psi
    return mean


def Hamiltonian_Mean(P, cc, psi):
    """Returns the mean of a Hamiltonian with a given state

    Args:
        P: pauli, Paulis of Hamiltonian
        cc: list[float], coefficients of Hamiltonian
        psi: numpy.array, state for mean

    Returns:
        numpy.float64, mean sum(c*<psi|P|psi>)
    """
    p = P.paulis()
    return sum(cc[i] * Mean(P.a_pauli(i), psi) for i in range(p))


def Var(P, psi):
    """Returns the variance of a single Pauli with a given state

    Args:
        P: pauli, Pauli for variance
        psi: numpy.array, state for variance

    Returns:
        numpy.float64, variance <psi|P^2|psi> - <psi|P|psi>^2
    """
    m = pauli_to_matrix(P)
    psi_dag = psi.conj().T
    var = (psi_dag @ m @ m @ psi) - (psi_dag @ m @ psi)**2
    return var.real


def Cov(P0, P1, psi):
    """Returns the covariance of two single Paulis with a given state.

    Args:
        P0: pauli, first Pauli for covariance
        P1: pauli, second Pauli for covariance
        psi: numpy.array, state for covariance

    Returns:
        numpy.float64, covariance <psi|P0P1|psi> - <psi|P0|psi><psi|P1|psi>
    """
    m0 = pauli_to_matrix(P0)
    m1 = pauli_to_matrix(P1)
    psi_dag = psi.conj().T
    cov = (psi_dag @ m0 @ m1 @ psi) - (psi_dag @ m0 @ psi) * (psi_dag @ m1 @ psi)
    return cov.real


def variance_graph(P, cc, psi):
    """Returns the graph of variances and covariances for a given Hamiltonian and ground state

    Args:
        P: pauli, set of Paulis in Hamiltonian
        cc: list[float], coefficients in Hamiltonian
        psi: numpy.array, ground state

    Returns:
        graph: variances and covariances of all Paulis with respect to ground state
    """
    p = P.paulis()
    mm = [pauli_to_matrix(P.a_pauli(i)) for i in range(p)]
    psi_dag = psi.conj().T
    cc1 = [psi_dag @ mm[i] @ psi for i in range(p)]
    cc2 = [psi_dag @ mm[i].conj().T @ psi for i in range(p)]
    return graph(np.array([[np.conj(cc[i0]) * cc[i1] * ((psi_dag @ mm[i0].conj().T @ mm[i1] @ psi) - cc2[i0] * cc1[i1]) for i1 in range(p)] for i0 in range(p)]))


def scale_variances(A, S):
    """Scales the entries in a variance graph with respect to number of measurements.

    Args:
        A: graph, variance matrix
        S: numpy.array, array for tracking number of measurements

    Returns:
        graph, scaled variance matrix
    """
    p = A.ord()
    S1 = S.copy()
    S1[range(p), range(p)] = [a if a != 0 else 1 for a in S1.diagonal()]
    s1 = 1 / S1.diagonal()
    return graph(S1 * A.adj * s1 * s1[:, None])


def commutation_graph(P):
    """Returns the commutation graph of a given Pauli.

    Args:
        P: pauli, Pauli to check for commutation relations

    Returns:
        graph: an edge is weighted 1 if the pair of Paulis commute
    """
    p = P.paulis()
    return graph(np.array([[1 - symplectic_inner_product(P.a_pauli(i0), P.a_pauli(i1)) for i1 in range(p)] for i0 in range(p)]))


def quditwise_commutation_graph(P):
    """Returns the quditwise commutation graph of a given Pauli.

    Args:
        P: pauli, Pauli to check for quditwise commutation relations

    Returns:
        graph: an edge is weighted 1 if the pair of Paulis quditwise commute
    """
    p = P.paulis()
    return graph(np.array([[1 - quditwise_inner_product(P.a_pauli(i0), P.a_pauli(i1)) for i1 in range(p)] for i0 in range(p)]))


def random_Ham(p, q, d):
    """Returns a random Hamiltonian with given number of Paulis, number of qudits, and Pauli weight

    Args:
        p (int): number of Paulis
        q (int): number of qudits
        d (int): max Pauli weight

    Returns:
        pauli: random set of Paulis satisfying input conditions
    """
    sss = []
    ssdict = {0: "I", 1: "Z", 2: "X", 3: "Y"}
    for i in range(p):
        rr = random.sample(range(q), d)
        sss.append("".join([ssdict[random.randint(0, 3)] if i1 in rr else "I" for i1 in range(q)]))
    return string_to_pauli(sss)


def print_Ham_string(P, cc):
    """Print list of Paulis in string form, together with coefficients

    Args:
        P: pauli, Pauli to be printed
        cc: list[int], coefficients for Hamiltonian

    Returns:
        None
    """
    X, Z = P.X, P.Z
    for i in range(P.paulis()):
        print(pauli_to_string(P.a_pauli(i)), end="")
        print('', cc[i])


def ground_state(P):
    """Returns the ground state of a given Hamiltonian

    Args:
        P: pauli, Paulis for Hamiltonian
        cc: list[int], coefficients for Hamiltonian

    Returns:
        numpy.array: eigenvector corresponding to lowest eigenvalue of Hamiltonian
    """
    m = P.matrix_form()  # sum(pauli_to_matrix(P.a_pauli(i)) * cc[i] for i in range(P.paulis()))

    m = m.toarray()
    print(scipy.linalg.ishermitian(np.around(m, 3)))
    val, vec = np.linalg.eig(m)
    val = np.real(val)
    vec = np.transpose(vec)

    tmp_index = val.argmin(axis=0)

    # print(val)
    # print(tmp_index)
    # print(val[tmp_index])
    gs = vec[tmp_index]
    gs = np.transpose(gs)
    gs = gs / np.linalg.norm(gs)

    if abs(min(val) - np.transpose(np.conjugate(gs)) @ m @ gs) > 10**-10:
        print("ERROR with the GS!!!")

    # print(gs)
    # print(np.linalg.norm(gs))
    # print(np.transpose(np.conjugate(gs))@m@gs)

    return gs

    # gval,gvec = scipy.sparse.linalg.eigsh(m,which='SA',k=1)
    # tmp_state = np.array([g for g in gvec[:,0]])
    # return tmp_state/np.linalg.norm(tmp_state)


# MEASUREMENT FUNCTIONS
# sample from distribution given by ground state and eigenstates of clique
#     optionally input a dictionary, which will be updated to track speed up future samples
def sample_(P, psi, aa, D={}):
    """Sample from distribution given by ground state and eigenstates of clique.

    Args:
        P: pauli, Paulis for Hamiltonian
        psi: numpy.array, ground state of Hamiltonian
        aa: list[int], clique to be measured
        D: dict, dictionary for storing pdf and negations for future samples

    Returns:
        list[int]: ith entry is +1/-1 for measurement outcome on ith element of aa
    """
    if str(aa) in D.keys():
        P1, pdf = D[str(aa)]
    else:
        P1 = P.copy()
        P1.delete_paulis_([i for i in range(P.paulis()) if i not in aa])
        C = diagonalize(P1)
        psi_diag = C.unitary() @ psi
        pdf = np.absolute(psi_diag * psi_diag.conj())
        P1 = act(P1, C)
        D[str(aa)] = (P1, pdf)
    p, q, phases, dims = P1.paulis(), P1.qudits(), P1.phases, P1.dims
    a1 = np.random.choice(np.prod(dims), p=pdf)
    bases_a1 = int_to_bases(a1, dims)
    ss = [(phases[i0] + sum((bases_a1[i1] * P1.Z[i0, i1] * P1.lcm) // P1.dims[i1] for i1 in range(q))) % P1.lcm for i0 in range(p)]
    return ss


# ESTIMATED PHYSICS FUNCTIONS
# Bayesian estimation of mean from samples
def bayes_Mean(xDict):
    """Bayesian estimation of mean from samples.

    Args:
        xDict (Dict): number of ++/+-/-+/-- outcomes for single Pauli

    Returns:
        float: Bayesian estimate of mean
    """
    x0, x1 = xDict[(1, 1)], xDict[(-1, -1)]
    return (x0 - x1) / (x0 + x1 + 2)


def bayes_Var(xDict):
    """Bayesian estimation of variance from samples.

    Args:
        xDict (Dict): number of ++/+-/-+/-- outcomes for single Pauli

    Returns:
        float: Bayesian variance of mean
    """
    lcm = int(np.sqrt(len(xDict)))
    alpha = [math.e**(2 * math.pi * 1j * i / lcm) for i in range(lcm)]
    alpha_conj = [math.e**(-2 * math.pi * 1j * i / lcm) for i in range(lcm)]
    s = sum(xDict[(i, i)] for i in range(lcm))
    return sum(alpha[i] * (alpha_conj[i] - alpha_conj[j]) * (xDict[(i, i)] + 1) * (xDict[(j, j)] + 1) / ((s + lcm) * (s + lcm + 1)) for i, j in itertools.product(range(lcm), repeat=2))


def bayes_Cov(xyDict, xDict, yDict):
    """Bayesian estimation of covariance from samples.

    Args:
        xyDict (Dict): number of ++/+-/-+/-- outcomes for pair of Paulis
        xDict (Dict): number of ++/+-/-+/-- outcomes for first Pauli
        yDict (Dict): number of ++/+-/-+/-- outcomes for second Pauli

    Returns:
        float: Bayesian estimate of covariance
    """
    lcm = int(np.sqrt(len(xyDict)))
    alpha = [math.e**(2 * math.pi * 1j * i / lcm) for i in range(lcm)]
    alpha_conj = [math.e**(-2 * math.pi * 1j * i / lcm) for i in range(lcm)]
    s = sum(xyDict[(i, j)] for i, j in itertools.product(range(lcm), repeat=2))
    return sum((alpha_conj[i0] - alpha_conj[j0]) * (alpha[i1] - alpha[j1]) * (xyDict[(i0, i1)] + 1) * (xyDict[(j0, j1)] + 1) / (2 * (s + lcm) * (s + lcm + 1)) for (i0, i1), (j0, j1) in itertools.product(itertools.product(range(lcm), repeat=2), repeat=2))
    # xy00,xy01,xy10,xy11 = xyDict[(1,1)],xyDict[(1,-1)],xyDict[(-1,1)],xyDict[(-1,-1)]
    # x0,x1 = xDict[(1,1)],xDict[(-1,-1)]
    # y0,y1 = yDict[(1,1)],yDict[(-1,-1)]
    # p00 = 4*((x0+1)*(y0+1))/((x0+x1+2)*(y0+y1+2))
    # p01 = 4*((x0+1)*(y1+1))/((x0+x1+2)*(y0+y1+2))
    # p10 = 4*((x1+1)*(y0+1))/((x0+x1+2)*(y0+y1+2))
    # p11 = 4*((x1+1)*(y1+1))/((x0+x1+2)*(y0+y1+2))
    # return 4*((xy00+p00)*(xy11+p11) - (xy01+p01)*(xy10+p10))/((xy00+xy01+xy10+xy11+4)*(xy00+xy01+xy10+xy11+5))


def bayes_variance_graph(X, cc):
    """
    Approximates the variance graph using Bayesian estimates.

    Args:
        X (np.ndarray of dict): Array for tracking measurement outcomes.
        cc (list of float): Coefficients of the Hamiltonian.

    Returns:
        np.ndarray: Variance graph calculated with Bayesian estimates.
    """
    p = len(cc)
    cc_conj = [np.conj(c) for c in cc]
    return graph(np.array([[cc_conj[i0] * cc[i0] * bayes_Var(X[i0, i0]) if i0 == i1 else cc_conj[i0] * cc[i1] * bayes_Cov(X[i0, i1], X[i0, i0], X[i1, i1]) for i1 in range(p)] for i0 in range(p)]))


def naive_Mean(xDict):
    """Naive estimation of mean from samples.

    Args:
        xDict (Dict): number of ++/+-/-+/-- outcomes for single Pauli

    Returns:
        float: Bayesian estimate of mean
    """
    x0, x1 = xDict[(1, 1)], xDict[(-1, -1)]
    if (x0 + x1) == 0:
        return 0
    return (x0 - x1) / (x0 + x1)


def naive_Var(xDict):
    """Naive estimation of variance from samples.

    Args:
        xDict (Dict): number of ++/+-/-+/-- outcomes for single Pauli

    Returns:
        float: Bayesian variance of mean
    """
    x0, x1 = xDict[(1, 1)], xDict[(-1, -1)]
    if (x0 + x1) == 0:
        return 2 / 3
    return 4 * (x0 * x1) / ((x0 + x1) * (x0 + x1))


def naive_Cov(xyDict, xDict, yDict):
    """Naive estimation of covariance from samples.

    Args:
        xyDict (Dict): number of ++/+-/-+/-- outcomes for pair of Paulis
        xDict (Dict): number of ++/+-/-+/-- outcomes for first Pauli
        yDict (Dict): number of ++/+-/-+/-- outcomes for second Pauli

    Returns:
        float: naive estimate of covariance
    """
    xy00, xy01, xy10, xy11 = xyDict[(1, 1)], xyDict[(1, -1)], xyDict[(-1, 1)], xyDict[(-1, -1)]
    x0, x1 = xDict[(1, 1)], xDict[(-1, -1)]
    y0, y1 = yDict[(1, 1)], yDict[(-1, -1)]
    if (xy00 + xy01 + xy10 + xy11) == 0:
        return 0
    return 4 * ((xy00) * (xy11) - (xy01) * (xy10)) / ((xy00 + xy01 + xy10 + xy11) * (xy00 + xy01 + xy10 + xy11))


def naive_variance_graph(X, cc):
    """Approximates the variance graph using naive estimates.

    Args:
        X (np.ndarray of dict): Array for tracking measurement outcomes.
        cc (list of float): Coefficients of the Hamiltonian.

    Returns:
        np.ndarray: Variance graph calculated with naive estimates.
    """
    p = len(cc)
    return graph(np.array([[(cc[i0]**2) * naive_Var(X[i0, i0]) if i0 == i1 else cc[i0] * cc[i1] * naive_Cov(X[i0, i1], X[i0, i0], X[i1, i1]) for i1 in range(p)] for i0 in range(p)]))


# SIMULATION ALGORITHMS
# convert from L,l notation to set of update steps
def Ll_updates(L, l, shots):
    """Convert from L,l notation to set of update steps.

    Args:
        L (int): number of sections into which shots should be split
        l (int): exponential scaling factor for size of sections
        shots (int): total number of shots required

    Returns:
        set[int]: set containing steps at which algorithm should update
    """
    r0_shots = shots / sum([(1 + l)**i for i in range(L)])
    shot_nums = [round(r0_shots * (1 + l)**i) for i in range(L - 1)]
    shot_nums.append(shots - sum(shot_nums))
    return set([0] + list(itertools.accumulate(shot_nums))[:-1])


def variance_estimate_(P, cc, psi, D, X, xxx):
    """Updates the variance matrix by sampling from pre-determined cliques.

    Args:
        P: pauli
            Paulis in Hamiltonian
        cc: list[int]
            coefficients in Hamiltonian
        psi: numpy.array
            ground state of Hamiltonian
        D: dict
            dictionary for storing pdf and negations for future samples
        X: numpy.array[dict]
            array of measurement outcome counts
        xxx: list[list[int]]
            list of cliques to-be-sampled

    Returns:
        tuple: (numpy.array[float], dict, numpy.array[dict]): (variance graph
        calculated with Bayesian estimates, (updated) dictionary for storing pdf and negations for future samples, (updated) array of measurement outcome counts)
    """
    p = P.paulis()
    index_set = set(range(p))
    for aa in xxx:
        aa1 = sorted(index_set.difference(aa))
        cc1 = sample_(P, psi, aa, D)
        for (a0, c0), (a1, c1) in itertools.product(zip(aa, cc1), repeat=2):
            X[a0, a1][(c0, c1)] += 1
    return bayes_variance_graph(X, cc).adj, D, X


def bucket_filling(P, cc, psi, shots, part_func, update_steps=set([]), repeats=(0, 1), full_simulation=False, general_commutation=True):
    """Partitions Hamiltonian and repeatedly samples cliques while minimizing total variance.

    Args:
        P: (pauli) - Paulis in Hamiltonian
        cc: (list{int}) - coefficients in Hamiltonian
        psi: (numpy.array) - ground state of Hamiltonian
        shots: (int) - number of samples to take
        part_func: (function) - function for determining partition
        update_steps: (set{int}) - steps at which variance graph should be updated
        repeats: (tuple{int}) - current iteration and total number of iterations
        full_simulation: (bool) - set True if full simulation is required
        general_commutation: (bool) - set True if general commutation is allowed

    Returns:
        (numpy.array{int}) - array containing number of times each pair of Paulis was measured together
        (numpy.array{dict}) - array of measurement outcome counts
        (list{list{int}}) - list of cliques which were sampled
    """
    p = P.paulis()
    X = np.array([[dict(zip([(i0, i1) for i0, i1 in itertools.product(range(P.lcm), repeat=2)], [0] * (P.lcm**2))) for a1 in range(p)] for a0 in range(p)])
    if general_commutation:
        CG = commutation_graph(P)
    else:
        CG = quditwise_commutation_graph(P)
    if part_func == weighted_vertex_covering_maximal_cliques:
        aaa = part_func(CG, cc=cc, k=3)
    else:
        aaa = part_func(CG)
    D = {}
    S = np.zeros((p, p), dtype=int)
    Ones = [np.ones((i, i), dtype=int) for i in range(p + 1)]
    index_set = set(range(p))
    xxx = []
    xxx1 = []
    S[range(p), range(p)] += np.array([1 for i in range(p)])
    for i0 in range(shots):
        if i0 == 0 or i0 in update_steps:
            V, D, X = variance_estimate_(P, cc, psi, D, X, xxx1)
            xxx1 = []
        S1 = S + Ones[p]
        s = 1 / (S.diagonal() | (S.diagonal() == 0))
        s1 = 1 / S1.diagonal()
        factor = p - np.count_nonzero(S.diagonal())
        S1[range(p), range(p)] = [a if a != 1 else -factor for a in S1.diagonal()]
        V1 = V * (S * s * s[:, None] - S1 * s1 * s1[:, None])
        V2 = 2 * V * (S * s * s[:, None] - S * s * s1[:, None])
        aaa, aaa1 = itertools.tee(aaa, 2)
        # aa = sorted(max(aaa1,key=lambda xx : np.abs(V1[xx][:,xx].sum()+V2[xx][:,list(index_set.difference(xx))].sum())))
        aa = sorted(random.sample(list(set([frozenset(aa1) for aa1 in aaa1])), 1)[0])
        xxx.append(aa)
        xxx1.append(aa)
        S[np.ix_(aa, aa)] += Ones[len(aa)]
        loading_bar([(i0, shots), repeats], scalings=[lambda x:x**(3 / 2)])

    S[range(p), range(p)] -= np.array([1 for i in range(p)])
    if full_simulation:
        for aa in xxx1:
            aa1 = sorted(index_set.difference(aa))
            cc1 = sample_(P, psi, aa, D)
            for (a0, c0), (a1, c1) in itertools.product(zip(aa, cc1), repeat=2):
                X[a0, a1][(c0, c1)] += 1
    else:
        X = None
    return S, X, xxx


def bucket_filling_mod(P,
                       cc,
                       psi,
                       shots,
                       part_func,
                       update_steps=set([]),
                       repeats=(0, 1),
                       full_simulation=False,
                       general_commutation=True,
                       best_possible=False):
    """Partitions Hamiltonian and repeatedly samples cliques while minimizing total variance.

    Args:
        P: pauli, Paulis in Hamiltonian
        cc: list of int, coefficients in Hamiltonian
        psi: numpy.array, ground state of Hamiltonian
        shots: int, number of samples to take
        part_func: function, function for determining partition
        update_steps: set of int, steps at which variance graph should be updated
        repeats: tuple of int, current iteration and total number of iterations
        full_simulation: bool, set True if full simulation is required
        general_commutation: bool, set True if general commutation is allowed

    Returns:
        numpy.array of int, array containing number of times each pair of Paulis was measured together
        numpy.array of dict, array of measurement outcome counts
        list of list of int, list of cliques which were sampled
    """
    if best_possible:
        vg = variance_graph(P, cc, psi)

    p = P.paulis()
    X = np.array([[dict(zip([(i0, i1) for i0, i1 in itertools.product(range(P.lcm), repeat=2)], [0] * (P.lcm**2))) for a1 in range(p)] for a0 in range(p)])
    if general_commutation:
        CG = commutation_graph(P)
    else:
        CG = quditwise_commutation_graph(P)
    if part_func == weighted_vertex_covering_maximal_cliques:
        aaa = part_func(CG, cc=cc, k=3)
    else:
        aaa = part_func(CG)
    D = {}
    S = np.zeros((p, p), dtype=int)
    Ones = [np.ones((i, i), dtype=int) for i in range(p + 1)]
    index_set = set(range(p))
    xxx = []
    xxx1 = []
    S[range(p), range(p)] += np.array([1 for i in range(p)])
    for i0 in range(shots):
        if i0 == 0 or i0 in update_steps:
            V, D, X = variance_estimate_(P, cc, psi, D, X, xxx1)
            if best_possible == True:
                V = vg.adj
            xxx1 = []
        S1 = S + Ones[p]
        s = 1 / (S.diagonal() | (S.diagonal() == 0))
        s1 = 1 / S1.diagonal()
        factor = p - np.count_nonzero(S.diagonal())
        S1[range(p), range(p)] = [a if a != 1 else -factor for a in S1.diagonal()]
        V1 = V * (S * s * s[:, None] - S1 * s1 * s1[:, None])
        V2 = 2 * V * (S * s * s[:, None] - S * s * s1[:, None])
        aaa, aaa1 = itertools.tee(aaa, 2)
        aa = sorted(max(aaa1, key=lambda xx: np.abs(V1[xx][:, xx].sum() + V2[xx][:, list(index_set.difference(xx))].sum())))
        # aa = sorted(random.sample(list(set([frozenset(aa1) for aa1 in aaa1])),1)[0])
        xxx.append(aa)
        xxx1.append(aa)
        S[np.ix_(aa, aa)] += Ones[len(aa)]
        loading_bar([(i0, shots), repeats], scalings=[lambda x:x**(3 / 2)])

    S[range(p), range(p)] -= np.array([1 for i in range(p)])
    if full_simulation:
        for aa in xxx1:
            aa1 = sorted(index_set.difference(aa))
            cc1 = sample_(P, psi, aa, D)
            for (a0, c0), (a1, c1) in itertools.product(zip(aa, cc1), repeat=2):
                X[a0, a1][(c0, c1)] += 1
    else:
        X = None
    return S, X, xxx


def equal_allocation_algorithm(P, cc, general_commutation=True):
    """
    Provides a dictionary which contains the circuits, Paulis, eigenvalues, et al. for every group in the partition
    Returns a dictionary indexed by the indices of the sets of Paulis (converetd to strings) containing:
    a circuit which diagonalizes these Paulis, alist of Paulis within the partition,
    a list of Paulis within the partition after applying diagonalization circuit, a list of coefficients corresponding to these Paulis,
    a list of tuples of measurement outcomes and eigenvalues (each measurement outcome is paired with a list of the corresponding eigenvalues of each Pauli
    the eigenvalues are stored as integers modulo the least common multiple of the dimensions)

    Args:
        P (pauli): Paulis in Hamiltonian
        cc (list[int]): coefficients in Hamiltonian
        general_commutation (bool): set True if general commutation is allowed

    Returns:
        dict: return dictionary.
    """
    if general_commutation:
        aaa = LDF(commutation_graph(P))
    else:
        aaa = LDF(quditwise_commutation_graph(P))
    measurement_dictionary = dict([])
    for aa in aaa:
        aa = sorted(aa)
        P1 = P.copy()
        P1.delete_paulis_([i for i in range(P.paulis()) if i not in aa])
        C = diagonalize(P1)
        P1 = act(P1, C)
        pauli_list = [P.a_pauli(a) for a in aa]
        diagonalized_pauli_list = [P1.a_pauli(i) for i in range(len(aa))]
        coefficient_list = [cc[a] for a in aa]
        eigenvalues_list = []
        for a1 in range(np.prod(P.dims)):
            eigenvalues_list.append([(P1.phases[i0] + sum((int_to_bases(a1, P1.dims)[i1] * P1.Z[i0, i1] * P1.lcm) // P1.dims[i1] for i1 in range(P1.qudits()))) % P1.lcm for i0 in range(P1.paulis())])
        measurement_dictionary[str(aa)] = (C, pauli_list, diagonalized_pauli_list, coefficient_list, eigenvalues_list)
    return measurement_dictionary


# determines the number of measurements required to reach desired error (as a proportion of the true mean)
#     note that this algorithm calculates the true ground state, mean, and variance
#     so it should be used to observe trends on small examples, but isn't realistic for large examples
def equal_allocation_measurements(P, cc, error, general_commutation=True):
    """
    Calculates the number of measurements required to reach a desired error
    (as a proportion of the true mean) for each clique in the commutation graph.

    Args:
        P: (pauli)     - Paulis in Hamiltonian
        cc: (list{int}) - coefficients in Hamiltonian
        error: (float)     - desired error as a proportion of the true mean
        general_commutation: (bool)      - set True if general commutation is allowed

    Returns:
        int: number of measurements each clique must be measured to reach this error
    """
    p = P.paulis()
    psi = ground_state(P, cc)
    true_Mean = Hamiltonian_Mean(P, cc, psi).real
    true_Variance_graph = variance_graph(P, cc, psi)
    if general_commutation:
        aaa = LDF(commutation_graph(P))
    else:
        aaa = LDF(quditwise_commutation_graph(P))
    S = np.zeros((p, p), dtype=int)
    Ones = [np.ones((i, i), dtype=int) for i in range(p + 1)]
    for aa in aaa:
        S[np.ix_(aa, aa)] += Ones[len(aa)]
    one_shot_error = np.sqrt(np.sum(scale_variances(true_Variance_graph, S).adj).real)
    return math.ceil((one_shot_error / (error * true_Mean))**2)


# EXTRA STUFF:
# Function to reconstruct expectation values (and errors)
def expt_rec(beta_exp,
             Ptot,
             Pe,
             Pm,
             cce,
             ccm,
             elmag,
             Michael_state=True,
             shots=1000,
             part_func=LDF,
             no_simulation=False,
             full_simulation=True,
             update_steps=set([10000, 100000]),
             printing=False,
             general_commutation=False):
    # beta value and coefficients of total Hamiltonian:
    beta = 10**beta_exp
    cct = ([beta * el for el in ccm]) + ([beta**-1 * el for el in cce])

    # parameters required below:
    pt, qt = Ptot.paulis(), Ptot.qudits()
    pe, qe = Pe.paulis(), Pe.qudits()
    pm, qm = Pm.paulis(), Pm.qudits()

    # import corresponding ground state
    if Michael_state:
        if elmag == "E":
            psi = np.load("./Expt_data/dms/E_basis/state_" + str(beta_exp) + ".npy")
        elif elmag == "B":
            psi = np.load("./Expt_data/dms/B_basis/state_" + str(beta_exp) + ".npy")
        if np.max(abs(np.array(psi @ psi) - np.array(psi))) < 10**-10:
            index_chosen = np.diag(psi).argmax()
            new_psi = psi[index_chosen]
            psi = new_psi / np.linalg.norm(new_psi)
        else:
            print("STATE IS NOT PURE!!!!")
    else:
        psi = ground_state(Ptot, cct)

    # if printing, print expected values:
    if printing:
        print("chosen beta: ", beta)
        print("total energy: ", Hamiltonian_Mean(Ptot, cct, psi).real)
        print("electric dimensionless energy: ", Hamiltonian_Mean(Pe, cce, psi).real)
        print("plaquette: ", -Hamiltonian_Mean(Pm, ccm, psi).real / 4)
        print()
        print()

    # determine experimental results:
    if not no_simulation:
        St, Xt, xxxt = bucket_filling(Ptot, cct, psi, shots, part_func, full_simulation=full_simulation, update_steps=update_steps, general_commutation=general_commutation)
        if printing:
            print("first round done :)")
        Se, Xe, xxxe = bucket_filling(Pe, cce, psi, shots, part_func, full_simulation=full_simulation, update_steps=update_steps, general_commutation=general_commutation)
        if printing:
            print("second round done :)")
        Sm, Xm, xxxm = bucket_filling(Pm, ccm, psi, shots, part_func, full_simulation=full_simulation, update_steps=update_steps, general_commutation=general_commutation)
        if printing:
            print("third round done - almost there! :)")
            print()
            print()

        # Getting the actual values:
        # total energy:
        tot_en_true = Hamiltonian_Mean(Ptot, cct, psi).real
        if printing:
            print('TRUE TOTAL ENERGY')
            print('True mean:', tot_en_true)
        if Xt is not None:
            tot_en_est = sum(cct[i0] * sum(Xt[i0, i0][i1, i1] * math.e**(2 * 1j * math.pi * i1 / Ptot.lcm) for i1 in range(Ptot.lcm)) / sum(Xt[i0, i0][i1, i1] for i1 in range(Ptot.lcm)) if sum(Xt[i0, i0][i1, i1] for i1 in range(Ptot.lcm)) > 0 else 0 for i0 in range(pt)).real
            if printing:
                print('Est. mean:', tot_en_est)
        tot_en_error_true = np.sqrt(np.sum(scale_variances(variance_graph(Ptot, cct, psi), St).adj).real)
        if printing:
            print('True error:', tot_en_error_true)
        if Xt is not None:
            tot_en_error_est = np.sqrt(np.sum(scale_variances(bayes_variance_graph(Xt, cct), St).adj)).real
            if printing:
                print('Est. error:', tot_en_error_est)
        if printing:
            print()
            print()

        # electric contribution:
        elec_true = Hamiltonian_Mean(Pe, cce, psi).real
        if printing:
            print('TRUE ELECTRIC (dimless) ENERGY')
            print('True mean:', elec_true)
        if Xe is not None:
            elec_est = sum(cce[i0] * sum(Xe[i0, i0][i1, i1] * math.e**(2 * 1j * math.pi * i1 / Pe.lcm) for i1 in range(Pe.lcm)) / sum(Xe[i0, i0][i1, i1] for i1 in range(Pe.lcm)) if sum(Xe[i0, i0][i1, i1] for i1 in range(Pe.lcm)) > 0 else 0 for i0 in range(pe)).real
            if printing:
                print('Est. mean:', elec_est)
        elec_error_true = np.sqrt(np.sum(scale_variances(variance_graph(Pe, cce, psi), Se).adj).real)
        if printing:
            print('True error:', elec_error_true)
        if Xe is not None:
            elec_error_est = np.sqrt(np.sum(scale_variances(bayes_variance_graph(Xe, cce), Se).adj)).real
            if printing:
                print('Est. error:', elec_error_est)
        if printing:
            print()
            print()

        # magnetic (plaquette) contribution
        plaq_true = -Hamiltonian_Mean(Pm, ccm, psi).real / 4
        if printing:
            print('TRUE PLAQUETTE')
            print('True mean:', plaq_true)
        if Xm is not None:
            plaq_est = -sum(ccm[i0] * sum(Xm[i0, i0][i1, i1] * math.e**(2 * 1j * math.pi * i1 / Pm.lcm) for i1 in range(Pm.lcm)) / sum(Xm[i0, i0][i1, i1] for i1 in range(Pm.lcm)) if sum(Xm[i0, i0][i1, i1] for i1 in range(Pm.lcm)) > 0 else 0 for i0 in range(pm)).real / 4
            if printing:
                print('Est. mean:', plaq_est)
        plaq_error_true = np.sqrt(np.sum(scale_variances(variance_graph(Pm, ccm, psi), Sm).adj).real) / 4
        if printing:
            print('True error:', plaq_error_true)
        if Xm is not None:
            plaq_error_est = np.sqrt(np.sum(scale_variances(bayes_variance_graph(Xm, ccm), Sm).adj)).real / 4
            if printing:
                print('Est. error:', plaq_error_est)
        if printing:
            print()
            print()

        # return
        if Xt is not None:
            return ([tot_en_true, tot_en_est, tot_en_error_true, tot_en_error_est],
                    [elec_true, elec_est, elec_error_true, elec_error_est],
                    [plaq_true, plaq_est, plaq_error_true, plaq_error_est])
        else:
            return ([tot_en_true, tot_en_error_true],
                    [elec_true, elec_error_true],
                    [plaq_true, plaq_error_true])
    else:
        return ([Hamiltonian_Mean(Ptot, cct, psi).real, "None", "None", "None"],
                [Hamiltonian_Mean(Pe, cce, psi).real, "None", "None", "None"],
                [-Hamiltonian_Mean(Pm, ccm, psi).real / 4, "None", "None", "None"])
