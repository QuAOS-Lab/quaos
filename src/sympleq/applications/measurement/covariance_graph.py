# GRAPHS

# a class for storing graphs as adjacency matrices
#     since we are dealing with covariance matrices with both vertex and edge weights,
#     this is a suitable format to capture that complexity

# Written by Andrew Jena, adapted for SympleQ by Rick Simon

import numpy as np
import itertools
import random
import networkx as nx
from sympleq.core.paulis import PauliSum


class graph:
    # Inputs:
    #     adj_mat - (numpy.array) - (weighted) adjacency matrix of graph
    #     dtype   - (numpy.dtype) - data type of graph weights
    def __init__(self, adj_mat=np.array([]), dtype=complex):
        self.adj = adj_mat.astype(dtype)

    # adds a vertex to self
    def add_vertex_(self, c=1):
        # Inputs:
        #     c - (float) - vertex weight
        if len(self.adj) == 0:
            self.adj = np.array([c])
        else:
            m0 = np.zeros((len(self.adj), 1))
            m1 = np.zeros((1, len(self.adj)))
            m2 = np.array([[c]])
            self.adj = np.block([[self.adj, m0], [m1, m2]])

    # weight a vertex
    def lade_vertex_(self, a, c):
        # Inputs:
        #     a - (int)   - vertex to be weighted
        #     c - (float) - vertex weight
        self.adj[a, a] = c

    # weight an edge
    def lade_edge_(self, a0, a1, c):
        # Inputs:
        #     a0 - (int)   - first vertex
        #     a1 - (int)   - second vertex
        #     c  - (float) - vertex weight
        self.adj[a0, a1] = c
        self.adj[a1, a0] = c

    # returns a set of the neighbors of a given vertex
    def neighbors(self, a) -> set[int]:
        # Inputs:
        #     a - (int) - vertex for which neighbors should be returned
        # Outputs:
        #     (list{int}) - set of neighbors of vertex a
        aa1 = set([])
        for i in range(self.ord()):
            if (a != i) and (self.adj[a, i] != 0):
                aa1.add(i)
        return aa1

    # returns list of all edges in self
    def edges(self) -> list[list[int]]:
        # Outputs:
        #     (list{list{int}}) - list of edges in self
        aaa = []
        for i0, i1 in itertools.combinations(range(self.ord()), 2):
            if i1 in self.neighbors(i0):
                aaa.append([i0, i1])
        return aaa

    # check whether a collection of vertices is a clique in self
    def clique(self, aa) -> bool:
        # Inputs:
        #     aa - (list{int}) - list of vertices to be checked for clique
        # Outputs:
        #     (bool) - True if aa is a clique in self; False otherwise
        for i0, i1 in itertools.combinations(aa, 2):
            if self.adj[i0, i1] == 0:
                return False
        return True

    # returns the degree of a given vertex
    def degree(self, a) -> int:
        # Inputs:
        #     a - (int) - vertex for which degree should be returned
        # Outputs:
        #     (int) - degree of vertex a
        return np.count_nonzero(self.adj[a, :])

    # returns the number of vertices in self
    def ord(self) -> int:
        # Outputs:
        #     (int) - number of vertices in self
        return self.adj.shape[0]

    # print adjacency matrix representation of self
    def print(self):
        for i0 in range(self.ord()):
            print('[', end=' ')
            for i1 in range(self.ord()):
                s = self.adj[i0, i1]
                if str(s)[0] == '-':
                    print(f'{self.adj[i0, i1]:.2f}', end=" ")
                else:
                    print(' ' + f'{self.adj[i0, i1]:.2f}', end=" ")
            print(']')

    # print self as a list of vertices together with their neighbors
    def print_neighbors(self):
        for i0 in range(self.ord()):
            print(i0, end=": ")
            for i1 in self.neighbors(i0):
                print(i1, end=" ")
            print()

    # return a deep copy of self
    def copy(self):
        # Outputs:
        #     (graph) - deep copy of self
        return graph(np.array([[self.adj[i0, i1] for i1 in range(self.ord())] for i0 in range(self.ord())]))


def nonempty_cliques(A) -> list[list[int]]:
    # Inputs:
    #     A - (graph) - graph for which all cliques should be found
    # Outputs:
    #     (list{list{int}}) - a list containing all non-empty cliques in A
    p = A.ord()
    aaa = set([frozenset([])])
    for i in range(p):
        iset = set([i])
        inter = A.neighbors(i)
        aaa |= set([frozenset(iset | (inter & aa)) for aa in aaa])
    aaa.remove(frozenset([]))
    return list([list(aa) for aa in aaa])


def all_maximal_cliques(A):
    # Inputs:
    #     A - (graph) - graph for which all cliques should be found
    # Outputs:
    #     (generator) - a generator over all maximal cliques in A
    p = A.ord()
    N = {}
    for i in range(p):
        N[i] = A.neighbors(i)
    nxG = nx.Graph()
    nxG.add_nodes_from([i for i in range(p)])
    nxG.add_edges_from([(i0, i1) for i0 in range(p) for i1 in N[i0]])
    return nx.algorithms.clique.find_cliques(nxG)


def weighted_vertex_covering_maximal_cliques(A, A1: graph | None = None, cc=None, k=1) -> list[list[int]]:
    # Inputs:
    #     A  - (graph)     - commutation graph for which covering should be found
    #     A1 - (graph)     - variance graph for which covering should be found
    #     cc - (list{int}) - coefficients of the Hamiltonian
    #     k  - (int)       - number of times each vertex should be covered
    # Outputs:
    #     (list{list{int}}) - a list containing cliques which cover A
    p = A.ord()
    if A1 is None and cc is None:
        return vertex_covering_maximal_cliques(A, k=k)
    elif A1 is None and cc is not None:
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
    elif A1 is not None and cc is None:
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
    else:
        raise NotImplementedError


def vertex_covering_maximal_cliques(A, k=1) -> list[list[int]]:
    # Inputs:
    #     A - (graph) - commutation graph for which covering should be found
    #     k - (int)   - number of times each vertex must be covered
    # Outputs:
    #     (list{list{int}}) - a list containing cliques which cover A
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


def post_process_cliques(A, aaa, k=1) -> list[list[int]]:
    # Inputs:
    #     A   - (graph)           - variance graph from which weights of cliques can be obtained
    #     aaa - (list{list{int}}) - a clique covering of the Hamiltonian
    #     k   - (int)             - number of times each vertex must be covered
    # Outputs:
    #     (list{list{int}}) - a list containing cliques which cover A
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


def commutation_graph(P: PauliSum) -> graph:
    # TODO Update with real symplectic product matrix, once mixed species supported
    spm = np.zeros([P.n_paulis(), P.n_paulis()], dtype=int)
    for i in range(P.n_paulis()):
        for j in range(i + 1, P.n_paulis()):
            P0 = P[i, :]
            P1 = P[j, :]
            xz = (P0.x_exp * P1.z_exp - P0.z_exp * P1.x_exp)
            spm[i, j] = np.sum(xz * np.array([P0.lcm] * len(P0.dimensions)) // P0.dimensions) % P0.lcm

    P_sym = spm + spm.T
    Spm_binary = P_sym.astype(bool)
    G = graph(np.array(1 - Spm_binary))
    return G


def quditwise_inner_product(PS1, PS2) -> np.bool_:
    if any(PS1.dimensions != PS2.dimensions):
        raise ValueError("Pauli strings must have the same dimensions for quditwise inner product.")
    X_1 = PS1.x_exp[:]
    X_2 = PS2.x_exp[:]
    Z_1 = PS1.z_exp[:]
    Z_2 = PS2.z_exp[:]
    dims = PS1.dimensions
    q_products = [np.sum(X_1[i] * Z_2[i] - Z_1[i] * X_2[i]) % dims[i] for i in range(len(dims))]
    return np.any(q_products)


def quditwise_commutation_graph(P) -> graph:
    p = P.n_paulis()
    adj_mat = np.array([[1 - quditwise_inner_product(P[i0, :], P[i1, :]) for i1 in range(p)] for i0 in range(p)])
    return graph(adj_mat)
