import numpy as np
import itertools
import re
import math 
    

# PAULIS

# a class for storing sets of Pauli operators as pairs of symplectic matrices
class pauli:
    def __init__(self,X,Z,dims=2,phases=0):
        # Inputs:
        #     X - (numpy.array) - X-part of Pauli in symplectic form with shape (p,q)
        #     Z - (numpy.array) - Z-part of Pauli in symplectic form with shape (p,q)
        if X.shape != Z.shape:
            raise Exception("X- and Z-parts must have same shape")
        if type(dims) == int:
            dims = np.array([dims]*X.shape[1])
        elif type(dims) == list:
            dims = np.array(dims)
        if type(phases) == int:
            phases = np.array([phases]*X.shape[0])
        elif type(phases) == list:
            phases = np.array(phases)
        self.X = X%dims
        self.Z = Z%dims
        self.lcm = np.lcm.reduce(dims)
        self.dims = dims
        self.phases = phases%self.lcm

    # check whether self has only X component
    def is_IX(self):
        # Outputs:
        #     (bool) - True if self has only X componenet, False otherwise
        return not np.any(self.Z)

    # check whether self has only Z component 
    def is_IZ(self):
        # Outputs:
        #     (bool) - True if self has only Z componenet, False otherwise
        return not np.any(self.X)

    # check whether the set of Paulis are pairwise commuting
    def is_commuting(self):
        # Outputs:
        #     (bool) - True if self is pairwise commuting set of Paulis
        p = self.paulis()
        PP = [self.a_pauli(i) for i in range(p)]
        return not any(symplectic_inner_product(PP[i0],PP[i1]) for i0,i1 in itertools.combinations(range(p),2))

    # check whether the set of Paulis are pairwise commuting on every qudit
    def is_quditwise_commuting(self):
        # Outputs:
        #     (bool) - True if self is pairwise quditwise commuting set of Paulis
        p = self.paulis()
        PP = [self.a_pauli(i) for i in range(p)]
        return not any(quditwise_inner_product(PP[i0],PP[i1]) for i0,i1 in itertools.combinations(range(p),2))

    # pull out the ath Pauli from self
    def a_pauli(self,a):
        # Inputs: 
        #     a - (int) - index of Pauli to be returned
        # Outputs:
        #     (pauli) - the ath Pauli in self
        return pauli(np.array([self.X[a,:]]),np.array([self.Z[a,:]]),self.dims,np.array([self.phases[a]]))

    # count the number of Paulis in self
    def paulis(self):
        # Output: (int)
        return self.X.shape[0]

    # count the number of qudits in self
    def qudits(self):
        # Outputs: (int)
        return self.X.shape[1]

    # delete Paulis indexed by aa
    def delete_paulis_(self,aa):
        # Inputs: 
        #     aa - (list of int)
        if type(aa) is int:
            self.X = np.delete(self.X,aa,axis=0)
            self.Z = np.delete(self.Z,aa,axis=0)
            self.phases = np.delete(self.phases,aa)
        else:
            for a in sorted(aa,reverse=True):
                self.X = np.delete(self.X,a,axis=0)
                self.Z = np.delete(self.Z,a,axis=0)
                self.phases = np.delete(self.phases,a)
        return self

    # return self after deletion of qudits indexed by aa
    def delete_qudits_(self,aa):
        # Inputs: 
        #     aa - (list of int)
        if type(aa) is int:
            self.X = np.delete(self.X,aa,axis=1)
            self.Z = np.delete(self.Z,aa,axis=1)
            self.dims = np.delete(self.dims,aa)
        else:
            for a in sorted(aa,reverse=True):
                self.X = np.delete(self.X,a,axis=1)
                self.Z = np.delete(self.Z,a,axis=1)
                self.dims = np.delete(self.dims,a)
        self.lcm = np.lcm.reduce(dims)

    # return deep copy of self
    def copy(self):
        # Outputs: (pauli)
        X = np.array([[self.X[i0,i1] for i1 in range(self.qudits())] for i0 in range(self.paulis())])
        Z = np.array([[self.Z[i0,i1] for i1 in range(self.qudits())] for i0 in range(self.paulis())])
        dims = np.array([self.dims[i] for i in range(self.qudits())])
        phases = np.array([self.phases[i] for i in range(self.paulis())])
        return pauli(X,Z,dims,phases)

    # # print string representation of self
    def print(self):
        sss,dims,phases = pauli_to_string(self)
        for ss in sss:
            print(ss)

    # print symplectic representation of self
    def print_symplectic(self):
        print(''.join(str(int(i)) for i in self.dims),''.join(str(int(i)) for i in self.dims))
        print('-'*self.qudits(),'-'*self.qudits())
        for i in range(self.paulis()):
            print(''.join(str(int(i1)) for i1 in self.X[i,:])+' '+''.join(str(int(i1)) for i1 in self.Z[i,:])+' '+str(self.phases[i])+'/'+str(self.lcm))


# convert a pauli object to its matrix representation
def pauli_to_matrix(P):
    # Inputs:
    #     P - (pauli) - must have shape (1,q)
    # Outputs:
    #     (scipy.sparse.csr_matrix) - matrix representation of input Pauli
    if P.paulis() != 1:
        raise Exception("Matrix can only be constructed for a single Pauli")
    X,Z,dims,phase = P.X[0],P.Z[0],P.dims,P.phases[0]
    return math.e**(phase*2*math.pi*1j/P.lcm)*tensor([XZ_mat(dims[i],X[i],Z[i]) for i in range(P.qudits())])

# convert a collection of strings (or single string) to a pauli object
def string_to_pauli(sss,dims=2,phases=0):
    # Inputs:
    #     sss - (list{str}) or (str) - string representation of Pauli
    # Outputs:
    #     (pauli) - Pauli corresponding to input string(s)
    if type(sss) is str:
        X = np.array([[re.split("x|z",s)[1] for s in sss.split()]],dtype=int)
        Z = np.array([[re.split("x|z",s)[2] for s in sss.split()]],dtype=int)
        return pauli(X,Z,dims,phases)
    else:
        X = np.array([[re.split('x|z',s)[1] for s in ss.split()] for ss in sss],dtype=int)
        Z = np.array([[re.split('x|z',s)[2] for s in ss.split()] for ss in sss],dtype=int)
        return pauli(X,Z,dims,phases)

# convert a pauli object to a collection of strings (or single string)
def pauli_to_string(P):
    # Inputs:
    #     P - (pauli) - Pauli to be stringified
    # Outputs:
    #     (str) or (list{str}) - string representation of Pauli
    X,Z = P.X,P.Z
    return [''.join('x'+str(X[i0,i1])+'z'+str(Z[i0,i1])+' ' for i1 in range(P.qudits()))[:-1] for i0 in range(P.paulis())],P.dims,P.phases

# the symplectic inner product of two pauli objects (each with a single Pauli)
def symplectic_inner_product(P0,P1):
    # Inputs:
    #     P0 - (pauli) - must have shape (1,q)
    #     P1 - (pauli) - must have shape (1,q)
    # Outputs:
    #     (bool) - symplectic inner product of Paulis
    if (P0.paulis() != 1) or (P1.paulis() != 1):
        raise Exception("Symplectic inner product only works with pair of single Paulis")
    if any(P0.dims-P1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    return bool(np.sum((P0.X*P1.Z-P0.Z*P1.X)*np.array([P0.lcm]*len(P0.dims))//P0.dims)%P0.lcm)

# the symplectic inner product of two pauli objects (each with a single Pauli)
def quditwise_inner_product(P0,P1):
    # Inputs:
    #     P0 - (pauli) - must have shape (1,q)
    #     P1 - (pauli) - must have shape (1,q)
    # Outputs:
    #     (bool) - quditwise inner product of Paulis
    if (P0.paulis() != 1) or (P1.paulis() != 1):
        raise Exception("Qubitwise inner product only works with pair of single Paulis")
    if any(P0.dims-P1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    return any(np.sum(P0.X[0,i]*P1.Z[0,i]-P0.Z[0,i]*P1.X[0,i])%P0.dims[i] for i in range(P0.qudits()))

# the product of two pauli objects
def pauli_product(P0,P1):
    # Inputs:
    #     P0 - (pauli) - must have shape (1,q)
    #     P1 - (pauli) - must have shape (1,q)
    # Outputs:
    #     (pauli) - product of Paulis
    if P0.paulis() != 1 or P1.paulis() != 1:
        raise Exception("Product can only be calculated for single Paulis")
    if any(P0.dims-P1.dims):
        raise Exception("Symplectic inner product only works if Paulis have same dimensions")
    return pauli(P0.X+P1.X,P0.Z+P1.Z,P0.dims,P0.phases+P1.phases)


def tensor(mm):
    # Inputs:
    #     mm - (list{scipy.sparse.csr_matrix}) - matrices to tensor
    # Outputs:
    #     (scipy.sparse.csr_matrix) - tensor product of matrices
    if len(mm) == 0:
        return matrix([])
    elif len(mm) == 1:
        return mm[0]
    else:
        return scipy.sparse.kron(mm[0],tensor(mm[1:]),format="csr")
