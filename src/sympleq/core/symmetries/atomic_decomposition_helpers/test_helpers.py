import numpy as np


def assert_certified_decomposition(F, p, Omega, blocks, T_global):
    """
    blocks: list of AtomicBlock (each has T_blk_amb)
    T_global: concatenation of block bases (columns), shape (2n,2n)
    """

    def mod_p(A):
        return A % p

    def rank_mod(A):
        # replace with your existing rank_mod for speed
        A = mod_p(A).astype(np.int64)
        # naive row-reduction placeholder:
        import sympy as sp
        return sp.Matrix(A.tolist()).rank()

    n2 = F.shape[0]
    assert Omega.shape == (n2, n2)

    # 1) Global basis is full rank and symplectic
    assert T_global.shape == (n2, n2)
    assert rank_mod(T_global) == n2, "T_global not invertible"
    G = mod_p(T_global.T @ Omega @ T_global)
    assert np.array_equal(G, mod_p(Omega)), "T_global not symplectic"

    # 2) Blocks span and are direct sum
    Tcat = np.concatenate([b.T_blk_amb for b in blocks], axis=1) if blocks else np.zeros((n2, 0), dtype=np.int64)
    assert Tcat.shape[1] == n2, "Block dims do not sum to full space"
    assert rank_mod(Tcat) == n2, "Blocks do not span full space"

    # 3) Each block is symplectic nondegenerate and F-invariant
    for i, b in enumerate(blocks):
        T = b.T_blk_amb
        k2 = T.shape[1]
        assert k2 % 2 == 0

        # nondegenerate restriction
        Om_blk = mod_p(T.T @ Omega @ T)
        assert rank_mod(Om_blk) == k2, f"Block {i} is degenerate"

        # invariance: span(F T) ⊆ span(T)
        FT = mod_p(F @ T)
        r1 = rank_mod(T)
        r2 = rank_mod(np.concatenate([T, FT], axis=1))
        assert r2 == r1, f"Block {i} not F-invariant"

    # 4) Pairwise symplectic orthogonality
    for i in range(len(blocks)):
        Ti = blocks[i].T_blk_amb
        for j in range(i + 1, len(blocks)):
            Tj = blocks[j].T_blk_amb
            X = mod_p(Ti.T @ Omega @ Tj)
            assert np.all(X == 0), f"Blocks {i} and {j} not symplectically orthogonal"
