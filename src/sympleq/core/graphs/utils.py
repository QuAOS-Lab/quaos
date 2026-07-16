import numpy as np


def qudit_coupling_graph(F: np.ndarray) -> dict[int, set[int]]:
    """
    Construct the qudit-coupling graph of a symplectic matrix F.

    Conventions:
      - F is a 2n x 2n integer matrix (symplectic over Z_d, dimension d implicit).
      - Tableau rows are row vectors p in Z_d^{2n}.
      - Clifford action is p' = p @ F.T.
      - Canonical generators:
          e_i^X = [e_i | 0]  -> image is column i of F
          e_i^Z = [0 | e_i]  -> image is column n+i of F
      - Support of a tableau row p = [x | z] is the set of qudits j
        such that (x_j, z_j) != (0, 0).

    The coupling graph has:
      - vertices: {0, 1, ..., n-1} (Python indexing for qudits 1..n),
      - an undirected edge {i, j} if the image of X_i or Z_i
        has nonzero support on qudit j (or vice versa).

    Returns:
      A dict mapping each qudit index i to a set of neighbouring qudits.
    """
    F = np.asarray(F)
    if F.ndim != 2 or F.shape[0] != F.shape[1]:
        raise ValueError("F must be a square 2D array.")

    m = F.shape[0]
    if m % 2 != 0:
        raise ValueError("F must have even dimension 2n x 2n.")
    n = m // 2

    # Adjacency as an undirected graph: neighbors[i] is a set of qudits j.
    neighbors: dict[int, set[int]] = {i: set() for i in range(n)}

    # Work over integers mod d in practice; here we just test for != 0.
    # For mixed-prime dims you’d reduce each column mod its local d_i first;
    # structurally, the logic is the same.
    for i in range(n):
        # Columns corresponding to X_i and Z_i images
        for col_idx in (i, n + i):
            col = F[:, col_idx]           # length 2n
            x_part = col[:n]
            z_part = col[n:]

            # Support: qudits where (x_j, z_j) != (0, 0)
            for j in range(n):
                if x_part[j] != 0 or z_part[j] != 0:
                    if j != i:
                        neighbors[i].add(j)
                        neighbors[j].add(i)  # ensure undirected

    return neighbors
