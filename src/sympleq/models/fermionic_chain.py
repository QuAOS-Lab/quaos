import numpy as np
from typing import Tuple
from sympleq.core.paulis import PauliSum


def _fermionic_chain_hamiltonian(
    n_qubits: int,
    J: float,
    V: float,
    D_vec: np.ndarray,
    periodic: bool = False,
    include_identity_shift: bool = False,
    tol: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build a Pauli-tableau representation of

        H = -J Σ_j (c_j^† c_{j+1} + h.c.) + Σ_j D_j n_j + V Σ_j n_j n_{j+1}

    using the identifications (as requested):
        c_j = (X_j + i Y_j)/2,  n_j = (Z_j + I)/2.

    Output:
      tableau : (M, 2n) uint8 array. Row m is [x_0..x_{n-1} | z_0..z_{n-1}]
      coeffs  : (M,) float array (real coefficients multiplying each Pauli operator)
      phases  : (M,) int array in {0,1,2,3} giving an overall prefactor (1j)**phase

    Convention for each row:
        P = (1j)**phase * Π_k X_k^{x_k} Z_k^{z_k}

    With this convention, Y = i X Z corresponds to (x=1,z=1,phase=1).

    Notes:
      * The hopping term simplifies to:  -(J/2) Σ_j (X_j X_{j+1} + Y_j Y_{j+1})
      * Terms proportional to identity (energy shifts) are included only if
        include_identity_shift=True.
    """
    n = int(n_qubits)
    if n <= 0:
        raise ValueError("n_qubits must be positive.")
    D_vec = np.asarray(D_vec, dtype=float).reshape(-1)
    if D_vec.shape[0] != n:
        raise ValueError(f"D_vec must have length n_qubits={n}, got {D_vec.shape[0]}.")

    # key: (phase, x_tuple, z_tuple) -> coefficient
    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], float] = {}

    def add_term(x: np.ndarray, z: np.ndarray, phase: int, coeff: float) -> None:
        phase &= 3
        key = (phase, tuple(int(v) for v in x), tuple(int(v) for v in z))
        terms[key] = terms.get(key, 0.0) + float(coeff)

    def zeros():
        return np.zeros(n, dtype=np.uint8), np.zeros(n, dtype=np.uint8)

    # neighbor list
    n_edges = n if periodic else (n - 1)

    # 1) Hopping:  -(J/2) Σ (XX + YY)
    for j in range(n_edges):
        k = (j + 1) % n

        x, z = zeros()
        x[j] = 1
        x[k] = 1
        add_term(x, z, phase=0, coeff=-(J / 2.0))

        # YY = (i^2) * (XZ ⊗ XZ) in our (phase, x|z) convention
        x, z = zeros()
        x[j] = 1
        z[j] = 1
        x[k] = 1
        z[k] = 1
        add_term(x, z, phase=2, coeff=-(J / 2.0))

    # 2) On-site: Σ D_j n_j = Σ (D_j/2) Z_j + (D_j/2) I
    for j, Dj in enumerate(D_vec):
        if Dj == 0:
            continue
        x, z = zeros()
        z[j] = 1
        add_term(x, z, phase=0, coeff=(Dj / 2.0))

        if include_identity_shift:
            x, z = zeros()
            add_term(x, z, phase=0, coeff=(Dj / 2.0))

    # 3) Interaction: V Σ n_j n_{j+1} = (V/4) Σ (ZZ + Z_j + Z_{j+1} + I)
    if V != 0:
        for j in range(n_edges):
            k = (j + 1) % n

            # ZZ
            x, z = zeros()
            z[j] = 1
            z[k] = 1
            add_term(x, z, phase=0, coeff=(V / 4.0))

            # Z_j
            x, z = zeros()
            z[j] = 1
            add_term(x, z, phase=0, coeff=(V / 4.0))

            # Z_{j+1}
            x, z = zeros()
            z[k] = 1
            add_term(x, z, phase=0, coeff=(V / 4.0))

            # I
            if include_identity_shift:
                x, z = zeros()
                add_term(x, z, phase=0, coeff=(V / 4.0))

    # prune zeros / tiny terms
    items = []
    for key, c in terms.items():
        if tol > 0.0 and abs(c) <= tol:
            continue
        if c != 0.0:
            items.append((key, c))

    # deterministic order: by phase, then x bits, then z bits
    items.sort(key=lambda kv: (kv[0][0], kv[0][1], kv[0][2]))

    M = len(items)
    tableau = np.zeros((M, 2 * n), dtype=np.uint8)
    coeffs = np.zeros(M, dtype=float)
    phases = np.zeros(M, dtype=np.int8)

    for i, ((phase, x_t, z_t), c) in enumerate(items):
        x = np.fromiter(x_t, count=n, dtype=np.uint8)
        z = np.fromiter(z_t, count=n, dtype=np.uint8)
        tableau[i, :n] = x
        tableau[i, n:] = z
        coeffs[i] = c
        phases[i] = phase

    return tableau, coeffs, phases


def fermionic_chain_hamiltonian(n: int, J: float, V: float, D_vec: np.ndarray, periodic: bool = False) -> 'PauliSum':
    tableau, coeffs, phases = _fermionic_chain_hamiltonian(n, J, V, D_vec, periodic)
    return PauliSum.from_tableau(tableau, weights=coeffs, phases=phases)
