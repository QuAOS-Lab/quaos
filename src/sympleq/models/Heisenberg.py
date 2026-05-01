import numpy as np

from sympleq.core.paulis import PauliString, PauliSum


def heisenberg_chain_hamiltonian(
    n_spins: int,
    j_x: float = 1.0,
    j_y: float = 1.0,
    j_z: float = 1.0,
    periodic: bool = False,
    all_to_all: bool = False,
) -> PauliSum:
    """
    Construct a spin-1/2 Heisenberg Hamiltonian with XX, YY, and ZZ interactions.
    """
    paulis = []
    weights = []
    dims = [2] * n_spins
    if all_to_all:
        pairs = [(i, j) for i in range(n_spins) for j in range(i + 1, n_spins)]
    else:
        pairs = [(i, i + 1) for i in range(n_spins - 1)]
    if periodic and not all_to_all and n_spins > 2:
        pairs.append((n_spins - 1, 0))

    for i, j in pairs:
        x = np.zeros(n_spins, dtype=int)
        z = np.zeros(n_spins, dtype=int)
        x[i] = 1
        x[j] = 1
        paulis.append(PauliString.from_exponents(x, z, dims))
        weights.append(j_x)

        x = np.zeros(n_spins, dtype=int)
        z = np.zeros(n_spins, dtype=int)
        x[i] = 1
        x[j] = 1
        z[i] = 1
        z[j] = 1
        paulis.append(PauliString.from_exponents(x, z, dims))
        weights.append(j_y)

        x = np.zeros(n_spins, dtype=int)
        z = np.zeros(n_spins, dtype=int)
        z[i] = 1
        z[j] = 1
        paulis.append(PauliString.from_exponents(x, z, dims))
        weights.append(j_z)

    return PauliSum.from_pauli_strings(paulis, weights=weights, phases=None)
