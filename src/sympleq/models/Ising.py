from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits import Gate
import numpy as np


def ising_chain_hamiltonian(n_spins: int,
                            J_zz: float | np.ndarray | list,
                            h_x: float | np.ndarray | list,
                            periodic: bool = False):
    """
    Constructs the Hamiltonian of the 1D Ising model in a transverse field.

    Parameters
    ----------
    n_spins : int
        The number of spins in the chain.
    J_zz : float or np.ndarray or list
        The Ising interaction strength between nearest-neighbour spins.
        If this is a matrix, it will be used as the interaction strength matrix
        should only have nearest neighbour terms
    h_x : float or np.ndarray or list
        The strength of the transverse field.
    periodic : bool, optional
        Whether the chain is periodic (default: False).

    Returns
    -------
    PauliSum
        The Hamiltonian as a PauliSum object.
    """
    paulis: list[PauliString] = []
    weights = []
    dims = [2 for _ in range(n_spins)]

    if isinstance(J_zz, np.ndarray) or isinstance(J_zz, list):
        J_zz = np.array(J_zz)  # Convert to numpy array if list
        if J_zz.shape != (n_spins,) and periodic:
            raise ValueError("J_zz must be a vector of size n_spins for periodic chain")
        elif J_zz.shape != (n_spins - 1,) and not periodic:
            raise ValueError("J_zz must be a vector of size n_spins-1 for closed chain")
        else:
            raise ValueError(f"J_zz shape, {J_zz.shape}, not compatible with boundary conditions."
                             "Should be (n_spins) for periodic chain or (n_spins-1) for closed chain.")
    elif isinstance(J_zz, float) or isinstance(J_zz, int):
        J_zz = J_zz * np.ones(n_spins if periodic else n_spins - 1)
    else:
        raise ValueError("J_zz must be a float or numpy array")

    if isinstance(h_x, np.ndarray):
        h_x = np.array(h_x)  # Convert to numpy array if list
        if h_x.shape != (n_spins,):
            raise ValueError("h_x must be a vector of size n_spins")
    elif isinstance(h_x, float) or isinstance(h_x, int):
        h_x = h_x * np.ones(n_spins)
    else:
        raise ValueError("h_x must be a float or numpy array")

    # this bit is just for typing
    if not isinstance(J_zz, np.ndarray):
        raise Exception("J_zz must be a numpy array")

    if not isinstance(h_x, np.ndarray):
        raise Exception("h_x must be a numpy array")

    # ZZ terms
    for i in range(n_spins - 1):
        zz = np.zeros(n_spins, dtype=int)
        zz[i] = 1
        zz[i + 1] = 1
        paulis.append(PauliString.from_exponents(np.zeros(n_spins, dtype=int), zz, dims))
        weights.append(J_zz[i])

    # Periodic ZZ term (last ↔ first spin)
    if periodic and n_spins > 2:
        zz = np.zeros(n_spins, dtype=int)
        zz[0] = 1
        zz[-1] = 1
        paulis.append(PauliString.from_exponents(np.zeros(n_spins, dtype=int), zz, dims))
        weights.append(J_zz[n_spins])

    # X terms (transverse field)
    for i in range(n_spins):
        x = np.zeros(n_spins, dtype=int)
        x[i] = 1
        paulis.append(PauliString.from_exponents(x, np.zeros(n_spins, dtype=int), dims))
        weights.append(h_x[i])

    return PauliSum.from_pauli_strings(paulis, weights=weights, phases=None)


def ising_2d_hamiltonian(n_x: int, n_y: int, J_zz: float, h_x: float, periodic: bool = False) -> PauliSum:
    """
    Constructs the Hamiltonian of a 2D Ising model with nearest-neighbor interactions
    and a transverse field.

    Parameters
    ----------
    n_x, n_y : int
        The number of spins in the x- and y-directions, respectively.
    J_zz : float
        The strength of the nearest-neighbor interactions.
    h_x : float
        The strength of the transverse field.
    periodic : bool, optional
        Whether the chain is periodic in both x- and y-directions (default: False).

    Returns
    -------
    PauliSum
        The Hamiltonian as a PauliSum object.
    """
    paulis = []
    weights = []
    n_spins = n_x * n_y
    dims = [2 for _ in range(n_spins)]

    def site_index(x, y):
        """Map 2D coordinates to 1D index in row-major order."""
        return y * n_x + x

    # ZZ terms (horizontal + vertical couplings)
    for x in range(n_x):
        for y in range(n_y):
            i = site_index(x, y)

            # horizontal coupling (x → x+1)
            if x + 1 < n_x or periodic:
                j = site_index((x + 1) % n_x, y)
                zz = np.zeros(n_spins, dtype=int)
                zz[i] = 1
                zz[j] = 1
                paulis.append(PauliString.from_exponents(np.zeros(n_spins, dtype=int), zz, dims))
                weights.append(J_zz)

            # vertical coupling (y → y+1)
            if y + 1 < n_y or periodic:
                j = site_index(x, (y + 1) % n_y)
                zz = np.zeros(n_spins, dtype=int)
                zz[i] = 1
                zz[j] = 1
                paulis.append(PauliString.from_exponents(np.zeros(n_spins, dtype=int), zz, dims))
                weights.append(J_zz)

    # X terms (transverse field)
    for i in range(n_spins):
        x = np.zeros(n_spins, dtype=int)
        x[i] = 1
        paulis.append(PauliString.from_exponents(x, np.zeros(n_spins, dtype=int), dims))
        weights.append(h_x)

    return PauliSum.from_pauli_strings(paulis, weights=weights, phases=None)


def ising_lower_triangular_hamiltonian(L: int, J_zz: float, h_x: float) -> PauliSum:
    """
    Construct a transverse-field Ising Hamiltonian on the lower-triangular half
    of an L x L square lattice (including the diagonal).

    Geometry:
      - Sites are integer coordinates (x, y) with 0 <= y <= x < L.
      - Nearest-neighbor ZZ couplings are added along +x and +y lattice edges
        whenever both endpoints are in the triangular region.
      - Uniform transverse X field is added on every site.

    Parameters
    ----------
    L : int
        Side length of the parent square. Must satisfy L >= 1.
    J_zz : float
        Coupling strength of nearest-neighbor ZZ terms.
    h_x : float
        Transverse-field strength for single-site X terms.

    Returns
    -------
    PauliSum
        The Hamiltonian as a PauliSum object on N = L*(L+1)/2 qubits.
    """
    if int(L) < 1:
        raise ValueError("L must be at least 1.")

    # Deterministic indexing of triangular sites in row-major order by x then y.
    coords: list[tuple[int, int]] = []
    for x in range(int(L)):
        for y in range(x + 1):
            coords.append((x, y))

    site_to_idx = {c: i for i, c in enumerate(coords)}
    n_spins = len(coords)
    dims = [2 for _ in range(n_spins)]

    paulis: list[PauliString] = []
    weights: list[float] = []
    z0 = np.zeros(n_spins, dtype=int)
    x0 = np.zeros(n_spins, dtype=int)

    # ZZ terms along lattice edges internal to the triangular domain.
    for (x, y), i in site_to_idx.items():
        for xn, yn in ((x + 1, y), (x, y + 1)):
            j = site_to_idx.get((xn, yn))
            if j is None:
                continue
            zz = z0.copy()
            zz[i] = 1
            zz[j] = 1
            paulis.append(PauliString.from_exponents(z0, zz, dims))
            weights.append(J_zz)

    # X terms (transverse field)
    for i in range(n_spins):
        x = x0.copy()
        x[i] = 1
        paulis.append(PauliString.from_exponents(x, z0, dims))
        weights.append(h_x)

    return PauliSum.from_pauli_strings(paulis, weights=weights, phases=None)


def heuristic_clifford_symmetry(n_spins: int, periodic: bool = False) -> Gate:
    A = np.zeros((n_spins, n_spins), dtype=int)
    B = np.ones((n_spins, n_spins), dtype=int)
    C = np.zeros((n_spins, n_spins), dtype=int)

    if periodic:
        A[0, 1] = 1
        A[1, 0] = 1
        for i in range(n_spins - 2):
            A[-1 - i, 2 + i] = 1
    else:
        # ones on anti-diagonal
        for i in range(n_spins):
            A[-i - 1, i] = 1

    F = np.block([[A, B], [C, A]])
    F_G = Gate('F', F,
               np.concatenate([np.zeros(n_spins, dtype=int), np.ones(n_spins, dtype=int)]))
    return F_G


def product_state_ising(
    N: int,
    kind: str = "x_plus",
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build a simple product state on an N-site qubit chain.

    Returns
    -------
    psi : (2**N,) complex ndarray
        Statevector in computational basis |0...0>,|0...1>,...,|1...1>.
    dims : (N,) int ndarray
        Dimensions (all 2's).
    """
    dims = np.full(N, 2, dtype=int)

    # Single-qubit kets
    ket0 = np.array([1.0, 0.0], dtype=np.complex128)  # |0> ~ |↑_z>
    ket1 = np.array([0.0, 1.0], dtype=np.complex128)  # |1> ~ |↓_z>
    ket_plus = (ket0 + ket1) / np.sqrt(2.0)
    ket_minus = (ket0 - ket1) / np.sqrt(2.0)

    if kind == "z_up":
        locals_ = [ket0 for _ in range(N)]
    elif kind == "z_neel":
        locals_ = [ket0 if i % 2 == 0 else ket1 for i in range(N)]
    elif kind == "x_plus":
        locals_ = [ket_plus for _ in range(N)]
    elif kind == "x_minus":
        locals_ = [ket_minus for _ in range(N)]
    else:
        raise ValueError(f"Unknown product_state kind '{kind}'.")

    psi = locals_[0]
    for k in range(1, N):
        psi = np.kron(psi, locals_[k])

    return psi, dims
