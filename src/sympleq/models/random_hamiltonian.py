import random
import numpy as np
from sympleq.core.paulis import PauliSum, PauliString
from sympleq.core.circuits import Gate
from sympleq.utils import int_to_bases, bases_to_int


def random_pauli_hamiltonian(num_paulis, qudit_dims, mode='rand'):
    """
    Generates a random Pauli Hamiltonian with the given number of Pauli operators plus their Hermitian conjugate pairs.

    Parameters:
        num_paulis (int): Number of Pauli operators to generate.
        qudit_dims (list): List of dimensions for each qudit.
        mode (str): 'rand' or 'uniform' - dictates form of weights in the PauliSum

    Returns:
        tuple: A random PauliSum
    """

    '''
    n_qudits = len(qudit_dims)
    pauli_strings = []
    coefficients = []

    for p in range(num_paulis):
        # np.random.randint(qudit_dims, size=n_qudits)
        x_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]
        # np.random.randint(qudit_dims, size=n_qudits)
        z_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]
        while np.all(np.array(x_exp) == 0) and np.all(np.array(z_exp) == 0):
            x_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]
            z_exp = [random.randint(0, qudit_dims[i] - 1) for i in range(n_qudits)]

        x_exp_H = np.zeros_like(x_exp)
        z_exp_H = np.zeros_like(z_exp)
        phase_factor = 1
        pauli_str = ''
        pauli_str_H = ''
    '''
    q2 = np.repeat(qudit_dims, 2)
    available_paulis = list(np.arange(int(np.prod(q2))))

    pauli_strings = []
    coefficients = []

    for _ in range(num_paulis):
        pauli_index = random.choice(available_paulis)
        available_paulis.remove(pauli_index)

        exponents = int_to_bases(pauli_index, q2)
        exponents_H = np.zeros_like(exponents)
        phase_factor = 1
        pauli_str = ' '
        pauli_str_H = ' '

        for j in range(len(qudit_dims)):
            r, s = int(exponents[2 * j]), int(exponents[2 * j + 1])
            pauli_str += f"x{r}z{s} "
            exponents_H[2 * j] = (-r) % qudit_dims[j]
            exponents_H[2 * j + 1] = (-s) % qudit_dims[j]
            pauli_str_H += f"x{exponents_H[2 * j]}z{exponents_H[2 * j + 1]} "

            omega = np.exp(2 * np.pi * 1j / qudit_dims[j])
            phase_factor *= omega**(r * s)

        pauli_strings.append(PauliString.from_string(pauli_str.strip(), dimensions=qudit_dims))
        if mode == 'rand' or mode == 'random':
            coeff = np.random.normal(0, 1) + 1j * np.random.normal(0, 1)
        elif mode == 'uniform' or mode == 'one':
            coeff = 1 + 0 * 1j
        elif mode[0:7] == 'randint':
            # mode is 'randint2', 'randint3', etc.
            d = int(mode[7:])
            coeff = np.random.randint(1, d + 1)

        if not np.array_equal(exponents, exponents_H):
            # random string not Hermitian, add conjugate pair
            conjugate_index = bases_to_int(exponents_H, q2)
            coefficients.append(coeff)
            coefficients.append(np.conj(coeff) * phase_factor)
            available_paulis.remove(conjugate_index)
            pauli_strings.append(PauliString.from_string(pauli_str_H.strip(), dimensions=qudit_dims))
        else:
            coefficients.append(coeff.real)

    rand_ham = PauliSum.from_pauli_strings(pauli_strings, weights=coefficients)
    # rand_ham.combine_equivalent_paulis()
    return rand_ham


def random_pauli_symmetry_hamiltonian(n_qudits: int, n_paulis: int, n_redundant=0,
                                      n_conditional=0, weight_mode='uniform', phase_mode='zero', shuffle=True):
    # 0: I, 1: X, 2: Z, 3: Y
    """
    Generate a random Pauli Hamiltonian with n_qudits qudits and n_paulis Pauli strings,
    with n_redundant redundant qubits and n_conditional conditional qubits.

    The Pauli strings are chosen randomly from the set of strings with the given
    structure. The weights are chosen randomly from the set of real numbers.

    Parameters
    ----------
    n_qudits : int
        The number of qudits in the Hamiltonian.
    n_paulis : int
        The number of Pauli strings in the Hamiltonian.
    n_redundant : int, optional
        The number of redundant qubits. Default is 0.
    n_conditional : int, optional
        The number of conditional qubits. Default is 0.
    weight_mode : str, optional
        The mode of the weights. Can be 'uniform' (default) or 'random'.
    phase_mode : str, optional
        The mode of the phases. Can be 'zero' (default) or 'random'.
    Returns
    -------
    P : PauliSum
        The random Pauli Hamiltonian.

    Examples
    --------
    >>> from sympleq.models.random_hamiltonian import random_pauli_symmetry_hamiltonian
    >>> random_pauli_symmetry_hamiltonian(2, 4)
    PauliSum of size 4x2 with 4 terms and 0 redundant or conditional qubits.
    """
    # TODO: Implementation for Qudits
    # TODO: Make sure that remaining paulis are always unique
    n_rest = n_qudits - n_redundant - n_conditional
    if n_paulis < 2 * n_rest:
        raise ValueError('Too few paulis for full basis with this number of independent qubits')
    # create general structure of the Pauli Hamiltonian before scrambling it with clifford gates
    # redundant qubits
    P = np.zeros((n_paulis, n_qudits), dtype=int)

    # conditional qubits
    for i in range(n_conditional):
        q = np.arange(n_redundant, n_redundant + n_conditional)[i]
        P[2 * n_rest + i, q] = 2  # Z
        for j in range(1, n_paulis - (2 * n_rest + i)):
            P[2 * n_rest + i + j, q] = np.random.choice([0, 2])  # I, Z

    # remaining qubits
    for i in range(n_qudits - (n_redundant + n_conditional)):
        q = np.arange(n_redundant + n_conditional, n_qudits)[i]
        P[i, q] = 1
        P[n_rest + i, q] = 2
        for j in range(1, n_rest - i):
            P[i + j, q] = np.random.choice([0, 1, 2, 3])
            P[n_rest + i + j, q] = np.random.choice([0, 1, 2, 3])

        for j in range(n_paulis - 2 * n_rest):
            P[2 * n_rest + j, q] = np.random.choice([0, 1, 2, 3])

    # Turn P-Matrix into PauliSum
    pauli_strings = ['' for _ in range(n_paulis)]
    for i in range(n_paulis):
        for j in range(n_qudits):
            if P[i, j] == 0:
                pauli_strings[i] += 'x0z0'
            elif P[i, j] == 1:
                pauli_strings[i] += 'x1z0'
            elif P[i, j] == 2:
                pauli_strings[i] += 'x0z1'
            elif P[i, j] == 3:
                pauli_strings[i] += 'x1z1'

            pauli_strings[i] += ' '
        pauli_strings[i] = pauli_strings[i].strip()

    if weight_mode == 'uniform':
        weights = np.ones(n_paulis, dtype=float)
    elif weight_mode == 'random':
        weights = np.random.rand(n_paulis)

    if phase_mode == 'zero':
        phases = np.zeros(n_paulis, dtype=int)
    elif phase_mode == 'random':
        phases = np.random.randint(0, 2, size=n_paulis, dtype=int)

    P = PauliSum.from_string(pauli_strings, dimensions=[2] * n_qudits, weights=weights, phases=phases)

    if shuffle:
        g = Gate.from_random(n_qudits, 2)
        P = g.act(P, tuple(np.arange(n_qudits)))

    return P


# def random_gate_symmetric_hamiltonian(G: Gate,
#                                       dimension: int,
#                                       qudit_indices: tuple[int, ...] | list[int],
#                                       n_qudits: int,
#                                       n_paulis: int | None = None,
#                                       weight_mode: str = 'uniform',
#                                       scrambled: bool = False):
#     """
#     Generate a random symmetric Hamiltonian from a gate G.

#     Parameters
#     ----------
#     G : Gate
#         The gate for which to generate the symmetric Hamiltonian.
#     n_qudits : int
#         The number of qudits in the resulting Hamiltonian. If None, it is set to G.dimension + 1.
#     n_paulis : int
#         The number of Pauli strings in the resulting Hamiltonian. If None, it is set to 2 * n_qudits.
#     weight_mode : str
#         Whether to use 'uniform' or 'random' weights in the Hamiltonian.

#     Returns
#     -------
#     P_sym : PauliSum
#         The symmetric Hamiltonian as a PauliSum.

#     Notes
#     -----
#     This function samples random Pauli strings and closes each one under the orbit of G, producing a sum that is
#     exactly invariant under G without needing to know the gate's global order. The weights are rounded to 10 decimal
#     places and Pauli strings with zero weight are removed.
#     """

#     if n_paulis is None:
#         n_paulis = 2 * n_qudits

#     if isinstance(qudit_indices, list):
#         qudit_indices = tuple(qudit_indices)

#     # Build full dimensions for the Hamiltonian, embedding the gate dimensions on the target indices.
#     dims = np.full(n_qudits, dimension, dtype=int)

#     rng = np.random.default_rng()

#     def _new_weight() -> float:
#         return float(rng.random()) if weight_mode == 'random' else 1.0

#     def _orbit(seed: PauliSum) -> list[PauliSum]:
#         """
#         Close the orbit of `seed` under G. Using the first repeat of (tableau, weight) as the stopping condition
#         guarantees an algebraically closed set without needing the global gate order.
#         """
#         seen = set()
#         orbit_terms = []
#         term = seed
#         while True:
#             key = (tuple(term.tableau[0]), complex(np.around(term.weights[0], decimals=12)))
#             if key in seen:
#                 break
#             seen.add(key)
#             orbit_terms.append(term)
#             term = G.act(term, qudit_indices).to_standard_form()
#         return orbit_terms

#     # Accumulate orbit-closed terms until we have at least n_paulis *non-zero* unique tableau rows.
#     # We combine coefficients on-the-fly to avoid producing an empty Hamiltonian due to cancellations.
#     coeff_by_row: dict[bytes, complex] = {}
#     row_by_key: dict[bytes, np.ndarray] = {}

#     def _add_coeff(row: np.ndarray, coeff: complex) -> None:
#         k = np.asarray(row, dtype=int).tobytes()
#         prev = coeff_by_row.get(k)
#         if prev is None:
#             coeff_by_row[k] = complex(coeff)
#             row_by_key[k] = np.asarray(row, dtype=int).copy()
#             return
#         new = prev + complex(coeff)
#         if abs(new) <= 1e-14:
#             # Exact (or near) cancellation; drop the term to keep the Hamiltonian compact.
#             coeff_by_row.pop(k, None)
#             row_by_key.pop(k, None)
#         else:
#             coeff_by_row[k] = new

#     max_attempts = 1000
#     attempts = 0
#     while len(coeff_by_row) < n_paulis:
#         attempts += 1
#         if attempts > max_attempts:
#             raise RuntimeError(
#                 f"Failed to generate a non-trivial symmetric Hamiltonian after {max_attempts} attempts "
#                 f"(have {len(coeff_by_row)} unique terms, want {n_paulis})."
#             )

#         seed = PauliSum.from_random(1, dims, rand_weights=False, rand_phases=False)
#         seed.weights = np.array([_new_weight()], dtype=complex)
#         seed = seed.to_standard_form()

#         for term in _orbit(seed):
#             _add_coeff(term.tableau[0], term.weights[0])

#     tableaus = np.vstack(list(row_by_key.values()))
#     weights = np.array(list(coeff_by_row.values()), dtype=np.complex128)
#     phases = np.zeros(weights.shape[0], dtype=int)

#     P_sym = PauliSum.from_tableau(tableaus, dimensions=dims, weights=weights, phases=phases)
#     P_sym.standardise()
#     P_sym.set_weights(np.around(P_sym.weights, decimals=10))
#     P_sym.remove_zero_weight_paulis()

#     if scrambled is True:
#         g = Gate.from_random(n_qudits, dims[0])
#         P_sym = g.act(P_sym, qudit_indices)

#     return P_sym


# --- GF(2) rank tracker for rows of length 2n (int dtype 0/1) ---
class GF2RankTracker:
    def __init__(self, m: int):
        self.m = int(m)
        self._piv = {}  # pivot_index -> row_bits

    @property
    def rank(self) -> int:
        return len(self._piv)

    def _reduce(self, x: int) -> int:
        for p in sorted(self._piv.keys(), reverse=True):
            if (x >> p) & 1:
                x ^= self._piv[p]
        return x

    def would_increase_rank(self, x: int) -> bool:
        return self._reduce(x) != 0

    def add(self, x: int) -> bool:
        x = self._reduce(x)
        if x == 0:
            return False
        p = x.bit_length() - 1
        # clean existing rows w.r.t new pivot
        for q in list(self._piv.keys()):
            if (self._piv[q] >> p) & 1:
                self._piv[q] ^= x
        self._piv[p] = x
        return True


def tableau_row_to_bits(v: np.ndarray) -> int:
    b = (np.asarray(v, dtype=np.uint8) & 1)
    packed = np.packbits(b, bitorder="little")
    return int.from_bytes(packed.tobytes(), byteorder="little", signed=False)


def gf2_rank_of_tableau(tableau: np.ndarray) -> int:
    if tableau.size == 0:
        return 0
    m = tableau.shape[1]
    r = GF2RankTracker(m)
    for row in tableau:
        r.add(tableau_row_to_bits(row))
    return r.rank


def tableau_basis_seeds(n_qudits: int) -> np.ndarray:
    """
    Return the 2n basis vectors in GF(2)^{2n} as a (2n, 2n) int array.
    First n are X_i, last n are Z_i in your [x|z] convention.
    """
    m = 2 * n_qudits
    B = np.zeros((m, m), dtype=np.int64)
    for i in range(m):
        B[i, i] = 1
    return B


def random_gate_symmetric_hamiltonian(G: Gate,
                                      dimension: int,
                                      qudit_indices: tuple[int, ...],
                                      n_qudits: int,
                                      n_paulis: int | None = None,
                                      weight_mode: str = 'uniform',
                                      scrambled: bool = False,
                                      # generation controls
                                      extra_orbit_budget: int = 10_000,   # how many extra random orbit seeds to add after basis orbits
                                      avoid_rounding: bool = True,        # recommended True for rank robustness
                                      ) -> PauliSum:
    """
    Symmetric Hamiltonian generator that is full-rank (rank 2n) by construction
    (before any combining/cancellation), using the 2n canonical basis seeds.

    Notes:
      - Returns >= n_paulis terms (because we add whole orbits).
      - Uses the library's projective action as-is: we keep term.weights/phases from G.act.
    """
    if n_qudits is None:
        n_qudits = len(qudit_indices)
    if n_paulis is None:
        n_paulis = 2 * n_qudits

    rng = np.random.default_rng()
    all_indices = tuple(range(n_qudits))

    def new_seed_weight() -> complex:
        if weight_mode == "random":
            # Use generic complex weights to reduce accidental cancellations.
            a = float(rng.random()) + 0.1
            b = float(rng.random()) + 0.1
            return complex(a, b)
        return complex(1.0)

    def orbit(seed_term: "PauliSum") -> list["PauliSum"]:
        """
        IMPORTANT: match your old working closure condition. In your representation
        the orbit may only close when both tableau and weight repeat.
        """
        seen = set()
        out = []
        term = seed_term.to_standard_form()
        while True:
            key = (tuple(term.tableau[0].tolist()), complex(np.around(term.weights[0], decimals=12)))
            if key in seen:
                break
            seen.add(key)
            out.append(term)
            term = G.act(pauli=term, qudits=qudit_indices).to_standard_form()
        return out

    # ---- Phase 1: add orbits for the 2n basis seeds ----
    basis_rows = tableau_basis_seeds(n_qudits)  # shape (2n, 2n)

    orbit_blocks = []  # each block: (tableau_rows, phases, weights)
    total_terms = 0

    for row in basis_rows:
        # Make a 1-term PauliSum for this basis tableau row
        w0 = new_seed_weight()
        seed = PauliSum.from_tableau(
            row[None, :],
            dimensions=[dimension] * n_qudits,
            weights=np.array([w0], dtype=complex),
            phases=np.array([0], dtype=int),
        ).to_standard_form()

        blk = orbit(seed)

        tab_rows = [t.tableau[0] for t in blk]
        phs = [int(t.phases[0]) for t in blk]
        wts = [complex(t.weights[0]) for t in blk]  # keep projective factors as represented

        orbit_blocks.append((np.vstack(tab_rows), np.array(phs, dtype=int), np.array(wts, dtype=complex)))
        total_terms += len(blk)

    # ---- Phase 2: add extra orbit blocks until we have >= n_paulis ----
    # (These can be dependent; they just fill out the Hamiltonian.)
    seeds_used = 0
    while total_terms < n_paulis and seeds_used < extra_orbit_budget:
        seeds_used += 1
        seed = PauliSum.from_random(1, [dimension] * n_qudits, rand_weights=False, rand_phases=False)
        seed.weights = np.array([new_seed_weight()], dtype=complex)
        seed = seed.to_standard_form()

        blk = orbit(seed)

        tab_rows = [t.tableau[0] for t in blk]
        phs = [int(t.phases[0]) for t in blk]
        wts = [complex(t.weights[0]) for t in blk]

        orbit_blocks.append((np.vstack(tab_rows), np.array(phs, dtype=int), np.array(wts, dtype=complex)))
        total_terms += len(blk)

    # ---- Build PauliSum from blocks ----
    all_tableau = np.vstack([b[0] for b in orbit_blocks])
    all_phases = np.concatenate([b[1] for b in orbit_blocks])
    all_weights = np.concatenate([b[2] for b in orbit_blocks])

    P_sym = PauliSum.from_tableau(all_tableau, dimensions=[dimension]
                                  * n_qudits, weights=all_weights, phases=all_phases)

    # Make combining safe: absorb phases into weights before combine if combine ignores phases.
    P_sym.phase_to_weight()

    # Combine/standardise. Avoid rounding here (rounding can induce exact cancellation).
    P_sym.combine_equivalent_paulis()
    P_sym.standardise()

    if not avoid_rounding:
        P_sym.set_weights(np.around(P_sym.weights, decimals=10))

    P_sym.remove_zero_weight_paulis()

    # If your downstream expects phases explicitly, you can move them back out:
    P_sym.weight_to_phase()
    P_sym.standardise()

    if scrambled:
        g = Gate.from_random(n_qudits, dimension)
        P_sym = g.act(P_sym, all_indices)

    return P_sym
