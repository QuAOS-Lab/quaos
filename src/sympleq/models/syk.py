import numpy as np
from itertools import combinations
from typing import Tuple, Optional
from sympleq.core.paulis import PauliSum
from .utils import _zeros_xz, _add_pauli_term, _finalize_terms, _pauli_mul


def _jw_majorana_operator(n_qubits: int, a: int) -> tuple[np.ndarray, np.ndarray, int]:
    r"""
    Jordan-Wigner Majorana operators on n_qubits fermionic modes (= n_qubits qubits).

    We define 2n Majoranas (a = 0,1,...,2n-1):
        gamma_{2j}   = (Π_{k<j} Z_k) X_j
        gamma_{2j+1} = (Π_{k<j} Z_k) Y_j

    In the convention Y = i X Z, gamma_{2j+1} has phase=1, x_j=1, z_{0..j}=1.
    """
    n = int(n_qubits)
    if a < 0 or a >= 2 * n:
        raise ValueError(f"Majorana index a={a} out of range for n_qubits={n} (need 0..{2 * n - 1}).")

    j, parity = divmod(a, 2)
    x, z = _zeros_xz(n)

    if j > 0:
        z[:j] = 1  # JW string

    x[j] = 1
    if parity == 0:
        # X_j
        phase = 0
    else:
        # Y_j = i X_j Z_j
        z[j] = 1
        phase = 1

    return x, z, phase


def _syk4_majorana_hamiltonian(
    n_qubits: int,
    J_scale: float = 1.0,
    seed: Optional[int] = None,
    couplings: Optional[np.ndarray] = None,
    standard_normalization: bool = True,
    include_constant: bool = False,
    tol: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""
    Build a qubit Pauli-tableau for the SYK_4 Majorana Hamiltonian on n_qubits fermionic modes:

        H = Σ_{a<b<c<d} J_{abcd} γ_a γ_b γ_c γ_d

    where there are N_M = 2*n_qubits Majoranas, mapped to qubits by Jordan-Wigner.

    Parameters
    ----------
    n_qubits : int
        Number of fermionic modes (also the number of qubits after JW mapping).
    J_scale : float
        Overall scale used when generating random couplings if couplings is None.
    seed : int | None
        RNG seed for generated couplings.
    couplings : np.ndarray | None
        Optional explicit antisymmetric coupling tensor values on ordered tuples.
        Expected shape (N_M, N_M, N_M, N_M). Only entries a<b<c<d are read.
    standard_normalization : bool
        If True and couplings is None, use
            std(J_abcd) = sqrt(3!)*J_scale / N_M^(3/2),
        a common SYK_4 normalization (up to convention factors).
        If False, use std(J_abcd) = J_scale.
    include_constant : bool
        Present for API symmetry (SYK_4 has no identity term here), ignored.
    tol : float
        Prune coefficients with |c| <= tol.

    Notes
    -----
    - With real J_abcd, the Hamiltonian is Hermitian.
    - In this tableau convention the resulting phases should be even (0 or 2).
    """
    n = int(n_qubits)
    if n <= 1:
        raise ValueError("n_qubits must be >= 2 for a nontrivial SYK_4 model.")
    N_M = 2 * n

    rng = np.random.default_rng(seed)

    if couplings is not None:
        J4 = np.asarray(couplings, dtype=float)
        if J4.shape != (N_M, N_M, N_M, N_M):
            raise ValueError(
                f"couplings must have shape ({N_M},{N_M},{N_M},{N_M}), got {J4.shape}."
            )

        def get_coupling(a, b, c, d):
            if not (a < b < c < d):
                raise ValueError(f"Couplings should only be read on ordered tuples a<b<c<d; got {(a, b, c, d)}.")
            return float(J4[a, b, c, d])
    else:
        if standard_normalization:
            std = (np.sqrt(6.0) * float(J_scale)) / (N_M ** 1.5)
        else:
            std = float(J_scale)

        # generate on the fly (saves memory vs full tensor)
        def get_coupling(a, b, c, d):
            return float(rng.normal(loc=0.0, scale=std))

    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], float] = {}

    # Precompute JW Majorana tableau reps
    gamma = [_jw_majorana_operator(n, a) for a in range(N_M)]

    for a, b, c, d in combinations(range(N_M), 4):
        J_abcd = get_coupling(a, b, c, d)
        if J_abcd == 0.0:
            continue

        x, z, p = gamma[a]
        x, z, p = _pauli_mul(x, z, p, *gamma[b])
        x, z, p = _pauli_mul(x, z, p, *gamma[c])
        x, z, p = _pauli_mul(x, z, p, *gamma[d])

        # For real SYK_4 couplings, phase should be 0 or 2 (Hermitian Pauli strings)
        # We do not hard-fail; just leave as-is in case of convention tweaks.
        _add_pauli_term(terms, x, z, phase=p, coeff=J_abcd)

    # include_constant exists only for interface symmetry (no-op here)
    _ = include_constant

    return _finalize_terms(n, terms, tol=tol)


def syk4_majorana_hamiltonian(
    n_qubits: int,
    J_scale: float = 1.0,
    seed: Optional[int] = None,
    couplings: Optional[np.ndarray] = None,
    standard_normalization: bool = True,
    tol: float = 0.0,
) -> 'PauliSum':
    tableau, coeffs, phases = _syk4_majorana_hamiltonian(
        n_qubits=n_qubits,
        J_scale=J_scale,
        seed=seed,
        couplings=couplings,
        standard_normalization=standard_normalization,
        tol=tol,
    )
    return PauliSum.from_tableau(tableau, weights=coeffs, phases=phases)
