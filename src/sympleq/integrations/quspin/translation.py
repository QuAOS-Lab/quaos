import numpy as np
from collections import defaultdict
from typing import Optional, Sequence, Tuple, Dict, Any
from quspin.basis import spin_basis_1d
from quspin.operators import hamiltonian


def quspin_hamiltonian_from_tableau(
    tableau: Any,
    *,
    weights: Optional[Sequence[complex]] = None,
    phases: Optional[Sequence[int | float]] = None,
    basis: Optional[spin_basis_1d] = None,
    dtype=np.complex128,
    phase_mode: str = "i_quarters",
    return_constant: bool = True,
) -> hamiltonian | Tuple[hamiltonian, complex]:
    """
    Convert a Pauli tableau representation of a Hamiltonian into a QuSpin hamiltonian.

    Supported inputs
    ----------------
    1) Concatenated tableau array of shape (m, 2n) with columns [X | Z], entries in {0,1}.
       - Provide `weights` (len m) and optional `phases` (len m).

    2) Dict with keys {"X": (m,n) array, "Z": (m,n) array, optional "weights", optional "phases"}.

    3) Minimal PauliSum-like: any object with attributes `x_exp`, `z_exp`, and optional
       `weights` / `phases` per term (each of shape (m, n) for exponents).

    Coefficients
    ------------
    The complex coefficient for term k is composed as:
        coeff_k = (weights[k] or 1) * phase_factor(phases[k])

    `phase_mode` selects how to interpret `phases`:
        - "i_quarters": integer phases in {0,1,2,3} meaning i**phase
        - "radians"   : float angle θ meaning exp(1j*θ)
        - "sign"      : integer phases in {0,1} meaning (-1)**phase
        - "none"      : ignore `phases` (treat as 0)

    Notes on Y and phases
    ---------------------
    We map (x_j, z_j) to site-operators as:
        (0,0) -> I (site omitted)
        (1,0) -> 'x'
        (0,1) -> 'z'
        (1,1) -> 'y'
    QuSpin's 'y' already includes the ±i factors inside σ_y. Therefore, **do not**
    add the canonical i^{∑_j x_j z_j} multiplier on top. If your tableau phase
    already encodes that canonical factor, keep `phases=None` (or set `phase_mode="none"`)
    and pass the full complex `weights` instead.

    Identity terms
    --------------
    Rows with all zeros (I⊗…⊗I) contribute a constant energy shift. We accumulate these
    in `const_shift` and (by default) return it alongside the QuSpin `hamiltonian`.

    Parameters
    ----------
    tableau : array-like, dict, or object
        See “Supported inputs”.
    weights : sequence of complex, optional
        Per-term coefficients. If omitted, defaults to 1.0 for all terms unless provided
        in `tableau` when it is a dict/object.
    phases : sequence of int|float, optional
        Per-term phases, interpreted by `phase_mode`.
    basis : quspin.basis.spin_basis_1d, optional
        If None, a default `spin_basis_1d(n)` (no symmetries) is constructed.
    dtype : numpy dtype
        dtype for the resulting Hamiltonian, default complex128.
    phase_mode : {"i_quarters","radians","sign","none"}
        How to interpret `phases` (see above).
    return_constant : bool
        If True, return a tuple (H, const_shift). Else return only H.

    Returns
    -------
    H : quspin.operators.hamiltonian
        The QuSpin Hamiltonian constructed from the tableau.
    const_shift : complex (optional)
        The constant energy shift from identity terms (if any).

    Example
    -------
    >>> # m=2 terms over n=3 qubits:
    >>> # term 0: X_0 Z_2  (weight 1.2)
    >>> # term 1: Y_1      (weight -0.5 i)
    >>> T = np.array([
    ...     [1,0,0, 0,0,1],   # x|z
    ...     [0,1,0, 0,1,0],
    ... ], dtype=int)
    >>> w = [1.2, -0.5j]
    >>> H, c = quspin_hamiltonian_from_tableau(T, weights=w, return_constant=True)
    >>> # Use H.eigsh(), H.dot(psi), etc.
    """
    # --- Normalize input to X, Z, weights, phases ---
    def _extract_XZ_w_phi(tab) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
        local_weights = None
        local_phases = None
        if isinstance(tab, dict):
            X = np.asarray(tab["X"], dtype=int)
            Z = np.asarray(tab["Z"], dtype=int)
            if "weights" in tab and weights is None:
                local_weights = np.asarray(tab["weights"], dtype=complex)
            if "phases" in tab and phases is None:
                local_phases = np.asarray(tab["phases"])
        elif hasattr(tab, "x_exp") and hasattr(tab, "z_exp"):
            X = np.asarray(tab.x_exp, dtype=int)
            Z = np.asarray(tab.z_exp, dtype=int)
            if hasattr(tab, "weights") and weights is None:
                local_weights = np.asarray(tab.weights, dtype=complex)
            if hasattr(tab, "phases") and phases is None:
                local_phases = np.asarray(tab.phases)
        else:
            T = np.asarray(tab)
            assert T.ndim == 2 and T.shape[1] % 2 == 0, \
                "Concatenated tableau must be shape (m, 2n)."
            n = T.shape[1] // 2
            X = T[:, :n].astype(int, copy=False)
            Z = T[:, n:].astype(int, copy=False)
        return X, Z, local_weights, local_phases

    X, Z, w_from_tab, phi_from_tab = _extract_XZ_w_phi(tableau)
    m, n = X.shape
    assert Z.shape == (m, n), "X and Z must have same shape (m, n)."

    # Resolve weights / phases order of precedence:
    if weights is not None:
        w = np.asarray(weights, dtype=complex)
    elif w_from_tab is not None:
        w = w_from_tab.astype(complex, copy=False)
    else:
        w = np.ones(m, dtype=complex)

    if phases is not None:
        phi = np.asarray(phases)
    else:
        phi = phi_from_tab

    # Phase factor interpreter
    def _phase_factor(p):
        if p is None:
            return 1.0
        if phase_mode == "none":
            return 1.0
        if phase_mode == "i_quarters":
            # integers in {0,1,2,3}: i**p
            return (1j) ** int(p)
        if phase_mode == "radians":
            # float angle theta: exp(i theta)
            return np.exp(1j * float(p))
        if phase_mode == "sign":
            # integers in {0,1}: (-1)**p
            return (-1) ** int(p)
        raise ValueError(f"Unknown phase_mode '{phase_mode}'.")

    # Build (op_str, sites) -> coefficient accumulator
    accumulator: Dict[Tuple[str, Tuple[int, ...]], complex] = defaultdict(complex)
    const_shift = 0.0 + 0.0j

    for k in range(m):
        x = X[k]
        z = Z[k]

        # Identify active sites and local operators
        sites = []
        ops = []
        for j in range(n):
            xj, zj = int(x[j]), int(z[j])
            if (xj | zj) == 0:
                continue
            if xj == 1 and zj == 0:
                ops.append("x")
                sites.append(j)
            elif xj == 0 and zj == 1:
                ops.append("z")
                sites.append(j)
            else:  # xj == 1 and zj == 1
                ops.append("y")
                sites.append(j)

        if len(sites) == 0:
            # Identity term: contributes a constant energy shift
            coeff = w[k] * (_phase_factor(phi[k]) if phi is not None else 1.0)
            const_shift += coeff
            continue

        op_str = "".join(ops)
        coeff = w[k] * (_phase_factor(phi[k]) if phi is not None else 1.0)
        accumulator[(op_str, tuple(sites))] += coeff

    # Prepare QuSpin static list
    static = []
    for (op_str, sites), coeff in accumulator.items():
        # QuSpin wants one entry per op_str with a site-coupling list:
        # [op_str, [[coeff, i, j, k, ...], ...]]
        static.append([op_str, [[complex(coeff), *sites]]])

    # Build basis if needed
    if basis is None:
        basis = spin_basis_1d(n)  # no symmetries by default

    H = hamiltonian(static, [], basis=basis, dtype=dtype)

    if return_constant:
        return H, complex(const_shift)
    return H
