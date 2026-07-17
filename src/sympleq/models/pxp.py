from __future__ import annotations

import numpy as np

from sympleq.core.paulis import PauliSum, PauliString


def pxp_model(
    n_qubits: int,
    coupling: float = 1.0,
    periodic: bool = False,
    normalized_projectors: bool = False,
) -> PauliSum:
    """
    Build the qubit PXP Hamiltonian.

    We use
        H = J * sum_i P_{i-1} X_i P_{i+1},
    with P = I - Z by default, or P = (I - Z)/2 if
    `normalized_projectors=True`.

    Parameters
    ----------
    n_qubits : int
        Number of qubits.
    coupling : float, optional
        Overall coupling J multiplying the Hamiltonian.
    periodic : bool, optional
        If True, use periodic boundary conditions; otherwise open boundaries.
    normalized_projectors : bool, optional
        If True, use projectors P=(I-Z)/2. If False, use P=(I-Z).

    Returns
    -------
    PauliSum
        PXP Hamiltonian as a PauliSum object.
    """
    n = int(n_qubits)
    if n < 3:
        raise ValueError("PXP model requires n_qubits >= 3.")

    dims = [2] * n
    paulis: list[PauliString] = []
    weights: list[float] = []

    prefactor = float(coupling) * (0.25 if normalized_projectors else 1.0)

    if periodic:
        centers = range(n)
    else:
        centers = range(1, n - 1)

    for i in centers:
        left = (i - 1) % n
        right = (i + 1) % n

        # +X_i
        x = np.zeros(n, dtype=int)
        z = np.zeros(n, dtype=int)
        x[i] = 1
        paulis.append(PauliString.from_exponents(x, z, dims))
        weights.append(prefactor)

        # -Z_{left} X_i
        x = np.zeros(n, dtype=int)
        z = np.zeros(n, dtype=int)
        x[i] = 1
        z[left] = 1
        paulis.append(PauliString.from_exponents(x, z, dims))
        weights.append(-prefactor)

        # -X_i Z_{right}
        x = np.zeros(n, dtype=int)
        z = np.zeros(n, dtype=int)
        x[i] = 1
        z[right] = 1
        paulis.append(PauliString.from_exponents(x, z, dims))
        weights.append(-prefactor)

        # +Z_{left} X_i Z_{right}
        x = np.zeros(n, dtype=int)
        z = np.zeros(n, dtype=int)
        x[i] = 1
        z[left] = 1
        z[right] = 1
        paulis.append(PauliString.from_exponents(x, z, dims))
        weights.append(prefactor)

    H = PauliSum.from_pauli_strings(paulis, weights=weights)
    H.combine_equivalent_paulis()
    H.remove_zero_weight_paulis()
    return H
