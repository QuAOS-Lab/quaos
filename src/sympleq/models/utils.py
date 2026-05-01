import numpy as np


def _zeros_xz(n_qubits: int) -> tuple[np.ndarray, np.ndarray]:
    return np.zeros(n_qubits, dtype=int), np.zeros(n_qubits, dtype=int)


def _add_pauli_term(
    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], complex],
    x: np.ndarray,
    z: np.ndarray,
    phase: int,
    coeff: complex,
) -> None:
    key = (int(phase) % 4, tuple(np.asarray(x, dtype=int)), tuple(np.asarray(z, dtype=int)))
    terms[key] = terms.get(key, 0.0) + coeff


def _finalize_terms(
    n_qubits: int,
    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], complex],
    tol: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = []
    coeffs = []
    phases = []

    for (phase, x, z), coeff in terms.items():
        if abs(coeff) <= tol:
            continue
        rows.append(np.array([*x, *z], dtype=int))
        coeffs.append(coeff)
        phases.append(phase)

    if not rows:
        return (
            np.zeros((0, 2 * n_qubits), dtype=int),
            np.zeros(0, dtype=float),
            np.zeros(0, dtype=int),
        )

    return np.vstack(rows), np.asarray(coeffs), np.asarray(phases, dtype=int)
