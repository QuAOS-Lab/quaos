import numpy as np


def _zeros_xz(n: int) -> tuple[np.ndarray, np.ndarray]:
    return np.zeros(n, dtype=np.uint8), np.zeros(n, dtype=np.uint8)


def _add_pauli_term(
    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], float],
    x: np.ndarray,
    z: np.ndarray,
    phase: int,
    coeff: float,
) -> None:
    phase &= 3
    key = (phase, tuple(int(v) for v in x), tuple(int(v) for v in z))
    terms[key] = terms.get(key, 0.0) + float(coeff)


def _finalize_terms(
    n: int,
    terms: dict[tuple[int, tuple[int, ...], tuple[int, ...]], float],
    tol: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    items = []
    for key, c in terms.items():
        if tol > 0.0 and abs(c) <= tol:
            continue
        if c != 0.0:
            items.append((key, c))

    # deterministic order
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


def _pauli_mul(
    x1: np.ndarray, z1: np.ndarray, p1: int,
    x2: np.ndarray, z2: np.ndarray, p2: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Multiply Paulis in the convention:
        P = (1j)**phase * Π_k X_k^{x_k} Z_k^{z_k}
    with x,z in {0,1}.

    Product rule:
      X^x Z^z  X^{x'} Z^{z'} = (-1)^{z·x'} X^{x+x'} Z^{z+z'}
    so phase increment gets +2*(z·x').
    """
    cross = int((np.dot(z1.astype(np.uint8), x2.astype(np.uint8)) % 2))
    phase = (int(p1) + int(p2) + 2 * cross) & 3
    x = (x1 ^ x2).astype(np.uint8)
    z = (z1 ^ z2).astype(np.uint8)
    return x, z, phase
