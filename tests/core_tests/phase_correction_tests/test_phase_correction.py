import numpy as np

from sympleq.core.circuits.phase_correction import (
    solve_phase_vector_h_from_residual as solve_from_circuits,
)
from sympleq.core.phase_correction import solve_phase_vector_h_from_residual
from sympleq.core.symmetries.phase_correction import (
    solve_phase_vector_h_from_residual as solve_from_symmetries,
)


def test_solve_phase_vector_qubit_residual():
    tableau = np.eye(2, dtype=int)
    delta = np.array([2, 0], dtype=int)

    h = solve_phase_vector_h_from_residual(tableau, delta, [2])

    assert np.array_equal(h % 4, np.array([2, 0]))


def test_solve_phase_vector_odd_prime_residual():
    tableau = np.eye(2, dtype=int)
    delta = np.array([1, 2], dtype=int)

    h = solve_phase_vector_h_from_residual(tableau, delta, [3])

    assert np.array_equal(h % 3, np.array([1, 2]))


def test_solve_phase_vector_mixed_dimensions_returns_none():
    tableau = np.zeros((1, 4), dtype=int)
    delta = np.zeros(1, dtype=int)

    h = solve_phase_vector_h_from_residual(tableau, delta, [2, 3])

    assert h is None


def test_phase_correction_compatibility_imports():
    assert solve_from_circuits is solve_phase_vector_h_from_residual
    assert solve_from_symmetries is solve_phase_vector_h_from_residual
