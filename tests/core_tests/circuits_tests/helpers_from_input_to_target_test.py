import numpy as np
import pytest

from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Circuit
from sympleq.core.circuits.helpers_from_input_to_target import (gf_rank,
                                                                independent_solution,
                                                                complete_basis,
                                                                map_tableau_to_target_tableau,
                                                                solve_mod_2p)
from sympleq.core.circuits.utils import symplectic_form, symplectic_product_matrix
from sympleq.models import random_hamiltonian


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("n_qudits", [10])
@pytest.mark.parametrize("num_pauli", [5, 10])
def test_independent_solution(dim: int, n_qudits: int, num_pauli: int):
    p = dim
    dimensions = [p] * n_qudits

    pl_sum = random_hamiltonian.random_pauli_hamiltonian(num_pauli, dimensions)
    target_pl_sum = Circuit.from_random(
        n_gates=10 * n_qudits**2,
        dimensions=dimensions,
    ).act(pl_sum)

    U = np.asarray(pl_sum.tableau, dtype=int) % p
    V = np.asarray(target_pl_sum.tableau, dtype=int) % p

    d = U.shape[1]
    Omega = symplectic_form(d // 2, p)

    found_independent_solution = False

    for e in np.eye(d, dtype=int):
        if gf_rank(np.vstack([U, e]), p) != gf_rank(U, p) + 1:
            continue

        A = (V @ Omega) % p
        b = (U @ Omega @ e) % p

        try:
            y = independent_solution(A, b, V, p)
        except ValueError:
            continue

        assert np.array_equal((A @ y) % p, b % p)
        assert gf_rank(np.vstack([V, y]), p) == gf_rank(V, p) + 1
        found_independent_solution = True
        break
    assert found_independent_solution


def test_independent_solution_no_solution():
    p = 2

    A = np.eye(2, dtype=int)
    b = np.array([0, 0], dtype=int)

    # V already spans GF(2)^2, so no y can be outside span(V).
    V = np.eye(2, dtype=int)

    with pytest.raises(ValueError, match="Target tableau already contain a complete basis"):
        independent_solution(A, b, V, p)


def test_independent_solution_raises_no_solution():
    p = 2

    A = np.array([
        [0, 0],
    ], dtype=int)

    b = np.array([1], dtype=int)

    V = np.array([
        [1, 0],
    ], dtype=int)

    with pytest.raises(ValueError):
        independent_solution(A, b, V, p)


def test_complete_basis_full_rank():
    p = 2
    n_qudits = 2
    dim = [p] * n_qudits

    generator_tableau = np.eye(2 * n_qudits, dtype=int)
    generator_ps = PauliSum.from_tableau(generator_tableau, dim)
    target_tableau = Circuit.from_random(
        n_gates=10 * n_qudits**2,
        dimensions=dim,
    ).act(generator_ps).tableau

    with pytest.raises(ValueError, match="Input and target tableaus already contain a complete basis"):
        complete_basis(generator_tableau, target_tableau, p)


@pytest.mark.parametrize("dim", [2, 3, 5])
@pytest.mark.parametrize("n_qudits", [3, 4])
def test_complete_basis_success(
    dim: int,
    n_qudits: int,
):
    num_pauli = n_qudits
    dimensions = [dim] * n_qudits

    pl_sum = random_hamiltonian.random_pauli_hamiltonian(num_pauli, dimensions)
    target_pl_sum = Circuit.from_random(
        n_gates=10 * n_qudits**2,
        dimensions=dimensions,
    ).act(pl_sum)

    input_tab = np.asarray(pl_sum.tableau, dtype=int) % dim
    target_tab = np.asarray(target_pl_sum.tableau, dtype=int) % dim

    U_full, V_full = complete_basis(input_tab, target_tab, dim)

    d = input_tab.shape[1]

    assert U_full.shape == (d, d)
    assert V_full.shape == (d, d)
    assert gf_rank(U_full, dim) == d
    assert gf_rank(V_full, dim) == d

    assert np.array_equal(
        symplectic_product_matrix(U_full, dim) % dim,
        symplectic_product_matrix(V_full, dim) % dim,
    )


@pytest.mark.parametrize("dim", [2, 3, 5])
@pytest.mark.parametrize("n_qudits", [3, 4])
@pytest.mark.parametrize("num_pauli", [2, 3, 30])
def test_map_tableau_to_target_tableau_success(
    dim: int,
    n_qudits: int,
    num_pauli: int
):
    dimensions = [dim] * n_qudits

    pl_sum = random_hamiltonian.random_pauli_hamiltonian(num_pauli, dimensions)
    target_pl_sum = Circuit.from_random(
        n_gates=10 * n_qudits**2,
        dimensions=dimensions,
    ).act(pl_sum)

    input_tab = np.asarray(pl_sum.tableau, dtype=int) % dim
    target_tab = np.asarray(target_pl_sum.tableau, dtype=int) % dim

    F = map_tableau_to_target_tableau(input_tab, target_tab, dim)

    assert F.shape == (input_tab.shape[1], input_tab.shape[1])
    assert np.array_equal((input_tab @ F) % dim, target_tab)


def test_map_tableau_to_target_tableau_failure():
    # Sends two independent vectors to the same target vector,
    # which is impossible for any invertible F.
    dim = 2
    input_tab = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
    ], dtype=int)

    target_tab = np.array([
        [1, 0, 0, 0],
        [1, 0, 0, 0],
    ], dtype=int)

    with pytest.raises(ValueError):
        map_tableau_to_target_tableau(input_tab, target_tab, dim)


@pytest.mark.parametrize("p", [3, 5, 7])
def test_solve_mod_2p_success(p: int):
    A = np.array([
        [1, 0, 2],
        [0, 1, 1],
        [1, 1, 0],
    ], dtype=int)

    h_expected = np.array([1, 2, 3], dtype=int) % (2 * p)
    delta = (A @ h_expected) % (2 * p)

    h_found = solve_mod_2p(A, delta, p)

    assert np.array_equal((A @ h_found) % (2 * p), delta)


@pytest.mark.parametrize("p", [3, 5, 7])
def test_solve_mod_2p_failure(p: int):
    A = np.array([
        [0, 0, 0],
    ], dtype=int)

    delta = np.array([1], dtype=int)

    with pytest.raises(ValueError):
        solve_mod_2p(A, delta, p)
