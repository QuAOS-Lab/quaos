from sympleq.core.circuits.find_symplectic import (
    map_single_pauli_string_to_target,
    find_symplectic_solution,
    find_symplectic_solution_extended,
    solve_gf2,
    map_pauli_sum_to_target_tableau,
    symplectic_from_pauli_permutation,
)
import numpy as np
import pytest
from sympleq.core.circuits.utils import transvection, transvection_matrix, symplectic_product_arrays
from sympleq.core.finite_field_solvers import get_linear_dependencies
from sympleq.models import random_hamiltonian
from sympleq.core.circuits import Circuit


class TestSymplecticSolver:

    def verify_solution_extended(self, u: np.ndarray, v: np.ndarray, w: np.ndarray,
                                 t_vectors: list | None = None) -> bool:
        """
        Verify that w satisfies all the extended symplectic conditions.

        Args:
            u, v, w: Binary vectors
            t_vectors: List of additional constraint vectors

        Returns:
            True if all conditions are satisfied
        """
        if t_vectors is None:
            t_vectors = []

        # Check primary conditions
        if not (symplectic_product_arrays(u, w) == 1 and symplectic_product_arrays(v, w) == 1):
            return False

        # Check additional conditions
        for t in t_vectors:
            if symplectic_product_arrays(t, w) != symplectic_product_arrays(t, v):
                return False

        return True

    def verify_solution(self, u: np.ndarray, v: np.ndarray, w: np.ndarray) -> bool:
        """
        Verify that w satisfies the symplectic conditions.

        Args:
            u, v, w: Binary vectors

        Returns:
            True if <u,w> = <v,w> = 1
        """
        return (symplectic_product_arrays(u, w) == 1 and symplectic_product_arrays(v, w) == 1)

    @staticmethod
    def assert_transvection_property(x, h):
        """
        For x, h in GF(2)^2n

        Z_h(x) = x if <x, h> = 0, and h if <x, h> = 1
        """
        if symplectic_product_arrays(x, h) == 0:
            assert np.all(transvection(h, x) == x)
            assert np.all((x.T @ transvection_matrix(h)) % 2 == x), f"\n{x}\n{h}\n{(x @ transvection_matrix(h)) % 2}"
        elif symplectic_product_arrays(x, h) == 1:
            assert np.all(transvection(h, x) == (h + x) % 2)
            assert np.all((x.T @ transvection_matrix(h)) % 2 == (h + x) % 2), (f"\n{(x + h) % 2}"
                                                                               f"\n{(x @ transvection_matrix(h)) % 2}")
        else:
            raise Exception('This test only works in GF(2)')

    def test_transvection_property(self):
        n = 10
        h = np.random.randint(2, size=2 * n)
        for _ in range(100):
            x = np.random.randint(2, size=2 * n)
            self.assert_transvection_property(x, h)

    def test_transvection_and_transvection_matrix(self):
        n = 10
        h = np.random.randint(2, size=2 * n)
        for _ in range(100):
            x = np.random.randint(2, size=2 * n)

            assert np.all(transvection(h, x) == (x @ transvection_matrix(h)) % 2)

    def test_map_single_pauli_string_to_target(self):
        n = 14
        p = 2
        for _ in range(1000):
            input_ps = np.random.randint(p, size=2 * n)
            target_ps = np.random.randint(p, size=2 * n)
            if np.array_equal(input_ps, np.zeros(2 * n)):
                input_ps[np.random.randint(2 * n)] = 1  # Ensure non-zero input
            if np.array_equal(target_ps, np.zeros(2 * n)):
                target_ps[np.random.randint(2 * n)] = 1  # Ensure non-zero target

            F_map = map_single_pauli_string_to_target(input_ps, target_ps)
            assert np.all((input_ps @ F_map) % p == target_ps), (f"\n{F_map}\n{input_ps}"
                                                                 f"\n{(input_ps @ F_map) % p}\n{target_ps}")

    def test_symplectic_from_pauli_permutation(self):
        tableau = np.eye(4, dtype=int)
        permutation = np.array([1, 0, 3, 2])

        F_map = symplectic_from_pauli_permutation(tableau, permutation, p=2)

        assert np.all((tableau @ F_map) % 2 == tableau[permutation])

    def test_symplectic_from_pauli_permutation_uses_inverse_by_default(self):
        tableau = np.array([[1, 0, 0, 0], [0, 0, 1, 0]], dtype=int)
        permutation = np.array([1, 0])

        with pytest.raises(ValueError, match="complete independent Pauli basis"):
            symplectic_from_pauli_permutation(tableau, permutation, p=2)

        F_map = symplectic_from_pauli_permutation(tableau, permutation, p=2, method="transvection")
        assert np.all((tableau @ F_map) % 2 == tableau[permutation])

    def test_find_w(self):
        """Test function that properly handles cases where no solution exists."""
        n = 10
        p = 2

        for i in range(1000):
            u = np.random.randint(p, size=2 * n)
            v = np.random.randint(p, size=2 * n)
            # ensure neither u and v are zero, as that would make the system consistent with w=0
            if np.array_equal(u, np.zeros(2 * n)):
                u[np.random.randint(2 * n)] = 1  # Ensure non-zero input
            if np.array_equal(v, np.zeros(2 * n)):
                v[np.random.randint(2 * n)] = 1  # Ensure non-zero target

            w = find_symplectic_solution(u, v)

            if w is not None:
                # If a solution is found, it must be valid
                assert self.verify_solution(u, v, w), (f"Solution verification failed for u={u},'\
                                                ' v={v}, w={w}: <u,w>={symplectic_product_arrays(u, w)},'\
                                                ' <v,w>={symplectic_product_arrays(v, w)}")
            else:
                # If no solution found, verify this is correct by checking if either vector is zero
                # or by attempting to solve and confirming inconsistency
                is_zero_case = np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n))
                if is_zero_case:
                    # Zero vector case - correctly identified as impossible
                    continue
                else:
                    # Non-zero case - verify the system is actually inconsistent
                    # by checking that the linear system has no solution
                    A = np.zeros((2, 2 * n), dtype=int)
                    A[0, :n] = u[n:]
                    A[0, n:] = u[:n]
                    A[1, :n] = v[n:]
                    A[1, n:] = v[:n]
                    b = np.array([1, 1])

                    # The solve_gf2 function should return None for inconsistent systems
                    result = solve_gf2(A, b)
                    assert result is None, f"System should be inconsistent but solver found solution for u={u}, v={v}"

    def test_extended_constraints(self):
        """Test the extended constraint functionality."""

        # Test case 1: Simple extended constraint
        u = np.array([1, 0, 0, 0])  # X_1
        v = np.array([0, 1, 0, 0])  # X_2
        t1 = np.array([0, 0, 1, 0])  # Z_1

        # Find w with additional constraint <t1, w> = <t1, v>
        w = find_symplectic_solution_extended(u, v, [t1])
        if w is not None:
            assert self.verify_solution_extended(u, v, w, [t1]), "Extended solution should be valid"

        # Test case 2: Multiple additional constraints
        u = np.array([1, 0, 0, 0])  # X_1
        v = np.array([0, 1, 0, 0])  # X_2
        t1 = np.array([0, 0, 1, 0])  # Z_1
        t2 = np.array([0, 0, 0, 1])  # Z_2

        w = find_symplectic_solution_extended(u, v, [t1, t2])

        # TODO: To add the random test the solver needs to return None when there is no solution. Currently it
        # raises an exception and breaks the test

        # Test case 4: Random extended constraints
        # n = 5
        # for _ in range(200):
        #     u = np.random.randint(2, size=2 * n)
        #     v = np.random.randint(2, size=2 * n)

        #     # Skip zero vectors
        #     if np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n)):
        #         continue

        #     # Add 1-2 random additional constraints
        #     num_constraints = np.random.randint(1, 3)
        #     t_vectors = [np.random.randint(2, size=2 * n) for _ in range(num_constraints)]

        #     w = find_symplectic_solution_extended(u, v, t_vectors)
        #     if w is not None:
        #         assert self.verify_solution_extended(u, v, w, t_vectors), "Failed verification"

    def test_map_pauli_sum_to_target(self):

        for _ in range(100):
            # choose random properties of the system
            n = np.random.randint(2, 5)  # , 50)  # Number of qudits
            allowed_dims = [2]  # , 3, 5, 7, 11]  # allowed dimensions
            dimensions = []  # dimensions
            for _ in range(n):
                dimensions.append(int(np.random.choice(allowed_dims)))
            m = int(np.random.randint(2, 2 * n - 1))  # Number of Paulis

            # define input hamiltonian
            pl_sum = random_hamiltonian.random_pauli_hamiltonian(m, dimensions)

            basis_indices, _ = get_linear_dependencies(pl_sum.tableau, pl_sum.lcm)
            pl_sum = pl_sum[basis_indices]

            # scramble input hamiltonian to get target
            C = Circuit.from_random(n_gates=10 * n**2, dimensions=dimensions)
            target_pl_sum = C.act(pl_sum)
            # target hamiltonian
            sym_sum = pl_sum.tableau
            target_sym_sum = target_pl_sum.tableau

            check_pl_sum = pl_sum.copy()
            check_pl_sum.combine_equivalent_paulis()

            if check_pl_sum.n_paulis() != pl_sum.n_paulis():
                continue  # Skip if not mappable

            F = map_pauli_sum_to_target_tableau(sym_sum, target_sym_sum, method="auto")

            # Verify the mapping
            mapped_sym_sum = (sym_sum @ F) % pl_sum.lcm
            assert np.array_equal(mapped_sym_sum, target_sym_sum), (
                "Mapping failed. The mapped Pauli sum is:\n"
                f"{mapped_sym_sum}\n"
                " while the target Pauli sum is:\n"
                f"{target_sym_sum}\n"
                "The matrix M is:\n"
                f"{F}\n"
            )

    def test_map_pauli_sum_to_target_complete_basis_fast_path(self):
        n = 2
        input_tab = np.eye(2 * n, dtype=int)
        F_expected = np.block([
            [np.zeros((n, n), dtype=int), np.eye(n, dtype=int)],
            [np.eye(n, dtype=int), np.zeros((n, n), dtype=int)],
        ])
        target_tab = input_tab @ F_expected % 2

        F = map_pauli_sum_to_target_tableau(input_tab, target_tab, p=2)

        assert np.array_equal(F, F_expected)
        assert np.array_equal((input_tab @ F) % 2, target_tab)

    def test_map_pauli_sum_to_target_default_inverse_requires_complete_basis(self):
        input_tab = np.array([[1, 0, 0, 0], [0, 0, 1, 0]], dtype=int)
        target_tab = input_tab[[1, 0]]

        with pytest.raises(ValueError, match="complete independent Pauli basis"):
            map_pauli_sum_to_target_tableau(input_tab, target_tab, p=2)

        F = map_pauli_sum_to_target_tableau(input_tab, target_tab, p=2, method="transvection")
        assert np.array_equal((input_tab @ F) % 2, target_tab)

    def test_map_pauli_sum_to_target_complete_basis_with_dependent_row(self):
        n = 2
        p = 3
        input_basis = np.eye(2 * n, dtype=int)
        dependent_row = (input_basis[0] + 2 * input_basis[2]) % p
        input_tab = np.vstack([input_basis, dependent_row])

        F_expected = np.block([
            [np.zeros((n, n), dtype=int), np.eye(n, dtype=int)],
            [-np.eye(n, dtype=int), np.zeros((n, n), dtype=int)],
        ]) % p
        target_tab = input_tab @ F_expected % p

        F = map_pauli_sum_to_target_tableau(input_tab, target_tab, p=p)

        assert np.array_equal(F, F_expected)
        assert np.array_equal((input_tab @ F) % p, target_tab)

    def test_map_pauli_sum_to_target_complete_basis_qudit_default_inverse(self):
        n = 2
        p = 5
        input_tab = np.eye(2 * n, dtype=int)
        F_expected = np.block([
            [np.zeros((n, n), dtype=int), np.eye(n, dtype=int)],
            [-np.eye(n, dtype=int), np.zeros((n, n), dtype=int)],
        ]) % p
        target_tab = input_tab @ F_expected % p

        F = map_pauli_sum_to_target_tableau(input_tab, target_tab, p=p)

        assert np.array_equal(F, F_expected)
        assert np.array_equal((input_tab @ F) % p, target_tab)
