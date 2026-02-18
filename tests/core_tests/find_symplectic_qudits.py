import numpy as np
import galois
from sympleq.core.circuits.utils import symplectic_product_arrays
from sympleq.core.circuits.find_symplectic_qudits import build_symplectic_for_transvection, \
    check_mappable_via_clifford, intermediate_transvection_solve, find_transvection_map, find_transvection_map_solve, \
    map_paulisum_to_target_paulisum
from sympleq.core.finite_field_solvers import get_linear_dependencies
from sympleq.models import random_hamiltonian
from sympleq.core.circuits import Circuit


class TestSymplecticSolverQudits:

    def verify_intermediate_ps(self, u: np.ndarray, v: np.ndarray, w: np.ndarray, p) -> bool:
        """
        Verify that w satisfies the symplectic conditions.

        Args:
            u: Input PauliSum
            v: Output PauliSum
            w: Intermediate PauliSum
            p: dimension; prime

        Returns:
            True if <u,w> != 0 and  <w,v> != 0
        """
        if not galois.is_prime(p):
            raise NotImplementedError(f"Prime dimension expected, got p={p}")

        return (symplectic_product_arrays(u, w, p) != 0 and symplectic_product_arrays(w, v, p) != 0)

    def test_build_intermediate_ps(self):
        """Test function that checks for intermediate solution."""
        n = 10
        p = 19
        if not galois.is_prime(p):
            raise NotImplementedError(f"Prime dimension expected, got p={p}")

        for _ in range(3000):
            u = np.random.randint(p, size=2 * n)
            v = np.random.randint(p, size=2 * n)
            is_zero_case = np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n))
            if is_zero_case:
                # Zero vector case - correctly identified as impossible
                continue

            elif symplectic_product_arrays(u, v, p) != 0:
                continue

            else:
                w = build_symplectic_for_transvection(u, v, p)

                if w is not None:
                    # If a solution is found, it must be valid
                    assert self.verify_intermediate_ps(u, v, w, p), (f"Solution verification failed for u={u},'\
                                                    ' v={v}, w={w}: <u,w>={symplectic_product_arrays(u, w, p)},'\
                                                    ' <v,w>={symplectic_product_arrays(w, v, p)}")

    def test_transvection_matrix_construct(self):
        n = 10
        p = 13
        u = np.random.randint(p, size=2 * n)
        v = np.random.randint(p, size=2 * n)

        if not galois.is_prime(p):
            raise NotImplementedError(f"Prime dimension expected, got p={p}")

        for _ in range(10000):
            is_zero_case = np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n))

            if is_zero_case:
                # Zero vector case - correctly identified as impossible
                continue

            F_h = find_transvection_map(u, v, p)
            assert (u @ F_h % p == v).all(), (f"Mapping failed for u={u},'\
                                                    ' v={v}, w={p}: ")

    def test_solve_for_intermediate_ps(self):
        n = 10
        for p in range(2, 20):
            if galois.is_prime(p):
                for _ in range(3000):
                    u = np.random.randint(p, size=2 * n)
                    v = np.random.randint(p, size=2 * n)

                    is_zero_case = np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n))
                    if is_zero_case:
                        # Zero vector case - correctly identified as impossible
                        continue
                    elif symplectic_product_arrays(u, v, p) != 0:
                        continue
                    else:
                        w = intermediate_transvection_solve(u, v, p)
                        assert symplectic_product_arrays(u, w, p) != 0 and symplectic_product_arrays(
                            w, v, p) != 0, f'{symplectic_product_arrays(u, w, p), symplectic_product_arrays(w, v, p)}'

    def test_transvection_matrix_solve(self):
        n = 10
        p = 2
        u = np.random.randint(p, size=2 * n)
        v = np.random.randint(p, size=2 * n)

        if not galois.is_prime(p):
            raise NotImplementedError(f"Prime dimension expected, got p={p}")

        for _ in range(10000):
            is_zero_case = np.array_equal(u, np.zeros(2 * n)) or np.array_equal(v, np.zeros(2 * n))

            if is_zero_case:
                continue

            F_h = find_transvection_map_solve(u, v, p)
            assert (u @ F_h % p == v).all(), (f"Mapping failed for u={u},'\
                                                    ' v={v}, w={p}: ")

    def test_map_paulisum_to_paulisum(self):
        n = 10  # number of qudits
        p = 7   # prime dimension
        m = 10  # number of paulis
        dimensions = [p] * n
        for _ in range(100):
            pl_sum = random_hamiltonian.random_pauli_hamiltonian(m, dimensions)

            basis_indices, _ = get_linear_dependencies(pl_sum.tableau, int(pl_sum.lcm))
            pl_sum = pl_sum[basis_indices]

            # scramble input Hamiltonian to get target
            C = Circuit.from_random(10 * n**2, dimensions=dimensions)
            target_pl_sum = C.act(pl_sum)

            input_tab = pl_sum.tableau
            output_tab = target_pl_sum.tableau

            if check_mappable_via_clifford(input_tab, output_tab, p):
                F_total = map_paulisum_to_target_paulisum(input_tab, output_tab, p)
                assert (input_tab @ F_total % p == output_tab).all(), f'could not map for {input_tab, output_tab}'
