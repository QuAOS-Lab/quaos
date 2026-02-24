from sympleq.core.finite_field_solvers import get_linear_dependencies
import numpy as np
import galois


class TestUtils():

    def single_linear_dependence_detection(self, qudit_dim, n_independent, dim, n_dep):
        GF = galois.GF(qudit_dim)

        # Step 1: Generate known linearly independent vectors
        basis = []
        while len(basis) < n_independent:
            candidate = GF.Random(dim)
            if basis:
                A = GF(basis)
                if np.linalg.matrix_rank(np.vstack([A, candidate])) > len(basis):
                    basis.append(candidate)
            else:
                basis.append(candidate)
        B = GF(basis)

        # Step 2: Generate dependent vectors as random combinations
        dep_rows = []
        dep_coeffs = []
        for _ in range(n_dep):
            coeff = GF.Random(n_independent)
            dep = coeff @ B
            dep_rows.append(dep)
            dep_coeffs.append(coeff)

        # Full matrix
        all_rows = np.vstack([B, dep_rows])

        # Step 3: Analyze
        pivots, deps = get_linear_dependencies(all_rows, qudit_dim)

        # Step 4: Check linear span equivalence (not just index identity)
        original_basis = B
        recovered_basis = GF(all_rows[pivots, :])

        rank_match = np.linalg.matrix_rank(np.array(original_basis, dtype=int)) == \
            np.linalg.matrix_rank(np.array(recovered_basis, dtype=int))

        assert rank_match, "Recovered basis does not match original in rank!"

        # Step 5: Check each dependent vector's decomposition
        for i, coeffs in deps.items():
            if coeffs is None:
                raise AssertionError(f"Dependency decomposition failed for row {i}")
            expected = GF.Zeros(dim)
            for j, c in coeffs:
                expected += GF(c) * all_rows[j]
            actual = all_rows[i]
            assert np.array_equal(expected, actual), f"Decomposition incorrect for row {i}"

    def test_linear_dependence_detection(self):

        for qudit_dim in [2, 7, 11, 13, 17]:
            for n_independent in [5, 10]:  # number of independent paulis
                for n_qudits in [10, 20]:
                    for n_dep in [3, 5]:  # number of dependent paulis

                        self.single_linear_dependence_detection(qudit_dim, n_independent, 2 * n_qudits, n_dep)
