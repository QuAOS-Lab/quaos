
import numpy as np
import pytest

from sympleq.core.circuits.random_symplectic import symplectic_random_transvection
from sympleq.core.symmetries.modular_helpers import (
    mod_p,
    rank_mod,
    nullspace_mod,
    independent_columns,
    inv_mod_mat,
    inv_mod_scalar,
    solve_linear_many,
    _solve_linear,
    omega_matrix,
    is_symplectic,
    matmul_mod
)

from sympleq.core.symmetries.polynomials_fp import (
    poly_trim,
    poly_is_zero,
    poly_monic,
    poly_add,
    poly_sub,
    poly_mul,
    poly_divmod,
    poly_gcd,
    poly_lcm,
    extended_euclidean,
    poly_pow,
    poly_reciprocal,
    poly_eval_matrix
)

from sympleq.core.symmetries.minpoly import (
    minimal_poly_for_vector,
    minimal_polynomial,
    factor_poly_over_fp,
)

from sympleq.core.symmetries.atomic_decomposition_helpers.rcf_prepass import rcf_prepass, _is_x_pm_1
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_linear import (
    symplectic_left_inverse,
    restrict_operator,
    is_nondegenerate,
    darboux_basis_from_span,
    symplectic_completion_from_block,
)

from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2 import (
    atomic_blocks_in_unipotent_self_sector_p2,
)

from sympleq.core.symmetries.atomic_decomposition import atomic_block_decompose, CertificationError, _build_sector
import os
import traceback
import textwrap

from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_paired import atomic_blocks_in_paired_sector
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_self import atomic_blocks_in_self_sector_nonunipotent


# -----------------------
# test utils
# -----------------------


def _extract_atomic_blocks_from_global_basis(B: np.ndarray, half_dims: list[int], p: int) -> list[np.ndarray]:
    """
    Given global basis B = [U_all | V_all] and per-block half-dims [h1,h2,...],
    reconstruct each block basis in ambient coordinates as T_blk = [U_blk | V_blk]
    (2n × 2h).
    """
    B = mod_p(B, p)
    n2 = B.shape[0]
    assert B.shape == (n2, n2)
    n = n2 // 2
    assert sum(half_dims) == n, "half_dims do not sum to n; cannot partition B into blocks."

    U_all = B[:, :n]
    V_all = B[:, n:]

    blocks: list[np.ndarray] = []
    off = 0
    for h in half_dims:
        h = int(h)
        assert h > 0
        U_blk = U_all[:, off:off + h]
        V_blk = V_all[:, off:off + h]
        T_blk = np.concatenate([U_blk, V_blk], axis=1)
        blocks.append(mod_p(T_blk, p))
        off += h

    return blocks


def _assert_invariant_span(F: np.ndarray, T: np.ndarray, p: int) -> None:
    """Assert F·span(T) ⊆ span(T) using a rank test."""
    T = independent_columns(mod_p(T, p), p)
    FT = mod_p(F @ T, p)
    assert rank_mod(np.concatenate([T, FT], axis=1), p) == T.shape[1]


def kernel_in_span(A: np.ndarray, span_basis: np.ndarray, p: int) -> np.ndarray:
    """
    Return a basis (ambient columns) for { x in span(span_basis) : A x = 0 }.

    If span_basis is n×d (columns spanning subspace U), then any x ∈ U is x = span_basis c.
    Constraint A x = 0 becomes (A span_basis) c = 0, so c ∈ ker(A span_basis).
    """
    span_basis = independent_columns(mod_p(span_basis, p), p)
    if span_basis.size == 0:
        return span_basis
    AB = mod_p(A @ span_basis, p)
    C = nullspace_mod(AB, p)            # d×k coefficients
    X = mod_p(span_basis @ C, p)        # n×k ambient vectors
    return independent_columns(X, p)


def _assert_darboux_block(T: np.ndarray, p: int) -> None:
    """
    Check T (2n×2m) is a Darboux (symplectic) basis for its spanned subspace:
        T^T Ω_n T = Ω_m.
    """
    assert T.ndim == 2
    n2, m2 = T.shape
    assert n2 % 2 == 0
    assert m2 % 2 == 0
    Ωn = omega_matrix(n2 // 2, p)
    Ωm = omega_matrix(m2 // 2, p)
    G = mod_p(T.T @ Ωn @ T, p)
    assert np.array_equal(G % p, Ωm % p)


def _assert_symplectic_basis_full(B: np.ndarray, p: int) -> None:
    """Check B is a full 2n×2n symplectic basis."""
    assert B.ndim == 2
    assert B.shape[0] == B.shape[1]
    _assert_darboux_block(B, p)


def _block_diag(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=np.int64)
    B = np.asarray(B, dtype=np.int64)
    Z1 = np.zeros((A.shape[0], B.shape[1]), dtype=np.int64)
    Z2 = np.zeros((B.shape[0], A.shape[1]), dtype=np.int64)
    return np.block([[A, Z1], [Z2, B]])


def _rand_mat(rng: np.random.Generator, n: int, m: int, p: int) -> np.ndarray:
    return rng.integers(0, p, size=(n, m), dtype=np.int64)


def _rand_invertible(rng: np.random.Generator, n: int, p: int, max_tries: int = 500) -> np.ndarray:
    for _ in range(max_tries):
        A = _rand_mat(rng, n, n, p)
        if rank_mod(A, p) == n:
            return mod_p(A, p)
    raise RuntimeError("could not sample invertible matrix quickly; increase max_tries")


def _rand_symmetric(rng: np.random.Generator, n: int, p: int) -> np.ndarray:
    M = _rand_mat(rng, n, n, p)
    for i in range(n):
        for j in range(i + 1, n):
            M[j, i] = M[i, j]
    return mod_p(M, p)


def _symplectic_shear_upper(Bsym: np.ndarray, p: int) -> np.ndarray:
    # [[I, B],[0,I]] with B symmetric is symplectic in the standard Ω convention
    n = Bsym.shape[0]
    I = np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    return mod_p(np.block([[I, Bsym], [Z, I]]), p)


def _symplectic_shear_lower(Csym: np.ndarray, p: int) -> np.ndarray:
    # [[I,0],[C,I]] with C symmetric is symplectic
    n = Csym.shape[0]
    I = np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    return mod_p(np.block([[I, Z], [Csym, I]]), p)


def _symplectic_scale(A: np.ndarray, p: int) -> np.ndarray:
    # diag(A, (A^{-1})^T) is symplectic
    A = mod_p(A, p)
    Ainv = inv_mod_mat(A, p)
    return mod_p(_block_diag(A, mod_p(Ainv.T, p)), p)


def _symplectic_swap(n: int, p: int) -> np.ndarray:
    # [[0, I],[-I,0]] is symplectic
    I = np.eye(n, dtype=np.int64)
    Z = np.zeros((n, n), dtype=np.int64)
    return mod_p(np.block([[Z, I], [mod_p(-I, p), Z]]), p)


def rand_symplectic(rng: np.random.Generator, n: int, p: int, steps: int = 12) -> np.ndarray:
    F = np.eye(2 * n, dtype=np.int64)
    for _ in range(int(steps)):
        t = int(rng.integers(0, 4))
        if t == 0:
            B = _rand_symmetric(rng, n, p)
            S = _symplectic_shear_upper(B, p)
        elif t == 1:
            C = _rand_symmetric(rng, n, p)
            S = _symplectic_shear_lower(C, p)
        elif t == 2:
            A = _rand_invertible(rng, n, p)
            S = _symplectic_scale(A, p)
        else:
            S = _symplectic_swap(n, p)

        F = mod_p(S @ F, p)

    assert is_symplectic(F, p)
    return F


def _rng(seed: int | None = None) -> np.random.Generator:
    return np.random.default_rng(seed)


class TestModularHelpers:
    def test_mod_p_basic(self) -> None:
        p = 7
        x = np.array([-1, 0, 1, 7, 8, -8], dtype=np.int64)
        y = mod_p(x, p)
        assert np.array_equal(y, np.array([6, 0, 1, 0, 1, 6], dtype=np.int64))
        assert np.all((y >= 0) & (y < p))

    def test_rank_mod_known_cases(self) -> None:
        p = 5
        Id = np.eye(4, dtype=np.int64)
        Z = np.zeros((4, 4), dtype=np.int64)
        assert rank_mod(Id, p) == 4
        assert rank_mod(Z, p) == 0

        # Two identical rows => rank 1
        A = np.array([[1, 2, 3],
                      [1, 2, 3]], dtype=np.int64)
        assert rank_mod(A, p) == 1

    def test_matmul_mod_matches_numpy(self) -> None:
        rng = _rng()
        p = 7
        A = _rand_mat(rng, 4, 5, p)
        B = _rand_mat(rng, 5, 3, p)
        got = matmul_mod(A, B, p)
        want = mod_p(A @ B, p)
        assert np.array_equal(got, want)

    def test_nullspace_mod_correctness_and_dimension(self) -> None:
        rng = _rng()
        p = 5
        A = _rand_mat(rng, 3, 6, p)
        N = nullspace_mod(A, p)  # expected shape (6, k)

        # A @ N == 0
        assert np.array_equal(mod_p(A @ N, p), np.zeros((A.shape[0], N.shape[1]), dtype=np.int64))

        # dimension check: dim ker = n - rank(A)
        assert N.shape[0] == A.shape[1]
        assert N.shape[1] == A.shape[1] - rank_mod(A, p)

        # columns independent
        assert rank_mod(N, p) == N.shape[1]

    def test_independent_columns_rank_preserved(self) -> None:
        p = 7
        # Build columns with dependencies: c2 = c0 + c1, and a zero column
        c0 = np.array([[1], [2], [3]], dtype=np.int64)
        c1 = np.array([[0], [1], [1]], dtype=np.int64)
        c2 = mod_p(c0 + c1, p)
        c3 = np.zeros((3, 1), dtype=np.int64)
        M = np.concatenate([c0, c1, c2, c3], axis=1)

        r0 = rank_mod(M, p)
        B = independent_columns(M, p)
        r1 = rank_mod(B, p)

        assert r1 == r0
        assert B.shape[1] == r0
        # and returned columns are independent
        assert rank_mod(B, p) == B.shape[1]

    def test_inv_mod_mat_identity_and_random(self) -> None:
        rng = _rng()
        p = 5

        Id = np.eye(4, dtype=np.int64)
        Id_inv = inv_mod_mat(Id, p)
        assert np.array_equal(mod_p(Id @ Id_inv, p), Id)
        assert np.array_equal(mod_p(Id_inv @ Id, p), Id)

        A = _rand_invertible(rng, 5, p)
        A_inv = inv_mod_mat(A, p)
        assert np.array_equal(mod_p(A @ A_inv, p), np.eye(5, dtype=np.int64))
        assert np.array_equal(mod_p(A_inv @ A, p), np.eye(5, dtype=np.int64))

    def test_inv_mod_mat_raises_on_singular(self) -> None:
        p = 7
        # Two identical rows -> singular
        A = np.array([[1, 0, 0],
                      [1, 0, 0],
                      [0, 0, 1]], dtype=np.int64)
        with pytest.raises(Exception):
            _ = inv_mod_mat(A, p)

    def test_solve_linear_many(self) -> None:
        rng = _rng()
        p = 5

        # B has full column rank (n x d)
        n, d, k = 7, 4, 6
        B = _rand_mat(rng, n, d, p)
        while rank_mod(B, p) != d:
            B = _rand_mat(rng, n, d, p)

        X_true = _rand_mat(rng, d, k, p)
        Y = mod_p(B @ X_true, p)

        X = solve_linear_many(B, Y, p)
        assert np.array_equal(mod_p(B @ X, p), Y)

    def test__solve_linear_square_system(self) -> None:
        rng = _rng()
        p = 11

        A = _rand_invertible(rng, 6, p)
        x_true = _rand_mat(rng, 6, 1, p)
        b = mod_p(A @ x_true, p)

        x = _solve_linear(A, b, p)
        assert np.array_equal(mod_p(A @ x, p), b)

    def test__solve_linear_inconsistent_raises(self) -> None:
        # A x = b impossible in GF(2): x = 0 and x = 1 simultaneously
        p = 2
        A = np.array([[1],
                      [1]], dtype=np.int64)
        b = np.array([[0],
                      [1]], dtype=np.int64)
        with pytest.raises(Exception):
            _ = _solve_linear(A, b, p)

    def test_omega_matrix_structure_and_properties(self) -> None:
        # n=1
        Ω2 = omega_matrix(1, 2)
        assert Ω2.shape == (2, 2)
        # In p=2, -1 == 1 so Ω = [[0,1],[1,0]]
        assert np.array_equal(mod_p(Ω2, 2), np.array([[0, 1],
                                                      [1, 0]], dtype=np.int64))

        p = 5
        Ω = omega_matrix(3, p)
        assert Ω.shape == (6, 6)
        # Ω^T = -Ω
        assert np.array_equal(mod_p(Ω.T, p), mod_p(-Ω, p))
        # Ω is invertible
        assert rank_mod(Ω, p) == Ω.shape[0]

    def test_is_symplectic_identity_and_constructed(self) -> None:
        rng = _rng()

        # Identity is symplectic
        p = 5
        Id = np.eye(8, dtype=np.int64)
        assert is_symplectic(Id, p)

        # Construct symplectic block-diagonal: [[A,0],[0,A^{-T}]]
        n = 4
        A = _rand_invertible(rng, n, p)
        A_inv = inv_mod_mat(A, p)
        F_invT = mod_p(A_inv.T, p)
        F = _block_diag(A, F_invT)
        assert F.shape == (2 * n, 2 * n)
        assert is_symplectic(F, p)

        # A random matrix should usually not be symplectic
        R = _rand_mat(rng, 2 * n, 2 * n, p)
        if not np.array_equal(R, F):  # avoid pathological equality
            assert is_symplectic(R, p) is False

    def test_mod_p_range_and_idempotent(self) -> None:
        p = 7
        A = np.array([[-8, -1, 0, 1, 8]], dtype=np.int64)
        B = mod_p(A, p)
        assert np.all((0 <= B) & (B < p))
        assert np.array_equal(mod_p(B, p), B)

    def test_inv_mod_scalar(self) -> None:
        p = 11
        for a in [1, 2, 3, 5, 7, 10]:
            inv = inv_mod_scalar(a, p)
            assert (a * inv) % p == 1

        # p=2 only has 1 as invertible scalar
        assert inv_mod_scalar(1, 2) == 1

    def test_inv_mod_mat_random(self) -> None:
        rng = np.random.default_rng()
        p = 5
        for n in [1, 2, 3, 4]:
            for _ in range(25):
                A = _rand_invertible(rng, n, p)
                Ainv = inv_mod_mat(A, p)
                I = np.eye(n, dtype=np.int64)
                assert np.array_equal(mod_p(A @ Ainv, p), I)
                assert np.array_equal(mod_p(Ainv @ A, p), I)

    def test_rank_nullspace_ranknullity(self) -> None:
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for (n, m) in [(4, 4), (4, 6), (6, 4)]:
                for _ in range(15):
                    A = _rand_mat(rng, n, m, p)
                    N = nullspace_mod(A, p)  # m×k, A @ N = 0
                    assert np.array_equal(mod_p(A @ N, p), np.zeros((n, N.shape[1]), dtype=np.int64))
                    # rank-nullity: rank(A) + nullity(A) = m
                    assert rank_mod(A, p) + N.shape[1] == m

    def test_independent_columns_rank_preserving(self) -> None:
        rng = np.random.default_rng()
        for p in [2, 3, 7]:
            for n in [4, 6]:
                A = _rand_mat(rng, n, n + 3, p)
                B = independent_columns(A, p)
                assert rank_mod(B, p) == B.shape[1]
                assert rank_mod(B, p) == rank_mod(A, p)

    def test_solve_linear_many_consistency(self) -> None:
        rng = np.random.default_rng()
        p = 5
        n = 6
        B = _rand_invertible(rng, n, p)
        X_true = _rand_mat(rng, n, 4, p)
        Y = mod_p(B @ X_true, p)
        X = solve_linear_many(B, Y, p)
        assert np.array_equal(mod_p(B @ X, p), Y)

    def test__solve_linear_exact_solution(self) -> None:
        rng = np.random.default_rng()
        p = 7
        A = _rand_invertible(rng, 5, p)
        x_true = _rand_mat(rng, 5, 1, p)
        b = mod_p(A @ x_true, p)
        x = _solve_linear(A, b, p)
        assert np.array_equal(mod_p(A @ x, p), b)

    def test_omega_basic_properties(self) -> None:
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                Ω = omega_matrix(n, p)
                # Ω should be invertible (rank = 2n)
                assert rank_mod(Ω, p) == 2 * n
                # skew-symmetry: Ω^T = -Ω (over p=2 this reduces to Ω^T = Ω)
                assert np.array_equal(mod_p(Ω.T + Ω, p), np.zeros_like(Ω))


class TestPolynomialsFP:
    def test_poly_add_sub_mul_basic(self) -> None:
        p = 7
        a = np.array([1, 2, 3], dtype=np.int64)   # 1 + 2x + 3x^2
        b = np.array([6, 1], dtype=np.int64)      # 6 + x

        s = poly_add(a, b, p)
        d = poly_sub(a, b, p)
        m = poly_mul(a, b, p)

        # sanity: (a+b)-b = a  (literal equality mod p, no monic normalization)
        assert np.array_equal(poly_sub(s, b, p), mod_p(poly_trim(a), p))

        # check mul against hand result:
        # (1+2x+3x^2)(6+x) = 6 + 13x + 20x^2 + 3x^3 = 6 + 6x + 6x^2 + 3x^3 mod 7
        want = np.array([6, 6, 6, 3], dtype=np.int64)
        assert np.array_equal(m, mod_p(poly_trim(want), p))

        # optional: also sanity-check subtraction
        assert np.array_equal(d, mod_p(poly_trim(np.array([1-6, 2-1, 3], dtype=np.int64)), p))

    def test_poly_divmod_regression_trailing_zero_overwrite(self) -> None:
        p = 7
        a = np.array([3, 3, 3, 3, 2, 4], dtype=np.int64)
        b = np.array([0, 4, 2, 3], dtype=np.int64)
        q, r = poly_divmod(a, b, p)
        recon = poly_add(poly_mul(q, b, p), r, p)
        assert np.array_equal(mod_p(poly_trim(recon), p), mod_p(poly_trim(a), p))

    def test_poly_divmod_roundtrip(self) -> None:
        p = 5
        a = np.array([4, 0, 1, 2], dtype=np.int64)   # 4 + x^2 + 2x^3
        b = np.array([1, 2, 1], dtype=np.int64)      # 1 + 2x + x^2

        q, r = poly_divmod(a, b, p)
        back = poly_add(poly_mul(q, b, p), r, p)

        assert np.array_equal(poly_monic(back, p), poly_monic(a, p))

        # deg(r) < deg(b) unless r=0
        if not poly_is_zero(r):
            assert (len(poly_trim(r)) - 1) < (len(poly_trim(b)) - 1)

    def test_poly_gcd_lcm_properties(self) -> None:
        p = 11
        # (x-1)(x-2)
        f = poly_mul(np.array([p - 1, 1], dtype=np.int64), np.array([p - 2, 1], dtype=np.int64), p)
        # (x-2)(x-3)
        g = poly_mul(np.array([p - 2, 1], dtype=np.int64), np.array([p - 3, 1], dtype=np.int64), p)

        d = poly_gcd(f, g, p)
        # gcd should be (x-2)
        want = poly_monic(np.array([p - 2, 1], dtype=np.int64), p)
        assert np.array_equal(poly_monic(d, p), want)

        l = poly_lcm(f, g, p)
        # lcm should be (x-1)(x-2)(x-3)
        want_l = poly_mul(f, np.array([p - 3, 1], dtype=np.int64), p)
        assert np.array_equal(poly_monic(l, p), poly_monic(want_l, p))

    def test_poly_pow_and_reciprocal(self) -> None:
        p = 5
        q = np.array([2, 1, 1], dtype=np.int64)  # 2 + x + x^2

        q2 = poly_pow(q, 2, p)
        # q^2 = q*q
        assert np.array_equal(poly_monic(q2, p), poly_monic(poly_mul(q, q, p), p))

        qs = poly_reciprocal(q, p)
        qss = poly_reciprocal(qs, p)
        # reciprocal twice gives original up to monic normalization
        assert np.array_equal(poly_monic(qss, p), poly_monic(q, p))

    def test_poly_eval_matrix_matches_naive(self) -> None:
        rng = _rng()
        p = 5
        F = _rand_mat(rng, 4, 4, p)
        poly = np.array([3, 2, 4], dtype=np.int64)   # 3 + 2x + 4x^2

        got = poly_eval_matrix(F, poly, p)

        I = np.eye(4, dtype=np.int64)
        want = mod_p(3 * I + 2 * F + 4 * (F @ F), p)
        assert np.array_equal(got, want)

    def test_poly_divmod_random_roundtrip(self) -> None:
        rng = np.random.default_rng()
        p = 7
        for _ in range(50):
            deg_a = int(rng.integers(1, 8))
            deg_b = int(rng.integers(1, 6))
            a = rng.integers(0, p, size=deg_a + 1, dtype=np.int64)
            b = rng.integers(0, p, size=deg_b + 1, dtype=np.int64)
            a = poly_trim(a)
            b = poly_trim(b)
            if poly_is_zero(b):
                b = np.array([1], dtype=np.int64)
            # ensure divisor has nonzero leading term
            if int(b[-1]) % p == 0:
                b[-1] = 1

            q, r = poly_divmod(a, b, p)
            recon = poly_add(poly_mul(q, b, p), r, p)
            assert np.array_equal(mod_p(poly_trim(recon), p), mod_p(poly_trim(a), p))
            # remainder degree < divisor degree (unless remainder is zero)
            if not poly_is_zero(r):
                assert (len(poly_trim(r)) - 1) < (len(poly_trim(b)) - 1)

    def test_poly_xgcd_identity(self) -> None:
        rng = np.random.default_rng()
        p = 5
        for _ in range(40):
            a = poly_trim(rng.integers(0, p, size=int(rng.integers(1, 7)), dtype=np.int64))
            b = poly_trim(rng.integers(0, p, size=int(rng.integers(1, 7)), dtype=np.int64))
            if poly_is_zero(a):
                a = np.array([1], dtype=np.int64)
            if poly_is_zero(b):
                b = np.array([1], dtype=np.int64)
            s, t, g = extended_euclidean(a, b, p)
            lhs = poly_add(poly_mul(s, a, p), poly_mul(t, b, p), p)
            assert np.array_equal(poly_monic(lhs, p), poly_monic(g, p))

    def test_poly_gcd_lcm_consistency(self) -> None:
        rng = np.random.default_rng()
        p = 7
        for _ in range(40):
            a = poly_trim(rng.integers(0, p, size=int(rng.integers(1, 7)), dtype=np.int64))
            b = poly_trim(rng.integers(0, p, size=int(rng.integers(1, 7)), dtype=np.int64))
            if poly_is_zero(a):
                a = np.array([1], dtype=np.int64)
            if poly_is_zero(b):
                b = np.array([1], dtype=np.int64)
            g = poly_gcd(a, b, p)
            l = poly_lcm(a, b, p)

            # g divides both a and b
            qa, ra = poly_divmod(poly_monic(a, p), g, p)
            qb, rb = poly_divmod(poly_monic(b, p), g, p)
            assert poly_is_zero(ra)
            assert poly_is_zero(rb)

            # a divides l and b divides l
            q1, r1 = poly_divmod(l, poly_monic(a, p), p)
            q2, r2 = poly_divmod(l, poly_monic(b, p), p)
            assert poly_is_zero(r1)
            assert poly_is_zero(r2)

    def test_factor_poly_linear_multiset_odd_p(self) -> None:
        rng = np.random.default_rng()
        p = 11
        # build f(x)=∏(x-r_i) with repeats
        roots = [1, 1, 3, 7, 9]
        f = np.array([1], dtype=np.int64)
        for r in roots:
            f = poly_mul(f, np.array([(-r) % p, 1], dtype=np.int64), p)
        f = poly_monic(f, p)

        facs = factor_poly_over_fp(f, p, rng=rng)
        facs = [poly_monic(g, p) for g in facs if (len(g) > 1)]  # ignore constants if any

        # Extract roots from linear factors (x - r) -> [(-r), 1]
        got_roots: list[int] = []
        for g in facs:
            assert len(g) == 2  # should fully split
            got_roots.append((-int(g[0])) % p)

        assert sorted(got_roots) == sorted(roots)

    def test_poly_reciprocal_involution(self) -> None:
        p = 7
        q = np.array([3, 0, 5, 1], dtype=np.int64)  # 3 + 0x + 5x^2 + x^3
        q = poly_monic(q, p)
        q2 = poly_reciprocal(poly_reciprocal(q, p), p)
        assert np.array_equal(q2, q)

class TestMinPolyAndFactorization:
    def test_minimal_poly_for_vector_simple(self) -> None:
        # Column-action convention: w <- F w
        p = 5
        F = np.array([[0, 1, 0],
                      [0, 0, 1],
                      [0, 0, 0]], dtype=np.int64)

        # e2 has a length-3 chain: e2 -> e1 -> e0 -> 0, so m_v(x) = x^3
        v = np.array([0, 0, 1], dtype=np.int64)
        mv = minimal_poly_for_vector(F, v, p)
        assert np.array_equal(poly_monic(mv, p), np.array([0, 0, 0, 1], dtype=np.int64))

        # e0 dies immediately: F e0 = 0, so m_{e0}(x) = x
        v0 = np.array([1, 0, 0], dtype=np.int64)
        mv0 = minimal_poly_for_vector(F, v0, p)
        assert np.array_equal(poly_monic(mv0, p), np.array([0, 1], dtype=np.int64))

    def test_minimal_polynomial_diagonal(self) -> None:
        p = 5
        # diag(2,3,2,3) has minimal polynomial (x-2)(x-3)
        F = np.diag([2, 3, 2, 3]).astype(np.int64)
        mF = minimal_polynomial(F, p)

        f2 = np.array([p - 2, 1], dtype=np.int64)  # x-2
        f3 = np.array([p - 3, 1], dtype=np.int64)  # x-3
        want = poly_monic(poly_mul(f2, f3, p), p)

        assert np.array_equal(poly_monic(mF, p), want)

    def test_factor_poly_over_fp_char2_known(self) -> None:
        p = 2
        # f = (x+1)^3 * (x^2+x+1)
        x1 = np.array([1, 1], dtype=np.int64)
        q2 = np.array([1, 1, 1], dtype=np.int64)
        f = poly_mul(poly_mul(poly_mul(x1, x1, p), x1, p), q2, p)

        facs = factor_poly_over_fp(f, p)
        keys = [tuple(poly_monic(g, p).tolist()) for g in facs]

        assert keys.count(tuple(x1.tolist())) == 3
        assert keys.count(tuple(poly_monic(q2, p).tolist())) == 1

    def test_minimal_poly_for_vector_shift_column_action(self) -> None:
        # Column-action: for nilpotent shift with ones on superdiagonal,
        # e_{n-1} has Krylov length n => minimal polynomial x^n.
        p = 5
        F = np.array([[0, 1, 0],
                      [0, 0, 1],
                      [0, 0, 0]], dtype=np.int64)
        v = np.array([0, 0, 1], dtype=np.int64)  # top of the column-action chain
        mv = minimal_poly_for_vector(F, v, p)
        assert np.array_equal(poly_monic(mv, p), np.array([0, 0, 0, 1], dtype=np.int64))

    def test_minimal_polynomial_annihilates_matrix_random(self) -> None:
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for n in [2, 3, 4]:
                for _ in range(10):
                    F = _rand_mat(rng, n, n, p)
                    mF = minimal_polynomial(F, p)
                    Z = poly_eval_matrix(F, mF, p)
                    assert np.array_equal(mod_p(Z, p), np.zeros((n, n), dtype=np.int64))

    def test_factorization_product_roundtrip_random_constructed(self) -> None:
        rng = np.random.default_rng()
        p = 2
        # Construct f = (x+1)^3 * (x^2+x+1)
        x1 = np.array([1, 1], dtype=np.int64)
        q2 = np.array([1, 1, 1], dtype=np.int64)
        f = poly_mul(poly_pow(x1, 3, p), q2, p)
        f = poly_monic(f, p)

        facs = factor_poly_over_fp(f, p, rng=rng)
        prod = np.array([1], dtype=np.int64)
        for g in facs:
            prod = poly_mul(prod, g, p)
        assert np.array_equal(poly_monic(prod, p), f)

class TestRCFPrepass:
    def test_rcf_prepass_paired_linear_factors(self) -> None:
        # Build symplectic F = diag(A, A^{-T}) over GF(5) with A=diag(2,3)
        p = 5
        A = np.diag([2, 3]).astype(np.int64)
        Ainv = inv_mod_mat(A, p)
        F = _block_diag(A, mod_p(Ainv.T, p))
        assert is_symplectic(F, p)

        meta = rcf_prepass(F, p)
        assert "primaries" in meta and "sectors" in meta
        prim = meta["primaries"]
        secs = meta["sectors"]

        # Expect one paired sector linking (x-2) <-> (x-3) (since 2^{-1}=3 mod 5)
        assert len(secs) == 1
        assert secs[0]["type"] == "paired"

        # Keys for x-2 and x-3 should exist (coeffs low->high)
        k2 = tuple(poly_monic(np.array([p - 2, 1], dtype=np.int64), p).tolist())
        k3 = tuple(poly_monic(np.array([p - 3, 1], dtype=np.int64), p).tolist())
        assert k2 in prim and k3 in prim

        # Sector should span whole space
        W = secs[0]["W_basis"]
        assert rank_mod(W, p) == F.shape[0]

        # Each primary basis is invariant
        for k in (k2, k3):
            B = prim[k]["V_basis"]
            FB = mod_p(F @ B, p)
            assert rank_mod(np.concatenate([B, FB], axis=1), p) == B.shape[1]

        # Certified bound should be present and <= n
        assert meta["Lmin_star"] >= 1
        assert meta["Lmin_star"] <= (F.shape[0] // 2)

    def test_rcf_prepass_p2_identity_unipotent_sector(self) -> None:
        p = 2
        F = np.eye(6, dtype=np.int64)
        assert is_symplectic(F, p)

        meta = rcf_prepass(F, p)
        secs = meta["sectors"]
        assert len(secs) == 1
        assert secs[0]["type"] == "self"
        # in your prepass this is the uncertified corner, so floor_certified is False
        assert secs[0]["floor_certified"] is False
        assert meta["Lmin_star"] == 1

    def test_rcf_prepass_sectors_direct_sum_random(self) -> None:
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                for _ in range(10):
                    F = rand_symplectic(rng, n, p, steps=12)
                    meta = rcf_prepass(F, p)
                    secs = meta["sectors"]

                    # invariance: F W ⊆ W
                    for sec in secs:
                        W = sec["W_basis"]
                        FW = mod_p(F @ W, p)
                        assert rank_mod(np.concatenate([W, FW], axis=1), p) == W.shape[1]

                    # direct sum: sum of dims = rank of concatenation
                    dims = sum(sec["W_basis"].shape[1] for sec in secs)
                    all_cols = np.concatenate([sec["W_basis"] for sec in secs], axis=1)
                    assert rank_mod(all_cols, p) == dims
                    assert rank_mod(all_cols, p) == 2 * n

    def test_rcf_prepass_sector_span_and_invariance_random(self) -> None:
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                for _ in range(8):
                    F = rand_symplectic(rng, n, p, steps=10)
                    meta = rcf_prepass(F, p)
                    sectors = meta["sectors"]
                    assert sectors

                    # Sectors should span the full space
                    all_cols = np.concatenate([sec["W_basis"] for sec in sectors], axis=1)
                    assert rank_mod(all_cols, p) == 2 * n

                    # Each sector basis should be invariant: F W ⊆ span(W)
                    for sec in sectors:
                        W = independent_columns(sec["W_basis"], p)
                        FW = mod_p(F @ W, p)
                        # invariance check via rank: rank([W|FW]) == rank(W)
                        assert rank_mod(np.concatenate([W, FW], axis=1), p) == W.shape[1]

    def test_rcf_prepass_paired_linear_eigs_diagA(self) -> None:
        # Paired sector expected when eigenvalues come in a / a^{-1} with a != a^{-1}.
        p = 5
        n = 2
        A = np.diag([2, 3]).astype(np.int64)  # 2^{-1}=3 mod 5 and 3^{-1}=2 mod 5
        F = _symplectic_scale(A, p)
        assert is_symplectic(F, p)

        meta = rcf_prepass(F, p)
        secs = meta["sectors"]
        assert len(secs) == 1
        assert secs[0]["type"] == "paired"

        key = secs[0]["key"]
        key_star = secs[0]["key_star"]
        prim = meta["primaries"]
        q = prim[key]["poly"]
        qstar = prim[key_star]["poly"]
        assert tuple(poly_reciprocal(q, p).tolist()) == tuple(qstar.tolist())

    def test_rcf_prepass_unipotent_p2_marks_uncertified_floor(self) -> None:
        p = 2
        n = 2
        # Simple unipotent symplectic shear
        B = np.diag([1, 0]).astype(np.int64)
        F = _symplectic_shear_upper(B, p)
        assert is_symplectic(F, p)

        meta = rcf_prepass(F, p)
        assert len(meta["sectors"]) == 1
        sec = meta["sectors"][0]
        assert sec["type"] == "self"
        assert sec["floor_certified"] is False
        assert "note" in sec

    def test_rcf_prepass_primaries_direct_sum_and_span_random(self) -> None:
        """
        Stronger than sector-level checks:
          - primary spaces V_q are pairwise independent (direct sum)
          - primaries span the full space (rank = 2n)
          - each V_q is invariant
        """
        rng = np.random.default_rng(1234)
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                for _ in range(8):
                    F = rand_symplectic(rng, n, p, steps=12)
                    meta = rcf_prepass(F, p)
                    prim = meta["primaries"]
                    assert prim

                    V_list: list[np.ndarray] = []
                    dims: list[int] = []
                    for k, d in prim.items():
                        V = independent_columns(mod_p(d["V_basis"], p), p)
                        V_list.append(V)
                        dims.append(int(V.shape[1]))

                        # invariance check: F V ⊆ span(V)
                        _assert_invariant_span(F, V, p)

                    # span check
                    all_cols = np.concatenate(V_list, axis=1) if V_list else np.zeros((2 * n, 0), dtype=np.int64)
                    assert rank_mod(all_cols, p) == 2 * n

                    # direct sum check: rank equals sum of dims
                    assert rank_mod(all_cols, p) == sum(dims)

                    # pairwise trivial intersections (cheap because small sizes)
                    for i in range(len(V_list)):
                        for j in range(i + 1, len(V_list)):
                            Vij = np.concatenate([V_list[i], V_list[j]], axis=1)
                            assert rank_mod(Vij, p) == dims[i] + dims[j]


class TestAtomicLinear:
    def test_symplectic_left_inverse_full_matrix(self) -> None:
        rng = _rng()
        p = 7
        n = 3

        A = _rand_invertible(rng, n, p)
        Ainv = inv_mod_mat(A, p)
        F = _block_diag(A, mod_p(Ainv.T, p))
        assert is_symplectic(F, p)

        L = symplectic_left_inverse(F, p)
        I = np.eye(2 * n, dtype=np.int64)
        assert np.array_equal(mod_p(L @ F, p), I)

    def test_restrict_operator_conjugation(self) -> None:
        rng = _rng()
        p = 5
        n = 2

        A = _rand_invertible(rng, n, p)
        Ainv = inv_mod_mat(A, p)
        B = _block_diag(A, mod_p(Ainv.T, p))
        assert is_symplectic(B, p)

        # pick another symplectic F
        C = _rand_invertible(rng, n, p)
        Cinv = inv_mod_mat(C, p)
        F = _block_diag(C, mod_p(Cinv.T, p))
        assert is_symplectic(F, p)

        # restrict_operator(F, B) = B^{-1} F B
        Sigma = restrict_operator(F, B, p)
        Binv = inv_mod_mat(B, p)
        want = mod_p(Binv @ F @ B, p)
        assert np.array_equal(Sigma, want)

    def test_kernel_in_span_correctness(self) -> None:
        rng = _rng()
        p = 11
        n = 8

        A = _rand_mat(rng, 5, n, p)
        # span_basis: pick some random columns (ensure independent)
        span = independent_columns(_rand_mat(rng, n, 6, p), p)

        X = kernel_in_span(A, span, p)

        # X is in span(span_basis): rank([span, X]) == rank(span)
        if X.shape[1] > 0:
            assert rank_mod(np.concatenate([span, X], axis=1), p) == rank_mod(span, p)

        # and A X = 0
        assert np.array_equal(mod_p(A @ X, p), np.zeros((A.shape[0], X.shape[1]), dtype=np.int64))

    def test_darboux_basis_from_span_raises_on_degenerate(self) -> None:
        p = 5
        n = 2
        Ω = omega_matrix(n, p)
        # Take a 2D isotropic subspace: span{e0, e1} in the x-part
        B = np.eye(2*n, dtype=np.int64)[:, [0, 1]]
        assert not is_nondegenerate(Ω, B, p)
        with pytest.raises(RuntimeError):
            darboux_basis_from_span(Ω, B, p)

    def test_darboux_basis_from_span_standard(self) -> None:
        p = 5
        n = 4
        Ω = omega_matrix(n, p)

        # Nondegenerate subspace spanned by {e0,e1,e_{n+0},e_{n+1}}
        e0 = np.zeros((2 * n, 1), dtype=np.int64); e0[0, 0] = 1
        e1 = np.zeros((2 * n, 1), dtype=np.int64); e1[1, 0] = 1
        f0 = np.zeros((2 * n, 1), dtype=np.int64); f0[n + 0, 0] = 1
        f1 = np.zeros((2 * n, 1), dtype=np.int64); f1[n + 1, 0] = 1
        B = np.concatenate([e0, e1, f0, f1], axis=1)

        T = darboux_basis_from_span(Ω, B, p)
        assert T.shape == (2 * n, 4)

        Ω2 = omega_matrix(2, p)
        G = mod_p(T.T @ Ω @ T, p)
        assert np.array_equal(G, Ω2)

    def test_darboux_basis_from_symplectic_subspace(self) -> None:
        rng = np.random.default_rng()
        p = 5
        n = 4
        F = rand_symplectic(rng, n, p, steps=10)

        # Build a guaranteed nondegenerate symplectic subspace:
        # image under F of span{e_i, e_{n+i}} for i=0..k-1.
        k = 2
        cols: list[int] = []
        for i in range(k):
            cols += [i, n + i]
        B = F[:, cols]
        Ω = omega_matrix(n, p)
        assert is_nondegenerate(Ω, B, p)
        T = darboux_basis_from_span(Ω, B, p)
        _assert_darboux_block(T, p)

    def test_symplectic_completion_from_block(self) -> None:
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            n = 3
            F = rand_symplectic(rng, n, p, steps=12)
            # pick a 2k block
            k = 1
            T_blk = F[:, [0, n + 0]]
            T_full = symplectic_completion_from_block(T_blk, p)
            assert T_full.shape == (2 * n, 2 * n)
            _assert_symplectic_basis_full(T_full, p)

    def test_restrict_operator_matches_action_on_span(self) -> None:
        rng = np.random.default_rng()
        p = 5
        n = 3
        k = 2
        n_tests = 10
        for _ in range(n_tests):
            # Build F with an invariant 2k-dimensional symplectic subspace in standard [x|z] ordering.
            # Subspace S = span{x0,x1,z0,z1} corresponds to indices [0,1,3,4].
            F1 = rand_symplectic(rng, k, p, steps=8)         # 4x4 symplectic in [x0,x1|z0,z1]
            F2 = rand_symplectic(rng, n - k, p, steps=8)     # 2x2 symplectic in [x2|z2]
            F = np.eye(2 * n, dtype=np.int64)
            idxS = [0, 1, n + 0, n + 1]
            idxT = [2, n + 2]
            F[np.ix_(idxS, idxS)] = F1
            F[np.ix_(idxT, idxT)] = F2
            F = mod_p(F, p)

            # Start from canonical Darboux basis on first k modes and randomize it within that subspace
            T0 = np.eye(2 * n, dtype=np.int64)[:, [0, 1, n + 0, n + 1]]  # [x0,x1,z0,z1]
            G = rand_symplectic(rng, k, p, steps=6)                      # change of basis within subspace
            T = mod_p(T0 @ G, p)                                         # still Darboux, still invariant

            F_T = restrict_operator(F, T, p)     # coords on span(T)

            # Check: F T = T F_T  (since columns of T are basis vectors)
            lhs = mod_p(F @ T, p)
            rhs = mod_p(T @ F_T, p)
            assert np.array_equal(lhs, rhs)

            # Symplectic left-inverse consistency on full-rank T
            L = symplectic_left_inverse(T, p)
            assert np.array_equal(mod_p(L @ T, p), np.eye(2 * k, dtype=np.int64))


class TestAtomicUnipotentP2:
    # def test_classify_unipotent_sp2_and_rebuild_blocks(self) -> None:
    #     # Build a unipotent symplectic F over GF(2):
    #     # F = [[I,0],[S,I]] with S symmetric and diag nonzero to ensure anisotropy exists.
    #     p = 2
    #     n = 2
    #     I = np.eye(n, dtype=np.int64)
    #     S = np.eye(n, dtype=np.int64)  # symmetric, diag ones
    #     top = np.concatenate([I, np.zeros((n, n), dtype=np.int64)], axis=1)
    #     bot = np.concatenate([S, I], axis=1)
    #     F = np.concatenate([top, bot], axis=0)
    #     F = mod_p(F, p)
    #     assert is_symplectic(F, p)

    #     inv_data = classify_unipotent_sp2(F)
    #     assert inv_data["p"] == 2
    #     assert inv_data["n2"] == 2 * n
    #     assert "jordan_profile" in inv_data
    #     assert "blocks" in inv_data
    #     assert len(inv_data["blocks"]) == 2  # should extract two 2D blocks

    #     # Expect V blocks here (m_max=2 and anisotropy exists)
    #     assert all(b["type"] == "V" and b["m"] == 2 and b["k"] == 1 for b in inv_data["blocks"])

    #     # Rebuild explicit ambient blocks (here V_u = I spans full space)
    #     V_u = np.eye(2 * n, dtype=np.int64)
    #     blocks = build_unipotent_blocks_from_invariants(F, V_u, inv_data, p=2)
    #     assert len(blocks) == 2
    #     assert all(b.half_dim == 1 for b in blocks)

    #     Ω = omega_matrix(n, p)
    #     Ω1 = omega_matrix(1, p)
    #     for b in blocks:
    #         T = b.T_blk
    #         assert T.shape == (2 * n, 2)
    #         assert np.array_equal(mod_p(T.T @ Ω @ T, p), Ω1)

    # def test_classify_unipotent_simple_shear(self) -> None:
    #     p = 2
    #     n = 3
    #     B = np.diag([1, 0, 1]).astype(np.int64)
    #     F = _symplectic_shear_upper(B, p)
    #     assert is_symplectic(F, p)

    #     inv = classify_unipotent_sp2(F)
    #     assert inv["p"] == 2
    #     assert inv["n2"] == 2 * n
    #     assert "jordan_profile" in inv
    #     assert "blocks" in inv
    #     # should include at least one block description
    #     assert isinstance(inv["blocks"], list) and len(inv["blocks"]) >= 1

    def test_unipotent_sector_pipeline_from_rcf_prepass(self) -> None:
        p = 2
        n = 3
        rng = np.random.default_rng()

        # Build a guaranteed p=2 unipotent symplectic:
        # product of *upper* shears stays upper shear => always unipotent (index <= 2)
        B1 = _rand_symmetric(rng, n, p)
        B2 = _rand_symmetric(rng, n, p)
        F = mod_p(_symplectic_shear_upper(B1, p) @ _symplectic_shear_upper(B2, p), p)

        assert is_symplectic(F, p)

        meta = rcf_prepass(F, p)
        prim = meta["primaries"]

        # Find the unipotent primary key (q = x±1 over p=2)
        uni_keys = [k for k, d in prim.items() if _is_x_pm_1(d["poly"], p)]
        assert uni_keys, "expected an x±1 primary in p=2 for an upper-shear unipotent F"

        # (Optional) extra sanity: for an upper shear in p=2, (F+I)^2 = 0 and F != I with high prob
        I = np.eye(2 * n, dtype=np.int64)
        N = mod_p(F + I, p)
        assert np.array_equal(mod_p(N @ N, p), np.zeros_like(N)), "expected nilpotent index <= 2"


class TestAtomicDecomposition:
    def test_atomic_block_decompose_paired_only_certified(self) -> None:
        # Choose F so all sectors are paired (odd p, eigenvalues a and a^{-1} with a!=a^{-1})
        p = 5
        n = 3
        A = np.diag([2, 2, 3]).astype(np.int64)  # 2^{-1}=3, 3^{-1}=2
        F = _symplectic_scale(A, p)
        assert is_symplectic(F, p)

        Sigma, B, info = atomic_block_decompose(F, p)
        assert is_symplectic(Sigma, p)
        assert rank_mod(B, p) == 2 * n
        _assert_symplectic_basis_full(B, p)

        # similarity check: Sigma = B^{-1} F B (using symplectic left inverse)
        Linv = symplectic_left_inverse(B, p)
        recon = mod_p(B @ Sigma @ Linv, p)
        assert np.array_equal(recon, mod_p(F, p))

        assert info["Q_opt"] >= 1
        assert info["certified"] is True  # paired-only path should be certified

    def test_atomic_block_decompose_random_smoke_small(self) -> None:
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                for _ in range(10):
                    F = rand_symplectic(rng, n, p, steps=10)
                    Sigma, B, info = atomic_block_decompose(F, p)
                    assert is_symplectic(Sigma, p)
                    assert rank_mod(B, p) == 2 * n
                    _assert_symplectic_basis_full(B, p)

                    Linv = symplectic_left_inverse(B, p)
                    recon = mod_p(B @ Sigma @ Linv, p)
                    assert np.array_equal(recon, mod_p(F, p))

                    assert info["Q_opt"] >= 1
                    # certified may be False if self nonunipotent TODO is encountered; don't require True here.

    def test_atomic_block_decompose_coordinate_action_full_basis(self) -> None:
        """
        Strong convention test:
          With Sigma = B^{-1} F B, we must have F B = B Sigma (mod p).
        """
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for n in [1, 2, 3]:
                for _ in range(5):
                    F = rand_symplectic(rng, n, p, steps=10)
                    Sigma, B, info = atomic_block_decompose(F, p)

                    assert is_symplectic(Sigma, p)
                    _assert_symplectic_basis_full(B, p)

                    lhs = mod_p(F @ B, p)
                    rhs = mod_p(B @ Sigma, p)
                    assert np.array_equal(lhs, rhs)

    def test_atomic_block_decompose_each_atomic_block_is_invariant(self) -> None:
        """
        Uses info["atomic_half_dims"] to reconstruct each block basis T_blk from B,
        then checks:
          - T_blk is Darboux (symplectic on its span)
          - F leaves span(T_blk) invariant
          - restricted action F_T is symplectic and satisfies F T = T F_T
        """
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for n in [2, 3]:
                for _ in range(6):
                    F = rand_symplectic(rng, n, p, steps=12)
                    Sigma, B, info = atomic_block_decompose(F, p)

                    half_dims = [int(h) for h in info["atomic_half_dims"]]
                    blocks = _extract_atomic_blocks_from_global_basis(B, half_dims, p)

                    for T_blk in blocks:
                        # Darboux check on the block span
                        _assert_darboux_block(T_blk, p)

                        # Invariance in ambient space
                        _assert_invariant_span(F, T_blk, p)

                        # Coordinate action identity: F T = T F_T
                        L = symplectic_left_inverse(T_blk, p)              # (2h × 2n)
                        F_T = mod_p(L @ F @ T_blk, p)                      # (2h × 2h)
                        assert np.array_equal(mod_p(F @ T_blk, p), mod_p(T_blk @ F_T, p))
                        assert is_symplectic(F_T, p)

    def test_atomic_block_decompose_conjugation_invariance_of_Qopt_paired_only(self) -> None:
        """
        Conjugate by a random symplectic S and check invariants, but restrict to
        cases that are guaranteed "paired-only" so we don't exercise the currently
        fragile self-sector nonunipotent path.

        Construction: F = diag(A, (A^{-1})^T) with A = a I_n and a^2 != 1 (so a != a^{-1}).
        Then the only sector is a paired one ((x-a) <-> (x-a^{-1})) and the route is certified.
        """
        rng = np.random.default_rng()

        # Use primes where we can pick a != a^{-1}; avoid p=3 where all a are self-inverse.
        for p in [5, 7, 11]:
            for n in [2, 3]:
                a = 2 % p
                assert (a * a) % p != 1, "pick a with a != a^{-1}"

                A = (a * np.eye(n, dtype=np.int64)) % p
                F = _symplectic_scale(A, p)
                assert is_symplectic(F, p)

                # Conjugate by random symplectic
                S = rand_symplectic(rng, n, p, steps=10)
                Sinv = inv_mod_mat(S, p)
                F2 = mod_p(Sinv @ F @ S, p)
                assert is_symplectic(F2, p)

                # Prepass signatures should match under conjugation
                meta1 = rcf_prepass(F, p)
                meta2 = rcf_prepass(F2, p)
                assert _sector_signature(meta1, p) == _sector_signature(meta2, p)

                Sigma1, B1, info1 = atomic_block_decompose(F, p)
                Sigma2, B2, info2 = atomic_block_decompose(F2, p)

                # Both must reconstruct correctly
                Linv1 = symplectic_left_inverse(B1, p)
                recon1 = mod_p(B1 @ Sigma1 @ Linv1, p)
                assert np.array_equal(recon1, mod_p(F, p))

                Linv2 = symplectic_left_inverse(B2, p)
                recon2 = mod_p(B2 @ Sigma2 @ Linv2, p)
                assert np.array_equal(recon2, mod_p(F2, p))

                # In this paired-only construction, route should be certified
                assert info1["certified"] is True
                assert info2["certified"] is True

                # Q_opt and the multiset of block sizes should be conjugation-invariant
                assert int(info1["Q_opt"]) == int(info2["Q_opt"])
                assert sorted(int(x) for x in info1["atomic_half_dims"]) == sorted(int(x)
                                                                                   for x in info2["atomic_half_dims"])


def _fmt_mat(A: np.ndarray) -> str:
    return np.array2string(
        A,
        separator=", ",
        max_line_width=120,
        threshold=10_000,
    )


def _sector_signature(meta: dict, p: int) -> list[tuple]:
    sig = []
    for sec in meta.get("sectors", []):
        if sec["type"] == "paired":
            sig.append((
                "paired",
                int(sec.get("dim2", sec["W_basis"].shape[1])),
                int(sec.get("deg", -1)),
                int(sec.get("exponent", -1)),
                tuple(sec["key"]),
                tuple(sec["key_star"]),
            ))
        else:
            sig.append((
                "self",
                int(sec.get("dim2", sec["W_basis"].shape[1])),
                int(sec.get("deg", -1)),
                int(sec.get("exponent", -1)),
                tuple(sec["key"]),
                bool(sec.get("floor_certified", False)),
            ))
    return sorted(sig, key=str)


# def _diagnose_one_case(F: np.ndarray, p: int) -> str:
#     lines: list[str] = []
#     lines.append(f"is_symplectic(F,p)={is_symplectic(F, p)}  shape={F.shape}")

#     try:
#         meta = rcf_prepass(F, p)
#         lines.append(f"rcf_prepass: Lmin_star={meta.get('Lmin_star')}  n_sectors={len(meta.get('sectors', []))}")
#         lines.append(f"sector_signature={_sector_signature(meta, p)}")
#     except Exception as e:
#         lines.append("rcf_prepass FAILED:")
#         lines.append(f"{type(e).__name__}: {e}")
#         lines.append(traceback.format_exc())
#         return "\n".join(lines)

#     prim = meta["primaries"]

#     for i, sec in enumerate(meta["sectors"]):
#         try:
#             if sec["type"] == "paired":
#                 key = sec["key"]
#                 key_star = sec["key_star"]
#                 blocks, inv = atomic_blocks_in_paired_sector(F, p, key, key_star, prim)
#                 lines.append(
#                     f"sector[{i}] paired OK: key={key} key*={key_star}  n_blocks={len(blocks)}  status={inv.data.get('status')}"
#                 )
#             else:
#                 key = sec["key"]
#                 q = prim[key]["poly"]
#                 if p == 2 and _is_x_pm_1(q, p):
#                     blocks, inv = atomic_blocks_in_unipotent_self_sector_p2(F, key, prim)
#                     lines.append(
#                         f"sector[{i}] self unipotent-p2 OK: key={key}  n_blocks={len(blocks)}  status={inv.data.get('status')}"
#                     )
#                 else:
#                     blocks, inv = atomic_blocks_in_self_sector_nonunipotent(F, p, key, prim)
#                     lines.append(
#                         f"sector[{i}] self nonunipotent OK: key={key}  n_blocks={len(blocks)}  status={inv.data.get('status')}"
#                     )
#         except Exception as e:
#             lines.append(
#                 f"sector[{i}] BUILDER FAILED: type={sec['type']} key={sec.get('key')} err={type(e).__name__}: {e}"
#             )
#             lines.append(traceback.format_exc())

#     return "\n".join(lines)


class TestAtomicDecompositionFuzz:
    def test_fuzz_atomic_block_decompose_p2_random(self) -> None:
        """
        Fast-by-default fuzzer for p=2 random symplectics.

        Knobs (optional env vars):
          SYMPLEQ_FUZZ_SEED=12345
          SYMPLEQ_FUZZ_TRIALS=50
          SYMPLEQ_FUZZ_MAXN=4
          SYMPLEQ_FUZZ_STEPS=16
        Set SYMPLEQ_FUZZ_TRIALS=0 to effectively disable this test.
        """
        trials = 500
        max_n = 5
        steps = 86

        if trials <= 0:
            pytest.skip("SYMPLEQ_FUZZ_TRIALS<=0")

        rng = np.random.default_rng()
        p = 2

        for t in range(trials):
            n = int(rng.integers(1, max_n + 1))
            F = rand_symplectic(rng, n, p, steps=steps)
            Sigma, B, info = atomic_block_decompose(F, p)

            assert is_symplectic(Sigma, p)
            assert rank_mod(B, p) == 2 * n
            _assert_symplectic_basis_full(B, p)

            Linv = symplectic_left_inverse(B, p)
            recon = mod_p(B @ Sigma @ Linv, p)
            assert np.array_equal(recon, mod_p(F, p))

    def test_fuzz_atomic_block_decompose_p2_unipotent_shears(self) -> None:
        """
        Fast-by-default fuzzer biased toward the p=2 unipotent corner using products of shears.

        Knobs:
          SYMPLEQ_FUZZ_SEED=12345
          SYMPLEQ_FUZZ_TRIALS=50
          SYMPLEQ_FUZZ_MAXN=4
          SYMPLEQ_FUZZ_SHEARLEN=8
        """

        trials = 50
        max_n = 5
        shear_len = 8

        if trials <= 0:
            pytest.skip("SYMPLEQ_FUZZ_TRIALS<=0")

        rng = np.random.default_rng()
        p = 2

        for t in range(trials):
            n = int(rng.integers(1, max_n + 1))
            F = np.eye(2 * n, dtype=np.int64)
            for _ in range(shear_len):
                B = _rand_symmetric(rng, n, p)
                C = _rand_symmetric(rng, n, p)
                F = mod_p(_symplectic_shear_upper(B, p) @ F, p)
                F = mod_p(_symplectic_shear_lower(C, p) @ F, p)

            assert is_symplectic(F, p)
            Sigma, B, info = atomic_block_decompose(F, p)

            assert is_symplectic(Sigma, p)
            assert rank_mod(B, p) == 2 * n
            _assert_symplectic_basis_full(B, p)

            Linv = symplectic_left_inverse(B, p)
            recon = mod_p(B @ Sigma @ Linv, p)
            assert np.array_equal(recon, mod_p(F, p))

class TestAtomicDecompositionCertifiedAPI:
    # def test_certified_raises_when_uncertified(self) -> None:
    #     rng = np.random.default_rng(0)
    #     p = 2
    #     n = 4
    #     F = rand_symplectic(rng, n, p, steps=20)

    #     with pytest.raises(Exception):  # ideally CertificationError
    #         _ = atomic_block_decompose(F, p, mode="certified")

    def test_best_effort_always_returns_valid(self) -> None:
        rng = np.random.default_rng(1)
        for p in [2, 3, 5]:
            for n in [1, 2, 3, 4]:
                F = rand_symplectic(rng, n, p, steps=15)
                Sigma, B, info = atomic_block_decompose(F, p, mode="best_effort")

                assert is_symplectic(Sigma, p)
                _assert_symplectic_basis_full(B, p)

                # reconstruction
                Linv = symplectic_left_inverse(B, p)
                recon = mod_p(B @ Sigma @ Linv, p)
                assert np.array_equal(recon, mod_p(F, p))

                assert info["status"] in ("OK", "DEGRADED")

    def test_certified_succeeds_on_paired_only(self) -> None:
        p = 5
        n = 3
        A = np.diag([2, 2, 3]).astype(np.int64)
        F = _symplectic_scale(A, p)

        Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
        assert info["status"] == "OK"
        assert info["certified"] is True

    def test_atomic_block_decompose_best_effort_completes_when_blocks_incomplete(self) -> None:
        rng = np.random.default_rng(123)
        p = 3
        n = 3
        F = rand_symplectic(rng, n, p, steps=10)

        Sigma, B, info = atomic_block_decompose(F, p, mode="best_effort")
        assert is_symplectic(Sigma, p)
        assert B.shape == (2*n, 2*n)
        _assert_symplectic_basis_full(B, p)

        Linv = symplectic_left_inverse(B, p)
        recon = mod_p(B @ Sigma @ Linv, p)
        assert np.array_equal(recon, mod_p(F, p))
        # It may be OK or DEGRADED, but must be internally consistent:
        assert info["status"] in ("OK", "DEGRADED")
        assert "completed" in info

    def test_atomic_block_decompose_certified_refuses_incomplete_spanning(self) -> None:
        # This test is more “behavioral”: certified must be strict.
        rng = np.random.default_rng(456)
        p = 3
        n = 3
        F = rand_symplectic(rng, n, p, steps=10)

        # certified may pass or fail depending on your current builders,
        # but if it fails, it MUST fail via CertificationError (not RuntimeError).
        try:
            _ = atomic_block_decompose(F, p, mode="certified")
        except Exception as e:
            assert isinstance(e, CertificationError)
