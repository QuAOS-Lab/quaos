from __future__ import annotations

import itertools

import numpy as np
import pytest

from sympleq.core.minimal_qudit_frame.helpers.mqf_extension import (
    GFpExtension,
    field_scalar_apply_to_vector,
)
from sympleq.core.symmetries.modular_helpers import mod_p, rank_mod


def _companion_multiplication_matrix(q: list[int], p: int) -> np.ndarray:
    """Matrix of multiplication by alpha on GF(p)[alpha]/q(alpha)."""
    q = [int(c) % p for c in q]
    d = len(q) - 1
    C = np.zeros((d, d), dtype=np.int64)
    for j in range(d - 1):
        C[j + 1, j] = 1
    for i in range(d):
        C[i, d - 1] = (-q[i]) % p
    return mod_p(C, p)


@pytest.mark.parametrize(
    "p,q",
    [
        (2, [1, 1, 1]),      # x^2+x+1
        (3, [1, 0, 1]),      # x^2+1
        (5, [1, 1, 1]),      # x^2+x+1
        (2, [1, 1, 0, 1]),   # x^3+x+1
    ],
)
def test_gfp_extension_field_axioms_and_trace_duality(p: int, q: list[int]) -> None:
    field = GFpExtension(q, p)
    elems = field.elements()
    nonzero = field.nonzero_elements()

    assert field.reduce_coeffs(q) == field.zero

    for a in elems:
        assert field.add(a, field.zero) == a
        assert field.mul(a, field.one) == a
        assert field.sub(a, a) == field.zero

    for a, b, c in itertools.product(elems, repeat=3):
        assert field.add(field.add(a, b), c) == field.add(a, field.add(b, c))
        assert field.mul(field.mul(a, b), c) == field.mul(a, field.mul(b, c))
        assert field.mul(a, field.add(b, c)) == field.add(field.mul(a, b), field.mul(a, c))

    for a in nonzero:
        assert field.mul(a, field.inv(a)) == field.one

    # Trace pairing must be nondegenerate, and from_trace_values must invert it.
    Tinv = field.trace_pairing_inverse()
    assert rank_mod(Tinv, p) == field.d
    for h in elems:
        values = [field.trace_to_base(field.mul(field.alpha_power(a), h)) for a in range(field.d)]
        assert field.from_trace_values(values) == h


@pytest.mark.parametrize(
    "p,q",
    [
        (2, [1, 1, 1]),
        (3, [1, 0, 1]),
        (5, [1, 1, 1]),
    ],
)
def test_reciprocal_conjugation_is_an_involution_for_self_reciprocal_q(p: int, q: list[int]) -> None:
    field = GFpExtension(q, p)
    assert field.conj(field.alpha) == field.alpha_inv()
    for a in field.elements():
        assert field.conj(field.conj(a)) == a
        assert field.conj(field.add(a, field.one)) == field.add(field.conj(a), field.one)
        assert field.conj(field.mul(a, field.alpha)) == field.mul(field.conj(a), field.alpha_inv())


@pytest.mark.parametrize(
    "p,q",
    [
        (2, [1, 1, 1]),
        (3, [1, 0, 1]),
        (5, [1, 1, 1]),
    ],
)
def test_field_scalar_apply_matches_polynomial_action(p: int, q: list[int]) -> None:
    field = GFpExtension(q, p)
    F = _companion_multiplication_matrix(q, p)

    for scalar in field.elements():
        for v_tuple in field.elements():
            v = field.coeff_vector(v_tuple)
            lhs = field_scalar_apply_to_vector(F, field, scalar, v, p)
            rhs = field.coeff_vector(field.mul(scalar, v_tuple))
            assert np.array_equal(lhs % p, rhs % p)
