"""
Focused tests for the p=2 self-reciprocal non-unipotent Hermitian sector.

These are intended to exercise q(x)=x^2+x+1 over GF(2), which is the smallest
non-unipotent self-reciprocal irreducible factor.
"""
from __future__ import annotations

import numpy as np
import pytest

from sympleq.core.minimal_qudit_frame import minimal_qudit_frame
from sympleq.core.minimal_qudit_frame.helpers.mqf_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import mod_p, inv_mod_mat, is_symplectic
from sympleq.core.circuits.random_symplectic import symplectic_random_transvection
from sympleq.core.minimal_qudit_frame import MinimalQuditFrameCertificationError


def companion_x2_x_1() -> np.ndarray:
    # Column-action companion for x^2+x+1: [[0,1],[1,1]].
    return np.array([[0, 1], [1, 1]], dtype=np.int64)


def symplectic_scale(A: np.ndarray, p: int = 2) -> np.ndarray:
    A = mod_p(A, p)
    Ainv = inv_mod_mat(A, p)
    Z = np.zeros_like(A)
    return mod_p(np.block([[A, Z], [Z, Ainv.T]]), p)


def direct_sum_grouped(F1: np.ndarray, F2: np.ndarray, p: int = 2) -> np.ndarray:
    n1 = F1.shape[0] // 2
    n2 = F2.shape[0] // 2
    A1, B1, C1, D1 = F1[:n1, :n1], F1[:n1, n1:], F1[n1:, :n1], F1[n1:, n1:]
    A2, B2, C2, D2 = F2[:n2, :n2], F2[:n2, n2:], F2[n2:, :n2], F2[n2:, n2:]
    Z12 = np.zeros((n1, n2), dtype=np.int64)
    Z21 = np.zeros((n2, n1), dtype=np.int64)
    A = np.block([[A1, Z12], [Z21, A2]])
    B = np.block([[B1, Z12], [Z21, B2]])
    C = np.block([[C1, Z12], [Z21, C2]])
    D = np.block([[D1, Z12], [Z21, D2]])
    return mod_p(np.block([[A, B], [C, D]]), p)


def test_p2_x2_x_1_single_sector_certifies() -> None:
    p = 2
    A = companion_x2_x_1()
    F = symplectic_scale(A, p)
    assert is_symplectic(F, p)

    Sigma, B, info = minimal_qudit_frame(F, p)
    verify_global_basis(F, B, Sigma, p)

    assert info["certified"] is True
    assert info["minimal_cost_certified"] is True
    assert sorted(int(h) for h in info["mqf_half_dims"]) == [1, 1]
    assert info["Q_opt"] == 1


def test_p2_x2_x_1_repeated_conjugated_certifies() -> None:
    p = 2
    A = companion_x2_x_1()
    F1 = symplectic_scale(A, p)
    F = direct_sum_grouped(F1, F1, p)
    assert is_symplectic(F, p)

    rng = np.random.default_rng(0)
    S = symplectic_random_transvection(F.shape[0] // 2, p, num_transvections=40, rng=rng)
    Fh = mod_p(inv_mod_mat(S, p) @ F @ S, p)
    assert is_symplectic(Fh, p)

    Sigma, B, info = minimal_qudit_frame(Fh, p)
    verify_global_basis(Fh, B, Sigma, p)

    assert info["certified"] is True
    assert info["minimal_cost_certified"] is True
    assert sorted(int(h) for h in info["mqf_half_dims"]) == [1, 1, 1, 1]
    assert info["Q_opt"] == 1
