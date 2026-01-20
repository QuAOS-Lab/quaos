from __future__ import annotations
import numpy as np
from .modular_helpers import mod_p, omega_matrix, inv_mod_mat, nullspace_mod, independent_columns


def symplectic_left_inverse(T: np.ndarray, p: int) -> np.ndarray:
    """
    For T with T^T Ω T invertible (typically Ω_k), return L such that L T = I on span(T).
    L = (T^T Ω T)^{-1} T^T Ω
    """
    n2 = T.shape[0]
    Ω = omega_matrix(n2 // 2, p)
    J = mod_p(T.T @ Ω @ T, p)
    Jinv = inv_mod_mat(J, p)
    return mod_p(Jinv @ T.T @ Ω, p)


def restrict_operator(F: np.ndarray, T: np.ndarray, p: int) -> np.ndarray:
    """
    Return F_restricted in the T-coordinates:  F_T = T^{-1} F T
    where T^{-1} means symplectic left inverse above.
    """
    L = symplectic_left_inverse(T, p)
    return mod_p(L @ F @ T, p)


def symplectic_orthogonal_complement_in_ambient(T: np.ndarray, p: int) -> np.ndarray:
    """
    Return a basis (columns) for W^⊥ in the ambient space, where W = span(T).
    W^⊥ = {x : T^T Ω x = 0}.
    """
    n2 = T.shape[0]
    Ω = omega_matrix(n2 // 2, p)
    A = mod_p(T.T @ Ω, p)
    N = nullspace_mod(A, p)
    return independent_columns(N, p)
