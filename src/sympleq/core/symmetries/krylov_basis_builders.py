import numpy as np
from .modular_helpers import mod_p, _solve_linear


def build_partner_in_span(K: np.ndarray, Omega: np.ndarray, span_basis: np.ndarray, p: int) -> np.ndarray:
    """
    Solve for z in span(span_basis) such that <K[:,b], z> = δ_{b, k-1}.
    """
    k = K.shape[1]
    # rows_b = (K[:,b])^T Omega span_basis  -> shape (k, d)
    A = mod_p(K.T @ Omega @ span_basis, p)
    b = np.zeros((k, 1), dtype=np.int64)
    b[-1, 0] = 1
    coeff = _solve_linear(A, b, p)     # (d,1)
    z = mod_p(span_basis @ coeff, p)   # (n2,1)
    return z
