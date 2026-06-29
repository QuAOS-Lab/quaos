import numpy as np

from sympleq.core.paulis._typing import HilbertOperator


def fidelity(rho: np.ndarray | HilbertOperator, sigma: np.ndarray | HilbertOperator) -> float:
    def _sparse_fidelity(rho: HilbertOperator, sigma: HilbertOperator) -> float:
        from scipy.sparse.linalg import eigsh
        lam, V = eigsh(rho)
        sqrt_lam = np.sqrt(np.maximum(lam, 0))
        inner = (V * sqrt_lam) .T.conj() @ sigma @ (V * sqrt_lam)
        # Eigendecompose M to get sqrt(M)
        mu = np.linalg.eigvalsh(inner)
        return float(np.real(np.sum(np.sqrt(np.maximum(mu, 0)))) ** 2)

    if isinstance(rho, HilbertOperator) and isinstance(sigma, HilbertOperator):
        return _sparse_fidelity(rho, sigma)

    from scipy.linalg import sqrtm
    sq_rho = sqrtm(rho)
    tmp_matrix = sq_rho @ sigma @ sq_rho
    return float(np.real(np.trace(sqrtm(tmp_matrix)) ** 2))
