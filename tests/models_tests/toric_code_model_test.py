import numpy as np
import pytest

from sympleq.models.toric_code import ToricCode


def test_toric_code_accepts_edge_resolved_gauge_coefficients():
    c_g = np.arange(1, 9, dtype=float)
    tc = ToricCode(Nx=2, Ny=2, c_x=0.0, c_z=0.0, c_g=c_g, periodic=True, d=2)

    terms, coeffs = tc.build_toric_code_hamiltonian()

    assert len(terms) == tc.n_qubits
    np.testing.assert_allclose(coeffs, c_g)


def test_toric_code_rejects_wrong_length_edge_resolved_gauge_coefficients():
    with pytest.raises(ValueError, match="Vector c_g must have length"):
        ToricCode(Nx=2, Ny=2, c_x=0.0, c_z=0.0, c_g=[1.0, 2.0], periodic=True, d=2)
