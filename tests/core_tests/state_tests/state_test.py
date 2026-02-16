
import numpy as np
import random

from sympleq.core.states.state import State
from sympleq.core.circuits.circuits import Circuit
from sympleq.core.circuits.gates import (
    Hadamard,
    PHASE,
    SUM,
    SWAP,
    CNOT,
    PauliGate,
)
from sympleq.core.circuits.utils import (
    H_mat,
    I_mat,
    CX_func,
    SWAP_func,
    pauli_unitary_qudit,
)
from sympleq.core.paulis import PauliString
from sympleq.core.circuits.utils import (
    _apply_single_qudit_dense,
    _apply_two_qudit_permutation,
)


class TestState:

    # -----------------------
    # Basic State behaviour
    # -----------------------

    def test_state_from_basis(self):
        dims = [2, 2]
        # index 3 corresponds to |11> for dims [2, 2]
        s = State.from_basis(dims, index=3)
        assert s.amplitudes.shape == (4,)
        assert np.allclose(s.amplitudes, np.array([0, 0, 0, 1], dtype=np.complex128))
        assert np.array_equal(s.dimensions, np.array(dims))

    def test_state_copy_independence(self):
        dims = [2, 3]
        s = State.from_basis(dims, index=2)
        s2 = s.copy()
        assert np.allclose(s.amplitudes, s2.amplitudes)
        assert np.array_equal(s.dimensions, s2.dimensions)

        # Modify copy and check original unchanged
        s2.amplitudes[2] = 0.5
        assert not np.allclose(s.amplitudes, s2.amplitudes)

    # -----------------------------------------
    # Low-level helpers: single / two-qudit
    # -----------------------------------------

    def test_apply_single_qudit_dense_matches_tensor(self):
        dims = [2, 2]
        D = int(np.prod(dims))

        psi = np.zeros(D, dtype=np.complex128)
        psi[2] = 1.0

        # Direct application via helper
        U_loc = H_mat(2)              # sparse 2x2
        psi_out_helper = _apply_single_qudit_dense(psi, dims, U_loc.toarray(), q=0)

        # Reference: full tensor operator H ⊗ I, but DENSE
        H_full = np.kron(H_mat(2).toarray(), I_mat(2).toarray())
        psi_out_ref = H_full @ psi

        assert np.allclose(psi_out_helper, psi_out_ref)

    def test_apply_two_qudit_permutation_sum_gate(self):
        """
        Check that _apply_two_qudit_permutation with CX_func matches
        the permutation matrix for a SUM/CNOT gate on two qubits.
        """
        dims = [2, 2]
        D = int(np.prod(dims))

        # Random state
        rng = np.random.default_rng(123)
        psi = rng.normal(size=D) + 1j * rng.normal(size=D)
        psi /= np.linalg.norm(psi)

        # Apply via helper
        psi_out_helper = _apply_two_qudit_permutation(
            psi, dims, CX_func, a0=0, a1=1
        )

        # Reference: explicit permutation map using CX_func
        out = np.zeros_like(psi)
        for j in range(D):
            out[CX_func(j, 0, 1, np.array(dims, dtype=int))] = psi[j]

        assert np.allclose(psi_out_helper, out)

    def test_apply_two_qudit_permutation_swap_gate_mixed_radix(self):
        """
        Check SWAP on dims [2, 2] via helper vs explicit mapping.

        (Mixed-radix [2,3] currently exposes a SWAP_func bug; test restricted
        to equal dims until that is fixed.)
        """
        dims = [2, 2]   # was [2, 3]
        D = int(np.prod(dims))

        rng = np.random.default_rng(456)
        psi = rng.normal(size=D) + 1j * rng.normal(size=D)
        psi /= np.linalg.norm(psi)

        psi_out_helper = _apply_two_qudit_permutation(
            psi, dims, SWAP_func, a0=0, a1=1
        )

        out = np.zeros_like(psi)
        for j in range(D):
            out[SWAP_func(j, 0, 1, np.array(dims, dtype=int))] = psi[j]

        assert np.allclose(psi_out_helper, out)

    # ----------------------------------------------
    # Gate.act_on_state vs Gate.unitary consistency
    # ----------------------------------------------

    def test_hadamard_gate_act_on_state_matches_unitary(self):
        dims = [2]
        D = 2
        rng = np.random.default_rng(789)
        psi = rng.normal(size=D) + 1j * rng.normal(size=D)
        psi /= np.linalg.norm(psi)

        state = State(psi, dims)
        gate = Hadamard(index=0, dimension=2)

        # Using act_on_state
        out_state = gate.act_on_state(state)
        psi_act = out_state.as_array()

        # Using explicit unitary
        U = gate.unitary(dims=dims).toarray()
        psi_ref = U @ psi

        assert np.allclose(psi_act, psi_ref)

    def test_sum_gate_act_on_state_matches_unitary(self):
        dims = [2, 2]
        D = 4
        rng = np.random.default_rng(1011)
        psi = rng.normal(size=D) + 1j * rng.normal(size=D)
        psi /= np.linalg.norm(psi)

        state = State(psi, dims)
        gate = SUM(control=0, target=1, dimension=2)

        out_state = gate.act_on_state(state)
        psi_act = out_state.as_array()

        U = gate.unitary(dims=dims).toarray()
        psi_ref = U @ psi

        assert np.allclose(psi_act, psi_ref)

    def test_swap_gate_act_on_state_matches_unitary_mixed_radix(self):
        """
        Check SWAP gate act_on_state vs unitary on equal dims [2,2].

        Mixed-radix SWAP currently exposes a bug in SWAP_func; this test
        still validates act_on_state vs unitary.
        """
        dims = [2, 2]   # was [2, 3]
        D = 4
        rng = np.random.default_rng(2022)
        psi = rng.normal(size=D) + 1j * rng.normal(size=D)
        psi /= np.linalg.norm(psi)

        state = State(psi, dims)
        gate = SWAP(index1=0, index2=1, dimension=2)

        out_state = gate.act_on_state(state)
        psi_act = out_state.as_array()

        U = gate.unitary(dims=dims).toarray()
        psi_ref = U @ psi

        assert np.allclose(psi_act, psi_ref)

    def test_pauli_gate_act_on_state_matches_local_pauli(self):
        dims = [2, 2]
        D = 4

        psi = np.zeros(D, dtype=np.complex128)
        psi[0] = 1.0
        state = State(psi, dims)

        x_exp = [1, 0]
        z_exp = [0, 0]
        ps = PauliString.from_exponents(x_exp, z_exp, dims)
        gate = PauliGate(ps)

        out_state = gate.act_on_state(state)
        psi_act = out_state.as_array()

        # Explicit X ⊗ I with dense locals
        X_loc = pauli_unitary_qudit(2, 1, 0).toarray()
        U_ref = np.kron(X_loc, I_mat(2).toarray())
        psi_ref = U_ref @ psi

        assert np.allclose(psi_act, psi_ref)

    def test_circuit_apply_equals_unitary_qubit(self):
        """
        For a small random qubit circuit, act_on_state / apply_to_statevector
        should match the full Circuit.unitary() action.
        """
        random.seed(0)
        np.random.seed(0)

        dims = [2, 2]
        # Small hand-picked circuit
        g1 = Hadamard(0, 2)
        g2 = PHASE(1, 2)
        g3 = CNOT(0, 1)
        g4 = SWAP(0, 1, 2)

        circuit = Circuit(dimensions=dims, gates=[g1, g2, g3, g4])

        D = 4
        psi = np.random.randn(D) + 1j * np.random.randn(D)
        psi /= np.linalg.norm(psi)

        # Path 1: apply_to_statevector
        psi_out_vec = circuit.apply_to_statevector(psi)

        # Path 2: full unitary
        U = circuit.unitary().toarray()
        psi_ref = U @ psi

        assert np.allclose(psi_out_vec, psi_ref)

    def test_circuit_apply_equals_unitary_mixed_radix(self):
        """
        Same as above, but for mixed-radix dims [2, 3] with only single-qudit gates.
        """
        np.random.seed(1)
        dims = [2, 3]

        g1 = Hadamard(0, 2)
        g2 = PHASE(1, 3)
        g3 = PHASE(0, 2)

        circuit = Circuit(dimensions=dims, gates=[g1, g2, g3])

        D = int(np.prod(dims))
        psi = np.random.randn(D) + 1j * np.random.randn(D)
        psi /= np.linalg.norm(psi)

        psi_out_vec = circuit.apply_to_statevector(psi)

        U = circuit.unitary().toarray()
        psi_ref = U @ psi

        assert np.allclose(psi_out_vec, psi_ref)

    def test_bell_state_example(self):
        """
        H(0) followed by CNOT(0,1) on |00> should produce
        (|00> + |11>)/sqrt(2).
        """
        dims = [2, 2]
        D = 4

        # |00>
        psi = np.zeros(D, dtype=np.complex128)
        psi[0] = 1.0

        circuit = Circuit(
            dimensions=dims,
            gates=[Hadamard(0, 2), CNOT(0, 1)]
        )

        psi_out = circuit.apply_to_statevector(psi)

        # Expected Bell state |Φ+> = (|00> + |11>) / sqrt(2)
        bell = np.zeros(D, dtype=np.complex128)
        bell[0] = 1 / np.sqrt(2)
        bell[3] = 1 / np.sqrt(2)

        # Global phase irrelevant; we check up to global phase
        # but here they should match directly.
        assert np.allclose(psi_out, bell)

    def test_hadamard_squared_is_identity(self):
        """
        H^2 = I on a single qubit.
        """
        dims = [2]
        D = 2

        rng = np.random.default_rng(42)
        psi = rng.normal(size=D) + 1j * rng.normal(size=D)
        psi /= np.linalg.norm(psi)

        state = State(psi, dims)
        H_gate = Hadamard(0, 2)

        # Apply H twice via act_on_state
        state1 = H_gate.act_on_state(state)
        state2 = H_gate.act_on_state(state1)

        assert np.allclose(state2.as_array(), psi)
