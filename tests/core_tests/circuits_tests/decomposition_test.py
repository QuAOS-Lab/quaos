import numpy as np

from sympleq.core.circuits import Circuit, GATES
from sympleq.core.circuits.utils import is_symplectic
from sympleq.core.circuits.circuits import GateSpec
from sympleq.core.circuits.gate_decomposition_to_circuit import (inv_gfp, mod_p, ensure_invertible_A_circuit, blocks,
                                                                 zeros, identity, _is_invertible_mod,
                                                                 synth_linear_A_to_gates, _emit_local_ops_for_D,
                                                                 synth_lower_from_symmetric,
                                                                 synth_upper_from_symmetric_via_H, _as_int_mod,
                                                                 decompose_symplectic_to_circuit,
                                                                 _compose_symp, _full_from_lower, gate_to_circuit)

from sympleq.core.paulis import PauliSum

# Convenience aliases for singleton gates
CX = GATES.CX
S = GATES.S


class TestDecomposition:
    def _full_from_upper(self, n, S, p):
        Id = np.eye(n, dtype=int)
        Z = np.zeros((n, n), dtype=int)
        return mod_p(np.block([[Id, S], [Z, Id]]), p)

    def random_invertible(self, n: int, rng: np.random.Generator, p: int) -> np.ndarray:
        while True:
            A = rng.integers(0, p, size=(n, n), dtype=int)
            try:
                _ = inv_gfp(A, p)
                return mod_p(A, p)
            except ValueError:
                pass

    def random_symmetric(self, n: int, rng: np.random.Generator, p: int) -> np.ndarray:
        M = rng.integers(0, p, size=(n, n), dtype=int)
        return mod_p(M + M.T, p)

    def build_F_from_LMR(self, n: int, A: np.ndarray, B_sym: np.ndarray, C_sym: np.ndarray, p: int) -> np.ndarray:
        """F = L M R with R=[I B;0 I], M=diag(A, (A^T)^{-1}), L=[I 0; C I]."""
        A = mod_p(A, p)
        B_sym = mod_p(B_sym, p)
        C_sym = mod_p(C_sym, p)
        A_invT = inv_gfp(A.T, p)
        Z = zeros((n, n))
        Id = identity(n)
        L = np.block([[Id, Z], [C_sym, Id]])
        M = np.block([[A, Z], [Z, A_invT]])
        R = np.block([[Id, B_sym], [Z, Id]])
        return mod_p(L @ M @ R, p)

    def random_symplectic(self, n: int, p: int, rng: np.random.Generator | None = None) -> np.ndarray:
        if rng is None:
            rng = np.random.default_rng()
        A = self.random_invertible(n, rng, p)
        B = self.random_symmetric(n, rng, p)
        C = self.random_symmetric(n, rng, p)
        F = self.build_F_from_LMR(n, A, B, C, p)
        assert is_symplectic(F, p)
        return F

    # ---------- preconditioner tests----------

    def preconditioner_noop_when_A_invertible(self, F: np.ndarray, p: int):
        n = F.shape[0] // 2
        from copy import deepcopy
        F0 = deepcopy(F)
        C_pre = ensure_invertible_A_circuit(F, p)
        if np.linalg.det(F[:n, :n]) % p != 0:  # or inv_gfp try/except
            # already invertible: either empty or something that keeps A invertible
            assert len(C_pre.gates) == 0
            F_prime = (C_pre.full_symplectic() @ F) % p
            assert np.array_equal(F_prime, F0 % p)

    def preconditioner_preserves_symplectic(self, F: np.ndarray, p: int):
        C_pre = ensure_invertible_A_circuit(F, p)
        F_prime = (C_pre.full_symplectic() @ F) % p
        assert is_symplectic(F_prime, p)

    def A_invertible_after_precondition(self, F: np.ndarray, p: int):
        n = F.shape[0] // 2
        C_pre = ensure_invertible_A_circuit(F, p)
        F_prime = (C_pre.full_symplectic() @ F) % p
        A = F_prime[:n, :n]
        _ = inv_gfp(A, p)  # should not raise

    def preconditioner_determinism(self, F: np.ndarray, p: int):
        C1 = ensure_invertible_A_circuit(F, p)
        C2 = ensure_invertible_A_circuit(F, p)
        assert [g.name for g in C1.gates] == [g.name for g in C2.gates]
        assert C1.qudit_indices == C2.qudit_indices

    def preconditioner_respects_depth(self, F: np.ndarray, p: int):
        C_pre = ensure_invertible_A_circuit(F, p, max_depth=4)
        assert len(C_pre.gates) <= 4

    def test_preconditioner_after_random_conjugation(self):
        num_trials = 10
        depth = 10
        n = 3
        rng = np.random.default_rng()

        for p in [2, 3, 5]:
            for _ in range(num_trials):
                F = self.random_symplectic(n, p, rng)
                # random circuit (uses only valid gates for dimension p)
                C = Circuit.from_random(n_gates=depth, dimensions=[p] * n)
                F_conj = (C.full_symplectic() @ F) % p
                assert is_symplectic(F_conj, p), "Conjugated F must remain symplectic"
                C_pre = ensure_invertible_A_circuit(F_conj, p)
                F_prime = (C_pre.full_symplectic() @ F_conj) % p
                A = F_prime[:n, :n]
                A_inv_ok = False
                try:
                    _ = inv_gfp(A, p)
                    A_inv_ok = True
                except ValueError:
                    A_inv_ok = False

                # --- assertions ---
                assert A_inv_ok, "A-block is not invertible after preconditioning"
                assert is_symplectic(F_prime, p), "Preconditioning broke symplecticity"

    def test_preconditioner(self) -> None:
        """
        Generate random symplectic F, build C_pre, and verify that
        A((C_pre)·F) is invertible (mod p). Also sanity-check symplecticity.
        """
        num_trials = 10

        rng = np.random.default_rng()

        n = 3
        for p in [2, 3, 5]:

            for t in range(num_trials):
                F = self.random_symplectic(n, p, rng)
                C_pre = ensure_invertible_A_circuit(F, p)
                F_prime = (C_pre.full_symplectic() @ F) % p  # left multiply for the matrix factorization step
                A, _, _, _ = blocks(F_prime, p)
                okA = _is_invertible_mod(A, p)
                okSp = is_symplectic(F_prime, p)
                if not okA or not okSp:
                    raise AssertionError(
                        f"[trial {t}] Preconditioner failed:\n"
                        f"  A invertible? {okA}\n"
                        f"  F' symplectic? {okSp}\n"
                        f"F:\n{F}\nF':\n{F_prime}\nA(F'):\n{A}"
                    )
                self.preconditioner_noop_when_A_invertible(F, p)
                self.A_invertible_after_precondition(F, p)
                self.preconditioner_preserves_symplectic(F, p)
                self.preconditioner_determinism(F, p)
                self.preconditioner_respects_depth(F, p)

    # M block tests
    def symp_of_gates(self, n, p, ops: list[GateSpec]):
        """Compose a list[GateSpec] into full symplectic matrix."""
        gates = [op[0] for op in ops]
        qudits = [op[1:] for op in ops]
        return Circuit.from_gates_and_qudits([p] * n, gates, qudits).composite_gate().symplectic % p

    def check_linear_from_A(self, n, p, A):
        A = mod_p(A, p)
        D = mod_p(inv_gfp(A, p).T, p)
        Z = np.zeros((n, n), dtype=int)
        return mod_p(np.block([[A, Z], [Z, D]]), p)

    def M_block_only(self, n, p, rng=np.random.default_rng(42)):
        A = self.random_invertible(n, rng, p)
        ops_M = synth_linear_A_to_gates(n, A, p)
        FM = self.symp_of_gates(n, p, ops_M)
        FM_exp = self.check_linear_from_A(n, p, A)
        assert np.array_equal(FM, FM_exp), f"\nA=\n{A}\nGot=\n{FM}\nExp=\n{FM_exp}"

    def local_Du_synthesis(self, p):
        for u in range(1, p):
            ops = _emit_local_ops_for_D(0, u, p)
            gates = [op[0] for op in ops]
            qudits = [op[1:] for op in ops]
            F = Circuit.from_gates_and_qudits([p], gates, qudits).composite_gate().symplectic % p
            target = np.array([[u, 0],
                               [0, pow(int(u), -1, p)]], dtype=int) % p
            assert np.array_equal(F, target), f"Failed D({u}) over GF({p})"

    def CX_direction_in_M(self, n: int, p: int, rng: np.random.Generator = np.random.default_rng(42)):
        # Try a few random (i, c, f) with i != c
        for _ in range(20):
            i = int(rng.integers(0, n))
            c = int(rng.integers(0, n))
            while i == c:
                c = int(rng.integers(0, n))
            f = int(rng.integers(1, p))   # 1..p-1

            A_exp = np.eye(n, dtype=int) % p
            A_exp[c, i] = (A_exp[c, i] + (-f) % p) % p  # NOTE: (c, i), not (i, c)
            ops = synth_linear_A_to_gates(n, A_exp, p)
            FM = self.symp_of_gates(n, p, ops)
            FM_exp = self.check_linear_from_A(n, p, A_exp)

            assert np.array_equal(FM, FM_exp), (
                f"Wrong CX mapping for (i={i}, c={c}, f={f}) over GF({p})\n"
                f"A_exp=\n{A_exp}\nGot=\n{FM}\nExp=\n{FM_exp}"
            )

    def local_Du_includes_identity_and_nontrivial(self, p):
        # identity case
        ops = _emit_local_ops_for_D(0, 1, p)
        gates = [op[0] for op in ops]
        qudits = [op[1:] for op in ops]
        M = Circuit.from_gates_and_qudits([p], gates, qudits).composite_gate().symplectic % p
        assert np.array_equal(M, np.eye(2, dtype=int) % p)

        # all nontrivial u
        for u in range(2, p):
            ops = _emit_local_ops_for_D(0, u, p)
            gates = [op[0] for op in ops]
            qudits = [op[1:] for op in ops]
            M = Circuit.from_gates_and_qudits([p], gates, qudits).composite_gate().symplectic % p
            inv_u = pow(u, -1, p)
            target = np.array([[u, 0], [0, inv_u]], dtype=int) % p
            assert np.array_equal(M, target), f"Failed D({u}) over GF({p})"

    def single_wire_scaling_M(self, p: int):
        n = 3
        for u in range(2, p):  # skip 0 and 1
            A = np.eye(n, dtype=int)
            A[1, 1] = u % p
            ops = synth_linear_A_to_gates(n, A, p)
            FM = self.symp_of_gates(n, p, ops)
            FM_exp = self.check_linear_from_A(n, p, A)
            assert np.array_equal(FM, FM_exp), f"Failed D({u}) on wire 1 over GF({p})"

    def A_block_of(self, gate, qudits, p):
        """Get the A block of the symplectic matrix for a gate acting on qudits."""
        F = Circuit.from_gates_and_qudits([p, p], [gate], [qudits]).composite_gate().symplectic % p
        return F[:2, :2]

    def one_col_add_in_M(self, n=3, p=5, c=0, i=2, f=1):
        A = np.eye(n, dtype=int) % p
        A_exp = A.copy()
        A_exp[:, i] = (A_exp[:, i] + f * A_exp[:, c]) % p  # col_i += f col_c

        # Build M ops for just that step:
        ops: list[GateSpec] = [(CX, i, c)] * (f % p)
        gates = [op[0] for op in ops]
        qudits = [op[1:] for op in ops]
        FM = Circuit.from_gates_and_qudits([p] * n, gates, qudits).composite_gate().symplectic % p
        assert np.array_equal(FM[:n, :n], A_exp)

    def test_m_block(self, rng=np.random.default_rng()):
        num_trials = 10
        n = 3
        for p in [2, 3, 5]:
            self.single_wire_scaling_M(p)
            self.local_Du_includes_identity_and_nontrivial(p)
            self.local_Du_synthesis(p)
            # test_SUM_direction_and_sign(n, p)

            # Expect A = [[1,1],[0,1]] for CX acting on (0, 1)
            A = self.A_block_of(CX, (0, 1), p)
            assert np.array_equal(A, np.array([[1, 0], [1, 1]])), A

            for _ in range(num_trials):
                self.CX_direction_in_M(n, p, rng)
                self.M_block_only(n, p, rng)
                self.one_col_add_in_M(n, p)

    def _rand_sym(self, n, p, rng):
        M = rng.integers(0, p, size=(n, n), dtype=int)
        return _as_int_mod(M + M.T, p)

    def _check_equal(self, F_got, F_exp, p, msg=""):
        G = _as_int_mod(F_got, p)
        E = _as_int_mod(F_exp, p)
        if not np.array_equal(G, E):
            raise AssertionError(f"{msg}\nGot=\n{G}\nExp=\n{E}")

    def test_L_block_unit(self):
        """Test L(C) = [[I,0],[C,I]] synthesis correctness."""
        num_trials = 10
        rng = np.random.default_rng()
        n = 3
        for p in [2, 3, 5]:
            for _ in range(num_trials):
                C_sym = self._rand_sym(n, p, rng)
                ops = synth_lower_from_symmetric(n, C_sym, p)
                F = _compose_symp(n, p, ops)
                F_exp = _full_from_lower(n, C_sym, p)
                self._check_equal(F, F_exp, p, f"[L] mismatch over GF({p})")

    def test_R_block_unit(self):
        """Test R(S) = [[I,S],[0,I]] synthesis correctness."""
        n = 3
        num_trials = 10

        rng = np.random.default_rng(42)
        for p in [2, 3, 5]:
            for _ in range(num_trials):
                S_sym = self._rand_sym(n, p, rng)
                ops = synth_upper_from_symmetric_via_H(n, S_sym, p)
                F = _compose_symp(n, p, ops)
                F_exp = self._full_from_upper(n, S_sym, p)
                self._check_equal(F, F_exp, p, f"[R] mismatch over GF({p})")

    def test_elementary_pairs(self):
        """Test elementary diagonal/off-diagonal increments for both L and R."""

        n = 3
        num_trials = 10
        rng = np.random.default_rng()
        for p in [2, 3, 5]:
            for _ in range(num_trials):
                # Diagonal (L)
                for i in range(n):
                    C = np.zeros((n, n), dtype=int)
                    C[i, i] = rng.integers(0, p)
                    Fl = _compose_symp(n, p, synth_lower_from_symmetric(n, C, p))
                    self._check_equal(Fl, _full_from_lower(n, C, p), p, f"[L diag] GF({p})")

                # Off-diagonal (L)
                for i in range(n):
                    for j in range(i + 1, n):
                        C = np.zeros((n, n), dtype=int)
                        t = rng.integers(0, p)
                        C[i, j] = C[j, i] = t % p
                        Fl = _compose_symp(n, p, synth_lower_from_symmetric(n, C, p))
                        self._check_equal(Fl, _full_from_lower(n, C, p), p, f"[L off] GF({p})")

                # Diagonal (R)
                for i in range(n):
                    S = np.zeros((n, n), dtype=int)
                    S[i, i] = rng.integers(0, p)
                    Fr = _compose_symp(n, p, synth_upper_from_symmetric_via_H(n, S, p))
                    self._check_equal(Fr, self._full_from_upper(n, S, p), p, f"[R diag] GF({p})")

    def test_full_decomposition_roundtrip(self):

        seed = 42
        n = 3
        num_trials = 10

        rng = np.random.default_rng(seed)
        for p in [2, 3, 5]:
            for t in range(num_trials):
                F = self.random_symplectic(n, p, rng=rng)
                C = decompose_symplectic_to_circuit(F, p, check=True)
                F_rec = C.composite_gate().symplectic % p
                self._check_equal(F_rec, F, p, f"[roundtrip] GF({p}) trial={t}")

    def test_full_decomposition_on_factored_instances(self):
        """
        Build F explicitly as L(A^{-1}B)·M(A)·R(CA^{-1}) (i.e., from random A, B_sym, C_sym),
        then check we recover it exactly via the decomposition.
        """

        num_trials = 10
        n = 3
        seed = 42

        rng = np.random.default_rng(seed)
        for p in [2, 3, 5]:
            for t in range(num_trials):
                # random invertible A and symmetric B,C
                A = self.random_invertible(n, rng, p)
                B_sym = self.random_symmetric(n, rng, p)
                C_sym = self.random_symmetric(n, rng, p)

                # compose a ground-truth F = L·M·R
                F_true = self.build_F_from_LMR(n, A, B_sym, C_sym, p)

                # Run decomposition
                C = decompose_symplectic_to_circuit(F_true, p, check=True)
                F_rec = C.composite_gate().symplectic % p

                self._check_equal(F_rec, F_true, p, f"[factored-instance] GF({p}) trial={t}")

    def test_phase_correction(self):

        num_trials = 10
        n = 3
        n_gates_in_C_in = 10

        for p in [2, 3, 5]:
            for _ in range(num_trials):
                C_in = Circuit.from_random(n_gates_in_C_in, dimensions=[p] * n)
                big_gate_in = C_in.composite_gate()

                C_out = gate_to_circuit(big_gate_in, dimensions=[p] * n)
                big_gate_out = C_out.composite_gate()

                assert np.array_equal(big_gate_in.symplectic % p, big_gate_out.symplectic % p)

                lcm = int(np.lcm.reduce([p] * n))
                mod_phase = 2 * lcm

                h_in = big_gate_in.phase_vector(p) % mod_phase
                h_out = big_gate_out.phase_vector(p) % mod_phase

                # If p is odd, Pauli can fix any even delta exactly
                if p % 2 == 1:
                    assert np.array_equal(h_in, h_out), f"Fail \nh_in={h_in} \nh_out={h_out}"
                else:
                    # For p=2, Pauli can't move phases mod 4; we can at least assert the
                    assert np.all((h_in - h_out) % 2 == 0), "p=2 residual odd parts cannot be Pauli-corrected"

                for _ in range(num_trials):
                    ps = PauliSum.from_random(2 * n, [p] * n)
                    all_qudits = tuple(range(n))
                    ps1 = big_gate_in.act(ps, all_qudits)
                    ps2 = big_gate_out.act(ps, all_qudits)
                    assert ps1.to_standard_form() == ps2.to_standard_form(), (f"Fail \n{ps1.to_standard_form()}"
                                                                              f"\n{ps2.to_standard_form()}")

    def test_gate_to_circuit_roundtrip_test(self):
        num_trials = 20
        n = 3
        for p in (2, 3, 5):
            for _ in range(num_trials):
                C_in = Circuit.from_random(12, dimensions=[p] * n)
                G_in = C_in.composite_gate()

                C_out = gate_to_circuit(G_in, dimensions=[p] * n)
                G_out = C_out.composite_gate()

                # Symplectic must match exactly mod p
                assert np.array_equal(G_in.symplectic % p, G_out.symplectic % p)

                # Phase behavior: exact for odd p, even-part match for p=2
                lcm = int(np.lcm.reduce([p] * n))
                MOD = 2 * lcm
                h_in = G_in.phase_vector(p) % MOD
                h_out = G_out.phase_vector(p) % MOD

                if p % 2 == 1:
                    assert np.array_equal(h_in, h_out)
                else:
                    # Qubits: Pauli can only move even phases
                    assert np.all((h_in - h_out) % 2 == 0)
