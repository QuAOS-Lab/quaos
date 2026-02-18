import numpy as np
import pytest

from sympleq.core.symmetries.atomic_decomposition import CertificationError
from sympleq.core.symmetries.block_decomposition import atomic_block_decompose, block_indexes
from sympleq.core.symmetries.modular_helpers import inv_mod_mat, mod_p, rank_mod
from sympleq.core.circuits.utils import is_symplectic
from sympleq.core.circuits.random_symplectic import symplectic_random_transvection
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_linear import restrict_operator
import pprint


def _omega(n: int, p: int) -> np.ndarray:
    Id = np.eye(n, dtype=int)
    Z = np.zeros((n, n), dtype=int)
    Ω = np.block([[Z, Id], [-Id, Z]]) % p
    return Ω


def direct_sum(F1, F2, p):
    """
    Symplectic direct sum in grouped ordering [x...|z...].

    If F1 = [[A1, B1], [C1, D1]] and F2 = [[A2, B2], [C2, D2]],
    this returns:
        [[diag(A1, A2), diag(B1, B2)],
         [diag(C1, C2), diag(D1, D2)]]
    which acts on [x1, x2, z1, z2].
    """
    if F1.ndim != 2 or F1.shape[0] != F1.shape[1]:
        raise ValueError("F1 must be square.")
    if F2.ndim != 2 or F2.shape[0] != F2.shape[1]:
        raise ValueError("F2 must be square.")
    if F1.shape[0] % 2 != 0 or F2.shape[0] % 2 != 0:
        raise ValueError("F1 and F2 must have even dimension (2n x 2n).")

    n1 = F1.shape[0] // 2
    n2 = F2.shape[0] // 2

    A1, B1 = F1[:n1, :n1], F1[:n1, n1:]
    C1, D1 = F1[n1:, :n1], F1[n1:, n1:]
    A2, B2 = F2[:n2, :n2], F2[:n2, n2:]
    C2, D2 = F2[n2:, :n2], F2[n2:, n2:]

    Z12 = np.zeros((n1, n2), dtype=int)
    Z21 = np.zeros((n2, n1), dtype=int)

    A = np.block([[A1, Z12], [Z21, A2]])
    B = np.block([[B1, Z12], [Z21, B2]])
    C = np.block([[C1, Z12], [Z21, C2]])
    D = np.block([[D1, Z12], [Z21, D2]])

    return np.block([[A, B], [C, D]]) % p


def _span_invariant(F: np.ndarray, T: np.ndarray, p: int) -> bool:
    # row action v -> vF corresponds to col action T -> F^T T
    A = mod_p(F.T @ T, p)
    return rank_mod(mod_p(np.concatenate([T, A], axis=1), p), p) == rank_mod(T, p)


def _is_nondegenerate(Omega: np.ndarray, T: np.ndarray, p: int) -> bool:
    # Nondegenerate on span(T) iff Gram has full rank
    G = mod_p(T.T @ Omega @ T, p)
    return rank_mod(G, p) == T.shape[1]


def _canonicalize_info(info: dict) -> tuple:
    # Keep only deterministic, basis-independent data
    invs = info.get("sector_invariants", [])
    # Each inv is AtomicInvariant; canonicalize using inv.data
    payload = []
    for inv in invs:
        d = inv.data
        # pick stable fields
        payload.append((
            inv.sector_type,
            inv.sector_key,
            d.get("deg"),
            d.get("exponent"),
            tuple(sorted((int(L), dd.get("mult"), dd.get("q1_count", None))
                         for L, dd in d.get("length_summary", {}).items()))
        ))
    payload.sort()
    half_dims = tuple(sorted(info.get("atomic_half_dims", [])))
    return (tuple(payload), half_dims, info.get("Q_opt"), info.get("Lmin_star"), info.get("certified"))


class TestAtomicBlocks:

    @pytest.mark.parametrize("p,n,seed", [(2, 6, 1), (3, 5, 2), (5, 4, 3)])
    def test_atomic_blocks_are_invariant_and_nondegenerate(self, p, n, seed):
        rng = np.random.default_rng(seed)
        n_transvections = 100
        # random symplectic to stress many sectors
        F = symplectic_random_transvection(n, p, num_transvections=n_transvections, rng=rng)
        assert is_symplectic(F, p)

        try:
            Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
        except CertificationError as e:
            import pprint
            pprint.pprint(e.info["failures"])
            raise
        assert info["certified"] is True

        Ω = _omega(n, p)

        # B should be symplectic and actually conjugate F -> Sigma
        assert is_symplectic(B, p)

        # Now verify each block inside Sigma is invariant + nondegenerate
        # If you have block_indexes(Sigma): iterate those blocks and extract T columns accordingly.
        from sympleq.core.symmetries.block_decomposition import block_indexes

        blocks = block_indexes(Sigma)
        for blk in blocks:
            # blk is list of qudit indices; build column selector for those qudits in [x|z] ordering
            cols = []
            for q in blk:
                cols.append(q)          # x_q
            for q in blk:
                cols.append(q + n)      # z_q
            cols = np.array(cols, dtype=int)

            T = np.eye(2 * n, dtype=int)[:, cols] % p  # basis for that block subspace
            assert _span_invariant(Sigma, T, p), "Block subspace not invariant under Sigma"
            assert _is_nondegenerate(Ω, T, p), "Block subspace degenerate under Ω"

    def test_blocks_symplectically_orthogonal_and_span(self):
        p, n = 2, 7
        rng = np.random.default_rng()
        n_transvections = 100
        F = symplectic_random_transvection(n, p, num_transvections=n_transvections, rng=rng)
        try:
            Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
        except CertificationError as e:
            import pprint
            pprint.pprint(e.info["failures"])
            raise
        assert info["certified"]

        Ω = _omega(n, p)
        from sympleq.core.symmetries.block_decomposition import block_indexes

        blocks = block_indexes(Sigma)

        # Build T_i for each block and check orthogonality
        Ts = []
        for blk in blocks:
            cols = np.array([*blk, *(q + n for q in blk)], dtype=int)
            Ts.append((np.eye(2 * n, dtype=int)[:, cols] % p))

        # pairwise orthogonality
        for i in range(len(Ts)):
            for j in range(i + 1, len(Ts)):
                Gij = mod_p(Ts[i].T @ Ω @ Ts[j], p)
                assert not Gij.any(), "Different blocks not symplectically orthogonal"

        # spanning: concatenation should have full rank 2n
        Tall = mod_p(np.concatenate(Ts, axis=1), p)
        assert rank_mod(Tall, p) == 2 * n

    def test_certified_mode_is_deterministic(self):
        p, n = 2, 7
        rng = np.random.default_rng()
        n_transvections = 100
        F = symplectic_random_transvection(n, p, num_transvections=n_transvections, rng=rng)

        try:
            Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
        except CertificationError as e:
            import pprint
            pprint.pprint(e.info["failures"])
            raise

        Sigma1, B1, info1 = atomic_block_decompose(F, p, mode="certified")
        Sigma2, B2, info2 = atomic_block_decompose(F, p, mode="certified")

        assert _canonicalize_info(info1) == _canonicalize_info(info2)

    def test_certified_never_degrades_p2_self_sector(self):
        p, n = 2, 8

        # pick a few Fs; if any triggers fallback, certified should raise
        for seed in [10, 11, 12, 13, 14]:
            rng = np.random.default_rng(seed)
            n_transvections = 100
            F = symplectic_random_transvection(n, p, num_transvections=n_transvections, rng=rng)
            try:
                Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
            except Exception:
                # this is acceptable: certified can fail
                continue

            # if certified returned, it must not have any degraded sector statuses
            for inv in info.get("sector_invariants", []):
                assert inv.data.get("status") == "OK"

    def test_certified_mode_failure_reports_sectors(self):
        p, n = 3, 6
        rng = np.random.default_rng(0)
        F = symplectic_random_transvection(n, p, num_transvections=50, rng=rng)

        try:
            atomic_block_decompose(F, p, mode="certified")
        except CertificationError as e:
            assert "failures" in e.info
            assert len(e.info["failures"]) > 0
            # ensure each failure reports the sector type
            assert all("sector_type" in f for f in e.info["failures"])


    @pytest.mark.parametrize("p,n,num_transv,seed", [
        (2, 4, 20, 0), (2, 4, 100, 1), (2, 6, 200, 2),
        (3, 4, 50, 3), (3, 6, 150, 4),
        (5, 4, 50, 5), (5, 6, 150, 6),
    ])
    def test_certified_stress_grid(self, p, n, num_transv, seed):
        rng = np.random.default_rng(seed)
        F = symplectic_random_transvection(n, p, num_transvections=num_transv, rng=rng)
        Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
        assert info["certified"]

    def test_certified_many_random(self):
        for p in [2, 3, 5]:
            for n in [4, 5, 6]:
                for seed in range(50):
                    rng = np.random.default_rng((p, n, seed))
                    F = symplectic_random_transvection(n, p, num_transvections=200, rng=rng)
                    try:
                        Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
                    except CertificationError as e:
                        print("\n=== CERTIFIED FAILURE ===")
                        print(f"p={p}, n={n}, seed={seed}")
                        pprint.pprint(e.info, width=120)
                        # optional: save reproducer
                        np.save(f"F_fail_p{p}_n{n}_seed{seed}.npy", F)
                        raise

    def test_conjugacy_invariance_of_invariants(self):
        p, n = 3, 6
        rng = np.random.default_rng(0)
        F = symplectic_random_transvection(n, p, num_transvections=120, rng=rng)
        S = symplectic_random_transvection(n, p, num_transvections=120, rng=rng)
        Sinv = inv_mod_mat(S, p)
        F2 = mod_p(Sinv @ F @ S, p)
        assert is_symplectic(F2, p)

        _, _, info1 = atomic_block_decompose(F, p, mode="certified")
        _, _, info2 = atomic_block_decompose(F2, p, mode="certified")

        sig1 = tuple(sorted(info1.get("atomic_half_dims", [])))
        sig2 = tuple(sorted(info2.get("atomic_half_dims", [])))
        assert sig1 == sig2

    def test_known_direct_sum_recovers_blocks(self):
        p, n1, n2 = 3, 2, 10
        rng = np.random.default_rng()
        n_tests = 100
        for seed in range(n_tests):
            rng = np.random.default_rng(seed)
            F1 = symplectic_random_transvection(n1, p, 60, rng)
            F2 = symplectic_random_transvection(n2, p, 60, rng)
            F = direct_sum(F1, F2, p)
            assert is_symplectic(F, p)

            # conjugate to hide the sum
            S = symplectic_random_transvection(n1 + n2, p, 120, rng)
            Sinv = inv_mod_mat(S, p)
            Fh = mod_p(Sinv @ F @ S, p)

            _, _, info = atomic_block_decompose(Fh, p, mode="certified")
            assert np.all(tuple(sorted(info["atomic_half_dims"])) <= tuple(sorted([n1, n2]))), f'block sizes should be at most the original blocks; got {info["atomic_half_dims"]} vs {n1, n2}'

    def test_each_atomic_block_is_indecomposable(self):
        p, n = 2, 7
        rng = np.random.default_rng()
        n_tests = 100
        for _ in range(n_tests):
            F = symplectic_random_transvection(n, p, 250, rng)

            Sigma, B, info = atomic_block_decompose(F, p, mode="certified")
            blocks = block_indexes(Sigma)

            for blk in blocks:
                cols = np.array([*blk, *(q+n for q in blk)], dtype=int)
                T = (np.eye(2 * n, dtype=int)[:, cols] % p)

                # restrict Sigma to that block
                Sig_blk = restrict_operator(Sigma, T, p)

                Sig2, B2, info2 = atomic_block_decompose(Sig_blk, p, mode="certified")
                blocks2 = block_indexes(Sig2)
                assert len(blocks2) == 1
            assert len(blocks2[0]) == len(blk)
