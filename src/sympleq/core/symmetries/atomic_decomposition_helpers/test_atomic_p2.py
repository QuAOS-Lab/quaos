# sympleq/core/symmetries/atomic_decomposition_helpers/test_atomic_p2.py
"""
Phase 4 -- independent verification of the p=2 atomic decomposition.

Three layers:
  * relation / distinction tests (review Table 1): conjugate presentations must
    give equal (and correct) qudit cost; the W/V dichotomy must separate the
    classes it should;
  * brute-force oracle comparison: pipeline cost AND the closed-form invariant
    lower bound must equal an exhaustive, theory-free ground truth for 2n <= 8
    (skipped when the invariant-subspace lattice is too large to brute force);
  * property tests: conjugation invariance, idempotence, independent
    verification, and that the implemented unipotent sectors certify minimal.

Runs under pytest, or standalone via ``python3 test_atomic_p2.py`` (the module
ships its own runner because pytest is not always available).
"""
from __future__ import annotations

import numpy as np

from sympleq.core.symmetries.atomic_decomposition import atomic_block_decompose
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import verify_global_basis
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2_generators import (
    canonical_unipotent_p2_block as B,
    direct_sum_unipotent_p2_blocks as DS,
    random_symplectic_conjugate_p2 as conj,
)
from sympleq.core.symmetries.atomic_decomposition_helpers._brute_force_oracle import (
    brute_force_cost,
    OracleInfeasible,
)
from sympleq.core.symmetries.modular_helpers import (
    mod_p,
    inv_mod_mat,
    omega_matrix,
    is_symplectic,
)

P = 2
ORACLE_CAP = 3000


class _Skip(Exception):
    """Raised to signal an intentionally skipped check (mapped to pytest.skip)."""


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _decompose(F):
    Sigma, Bm, info = atomic_block_decompose(F, P)
    verify_global_basis(F, Bm, Sigma, P)          # independent symplectic check
    return Sigma, Bm, info


def _cost(info) -> int:
    return int(info["cost_certificate"]["qudit_cost"])


def _lower_bound(info):
    return info["cost_certificate"]["lower_bound"]


def _blocks_meta(info):
    out = []
    for inv in info.get("sector_invariants", []):
        for b in (getattr(inv, "data", {}) or {}).get("blocks", []):
            out.append((str(b["type"]), int(b["L"]), int(b["half_dim"])))
    return tuple(sorted(out))


def _fingerprint(info):
    """
    Conjugation-invariant signature of a decomposition.

    Note: the per-block alpha/beta refinement (``V`` vs ``V_beta``) is NOT
    individually conjugation-invariant -- only the *total Arf* (sum of beta mod 2)
    within each (family, length) class is (review relation R2:
    ``V_alpha(2k)^2 ~= V(2k)^2``).  So we key on the coarse family (V/W), the
    length multiset, the half-dimensions, and the Arf per (family, length).
    """
    from collections import Counter

    fam_count: Counter = Counter()
    arf: Counter = Counter()
    for (t, L, _hd) in _blocks_meta(info):
        fam = "V" if t.startswith("V") else "W"
        beta = 1 if t.endswith("_beta") else 0
        fam_count[(fam, L)] += 1
        arf[(fam, L)] ^= beta
    return (
        _cost(info),
        tuple(sorted(info.get("atomic_half_dims", []))),
        tuple(sorted(fam_count.items())),
        tuple(sorted(arf.items())),
    )


# --------------------------------------------------------------------------- #
# 1. Relation / distinction tests  (review Table 1)
# --------------------------------------------------------------------------- #
def test_R1_V3_equals_W_perp_V_k1():
    """V(2)^3  ~=  W(2) (+) V(2),  both qudit cost 1 (2n = 6)."""
    a = DS([("V_beta", 2, 0)] * 3)
    b = DS([("W", 2, None), ("V_beta", 2, 0)])
    _, _, ia = _decompose(a)
    _, _, ib = _decompose(b)
    assert _cost(ia) == 1, _cost(ia)
    assert _cost(ib) == 1, _cost(ib)
    assert ia["certified_minimal"] and ib["certified_minimal"]


def test_R1_V3_equals_W_perp_V_k2_pipeline():
    """V(4)^3  ~=  W(4) (+) V(4),  both qudit cost 2 (2n = 12, oracle out of range)."""
    a = DS([("V_beta", 4, 0)] * 3)
    b = DS([("W", 4, None), ("V_beta", 4, 0)])
    _, _, ia = _decompose(a)
    _, _, ib = _decompose(b)
    assert _cost(ia) == 2 and _cost(ib) == 2
    assert ia["certified_minimal"] and ib["certified_minimal"]


def test_R2_Valpha2_sq_equals_V2_sq():
    """V_alpha(2)^2  ~=  V(2)^2,  both cost 1."""
    a = DS([("V_beta", 2, 1)] * 2)
    b = DS([("V_beta", 2, 0)] * 2)
    _, _, ia = _decompose(a)
    _, _, ib = _decompose(b)
    assert _cost(ia) == 1 and _cost(ib) == 1


def test_D1_W2_vs_V2sq_distinguished():
    """W(2) has cost 2; V(2)^2 has cost 1 -- separated only by the top form b_2."""
    _, _, iw = _decompose(B("W", 2))
    _, _, iv = _decompose(DS([("V_beta", 2, 0)] * 2))
    assert _cost(iw) == 2 and iw["certified_minimal"]
    assert _cost(iv) == 1 and iv["certified_minimal"]
    assert _cost(iw) != _cost(iv)


def test_D1_W4_vs_V4sq_distinguished():
    """W(4) cost 4 vs V(4)^2 cost 2 (2n = 8)."""
    _, _, iw = _decompose(B("W", 4))
    _, _, iv = _decompose(DS([("V_beta", 4, 0)] * 2))
    assert _cost(iw) == 4 and _cost(iv) == 2


def test_D2_W3_vs_Wbeta3_same_cost():
    """W(3) and W_beta(3) have the same cost 3 (odd length -> always W-type)."""
    _, _, iw = _decompose(B("W", 3))
    _, _, iwb = _decompose(B("W_beta", 3, 1))
    assert _cost(iw) == 3 and _cost(iwb) == 3


def test_pure_W_not_split_into_V():
    """A pure W(2k) (alternating top form) must NOT be re-decomposed into V's."""
    for k in (1, 2):
        _, _, info = _decompose(B("W", 2 * k))
        assert _cost(info) == 2 * k
        assert ("W", 2 * k, 2 * k) in _blocks_meta(info)


# --------------------------------------------------------------------------- #
# 2. Brute-force oracle comparison (ground truth, 2n <= 8)
# --------------------------------------------------------------------------- #
def _oracle_battery():
    return {
        "I2": np.eye(2, dtype=np.int64),
        "I4": np.eye(4, dtype=np.int64),
        "Vb2_0": B("V_beta", 2, 0),
        "Vb2_1": B("V_beta", 2, 1),
        "Vb4_0": B("V_beta", 4, 0),
        "Vb4_1": B("V_beta", 4, 1),
        "W2": B("W", 2),
        "W3": B("W", 3),
        "W4": B("W", 4),
        "Wb3": B("W_beta", 3, 1),
        "W2+V2": DS([("W", 2, None), ("V_beta", 2, 0)]),
        "V2^3": DS([("V_beta", 2, 0)] * 3),
        "V4+V2": DS([("V_beta", 4, 0), ("V_beta", 2, 0)]),
        "V4+V4": DS([("V_beta", 4, 0)] * 2),
    }


def test_oracle_matches_pipeline_cost_and_lower_bound():
    """Pipeline cost AND the closed-form invariant lower bound equal the exact
    brute-force optimum on every feasible case."""
    checked = 0
    for name, F in _oracle_battery().items():
        try:
            truth = brute_force_cost(F, P, cap=ORACLE_CAP)
        except OracleInfeasible:
            continue  # lattice too large to brute force; skip this case
        _, _, info = _decompose(F)
        assert _cost(info) == truth, f"{name}: pipeline {_cost(info)} != oracle {truth}"
        assert _lower_bound(info) == truth, f"{name}: lower_bound {_lower_bound(info)} != oracle {truth}"
        checked += 1
    assert checked >= 8, f"too few oracle cases were feasible ({checked})"


def test_oracle_matches_on_random_conjugates():
    """Conjugates have the same exact cost; pipeline and oracle still agree."""
    bases = {
        "W2+V2": DS([("W", 2, None), ("V_beta", 2, 0)]),
        "V2^3": DS([("V_beta", 2, 0)] * 3),
        "W3": B("W", 3),
        "V4+V2": DS([("V_beta", 4, 0), ("V_beta", 2, 0)]),
    }
    checked = 0
    for name, F0 in bases.items():
        try:
            truth = brute_force_cost(F0, P, cap=ORACLE_CAP)
        except OracleInfeasible:
            continue
        for seed in (1, 2, 3):
            Fc, _ = conj(F0, seed=seed)
            try:
                tc = brute_force_cost(Fc, P, cap=ORACLE_CAP)
            except OracleInfeasible:
                continue
            assert tc == truth, f"{name} seed{seed}: oracle not conjugation-invariant"
            _, _, info = _decompose(Fc)
            assert _cost(info) == truth, f"{name} seed{seed}: pipeline {_cost(info)} != oracle {truth}"
            checked += 1
    assert checked >= 6, f"too few conjugate cases checked ({checked})"


# --------------------------------------------------------------------------- #
# 3. Property tests
# --------------------------------------------------------------------------- #
def _property_battery():
    return {
        "W2": B("W", 2),
        "W3": B("W", 3),
        "W4": B("W", 4),
        "Vb4_0": B("V_beta", 4, 0),
        "Vb4_1": B("V_beta", 4, 1),
        "Wb3": B("W_beta", 3, 1),
        "W2+V2": DS([("W", 2, None), ("V_beta", 2, 0)]),
        "V2^3": DS([("V_beta", 2, 0)] * 3),
        "WV4": DS([("W", 4, None), ("V_beta", 4, 0)]),
        "mixed": DS([("W", 3, None), ("V_beta", 4, 0), ("W", 2, None)]),
    }


def test_conjugation_invariance_of_fingerprint():
    for name, F0 in _property_battery().items():
        if F0.shape[0] // 2 < 2:
            continue
        _, _, base = _decompose(F0)
        fp = _fingerprint(base)
        for seed in (1, 2, 4, 7):
            Fc, _ = conj(F0, seed=seed)
            _, _, info = _decompose(Fc)
            assert _fingerprint(info) == fp, f"{name} seed{seed}: fingerprint changed {fp} -> {_fingerprint(info)}"


def test_idempotence():
    """Re-decomposing the block-diagonal Sigma yields the same cost."""
    for name, F0 in _property_battery().items():
        Sigma, _, info = _decompose(F0)
        _, _, info2 = _decompose(Sigma)
        assert _cost(info2) == _cost(info), f"{name}: cost changed under re-decomposition"


def test_implemented_sectors_certify_minimal():
    """Every unipotent case is built by the deterministic family, so it must be
    certified AND certified minimal (lower_bound == attained); a fired rad b_L
    assertion or a non-additive kernel profile would flip these to False."""
    for name, F0 in _property_battery().items():
        _, _, info = _decompose(F0)
        cc = info["cost_certificate"]
        assert info["certified"], f"{name}: not certified"
        assert info["certified_minimal"], f"{name}: not certified minimal"
        assert cc["lower_bound"] == cc["qudit_cost"] == cc["attained"], f"{name}: {cc}"


# --------------------------------------------------------------------------- #
# 4. Odd-p coverage (paired & odd-p self sectors) via the same oracle
# --------------------------------------------------------------------------- #
# No odd-p block generators exist in the package, so we build a few symplectic
# fixtures directly: companion/paired blocks, a U/V-interleaving symplectic
# direct sum, and a symplectic conjugator from products of transvections.

def _companion(coeffs_low_to_high, p):
    d = len(coeffs_low_to_high)
    C = np.zeros((d, d), dtype=np.int64)
    for i in range(1, d):
        C[i, i - 1] = 1
    for i in range(d):
        C[i, d - 1] = (-int(coeffs_low_to_high[i])) % p
    return C


def _paired_block(A, p):
    """diag(A, A^{-T}) -- symplectic for any invertible A (a paired sector when
    A's minimal polynomial q satisfies q != q*)."""
    A = mod_p(A, p)
    d = A.shape[0]
    B = mod_p(inv_mod_mat(A, p).T, p)
    F = np.zeros((2 * d, 2 * d), dtype=np.int64)
    F[:d, :d] = A
    F[d:, d:] = B
    return F


def _sym_dsum(blocks, p):
    """Symplectic direct sum: re-orders each 2d_i block into the global
    (U_1..U_k | V_1..V_k) coordinate layout so the result is symplectic."""
    ds = [b.shape[0] // 2 for b in blocks]
    n = sum(ds)
    F = np.zeros((2 * n, 2 * n), dtype=np.int64)
    off = 0
    for b, d in zip(blocks, ds):
        F[off:off + d, off:off + d] = b[:d, :d]
        F[off:off + d, n + off:n + off + d] = b[:d, d:]
        F[n + off:n + off + d, off:off + d] = b[d:, :d]
        F[n + off:n + off + d, n + off:n + off + d] = b[d:, d:]
        off += d
    return mod_p(F, p)


def _transvection(v, c, n2, p):
    Om = omega_matrix(n2 // 2, p)
    v = v.reshape(-1, 1).astype(np.int64)
    return mod_p(np.eye(n2, dtype=np.int64) - (c * (v @ v.T) @ Om), p)


def _random_symplectic(n2, p, seed, steps=6):
    rng = np.random.default_rng(seed)
    S = np.eye(n2, dtype=np.int64)
    for _ in range(steps):
        v = rng.integers(0, p, size=n2).astype(np.int64)
        if not np.any(v % p):
            continue
        S = mod_p(_transvection(v, int(rng.integers(1, p)), n2, p) @ S, p)
    return S


def _conjugate(F, S, p):
    return mod_p(inv_mod_mat(S, p) @ mod_p(F, p) @ S, p)


def _I(n2):
    return np.eye(n2, dtype=np.int64)


def _odd_p_battery():
    p3_paired_d2 = _paired_block(_companion([2, 1], 3), 3)   # q = x^2+x+2, q != q*
    p5_paired_d1 = _paired_block(np.array([[2]], dtype=np.int64), 5)  # diag(2,3)
    return [
        # (name, F, p, hand_known_cost or None)
        ("p3_I2", _I(2), 3, 1),
        ("p3_I4", _I(4), 3, 1),
        ("p3_transvection", np.array([[1, 1], [0, 1]], dtype=np.int64), 3, 1),
        ("p3_paired_d2", p3_paired_d2, 3, 2),                       # q != q* -> cost deg = 2
        ("p3_paired_d2+J2", _sym_dsum([p3_paired_d2, np.array([[1, 1], [0, 1]], dtype=np.int64)], 3), 3, None),
        ("p3_J2+J2", _sym_dsum([np.array([[1, 1], [0, 1]], dtype=np.int64)] * 2, 3), 3, None),
        ("p5_I2", _I(2), 5, 1),
        ("p5_paired_d1", p5_paired_d1, 5, 1),                       # diag(2,3), q != q* -> cost 1
    ]


def test_oracle_self_validation_odd_p():
    """The oracle reproduces hand-known costs over GF(3)/GF(5) (validates the
    ground truth itself for odd p, independent of the pipeline)."""
    checked = 0
    for name, F, p, want in _odd_p_battery():
        if want is None:
            continue
        try:
            got = brute_force_cost(F, p, cap=ORACLE_CAP)
        except OracleInfeasible:
            continue
        assert got == want, f"{name}: oracle {got} != hand-known {want}"
        checked += 1
    assert checked >= 5, f"too few odd-p oracle self-checks ({checked})"


def test_oracle_matches_pipeline_odd_p():
    """Pipeline cost AND invariant lower bound equal the oracle over GF(3)/GF(5)
    -- the paired and odd-p self sectors' first ground-truth coverage."""
    checked = 0
    for name, F, p, _want in _odd_p_battery():
        assert is_symplectic(F, p), f"{name}: fixture not symplectic"
        try:
            truth = brute_force_cost(F, p, cap=ORACLE_CAP)
        except OracleInfeasible:
            continue
        Sigma, Bm, info = atomic_block_decompose(F, p)
        verify_global_basis(F, Bm, Sigma, p)
        cc = info["cost_certificate"]
        assert int(cc["qudit_cost"]) == truth, f"{name}: pipeline {cc['qudit_cost']} != oracle {truth}"
        assert cc["lower_bound"] == truth, f"{name}: lower_bound {cc['lower_bound']} != oracle {truth}"
        checked += 1
    assert checked >= 6, f"too few odd-p oracle comparisons ({checked})"


def test_oracle_odd_p_conjugation_invariance():
    """Exact cost is invariant under random symplectic conjugation over GF(p),
    and the pipeline tracks it (exercises the paired builder in skew bases)."""
    cases = [
        ("p3_paired_d2", _paired_block(_companion([2, 1], 3), 3), 3),
        ("p3_J2+J2", _sym_dsum([np.array([[1, 1], [0, 1]], dtype=np.int64)] * 2, 3), 3),
        ("p5_paired_d1", _paired_block(np.array([[2]], dtype=np.int64), 5), 5),
    ]
    checked = 0
    for name, F0, p in cases:
        try:
            truth = brute_force_cost(F0, p, cap=ORACLE_CAP)
        except OracleInfeasible:
            continue
        for seed in (1, 2, 3):
            S = _random_symplectic(F0.shape[0], p, seed)
            Fc = _conjugate(F0, S, p)
            assert is_symplectic(Fc, p)
            try:
                tc = brute_force_cost(Fc, p, cap=ORACLE_CAP)
            except OracleInfeasible:
                continue
            assert tc == truth, f"{name} seed{seed}: oracle not conjugation-invariant ({tc} != {truth})"
            Sigma, Bm, info = atomic_block_decompose(Fc, p)
            verify_global_basis(Fc, Bm, Sigma, p)
            assert int(info["cost_certificate"]["qudit_cost"]) == truth, f"{name} seed{seed}: pipeline != oracle"
            checked += 1
    assert checked >= 6, f"too few odd-p conjugate checks ({checked})"


# --------------------------------------------------------------------------- #
# 5. Self / unitary sectors (q self-reciprocal, q != x±1) vs the oracle
# --------------------------------------------------------------------------- #
# diag(A, A^{-T}) with A = companion(self-reciprocal irreducible q) is symplectic
# and routes to a *self* sector: the p=2 Hermitian builder
# (atomic_self_p2_unitary) for p=2, or the odd-p Witt builder (atomic_self)
# otherwise. Both achieve the minimal cost on generic bases; the odd-p builder
# additionally has a known basis-dependent non-minimality on the symmetric
# block-diagonal presentation, characterized by the last test below.

def _self_block(coeffs_low_to_high, p):
    """A = companion(q), F = diag(A, A^{-T}); self sector when q is self-reciprocal."""
    return _paired_block(_companion(coeffs_low_to_high, p), p)


def test_oracle_self_p2_unitary():
    """p=2 Hermitian self-sector (atomic_self_p2_unitary): pipeline cost, lower
    bound and oracle agree and the sector certifies minimal, including under
    symplectic conjugation."""
    F = _self_block([1, 1], 2)  # q = x^2 + x + 1 (self-reciprocal, irreducible)
    assert is_symplectic(F, 2)
    truth = brute_force_cost(F, 2, cap=ORACLE_CAP)
    Sigma, Bm, info = atomic_block_decompose(F, 2)
    verify_global_basis(F, Bm, Sigma, 2)
    cc = info["cost_certificate"]
    assert int(cc["qudit_cost"]) == truth, f"unitary cost {cc['qudit_cost']} != oracle {truth}"
    assert cc["lower_bound"] == truth and info["certified_minimal"]
    for seed in (1, 2, 3):
        Fc = _conjugate(F, _random_symplectic(F.shape[0], 2, seed), 2)
        Sigma2, Bm2, info2 = atomic_block_decompose(Fc, 2)
        verify_global_basis(Fc, Bm2, Sigma2, 2)
        assert int(info2["cost_certificate"]["qudit_cost"]) == truth


def test_oracle_self_witt_generic_basis():
    """Odd-p self-reciprocal self-sector (atomic_self / Witt path) on generic
    bases: random symplectic conjugates achieve the minimal cost (== oracle ==
    lower bound) and certify minimal. This is the builder's main, optimal path."""
    cases = [
        ("p3_x2+1", [1, 0], 3),
        ("p5_x2+x+1", [1, 1], 5),
        ("p3_deg4_cyclotomic5", [1, 1, 1, 1], 3),
    ]
    checked = 0
    for name, coeffs, p in cases:
        F0 = _self_block(coeffs, p)
        try:
            truth = brute_force_cost(F0, p, cap=ORACLE_CAP)
        except OracleInfeasible:
            continue
        for seed in (1, 2, 3):
            Fc = _conjugate(F0, _random_symplectic(F0.shape[0], p, seed), p)
            assert is_symplectic(Fc, p)
            Sigma, Bm, info = atomic_block_decompose(Fc, p)
            verify_global_basis(Fc, Bm, Sigma, p)
            cc = info["cost_certificate"]
            assert int(cc["qudit_cost"]) == truth, f"{name} seed{seed}: cost {cc['qudit_cost']} != oracle {truth}"
            assert cc["lower_bound"] == truth and info["certified_minimal"], f"{name} seed{seed}: not certified minimal"
            checked += 1
    assert checked >= 6, f"too few generic-basis self-sector checks ({checked})"


def test_self_witt_special_basis_known_nonminimal():
    """Characterization of a KNOWN limitation: on the symmetric block-diagonal
    presentation diag(A, A^{-T}), the odd-p Witt self-builder returns a VALID but
    non-minimal decomposition (cost == deg(q) == 2 * oracle). The invariant lower
    bound still equals the oracle, and certified_minimal is correctly False --
    i.e. the gap is reported, never silently certified. (Generic conjugates of the
    same matrix attain the minimum: see test_oracle_self_witt_generic_basis.)
    """
    cases = [("p3_x2+1", [1, 0], 3), ("p5_x2+x+1", [1, 1], 5)]
    for name, coeffs, p in cases:
        F = _self_block(coeffs, p)
        assert is_symplectic(F, p)
        truth = brute_force_cost(F, p, cap=ORACLE_CAP)
        Sigma, Bm, info = atomic_block_decompose(F, p)
        verify_global_basis(F, Bm, Sigma, p)  # decomposition is still valid
        cc = info["cost_certificate"]
        deg = len(coeffs)
        assert cc["lower_bound"] == truth, f"{name}: lower_bound {cc['lower_bound']} != oracle {truth}"
        assert int(cc["qudit_cost"]) == deg == 2 * truth, f"{name}: expected known cost {deg}, got {cc['qudit_cost']}"
        assert not info["certified_minimal"], f"{name}: non-minimal result was wrongly certified"


# --------------------------------------------------------------------------- #
# standalone runner (pytest not required)
# --------------------------------------------------------------------------- #
def _run_all():
    import time
    tests = sorted(k for k, v in globals().items() if k.startswith("test_") and callable(v))
    passed = failed = skipped = 0
    for name in tests:
        t = time.time()
        try:
            globals()[name]()
            print(f"  PASS  {name}  ({time.time() - t:.2f}s)")
            passed += 1
        except _Skip as e:
            print(f"  SKIP  {name}: {e}")
            skipped += 1
        except AssertionError as e:
            print(f"  FAIL  {name}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR {name}: {type(e).__name__}: {e}")
            failed += 1
    print(f"\n{passed} passed, {failed} failed, {skipped} skipped, out of {len(tests)}")
    return failed == 0


if __name__ == "__main__":
    import sys
    sys.exit(0 if _run_all() else 1)
