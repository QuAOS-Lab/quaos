import os
import pprint

import numpy as np
import pytest

from sympleq.core.circuits.random_symplectic import symplectic_random_transvection
from sympleq.core.symmetries.atomic_decomposition import CertificationError, atomic_block_decompose
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import (
    inv_mod_mat,
    is_symplectic,
    mod_p,
    rank_mod,
)


def _env_int(name: str, default: int) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return int(default)
    try:
        return int(raw)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer, got {raw!r}") from exc


def _rand_symplectic(rng: np.random.Generator, n: int, p: int, steps: int) -> np.ndarray:
    F = symplectic_random_transvection(n, p, num_transvections=int(steps), rng=rng)
    F = mod_p(F, p)
    assert is_symplectic(F, p)
    return F


def direct_sum_grouped_blocks(blocks: list[np.ndarray], p: int) -> np.ndarray:
    """
    Symplectic direct sum in grouped [x_1,...,x_n | z_1,...,z_n] ordering.

    If block i has shape 2n_i × 2n_i and is written as

        [[A_i, B_i],
         [C_i, D_i]],

    this returns

        [[diag(A_i), diag(B_i)],
         [diag(C_i), diag(D_i)]].
    """
    if not blocks:
        raise ValueError("Need at least one block.")

    ns: list[int] = []
    parts: list[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = []

    for F in blocks:
        F = mod_p(np.asarray(F, dtype=np.int64), p)
        if F.ndim != 2 or F.shape[0] != F.shape[1] or F.shape[0] % 2 != 0:
            raise ValueError(f"Each block must be square with even dimension; got {F.shape}.")
        assert is_symplectic(F, p)

        n = F.shape[0] // 2
        ns.append(n)
        A, B = F[:n, :n], F[:n, n:]
        C, D = F[n:, :n], F[n:, n:]
        parts.append((A, B, C, D))

    total_n = sum(ns)

    def block_diag(mats: list[np.ndarray]) -> np.ndarray:
        out = np.zeros((total_n, total_n), dtype=np.int64)
        off = 0
        for M in mats:
            nr, nc = M.shape
            assert nr == nc
            out[off:off + nr, off:off + nc] = M
            off += nr
        return mod_p(out, p)

    A_big = block_diag([x[0] for x in parts])
    B_big = block_diag([x[1] for x in parts])
    C_big = block_diag([x[2] for x in parts])
    D_big = block_diag([x[3] for x in parts])

    F_big = mod_p(np.block([[A_big, B_big], [C_big, D_big]]), p)
    assert F_big.shape == (2 * total_n, 2 * total_n)
    assert is_symplectic(F_big, p)
    return F_big


def _random_block_sizes(
    rng: np.random.Generator,
    *,
    q_max: int,
    max_blocks: int,
    max_total_n: int,
) -> list[int]:
    """
    Randomly choose a nonempty list of block sizes q_i with q_i <= q_max and
    sum(q_i) <= max_total_n.
    """
    n_blocks = int(rng.integers(1, max_blocks + 1))
    sizes: list[int] = []

    remaining = int(max_total_n)
    for _ in range(n_blocks):
        if remaining <= 0:
            break
        q_hi = min(int(q_max), remaining)
        q = int(rng.integers(1, q_hi + 1))
        sizes.append(q)
        remaining -= q

    if not sizes:
        sizes = [1]
    return sizes


def _make_random_bounded_cost_instance(
    rng: np.random.Generator,
    *,
    p: int,
    block_sizes: list[int],
    block_steps_min: int = 15,
    block_steps_max: int = 70,
    scramble_steps_min: int = 30,
    scramble_steps_max: int = 140,
) -> tuple[np.ndarray, np.ndarray, int, dict]:
    """
    Return (F_direct, F_scrambled, Q_max, metadata).
    """
    blocks: list[np.ndarray] = []
    block_steps: list[int] = []

    for q in block_sizes:
        steps = int(rng.integers(block_steps_min, block_steps_max + 1))
        block_steps.append(steps)
        blocks.append(_rand_symplectic(rng, int(q), p, steps=steps))

    F_direct = direct_sum_grouped_blocks(blocks, p)
    n_total = F_direct.shape[0] // 2
    Q_max = max(int(q) for q in block_sizes)

    scramble_steps = int(rng.integers(scramble_steps_min, scramble_steps_max + 1))
    S = _rand_symplectic(rng, n_total, p, steps=scramble_steps)
    F_scrambled = mod_p(inv_mod_mat(S, p) @ F_direct @ S, p)

    assert is_symplectic(F_direct, p)
    assert is_symplectic(F_scrambled, p)

    meta = {
        "p": int(p),
        "block_sizes": [int(q) for q in block_sizes],
        "Q_max": int(Q_max),
        "n_total": int(n_total),
        "block_steps": block_steps,
        "scramble_steps": int(scramble_steps),
    }
    return F_direct, F_scrambled, Q_max, meta


def _assert_bounded_cost_decomposition(
    F: np.ndarray,
    p: int,
    Q_max: int,
    *,
    trial_meta: dict,
) -> None:
    try:
        Sigma, B, info = atomic_block_decompose(F, p)
    except CertificationError as exc:
        pytest.fail(
            "Certified decomposition failed during randomized bounded-cost test.\n"
            f"trial_meta={pprint.pformat(trial_meta, width=120)}\n"
            f"CertificationError.info={pprint.pformat(exc.info, width=140)}"
        )

    verify_global_basis(F, B, Sigma, p)

    assert info["certified"] is True, pprint.pformat({"trial_meta": trial_meta, "info": info}, width=140)
    assert rank_mod(B, p) == F.shape[0]
    assert sum(int(h) for h in info["atomic_half_dims"]) == F.shape[0] // 2

    Q_cost = int(info.get("Q_att", info.get("attained_qudit_cost", info["qudit_cost"])))
    assert Q_cost == max(int(h) for h in info["atomic_half_dims"])
    if info.get("minimal_cost_certified", False):
        assert info.get("Q_opt") == Q_cost
    else:
        assert info.get("Q_opt") is None
    # assert Q_cost <= int(Q_max), (
    #     "Recovered qudit cost exceeds the known explicit block construction.\n"
    #     f"Q_cost={Q_cost}, Q_max={Q_max}\n"
    #     f"atomic_half_dims={info['atomic_half_dims']}\n"
    #     f"trial_meta={pprint.pformat(trial_meta, width=120)}"
    # )
    if Q_cost > int(Q_max):
        sector_payload = []
        for inv in info.get("sector_invariants", []):
            d = getattr(inv, "data", {})
            sector_payload.append(
                {
                    "sector_type": getattr(inv, "sector_type", None),
                    "sector_key": getattr(inv, "sector_key", None),
                    "status": d.get("status"),
                    "block_half_dims": d.get("block_half_dims"),
                    "block_records": d.get("block_records"),
                    "p2_unipotent": d.get("p2_unipotent"),
                    "unitary": d.get("unitary"),
                    "length_summary": d.get("length_summary"),
                    "cost_certificate": d.get("cost_certificate"),
                }
            )

        pytest.fail(
            "Recovered qudit cost exceeds the known explicit block construction.\n"
            f"Q_cost={Q_cost}, Q_max={Q_max}\n"
            f"atomic_half_dims={info['atomic_half_dims']}\n"
            f"trial_meta={pprint.pformat(trial_meta, width=120)}\n"
            f"sector_payload={pprint.pformat(sector_payload, width=140)}"
        )

    if info.get("minimal_cost_certified", False):
        cert = info.get("cost_certificate", {})
        assert int(cert.get("attained")) == Q_cost, pprint.pformat({"trial_meta": trial_meta, "cert": cert}, width=140)
        assert int(cert.get("qudit_cost", Q_cost)) == Q_cost
        assert int(cert.get("lower_bound")) == Q_cost
        assert cert.get("certified_minimal") is True


class TestAtomicDecompositionRandomBoundedCost:
    def test_many_random_bounded_cost_instances(self) -> None:
        """
        Non-deterministic stress test.

        Every run samples a fresh random collection of bounded-cost decomposable
        symplectic matrices.
        """
        trials = _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_TRIALS", 20)
        q_max_limit = _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_MAX_Q", 4)
        max_blocks = _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_MAX_BLOCKS", 8)
        max_total_n = _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_MAX_TOTAL_N", 18)

        ss = np.random.SeedSequence(_env_int("SYMPLEQ_BOUNDED_COST_RANDOM_SEED", 123456))
        rng = np.random.default_rng(ss)

        primes = [2]

        for trial in range(int(trials)):
            p = int(rng.choice(primes))
            q_max = int(rng.integers(1, q_max_limit + 1))
            block_sizes = _random_block_sizes(
                rng,
                q_max=q_max,
                max_blocks=max_blocks,
                max_total_n=max_total_n,
            )

            _F_direct, F_scrambled, Q_max, meta = _make_random_bounded_cost_instance(
                rng,
                p=p,
                block_sizes=block_sizes,
            )
            meta.update(
                {
                    "trial": int(trial),
                    "seedsequence_entropy": ss.entropy,
                    "q_max_limit": int(q_max_limit),
                    "max_blocks": int(max_blocks),
                    "max_total_n": int(max_total_n),
                }
            )

            _assert_bounded_cost_decomposition(
                F_scrambled,
                p,
                Q_max,
                trial_meta=meta,
            )

    def test_many_random_bounded_cost_instances_unscrambled_sanity(self) -> None:
        """
        Smaller non-deterministic sanity check without scrambling.

        This helps distinguish direct-sum construction failures from failures
        caused by the random conjugating symplectic.
        """
        trials = max(5, _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_TRIALS", 20) // 10)
        q_max_limit = _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_MAX_Q", 4)
        max_blocks = _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_MAX_BLOCKS", 8)
        max_total_n = _env_int("SYMPLEQ_BOUNDED_COST_RANDOM_MAX_TOTAL_N", 18)

        ss = np.random.SeedSequence(_env_int("SYMPLEQ_BOUNDED_COST_RANDOM_SEED", 123456))
        rng = np.random.default_rng(ss)
        primes = [2]

        for trial in range(int(trials)):
            p = int(rng.choice(primes))
            q_max = int(rng.integers(1, q_max_limit + 1))
            block_sizes = _random_block_sizes(
                rng,
                q_max=q_max,
                max_blocks=max_blocks,
                max_total_n=max_total_n,
            )

            F_direct, _F_scrambled, Q_max, meta = _make_random_bounded_cost_instance(
                rng,
                p=p,
                block_sizes=block_sizes,
                scramble_steps_min=1,
                scramble_steps_max=1,
            )
            meta.update(
                {
                    "trial": int(trial),
                    "seedsequence_entropy": ss.entropy,
                    "unscrambled": True,
                    "q_max_limit": int(q_max_limit),
                    "max_blocks": int(max_blocks),
                    "max_total_n": int(max_total_n),
                }
            )

            _assert_bounded_cost_decomposition(
                F_direct,
                p,
                Q_max,
                trial_meta=meta,
            )

def test_scrambled_direct_sum_of_one_qudit_blocks_has_cost_one() -> None:
    p = 2
    rng = np.random.default_rng(12346)

    block_sizes = [1] * 8
    F_direct, F_scrambled, Q_max, meta = _make_random_bounded_cost_instance(
        rng,
        p=p,
        block_sizes=block_sizes,
        block_steps_min=5,
        block_steps_max=60,
        scramble_steps_min=80,
        scramble_steps_max=120,
    )

    assert Q_max == 1
    assert is_symplectic(F_direct, p)
    assert is_symplectic(F_scrambled, p)

    Sigma, B, info = atomic_block_decompose(F_scrambled, p)
    verify_global_basis(F_scrambled, B, Sigma, p)

    assert int(info["Q_opt"]) == 1
    assert sorted(int(h) for h in info["atomic_half_dims"]) == [1] * 8
    assert info["minimal_cost_certified"] is True