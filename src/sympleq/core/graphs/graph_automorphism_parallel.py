from __future__ import annotations

from typing import Any, cast
import os
import queue
import multiprocessing as mp

import numpy as np

from .graph_automorphism_search import (
    PreparedGASearch,
    prepare_clifford_ga_search,
    clifford_ga_search_from_prepared,
)
from .graph_automorphism_kernels import _consistent_numba, _consistent_bitset, _build_bitrows_binary


# Globals inherited under fork; avoids pickling large contexts.
_PAR_PREPARED: PreparedGASearch | None = None
_PAR_KWARGS: dict[str, Any] | None = None
_PAR_BASE_SEED: int = 0


def _warmup_numba() -> None:
    """Compile numba kernels before forking so workers avoid compile latency (Linux fork)."""
    S = np.zeros((2, 2), dtype=np.int64)
    phi = -np.ones(2, dtype=np.int64)
    stack = np.zeros(2, dtype=np.int64)
    _ = _consistent_numba(S, phi, stack, 0, 0, 0)
    bits, _ = _build_bitrows_binary(S)
    _ = _consistent_bitset(bits, phi, stack, 0, 0, 0)


def _worker(task_q: Any, result_q: Any, found_event: Any) -> None:
    global _PAR_PREPARED, _PAR_KWARGS, _PAR_BASE_SEED

    prepared = _PAR_PREPARED
    kwargs = _PAR_KWARGS or {}
    base_seed = int(_PAR_BASE_SEED)

    if prepared is None:
        return

    # Each worker loops over multiple restarts to amortize JIT and preprocessing.
    while True:
        try:
            idx = task_q.get_nowait()
        except queue.Empty:
            return

        try:
            if found_event.is_set():
                return
        except Exception:
            pass

        seed = base_seed + int(idx) * 1_000_003
        out = clifford_ga_search_from_prepared(
            prepared,
            k_wanted=1,
            random_seed=seed,
            shuffle_domain_order=True,
            progress=False,
            stop_event=found_event,
            # check the stop flag fairly often; this is for parallel early-exit
            stop_check_every=2048,
            **kwargs,
        )
        if out:
            try:
                result_q.put(out[0])
            except Exception:
                pass
            try:
                found_event.set()
            except Exception:
                pass
            return


def clifford_graph_automorphism_search_random_restarts(
    pauli_sum,
    k_wanted: int = 1,
    n_restarts: int = 8,
    n_jobs: int | None = None,
    base_seed: int = 0,
    start_method: str = "fork",
    warmup: bool = True,
    # parameters forwarded to prepare/search
    dynamic_refine_every: int = 0,
    extra_column_invariants: str = "lc",
    p2_bitset: str | bool = "auto",
    color_mode: str = "wl",
    max_wl_rounds: int = 10,
) -> list:
    """Parallel random-restart search for a *single* symmetry.

    This runs multiple independent restarts with different tie-break seeds in parallel.
    The first worker to find a valid symmetry stops the others.

    Notes
    -----
    - Designed for Linux: uses start_method='fork' so workers inherit the prepared context
      without pickling/copying large arrays.
    - If fork is not available, we fall back to sequential restarts.
    """
    if k_wanted != 1:
        raise NotImplementedError("Parallel random restarts currently implemented for k_wanted=1 only.")

    n_restarts = int(max(1, n_restarts))

    # Precompute the expensive invariants once.
    prepared = prepare_clifford_ga_search(
        pauli_sum,
        extra_column_invariants=extra_column_invariants,
        p2_bitset=p2_bitset,
        color_mode=color_mode,
        max_wl_rounds=max_wl_rounds,
    )

    # If only one restart requested, just run once.
    if n_restarts == 1:
        return clifford_ga_search_from_prepared(
            prepared,
            k_wanted=1,
            random_seed=int(base_seed),
            shuffle_domain_order=True,
            progress=False,
            dynamic_refine_every=dynamic_refine_every,
        )

    if n_jobs is None:
        n_jobs = min(n_restarts, max(1, (os.cpu_count() or 1)))
    n_jobs = int(max(1, n_jobs))

    # If no parallelism desired, do serial restarts (still useful).
    if n_jobs == 1:
        for r in range(n_restarts):
            out = clifford_ga_search_from_prepared(
                prepared,
                k_wanted=1,
                random_seed=int(base_seed) + r * 1_000_003,
                shuffle_domain_order=True,
                progress=False,
                dynamic_refine_every=dynamic_refine_every,
            )
            if out:
                return out
        return []

    # Try to get the requested multiprocessing context.
    try:
        ctx = mp.get_context(start_method)
    except Exception:
        # Fallback to serial if context isn't available.
        for r in range(n_restarts):
            out = clifford_ga_search_from_prepared(
                prepared,
                k_wanted=1,
                random_seed=int(base_seed) + r * 1_000_003,
                shuffle_domain_order=True,
                progress=False,
                dynamic_refine_every=dynamic_refine_every,
            )
            if out:
                return out
        return []

    if warmup and start_method == "fork":
        _warmup_numba()

    # Set globals for forked workers.
    global _PAR_PREPARED, _PAR_KWARGS, _PAR_BASE_SEED
    _PAR_PREPARED = prepared
    _PAR_KWARGS = dict(dynamic_refine_every=dynamic_refine_every)
    _PAR_BASE_SEED = int(base_seed)
    ctx_any = cast(Any, ctx)
    found_event = ctx_any.Event()
    task_q = ctx_any.Queue()
    result_q = ctx_any.Queue(maxsize=1)

    for r in range(n_restarts):
        task_q.put(r)

    procs = [ctx_any.Process(target=_worker, args=(task_q, result_q, found_event)) for _ in range(n_jobs)]
    for p in procs:
        p.daemon = True
        p.start()

    gate = None
    try:
        # Wait for a result while workers run; poll to allow quick exit.
        while True:
            try:
                gate = result_q.get(timeout=0.2)
                break
            except queue.Empty:
                if not any(p.is_alive() for p in procs):
                    break
                continue
    finally:
        try:
            found_event.set()
        except Exception:
            pass
        for p in procs:
            if p.is_alive():
                try:
                    p.terminate()
                except Exception:
                    pass
        for p in procs:
            try:
                p.join(timeout=0.2)
            except Exception:
                pass

    return [gate] if gate is not None else []
