from __future__ import annotations

from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import galois

from sympleq.core.graphs.graph_coloring import _build_base_partition
from sympleq.core.finite_field_solvers import get_linear_dependencies, _select_row_basis_indices
from sympleq.core.paulis import PauliSum
from sympleq.core.circuits import Gate

from .graph_automorphism_kernels import _ConsistencyChecker

# gf2 inverse may live either in this package (older layout) or in the shared solvers
try:  # pragma: no cover
    from sympleq.core.finite_field_solvers import gf2_inv  # type: ignore
except Exception:  # pragma: no cover
    from .graph_automorphism_gf2 import gf2_inv  # type: ignore
from .graph_automorphism_hashing import coeff_ids, choose_base_anchors, compute_anchor_hash
from .graph_automorphism_leaf import LeafContext, check_leaf
from .graph_automorphism_code import compute_induced_completion_matrix_gf2


@dataclass
class PreparedGASearch:
    # problem basics
    pauli: PauliSum
    p: int
    n: int
    labels: list[int]

    # invariants / constraints
    S_mod: np.ndarray
    G: galois.FieldArray
    G_mod2: np.ndarray | None
    basis_order: list[int]
    coeffs: np.ndarray

    # base partition for vertex pruning
    base_colors: np.ndarray
    base_classes: dict[int, list[int]]

    # basis columns (code/matroid information set) as indices in [0..n)
    B_cols: np.ndarray
    is_basis: np.ndarray

    # GF(2) column keys for induced-completion pruning (only if p==2)
    col_keys: list[bytes] | None
    col_bucket: dict[bytes, list[int]] | None

    # fast consistency checker
    consistency: _ConsistencyChecker

    # leaf verification context
    leaf_ctx: LeafContext


def prepare_clifford_ga_search(
    pauli_sum: PauliSum,
    *,
    dynamic_refine_every: int = 0,
    extra_column_invariants: str = "lc",
    p2_bitset: str | bool = "auto",
    color_mode: str = "wl",
    max_wl_rounds: int = 10,
) -> PreparedGASearch:
    """Precompute all data independent of restart seeds.

    This is the expensive part; random-restart runs should call this once per process.
    """
    # canonical gauge
    pauli = pauli_sum.copy()
    pauli.weight_to_phase()

    independent_labels, dependencies = get_linear_dependencies(pauli.tableau, 2)
    labels = sorted(set(independent_labels) | set(dependencies.keys()))

    S_mod = pauli.symplectic_product_matrix()
    G, basis_order = pauli.matroid()
    coeffs = np.asarray(pauli.weights)

    # only uniform dimension supported here (matches your current code)
    if not np.all([pauli.dimensions[i] == pauli.dimensions[0] for i in range(1, len(pauli.dimensions))]):
        raise ValueError(
            "All qudits must have same dimension for now. The key things to fix are: "
            "_gf_solve_one_solution, and the symplectic_solver for F."
        )

    p = int(pauli.lcm)
    n = len(labels)

    # basis columns as indices into the current "labels" order
    lab_to_idx = {lab: i for i, lab in enumerate(labels)}
    B_cols = np.array([lab_to_idx[b] for b in basis_order], dtype=np.int64)
    is_basis = np.zeros(n, dtype=bool)
    is_basis[B_cols] = True

    G_mod2: np.ndarray | None = None
    if p == 2:
        G_mod2 = (np.asarray(G, dtype=np.uint8) & 1)

    # Precompute per-column packed keys for GF(2) induced completion
    col_keys: list[bytes] | None = None
    col_bucket: dict[bytes, list[int]] | None = None
    if G_mod2 is not None:
        # packbits produces a compact representation; bytes are hashable
        col_keys = [np.packbits(G_mod2[:, j], bitorder="little").tobytes() for j in range(n)]
        col_bucket = {}
        for j, k0 in enumerate(col_keys):
            col_bucket.setdefault(k0, []).append(j)

    # ---- extra column invariants for base partition (optional) ----
    col_invariants = None
    if extra_column_invariants != "none":
        toks = {t.strip().lower() for t in extra_column_invariants.replace(",", "+").split("+") if t.strip()}
        toks.discard("none")

        feats: list[np.ndarray] = []

        if "hist" in toks:
            # Heuristic. Useful for ordering / early pruning but not provably complete.
            inv_hist = np.zeros((n, min(p, 16)), dtype=np.int64)
            if p == 2 and G_mod2 is not None:
                G_for_hist = G_mod2
            else:
                G_for_hist = (np.asarray(G, dtype=int) % p)
            for j in range(n):
                col = np.asarray(G_for_hist[:, j], dtype=int)
                cnt = np.bincount(col, minlength=p)
                inv_hist[j, :min(p, 16)] = cnt[:min(p, 16)]
            feats.append(inv_hist)

        if ("lc" in toks) or ("loop_coloop" in toks) or ("loops_coloops" in toks):
            # loop: column is zero
            if p == 2 and G_mod2 is not None:
                is_loop = np.all(G_mod2 == 0, axis=0)
            else:
                G_int = (np.asarray(G, dtype=int) % p)
                is_loop = np.all(G_int == 0, axis=0)

            # coloop: relative to chosen basis_order; computed via basis coordinates
            is_coloop = np.zeros(n, dtype=bool)
            try:
                lab_to_idx = {lab: i for i, lab in enumerate(labels)}
                B_cols = np.array([lab_to_idx[b] for b in basis_order], dtype=int)
                B_mask = np.zeros(n, dtype=bool)
                B_mask[B_cols] = True
                nonB = np.where(~B_mask)[0]

                if p == 2 and G_mod2 is not None:
                    C = G_mod2[:, B_cols]
                    C_inv = gf2_inv(C)
                    X = (C_inv @ G_mod2) & 1
                    basis_row_used = np.any(X[:, nonB] != 0, axis=1) if nonB.size else np.zeros(X.shape[0], dtype=bool)
                else:
                    C = G[:, B_cols]
                    U = np.linalg.inv(C)
                    X = U @ G
                    X_np = np.asarray(X, dtype=int) % p
                    basis_row_used = np.any(X_np[:, nonB] != 0, axis=1) if nonB.size else np.zeros(X_np.shape[0], dtype=bool)

                is_coloop[B_cols] = ~basis_row_used
            except Exception:
                # preserve correctness: no pruning if anything goes wrong
                pass

            inv_lc = np.stack([is_loop.astype(np.int64), is_coloop.astype(np.int64)], axis=1)
            feats.append(inv_lc)

        unknown = toks - {"hist", "lc", "loop_coloop", "loops_coloops"}
        if unknown:
            raise ValueError("extra_column_invariants must be 'none', 'hist', 'lc', or 'hist+lc'.")

        if feats:
            col_invariants = np.hstack(feats) if len(feats) > 1 else feats[0]

    base_colors, base_classes = _build_base_partition(
        S_mod,
        p,
        coeffs=coeffs if coeffs is not None else None,
        col_invariants=col_invariants if color_mode == "wl" else None,
        max_rounds=max_wl_rounds,
        color_mode=color_mode,
    )

    use_bitset = (p == 2 and (p2_bitset is True or (p2_bitset == "auto" and n <= 256)))
    consistency = _ConsistencyChecker(S_mod, use_bitset)

    # ---- leaf precomputation ----
    pauli_standard = pauli.to_standard_form()
    pauli_standard.weight_to_phase()

    ref_tableau = pauli_standard.tableau.astype(int, copy=False)
    ref_phases = np.asarray(pauli_standard.phases, dtype=int)
    ref_weights = np.asarray(pauli_standard.weights)

    base_tableau = pauli.tableau.astype(int, copy=False)
    base_weights = np.asarray(pauli.weights)
    base_phases = np.asarray(pauli.phases, dtype=int)

    basis_indices = np.asarray(independent_labels, dtype=int)
    basis_source_ps = pauli[basis_indices]

    # precompute inverse if full-rank (rank == 2*n_qudits)
    basis_src_inv_gf2: np.ndarray | None = None
    basis_src_inv_gfp: galois.FieldArray | None = None
    basis_src = np.asarray(basis_source_ps.tableau, dtype=int)
    if basis_src.shape[0] == basis_src.shape[1]:
        if p == 2:
            try:
                basis_src_inv_gf2 = gf2_inv(basis_src & 1)
            except np.linalg.LinAlgError:
                basis_src_inv_gf2 = None
        else:
            GF = galois.GF(int(p))
            try:
                basis_src_inv_gfp = cast(galois.FieldArray, np.linalg.inv(GF(basis_src % p)))
            except np.linalg.LinAlgError:
                basis_src_inv_gfp = None

    dims_array = np.asarray(pauli.dimensions, dtype=int)
    row_basis_cache: dict[str, np.ndarray] = {}
    if dims_array.size and np.all(dims_array == dims_array[0]):
        p_uni = int(dims_array[0])
        if p_uni == 2:
            row_basis_cache["gf2"] = _select_row_basis_indices(base_tableau % 2, 2, base_tableau.shape[1])
        else:
            row_basis_cache["gfp"] = _select_row_basis_indices(base_tableau % p_uni, p_uni, base_tableau.shape[1])

    identity_perm = np.arange(pauli.n_paulis(), dtype=np.int64)

    leaf_ctx = LeafContext(
        p=p,
        two_lcm=2 * int(pauli.lcm),
        n_qudits=pauli.n_qudits(),
        identity_perm=identity_perm,
        S_mod=S_mod,
        G=G,
        G_mod2=G_mod2,
        basis_order=basis_order,
        labels=labels,
        pauli_sum=pauli,
        ref_tableau=ref_tableau,
        ref_phases=ref_phases,
        ref_weights=ref_weights,
        base_tableau=base_tableau,
        base_weights=base_weights,
        base_phases=base_phases,
        basis_indices=basis_indices,
        basis_source_ps=basis_source_ps,
        basis_src_inv_gf2=basis_src_inv_gf2,
        basis_src_inv_gfp=basis_src_inv_gfp,
        row_basis_cache=row_basis_cache,
    )

    return PreparedGASearch(
        pauli=pauli,
        p=p,
        n=n,
        labels=labels,
        S_mod=S_mod,
        G=G,
        G_mod2=G_mod2,
        basis_order=basis_order,
        coeffs=coeffs,
        base_colors=base_colors,
        base_classes=base_classes,
        B_cols=B_cols,
        is_basis=is_basis,
        col_keys=col_keys,
        col_bucket=col_bucket,
        consistency=consistency,
        leaf_ctx=leaf_ctx,
    )


def clifford_ga_search_from_prepared(
    prepared: PreparedGASearch,
    *,
    k_wanted: int,
    dynamic_refine_every: int = 0,
    random_seed: int = 0,
    shuffle_domain_order: bool = True,
    progress: bool = False,
    progress_every: int = 2048,
    stop_event: Any | None = None,
    stop_check_every: int = 4096,
    # --- new toggles ---
    use_basis_first_ordering: bool = False,
    use_code_induced_completion: bool = False,
) -> list[Gate]:
    """Run the DFS search using a prepared context.

    random_seed affects only tie-breaking / candidate ordering.
    """
    pauli = prepared.pauli
    p = prepared.p
    n = prepared.n
    S_mod = prepared.S_mod
    coeffs = prepared.coeffs
    base_colors = prepared.base_colors
    base_classes = prepared.base_classes
    leaf_ctx = prepared.leaf_ctx
    B_cols = prepared.B_cols
    is_basis = prepared.is_basis
    k_basis = int(B_cols.size)
    col_keys = prepared.col_keys
    col_bucket = prepared.col_bucket

    rng = np.random.default_rng(int(random_seed))

    # anchor-based ordering keys (restart-dependent)
    coeff_id = coeff_ids(coeffs)
    base_anchors = choose_base_anchors(base_classes, max_anchors=32, rng=rng)
    anchors = np.array(base_anchors, dtype=np.int64)
    key_hash = compute_anchor_hash(S_mod, anchors, base_colors, coeff_id, seed=int(random_seed))

    # domain order (restart-dependent tie-breaking)
    base_order = sorted(base_classes.keys(), key=lambda c: -len(base_classes[c]))
    domain_buckets = [list(base_classes[c]) for c in base_order]
    if shuffle_domain_order:
        for b in domain_buckets:
            rng.shuffle(b)
    domain_order = [i for bucket in domain_buckets for i in bucket]

    # remaining-candidate counts for MRV
    rem_counts: dict[Any, int] = {}
    if coeffs is None:
        rem_counts.update({c: len(base_classes[c]) for c in base_classes})
    else:
        for idx in range(n):
            key = (int(base_colors[idx]), coeffs[idx])
            rem_counts[key] = rem_counts.get(key, 0) + 1

    # state
    phi = -np.ones(n, dtype=np.int64)
    used = np.zeros(n, dtype=bool)

    mapped_stack = np.empty(n, dtype=np.int64)
    mapped_len = 0

    results: list[Gate] = []
    steps = 0

    consistency = prepared.consistency

    def _bucket_key(idx: int) -> Any:
        if coeffs is None:
            return int(base_colors[idx])
        return (int(base_colors[idx]), coeffs[idx])

    def _dec_count(y_idx: int) -> None:
        if coeffs is None:
            rem_counts[int(base_colors[y_idx])] -= 1
        else:
            key = (int(base_colors[y_idx]), coeffs[y_idx])
            rem_counts[key] -= 1

    def _inc_count(y_idx: int) -> None:
        if coeffs is None:
            rem_counts[int(base_colors[y_idx])] += 1
        else:
            key = (int(base_colors[y_idx]), coeffs[y_idx])
            rem_counts[key] += 1

    basis_mapped_count = 0

    # If active, code_C is the induced-completion matrix C = G[:, pi(B)] over GF(2).
    code_C: np.ndarray | None = None
    code_C_u16: np.ndarray | None = None

    def select_next() -> int:
        """Select next domain vertex.

        If use_basis_first_ordering is enabled, prioritize unmapped basis vertices
        until the basis is fully mapped.
        """
        best_i, best_rem = -1, 10**9

        if use_basis_first_ordering and basis_mapped_count < k_basis:
            for i in domain_order:
                if phi[i] >= 0 or (not bool(is_basis[i])):
                    continue
                if coeffs is None:
                    rem = rem_counts[int(base_colors[i])]
                else:
                    rem = rem_counts.get((int(base_colors[i]), coeffs[i]), 0)
                if rem < best_rem:
                    best_i, best_rem = i, rem
                    if rem <= 1:
                        break
            if best_i >= 0:
                return best_i

        # fallback: original MRV
        best_i, best_rem = -1, 10**9
        for i in domain_order:
            if phi[i] >= 0:
                continue
            if coeffs is None:
                rem = rem_counts[int(base_colors[i])]
            else:
                rem = rem_counts.get((int(base_colors[i]), coeffs[i]), 0)
            if rem < best_rem:
                best_i, best_rem = i, rem
                if rem <= 1:
                    break
        return best_i

    def _code_target_key(i: int) -> bytes:
        """Compute the target column key for i under induced completion.

        Requires code_C_u16 to be set and prepared.G_mod2 to be available.
        """
        assert code_C_u16 is not None
        assert prepared.G_mod2 is not None
        col = prepared.G_mod2[:, i].astype(np.uint16, copy=False)
        t = (code_C_u16 @ col) & 1
        return np.packbits(t.astype(np.uint8, copy=False), bitorder="little").tobytes()

    def _try_induced_complete() -> Gate | None:
        """Attempt to complete pi deterministically using g_{pi(i)} = C g_i.

        If every induced target has a unique image (after coefficient filtering) and
        agrees with current partial assignments, build the full permutation and
        run the full leaf verification. Returns a Gate if successful.

        If completion is ambiguous (duplicate targets), returns None and the caller
        should continue with DFS, using code-based candidate restriction.
        """
        if not use_code_induced_completion:
            return None
        if prepared.G_mod2 is None or col_keys is None or col_bucket is None or code_C_u16 is None:
            return None

        pi_full = -np.ones(n, dtype=np.int64)
        used_full = np.zeros(n, dtype=bool)

        # seed with current assignments
        for t in range(mapped_len):
            i0 = int(mapped_stack[t])
            y0 = int(phi[i0])
            if y0 < 0:
                continue
            pi_full[i0] = y0
            used_full[y0] = True

        # try to fill remaining entries
        for i0 in range(n):
            tgt = _code_target_key(i0)
            # candidates are all columns with matching key
            cand = list(col_bucket.get(tgt, []))
            if coeffs is not None:
                cand = [j for j in cand if coeffs[j] == coeffs[i0]]

            if len(cand) != 1:
                return None
            y0 = int(cand[0])
            if pi_full[i0] >= 0 and pi_full[i0] != y0:
                return None
            if pi_full[i0] < 0:
                if used_full[y0]:
                    return None
                pi_full[i0] = y0
                used_full[y0] = True

        # full candidate permutation constructed; run full verification
        return check_leaf(pi_full, leaf_ctx)

    # optional dynamic refine: update ordering keys by individualizing mapped vertices
    def dynamic_refine() -> tuple[np.ndarray, np.ndarray]:
        nonlocal anchors
        if dynamic_refine_every <= 0:
            return anchors, key_hash
        mapped = mapped_stack[:mapped_len]
        K = 16
        extra = mapped[:K]
        anchors_new = np.unique(np.concatenate([anchors, extra]))
        Amax = 64
        if anchors_new.size > Amax:
            anchors_new = anchors_new[:Amax]
        anchors = anchors_new
        return anchors, compute_anchor_hash(S_mod, anchors, base_colors, coeff_id, seed=int(random_seed))

    @dataclass
    class _DFSFrame:
        i: int
        bi: int
        mapped_len: int
        candidate: list[int]
        idx: int = 0
        assigned_y: int = -1

    def _undo_assignment(frame: _DFSFrame) -> None:
        nonlocal mapped_len, basis_mapped_count, code_C, code_C_u16
        y = int(frame.assigned_y)
        if y < 0:
            return
        mapped_len -= 1
        phi[frame.i] = -1
        used[y] = False
        _inc_count(y)
        if bool(is_basis[frame.i]):
            basis_mapped_count -= 1
            # changing basis mapping invalidates induced-completion matrix
            code_C = None
            code_C_u16 = None
        frame.assigned_y = -1

    def _make_frame() -> _DFSFrame | None:
        nonlocal steps, key_hash
        if len(results) >= k_wanted:
            return None

        if dynamic_refine_every and (steps % dynamic_refine_every == 0):
            _, key_hash = dynamic_refine()
        steps += 1

        i = int(select_next())
        if i < 0:
            return None
        bi = int(base_colors[i])
        frame_mapped_len = mapped_len

        candidate = [y for y in base_classes[bi] if not used[y]]
        if coeffs is not None:
            candidate = [y for y in candidate if coeffs[i] == coeffs[y]]

        # If induced completion is active, restrict candidates by code target.
        if use_code_induced_completion and (code_C_u16 is not None) and (col_keys is not None):
            tgt = _code_target_key(i)
            candidate = [y for y in candidate if col_keys[y] == tgt]

        cand = np.array(candidate, dtype=np.int64)
        if cand.size:
            cand = cand[np.argsort(key_hash[cand], kind="mergesort")]
            candidate = cand.tolist()

        return _DFSFrame(i=i, bi=bi, mapped_len=int(frame_mapped_len), candidate=candidate)

    # progress bar (optional)
    progress_bar = None
    pending_progress = 0
    leaves_checked = 0
    if progress:
        try:
            from tqdm.auto import tqdm

            progress_bar = tqdm(
                total=None,
                desc="Clifford automorphism search",
                unit="node",
                leave=False,
                dynamic_ncols=True,
                mininterval=0.2,
            )
            progress_bar.set_postfix(found=0, leaves=0, refresh=False)
        except Exception:
            progress_bar = None

    stack: list[_DFSFrame] = []
    loop_iters = 0

    try:
        while True:
            loop_iters += 1

            if stop_event is not None and (loop_iters % int(stop_check_every) == 0):
                try:
                    if stop_event.is_set():
                        break
                except Exception:
                    # if a non-Event is passed, ignore
                    pass

            if progress_bar is not None:
                pending_progress += 1
                if pending_progress >= progress_every:
                    progress_bar.update(pending_progress)
                    pending_progress = 0
                    progress_bar.set_postfix(found=len(results), leaves=leaves_checked, refresh=False)

            if len(results) >= k_wanted:
                break

            # leaf
            if mapped_len == n:
                leaves_checked += 1
                pi = phi.copy()
                leaf = check_leaf(pi, leaf_ctx)
                if leaf is not None:
                    results.append(leaf)
                    break  # preserve behavior: stop after first symmetry
                if not stack:
                    break
                _undo_assignment(stack[-1])
                continue

            # ensure a frame
            if not stack or stack[-1].assigned_y >= 0:
                fr = _make_frame()
                if fr is None:
                    if not stack:
                        break
                    _undo_assignment(stack[-1])
                    continue
                stack.append(fr)

            frame = stack[-1]

            assigned = False
            while frame.idx < len(frame.candidate):
                y = int(frame.candidate[frame.idx])
                frame.idx += 1
                if used[y]:
                    continue
                if not consistency(phi, mapped_stack, frame.mapped_len, frame.i, y):
                    continue

                # If induced completion is active, reject y that violates the induced target.
                if use_code_induced_completion and (code_C_u16 is not None) and (col_keys is not None):
                    if col_keys[y] != _code_target_key(frame.i):
                        continue

                phi[frame.i] = y
                used[y] = True
                _dec_count(y)
                frame.assigned_y = y

                # update basis bookkeeping and maybe activate induced completion
                if bool(is_basis[frame.i]):
                    basis_mapped_count += 1
                    # basis mapping changed => reset induced completion until re-activated
                    code_C = None
                    code_C_u16 = None

                # activate induced completion once basis is fully mapped
                if (
                    use_code_induced_completion
                    and code_C is None
                    and prepared.G_mod2 is not None
                    and basis_mapped_count == k_basis
                ):
                    C = compute_induced_completion_matrix_gf2(prepared.G_mod2, B_cols, phi)
                    if C is None:
                        # invalid basis image
                        phi[frame.i] = -1
                        used[y] = False
                        _inc_count(y)
                        frame.assigned_y = -1
                        if bool(is_basis[frame.i]):
                            basis_mapped_count -= 1
                        continue
                    code_C = C
                    code_C_u16 = code_C.astype(np.uint16, copy=False)

                    # Check already mapped columns are consistent with the induced rule.
                    if col_keys is not None:
                        ok = True
                        for t in range(frame.mapped_len):
                            j = int(mapped_stack[t])
                            yj = int(phi[j])
                            if yj < 0:
                                continue
                            if col_keys[yj] != _code_target_key(j):
                                ok = False
                                break
                        if not ok:
                            # rollback
                            code_C = None
                            code_C_u16 = None
                            phi[frame.i] = -1
                            used[y] = False
                            _inc_count(y)
                            frame.assigned_y = -1
                            if bool(is_basis[frame.i]):
                                basis_mapped_count -= 1
                            continue

                    # Try deterministic completion; if it succeeds, we are done.
                    gate = _try_induced_complete()
                    if gate is not None:
                        results.append(gate)
                        assigned = True
                        break
                mapped_stack[mapped_len] = frame.i
                mapped_len += 1
                assigned = True
                break

            if assigned:
                continue

            # no candidates
            stack.pop()
            if not stack:
                break
            _undo_assignment(stack[-1])
    finally:
        if progress_bar is not None:
            if pending_progress:
                progress_bar.update(pending_progress)
            progress_bar.set_postfix(found=len(results), leaves=leaves_checked, refresh=False)
            progress_bar.close()

    return results[:k_wanted]


def clifford_graph_automorphism_search(
    pauli_sum: PauliSum,
    k_wanted: int,
    dynamic_refine_every: int = 0,
    extra_column_invariants: str = "lc",
    p2_bitset: str | bool = "auto",
    color_mode: str = "wl",
    max_wl_rounds: int = 10,
    progress: bool = False,
    progress_every: int = 2048,
    # --- new for random restarts / parallel ---
    random_seed: int = 0,
    shuffle_domain_order: bool = False,
    stop_event: Any | None = None,
    stop_check_every: int = 4096,
    # --- toggles ---
    use_basis_first_ordering: bool = False,
    use_code_induced_completion: bool = False,
) -> list[Gate]:
    """Convenience wrapper: prepares and runs a single search."""
    prepared = prepare_clifford_ga_search(
        pauli_sum,
        dynamic_refine_every=dynamic_refine_every,
        extra_column_invariants=extra_column_invariants,
        p2_bitset=p2_bitset,
        color_mode=color_mode,
        max_wl_rounds=max_wl_rounds,
    )
    return clifford_ga_search_from_prepared(
        prepared,
        k_wanted=k_wanted,
        dynamic_refine_every=dynamic_refine_every,
        random_seed=random_seed,
        shuffle_domain_order=shuffle_domain_order,
        progress=progress,
        progress_every=progress_every,
        stop_event=stop_event,
        stop_check_every=stop_check_every,
        use_basis_first_ordering=use_basis_first_ordering,
        use_code_induced_completion=use_code_induced_completion,
    )
