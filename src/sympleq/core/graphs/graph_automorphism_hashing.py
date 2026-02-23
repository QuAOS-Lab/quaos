from __future__ import annotations

import numpy as np


def coeff_ids(coeffs: np.ndarray | None) -> np.ndarray | None:
    if coeffs is None:
        return None
    _, ids = np.unique(coeffs, return_inverse=True)
    return ids.astype(np.int64)


def choose_base_anchors(
    base_classes: dict[int, list[int]],
    max_anchors: int,
    rng: np.random.Generator | None = None,
) -> list[int]:
    """Choose anchors from small colour classes first.

    Anchors are only used for candidate ordering (hash keys), so any choice is ok.
    Using small classes is typically more discriminative.
    """
    anchors: list[int] = []
    for c in sorted(base_classes.keys(), key=lambda cc: len(base_classes[cc])):
        cls = base_classes[c]
        if not cls:
            continue
        if rng is None:
            anchors.append(int(cls[0]))
        else:
            anchors.append(int(cls[int(rng.integers(0, len(cls)))]))
        if len(anchors) >= max_anchors:
            break
    return anchors


def compute_anchor_hash(
    S_mod: np.ndarray,
    anchors: np.ndarray,
    base_colors: np.ndarray,
    coeff_id: np.ndarray | None,
    seed: int = 0,
) -> np.ndarray:
    """Compute a 64-bit hash per vertex used only for ordering.

    Hash includes (base_color, coeff_id?, S[v,anchors], S[anchors,v]). Collisions are fine.
    Complexity: O(M * |anchors|).
    """
    M = S_mod.shape[0]
    A = anchors.size
    width = 1 + (1 if coeff_id is not None else 0) + 2 * A

    feats = np.empty((M, width), dtype=np.int16)
    col = 0
    feats[:, col] = base_colors.astype(np.int16, copy=False)
    col += 1
    if coeff_id is not None:
        feats[:, col] = coeff_id.astype(np.int16, copy=False)
        col += 1
    feats[:, col:col + A] = S_mod[:, anchors].astype(np.int16, copy=False)
    col += A
    feats[:, col:col + A] = S_mod[anchors, :].T.astype(np.int16, copy=False)

    rng = np.random.default_rng(seed)
    w = rng.integers(1, np.iinfo(np.uint64).max, size=width, dtype=np.uint64)
    h = (feats.astype(np.uint64) * w).sum(axis=1, dtype=np.uint64)
    return h
