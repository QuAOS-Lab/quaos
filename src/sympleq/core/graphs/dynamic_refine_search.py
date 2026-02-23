import numpy as np


def dynamic_refine_individualized(
    S_mod: np.ndarray,
    p: int,
    base_colors: np.ndarray,
    coeffs: np.ndarray | None,
    phi: np.ndarray,          # length M, -1 for unmapped, else image index
    used: np.ndarray,         # length M, True if used as image
    max_rounds: int = 2,
    individualize: str = "domain+image",   # "domain" | "image" | "domain+image"
) -> np.ndarray:
    """
    Return a state-dependent color vector (length M) for *ordering candidates only*.

    Idea:
      - Start from base_colors (already includes WL + coeff labels if you built it that way).
      - 'Individualize' already-mapped domain vertices and/or already-used image vertices
        by assigning them unique fresh colors.
      - Run a couple of 1-WL refinement rounds on S_mod using these individualized colors.

    This does NOT affect correctness: it's purely a heuristic ordering.
    """

    M = S_mod.shape[0]
    colors = np.asarray(base_colors, dtype=np.int64).copy()

    # Optional: bake coeffs into initial colors if base_colors didn't already do it.
    # If your base_colors already includes coeff partitioning, you can delete this block.
    if coeffs is not None:
        # Make a stable combined label: (color, coeff_id)
        # Use factor large enough to avoid collisions; M is safe.
        # If coeffs are floats/complex, you should pre-hash them elsewhere.
        try:
            coeff_view = coeffs
            # map arbitrary coeff objects -> integers deterministically
            _, coeff_ids = np.unique(coeff_view, return_inverse=True)
            colors = colors * (M + 1) + coeff_ids.astype(np.int64)
        except Exception:
            pass

    # Individualize mapped/used vertices: give them unique colors so WL “sees” them
    next_color = int(colors.max()) + 1

    if individualize in ("domain", "domain+image"):
        dom = np.where(phi >= 0)[0]
        for v in dom:
            colors[v] = next_color
            next_color += 1

    if individualize in ("image", "domain+image"):
        img = np.where(used)[0]
        # Don't reuse colors already assigned to domain individualized vertices
        for v in img:
            colors[v] = next_color
            next_color += 1

    # 1-WL refinement rounds on complete directed edge-colored graph with edge labels S_mod
    # We update each vertex color by hashing the multiset of (edge_color, neighbor_color).
    # For a complete graph this is O(M^2) per round, but you only do 1–2 rounds.
    for _ in range(max_rounds):
        # Build signatures: for each i, count neighbor colors per edge label.
        # We'll do it with a stable hashing trick based on sorting pairs.
        sigs = []
        for i in range(M):
            # pair (S_ij, colors[j]) for all j
            pair = np.stack([S_mod[i].astype(np.int64), colors], axis=1)
            # sort rows lexicographically to get a canonical multiset representation
            pair = pair[np.lexsort((pair[:, 1], pair[:, 0]))]
            sigs.append(pair.reshape(-1))  # flatten

        # Compress signatures to new colors (stable)
        # Convert list of 1D arrays into a single 2D array with object dtype keys:
        keys = [tuple(s.tolist()) for s in sigs]
        _, new_colors = np.unique(keys, return_inverse=True)
        new_colors = new_colors.astype(np.int64)

        # If nothing changes, stop early
        if np.array_equal(new_colors, colors):
            break
        colors = new_colors

    return colors

