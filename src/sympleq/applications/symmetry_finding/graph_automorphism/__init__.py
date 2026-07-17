from __future__ import annotations

from typing import Any

__all__ = [
    "find_graph_automorphism_clifford_symmetries",
    "find_igraph_clifford_symmetries",
]


def __getattr__(name: str) -> Any:
    if name == "find_igraph_clifford_symmetries":
        from .igraph_automorphism import find_igraph_clifford_symmetries

        return find_igraph_clifford_symmetries
    if name == "find_graph_automorphism_clifford_symmetries":
        from .graph_automorphism_search import clifford_graph_automorphism_search

        return clifford_graph_automorphism_search
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
