import numpy as np

from sympleq.core.graphs.graph_automorphism_circuits import (
    extract_fundamental_relations_from_matroid,
)


def test_extract_fundamental_relations_from_matroid_gf2_nontrivial_basis():
    G = np.array(
        [
            [0, 1, 1, 1],
            [1, 0, 1, 0],
        ],
        dtype=int,
    )

    relations = extract_fundamental_relations_from_matroid(
        G,
        np.array([1, 2]),
        p=2,
    )

    assert [(rel.dependent_index, rel.indices, rel.coefficients) for rel in relations] == [
        (0, (0, 1, 2), (1, 1, 1)),
        (3, (1, 3), (1, 1)),
    ]
