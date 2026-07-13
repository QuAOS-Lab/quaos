from __future__ import annotations

import numpy as np

from sympleq.core.symmetries.atomic_decomposition import atomic_block_decompose
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_types import SectorCostCertificate
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_unipotent_p2_generators import (
    canonical_W_block_p2,
    direct_sum_unipotent_p2_blocks,
)
from sympleq.core.symmetries.atomic_decomposition_helpers.atomic_verify import verify_global_basis
from sympleq.core.symmetries.modular_helpers import is_symplectic, mod_p


def _p2_payload(info: dict) -> dict:
    for inv in info.get("sector_invariants", []):
        data = getattr(inv, "data", {})
        if isinstance(data, dict) and "p2_unipotent" in data:
            return data["p2_unipotent"]
    raise AssertionError("No p=2 unipotent payload found")


def test_phase_c_sector_certificates_have_typed_and_legacy_views() -> None:
    F = canonical_W_block_p2(3)
    Sigma, B, info = atomic_block_decompose(F, 2)
    verify_global_basis(F, B, Sigma, 2)

    cert = info["cost_certificate"]
    assert isinstance(cert["attained"], int)
    assert cert["attained"] == info["Q_att"] == info["Q_opt"]
    assert cert["certified_minimal"] is True

    for inv in info.get("sector_invariants", []):
        data = getattr(inv, "data", {})
        assert data.get("status") == "OK"
        typed = data.get("sector_cost_certificate")
        legacy = data.get("cost_certificate")
        assert isinstance(typed, SectorCostCertificate)
        assert isinstance(legacy, dict)
        assert typed.as_dict() == legacy
        assert legacy["extraction_attained"] is True
        assert legacy["lower_bound_complete"] is True
        assert legacy["complete"] is True
        assert legacy["certified_minimal_sector"] is True


def test_phase_c_p2_unipotent_uses_quotient_normal_form_not_enumeration() -> None:
    F = direct_sum_unipotent_p2_blocks(
        [
            ("W", 3, None),
            ("V_beta", 4, 1),
            ("W_beta", 5, 1),
            ("V_beta", 2, 0),
        ]
    )
    F = mod_p(F, 2)
    assert is_symplectic(F, 2)

    Sigma, B, info = atomic_block_decompose(F, 2)
    verify_global_basis(F, B, Sigma, 2)

    payload = _p2_payload(info)
    assert payload["classification_complete"] is True
    assert payload["used_best_effort_sweep"] is False
    assert payload["candidate_enumeration_used"] is False
    assert payload["certified_extraction_policy"] == "quotient_top_normal_form_no_enumeration"
    assert sorted(int(h) for h in info["atomic_half_dims"]) == sorted([3, 2, 5, 1])
    assert info["Q_opt"] == 5
