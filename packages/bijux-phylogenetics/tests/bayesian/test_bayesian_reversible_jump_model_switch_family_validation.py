from __future__ import annotations

import pytest

from bijux_phylogenetics.bayesian import (
    list_reversible_jump_model_switch_families,
    validate_reversible_jump_model_switch_family,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_reversible_jump_model_switch_family_validation_accepts_declared_family() -> (
    None
):
    assert list_reversible_jump_model_switch_families() == (
        "nucleotide-substitution-model",
    )
    assert (
        validate_reversible_jump_model_switch_family("nucleotide-substitution-model")
        == "nucleotide-substitution-model"
    )


def test_reversible_jump_model_switch_family_validation_rejects_unknown_family() -> (
    None
):
    with pytest.raises(
        PhylogeneticsError,
        match="requires one declared model family",
    ):
        validate_reversible_jump_model_switch_family("continuous-trait-model")
