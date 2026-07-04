from __future__ import annotations

import pytest

from bijux_phylogenetics.parity import validate_reference_parity_examples


def test_validate_reference_parity_examples_supports_method_selection() -> None:
    report = validate_reference_parity_examples(methods=["pgls"])

    assert report.all_passed is True
    assert report.covered_methods == ["pgls"]
    assert report.case_count == 3
    assert {observation.method for observation in report.observations} == {"pgls"}


def test_validate_reference_parity_examples_rejects_unknown_methods() -> None:
    with pytest.raises(ValueError, match="unknown reference parity methods"):
        validate_reference_parity_examples(methods=["not-a-method"])
