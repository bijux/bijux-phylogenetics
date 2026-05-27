from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    load_sankoff_cost_matrix,
    validate_sankoff_cost_matrix,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_sankoff_validation_surface() -> None:
    assert parsimony_api.validate_sankoff_cost_matrix is validate_sankoff_cost_matrix


def test_validate_sankoff_cost_matrix_reports_nonzero_diagonal_warnings() -> None:
    report = validate_sankoff_cost_matrix(
        fixture("sankoff_diagonal_nonzero_cost_matrix.tsv"),
        observed_states=["red", "green", "blue"],
    )

    assert [warning.code for warning in report.validation_warnings] == [
        "parsimony_cost_matrix_diagonal_nonzero"
    ]
    assert report.validation_warnings[0].details["states"] == ["red"]


def test_validate_sankoff_cost_matrix_reports_unused_state_warnings() -> None:
    report = validate_sankoff_cost_matrix(
        fixture("sankoff_unused_state_cost_matrix.tsv"),
        observed_states=["red", "green", "blue"],
    )

    assert [warning.code for warning in report.validation_warnings] == [
        "parsimony_cost_matrix_unused_states"
    ]
    assert report.validation_warnings[0].details["unused_states"] == ["yellow"]


def test_validate_sankoff_cost_matrix_rejects_asymmetry_when_forbidden() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        validate_sankoff_cost_matrix(
            fixture("sankoff_asymmetric_cost_matrix.tsv"),
            observed_states=["red", "green", "blue"],
        )

    assert error_info.value.code == "parsimony_cost_matrix_asymmetric"
    assert error_info.value.details["asymmetric_pairs"] == [
        {
            "from_state": "red",
            "to_state": "green",
            "forward_cost": 1.0,
            "reverse_cost": 2.0,
        }
    ]


def test_validate_sankoff_cost_matrix_allows_asymmetry_when_requested() -> None:
    report = validate_sankoff_cost_matrix(
        fixture("sankoff_asymmetric_cost_matrix.tsv"),
        observed_states=["red", "green", "blue"],
        allow_asymmetric_costs=True,
    )

    assert report.states == ["red", "green", "blue"]
    assert report.validation_warnings == []


def test_load_sankoff_cost_matrix_preserves_validation_warning_surface() -> None:
    report = load_sankoff_cost_matrix(
        fixture("sankoff_diagonal_nonzero_cost_matrix.tsv"),
        observed_states=["red", "green", "blue"],
    )

    assert [warning.code for warning in report.validation_warnings] == [
        "parsimony_cost_matrix_diagonal_nonzero"
    ]
