from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    load_sankoff_cost_matrix,
    score_sankoff,
    validate_sankoff_cost_matrix,
    write_sankoff_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_sankoff_surface() -> None:
    assert parsimony_api.load_sankoff_cost_matrix is load_sankoff_cost_matrix
    assert parsimony_api.validate_sankoff_cost_matrix is validate_sankoff_cost_matrix
    assert parsimony_api.score_sankoff is score_sankoff
    assert parsimony_api.write_sankoff_artifacts is write_sankoff_artifacts


def test_score_sankoff_matches_hand_computed_five_taxon_fixture() -> None:
    report = score_sankoff(
        fixture("sankoff_tree_5_taxa.nwk"),
        fixture("sankoff_character_matrix.tsv"),
        fixture("sankoff_cost_matrix.tsv"),
    )

    assert report.total_cost == 3.0
    assert [(row.character_id, row.minimum_cost) for row in report.step_rows] == [
        ("char01_color", 3.0),
    ]
    node_costs = {
        (row.character_id, row.node, row.state): (row.cost, row.is_optimal_state)
        for row in report.node_cost_rows
    }
    assert node_costs[("char01_color", "A|B", "red")] == (1.0, True)
    assert node_costs[("char01_color", "A|B", "green")] == (3.0, False)
    assert node_costs[("char01_color", "A|B", "blue")] == (1.0, True)
    assert node_costs[("char01_color", "A|B|C", "red")] == (3.0, False)
    assert node_costs[("char01_color", "A|B|C", "green")] == (2.0, True)
    assert node_costs[("char01_color", "A|B|C", "blue")] == (2.0, True)
    assert node_costs[("char01_color", "D|E", "red")] == (3.0, False)
    assert node_costs[("char01_color", "D|E", "green")] == (1.0, True)
    assert node_costs[("char01_color", "D|E", "blue")] == (1.0, True)
    assert node_costs[("char01_color", "A|B|C|D|E", "red")] == (5.0, False)
    assert node_costs[("char01_color", "A|B|C|D|E", "green")] == (3.0, True)
    assert node_costs[("char01_color", "A|B|C|D|E", "blue")] == (3.0, True)

    selections = {
        (row.character_id, row.node): (row.optimal_states, row.tie_states)
        for row in report.node_selection_rows
    }
    assert selections[("char01_color", "A|B")] == (
        ["red", "blue"],
        ["red", "blue"],
    )
    assert selections[("char01_color", "A|B|C")] == (
        ["green", "blue"],
        ["green", "blue"],
    )
    assert selections[("char01_color", "D|E")] == (
        ["green", "blue"],
        ["green", "blue"],
    )
    assert selections[("char01_color", "A|B|C|D|E")] == (
        ["green", "blue"],
        ["green", "blue"],
    )


def test_load_sankoff_cost_matrix_rejects_negative_costs() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        load_sankoff_cost_matrix(fixture("sankoff_negative_cost_matrix.tsv"))

    assert error_info.value.code == "parsimony_cost_matrix_negative_cost"
    assert error_info.value.details["row_state"] == "red"
    assert error_info.value.details["column_state"] == "blue"


def test_load_sankoff_cost_matrix_rejects_non_square_matrices() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        load_sankoff_cost_matrix(fixture("sankoff_non_square_cost_matrix.tsv"))

    assert error_info.value.code == "parsimony_cost_matrix_not_square"
    assert error_info.value.details["row_count"] == 2
    assert error_info.value.details["column_count"] == 3


def test_load_sankoff_cost_matrix_rejects_inconsistent_labels() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        load_sankoff_cost_matrix(fixture("sankoff_inconsistent_labels_cost_matrix.tsv"))

    assert error_info.value.code == "parsimony_cost_matrix_inconsistent_labels"
    assert error_info.value.details["column_labels"] == ["red", "green", "blue"]


def test_score_sankoff_rejects_character_states_missing_from_cost_matrix() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        score_sankoff(
            fixture("sankoff_tree_5_taxa.nwk"),
            fixture("sankoff_missing_state_character_matrix.tsv"),
            fixture("sankoff_cost_matrix.tsv"),
        )

    assert error_info.value.code == "parsimony_cost_matrix_missing_states"
    assert error_info.value.details["missing_states"] == ["yellow"]


def test_score_sankoff_rejects_asymmetric_cost_matrices_before_scoring() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        score_sankoff(
            fixture("sankoff_tree_5_taxa.nwk"),
            fixture("sankoff_character_matrix.tsv"),
            fixture("sankoff_asymmetric_cost_matrix.tsv"),
        )

    assert error_info.value.code == "parsimony_cost_matrix_asymmetric"


def test_score_sankoff_records_diagonal_and_unused_state_warnings() -> None:
    diagonal_report = score_sankoff(
        fixture("sankoff_tree_5_taxa.nwk"),
        fixture("sankoff_character_matrix.tsv"),
        fixture("sankoff_diagonal_nonzero_cost_matrix.tsv"),
    )
    unused_state_report = score_sankoff(
        fixture("sankoff_tree_5_taxa.nwk"),
        fixture("sankoff_character_matrix.tsv"),
        fixture("sankoff_unused_state_cost_matrix.tsv"),
    )

    assert [warning.code for warning in diagonal_report.validation_warnings] == [
        "parsimony_cost_matrix_diagonal_nonzero"
    ]
    assert [warning.code for warning in unused_state_report.validation_warnings] == [
        "parsimony_cost_matrix_unused_states"
    ]


def test_write_sankoff_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = score_sankoff(
        fixture("sankoff_tree_5_taxa.nwk"),
        fixture("sankoff_character_matrix.tsv"),
        fixture("sankoff_cost_matrix.tsv"),
    )

    outputs = write_sankoff_artifacts(tmp_path / "sankoff-run", report)

    assert set(outputs) == {
        "steps_path",
        "node_costs_path",
        "node_selection_path",
        "run_json_path",
    }
    assert (
        outputs["steps_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tminimum_cost\tobserved_states\tmatrix_states\tcharacter_weight\tweighted_score\n"
        )
    )
    assert (
        outputs["node_costs_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tnode\tnode_name\tdescendant_taxa\tstate\tcost\tis_optimal_state\n"
        )
    )
    assert (
        outputs["node_selection_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tnode\tnode_name\tdescendant_taxa\toptimal_states\ttie_states\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "sankoff"
    assert payload["total_cost"] == 3.0
    assert payload["total_weighted_score"] == 3.0
    assert payload["validation_warnings"] == []
