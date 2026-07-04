from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    load_parsimony_character_weights,
    tree_length,
    write_parsimony_tree_length_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_tree_length_surface() -> None:
    assert parsimony_api.tree_length is tree_length
    assert (
        parsimony_api.load_parsimony_character_weights
        is load_parsimony_character_weights
    )
    assert (
        parsimony_api.write_parsimony_tree_length_artifacts
        is write_parsimony_tree_length_artifacts
    )


def test_tree_length_matches_hand_computed_unweighted_fitch_fixture() -> None:
    report = tree_length(
        fixture("fitch_tree.nwk"),
        fixture("fitch_binary_matrix.tsv"),
        method="fitch",
    )

    assert report.algorithm == "parsimony-tree-length"
    assert report.method == "fitch"
    assert report.raw_total_score == 2.0
    assert report.total_score == 2.0
    assert [
        (
            row.character_id,
            row.raw_score,
            row.character_weight,
            row.weighted_score,
        )
        for row in report.step_rows
    ] == [
        ("char01_split", 1.0, 1.0, 1.0),
        ("char02_left_conflict", 1.0, 1.0, 1.0),
    ]


def test_tree_length_matches_hand_computed_weighted_fitch_fixture() -> None:
    report = tree_length(
        fixture("fitch_tree.nwk"),
        fixture("fitch_binary_matrix.tsv"),
        method="fitch",
        character_weights=fixture("fitch_character_weights.tsv"),
    )

    assert report.raw_total_score == 2.0
    assert report.total_score == 3.5
    assert [
        (
            row.character_id,
            row.raw_score,
            row.character_weight,
            row.weighted_score,
        )
        for row in report.step_rows
    ] == [
        ("char01_split", 1.0, 1.5, 1.5),
        ("char02_left_conflict", 1.0, 2.0, 2.0),
    ]


def test_tree_length_matches_hand_computed_weighted_wagner_fixture() -> None:
    report = tree_length(
        fixture("fitch_tree.nwk"),
        fixture("wagner_ordinal_matrix.tsv"),
        method="wagner",
        character_weights=fixture("wagner_character_weights.tsv"),
    )

    assert report.raw_total_score == 6.0
    assert report.total_score == 7.5
    assert [
        (
            row.character_id,
            row.raw_score,
            row.character_weight,
            row.weighted_score,
        )
        for row in report.step_rows
    ] == [
        ("char01_gradient", 3.0, 0.5, 1.5),
        ("char02_distance_gap", 3.0, 2.0, 6.0),
    ]


def test_tree_length_requires_sankoff_cost_matrix() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        tree_length(
            fixture("sankoff_tree_5_taxa.nwk"),
            fixture("sankoff_character_matrix.tsv"),
            method="sankoff",
        )

    assert error_info.value.code == "parsimony_tree_length_cost_matrix_required"


def test_tree_length_rejects_asymmetric_sankoff_cost_matrix_by_default() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        tree_length(
            fixture("sankoff_tree_5_taxa.nwk"),
            fixture("sankoff_character_matrix.tsv"),
            method="sankoff",
            cost_matrix=fixture("sankoff_asymmetric_cost_matrix.tsv"),
        )

    assert error_info.value.code == "parsimony_cost_matrix_asymmetric"


def test_tree_length_allows_asymmetric_sankoff_cost_matrix_when_requested() -> None:
    report = tree_length(
        fixture("sankoff_tree_5_taxa.nwk"),
        fixture("sankoff_character_matrix.tsv"),
        method="sankoff",
        cost_matrix=fixture("sankoff_asymmetric_cost_matrix.tsv"),
        allow_asymmetric_costs=True,
    )

    assert report.method == "sankoff"
    assert report.raw_total_score >= 0.0


def test_tree_length_rejects_missing_character_weights() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        tree_length(
            fixture("fitch_tree.nwk"),
            fixture("fitch_binary_matrix.tsv"),
            method="fitch",
            character_weights=fixture("parsimony_missing_character_weight.tsv"),
        )

    assert error_info.value.code == "parsimony_character_weight_missing_character"
    assert error_info.value.details["missing_characters"] == ["char02_left_conflict"]


def test_load_parsimony_character_weights_rejects_nonnumeric_values() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        load_parsimony_character_weights(
            fixture("parsimony_invalid_character_weight.tsv")
        )

    assert error_info.value.code == "parsimony_character_weight_invalid_value"
    assert error_info.value.details["character_id"] == "char01_split"


def test_write_parsimony_tree_length_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = tree_length(
        fixture("fitch_tree.nwk"),
        fixture("fitch_binary_matrix.tsv"),
        method="fitch",
        character_weights=fixture("fitch_character_weights.tsv"),
    )

    outputs = write_parsimony_tree_length_artifacts(
        tmp_path / "tree-length-run", report
    )

    assert set(outputs) == {"scores_path", "run_json_path"}
    assert (
        outputs["scores_path"]
        .read_text(encoding="utf-8")
        .startswith("character_id\traw_score\tcharacter_weight\tweighted_score\n")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-tree-length"
    assert payload["method"] == "fitch"
    assert payload["raw_total_score"] == 2.0
    assert payload["total_score"] == 3.5
