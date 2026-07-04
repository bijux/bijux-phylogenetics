from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    consistency_index,
    write_parsimony_consistency_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_consistency_index_surface() -> None:
    assert parsimony_api.consistency_index is consistency_index
    assert (
        parsimony_api.write_parsimony_consistency_artifacts
        is write_parsimony_consistency_artifacts
    )


def test_consistency_index_matches_hand_computed_fitch_fixture() -> None:
    report = consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("consistency_index_matrix.tsv"),
        method="fitch",
    )

    assert report.algorithm == "parsimony-consistency-index"
    assert report.method == "fitch"
    assert report.included_character_count == 3
    assert report.excluded_character_count == 1
    assert report.minimum_possible_steps_total == 3.0
    assert report.observed_steps_total == 4.0
    assert report.consistency_index == 0.75
    assert report.undefined_reason is None
    assert [
        (
            row.character_id,
            row.character_kind,
            row.minimum_possible_steps,
            row.observed_steps,
            row.consistency_index,
            row.undefined_reason,
        )
        for row in report.character_rows
    ] == [
        ("char01_constant", "constant", 0.0, 0.0, None, "constant_character"),
        ("char02_singleton", "parsimony-uninformative", 1.0, 1.0, 1.0, None),
        ("char03_split", "parsimony-informative", 1.0, 1.0, 1.0, None),
        ("char04_homoplastic", "parsimony-informative", 1.0, 2.0, 0.5, None),
    ]


def test_consistency_index_reports_undefined_aggregate_for_all_constant_matrix() -> (
    None
):
    report = consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("consistency_index_constant_matrix.tsv"),
        method="fitch",
    )

    assert report.included_character_count == 0
    assert report.excluded_character_count == 2
    assert report.minimum_possible_steps_total == 0.0
    assert report.observed_steps_total == 0.0
    assert report.consistency_index is None
    assert report.undefined_reason == "no_variable_characters"


def test_consistency_index_supports_wagner_minimum_possible_step_policy() -> None:
    report = consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("wagner_ordinal_matrix.tsv"),
        method="wagner",
    )

    assert report.minimum_possible_steps_total == 6.0
    assert report.observed_steps_total == 6.0
    assert report.consistency_index == 1.0
    assert [
        (row.character_id, row.minimum_possible_steps, row.observed_steps)
        for row in report.character_rows
    ] == [
        ("char01_gradient", 3.0, 3.0),
        ("char02_distance_gap", 3.0, 3.0),
    ]


def test_consistency_index_rejects_sankoff_methods_until_step_matrix_minima_are_owned() -> (
    None
):
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        consistency_index(
            fixture("sankoff_tree_5_taxa.nwk"),
            fixture("sankoff_character_matrix.tsv"),
            method="sankoff",
        )

    assert error_info.value.code == "parsimony_consistency_index_method_unsupported"


def test_write_parsimony_consistency_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("consistency_index_matrix.tsv"),
        method="fitch",
    )

    outputs = write_parsimony_consistency_artifacts(
        tmp_path / "consistency-run", report
    )

    assert set(outputs) == {"indices_path", "run_json_path"}
    assert (
        outputs["indices_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tcharacter_kind\tobserved_states\tminimum_possible_steps\tobserved_steps\tconsistency_index\tundefined_reason\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-consistency-index"
    assert payload["consistency_index"] == 0.75
