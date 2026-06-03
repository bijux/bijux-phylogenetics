from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    retention_index,
    write_parsimony_retention_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_retention_index_surface() -> None:
    assert parsimony_api.retention_index is retention_index
    assert (
        parsimony_api.write_parsimony_retention_artifacts
        is write_parsimony_retention_artifacts
    )


def test_retention_index_matches_hand_computed_fitch_fixture() -> None:
    report = retention_index(
        fixture("fitch_tree.nwk"),
        fixture("retention_index_matrix.tsv"),
        method="fitch",
    )

    assert report.algorithm == "parsimony-retention-index"
    assert report.method == "fitch"
    assert report.included_character_count == 2
    assert report.excluded_character_count == 2
    assert report.minimum_possible_steps_total == 2.0
    assert report.maximum_possible_steps_total == 4.0
    assert report.observed_steps_total == 3.0
    assert report.retention_index == 0.5
    assert report.undefined_reason is None
    assert [
        (
            row.character_id,
            row.character_kind,
            row.minimum_possible_steps,
            row.maximum_possible_steps,
            row.observed_steps,
            row.retention_index,
            row.undefined_reason,
        )
        for row in report.character_rows
    ] == [
        ("char01_constant", "constant", 0.0, 0.0, 0.0, None, "zero_range_character"),
        (
            "char02_singleton",
            "parsimony-uninformative",
            1.0,
            1.0,
            1.0,
            None,
            "zero_range_character",
        ),
        ("char03_split", "parsimony-informative", 1.0, 2.0, 1.0, 1.0, None),
        ("char04_homoplastic", "parsimony-informative", 1.0, 2.0, 2.0, 0.0, None),
    ]


def test_retention_index_reports_undefined_aggregate_for_all_zero_range_matrix() -> (
    None
):
    report = retention_index(
        fixture("fitch_tree.nwk"),
        fixture("retention_index_constant_matrix.tsv"),
        method="fitch",
    )

    assert report.included_character_count == 0
    assert report.excluded_character_count == 2
    assert report.minimum_possible_steps_total == 0.0
    assert report.maximum_possible_steps_total == 0.0
    assert report.observed_steps_total == 0.0
    assert report.retention_index is None
    assert report.undefined_reason == "no_defined_retention_characters"


def test_retention_index_rejects_methods_without_owned_maximum_step_contracts() -> None:
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        retention_index(
            fixture("fitch_tree.nwk"),
            fixture("wagner_ordinal_matrix.tsv"),
            method="wagner",
        )

    assert error_info.value.code == "parsimony_retention_index_method_unsupported"


def test_write_parsimony_retention_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = retention_index(
        fixture("fitch_tree.nwk"),
        fixture("retention_index_matrix.tsv"),
        method="fitch",
    )

    outputs = write_parsimony_retention_artifacts(tmp_path / "retention-run", report)

    assert set(outputs) == {"indices_path", "run_json_path"}
    assert (
        outputs["indices_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "character_id\tcharacter_kind\tobserved_states\tminimum_possible_steps\tmaximum_possible_steps\tobserved_steps\tretention_index\tundefined_reason\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-retention-index"
    assert payload["retention_index"] == 0.5
