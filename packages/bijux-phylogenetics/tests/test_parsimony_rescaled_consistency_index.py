from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.parsimony as parsimony_api
from bijux_phylogenetics.parsimony import (
    consistency_index,
    rescaled_consistency_index,
    retention_index,
    write_parsimony_rescaled_consistency_artifacts,
)
from bijux_phylogenetics.runtime.errors import ParsimonyAnalysisError

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_parsimony_gateway_exports_rescaled_consistency_surface() -> None:
    assert parsimony_api.rescaled_consistency_index is rescaled_consistency_index
    assert (
        parsimony_api.write_parsimony_rescaled_consistency_artifacts
        is write_parsimony_rescaled_consistency_artifacts
    )


def test_rescaled_consistency_index_matches_ci_times_ri_from_tested_surfaces() -> None:
    rc_report = rescaled_consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("rescaled_consistency_index_matrix.tsv"),
        method="fitch",
    )
    ci_report = consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("rescaled_consistency_index_matrix.tsv"),
        method="fitch",
    )
    ri_report = retention_index(
        fixture("fitch_tree.nwk"),
        fixture("rescaled_consistency_index_matrix.tsv"),
        method="fitch",
    )

    assert rc_report.algorithm == "parsimony-rescaled-consistency-index"
    assert rc_report.method == "fitch"
    assert rc_report.ci == ci_report.consistency_index == 0.75
    assert rc_report.ri == ri_report.retention_index == 0.5
    assert rc_report.rc == 0.375
    assert rc_report.rc == rc_report.ci * rc_report.ri
    assert [
        (row.character_id, row.ci, row.ri, row.rc, row.undefined_reason)
        for row in rc_report.character_rows
    ] == [
        (
            "char01_constant",
            None,
            None,
            None,
            "constant_character|zero_range_character",
        ),
        ("char02_singleton", 1.0, None, None, "zero_range_character"),
        ("char03_split", 1.0, 1.0, 1.0, None),
        ("char04_homoplastic", 0.5, 0.0, 0.0, None),
    ]


def test_rescaled_consistency_index_reports_undefined_aggregate_when_ci_and_ri_are_undefined() -> (
    None
):
    report = rescaled_consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("rescaled_consistency_index_constant_matrix.tsv"),
        method="fitch",
    )

    assert report.ci is None
    assert report.ri is None
    assert report.rc is None
    assert (
        report.undefined_reason
        == "no_variable_characters|no_defined_retention_characters"
    )


def test_rescaled_consistency_index_rejects_methods_without_common_ci_ri_support() -> (
    None
):
    with pytest.raises(ParsimonyAnalysisError) as error_info:
        rescaled_consistency_index(
            fixture("fitch_tree.nwk"),
            fixture("wagner_ordinal_matrix.tsv"),
            method="wagner",
        )

    assert error_info.value.code == "parsimony_rescaled_consistency_method_unsupported"


def test_write_parsimony_rescaled_consistency_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    report = rescaled_consistency_index(
        fixture("fitch_tree.nwk"),
        fixture("rescaled_consistency_index_matrix.tsv"),
        method="fitch",
    )

    outputs = write_parsimony_rescaled_consistency_artifacts(
        tmp_path / "rc-run", report
    )

    assert set(outputs) == {"indices_path", "run_json_path"}
    assert (
        outputs["indices_path"]
        .read_text(encoding="utf-8")
        .startswith("character_id\tci\tri\trc\tundefined_reason\n")
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "parsimony-rescaled-consistency-index"
    assert payload["ci"] == 0.75
    assert payload["ri"] == 0.5
    assert payload["rc"] == 0.375
