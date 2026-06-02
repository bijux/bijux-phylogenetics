from __future__ import annotations

import json
from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    summarize_likelihood_wrapper_correspondence,
    write_likelihood_wrapper_correspondence_artifacts,
)

pytestmark = pytest.mark.slow


def test_public_likelihood_gateway_exports_wrapper_correspondence_surface() -> None:
    assert (
        likelihood_api.summarize_likelihood_wrapper_correspondence
        is summarize_likelihood_wrapper_correspondence
    )
    assert (
        likelihood_api.write_likelihood_wrapper_correspondence_artifacts
        is write_likelihood_wrapper_correspondence_artifacts
    )


def test_likelihood_wrapper_correspondence_classifies_native_and_wrapper_cases(
    tmp_path: Path,
) -> None:
    report = summarize_likelihood_wrapper_correspondence()

    assert report.case_count == 4
    assert report.supported_case_count == 3
    assert report.exact_match_case_count == 1
    assert report.tolerance_match_case_count == 1
    assert report.expected_model_assumption_difference_case_count == 1
    assert report.unsupported_case_count == 1
    assert report.native_bug_case_count == 0
    assert report.blocking_case_count == 0
    assert report.all_supported_cases_clear is True

    observations = {row.case_id: row for row in report.observations}
    exact = observations["iqtree-small-topology-reference"]
    assert exact.status == "exact-match"
    assert exact.observed_output["topology_equal"] is True
    assert exact.observed_output["robinson_foulds_distance"] == 0

    tolerance = observations["iqtree-small-rounded-branch-reference"]
    assert tolerance.status == "tolerance-match"
    assert tolerance.observed_output["topology_equal"] is True
    assert tolerance.observed_output["branch_score_distance"] <= tolerance.tolerance

    fasttree = observations["fasttree-small-approximate-reference"]
    assert fasttree.status == "expected-model-assumption-difference"
    assert fasttree.observed_output["topology_equal"] is True
    assert fasttree.expected_output["support_labels_present"] is True

    unsupported = observations["raxml-style-small-reference"]
    assert unsupported.status == "unsupported-case"
    assert unsupported.supported is False
    assert unsupported.blocking is False

    outputs = write_likelihood_wrapper_correspondence_artifacts(
        tmp_path / "likelihood-wrapper-correspondence",
        report,
    )
    assert set(outputs) == {
        "summary_path",
        "observations_path",
        "run_json_path",
    }
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith("status\tcase_count\tblocking_case_count\tcase_ids\n")
    )
    assert (
        outputs["observations_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "case_id\twrapper_engine\tnative_surface\twrapper_surface\tcomparison_policy\tstatus\tsupported\tblocking\ttolerance\trationale\tinput_fixtures\texpected_output\tobserved_output\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["case_count"] == 4
    assert payload["native_bug_case_count"] == 0
    assert payload["summary_rows"][0]["status"] == "exact-match"
    assert {row["case_id"] for row in payload["observations"]} == set(observations)
