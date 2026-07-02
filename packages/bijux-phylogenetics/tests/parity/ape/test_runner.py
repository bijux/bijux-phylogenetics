from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.parity import run_ape_parity_cases
from tests.support.fake_reference_parity import fake_ape_rscript


@pytest.mark.slow
def test_run_ape_parity_cases_passes_representative_shards_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")

    report = run_ape_parity_cases(
        shard_ids=[
            "load_tree+extract_tree_clades",
            "extract_tree_clade_by_node_id",
            "find_tree_mrca",
            "assess_tree_monophyly",
            "summarize_brownian_covariance",
            "compute_tree_branching_times",
            "compute_diversification_gamma_statistic",
            "translate_coding_alignment",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is True
    assert report.case_count == 45
    assert report.passed_case_count == 45
    assert report.failed_case_count == 0
    assert report.skipped_case_count == 0
    assert [row.function_name for row in report.summary_rows] == [
        "ape::branching.times",
        "ape::extract.clade",
        "ape::gammaStat",
        "ape::getMRCA",
        "ape::is.monophyletic",
        "ape::read.tree",
        "ape::trans",
        "ape::vcv.phylo",
    ]
    assert all(observation.r_version == "4.6.0" for observation in report.observations)
    assert all(
        observation.ape_version == "5.0.0" for observation in report.observations
    )
    assert all(
        observation.reproducible_artifact_root is None
        for observation in report.observations
    )
    internal_label_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "read-tree-internal-node-labels"
    )
    assert internal_label_case.reference_summary is not None
    assert internal_label_case.reference_summary["tree_count"] == 1
    quoted_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "read-tree-quoted-taxon-labels"
    )
    assert quoted_case.reference_summary is not None
    assert set(quoted_case.reference_summary["tip_labels"]) == {
        "A.B-1",
        "Homo sapiens",
        "Mus musculus",
    }
    extract_clade_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "extract-clade-mammals"
    )
    assert extract_clade_case.reference_summary is not None
    assert extract_clade_case.reference_summary["requested_node_id"] == 6
    assert extract_clade_case.reference_summary["matched_node_id"] == 6
    assert extract_clade_case.reference_summary["matched_node_name"] == "Mammals"
    get_mrca_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "get-mrca-balanced-duplicate-request"
    )
    assert get_mrca_case.reference_summary is not None
    assert get_mrca_case.reference_summary["duplicate_requested_taxa"] == ["A"]
    assert get_mrca_case.reference_summary["matched_node_id"] == 6
    monophyly_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "is-monophyletic-rooted-three-tip-reroot-true"
    )
    assert monophyly_case.reference_summary is not None
    assert monophyly_case.reference_summary["monophyletic"] is True
    assert monophyly_case.reference_summary["complementary_clade_used"] is True
    assert monophyly_case.reference_summary["matched_extra_taxa"] == ["D"]
    vcv_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "vcv-zero-branch-singular"
    )
    assert vcv_case.reference_summary is not None
    assert vcv_case.reference_summary["singular"] is True
    assert vcv_case.reference_summary["near_singular"] is True
    assert vcv_case.reference_summary["matrix_dimension"] == 4
    branching_times_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "branching-times-zero-internal-branch"
    )
    assert branching_times_case.reference_summary is not None
    assert branching_times_case.reference_summary["internal_node_count"] == 3
    assert branching_times_case.reference_summary["root_age"] == pytest.approx(0.3)
    assert branching_times_case.reference_summary["zero_branch_length_count"] == 1
    gamma_stat_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "gamma-stat-zero-internal-branch"
    )
    assert gamma_stat_case.reference_summary is not None
    assert gamma_stat_case.reference_summary["tip_count"] == 4
    assert gamma_stat_case.reference_summary["bifurcating"] is True
    assert gamma_stat_case.reference_summary["gamma_statistic"] == pytest.approx(
        -0.979795897113271,
        abs=1e-12,
    )
    translation_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "dna-translation-terminal-stop"
    )
    assert translation_case.reference_summary is not None
    assert translation_case.reference_summary["stop_codon_count"] == 1
    frame_error_case = next(
        observation
        for observation in report.observations
        if observation.case_id == "dna-translation-frame-error-truncation"
    )
    assert frame_error_case.reference_summary is not None
    assert frame_error_case.reference_summary["dropped_trailing_nucleotide_count"] == 2
    assert frame_error_case.reference_summary["warning_count"] == 1


def test_run_ape_parity_cases_records_branch_length_failure_for_tree_structure_mismatch(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(
        tmp_path / "fake-ape-rscript",
        normalized_tree_overrides={
            "read-tree-balanced-rooted-ultrametric": "(A:0.2,B:0.2,(C:0.1,D:0.1):0.1);\n"
        },
    )

    report = run_ape_parity_cases(
        case_ids=["read-tree-balanced-rooted-ultrametric"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is False
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason is not None
    assert "branch lengths differ" in observation.mismatch_reason
    assert observation.reproducible_artifact_root is not None
    artifact_root = observation.reproducible_artifact_root
    assert artifact_root.exists()
    comparison_payload = json.loads(
        (artifact_root / "comparison.json").read_text(encoding="utf-8")
    )
    assert "branch lengths differ" in comparison_payload["mismatch_reason"]
    observed_summary = json.loads(
        (artifact_root / "reference-summary.observed.json").read_text(encoding="utf-8")
    )
    assert observed_summary["rooted"] is True
    assert (artifact_root / "bijux-normalized.txt").exists()


def test_run_ape_parity_cases_records_covariance_tables_for_vcv_summary_mismatch(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(
        tmp_path / "fake-ape-rscript",
        summary_overrides={"matrix_rank": 99},
    )

    report = run_ape_parity_cases(
        case_ids=["vcv-rooted-ultrametric"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is False
    observation = report.observations[0]
    assert observation.status == "failed"
    assert observation.mismatch_reason == "summary_mismatch"
    assert observation.reproducible_artifact_root is not None
    artifact_root = observation.reproducible_artifact_root
    assert (artifact_root / "reference-covariance-matrix.tsv").exists()
    assert (artifact_root / "reference-covariance-long.tsv").exists()
    assert (artifact_root / "reference-rows.observed.tsv").exists()
    assert (artifact_root / "bijux-rows.tsv").exists()


def test_run_ape_parity_cases_passes_expected_rooting_errors_against_fake_reference_runner(
    tmp_path: Path,
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")

    report = run_ape_parity_cases(
        case_ids=[
            "root-tree-missing-outgroup",
            "root-tree-non-monophyletic-outgroup",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "ape-parity-failures",
    )

    assert report.all_passed is True
    assert [observation.status for observation in report.observations] == [
        "passed",
        "passed",
    ]


def test_run_ape_parity_cases_marks_missing_rscript_as_skipped(tmp_path: Path) -> None:
    report = run_ape_parity_cases(
        case_ids=["read-tree-balanced-rooted-ultrametric"],
        rscript_executable=str(tmp_path / "missing-rscript"),
        failure_root=tmp_path / "ape-parity-failures",
    )

    observation = report.observations[0]
    assert observation.status == "skipped"
    assert observation.mismatch_reason == "rscript_unavailable"
    assert observation.reproducible_artifact_root is not None
    assert observation.reproducible_artifact_root.exists()
