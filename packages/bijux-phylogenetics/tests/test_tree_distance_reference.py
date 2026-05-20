from __future__ import annotations

from bijux_phylogenetics.compare.reference import (
    validate_tree_distance_reference_examples,
)


def test_validate_tree_distance_reference_examples_reports_passing_cases() -> None:
    report = validate_tree_distance_reference_examples()

    assert report.case_count == 27
    assert report.external_case_count == 25
    assert report.policy_case_count == 2
    assert report.all_passed is True
    assert {observation.category for observation in report.observations} == {
        "external-reference",
        "overlap-policy",
    }


def test_validate_tree_distance_reference_examples_covers_polytomy_and_overlap() -> (
    None
):
    report = validate_tree_distance_reference_examples()
    observations = {row.case_id: row for row in report.observations}

    branch_score = observations["branch-score-polytomy-overlap-pruned"]
    rooted_rf = observations["rf-rooted-polytomy-binary"]
    unrooted_rf = observations["rf-unrooted-polytomy-star"]

    assert branch_score.passed is True
    assert branch_score.observed_metrics["same_taxon_set"] is False
    assert branch_score.observed_metrics["branch_score_distance"] == 0.45825756949558405
    assert branch_score.observed_metrics["left_only_split_count"] == 1
    assert rooted_rf.observed_metrics["robinson_foulds_distance"] == 3
    assert rooted_rf.observed_metrics["normalized_robinson_foulds"] == 1.0
    assert unrooted_rf.observed_metrics["robinson_foulds_distance"] == 0
    assert unrooted_rf.observed_metrics["normalized_robinson_foulds"] == 0.0


def test_validate_tree_distance_reference_examples_keeps_overlap_policy_explicit() -> (
    None
):
    report = validate_tree_distance_reference_examples()
    observations = {row.case_id: row for row in report.observations}

    rf_policy = observations["rf-require-identical-overlap-policy"]
    branch_score_policy = observations["branch-score-require-identical-overlap-policy"]

    assert rf_policy.passed is True
    assert branch_score_policy.passed is True
    assert "identical taxon sets" in rf_policy.observed_metrics["error"]
    assert "identical taxon sets" in branch_score_policy.observed_metrics["error"]
