from __future__ import annotations

from bijux_phylogenetics.bayesian.wrapper_correspondence import (
    BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES,
    BayesianWrapperCorrespondenceReport,
    summarize_bayesian_wrapper_correspondence,
)


def test_summarize_bayesian_wrapper_correspondence_reports_governed_statuses() -> None:
    report = summarize_bayesian_wrapper_correspondence()

    assert isinstance(report, BayesianWrapperCorrespondenceReport)
    assert report.case_count == 6
    assert report.supported_case_count == 5
    assert report.exact_match_case_count == 3
    assert report.tolerance_match_case_count == 1
    assert report.expected_model_assumption_difference_case_count == 1
    assert report.unsupported_case_count == 1
    assert report.native_bug_case_count == 0
    assert report.blocking_case_count == 0
    assert report.all_supported_cases_clear is True
    assert [row.status for row in report.summary_rows] == list(
        BAYESIAN_WRAPPER_CORRESPONDENCE_STATUSES
    )


def test_summarize_bayesian_wrapper_correspondence_separates_wrapper_outcomes() -> None:
    report = summarize_bayesian_wrapper_correspondence()
    by_case = {observation.case_id: observation for observation in report.observations}

    beast_log = by_case["beast-log-parameter-summaries-strict-yule"]
    beast_consensus = by_case["beast-consensus-tree-strict-yule"]
    beast_mcc = by_case["beast-maximum-clade-credibility-tree-strict-yule"]
    mrbayes_topology = by_case["mrbayes-consensus-topology-partitioned-analysis"]
    mrbayes_branch_lengths = by_case[
        "mrbayes-consensus-branch-length-semantics-partitioned-analysis"
    ]
    revbayes = by_case["revbayes-governed-posterior-corpus"]

    assert beast_log.status == "tolerance-match"
    assert beast_log.tolerance == 1e-9
    assert "currently matches the cached reference exactly" in beast_log.rationale
    assert (
        beast_log.observed_output["posterior.mean"]
        == beast_log.expected_output["posterior.mean"]
    )

    assert beast_consensus.status == "exact-match"
    assert (
        beast_consensus.observed_output["consensus_newick"]
        == beast_consensus.expected_output["consensus_newick"]
    )
    assert beast_consensus.observed_output["annotated_node_count"] == 2

    assert beast_mcc.status == "exact-match"
    assert beast_mcc.observed_output["selected_tree_index"] == 13
    assert (
        beast_mcc.observed_output["mcc_newick"]
        == beast_mcc.expected_output["mcc_newick"]
    )

    assert mrbayes_topology.status == "exact-match"
    assert mrbayes_topology.observed_output["rooted_robinson_foulds_distance"] == 0
    assert mrbayes_topology.observed_output["topology_equal"] is True

    assert mrbayes_branch_lengths.status == "expected-model-assumption-difference"
    assert (
        mrbayes_branch_lengths.observed_output["rooted_robinson_foulds_distance"] == 0
    )
    assert (
        mrbayes_branch_lengths.observed_output["same_topology_different_branch_lengths"]
        is True
    )
    assert float(mrbayes_branch_lengths.observed_output["branch_score_distance"]) > 0.0
    assert (
        mrbayes_branch_lengths.observed_output["native_consensus_newick"]
        != mrbayes_branch_lengths.observed_output["wrapper_consensus_newick"]
    )

    assert revbayes.status == "unsupported-case"
    assert revbayes.supported is False
    assert revbayes.input_fixtures == []
