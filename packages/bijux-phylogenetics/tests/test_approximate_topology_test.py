from __future__ import annotations

import json
import math
from pathlib import Path

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_nucleotide_approximate_topology_test_from_alignment,
    validate_approximate_topology_test_replicate_count,
    write_approximate_topology_test_artifacts,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def write_candidate_tree_set(path: Path, trees: list[str]) -> Path:
    path.write_text("".join(f"{tree}\n" for tree in trees), encoding="utf-8")
    return path


def test_package_likelihood_gateway_exports_approximate_topology_test_surface() -> None:
    assert (
        likelihood_api.evaluate_nucleotide_approximate_topology_test_from_alignment
        is evaluate_nucleotide_approximate_topology_test_from_alignment
    )
    assert (
        likelihood_api.write_approximate_topology_test_artifacts
        is write_approximate_topology_test_artifacts
    )
    assert (
        likelihood_api.validate_approximate_topology_test_replicate_count
        is validate_approximate_topology_test_replicate_count
    )


def test_approximate_topology_test_reports_observed_delta_distribution_and_caution(
    tmp_path: Path,
) -> None:
    tree_set_path = write_candidate_tree_set(
        tmp_path / "candidate-trees.nwk",
        [
            "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);",
        ],
    )

    report = evaluate_nucleotide_approximate_topology_test_from_alignment(
        tree_set_path,
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
        resampling_replicate_count=8,
        resampling_seed=7,
    )

    assert report.algorithm == "nucleotide-approximate-topology-test"
    assert report.model_name == "JC69"
    assert report.tree_count == 2
    assert report.site_count == 10
    assert report.pattern_count == 6
    assert report.resampling_method == "site-resampling-with-replacement"
    assert report.resampling_replicate_count == 8
    assert report.resampling_seed == 7
    assert (
        report.caution_label
        == "site-resampled p-like statistics are approximate ranking aids and must not be interpreted as AU/SH-style p-values"
    )
    assert report.observed_best_tree_id == "candidate-tree-2"
    assert [row.candidate_tree_id for row in report.summary_rows] == [
        "candidate-tree-1",
        "candidate-tree-2",
    ]
    assert len(report.resampling_rows) == 16

    loser = next(row for row in report.summary_rows if not row.observed_best_tree)
    loser_distribution = [
        row
        for row in report.resampling_rows
        if row.candidate_tree_id == loser.candidate_tree_id
    ]
    winner = next(row for row in report.summary_rows if row.observed_best_tree)
    expected_p_like = sum(
        1.0
        for row in loser_distribution
        if row.candidate_matches_or_beats_observed_best
    ) / float(report.resampling_replicate_count)

    assert math.isclose(
        loser.observed_delta_log_likelihood,
        winner.observed_log_likelihood - loser.observed_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        loser.p_like_statistic, expected_p_like, rel_tol=0.0, abs_tol=1e-12
    )
    assert loser.p_like_statistic == 0.0
    assert loser.caution_label == report.caution_label


def test_write_approximate_topology_test_artifacts_materializes_governed_output_family(
    tmp_path: Path,
) -> None:
    tree_set_path = write_candidate_tree_set(
        tmp_path / "candidate-trees.nwk",
        [
            "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);",
            "(((A:0.1,B:0.1):0.1,C:0.1):0.1,D:0.1);",
        ],
    )
    report = evaluate_nucleotide_approximate_topology_test_from_alignment(
        tree_set_path,
        fixture("alignments", "jc69_site_pattern_alignment.fasta"),
        model_name="jc69",
        resampling_replicate_count=8,
        resampling_seed=7,
    )

    outputs = write_approximate_topology_test_artifacts(
        tmp_path / "approximate-topology-test-run",
        report,
    )

    assert set(outputs) == {
        "candidate_tree_path",
        "summary_path",
        "resampling_path",
        "run_json_path",
    }
    assert (
        outputs["summary_path"]
        .read_text(encoding="utf-8")
        .startswith(
            "candidate_tree_id\tcandidate_tree_label\tobserved_log_likelihood\tobserved_delta_log_likelihood\tobserved_best_tree\tresampling_win_count\tresampling_frequency\tp_like_statistic\tresampling_mean_delta_log_likelihood\tresampling_min_delta_log_likelihood\tresampling_max_delta_log_likelihood\tcaution_label\ttree_newick\n"
        )
    )
    payload = json.loads(outputs["run_json_path"].read_text(encoding="utf-8"))
    assert payload["algorithm"] == "nucleotide-approximate-topology-test"
    assert payload["observed_best_tree_id"] == "candidate-tree-2"
    assert payload["resampling_replicate_count"] == 8
    assert len(payload["summary_rows"]) == 2
    assert len(payload["resampling_rows"]) == 16


def test_validate_approximate_topology_test_rejects_zero_replicates() -> None:
    try:
        validate_approximate_topology_test_replicate_count(0)
    except ValueError as error:
        assert str(error) == "resampling_replicate_count must be at least one"
    else:
        raise AssertionError(
            "approximate topology test must require at least one replicate"
        )
