from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import search_nucleotide_likelihood_nni_from_alignment
from bijux_phylogenetics.phylo.likelihood.nni_search import (
    validate_nucleotide_likelihood_nni_improvement_policy,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_nni_first_improvement_stops_after_first_accepted_candidate() -> (
    None
):
    best_report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        improvement_policy="best-improvement",
        upper_branch_length_bound=1.0,
    )
    first_report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        improvement_policy="first-improvement",
        upper_branch_length_bound=1.0,
    )

    assert first_report.improvement_policy == "first-improvement"
    assert best_report.improvement_policy == "best-improvement"
    assert first_report.accepted_move_count == best_report.accepted_move_count == 2
    assert first_report.final_tree_newick == best_report.final_tree_newick
    assert math.isclose(
        first_report.final_log_likelihood,
        best_report.final_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert len([row for row in best_report.candidate_rows if row.iteration == 1]) == 4
    assert len([row for row in first_report.candidate_rows if row.iteration == 1]) == 4
    assert len([row for row in best_report.candidate_rows if row.iteration == 2]) == 4
    assert len([row for row in first_report.candidate_rows if row.iteration == 2]) == 2
    assert len([row for row in best_report.candidate_rows if row.iteration == 3]) == 4
    assert len([row for row in first_report.candidate_rows if row.iteration == 3]) == 4
    assert first_report.candidate_rows[-1].iteration == 3


def test_likelihood_nni_first_improvement_preserves_deterministic_candidate_order() -> (
    None
):
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        improvement_policy="first-improvement",
        upper_branch_length_bound=1.0,
    )

    selected_rows = [row for row in report.candidate_rows if row.selected_best_move]

    assert [row.candidate_order for row in selected_rows] == [4, 2]
    assert [row.iteration for row in selected_rows] == [1, 2]
    assert selected_rows[0].pivot_branch_id == "A|C"
    assert selected_rows[1].pivot_branch_id == "A|B|C"
    assert all(
        row.iteration != 2 or row.candidate_order <= 2
        for row in report.candidate_rows
    )
    assert report.trace_rows[-1].stopping_reason == "no-improving-neighbor"


def test_likelihood_nni_improvement_policy_validation_rejects_unknown_name() -> None:
    try:
        validate_nucleotide_likelihood_nni_improvement_policy("greedy")
    except ValueError as error:
        assert "improvement_policy must be one of" in str(error)
    else:
        raise AssertionError("unknown NNI improvement policy should fail")
