from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_nni_best_improvement_records_full_candidate_ledger() -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    iteration_rows = {
        iteration: [row for row in report.candidate_rows if row.iteration == iteration]
        for iteration in {row.iteration for row in report.candidate_rows}
    }

    assert report.accepted_move_count == 2
    assert len(report.candidate_rows) == 12
    assert sorted(iteration_rows) == [1, 2, 3]
    assert [len(iteration_rows[index]) for index in [1, 2, 3]] == [4, 4, 4]
    assert [
        [row.candidate_order for row in iteration_rows[index]] for index in [1, 2, 3]
    ] == [[1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]]

    selected_rows = [row for row in report.candidate_rows if row.selected_best_move]
    assert len(selected_rows) == report.accepted_move_count
    assert [row.iteration for row in selected_rows] == [1, 2]
    assert selected_rows[0].pivot_branch_id == "A|C"
    assert selected_rows[0].sibling_clade_id == "B"
    assert selected_rows[0].exchanged_clade_id == "C"
    assert math.isclose(
        selected_rows[0].log_likelihood_delta,
        17.499717347003227,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert (
        selected_rows[0].candidate_tree_newick == report.trace_rows[1].tree_after_newick
    )
    assert selected_rows[1].pivot_branch_id == "A|B|C"
    assert selected_rows[1].sibling_clade_id == "D"
    assert selected_rows[1].exchanged_clade_id == "A|B"
    assert math.isclose(
        selected_rows[1].log_likelihood_delta,
        2.8075503046659165,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert (
        selected_rows[1].candidate_tree_newick == report.trace_rows[2].tree_after_newick
    )


def test_likelihood_nni_best_improvement_marks_local_optimum_iteration_without_choice() -> (
    None
):
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    local_optimum_rows = [row for row in report.candidate_rows if row.iteration == 3]

    assert local_optimum_rows
    assert all(row.selected_best_move is False for row in local_optimum_rows)
    assert all(row.improving_move is False for row in local_optimum_rows)
    assert all(row.log_likelihood_delta <= 0.0 for row in local_optimum_rows)
    assert report.trace_rows[-1].stopping_reason == "no-improving-neighbor"
