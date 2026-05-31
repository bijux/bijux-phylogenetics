from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_multi_start_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_multi_start_reports_per_run_iteration_counts() -> None:
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    assert [row.search_iteration_count for row in report.run_summaries] == [3, 3, 3]
    assert all(
        row.search_iteration_count >= row.accepted_move_count
        for row in report.run_summaries
    )
