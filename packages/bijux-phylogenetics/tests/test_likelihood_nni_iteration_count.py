from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import search_nucleotide_likelihood_nni_from_alignment

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_nni_search_reports_full_iteration_count() -> None:
    report = search_nucleotide_likelihood_nni_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        upper_branch_length_bound=1.0,
    )

    assert report.accepted_move_count == 2
    assert report.iteration_count == 3
    assert report.iteration_count == max(
        row.iteration for row in report.candidate_rows
    )
    assert report.iteration_count == report.accepted_move_count + 1
    assert report.trace_rows[-1].stopping_reason == "no-improving-neighbor"
