from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_multi_start_from_alignment,
)

pytestmark = pytest.mark.slow

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_multi_start_reports_final_likelihood_ranks() -> None:
    report = search_nucleotide_likelihood_multi_start_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
        start_tree_count=3,
        start_tree_seed=17,
        upper_branch_length_bound=1.0,
    )

    assert report.algorithm == "nucleotide-likelihood-multi-start-search"
    assert report.best_run_source_label == "input-tree"
    assert math.isclose(
        report.best_final_log_likelihood,
        -34.13524969797671,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [row.start_tree_source_label for row in report.run_summaries] == [
        "input-tree",
        "random-tree-seed-17",
        "random-tree-seed-18",
    ]
    assert [row.final_likelihood_rank for row in report.run_summaries] == [1, 2, 3]
    assert [row.best_run for row in report.run_summaries] == [True, False, False]
