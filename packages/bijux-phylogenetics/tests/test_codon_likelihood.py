from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood.codon import (
    evaluate_codon_ctmc_tree_likelihood_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_codon_ctmc_likelihood_returns_finite_report() -> None:
    report = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 2
    assert report.pattern_count == 2
    assert report.compression_used is False
    assert report.state_count == 61
    assert report.genetic_code_id == 1
    assert report.genetic_code_name == "Standard"
    assert report.codon_frequency_source == "uniform"
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert math.isfinite(report.log_likelihood)
