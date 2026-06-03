from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood.f81 as f81_likelihood
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_f81_likelihood_rejects_non_normalized_base_frequencies() -> None:
    with pytest.raises(
        InvalidAlignmentError,
        match="F81 likelihood base frequencies must sum to one within tolerance",
    ):
        f81_likelihood.evaluate_f81_tree_likelihood_from_alignment(
            fixture("trees", "f81_likelihood_tree_2_taxa.nwk"),
            fixture("alignments", "f81_likelihood_alignment_2_taxa.fasta"),
            base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.4},
        )
