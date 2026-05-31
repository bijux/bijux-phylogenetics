from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.phylo.likelihood import (
    build_nucleotide_likelihood_starting_tree_pool_from_alignment,
    search_nucleotide_likelihood_multi_start,
    validate_nucleotide_likelihood_starting_tree_from_alignment,
)
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidBranchLengthError,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_starting_tree_validation_accepts_binary_root_representation() -> None:
    validate_nucleotide_likelihood_starting_tree_from_alignment(
        fixture("trees", "jc69_likelihood_nni_start_tree_4_taxa.nwk"),
        fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
        model_name="jc69",
    )


def test_likelihood_starting_tree_validation_rejects_duplicate_tip_labels() -> None:
    with pytest.raises(
        AlignmentTaxonMismatchError,
        match="uniquely named tree tips",
    ):
        validate_nucleotide_likelihood_starting_tree_from_alignment(
            fixture("trees", "example_tree_duplicate.nwk"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
        )


def test_likelihood_starting_tree_validation_rejects_missing_alignment_taxa() -> None:
    with pytest.raises(
        AlignmentTaxonMismatchError,
        match="alignment-only taxa: D",
    ):
        validate_nucleotide_likelihood_starting_tree_from_alignment(
            loads_newick("((A:0.1,B:0.1):0.2,C:0.3);"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
        )


def test_likelihood_starting_tree_validation_rejects_negative_branch_lengths() -> None:
    with pytest.raises(
        InvalidBranchLengthError,
        match="does not accept negative branch lengths",
    ):
        validate_nucleotide_likelihood_starting_tree_from_alignment(
            loads_newick("(((A:-0.1,B:0.2):0.3,C:0.4):0.5,D:0.6);"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
        )


def test_likelihood_starting_tree_validation_rejects_missing_branch_lengths() -> None:
    with pytest.raises(
        InvalidBranchLengthError,
        match="requires explicit branch lengths on every edge",
    ):
        validate_nucleotide_likelihood_starting_tree_from_alignment(
            loads_newick("(((A,B):0.3,C:0.4):0.5,D:0.6);"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
        )


def test_likelihood_starting_tree_validation_rejects_nonbinary_topology() -> None:
    with pytest.raises(
        ValueError,
        match="requires a strictly binary tree",
    ):
        validate_nucleotide_likelihood_starting_tree_from_alignment(
            fixture("trees", "example_tree_polytomy.nwk"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
        )


def test_likelihood_starting_tree_pool_rejects_unrooted_representation() -> None:
    with pytest.raises(
        ValueError,
        match="requires a rooted binary tree",
    ):
        build_nucleotide_likelihood_starting_tree_pool_from_alignment(
            fixture("trees", "example_tree_unrooted.nwk"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
            random_start_tree_count=1,
            random_start_tree_seed=17,
        )


def test_likelihood_multi_start_search_rejects_invalid_start_tree_before_search() -> None:
    with pytest.raises(
        AlignmentTaxonMismatchError,
        match="alignment-only taxa: D",
    ):
        search_nucleotide_likelihood_multi_start(
            loads_newick("((A:0.1,B:0.1):0.2,C:0.3);"),
            fixture("alignments", "jc69_likelihood_nni_alignment_4_taxa.fasta"),
            model_name="jc69",
            start_tree_count=2,
            start_tree_seed=17,
        )
