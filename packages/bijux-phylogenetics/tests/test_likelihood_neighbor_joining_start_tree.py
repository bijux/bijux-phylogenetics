from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodNeighborJoiningStartTreeReport,
    build_likelihood_neighbor_joining_tree,
    build_likelihood_neighbor_joining_tree_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_likelihood_gateway_exports_neighbor_joining_start_tree_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodNeighborJoiningStartTreeReport
        is NucleotideLikelihoodNeighborJoiningStartTreeReport
    )
    assert (
        likelihood_api.build_likelihood_neighbor_joining_tree
        is build_likelihood_neighbor_joining_tree
    )
    assert (
        likelihood_api.build_likelihood_neighbor_joining_tree_from_alignment
        is build_likelihood_neighbor_joining_tree_from_alignment
    )


def test_build_likelihood_neighbor_joining_tree_reports_alignment_derived_nj() -> None:
    tree, report = build_likelihood_neighbor_joining_tree_from_alignment(
        fixture("alignments", "example_alignment_distance.fasta")
    )

    assert tree.rooted is True
    assert (
        dumps_newick(tree)
        == "((A:0.0625,B:0.0625)Inner1:0,(C:0.0625,D:0.0625)Inner2:0.4375);"
    )
    assert report.algorithm == "nucleotide-likelihood-neighbor-joining-start-tree"
    assert report.distance_method == "neighbor-joining"
    assert report.distance_model == "p-distance"
    assert report.gap_handling == "pairwise-deletion"
    assert report.ambiguity_policy == "ignore"
    assert report.alignment_path == str(
        fixture("alignments", "example_alignment_distance.fasta")
    )
    assert report.taxon_count == 4
    assert report.site_count == 8
    assert report.pair_count == 6
    assert report.distance_warnings == []
    assert report.rooting_policy == "largest-root-clade"
    assert (
        report.unrooted_tree_newick
        == "((A:0.0625,B:0.0625)Inner1:0.4375,C:0.0625,D:0.0625)Inner2;"
    )
    assert report.rooted_tree_newick == dumps_newick(tree)


def test_build_likelihood_neighbor_joining_tree_roots_two_taxon_tree_without_rewriting_distances() -> (
    None
):
    tree, report = build_likelihood_neighbor_joining_tree_from_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )

    assert tree.rooted is True
    assert len(tree.root.children) == 2
    assert report.taxon_count == 2
    assert report.pair_count == 1
    assert report.rooted_tree_newick == dumps_newick(tree)
