from __future__ import annotations

from pathlib import Path

import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    GeneTreeConflictArtifactReport,
    GeneTreeConflictQuartetSummary,
    GeneTreeConflictReferenceTree,
    GeneTreeConflictSummaryReport,
    summarize_gene_tree_conflicts,
    write_gene_tree_conflict_artifacts,
    write_gene_tree_conflict_quartet_table,
    write_gene_tree_conflict_summary_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_tree_gateway_exports_gene_tree_conflict_surface() -> None:
    assert trees_api.GeneTreeConflictArtifactReport is GeneTreeConflictArtifactReport
    assert trees_api.GeneTreeConflictQuartetSummary is GeneTreeConflictQuartetSummary
    assert trees_api.GeneTreeConflictReferenceTree is GeneTreeConflictReferenceTree
    assert trees_api.GeneTreeConflictSummaryReport is GeneTreeConflictSummaryReport
    assert trees_api.summarize_gene_tree_conflicts is summarize_gene_tree_conflicts
    assert (
        trees_api.write_gene_tree_conflict_artifacts
        is write_gene_tree_conflict_artifacts
    )
    assert (
        trees_api.write_gene_tree_conflict_quartet_table
        is write_gene_tree_conflict_quartet_table
    )
    assert (
        trees_api.write_gene_tree_conflict_summary_table
        is write_gene_tree_conflict_summary_table
    )


def test_summarize_gene_tree_conflicts_reports_exact_known_conflicting_clades() -> None:
    report = summarize_gene_tree_conflicts(
        fixture("example_tree_set_left.nwk"),
        credibility_threshold=0.3,
    )

    assert report.tree_count == 3
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.reference_tree.selection_method == (
        "dominant-rooted-topology-representative"
    )
    assert report.reference_tree.frequency == 0.666666666666667
    assert report.reference_tree.newick == "((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);"
    assert [
        (row.clade, row.frequency) for row in report.clade_frequencies.clade_frequencies
    ] == [
        ("A|B", 0.666666666666667),
        ("A|C", 0.333333333333333),
        ("B|D", 0.333333333333333),
        ("C|D", 0.666666666666667),
    ]
    assert report.quartet_concordance.branch_count == 1
    assert report.quartet_concordance.total_quartet_count == 3
    assert report.quartet_concordance.concordant_quartet_count == 2
    assert report.quartet_concordance.discordant_first_quartet_count == 1
    assert report.quartet_concordance.discordant_second_quartet_count == 0
    assert report.quartet_concordance.uninformative_quartet_count == 0
    assert report.rogue_taxa.rows[0].taxon == "A"
    assert report.clade_conflicts.conflict_count == 4
    assert [
        (row.left_clade, row.right_clade, row.combined_frequency)
        for row in report.clade_conflicts.conflicts
    ] == [
        ("A|B", "A|C", 1.0),
        ("A|B", "B|D", 1.0),
        ("A|C", "C|D", 1.0),
        ("B|D", "C|D", 1.0),
    ]
