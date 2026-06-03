from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
import bijux_phylogenetics.trees as trees_api
from bijux_phylogenetics.trees import (
    GeneTreeQuartetConcordanceReport,
    GeneTreeQuartetConcordanceRow,
    compute_gene_tree_quartet_concordance_factors,
    write_gene_tree_quartet_concordance_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_package_tree_gateway_exports_quartet_concordance_surface() -> None:
    assert trees_api.GeneTreeQuartetConcordanceRow is GeneTreeQuartetConcordanceRow
    assert (
        trees_api.GeneTreeQuartetConcordanceReport is GeneTreeQuartetConcordanceReport
    )
    assert (
        trees_api.compute_gene_tree_quartet_concordance_factors
        is compute_gene_tree_quartet_concordance_factors
    )
    assert (
        trees_api.write_gene_tree_quartet_concordance_table
        is write_gene_tree_quartet_concordance_table
    )


def test_compute_gene_tree_quartet_concordance_matches_hand_counted_fixture() -> None:
    report = compute_gene_tree_quartet_concordance_factors(
        fixture("quartet_concordance_species_tree_4_taxa.nwk"),
        fixture("quartet_concordance_gene_trees_4_taxa.nwk"),
    )

    assert report.tree_count == 5
    assert report.shared_taxa == ["A", "B", "C", "D"]
    assert report.branch_count == 1
    assert report.total_quartet_count == 5
    assert report.concordant_quartet_count == 2
    assert report.discordant_first_quartet_count == 1
    assert report.discordant_second_quartet_count == 1
    assert report.uninformative_quartet_count == 1
    assert report.informative_quartet_count == 4
    assert report.rows == [
        GeneTreeQuartetConcordanceRow(
            branch_id="A|B::C|D",
            left_taxa=["A", "B"],
            right_taxa=["C", "D"],
            quartet_count_per_tree=1,
            concordant_quartet_count=2,
            discordant_first_quartet_count=1,
            discordant_second_quartet_count=1,
            uninformative_quartet_count=1,
            informative_quartet_count=4,
            concordance_factor=0.5,
            concordant_frequency=0.4,
            discordant_first_frequency=0.2,
            discordant_second_frequency=0.2,
            uninformative_frequency=0.2,
        )
    ]


def test_compute_gene_tree_quartet_concordance_handles_multiple_quartets_per_branch() -> (
    None
):
    report = compute_gene_tree_quartet_concordance_factors(
        fixture("quartet_concordance_species_tree_5_taxa.nwk"),
        fixture("quartet_concordance_gene_trees_5_taxa.nwk"),
    )

    assert report.branch_count == 2
    assert report.total_quartet_count == 12
    assert report.concordant_quartet_count == 12
    assert report.discordant_first_quartet_count == 0
    assert report.discordant_second_quartet_count == 0
    assert report.uninformative_quartet_count == 0
    assert report.informative_quartet_count == 12
    assert [
        (
            row.branch_id,
            row.left_taxa,
            row.right_taxa,
            row.quartet_count_per_tree,
            row.concordant_quartet_count,
            row.concordance_factor,
        )
        for row in report.rows
    ] == [
        ("A|B::C|D|E", ["A", "B"], ["C", "D", "E"], 3, 6, 1.0),
        ("D|E::A|B|C", ["D", "E"], ["A", "B", "C"], 3, 6, 1.0),
    ]


def test_compute_gene_tree_quartet_concordance_requires_exact_taxon_match(
    tmp_path: Path,
) -> None:
    species_tree = tmp_path / "quartet-concordance-species-tree.nwk"
    gene_tree_set = tmp_path / "quartet-concordance-gene-tree-set.nwk"
    species_tree.write_text("((A,B),(C,D));\n", encoding="utf-8")
    gene_tree_set.write_text("((A,B),(C,E));\n", encoding="utf-8")

    with pytest.raises(InvalidAlignmentError, match="exact same taxon set"):
        compute_gene_tree_quartet_concordance_factors(
            species_tree,
            gene_tree_set,
        )
