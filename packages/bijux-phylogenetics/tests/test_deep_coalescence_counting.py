from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.compare as compare_api
from bijux_phylogenetics.compare import (
    DeepCoalescenceBranchRow,
    DeepCoalescenceReport,
    DeepCoalescenceTaxonMapRow,
    compare_tree_paths,
    count_deep_coalescences,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_compare_gateway_exports_deep_coalescence_surface() -> None:
    assert compare_api.DeepCoalescenceBranchRow is DeepCoalescenceBranchRow
    assert compare_api.DeepCoalescenceReport is DeepCoalescenceReport
    assert compare_api.DeepCoalescenceTaxonMapRow is DeepCoalescenceTaxonMapRow
    assert compare_api.count_deep_coalescences is count_deep_coalescences


def test_count_deep_coalescences_matches_hand_checked_extra_lineage_example() -> None:
    report = count_deep_coalescences(
        fixture("trees", "deep_coalescence_species_tree_3_taxa.nwk"),
        fixture("trees", "deep_coalescence_gene_tree_4_tips.nwk"),
        taxon_map_path=fixture(
            "metadata",
            "deep_coalescence_gene_taxon_map_4_tips.tsv",
        ),
    )

    assert report.observed_species_taxa == ["A", "B", "C"]
    assert report.species_only_taxa == []
    assert report.gene_tip_count == 4
    assert report.deep_coalescence_total == 2
    assert [(row.gene_taxon, row.species_taxon) for row in report.mapping_rows] == [
        ("A__1", "A"),
        ("A__2", "A"),
        ("B__1", "B"),
        ("C__1", "C"),
    ]
    branch_rows = {row.species_branch: row for row in report.branch_rows}
    assert (
        branch_rows["A"].branch_role,
        branch_rows["A"].lineage_count_entering,
        branch_rows["A"].coalescent_event_count,
        branch_rows["A"].lineage_count_exiting,
        branch_rows["A"].extra_lineage_count,
    ) == ("tip-branch", 2, 0, 2, 1)
    assert (
        branch_rows["A|B"].branch_role,
        branch_rows["A|B"].lineage_count_entering,
        branch_rows["A|B"].coalescent_event_count,
        branch_rows["A|B"].lineage_count_exiting,
        branch_rows["A|B"].extra_lineage_count,
    ) == ("internal-branch", 3, 1, 2, 1)
    assert (
        branch_rows["A|B|C"].branch_role,
        branch_rows["A|B|C"].lineage_count_entering,
        branch_rows["A|B|C"].coalescent_event_count,
        branch_rows["A|B|C"].lineage_count_exiting,
        branch_rows["A|B|C"].included_in_deep_coalescence_total,
    ) == ("root-population", 3, 2, 1, False)


def test_count_deep_coalescences_is_not_robinson_foulds_distance() -> None:
    species_tree_path = fixture("trees", "deep_coalescence_species_tree_3_taxa.nwk")
    gene_tree_path = fixture(
        "trees",
        "deep_coalescence_discordant_gene_tree_3_taxa.nwk",
    )

    report = count_deep_coalescences(species_tree_path, gene_tree_path)
    comparison = compare_tree_paths(
        species_tree_path,
        gene_tree_path,
        rf_mode="rooted",
        taxon_overlap_policy="require-identical",
    )

    assert report.deep_coalescence_total == 1
    assert comparison.robinson_foulds_distance == 2


def test_count_deep_coalescences_requires_taxon_map_for_nonmatching_gene_tips() -> None:
    with pytest.raises(
        ValueError,
        match="requires --taxon-map when gene tips do not exactly match species-tree taxa",
    ):
        count_deep_coalescences(
            fixture("trees", "deep_coalescence_species_tree_3_taxa.nwk"),
            fixture("trees", "deep_coalescence_gene_tree_4_tips.nwk"),
        )
