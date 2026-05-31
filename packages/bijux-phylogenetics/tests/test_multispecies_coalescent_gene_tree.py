from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.simulation import simulate_multispecies_coalescent_gene_tree

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_simulate_multispecies_coalescent_gene_tree_records_species_tree_events() -> (
    None
):
    gene_tree, report = simulate_multispecies_coalescent_gene_tree(
        fixture("trees", "multispecies_coalescent_species_tree_3_taxa.nwk"),
        sample_count_table_path=fixture(
            "metadata", "multispecies_coalescent_sample_counts_3_taxa.tsv"
        ),
        population_size_table_path=fixture(
            "metadata", "multispecies_coalescent_population_sizes_3_taxa.tsv"
        ),
        seed=7,
    )

    assert report.model == "multispecies-coalescent"
    assert report.species_tip_count == 3
    assert report.gene_tip_count == 4
    assert report.deep_coalescence_total == 1
    assert sorted(gene_tree.tip_names) == ["A__1", "A__2", "B__1", "C__1"]
    assert [(row.species_taxon, row.sample_count) for row in report.sample_rows] == [
        ("A", 2),
        ("B", 1),
        ("C", 1),
    ]

    branch_rows = {row.species_branch: row for row in report.branch_rows}
    assert (
        branch_rows["A"].lineage_count_entering,
        branch_rows["A"].coalescent_event_count,
        branch_rows["A"].lineage_count_exiting,
        branch_rows["A"].extra_lineage_count,
    ) == (2, 1, 1, 0)
    assert (
        branch_rows["A|B"].branch_role,
        branch_rows["A|B"].lineage_count_entering,
        branch_rows["A|B"].coalescent_event_count,
        branch_rows["A|B"].lineage_count_exiting,
        branch_rows["A|B"].extra_lineage_count,
    ) == ("internal-branch", 2, 0, 2, 1)
    assert (
        branch_rows["A|B|C"].branch_role,
        branch_rows["A|B|C"].lineage_count_entering,
        branch_rows["A|B|C"].coalescent_event_count,
        branch_rows["A|B|C"].lineage_count_exiting,
        branch_rows["A|B|C"].included_in_deep_coalescence_total,
    ) == ("root-population", 3, 2, 1, False)
    assert [row.species_branch for row in report.event_rows] == [
        "A",
        "A|B|C",
        "A|B|C",
    ]


def test_simulate_multispecies_coalescent_gene_tree_rejects_non_ultrametric_tree(
    tmp_path: Path,
) -> None:
    species_tree_path = tmp_path / "non-ultrametric-species-tree.nwk"
    species_tree_path.write_text("((A:1,B:2):1,C:2);", encoding="utf-8")

    with pytest.raises(
        ValueError,
        match="requires an ultrametric species tree",
    ):
        simulate_multispecies_coalescent_gene_tree(species_tree_path)
