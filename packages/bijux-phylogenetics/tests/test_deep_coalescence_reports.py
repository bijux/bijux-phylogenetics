from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    write_deep_coalescence_branch_table,
    write_deep_coalescence_taxon_map_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_deep_coalescence_branch_table_emits_species_branch_ledger(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "deep-coalescence.tsv"

    write_deep_coalescence_branch_table(
        output_path,
        fixture("trees", "deep_coalescence_species_tree_3_taxa.nwk"),
        fixture("trees", "deep_coalescence_gene_tree_4_tips.nwk"),
        taxon_map_path=fixture(
            "metadata",
            "deep_coalescence_gene_taxon_map_4_tips.tsv",
        ),
    )

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "species_branch\tbranch_role\tdescendant_species\tlineage_count_entering\tcoalescent_event_count\tlineage_count_exiting\textra_lineage_count\tincluded_in_deep_coalescence_total\tdeep_coalescence_total\tgene_tip_count",
        "A|B|C\troot-population\tA|B|C\t3\t2\t1\t0\tfalse\t2\t4",
        "C\ttip-branch\tC\t1\t0\t1\t0\ttrue\t2\t4",
        "A|B\tinternal-branch\tA|B\t3\t1\t2\t1\ttrue\t2\t4",
        "B\ttip-branch\tB\t1\t0\t1\t0\ttrue\t2\t4",
        "A\ttip-branch\tA\t2\t0\t2\t1\ttrue\t2\t4",
    ]


def test_write_deep_coalescence_taxon_map_table_emits_resolved_mapping_rows(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "deep-coalescence-taxon-map.tsv"

    write_deep_coalescence_taxon_map_table(
        output_path,
        fixture("trees", "deep_coalescence_species_tree_3_taxa.nwk"),
        fixture("trees", "deep_coalescence_gene_tree_4_tips.nwk"),
        taxon_map_path=fixture(
            "metadata",
            "deep_coalescence_gene_taxon_map_4_tips.tsv",
        ),
    )

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "gene_taxon\tspecies_taxon",
        "A__1\tA",
        "A__2\tA",
        "B__1\tB",
        "C__1\tC",
    ]
