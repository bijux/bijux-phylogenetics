from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.topology import (
    write_duplication_loss_transfer_event_table,
    write_duplication_loss_transfer_taxon_map_table,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_write_dlt_event_table_emits_reconciliation_ledger(tmp_path: Path) -> None:
    output_path = tmp_path / "duplication-loss-transfer.tsv"

    write_duplication_loss_transfer_event_table(
        output_path,
        fixture("trees", "duplication_loss_transfer_species_tree_4_taxa.nwk"),
        fixture("trees", "duplication_loss_transfer_gene_tree_4_tips.nwk"),
        taxon_map_path=fixture(
            "metadata",
            "duplication_loss_transfer_gene_taxon_map_4_tips.tsv",
        ),
    )

    assert output_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "gene_node\tgene_node_name\tdescendant_gene_tips\tevent_type\tmapped_species_branch\tmapped_descendant_species\tleft_child_gene_node\tright_child_gene_node\tleft_child_species_branch\tright_child_species_branch\ttransferred_child_side\ttransfer_recipient_branch\tloss_branches\tevent_cost\treconciliation_score\tduplication_cost\tloss_cost\ttransfer_cost",
        "origin\t\tA__1|A__2|C__1|D__1\torigin\tA|B|C|D\tA|B|C|D\troot:clade:A__1|A__2|C__1|D__1\t\tC|D\t\t\tC|D\tA|B\t1.0\t6.0\t2.0\t1.0\t3.0",
        "root:clade:A__1|A__2|C__1|D__1\t\tA__1|A__2|C__1|D__1\ttransfer\tC|D\tC|D\troot:clade:A__1|A__2|C__1|D__1/clade:C__1|D__1\troot:clade:A__1|A__2|C__1|D__1/clade:A__1|A__2\tC|D\tA\tright\tA\t\t3.0\t6.0\t2.0\t1.0\t3.0",
        "root:clade:A__1|A__2|C__1|D__1/clade:C__1|D__1\t\tC__1|D__1\tspeciation\tC|D\tC|D\troot:clade:A__1|A__2|C__1|D__1/clade:C__1|D__1/taxon:C__1\troot:clade:A__1|A__2|C__1|D__1/clade:C__1|D__1/taxon:D__1\tC\tD\t\t\t\t0.0\t6.0\t2.0\t1.0\t3.0",
    ]


def test_write_dlt_taxon_map_table_emits_association_rows(tmp_path: Path) -> None:
    output_path = tmp_path / "duplication-loss-transfer-taxon-map.tsv"

    write_duplication_loss_transfer_taxon_map_table(
        output_path,
        fixture("trees", "duplication_loss_transfer_species_tree_4_taxa.nwk"),
        fixture("trees", "duplication_loss_transfer_gene_tree_4_tips.nwk"),
        taxon_map_path=fixture(
            "metadata",
            "duplication_loss_transfer_gene_taxon_map_4_tips.tsv",
        ),
    )

    assert output_path.read_text(encoding="utf-8").splitlines() == [
        "gene_taxon\tspecies_taxon",
        "A__1\tA",
        "A__2\tA",
        "C__1\tC",
        "D__1\tD",
    ]
