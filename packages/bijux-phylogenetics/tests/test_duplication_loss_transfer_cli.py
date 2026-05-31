from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_cli_compare_duplication_loss_transfer_reports_metrics_and_writes_ledgers(
    tmp_path: Path,
    capsys,
) -> None:
    event_table_path = tmp_path / "duplication-loss-transfer.tsv"
    mapping_table_path = tmp_path / "duplication-loss-transfer-taxon-map.tsv"

    exit_code = main(
        [
            "compare",
            "duplication-loss-transfer",
            str(fixture("trees", "duplication_loss_transfer_species_tree_4_taxa.nwk")),
            str(fixture("trees", "duplication_loss_transfer_gene_tree_4_tips.nwk")),
            "--taxon-map",
            str(
                fixture(
                    "metadata",
                    "duplication_loss_transfer_gene_taxon_map_4_tips.tsv",
                )
            ),
            "--out",
            str(event_table_path),
            "--mapping-out",
            str(mapping_table_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["observed_species_taxa"] == 3
    assert payload["metrics"]["species_only_taxa"] == 1
    assert payload["metrics"]["gene_tip_count"] == 4
    assert payload["metrics"]["reconciliation_score"] == 6.0
    assert payload["metrics"]["duplication_event_count"] == 1
    assert payload["metrics"]["loss_event_count"] == 1
    assert payload["metrics"]["transfer_event_count"] == 1
    assert payload["metrics"]["speciation_event_count"] == 1
    assert payload["data"]["event_rows"][1]["event_type"] == "transfer"
    assert payload["data"]["event_rows"][1]["transfer_recipient_branch"] == "A"
    assert payload["outputs"] == [
        str(event_table_path),
        str(mapping_table_path),
    ]
    assert event_table_path.read_text(encoding="utf-8").startswith(
        "gene_node\tgene_node_name\tdescendant_gene_tips\tevent_type\tmapped_species_branch\t"
    )
    assert mapping_table_path.read_text(encoding="utf-8") == (
        "gene_taxon\tspecies_taxon\nA__1\tA\nA__2\tA\nC__1\tC\nD__1\tD\n"
    )
