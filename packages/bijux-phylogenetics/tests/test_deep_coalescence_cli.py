from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).resolve().parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_cli_compare_deep_coalescence_reports_metrics_and_writes_ledgers(
    tmp_path: Path,
    capsys,
) -> None:
    branch_table_path = tmp_path / "deep-coalescence.tsv"
    mapping_table_path = tmp_path / "deep-coalescence-taxon-map.tsv"

    exit_code = main(
        [
            "compare",
            "deep-coalescence",
            str(fixture("trees", "deep_coalescence_species_tree_3_taxa.nwk")),
            str(fixture("trees", "deep_coalescence_gene_tree_4_tips.nwk")),
            "--taxon-map",
            str(fixture("metadata", "deep_coalescence_gene_taxon_map_4_tips.tsv")),
            "--out",
            str(branch_table_path),
            "--mapping-out",
            str(mapping_table_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["observed_species_taxa"] == 3
    assert payload["metrics"]["species_only_taxa"] == 0
    assert payload["metrics"]["gene_tip_count"] == 4
    assert payload["metrics"]["deep_coalescence_total"] == 2
    assert payload["data"]["branch_rows"][2]["species_branch"] == "A|B"
    assert payload["data"]["branch_rows"][2]["extra_lineage_count"] == 1
    assert payload["outputs"] == [
        str(branch_table_path),
        str(mapping_table_path),
    ]
    assert branch_table_path.read_text(encoding="utf-8").startswith(
        "species_branch\tbranch_role\tdescendant_species\tlineage_count_entering\t"
    )
    assert mapping_table_path.read_text(encoding="utf-8") == (
        "gene_taxon\tspecies_taxon\nA__1\tA\nA__2\tA\nB__1\tB\nC__1\tC\n"
    )
