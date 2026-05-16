from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_compare_influence_reports_ranked_taxa_and_table_output(
    capsys, tmp_path: Path
) -> None:
    table_path = tmp_path / "taxon-influence.tsv"

    exit_code = main(
        [
            "compare",
            "influence",
            str(fixture("example_tree_taxon_influence_left.nwk")),
            str(fixture("example_tree_taxon_influence_right.nwk")),
            "--out",
            str(table_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["shared_taxa"] == 5
    assert payload["metrics"]["top_influential_taxon"] == "C"
    assert payload["metrics"]["taxa_with_topology_change"] >= 1
    assert payload["metrics"]["taxa_with_support_change"] >= 1
    assert payload["outputs"] == [str(table_path)]
    assert payload["data"]["rows"][0]["taxon"] == "C"
