from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def metadata_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_cli_topology_clades_reports_clade_rows_and_metadata_columns(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "clades.tsv"
    exit_code = main(
        [
            "topology",
            "clades",
            str(tree_fixture("example_tree_support_conflict_left.nwk")),
            "--metadata",
            str(metadata_fixture("example_metadata.tsv")),
            "--metadata-column",
            "species",
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 1
    assert payload["metrics"]["clade_count"] == 7
    assert payload["metrics"]["metadata_column_count"] == 1
    assert payload["outputs"] == [str(output)]
    assert payload["data"]["metadata_columns"] == ["species"]
    assert any(row["clade_id"] == "A|B" for row in payload["data"]["rows"])
    assert output.read_text(encoding="utf-8").startswith(
        "source_path\ttree_index\tnode_kind\tclade_id\tnode_label\t"
    )


def test_cli_tree_set_clades_reports_per_tree_clade_rows(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "tree-set-clades.tsv"
    exit_code = main(
        [
            "tree-set",
            "clades",
            str(tree_fixture("example_tree_set_left.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["clade_count"] == 21
    assert payload["metrics"]["metadata_column_count"] == 0
    assert payload["outputs"] == [str(output)]
    assert (
        output.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("source_path\ttree_index\tnode_kind\tclade_id\t")
    )
