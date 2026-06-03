from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_tree_set_compatibility_graph_writes_expected_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "clade-compatibility"

    exit_code = main(
        [
            "tree-set",
            "compatibility-graph",
            str(fixture("clade_compatibility_tree_set.nwk")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["tree_count"] == 2
    assert payload["metrics"]["shared_taxon_count"] == 4
    assert payload["metrics"]["node_count"] == 4
    assert payload["metrics"]["edge_count"] == 6
    assert payload["metrics"]["compatible_edge_count"] == 2
    assert payload["metrics"]["conflict_edge_count"] == 4
    assert sorted(Path(path).name for path in payload["outputs"]) == [
        "clade-compatibility-edges.tsv",
        "clade-compatibility-nodes.tsv",
        "clade-compatibility.dot",
    ]

    node_lines = (
        (out_dir / "clade-compatibility-nodes.tsv")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    edge_lines = (
        (out_dir / "clade-compatibility-edges.tsv")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    dot_text = (out_dir / "clade-compatibility.dot").read_text(encoding="utf-8")

    assert node_lines[0] == (
        "clade\ttree_count\tfrequency\tcompatible_neighbor_count\tconflict_neighbor_count"
    )
    assert edge_lines[0] == (
        "left_clade\tright_clade\tcompatibility_relation\tcompatibility_reason\t"
        "left_tree_count\tright_tree_count\tleft_frequency\tright_frequency"
    )
    assert dot_text.startswith("graph clade_compatibility {\n")
