from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_topology_shape_reports_tree_metrics(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree-shape.tsv"
    exit_code = main(
        [
            "topology",
            "shape",
            str(tree_fixture("example_tree_ladderized.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 1
    assert payload["metrics"]["cherry_count"] == 1
    assert payload["metrics"]["sackin_imbalance_index"] == 9
    assert payload["metrics"]["colless_imbalance_index"] == 3.0
    assert payload["metrics"]["tree_height_edges"] == 3
    assert payload["metrics"]["imbalance_summary"] == "ladderized"
    assert payload["outputs"] == [str(output)]
    assert payload["data"]["rows"][0]["comb_like"] is True
    assert output.read_text(encoding="utf-8").startswith(
        "source_path\ttree_index\trooted\t"
    )


def test_cli_tree_set_shape_reports_aggregate_counts(tmp_path: Path, capsys) -> None:
    output = tmp_path / "tree-set-shape.tsv"
    exit_code = main(
        [
            "tree-set",
            "shape",
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
    assert payload["metrics"]["balanced_tree_count"] == 3
    assert payload["metrics"]["ladderized_tree_count"] == 0
    assert payload["metrics"]["star_like_tree_count"] == 0
    assert payload["metrics"]["comb_like_tree_count"] == 0
    assert payload["outputs"] == [str(output)]
    assert payload["data"]["aggregate"]["mean_cherry_count"] == 2.0
    assert (
        output.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("source_path\ttree_index\trooted\t")
    )
