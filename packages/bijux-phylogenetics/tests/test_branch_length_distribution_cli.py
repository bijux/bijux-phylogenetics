from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_topology_branch_lengths_reports_summary_and_output(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "branch-lengths.tsv"
    exit_code = main(
        [
            "topology",
            "branch-lengths",
            str(tree_fixture("example_tree_long_branch.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 1
    assert payload["metrics"]["branch_count"] == 6
    assert payload["metrics"]["zero_length_branch_count"] == 0
    assert payload["metrics"]["negative_branch_count"] == 0
    assert payload["metrics"]["long_outlier_count"] == 1
    assert payload["metrics"]["median_branch_length"] == 0.1
    assert payload["outputs"] == [str(output)]
    assert any(row["long_outlier"] for row in payload["data"]["rows"])
    assert output.read_text(encoding="utf-8").startswith(
        "source_path\ttree_index\tnode\tbranch_type\t"
    )


def test_cli_tree_set_branch_lengths_reports_aggregate_counts(
    tmp_path: Path, capsys
) -> None:
    output = tmp_path / "tree-set-branch-lengths.tsv"
    exit_code = main(
        [
            "tree-set",
            "branch-lengths",
            str(tree_fixture("example_tree_set_branch_lengths.nwk")),
            "--out",
            str(output),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["branch_count"] == 18
    assert payload["metrics"]["zero_length_branch_count"] == 3
    assert payload["metrics"]["negative_branch_count"] == 0
    assert payload["metrics"]["long_outlier_count"] == 1
    assert payload["outputs"] == [str(output)]
    assert payload["data"]["aggregate"]["maximum_branch_length"] == 1.0
    assert (
        output.read_text(encoding="utf-8")
        .splitlines()[0]
        .startswith("source_path\ttree_index\tnode\tbranch_type\t")
    )
