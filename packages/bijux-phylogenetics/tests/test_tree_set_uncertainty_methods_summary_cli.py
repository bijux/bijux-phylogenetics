from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "trees"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_cli_tree_set_methods_summary_writes_metrics(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "tree-set-uncertainty-methods-summary.md"

    exit_code = main(
        [
            "tree-set",
            "methods-summary",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["warning_count"] >= 1
    assert payload["metrics"]["tree_count"] == 3
    assert payload["metrics"]["rooted_topology_count"] == 2
    assert payload["metrics"]["topology_cluster_count"] == 2
    assert payload["metrics"]["unstable_taxon_count"] == 4
    assert payload["metrics"]["multimodal"] is True
    assert output_path.exists()
