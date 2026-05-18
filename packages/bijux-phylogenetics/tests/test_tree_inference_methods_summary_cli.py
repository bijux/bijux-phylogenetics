from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "expected" / "fasta_to_tree"


def workflow_fixture(dataset: str) -> Path:
    return FIXTURES / dataset / f"{dataset}.manifest.json"


def test_cli_report_tree_inference_methods_summary_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "tree-inference-methods-summary.md"

    exit_code = main(
        [
            "report",
            "tree-inference-methods-summary",
            str(workflow_fixture("pleistocene-bear-cytb-fragments")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["warning_count"] >= 1
    assert payload["metrics"]["selected_model"] == "HKY+F"
    assert payload["metrics"]["bootstrap_replicates"] == 1000
    assert payload["metrics"]["trimmed_alignment_length"] == 1140
    assert payload["metrics"]["supported_node_count"] == 2
    assert output_path.exists()
