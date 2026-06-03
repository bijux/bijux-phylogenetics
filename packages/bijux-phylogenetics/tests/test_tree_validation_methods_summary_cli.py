from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_report_tree_validation_methods_summary_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "tree-validation-methods-summary.md"

    exit_code = main(
        [
            "report",
            "tree-validation-methods-summary",
            str(tree_fixture("example_tree_unrooted.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["warning_count"] >= 1
    assert payload["metrics"]["blocked_context_count"] == 2
    assert payload["metrics"]["repair_item_count"] == 0
    assert payload["metrics"]["validity_decision"] == "valid_with_warnings"
    assert output_path.exists()
