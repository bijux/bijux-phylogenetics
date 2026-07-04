from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def tree_fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_report_tree_package_writes_review_directory(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "tree-report-package"
    exit_code = main(
        [
            "report",
            "tree-package",
            str(tree_fixture("example_tree_support_left.nwk")),
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["tip_count"] == 4
    assert payload["metrics"]["supported_branch_count"] == 3
    assert payload["metrics"]["rendered_support_count"] == 2
    assert payload["metrics"]["methods_warning_count"] == 0
    assert payload["metrics"]["reviewer_audit_item_count"] == 5
    assert payload["metrics"]["method_tier"] == "advisory"
    assert payload["metrics"]["method_inference_mode"] == "review-only"
    assert len(payload["outputs"]) == 8
    assert (output_dir / "tree-report.html").exists()
    assert (output_dir / "tree-image.svg").exists()
    assert (output_dir / "tree-validation-methods-summary.md").exists()
    assert (output_dir / "reviewer-audit-checklist.tsv").exists()
    assert (output_dir / "support-table.tsv").exists()
    assert (output_dir / "clade-table.tsv").exists()
    assert (output_dir / "branch-stats.tsv").exists()
