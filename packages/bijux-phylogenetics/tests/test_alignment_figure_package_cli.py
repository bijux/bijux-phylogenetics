from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "alignments" / name


def test_cli_report_alignment_package_writes_review_directory(
    tmp_path: Path, capsys
) -> None:
    output_dir = tmp_path / "alignment-figure-package"
    exit_code = main(
        [
            "report",
            "alignment-package",
            str(fixture("example_alignment.fasta")),
            "--out-dir",
            str(output_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["publication_ready"] is True
    assert payload["metrics"]["quality_score"] == 100.0
    assert payload["metrics"]["suspicious_alignment"] is False
    assert payload["metrics"]["heatmap_row_count"] == 4
    assert payload["metrics"]["heatmap_bin_count"] == 8
    assert payload["metrics"]["plotted_window_count"] == 1
    assert payload["metrics"]["plotted_sequence_count"] == 4
    assert payload["metrics"]["reviewer_audit_item_count"] == 5
    assert len(payload["outputs"]) == 12
    assert (output_dir / "alignment-missingness-heatmap.svg").exists()
    assert (output_dir / "alignment-site-quality-summary.svg").exists()
    assert (output_dir / "alignment-sequence-quality-panel.svg").exists()
    assert (output_dir / "alignment-quality-review.html").exists()
    assert (output_dir / "alignment-quality-package.manifest.json").exists()
    assert (output_dir / "figure-reproducibility.manifest.json").exists()
    assert (output_dir / "reviewer-audit-checklist.tsv").exists()
