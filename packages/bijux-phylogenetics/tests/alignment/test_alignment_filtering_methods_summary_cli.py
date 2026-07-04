from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    direct = Path(__file__).parent / "fixtures" / "alignments" / name
    if direct.exists():
        return direct
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_cli_report_alignment_filtering_methods_summary_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "alignment-filtering-methods-summary.md"

    exit_code = main(
        [
            "report",
            "alignment-filtering-methods-summary",
            str(fixture("example_alignment_filtering.fasta")),
            "--profile",
            "moderate",
            "--group-table",
            str(fixture("example_alignment_groups.tsv")),
            "--group-column",
            "region",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["warning_count"] >= 1
    assert payload["metrics"]["removed_site_count"] == 4
    assert payload["metrics"]["removed_sequence_count"] == 1
    assert payload["metrics"]["retained_sequence_count"] == 3
    assert payload["metrics"]["retained_alignment_length"] == 8
    assert payload["metrics"]["profile_name"] == "moderate"
    assert output_path.exists()
