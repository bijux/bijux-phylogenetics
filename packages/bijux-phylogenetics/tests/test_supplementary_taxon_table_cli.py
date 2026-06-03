from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_cli_report_supplementary_taxon_table_writes_table(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-taxa.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-taxon-table",
            "--tree",
            str(fixture("example_taxon_workflow_tree.nwk")),
            "--metadata",
            str(fixture("example_taxon_workflow_metadata.csv")),
            "--traits",
            str(fixture("example_taxon_workflow_traits.csv")),
            "--alignment",
            str(fixture("example_taxon_workflow_alignment.fasta")),
            "--filtered-alignment",
            str(fixture("example_taxon_workflow_filtered_alignment.fasta")),
            "--inference-tree",
            str(fixture("example_taxon_workflow_inference.nwk")),
            "--reported-taxa",
            str(fixture("example_taxon_workflow_reported.csv")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert output_path.exists()
    assert payload["metrics"]["row_count"] == 4
    assert payload["metrics"]["analysis_included_count"] == 2
    assert payload["metrics"]["analysis_excluded_count"] == 2
    assert payload["metrics"]["reporting_retained_count"] == 1
    assert payload["metrics"]["reporting_dropped_count"] == 3
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["taxon"] == "A"
    assert "metadata_group" in rows[0]
