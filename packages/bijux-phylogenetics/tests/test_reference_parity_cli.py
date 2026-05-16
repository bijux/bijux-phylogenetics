from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def test_reference_parity_cli_reports_pass_and_writes_tables(
    tmp_path: Path,
    capsys,
) -> None:
    summary_path = tmp_path / "reference-parity-summary.tsv"
    observation_path = tmp_path / "reference-parity-observations.tsv"
    exit_code = main(
        [
            "parity",
            "--summary-out",
            str(summary_path),
            "--observations-out",
            str(observation_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)
    assert exit_code == 0
    assert payload["metrics"]["all_passed"] is True
    assert payload["metrics"]["case_count"] == 14
    assert payload["metrics"]["method_count"] == 10
    assert payload["metrics"]["extended"] is False
    assert Path(payload["data"]["summary_table"]).exists()
    assert Path(payload["data"]["observation_table"]).exists()
