from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from tests.support.fake_reference_parity import fake_ape_rscript


@pytest.mark.slow
def test_parity_cli_runs_live_ape_harness_and_writes_tables(
    tmp_path: Path, capsys
) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")
    summary_path = tmp_path / "ape-parity-summary.tsv"
    observation_path = tmp_path / "ape-parity-observations.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "ape-live",
            "--ape-rscript-executable",
            str(rscript),
            "--summary-out",
            str(summary_path),
            "--observations-out",
            str(observation_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["reference_source"] == "ape-live"
    assert payload["metrics"]["case_count"] == 180
    assert payload["metrics"]["function_count"] == 28
    assert payload["metrics"]["skipped_case_count"] == 0
    assert summary_path.exists()
    assert observation_path.exists()


def test_parity_cli_restricts_live_ape_cases(tmp_path: Path, capsys) -> None:
    rscript = fake_ape_rscript(tmp_path / "fake-ape-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "ape-live",
            "--ape-rscript-executable",
            str(rscript),
            "--ape-case",
            "get-mrca-balanced-two-tip",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    report = payload["data"]["report"]
    assert report["observations"][0]["case_id"] == "get-mrca-balanced-two-tip"
