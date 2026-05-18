from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from tests.support.fake_geiger_parity import fake_geiger_rscript

pytestmark = pytest.mark.slow


def test_parity_cli_runs_live_geiger_harness_and_writes_tables(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    summary_path = tmp_path / "geiger-parity-summary.tsv"
    observation_path = tmp_path / "geiger-parity-observations.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
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
    assert payload["metrics"]["reference_source"] == "geiger-live"
    assert payload["metrics"]["case_count"] == 4
    assert payload["metrics"]["function_count"] == 3
    assert payload["metrics"]["skipped_case_count"] == 0
    assert summary_path.exists()
    assert observation_path.exists()


def test_parity_cli_restricts_live_geiger_cases(tmp_path: Path, capsys) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitcontinuous-eb-early-burst-rate-recovery",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitcontinuous-eb-early-burst-rate-recovery"
    assert observation["model_name"] == "EB"
