from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from tests.support.fake_phytools_parity import fake_phytools_rscript

pytestmark = pytest.mark.slow


def test_parity_cli_runs_live_phytools_harness_and_writes_tables(
    tmp_path: Path, capsys
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")
    summary_path = tmp_path / "phytools-parity-summary.tsv"
    observation_path = tmp_path / "phytools-parity-observations.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "phytools-live",
            "--phytools-rscript-executable",
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
    assert payload["metrics"]["reference_source"] == "phytools-live"
    assert payload["metrics"]["case_count"] == 65
    assert payload["metrics"]["function_count"] == 21
    assert payload["metrics"]["skipped_case_count"] == 0
    assert summary_path.exists()
    assert observation_path.exists()


def test_parity_cli_restricts_live_phytools_cases(tmp_path: Path, capsys) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "phytools-live",
            "--phytools-rscript-executable",
            str(rscript),
            "--phytools-case",
            "simmap-er-binary-twenty-four-taxa",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    report = payload["data"]["report"]
    assert report["observations"][0]["case_id"] == "simmap-er-binary-twenty-four-taxa"


def test_parity_cli_restricts_live_phytools_sim_history_root_prior_case(
    tmp_path: Path, capsys
) -> None:
    rscript = fake_phytools_rscript(tmp_path / "fake-phytools-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "phytools-live",
            "--phytools-rscript-executable",
            str(rscript),
            "--phytools-case",
            "sim-history-binary-root-prior-example-tree",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "sim-history-binary-root-prior-example-tree"
    assert observation["bijux_summary"]["fixed_root_state"] is None
    assert observation["bijux_summary"]["root_prior_probabilities"] == {
        "0": 0.2,
        "1": 0.8,
    }
