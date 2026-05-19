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
    triage_path = tmp_path / "geiger-optimizer-triage.tsv"

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
            "--optimizer-triage-out",
            str(triage_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["reference_source"] == "geiger-live"
    assert payload["metrics"]["case_count"] == 48
    assert payload["metrics"]["function_count"] == 17
    assert payload["metrics"]["skipped_case_count"] == 0
    assert summary_path.exists()
    assert observation_path.exists()
    assert triage_path.exists()
    assert payload["data"]["optimizer_triage_table"] == str(triage_path)
    assert len(payload["data"]["report"]["optimizer_triage_rows"]) == payload["metrics"][
        "case_count"
    ]


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
            "fitcontinuous-lambda-weak-signal-review",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert observation["model_name"] == "lambda"


def test_parity_cli_writes_geiger_optimizer_triage_table_for_single_case(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    triage_path = tmp_path / "geiger-optimizer-triage.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitcontinuous-lambda-weak-signal-review",
            "--optimizer-triage-out",
            str(triage_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert triage_path.exists()
    assert payload["data"]["optimizer_triage_table"] == str(triage_path)
    triage_row = payload["data"]["report"]["optimizer_triage_rows"][0]
    assert triage_row["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert triage_row["mismatch_type"] == "no_algorithm_mismatch"


def test_parity_cli_writes_geiger_boundary_warning_table_for_single_case(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    boundary_path = tmp_path / "geiger-boundary-warning.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitcontinuous-lambda-weak-signal-review",
            "--boundary-warning-out",
            str(boundary_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert boundary_path.exists()
    assert payload["data"]["boundary_warning_table"] == str(boundary_path)
    row = payload["data"]["report"]["boundary_warning_rows"][0]
    assert row["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert row["reference_hit_lower_boundary"] is True
    assert row["bijux_stable_conclusion_supported"] is False


def test_parity_cli_writes_geiger_model_confidence_table_for_single_case(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    confidence_path = tmp_path / "geiger-model-confidence.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitcontinuous-model-comparison-brownian-review",
            "--model-confidence-out",
            str(confidence_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert confidence_path.exists()
    assert payload["data"]["model_confidence_table"] == str(confidence_path)
    row = payload["data"]["report"]["model_confidence_rows"][0]
    assert row["case_id"] == "fitcontinuous-model-comparison-brownian-review"
    assert row["weight_basis"] == "AICc"
    assert row["reference_best_model"] == "brownian"


def test_parity_cli_writes_geiger_parameterization_registry_for_single_case(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    registry_path = tmp_path / "geiger-parameterization-registry.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitcontinuous-eb-early-burst-rate-recovery",
            "--parameterization-registry-out",
            str(registry_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert registry_path.exists()
    assert payload["data"]["parameterization_registry_table"] == str(registry_path)
    registry_row = payload["data"]["report"]["parameterization_registry_rows"][0]
    assert registry_row["case_id"] == "fitcontinuous-eb-early-burst-rate-recovery"
    assert registry_row["expected_divergence_kind"] == (
        "continuous-early-burst-sign-and-bound-convention"
    )


def test_parity_cli_writes_geiger_likelihood_policy_table_for_single_case(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    policy_path = tmp_path / "geiger-likelihood-policy.tsv"

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitcontinuous-bm-example-tree",
            "--likelihood-policy-out",
            str(policy_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert policy_path.exists()
    assert payload["data"]["likelihood_policy_table"] == str(policy_path)
    row = payload["data"]["report"]["likelihood_policy_rows"][0]
    assert row["case_id"] == "fitcontinuous-bm-example-tree"
    assert row["case_level_raw_log_likelihood_comparable"] is True


def test_parity_cli_restricts_live_geiger_discrete_cases(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitdiscrete-er-binary-twenty-four-taxa",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitdiscrete-er-binary-twenty-four-taxa"
    assert observation["model_name"] == "ER"


def test_parity_cli_restricts_live_geiger_symmetric_discrete_cases(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitdiscrete-sym-three-state-twenty-four-taxa",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitdiscrete-sym-three-state-twenty-four-taxa"
    assert observation["model_name"] == "SYM"


def test_parity_cli_restricts_live_geiger_lambda_discrete_cases(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitdiscrete-lambda-weak-signal-review",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitdiscrete-lambda-weak-signal-review"
    assert observation["model_name"] == "ER"


def test_parity_cli_restricts_live_geiger_kappa_discrete_cases(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitdiscrete-kappa-strong-signal-review",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitdiscrete-kappa-strong-signal-review"
    assert observation["model_name"] == "ER"


def test_parity_cli_restricts_live_geiger_delta_discrete_cases(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitdiscrete-delta-late-change-boundary-review",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitdiscrete-delta-late-change-boundary-review"
    assert observation["model_name"] == "ER"


def test_parity_cli_restricts_live_geiger_early_burst_discrete_cases(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitdiscrete-early-burst-early-change-review",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitdiscrete-early-burst-early-change-review"
    assert observation["model_name"] == "ER"


def test_parity_cli_restricts_live_geiger_ard_discrete_cases(
    tmp_path: Path,
    capsys,
) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")

    exit_code = main(
        [
            "parity",
            "--reference-source",
            "geiger-live",
            "--geiger-rscript-executable",
            str(rscript),
            "--geiger-case",
            "fitdiscrete-ard-binary-twenty-four-taxa",
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["case_count"] == 1
    observation = payload["data"]["report"]["observations"][0]
    assert observation["case_id"] == "fitdiscrete-ard-binary-twenty-four-taxa"
    assert observation["model_name"] == "ARD"
