from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from bijux_phylogenetics.parity import (
    run_geiger_parity_cases,
    write_geiger_boundary_warning_table,
    write_geiger_likelihood_policy_table,
    write_geiger_model_confidence_table,
    write_geiger_optimizer_triage_table,
    write_geiger_parameterization_registry_table,
    write_geiger_parity_observation_table,
    write_geiger_parity_summary_table,
)
from tests.support.fake_geiger_parity import fake_geiger_rscript
from tests.support.geiger_fitcontinuous_brownian_reference import (
    GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_lambda_reference import (
    GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_model_comparison_reference import (
    GEIGER_FITCONTINUOUS_MODEL_COMPARISON_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitcontinuous_ou_reference import (
    GEIGER_FITCONTINUOUS_OU_REFERENCE_PAYLOADS,
)
from tests.support.geiger_fitdiscrete_er_reference import (
    GEIGER_FITDISCRETE_ER_REFERENCE_PAYLOADS,
)

pytestmark = pytest.mark.slow


def test_write_geiger_optimizer_triage_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    triage_path = tmp_path / "geiger-optimizer-triage.tsv"

    write_geiger_optimizer_triage_table(triage_path, report)

    with triage_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert rows[0]["mismatch_type"] == "no_algorithm_mismatch"


def test_write_geiger_boundary_warning_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    boundary_path = tmp_path / "geiger-boundary-warning.tsv"

    write_geiger_boundary_warning_table(boundary_path, report)

    with boundary_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert rows[0]["reference_hit_lower_boundary"] == "True"
    assert "weak_phylogenetic_signal" in rows[0]["bijux_boundary_warning_kinds"]


def test_write_geiger_likelihood_policy_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(tmp_path / "fake-geiger-rscript")
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-bm-example-tree"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    policy_path = tmp_path / "geiger-likelihood-policy.tsv"

    write_geiger_likelihood_policy_table(policy_path, report)

    with policy_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-bm-example-tree"
    assert rows[0]["case_level_raw_log_likelihood_comparable"] == "True"


def test_write_geiger_model_confidence_table_writes_rows(tmp_path: Path) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_MODEL_COMPARISON_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-model-comparison-brownian-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    confidence_path = tmp_path / "geiger-model-confidence.tsv"

    write_geiger_model_confidence_table(confidence_path, report)

    with confidence_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 7
    assert rows[0]["case_id"] == "fitcontinuous-model-comparison-brownian-review"
    assert rows[0]["weight_basis"] == "AICc"


def test_write_geiger_parameterization_registry_table_writes_rows(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=GEIGER_FITCONTINUOUS_LAMBDA_REFERENCE_PAYLOADS,
    )
    report = run_geiger_parity_cases(
        case_ids=["fitcontinuous-lambda-weak-signal-review"],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    registry_path = tmp_path / "geiger-parameterization-registry.tsv"

    write_geiger_parameterization_registry_table(registry_path, report)

    with registry_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["case_id"] == "fitcontinuous-lambda-weak-signal-review"
    assert rows[0]["reference_parameter_name"] == "lambda"


def test_write_geiger_parity_tables_writes_summary_and_observations(
    tmp_path: Path,
) -> None:
    rscript = fake_geiger_rscript(
        tmp_path / "fake-geiger-rscript",
        reference_payloads=(
            GEIGER_FITCONTINUOUS_BROWNIAN_REFERENCE_PAYLOADS
            | GEIGER_FITCONTINUOUS_OU_REFERENCE_PAYLOADS
            | GEIGER_FITDISCRETE_ER_REFERENCE_PAYLOADS
        ),
    )
    report = run_geiger_parity_cases(
        shard_ids=[
            "geiger-fitcontinuous-model-bm",
            "geiger-fitcontinuous-model-ou",
            "geiger-fitdiscrete-model-er",
        ],
        rscript_executable=str(rscript),
        failure_root=tmp_path / "geiger-parity-failures",
    )
    summary_path = tmp_path / "geiger-parity-summary.tsv"
    observation_path = tmp_path / "geiger-parity-observations.tsv"
    triage_path = tmp_path / "geiger-optimizer-triage.tsv"

    write_geiger_parity_summary_table(summary_path, report)
    write_geiger_parity_observation_table(observation_path, report)
    write_geiger_optimizer_triage_table(triage_path, report)

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    assert summary_rows[0].startswith("function_name\tcase_count")
    assert any("geiger::fitContinuous(model='OU')" in row for row in summary_rows[1:])
    with observation_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 10
    assert rows[0]["model_name"] in {
        "BM",
        "OU",
        "ER",
    }
    optimizer_settings = json.loads(rows[0]["optimizer_settings"])
    assert "reference_control_policy" in optimizer_settings
    with triage_path.open(encoding="utf-8", newline="") as handle:
        triage_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(triage_rows) == len(rows)
    assert triage_rows[0]["mismatch_type"] in {
        "no_algorithm_mismatch",
        "parameter_surface_not_applicable",
    }
