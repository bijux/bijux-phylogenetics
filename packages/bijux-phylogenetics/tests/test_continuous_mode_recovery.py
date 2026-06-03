from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.comparative.continuous_mode_recovery as continuous_mode_recovery_api
from bijux_phylogenetics.comparative.continuous_mode_recovery import (
    geiger_fitcontinuous_recovery_reference_payload,
    run_continuous_mode_recovery,
    write_continuous_mode_recovery_execution_table,
    write_continuous_mode_recovery_model_choice_table,
    write_continuous_mode_recovery_parameter_comparison_table,
    write_continuous_mode_recovery_parameter_table,
    write_continuous_mode_recovery_summary_table,
    write_continuous_mode_recovery_warning_table,
    write_geiger_fitcontinuous_recovery_reference_payload_table,
)
from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    load_continuous_mode_recovery_panel_dataset,
    run_continuous_mode_recovery_panel_workflow,
)


@pytest.mark.slow
def test_run_continuous_mode_recovery_governs_packaged_bijux_and_geiger_review_cases(
    tmp_path: Path,
) -> None:
    dataset = load_continuous_mode_recovery_panel_dataset()
    workflow_report = run_continuous_mode_recovery_panel_workflow()
    scenarios = [case.scenario for case in workflow_report.recovery_report.case_reports]

    report = run_continuous_mode_recovery(
        dataset.default_tree_path, scenarios, artifacts_root=tmp_path / "artifacts"
    )

    case_by_id = {case.scenario.case_id: case for case in report.case_reports}
    assert set(case_by_id) == {
        "brownian-sigma-recovery",
        "ou-parameter-recovery",
        "early-burst-rate-recovery",
        "weak-ou-identifiability",
        "lambda-transformed-branch-review",
        "kappa-transformed-branch-review",
        "delta-transformed-branch-review",
    }
    assert case_by_id["brownian-sigma-recovery"].selected_model == "brownian"
    assert case_by_id["brownian-sigma-recovery"].geiger_selected_model == "brownian"
    assert case_by_id["ou-parameter-recovery"].selected_model == "ornstein-uhlenbeck"
    assert case_by_id["ou-parameter-recovery"].geiger_selected_model == "brownian"
    assert case_by_id["early-burst-rate-recovery"].selected_model == "early-burst"
    assert (
        case_by_id["early-burst-rate-recovery"].geiger_selected_model == "early-burst"
    )
    assert case_by_id["weak-ou-identifiability"].selected_model == "brownian"
    assert case_by_id["weak-ou-identifiability"].geiger_selected_model == "brownian"
    assert case_by_id["lambda-transformed-branch-review"].selected_model == (
        "ornstein-uhlenbeck"
    )
    assert case_by_id["lambda-transformed-branch-review"].geiger_selected_model == (
        "ornstein-uhlenbeck"
    )
    assert case_by_id["kappa-transformed-branch-review"].selected_model == (
        "pagel-kappa"
    )
    assert case_by_id["kappa-transformed-branch-review"].geiger_selected_model == (
        "pagel-kappa"
    )
    assert case_by_id["delta-transformed-branch-review"].selected_model == (
        "pagel-kappa"
    )
    assert case_by_id["delta-transformed-branch-review"].geiger_selected_model == (
        "pagel-kappa"
    )
    ou_theta_rows = [
        row
        for row in case_by_id["ou-parameter-recovery"].parameter_rows
        if row.parameter == "theta"
    ]
    assert len(ou_theta_rows) == 2
    assert {row.recovery_engine for row in ou_theta_rows} == {"bijux", "geiger"}
    assert all(row.within_tolerance for row in ou_theta_rows)
    assert all(row.true_value == 2.8 for row in ou_theta_rows)
    assert case_by_id["weak-ou-identifiability"].parameter_rows == []
    assert case_by_id["weak-ou-identifiability"].expected_warning_kinds_present is True
    assert (
        case_by_id["lambda-transformed-branch-review"].expected_warning_kinds_present
        is True
    )
    assert (
        case_by_id["kappa-transformed-branch-review"].expected_warning_kinds_present
        is True
    )
    assert (
        case_by_id["delta-transformed-branch-review"].expected_warning_kinds_present
        is True
    )
    assert all(
        row.case_id != "weak-ou-identifiability"
        for row in case_by_id["weak-ou-identifiability"].parameter_comparison_rows
    )
    assert case_by_id["brownian-sigma-recovery"].traits_path is not None


@pytest.mark.slow
def test_continuous_mode_recovery_writers_emit_paired_benchmark_ledgers(
    tmp_path: Path,
) -> None:
    report = run_continuous_mode_recovery_panel_workflow().recovery_report

    summary_path = write_continuous_mode_recovery_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    parameter_path = write_continuous_mode_recovery_parameter_table(
        tmp_path / "parameter-recovery.tsv",
        report,
    )
    parameter_comparison_path = (
        write_continuous_mode_recovery_parameter_comparison_table(
            tmp_path / "parameter-comparison.tsv",
            report,
        )
    )
    model_choice_path = write_continuous_mode_recovery_model_choice_table(
        tmp_path / "model-choice.tsv",
        report,
    )
    execution_path = write_continuous_mode_recovery_execution_table(
        tmp_path / "execution-review.tsv",
        report,
    )
    warning_path = write_continuous_mode_recovery_warning_table(
        tmp_path / "warning-review.tsv",
        report,
    )
    geiger_payload_path = write_geiger_fitcontinuous_recovery_reference_payload_table(
        tmp_path / "geiger-reference.tsv",
        report,
    )

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    parameter_rows = parameter_path.read_text(encoding="utf-8").splitlines()
    parameter_comparison_rows = parameter_comparison_path.read_text(
        encoding="utf-8"
    ).splitlines()
    model_choice_rows = model_choice_path.read_text(encoding="utf-8").splitlines()
    execution_rows = execution_path.read_text(encoding="utf-8").splitlines()
    warning_rows = warning_path.read_text(encoding="utf-8").splitlines()
    geiger_payload_rows = geiger_payload_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith("case_id\tlabel\ttree_path\tgenerating_model")
    assert any(
        row.startswith("lambda-transformed-branch-review\t") for row in summary_rows[1:]
    )
    assert parameter_rows[0].startswith(
        "case_id\tgenerating_model\trecovery_engine\tfitted_model"
    )
    assert any(
        row.startswith(
            "ou-parameter-recovery\tornstein-uhlenbeck\tgeiger\tornstein-uhlenbeck\t"
        )
        for row in parameter_rows[1:]
    )
    assert parameter_comparison_rows[0].startswith(
        "case_id\tgenerating_model\tparameter\ttrue_value"
    )
    assert any(
        row.startswith("ou-parameter-recovery\tornstein-uhlenbeck\ttheta\t2.8\t")
        for row in parameter_comparison_rows[1:]
    )
    assert any(
        row.startswith("early-burst-rate-recovery\tearly-burst\trate_change\t")
        for row in parameter_comparison_rows[1:]
    )
    assert model_choice_rows[0].startswith(
        "case_id\tgenerating_model\trecovery_engine\texpected_selected_model"
    )
    assert any(
        row.startswith(
            "delta-transformed-branch-review\tpagel-delta\tgeiger\t\tpagel-kappa\t"
        )
        for row in model_choice_rows[1:]
    )
    assert execution_rows[0].startswith(
        "case_id\trecovery_engine\toperation\tfitted_model"
    )
    assert any(
        row.startswith(
            "kappa-transformed-branch-review\tgeiger\tmodel-comparison\tcandidate-set\tok\tpagel-kappa\t"
        )
        for row in execution_rows[1:]
    )
    assert warning_rows[0] == "case_id\trecovery_engine\tfitted_model\tkind\tmessage"
    assert any(
        row.startswith(
            "weak-ou-identifiability\tbijux\tornstein-uhlenbeck\tflat_likelihood\t"
        )
        for row in warning_rows[1:]
    )
    assert (
        geiger_payload_rows[0] == "case_id\tfit_summary_json\tcomparison_summary_json"
    )
    assert any(
        row.startswith("brownian-sigma-recovery\t") for row in geiger_payload_rows[1:]
    )


def test_geiger_fitcontinuous_recovery_reference_payload_exposes_governed_case_summary() -> (
    None
):
    payload = geiger_fitcontinuous_recovery_reference_payload("brownian-sigma-recovery")
    assert payload["fit_summary"]["model_name"] == "BM"
    assert payload["comparison_summary"]["selected_model"] == "brownian"


def test_public_runtime_exports_include_continuous_mode_recovery_surface() -> None:
    assert (
        continuous_mode_recovery_api.run_continuous_mode_recovery
        is run_continuous_mode_recovery
    )
    assert (
        continuous_mode_recovery_api.write_continuous_mode_recovery_summary_table
        is write_continuous_mode_recovery_summary_table
    )
    assert (
        continuous_mode_recovery_api.write_continuous_mode_recovery_parameter_table
        is write_continuous_mode_recovery_parameter_table
    )
    assert (
        continuous_mode_recovery_api.write_continuous_mode_recovery_parameter_comparison_table
        is write_continuous_mode_recovery_parameter_comparison_table
    )
    assert (
        continuous_mode_recovery_api.write_continuous_mode_recovery_model_choice_table
        is write_continuous_mode_recovery_model_choice_table
    )
    assert (
        continuous_mode_recovery_api.write_continuous_mode_recovery_execution_table
        is write_continuous_mode_recovery_execution_table
    )
    assert (
        continuous_mode_recovery_api.write_continuous_mode_recovery_warning_table
        is write_continuous_mode_recovery_warning_table
    )
    assert (
        continuous_mode_recovery_api.geiger_fitcontinuous_recovery_reference_payload
        is geiger_fitcontinuous_recovery_reference_payload
    )
    assert (
        continuous_mode_recovery_api.write_geiger_fitcontinuous_recovery_reference_payload_table
        is write_geiger_fitcontinuous_recovery_reference_payload_table
    )
