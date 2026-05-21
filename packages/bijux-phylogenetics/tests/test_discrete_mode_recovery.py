from __future__ import annotations

from pathlib import Path

import pytest

import bijux_phylogenetics.comparative.discrete_mode_recovery as discrete_mode_recovery_api
from bijux_phylogenetics.comparative.discrete_mode_recovery import (
    geiger_fitdiscrete_recovery_reference_payload,
    run_discrete_mode_recovery,
    write_discrete_mode_recovery_execution_table,
    write_discrete_mode_recovery_model_choice_table,
    write_discrete_mode_recovery_parameter_comparison_table,
    write_discrete_mode_recovery_parameter_table,
    write_discrete_mode_recovery_rate_comparison_table,
    write_discrete_mode_recovery_rate_table,
    write_discrete_mode_recovery_summary_table,
    write_discrete_mode_recovery_warning_table,
    write_geiger_fitdiscrete_recovery_reference_payload_table,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    load_discrete_mode_recovery_panel_dataset,
    run_discrete_mode_recovery_panel_workflow,
)


@pytest.mark.slow
def test_run_discrete_mode_recovery_governs_packaged_bijux_and_geiger_review_cases(
    tmp_path: Path,
) -> None:
    dataset = load_discrete_mode_recovery_panel_dataset()
    workflow_report = run_discrete_mode_recovery_panel_workflow()
    scenarios = [case.scenario for case in workflow_report.recovery_report.case_reports]

    report = run_discrete_mode_recovery(
        dataset.default_tree_path,
        scenarios,
        artifacts_root=tmp_path / "artifacts",
    )

    case_by_id = {case.scenario.case_id: case for case in report.case_reports}
    assert set(case_by_id) == {
        "er-three-state-rate-recovery",
        "sym-three-state-rate-recovery",
        "ard-three-state-weak-identification-review",
        "ard-five-state-overparameterized-review",
    }
    assert case_by_id["er-three-state-rate-recovery"].selected_model == "equal-rates"
    assert case_by_id["sym-three-state-rate-recovery"].selected_model == "symmetric"
    assert (
        case_by_id[
            "ard-three-state-weak-identification-review"
        ].selection_matches_expectation
        is None
    )
    assert (
        case_by_id[
            "ard-three-state-weak-identification-review"
        ].expected_warning_kinds_present
        is True
    )
    assert any(
        row.kind == "equal_rates_preferred"
        for row in case_by_id["ard-three-state-weak-identification-review"].warning_rows
        if row.recovery_engine == "bijux"
    )
    assert (
        case_by_id[
            "ard-five-state-overparameterized-review"
        ].overparameterized_review_matches_expectation
        is True
    )
    assert (
        case_by_id[
            "ard-five-state-overparameterized-review"
        ].expected_warning_kinds_present
        is True
    )
    assert any(
        row.kind == "overparameterized"
        for row in case_by_id["ard-five-state-overparameterized-review"].warning_rows
        if row.recovery_engine == "bijux"
    )
    assert all(
        row.tolerance is not None
        for row in case_by_id["er-three-state-rate-recovery"].rate_rows
    )
    assert any(
        row.tolerance is None
        for row in case_by_id["ard-three-state-weak-identification-review"].rate_rows
    )
    assert case_by_id["er-three-state-rate-recovery"].traits_path is not None


@pytest.mark.slow
def test_discrete_mode_recovery_writers_emit_paired_benchmark_ledgers(
    tmp_path: Path,
) -> None:
    report = run_discrete_mode_recovery_panel_workflow().recovery_report

    summary_path = write_discrete_mode_recovery_summary_table(
        tmp_path / "summary.tsv",
        report,
    )
    parameter_path = write_discrete_mode_recovery_parameter_table(
        tmp_path / "parameter-recovery.tsv",
        report,
    )
    parameter_comparison_path = write_discrete_mode_recovery_parameter_comparison_table(
        tmp_path / "parameter-comparison.tsv",
        report,
    )
    rate_path = write_discrete_mode_recovery_rate_table(
        tmp_path / "rate-recovery.tsv",
        report,
    )
    rate_comparison_path = write_discrete_mode_recovery_rate_comparison_table(
        tmp_path / "rate-comparison.tsv",
        report,
    )
    model_choice_path = write_discrete_mode_recovery_model_choice_table(
        tmp_path / "model-choice.tsv",
        report,
    )
    execution_path = write_discrete_mode_recovery_execution_table(
        tmp_path / "execution-review.tsv",
        report,
    )
    warning_path = write_discrete_mode_recovery_warning_table(
        tmp_path / "warning-review.tsv",
        report,
    )
    geiger_payload_path = write_geiger_fitdiscrete_recovery_reference_payload_table(
        tmp_path / "geiger-reference.tsv",
        report,
    )

    summary_rows = summary_path.read_text(encoding="utf-8").splitlines()
    parameter_rows = parameter_path.read_text(encoding="utf-8").splitlines()
    parameter_comparison_rows = parameter_comparison_path.read_text(
        encoding="utf-8"
    ).splitlines()
    rate_rows = rate_path.read_text(encoding="utf-8").splitlines()
    rate_comparison_rows = rate_comparison_path.read_text(encoding="utf-8").splitlines()
    model_choice_rows = model_choice_path.read_text(encoding="utf-8").splitlines()
    execution_rows = execution_path.read_text(encoding="utf-8").splitlines()
    warning_rows = warning_path.read_text(encoding="utf-8").splitlines()
    geiger_payload_rows = geiger_payload_path.read_text(encoding="utf-8").splitlines()

    assert summary_rows[0].startswith("case_id\tlabel\ttree_path\tgenerating_model")
    assert any(
        row.startswith("ard-three-state-weak-identification-review\t")
        for row in summary_rows[1:]
    )
    assert parameter_rows[0].startswith(
        "case_id\tgenerating_model\ttransform\trecovery_engine\tfitted_model"
    )
    assert parameter_comparison_rows[0].startswith(
        "case_id\tgenerating_model\ttransform\tparameter\ttrue_value"
    )
    assert rate_rows[0].startswith(
        "case_id\tgenerating_model\ttransform\trecovery_engine\tfitted_model"
    )
    assert any(
        row.startswith(
            (
                "sym-three-state-rate-recovery\tsymmetric\t\tgeiger\tall-rates-different\t",
                "sym-three-state-rate-recovery\tsymmetric\t\tgeiger\tsymmetric\t",
            )
        )
        for row in rate_rows[1:]
    )
    assert rate_comparison_rows[0].startswith(
        "case_id\tgenerating_model\ttransform\tsource_state\ttarget_state\ttrue_rate"
    )
    assert any(
        row.startswith("er-three-state-rate-recovery\tequal-rates\t\tcentral\tnorth\t")
        for row in rate_comparison_rows[1:]
    )
    assert model_choice_rows[0].startswith(
        "case_id\tgenerating_model\ttransform\trecovery_engine\texpected_selected_model"
    )
    assert any(
        row.startswith(
            "ard-five-state-overparameterized-review\tall-rates-different\t\tgeiger\t\t"
        )
        for row in model_choice_rows[1:]
    )
    assert execution_rows[0].startswith(
        "case_id\trecovery_engine\toperation\tfitted_model"
    )
    assert any(
        row.startswith(
            "ard-three-state-weak-identification-review\tbijux\tfit\tall-rates-different\tok\t"
        )
        for row in execution_rows[1:]
    )
    assert warning_rows[0] == "case_id\trecovery_engine\tfitted_model\tkind\tmessage"
    assert any(
        row.startswith(
            "ard-three-state-weak-identification-review\tbijux\tall-rates-different\tequal_rates_preferred\t"
        )
        for row in warning_rows[1:]
    )
    assert (
        geiger_payload_rows[0] == "case_id\tfit_summary_json\tcomparison_summary_json"
    )
    assert any(
        row.startswith("er-three-state-rate-recovery\t")
        for row in geiger_payload_rows[1:]
    )


def test_geiger_fitdiscrete_recovery_reference_payload_exposes_governed_case_summary() -> (
    None
):
    payload = geiger_fitdiscrete_recovery_reference_payload(
        "er-three-state-rate-recovery"
    )
    assert payload["fit_summary"]["model_name"] == "ER"
    assert payload["comparison_summary"]["selected_model"] in {
        "equal-rates",
        "symmetric",
        "all-rates-different",
    }


def test_public_runtime_exports_include_discrete_mode_recovery_surface() -> None:
    assert (
        discrete_mode_recovery_api.run_discrete_mode_recovery
        is run_discrete_mode_recovery
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_summary_table
        is write_discrete_mode_recovery_summary_table
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_rate_table
        is write_discrete_mode_recovery_rate_table
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_rate_comparison_table
        is write_discrete_mode_recovery_rate_comparison_table
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_model_choice_table
        is write_discrete_mode_recovery_model_choice_table
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_parameter_table
        is write_discrete_mode_recovery_parameter_table
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_parameter_comparison_table
        is write_discrete_mode_recovery_parameter_comparison_table
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_execution_table
        is write_discrete_mode_recovery_execution_table
    )
    assert (
        discrete_mode_recovery_api.write_discrete_mode_recovery_warning_table
        is write_discrete_mode_recovery_warning_table
    )
    assert (
        discrete_mode_recovery_api.geiger_fitdiscrete_recovery_reference_payload
        is geiger_fitdiscrete_recovery_reference_payload
    )
    assert (
        discrete_mode_recovery_api.write_geiger_fitdiscrete_recovery_reference_payload_table
        is write_geiger_fitdiscrete_recovery_reference_payload_table
    )
