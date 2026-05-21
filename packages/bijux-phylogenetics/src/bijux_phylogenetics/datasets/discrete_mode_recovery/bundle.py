from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.comparative.discrete_mode_recovery import (
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
from bijux_phylogenetics.simulation import write_discrete_trait_table

from .models import (
    DiscreteModeRecoveryPanelWorkflowBundle,
    DiscreteModeRecoveryPanelWorkflowReport,
)


def write_discrete_mode_recovery_panel_workflow_bundle(
    output_root: Path,
    report: DiscreteModeRecoveryPanelWorkflowReport,
) -> DiscreteModeRecoveryPanelWorkflowBundle:
    """Write the reviewer-facing recovery outputs for the packaged panel."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    selection_review_case_count = sum(
        1
        for case in report.recovery_report.case_reports
        if case.selection_matches_expectation is not None
    )
    selection_match_count = sum(
        1
        for case in report.recovery_report.case_reports
        if case.selection_matches_expectation is True
    )
    geiger_selection_match_count = sum(
        1
        for case in report.recovery_report.case_reports
        if case.geiger_selection_matches_expectation is True
    )
    parameter_rows = [
        row
        for case in report.recovery_report.case_reports
        for row in case.parameter_rows
    ]
    governed_parameter_rows = [
        row for row in parameter_rows if row.tolerance is not None
    ]
    parameter_comparison_rows = [
        row
        for case in report.recovery_report.case_reports
        for row in case.parameter_comparison_rows
    ]
    parameter_pass_count = sum(
        1 for row in governed_parameter_rows if row.within_tolerance is True
    )
    rate_rows = [
        row for case in report.recovery_report.case_reports for row in case.rate_rows
    ]
    governed_rate_rows = [row for row in rate_rows if row.tolerance is not None]
    rate_comparison_rows = [
        row
        for case in report.recovery_report.case_reports
        for row in case.rate_comparison_rows
    ]
    governed_rate_comparison_rows = [
        row for row in rate_comparison_rows if row.tolerance is not None
    ]
    rate_pass_count = sum(
        1 for row in governed_rate_rows if row.within_tolerance is True
    )
    expected_warning_case_count = sum(
        1
        for case in report.recovery_report.case_reports
        if case.scenario.expected_warning_kinds
    )
    expected_warning_present_count = sum(
        1
        for case in report.recovery_report.case_reports
        if case.scenario.expected_warning_kinds and case.expected_warning_kinds_present
    )
    parameter_closer_to_truth_count_bijux = sum(
        1 for row in parameter_comparison_rows if row.closer_engine == "bijux"
    )
    parameter_closer_to_truth_count_geiger = sum(
        1 for row in parameter_comparison_rows if row.closer_engine == "geiger"
    )
    rate_closer_to_truth_count_bijux = sum(
        1 for row in rate_comparison_rows if row.closer_engine == "bijux"
    )
    rate_closer_to_truth_count_geiger = sum(
        1 for row in rate_comparison_rows if row.closer_engine == "geiger"
    )
    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report=report,
        selection_review_case_count=selection_review_case_count,
        selection_match_count=selection_match_count,
        geiger_selection_match_count=geiger_selection_match_count,
        parameter_pass_count=parameter_pass_count,
        governed_parameter_row_count=len(governed_parameter_rows),
        parameter_row_count=len(parameter_rows),
        parameter_comparison_row_count=len(parameter_comparison_rows),
        parameter_closer_to_truth_count_bijux=parameter_closer_to_truth_count_bijux,
        parameter_closer_to_truth_count_geiger=parameter_closer_to_truth_count_geiger,
        rate_pass_count=rate_pass_count,
        governed_rate_row_count=len(governed_rate_rows),
        rate_row_count=len(rate_rows),
        governed_rate_comparison_row_count=len(governed_rate_comparison_rows),
        rate_comparison_row_count=len(rate_comparison_rows),
        rate_closer_to_truth_count_bijux=rate_closer_to_truth_count_bijux,
        rate_closer_to_truth_count_geiger=rate_closer_to_truth_count_geiger,
        expected_warning_case_count=expected_warning_case_count,
        expected_warning_present_count=expected_warning_present_count,
    )
    recovery_summary_path = write_discrete_mode_recovery_summary_table(
        output_root / "recovery-summary.tsv",
        report.recovery_report,
    )
    parameter_recovery_path = write_discrete_mode_recovery_parameter_table(
        output_root / "parameter-recovery.tsv",
        report.recovery_report,
    )
    parameter_comparison_path = write_discrete_mode_recovery_parameter_comparison_table(
        output_root / "parameter-comparison.tsv",
        report.recovery_report,
    )
    rate_recovery_path = write_discrete_mode_recovery_rate_table(
        output_root / "rate-recovery.tsv",
        report.recovery_report,
    )
    rate_comparison_path = write_discrete_mode_recovery_rate_comparison_table(
        output_root / "rate-comparison.tsv",
        report.recovery_report,
    )
    model_choice_path = write_discrete_mode_recovery_model_choice_table(
        output_root / "model-choice.tsv",
        report.recovery_report,
    )
    execution_review_path = write_discrete_mode_recovery_execution_table(
        output_root / "execution-review.tsv",
        report.recovery_report,
    )
    warning_review_path = write_discrete_mode_recovery_warning_table(
        output_root / "warning-review.tsv",
        report.recovery_report,
    )
    geiger_reference_path = write_geiger_fitdiscrete_recovery_reference_payload_table(
        output_root / "geiger-reference.tsv",
        report.recovery_report,
    )
    simulated_traits_root = output_root / "simulated-traits"
    simulated_traits_root.mkdir(parents=True, exist_ok=True)
    for case in report.recovery_report.case_reports:
        write_discrete_trait_table(
            simulated_traits_root / f"{case.scenario.case_id}.tsv",
            case.simulation,
        )
    return DiscreteModeRecoveryPanelWorkflowBundle(
        output_root=output_root,
        selection_review_case_count=selection_review_case_count,
        selection_match_count=selection_match_count,
        geiger_selection_match_count=geiger_selection_match_count,
        parameter_pass_count=parameter_pass_count,
        governed_parameter_row_count=len(governed_parameter_rows),
        parameter_row_count=len(parameter_rows),
        parameter_comparison_row_count=len(parameter_comparison_rows),
        parameter_closer_to_truth_count_bijux=parameter_closer_to_truth_count_bijux,
        parameter_closer_to_truth_count_geiger=parameter_closer_to_truth_count_geiger,
        rate_pass_count=rate_pass_count,
        governed_rate_row_count=len(governed_rate_rows),
        rate_row_count=len(rate_rows),
        governed_rate_comparison_row_count=len(governed_rate_comparison_rows),
        rate_comparison_row_count=len(rate_comparison_rows),
        rate_closer_to_truth_count_bijux=rate_closer_to_truth_count_bijux,
        rate_closer_to_truth_count_geiger=rate_closer_to_truth_count_geiger,
        expected_warning_case_count=expected_warning_case_count,
        expected_warning_present_count=expected_warning_present_count,
        workflow_summary_path=workflow_summary_path,
        recovery_summary_path=recovery_summary_path,
        parameter_recovery_path=parameter_recovery_path,
        parameter_comparison_path=parameter_comparison_path,
        rate_recovery_path=rate_recovery_path,
        rate_comparison_path=rate_comparison_path,
        model_choice_path=model_choice_path,
        execution_review_path=execution_review_path,
        warning_review_path=warning_review_path,
        geiger_reference_path=geiger_reference_path,
        simulated_traits_root=simulated_traits_root,
    )


def _write_workflow_summary_table(
    path: Path,
    *,
    report: DiscreteModeRecoveryPanelWorkflowReport,
    selection_review_case_count: int,
    selection_match_count: int,
    geiger_selection_match_count: int,
    parameter_pass_count: int,
    governed_parameter_row_count: int,
    parameter_row_count: int,
    parameter_comparison_row_count: int,
    parameter_closer_to_truth_count_bijux: int,
    parameter_closer_to_truth_count_geiger: int,
    rate_pass_count: int,
    governed_rate_row_count: int,
    rate_row_count: int,
    governed_rate_comparison_row_count: int,
    rate_comparison_row_count: int,
    rate_closer_to_truth_count_bijux: int,
    rate_closer_to_truth_count_geiger: int,
    expected_warning_case_count: int,
    expected_warning_present_count: int,
) -> Path:
    rows = [
        "\t".join(
            [
                "dataset_id",
                "taxon_count",
                "tree_count",
                "case_count",
                "selection_review_case_count",
                "selection_match_count",
                "geiger_selection_match_count",
                "parameter_pass_count",
                "governed_parameter_row_count",
                "parameter_row_count",
                "parameter_comparison_row_count",
                "parameter_closer_to_truth_count_bijux",
                "parameter_closer_to_truth_count_geiger",
                "rate_pass_count",
                "governed_rate_row_count",
                "rate_row_count",
                "governed_rate_comparison_row_count",
                "rate_comparison_row_count",
                "rate_closer_to_truth_count_bijux",
                "rate_closer_to_truth_count_geiger",
                "expected_warning_case_count",
                "expected_warning_present_count",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.taxon_count),
                str(report.dataset.tree_count),
                str(report.dataset.case_count),
                str(selection_review_case_count),
                str(selection_match_count),
                str(geiger_selection_match_count),
                str(parameter_pass_count),
                str(governed_parameter_row_count),
                str(parameter_row_count),
                str(parameter_comparison_row_count),
                str(parameter_closer_to_truth_count_bijux),
                str(parameter_closer_to_truth_count_geiger),
                str(rate_pass_count),
                str(governed_rate_row_count),
                str(rate_row_count),
                str(governed_rate_comparison_row_count),
                str(rate_comparison_row_count),
                str(rate_closer_to_truth_count_bijux),
                str(rate_closer_to_truth_count_geiger),
                str(expected_warning_case_count),
                str(expected_warning_present_count),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path
