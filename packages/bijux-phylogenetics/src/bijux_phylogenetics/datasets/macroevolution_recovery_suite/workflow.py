from __future__ import annotations

from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    ContinuousModeRecoveryPanelWorkflowReport,
    run_continuous_mode_recovery_panel_workflow,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    DiscreteModeRecoveryPanelWorkflowReport,
    run_discrete_mode_recovery_panel_workflow,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    KnownAnswerReferenceWorkflowReport,
    run_known_answer_reference_workflow,
)
from bijux_phylogenetics.simulation import (
    validate_geiger_sim_char_reference_examples,
)

from .models import (
    MacroevolutionRecoverySuiteComponentSummary,
    MacroevolutionRecoverySuiteWorkflowReport,
)
from .panel import load_macroevolution_recovery_suite_dataset


def run_macroevolution_recovery_suite_workflow() -> (
    MacroevolutionRecoverySuiteWorkflowReport
):
    """Run the governed macroevolution recovery suite across all component panels."""
    dataset = load_macroevolution_recovery_suite_dataset()
    continuous_panel_workflow = run_continuous_mode_recovery_panel_workflow()
    discrete_panel_workflow = run_discrete_mode_recovery_panel_workflow()
    known_answer_panel_workflow = run_known_answer_reference_workflow()
    sim_char_report = validate_geiger_sim_char_reference_examples()
    return MacroevolutionRecoverySuiteWorkflowReport(
        dataset=dataset,
        continuous_component=_summarize_continuous_component(
            continuous_panel_workflow,
        ),
        discrete_component=_summarize_discrete_component(
            discrete_panel_workflow,
        ),
        known_answer_component=_summarize_known_answer_component(
            known_answer_panel_workflow,
        ),
        continuous_panel_workflow=continuous_panel_workflow,
        discrete_panel_workflow=discrete_panel_workflow,
        known_answer_panel_workflow=known_answer_panel_workflow,
        sim_char_case_count=sim_char_report.case_count,
        sim_char_all_passed=sim_char_report.all_passed,
    )


def _summarize_continuous_component(
    workflow_report: ContinuousModeRecoveryPanelWorkflowReport,
) -> MacroevolutionRecoverySuiteComponentSummary:
    case_reports = workflow_report.recovery_report.case_reports
    selection_review_case_count = sum(
        1 for case in case_reports if case.selection_matches_expectation is not None
    )
    selection_match_count = sum(
        1 for case in case_reports if case.selection_matches_expectation is True
    )
    geiger_selection_match_count = sum(
        1 for case in case_reports if case.geiger_selection_matches_expectation is True
    )
    parameter_rows = [row for case in case_reports for row in case.parameter_rows]
    parameter_comparison_rows = [
        row for case in case_reports for row in case.parameter_comparison_rows
    ]
    expected_warning_case_count = sum(
        1 for case in case_reports if case.scenario.expected_warning_kinds
    )
    expected_warning_present_count = sum(
        1
        for case in case_reports
        if case.scenario.expected_warning_kinds and case.expected_warning_kinds_present
    )
    return MacroevolutionRecoverySuiteComponentSummary(
        dataset_id=workflow_report.dataset.dataset_id,
        label=workflow_report.dataset.label,
        expected_output_root=workflow_report.dataset.reference_output_root,
        case_count=workflow_report.dataset.case_count,
        taxon_count=workflow_report.dataset.taxon_count,
        selection_review_case_count=selection_review_case_count,
        selection_match_count=selection_match_count,
        geiger_selection_match_count=geiger_selection_match_count,
        governed_value_pass_count=sum(
            1 for row in parameter_rows if row.within_tolerance
        ),
        governed_value_row_count=len(parameter_rows),
        governed_comparison_row_count=len(parameter_comparison_rows),
        expected_warning_case_count=expected_warning_case_count,
        expected_warning_present_count=expected_warning_present_count,
        truth_threshold_pass_count=0,
        truth_threshold_row_count=0,
    )


def _summarize_discrete_component(
    workflow_report: DiscreteModeRecoveryPanelWorkflowReport,
) -> MacroevolutionRecoverySuiteComponentSummary:
    case_reports = workflow_report.recovery_report.case_reports
    selection_review_case_count = sum(
        1 for case in case_reports if case.selection_matches_expectation is not None
    )
    selection_match_count = sum(
        1 for case in case_reports if case.selection_matches_expectation is True
    )
    geiger_selection_match_count = sum(
        1 for case in case_reports if case.geiger_selection_matches_expectation is True
    )
    governed_rate_rows = [
        row
        for case in case_reports
        for row in case.rate_rows
        if row.tolerance is not None
    ]
    governed_rate_comparison_rows = [
        row
        for case in case_reports
        for row in case.rate_comparison_rows
        if row.tolerance is not None
    ]
    expected_warning_case_count = sum(
        1 for case in case_reports if case.scenario.expected_warning_kinds
    )
    expected_warning_present_count = sum(
        1
        for case in case_reports
        if case.scenario.expected_warning_kinds and case.expected_warning_kinds_present
    )
    return MacroevolutionRecoverySuiteComponentSummary(
        dataset_id=workflow_report.dataset.dataset_id,
        label=workflow_report.dataset.label,
        expected_output_root=workflow_report.dataset.reference_output_root,
        case_count=workflow_report.dataset.case_count,
        taxon_count=workflow_report.dataset.taxon_count,
        selection_review_case_count=selection_review_case_count,
        selection_match_count=selection_match_count,
        geiger_selection_match_count=geiger_selection_match_count,
        governed_value_pass_count=sum(
            1 for row in governed_rate_rows if row.within_tolerance is True
        ),
        governed_value_row_count=len(governed_rate_rows),
        governed_comparison_row_count=len(governed_rate_comparison_rows),
        expected_warning_case_count=expected_warning_case_count,
        expected_warning_present_count=expected_warning_present_count,
        truth_threshold_pass_count=0,
        truth_threshold_row_count=0,
    )


def _summarize_known_answer_component(
    workflow_report: KnownAnswerReferenceWorkflowReport,
) -> MacroevolutionRecoverySuiteComponentSummary:
    threshold_rows = workflow_report.threshold_evaluation_rows
    return MacroevolutionRecoverySuiteComponentSummary(
        dataset_id=workflow_report.dataset.dataset_id,
        label=workflow_report.dataset.label,
        expected_output_root=workflow_report.dataset.reference_output_root,
        case_count=len(threshold_rows),
        taxon_count=workflow_report.dataset.taxon_count,
        selection_review_case_count=0,
        selection_match_count=0,
        geiger_selection_match_count=0,
        governed_value_pass_count=0,
        governed_value_row_count=0,
        governed_comparison_row_count=0,
        expected_warning_case_count=0,
        expected_warning_present_count=0,
        truth_threshold_pass_count=sum(1 for row in threshold_rows if row.passed),
        truth_threshold_row_count=len(threshold_rows),
    )
