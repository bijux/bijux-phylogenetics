from __future__ import annotations

import csv
from pathlib import Path

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
    """Aggregate the governed macroevolution recovery evidence into one suite view."""
    dataset = load_macroevolution_recovery_suite_dataset()
    sim_char_report = validate_geiger_sim_char_reference_examples()
    return MacroevolutionRecoverySuiteWorkflowReport(
        dataset=dataset,
        continuous_component=_load_continuous_component_summary(
            dataset.continuous_panel.reference_output_root,
            label=dataset.continuous_panel.label,
        ),
        discrete_component=_load_discrete_component_summary(
            dataset.discrete_panel.reference_output_root,
            label=dataset.discrete_panel.label,
        ),
        known_answer_component=_load_known_answer_component_summary(
            dataset.known_answer_panel.reference_output_root,
            label=dataset.known_answer_panel.label,
        ),
        sim_char_case_count=sim_char_report.case_count,
        sim_char_all_passed=sim_char_report.all_passed,
    )


def _load_continuous_component_summary(
    expected_output_root: Path,
    *,
    label: str,
) -> MacroevolutionRecoverySuiteComponentSummary:
    row = _read_workflow_summary_row(expected_output_root / "workflow-summary.tsv")
    return MacroevolutionRecoverySuiteComponentSummary(
        dataset_id=row["dataset_id"],
        label=label,
        expected_output_root=expected_output_root,
        case_count=int(row["case_count"]),
        taxon_count=int(row["taxon_count"]),
        selection_review_case_count=int(row["selection_review_case_count"]),
        selection_match_count=int(row["selection_match_count"]),
        geiger_selection_match_count=int(row["geiger_selection_match_count"]),
        governed_value_pass_count=int(row["parameter_pass_count"]),
        governed_value_row_count=int(row["parameter_row_count"]),
        governed_comparison_row_count=int(row["parameter_comparison_row_count"]),
        expected_warning_case_count=int(row["expected_warning_case_count"]),
        expected_warning_present_count=int(row["expected_warning_present_count"]),
        truth_threshold_pass_count=0,
        truth_threshold_row_count=0,
    )


def _load_discrete_component_summary(
    expected_output_root: Path,
    *,
    label: str,
) -> MacroevolutionRecoverySuiteComponentSummary:
    row = _read_workflow_summary_row(expected_output_root / "workflow-summary.tsv")
    return MacroevolutionRecoverySuiteComponentSummary(
        dataset_id=row["dataset_id"],
        label=label,
        expected_output_root=expected_output_root,
        case_count=int(row["case_count"]),
        taxon_count=int(row["taxon_count"]),
        selection_review_case_count=int(row["selection_review_case_count"]),
        selection_match_count=int(row["selection_match_count"]),
        geiger_selection_match_count=int(row["geiger_selection_match_count"]),
        governed_value_pass_count=int(row["rate_pass_count"]),
        governed_value_row_count=int(row["governed_rate_row_count"]),
        governed_comparison_row_count=int(row["governed_rate_comparison_row_count"]),
        expected_warning_case_count=int(row["expected_warning_case_count"]),
        expected_warning_present_count=int(row["expected_warning_present_count"]),
        truth_threshold_pass_count=0,
        truth_threshold_row_count=0,
    )


def _load_known_answer_component_summary(
    expected_output_root: Path,
    *,
    label: str,
) -> MacroevolutionRecoverySuiteComponentSummary:
    row = _read_workflow_summary_row(expected_output_root / "workflow-summary.tsv")
    return MacroevolutionRecoverySuiteComponentSummary(
        dataset_id=row["dataset_id"],
        label=label,
        expected_output_root=expected_output_root,
        case_count=int(row["threshold_row_count"]),
        taxon_count=int(row["taxon_count"]),
        selection_review_case_count=0,
        selection_match_count=0,
        geiger_selection_match_count=0,
        governed_value_pass_count=0,
        governed_value_row_count=0,
        governed_comparison_row_count=0,
        expected_warning_case_count=0,
        expected_warning_present_count=0,
        truth_threshold_pass_count=int(row["threshold_pass_count"]),
        truth_threshold_row_count=int(row["threshold_row_count"]),
    )


def _read_workflow_summary_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader, None)
    if row is None:
        raise ValueError(f"workflow summary is empty: {path}")
    return {key: value for key, value in row.items() if key is not None and value is not None}
