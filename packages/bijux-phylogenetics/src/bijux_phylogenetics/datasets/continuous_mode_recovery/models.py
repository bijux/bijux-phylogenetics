from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.continuous_mode_recovery import (
    ContinuousModeRecoveryReport,
)


@dataclass(slots=True)
class ContinuousModeRecoveryPanelDataset:
    """Packaged deterministic recovery panel with governed strong and weak cases."""

    dataset_id: str
    label: str
    dataset_root: Path
    default_tree_path: Path
    reference_tree_paths: list[Path]
    simulation_cases_path: Path
    reference_output_root: Path
    taxon_count: int
    tree_count: int
    case_count: int
    source_summary: str


@dataclass(slots=True)
class ContinuousModeRecoveryPanelExportResult:
    """Materialized copy of the packaged continuous-mode recovery panel."""

    output_root: Path
    readme_path: Path
    default_tree_path: Path
    reference_tree_root: Path
    simulation_cases_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class ContinuousModeRecoveryPanelWorkflowReport:
    """One recovery workflow run over the packaged continuous-mode panel."""

    dataset: ContinuousModeRecoveryPanelDataset
    recovery_report: ContinuousModeRecoveryReport


@dataclass(slots=True)
class ContinuousModeRecoveryPanelWorkflowBundle:
    """Written reviewer-facing outputs for the packaged continuous-mode panel."""

    output_root: Path
    selection_review_case_count: int
    selection_match_count: int
    geiger_selection_match_count: int
    parameter_pass_count: int
    parameter_row_count: int
    parameter_comparison_row_count: int
    parameter_closer_to_truth_count_bijux: int
    parameter_closer_to_truth_count_geiger: int
    expected_warning_case_count: int
    expected_warning_present_count: int
    workflow_summary_path: Path
    recovery_summary_path: Path
    parameter_recovery_path: Path
    parameter_comparison_path: Path
    model_choice_path: Path
    execution_review_path: Path
    warning_review_path: Path
    geiger_reference_path: Path
    simulated_traits_root: Path


@dataclass(slots=True)
class ContinuousModeRecoveryPanelDemoResult:
    """Dataset export plus recovery workflow outputs for the public demo."""

    output_root: Path
    dataset: ContinuousModeRecoveryPanelDataset
    dataset_export: ContinuousModeRecoveryPanelExportResult
    workflow_bundle: ContinuousModeRecoveryPanelWorkflowBundle
    overview_path: Path
