from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.discrete_mode_recovery import (
    DiscreteModeRecoveryReport,
)


@dataclass(slots=True)
class DiscreteModeRecoveryPanelDataset:
    """Packaged deterministic recovery panel with stable and review discrete cases."""

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
class DiscreteModeRecoveryPanelExportResult:
    """Materialized copy of the packaged discrete-mode recovery panel."""

    output_root: Path
    readme_path: Path
    default_tree_path: Path
    reference_tree_root: Path
    simulation_cases_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class DiscreteModeRecoveryPanelWorkflowReport:
    """One recovery workflow run over the packaged discrete-mode panel."""

    dataset: DiscreteModeRecoveryPanelDataset
    recovery_report: DiscreteModeRecoveryReport


@dataclass(slots=True)
class DiscreteModeRecoveryPanelWorkflowBundle:
    """Written reviewer-facing outputs for the packaged discrete-mode panel."""

    output_root: Path
    selection_review_case_count: int
    selection_match_count: int
    geiger_selection_match_count: int
    parameter_pass_count: int
    governed_parameter_row_count: int
    parameter_row_count: int
    parameter_comparison_row_count: int
    parameter_closer_to_truth_count_bijux: int
    parameter_closer_to_truth_count_geiger: int
    rate_pass_count: int
    governed_rate_row_count: int
    rate_row_count: int
    governed_rate_comparison_row_count: int
    rate_comparison_row_count: int
    rate_closer_to_truth_count_bijux: int
    rate_closer_to_truth_count_geiger: int
    expected_warning_case_count: int
    expected_warning_present_count: int
    workflow_summary_path: Path
    recovery_summary_path: Path
    parameter_recovery_path: Path
    parameter_comparison_path: Path
    rate_recovery_path: Path
    rate_comparison_path: Path
    model_choice_path: Path
    execution_review_path: Path
    warning_review_path: Path
    geiger_reference_path: Path
    simulated_traits_root: Path


@dataclass(slots=True)
class DiscreteModeRecoveryPanelDemoResult:
    """Dataset export plus recovery workflow outputs for the public demo."""

    output_root: Path
    dataset: DiscreteModeRecoveryPanelDataset
    dataset_export: DiscreteModeRecoveryPanelExportResult
    workflow_bundle: DiscreteModeRecoveryPanelWorkflowBundle
    overview_path: Path
