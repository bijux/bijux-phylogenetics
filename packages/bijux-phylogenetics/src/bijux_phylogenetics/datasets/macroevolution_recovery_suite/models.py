from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    ContinuousModeRecoveryPanelDataset,
    ContinuousModeRecoveryPanelExportResult,
    ContinuousModeRecoveryPanelWorkflowReport,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    DiscreteModeRecoveryPanelDataset,
    DiscreteModeRecoveryPanelExportResult,
    DiscreteModeRecoveryPanelWorkflowReport,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    KnownAnswerReferenceDataset,
    KnownAnswerReferenceExportResult,
    KnownAnswerReferenceWorkflowReport,
)


@dataclass(slots=True)
class MacroevolutionRecoverySuiteDataset:
    """Packaged suite that unifies governed macroevolution recovery surfaces."""

    dataset_id: str
    label: str
    dataset_root: Path
    reference_output_root: Path
    continuous_panel: ContinuousModeRecoveryPanelDataset
    discrete_panel: DiscreteModeRecoveryPanelDataset
    known_answer_panel: KnownAnswerReferenceDataset
    component_count: int
    geiger_component_count: int
    max_taxon_count: int
    total_recovery_case_count: int
    geiger_recovery_case_count: int
    truth_threshold_row_count: int
    source_summary: str


@dataclass(slots=True)
class MacroevolutionRecoverySuiteExportResult:
    """Materialized copy of the packaged macroevolution recovery suite."""

    output_root: Path
    readme_path: Path
    component_root: Path
    continuous_panel_export: ContinuousModeRecoveryPanelExportResult
    discrete_panel_export: DiscreteModeRecoveryPanelExportResult
    known_answer_panel_export: KnownAnswerReferenceExportResult
    expected_output_root: Path


@dataclass(slots=True)
class MacroevolutionRecoverySuiteComponentSummary:
    """One governed component surface folded into the suite-level summary."""

    dataset_id: str
    label: str
    expected_output_root: Path
    case_count: int
    taxon_count: int
    selection_review_case_count: int
    selection_match_count: int
    geiger_selection_match_count: int
    governed_value_pass_count: int
    governed_value_row_count: int
    governed_comparison_row_count: int
    expected_warning_case_count: int
    expected_warning_present_count: int
    truth_threshold_pass_count: int
    truth_threshold_row_count: int


@dataclass(slots=True)
class MacroevolutionRecoverySuiteWorkflowReport:
    """One governed aggregate over the package-owned recovery component surfaces."""

    dataset: MacroevolutionRecoverySuiteDataset
    continuous_component: MacroevolutionRecoverySuiteComponentSummary
    discrete_component: MacroevolutionRecoverySuiteComponentSummary
    known_answer_component: MacroevolutionRecoverySuiteComponentSummary
    continuous_panel_workflow: ContinuousModeRecoveryPanelWorkflowReport
    discrete_panel_workflow: DiscreteModeRecoveryPanelWorkflowReport
    known_answer_panel_workflow: KnownAnswerReferenceWorkflowReport
    sim_char_case_count: int
    sim_char_all_passed: bool


@dataclass(slots=True)
class MacroevolutionRecoverySuiteWorkflowBundle:
    """Written suite-level outputs plus bundled component recovery evidence."""

    output_root: Path
    component_root: Path
    continuous_component_root: Path
    discrete_component_root: Path
    known_answer_component_root: Path
    component_count: int
    geiger_component_count: int
    total_recovery_case_count: int
    geiger_recovery_case_count: int
    max_taxon_count: int
    selection_review_case_count: int
    selection_match_count: int
    geiger_selection_match_count: int
    governed_value_pass_count: int
    governed_value_row_count: int
    governed_comparison_row_count: int
    expected_warning_case_count: int
    expected_warning_present_count: int
    truth_threshold_pass_count: int
    truth_threshold_row_count: int
    sim_char_case_count: int
    sim_char_all_passed: bool
    requirement_pass_count: int
    requirement_row_count: int
    workflow_summary_path: Path
    component_summary_path: Path
    requirement_summary_path: Path
    sim_char_summary_path: Path


@dataclass(slots=True)
class MacroevolutionRecoverySuiteDemoResult:
    """Dataset export plus suite workflow outputs for the public demo."""

    output_root: Path
    dataset: MacroevolutionRecoverySuiteDataset
    dataset_export: MacroevolutionRecoverySuiteExportResult
    workflow_bundle: MacroevolutionRecoverySuiteWorkflowBundle
    overview_path: Path
