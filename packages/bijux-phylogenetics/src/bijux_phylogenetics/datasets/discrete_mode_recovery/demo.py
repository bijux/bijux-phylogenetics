from __future__ import annotations

from pathlib import Path
import shutil

from .bundle import write_discrete_mode_recovery_panel_workflow_bundle
from .export import export_discrete_mode_recovery_panel_dataset
from .models import (
    DiscreteModeRecoveryPanelDemoResult,
    DiscreteModeRecoveryPanelWorkflowBundle,
    DiscreteModeRecoveryPanelWorkflowReport,
)
from .workflow import run_discrete_mode_recovery_panel_workflow


def run_discrete_mode_recovery_panel_demo(
    output_root: Path,
) -> DiscreteModeRecoveryPanelDemoResult:
    """Materialize the packaged panel and rerun the recovery outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    workflow_report = run_discrete_mode_recovery_panel_workflow()
    dataset_export = export_discrete_mode_recovery_panel_dataset(
        output_root / "dataset"
    )
    workflow_bundle = write_discrete_mode_recovery_panel_workflow_bundle(
        output_root / "workflow",
        workflow_report,
    )
    overview_path = _write_overview(
        output_root / "overview.md",
        workflow_report,
        workflow_bundle,
    )
    return DiscreteModeRecoveryPanelDemoResult(
        output_root=output_root,
        dataset=workflow_report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_overview(
    path: Path,
    report: DiscreteModeRecoveryPanelWorkflowReport,
    bundle: DiscreteModeRecoveryPanelWorkflowBundle,
) -> Path:
    lines = [
        "# Discrete Trait-Model Recovery Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- governed trees: `{report.dataset.tree_count}`",
        f"- largest taxon count: `{report.dataset.taxon_count}`",
        f"- recovery cases: `{report.dataset.case_count}`",
        f"- selection review cases: `{bundle.selection_review_case_count}`",
        f"- Bijux model-selection matches expectation: `{bundle.selection_match_count}`",
        f"- geiger model-selection matches expectation: `{bundle.geiger_selection_match_count}`",
        f"- transform parameters within tolerance: `{bundle.parameter_pass_count}/{bundle.governed_parameter_row_count}`",
        f"- all transform-parameter recovery rows: `{bundle.parameter_row_count}`",
        f"- paired transform-parameter comparisons: `{bundle.parameter_comparison_row_count}`",
        f"- transform parameters closer to truth in Bijux: `{bundle.parameter_closer_to_truth_count_bijux}`",
        f"- transform parameters closer to truth in geiger: `{bundle.parameter_closer_to_truth_count_geiger}`",
        f"- rate recoveries within tolerance: `{bundle.rate_pass_count}/{bundle.governed_rate_row_count}`",
        f"- all rate recovery rows: `{bundle.rate_row_count}`",
        f"- governed paired rate comparisons: `{bundle.governed_rate_comparison_row_count}`",
        f"- all paired rate comparisons: `{bundle.rate_comparison_row_count}`",
        f"- transition rates closer to truth in Bijux: `{bundle.rate_closer_to_truth_count_bijux}`",
        f"- transition rates closer to truth in geiger: `{bundle.rate_closer_to_truth_count_geiger}`",
        f"- expected warning cases satisfied: `{bundle.expected_warning_present_count}/{bundle.expected_warning_case_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- recovery summary: `{bundle.recovery_summary_path.name}`",
        f"- transform-parameter recovery ledger: `{bundle.parameter_recovery_path.name}`",
        f"- transform-parameter comparison ledger: `{bundle.parameter_comparison_path.name}`",
        f"- rate recovery ledger: `{bundle.rate_recovery_path.name}`",
        f"- rate comparison ledger: `{bundle.rate_comparison_path.name}`",
        f"- model-choice ledger: `{bundle.model_choice_path.name}`",
        f"- execution review ledger: `{bundle.execution_review_path.name}`",
        f"- warning ledger: `{bundle.warning_review_path.name}`",
        f"- stored geiger reference ledger: `{bundle.geiger_reference_path.name}`",
        f"- simulated traits directory: `{bundle.simulated_traits_root.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
