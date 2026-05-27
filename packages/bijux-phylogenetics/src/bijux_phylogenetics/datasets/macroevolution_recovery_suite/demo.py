from __future__ import annotations

from pathlib import Path
import shutil

from .bundle import write_macroevolution_recovery_suite_workflow_bundle
from .export import export_macroevolution_recovery_suite_dataset
from .models import (
    MacroevolutionRecoverySuiteDemoResult,
    MacroevolutionRecoverySuiteWorkflowBundle,
    MacroevolutionRecoverySuiteWorkflowReport,
)
from .workflow import run_macroevolution_recovery_suite_workflow


def run_macroevolution_recovery_suite_demo(
    output_root: Path,
) -> MacroevolutionRecoverySuiteDemoResult:
    """Materialize the packaged suite and rerun the governed recovery outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    report = run_macroevolution_recovery_suite_workflow()
    dataset_export = export_macroevolution_recovery_suite_dataset(
        output_root / "dataset"
    )
    workflow_bundle = write_macroevolution_recovery_suite_workflow_bundle(
        output_root / "workflow",
        report,
    )
    overview_path = _write_overview(
        output_root / "overview.md",
        report,
        workflow_bundle,
    )
    return MacroevolutionRecoverySuiteDemoResult(
        output_root=output_root,
        dataset=report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_overview(
    path: Path,
    report: MacroevolutionRecoverySuiteWorkflowReport,
    bundle: MacroevolutionRecoverySuiteWorkflowBundle,
) -> Path:
    lines = [
        "# Macroevolution Recovery Suite Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- component panels: `{bundle.component_count}`",
        f"- geiger-backed recovery panels: `{bundle.geiger_component_count}`",
        f"- governed recovery evaluations: `{bundle.total_recovery_case_count}`",
        f"- geiger-backed recovery evaluations: `{bundle.geiger_recovery_case_count}`",
        f"- maximum governed taxon count: `{bundle.max_taxon_count}`",
        f"- geiger selection matches: `{bundle.geiger_selection_match_count}`",
        f"- governed recovery values within tolerance: `{bundle.governed_value_pass_count}/{bundle.governed_value_row_count}`",
        f"- known-answer truth thresholds passed: `{bundle.truth_threshold_pass_count}/{bundle.truth_threshold_row_count}`",
        f"- sim.char parity passed: `{str(bundle.sim_char_all_passed).lower()}` over `{bundle.sim_char_case_count}` cases",
        f"- suite requirements passed: `{bundle.requirement_pass_count}/{bundle.requirement_row_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- component summary: `{bundle.component_summary_path.name}`",
        f"- requirement summary: `{bundle.requirement_summary_path.name}`",
        f"- sim.char summary: `{bundle.sim_char_summary_path.name}`",
        f"- component bundle root: `{bundle.component_root.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
