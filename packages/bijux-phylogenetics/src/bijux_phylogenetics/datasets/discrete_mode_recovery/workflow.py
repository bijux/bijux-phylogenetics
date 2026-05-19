from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.comparative.discrete_mode_recovery import (
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
from bijux_phylogenetics.simulation import write_discrete_trait_table

from .models import (
    DiscreteModeRecoveryPanelDemoResult,
    DiscreteModeRecoveryPanelExportResult,
    DiscreteModeRecoveryPanelWorkflowBundle,
    DiscreteModeRecoveryPanelWorkflowReport,
)
from .panel import DEFAULT_TREE_FILE, load_discrete_mode_recovery_panel_dataset
from .scenarios import load_discrete_mode_recovery_panel_scenarios


def export_discrete_mode_recovery_panel_dataset(
    destination: Path,
) -> DiscreteModeRecoveryPanelExportResult:
    """Copy the packaged discrete-mode recovery panel and reference outputs."""
    dataset = load_discrete_mode_recovery_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = shutil.copy2(
        dataset.dataset_root / "README.md",
        destination / "README.md",
    )
    reference_tree_root = destination / "trees"
    shutil.copytree(dataset.dataset_root / "trees", reference_tree_root)
    simulation_cases_path = shutil.copy2(
        dataset.simulation_cases_path,
        destination / "simulation-cases.tsv",
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return DiscreteModeRecoveryPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        default_tree_path=reference_tree_root / Path(DEFAULT_TREE_FILE).name,
        reference_tree_root=reference_tree_root,
        simulation_cases_path=Path(simulation_cases_path),
        expected_output_root=expected_output_root,
    )


def run_discrete_mode_recovery_panel_workflow() -> DiscreteModeRecoveryPanelWorkflowReport:
    """Run the governed recovery workflow over the packaged discrete-mode panel."""
    dataset = load_discrete_mode_recovery_panel_dataset()
    recovery_report = run_discrete_mode_recovery(
        dataset.default_tree_path,
        load_discrete_mode_recovery_panel_scenarios(dataset),
    )
    return DiscreteModeRecoveryPanelWorkflowReport(
        dataset=dataset,
        recovery_report=recovery_report,
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
    governed_parameter_rows = [row for row in parameter_rows if row.tolerance is not None]
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
        parameter_closer_to_truth_count_bijux=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "bijux"
        ),
        parameter_closer_to_truth_count_geiger=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "geiger"
        ),
        rate_pass_count=rate_pass_count,
        governed_rate_row_count=len(governed_rate_rows),
        rate_row_count=len(rate_rows),
        governed_rate_comparison_row_count=len(governed_rate_comparison_rows),
        rate_comparison_row_count=len(rate_comparison_rows),
        rate_closer_to_truth_count_bijux=sum(
            1 for row in rate_comparison_rows if row.closer_engine == "bijux"
        ),
        rate_closer_to_truth_count_geiger=sum(
            1 for row in rate_comparison_rows if row.closer_engine == "geiger"
        ),
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
        parameter_closer_to_truth_count_bijux=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "bijux"
        ),
        parameter_closer_to_truth_count_geiger=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "geiger"
        ),
        rate_pass_count=rate_pass_count,
        governed_rate_row_count=len(governed_rate_rows),
        rate_row_count=len(rate_rows),
        governed_rate_comparison_row_count=len(governed_rate_comparison_rows),
        rate_comparison_row_count=len(rate_comparison_rows),
        rate_closer_to_truth_count_bijux=sum(
            1 for row in rate_comparison_rows if row.closer_engine == "bijux"
        ),
        rate_closer_to_truth_count_geiger=sum(
            1 for row in rate_comparison_rows if row.closer_engine == "geiger"
        ),
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


def run_discrete_mode_recovery_panel_demo(
    output_root: Path,
) -> DiscreteModeRecoveryPanelDemoResult:
    """Materialize the packaged panel and rerun the recovery outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    workflow_report = run_discrete_mode_recovery_panel_workflow()
    dataset_export = export_discrete_mode_recovery_panel_dataset(output_root / "dataset")
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
