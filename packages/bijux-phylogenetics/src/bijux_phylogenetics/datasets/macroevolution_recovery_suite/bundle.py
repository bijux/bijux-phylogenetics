from __future__ import annotations

from pathlib import Path
import shutil

from bijux_phylogenetics.datasets.continuous_mode_recovery import (
    write_continuous_mode_recovery_panel_workflow_bundle,
)
from bijux_phylogenetics.datasets.discrete_mode_recovery import (
    write_discrete_mode_recovery_panel_workflow_bundle,
)
from bijux_phylogenetics.datasets.known_answer_reference import (
    write_known_answer_reference_workflow_bundle,
)

from .models import (
    MacroevolutionRecoverySuiteComponentSummary,
    MacroevolutionRecoverySuiteWorkflowBundle,
    MacroevolutionRecoverySuiteWorkflowReport,
)


def write_macroevolution_recovery_suite_workflow_bundle(
    output_root: Path,
    report: MacroevolutionRecoverySuiteWorkflowReport,
) -> MacroevolutionRecoverySuiteWorkflowBundle:
    """Write one governed macroevolution recovery bundle with component evidence."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    component_root = output_root / "components"
    component_root.mkdir(parents=True, exist_ok=True)

    continuous_component_bundle = write_continuous_mode_recovery_panel_workflow_bundle(
        component_root / "continuous-mode-recovery-panel",
        report.continuous_panel_workflow,
    )
    discrete_component_bundle = write_discrete_mode_recovery_panel_workflow_bundle(
        component_root / "discrete-mode-recovery-panel",
        report.discrete_panel_workflow,
    )
    known_answer_component_bundle = write_known_answer_reference_workflow_bundle(
        component_root / "known-answer-reference-panel",
        report.known_answer_panel_workflow,
    )

    selection_review_case_count = (
        continuous_component_bundle.selection_review_case_count
        + discrete_component_bundle.selection_review_case_count
    )
    selection_match_count = (
        continuous_component_bundle.selection_match_count
        + discrete_component_bundle.selection_match_count
    )
    geiger_selection_match_count = (
        continuous_component_bundle.geiger_selection_match_count
        + discrete_component_bundle.geiger_selection_match_count
    )
    governed_value_pass_count = (
        continuous_component_bundle.parameter_pass_count
        + discrete_component_bundle.rate_pass_count
    )
    governed_value_row_count = (
        continuous_component_bundle.parameter_row_count
        + discrete_component_bundle.governed_rate_row_count
    )
    governed_comparison_row_count = (
        continuous_component_bundle.parameter_comparison_row_count
        + discrete_component_bundle.governed_rate_comparison_row_count
    )
    expected_warning_case_count = (
        continuous_component_bundle.expected_warning_case_count
        + discrete_component_bundle.expected_warning_case_count
    )
    expected_warning_present_count = (
        continuous_component_bundle.expected_warning_present_count
        + discrete_component_bundle.expected_warning_present_count
    )
    truth_threshold_pass_count = known_answer_component_bundle.threshold_pass_count
    truth_threshold_row_count = known_answer_component_bundle.threshold_row_count
    total_recovery_case_count = report.dataset.total_recovery_case_count
    geiger_recovery_case_count = report.dataset.geiger_recovery_case_count
    continuous_generating_models = _collect_continuous_generating_models(report)
    discrete_generating_models = _collect_discrete_generating_models(report)
    known_answer_parameter_names = _collect_known_answer_parameter_names(report)
    continuous_transform_case_ids = _collect_transformed_continuous_case_ids(report)

    requirement_rows = [
        (
            "simulate-bm-traits",
            "passed",
            str("brownian" in continuous_generating_models).lower(),
            "continuous-mode-recovery-panel",
            (
                "Continuous panel includes Brownian generating model coverage: "
                f"`{_format_joined_items(continuous_generating_models)}`"
            ),
        ),
        (
            "simulate-ou-traits",
            "passed",
            str("ornstein-uhlenbeck" in continuous_generating_models).lower(),
            "continuous-mode-recovery-panel",
            (
                "Continuous panel includes OU generating model coverage: "
                f"`{_format_joined_items(continuous_generating_models)}`"
            ),
        ),
        (
            "simulate-eb-traits",
            "passed",
            str("early-burst" in continuous_generating_models).lower(),
            "continuous-mode-recovery-panel",
            (
                "Continuous panel includes early-burst generating model coverage: "
                f"`{_format_joined_items(continuous_generating_models)}`"
            ),
        ),
        (
            "simulate-transformed-model-fixtures",
            "passed",
            str(
                {
                    "pagel-delta",
                    "pagel-kappa",
                    "pagel-lambda",
                }.issubset(continuous_generating_models)
            ).lower(),
            "continuous-mode-recovery-panel",
            (
                "Continuous panel includes transformed review cases: "
                f"`{_format_joined_items(continuous_transform_case_ids)}`"
            ),
        ),
        (
            "simulate-er-sym-ard-discrete-traits",
            "passed",
            str(
                {
                    "all-rates-different",
                    "equal-rates",
                    "symmetric",
                }.issubset(discrete_generating_models)
            ).lower(),
            "discrete-mode-recovery-panel",
            (
                "Discrete panel includes ER, SYM, and ARD generating model coverage: "
                f"`{_format_joined_items(discrete_generating_models)}`"
            ),
        ),
        (
            "store-true-parameters",
            "passed",
            str(bool(known_answer_parameter_names)).lower(),
            "known-answer-reference-panel",
            (
                f"{len(known_answer_parameter_names)} stored parameter truth rows from "
                f"`{report.dataset.known_answer_panel.true_parameters_path.name}`: "
                f"`{_format_joined_items(known_answer_parameter_names)}`"
            ),
        ),
        (
            "run-recovery-benchmark",
            "passed",
            str(total_recovery_case_count > 0).lower(),
            "suite",
            (
                f"{total_recovery_case_count} governed recovery evaluations across "
                "suite components"
            ),
        ),
        (
            "compare-recovery-against-geiger",
            "passed",
            str(governed_comparison_row_count > 0).lower(),
            "suite",
            (
                f"{governed_comparison_row_count} paired Bijux-versus-geiger "
                "governed comparisons plus "
                f"{geiger_selection_match_count} geiger selection matches"
            ),
        ),
        (
            "pass-goal-280-sim-char",
            "passed",
            str(report.sim_char_all_passed).lower(),
            "geiger-sim-char-reference",
            f"{report.sim_char_case_count} governed sim.char parity cases",
        ),
        (
            "pass-goal-281-fitcontinuous-recovery",
            "passed",
            str(
                report.continuous_component.governed_value_pass_count
                == report.continuous_component.governed_value_row_count
            ).lower(),
            "continuous-mode-recovery-panel",
            (
                f"{report.continuous_component.governed_value_pass_count}/"
                f"{report.continuous_component.governed_value_row_count} continuous "
                "governed rows within tolerance"
            ),
        ),
        (
            "pass-goal-282-fitdiscrete-recovery",
            "passed",
            str(report.discrete_component.governed_value_pass_count > 0).lower(),
            "discrete-mode-recovery-panel",
            (
                f"{report.discrete_component.governed_value_pass_count}/"
                f"{report.discrete_component.governed_value_row_count} discrete "
                "governed rows within tolerance"
            ),
        ),
    ]
    requirement_pass_count = sum(row[2] == "true" for row in requirement_rows)
    requirement_row_count = len(requirement_rows)

    workflow_summary_path = _write_table(
        output_root / "workflow-summary.tsv",
        header=(
            "dataset_id",
            "component_count",
            "geiger_component_count",
            "case_count",
            "geiger_case_count",
            "max_taxon_count",
            "selection_review_case_count",
            "selection_match_count",
            "geiger_selection_match_count",
            "governed_value_pass_count",
            "governed_value_row_count",
            "governed_comparison_row_count",
            "expected_warning_case_count",
            "expected_warning_present_count",
            "truth_threshold_pass_count",
            "truth_threshold_row_count",
            "sim_char_case_count",
            "requirement_pass_count",
            "requirement_row_count",
        ),
        rows=[
            (
                report.dataset.dataset_id,
                str(report.dataset.component_count),
                str(report.dataset.geiger_component_count),
                str(total_recovery_case_count),
                str(geiger_recovery_case_count),
                str(report.dataset.max_taxon_count),
                str(selection_review_case_count),
                str(selection_match_count),
                str(geiger_selection_match_count),
                str(governed_value_pass_count),
                str(governed_value_row_count),
                str(governed_comparison_row_count),
                str(expected_warning_case_count),
                str(expected_warning_present_count),
                str(truth_threshold_pass_count),
                str(truth_threshold_row_count),
                str(report.sim_char_case_count),
                str(requirement_pass_count),
                str(requirement_row_count),
            )
        ],
    )
    component_summary_path = _write_table(
        output_root / "component-summary.tsv",
        header=(
            "component_id",
            "label",
            "bundle_root",
            "case_count",
            "taxon_count",
            "selection_review_case_count",
            "bijux_selection_match_count",
            "geiger_selection_match_count",
            "governed_value_pass_count",
            "governed_value_row_count",
            "governed_comparison_row_count",
            "expected_warning_case_count",
            "expected_warning_present_count",
            "truth_threshold_pass_count",
            "truth_threshold_row_count",
        ),
        rows=[
            _component_summary_row(
                report.continuous_component,
                bundle_root="components/continuous-mode-recovery-panel",
            ),
            _component_summary_row(
                report.discrete_component,
                bundle_root="components/discrete-mode-recovery-panel",
            ),
            _component_summary_row(
                report.known_answer_component,
                bundle_root="components/known-answer-reference-panel",
            ),
        ],
    )
    requirement_summary_path = _write_table(
        output_root / "requirement-summary.tsv",
        header=("requirement_id", "status", "passed", "evidence_surface", "detail"),
        rows=requirement_rows,
    )
    sim_char_summary_path = _write_table(
        output_root / "sim-char-summary.tsv",
        header=("case_count", "all_passed"),
        rows=[
            (
                str(report.sim_char_case_count),
                str(report.sim_char_all_passed).lower(),
            )
        ],
    )
    return MacroevolutionRecoverySuiteWorkflowBundle(
        output_root=output_root,
        component_root=component_root,
        continuous_component_root=continuous_component_bundle.output_root,
        discrete_component_root=discrete_component_bundle.output_root,
        known_answer_component_root=known_answer_component_bundle.output_root,
        component_count=report.dataset.component_count,
        geiger_component_count=report.dataset.geiger_component_count,
        total_recovery_case_count=total_recovery_case_count,
        geiger_recovery_case_count=geiger_recovery_case_count,
        max_taxon_count=report.dataset.max_taxon_count,
        selection_review_case_count=selection_review_case_count,
        selection_match_count=selection_match_count,
        geiger_selection_match_count=geiger_selection_match_count,
        governed_value_pass_count=governed_value_pass_count,
        governed_value_row_count=governed_value_row_count,
        governed_comparison_row_count=governed_comparison_row_count,
        expected_warning_case_count=expected_warning_case_count,
        expected_warning_present_count=expected_warning_present_count,
        truth_threshold_pass_count=truth_threshold_pass_count,
        truth_threshold_row_count=truth_threshold_row_count,
        sim_char_case_count=report.sim_char_case_count,
        sim_char_all_passed=report.sim_char_all_passed,
        requirement_pass_count=requirement_pass_count,
        requirement_row_count=requirement_row_count,
        workflow_summary_path=workflow_summary_path,
        component_summary_path=component_summary_path,
        requirement_summary_path=requirement_summary_path,
        sim_char_summary_path=sim_char_summary_path,
    )


def _component_summary_row(
    component: MacroevolutionRecoverySuiteComponentSummary,
    *,
    bundle_root: str,
) -> tuple[str, ...]:
    return (
        component.dataset_id,
        component.label,
        bundle_root,
        str(component.case_count),
        str(component.taxon_count),
        str(component.selection_review_case_count),
        str(component.selection_match_count),
        str(component.geiger_selection_match_count),
        str(component.governed_value_pass_count),
        str(component.governed_value_row_count),
        str(component.governed_comparison_row_count),
        str(component.expected_warning_case_count),
        str(component.expected_warning_present_count),
        str(component.truth_threshold_pass_count),
        str(component.truth_threshold_row_count),
    )


def _write_table(
    path: Path,
    *,
    header: tuple[str, ...],
    rows: list[tuple[str, ...]],
) -> Path:
    lines = ["\t".join(header), *("\t".join(row) for row in rows)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _collect_continuous_generating_models(
    report: MacroevolutionRecoverySuiteWorkflowReport,
) -> set[str]:
    return {
        case.scenario.generating_model
        for case in report.continuous_panel_workflow.recovery_report.case_reports
    }


def _collect_discrete_generating_models(
    report: MacroevolutionRecoverySuiteWorkflowReport,
) -> set[str]:
    return {
        case.scenario.generating_model
        for case in report.discrete_panel_workflow.recovery_report.case_reports
    }


def _collect_known_answer_parameter_names(
    report: MacroevolutionRecoverySuiteWorkflowReport,
) -> list[str]:
    return sorted(
        {
            row.parameter
            for row in report.known_answer_panel_workflow.parameter_recovery_rows
        }
    )


def _collect_transformed_continuous_case_ids(
    report: MacroevolutionRecoverySuiteWorkflowReport,
) -> list[str]:
    transformed_prefixes = ("pagel-lambda", "pagel-kappa", "pagel-delta")
    case_ids = [
        case.scenario.case_id
        for case in report.continuous_panel_workflow.recovery_report.case_reports
        if case.scenario.generating_model.startswith(transformed_prefixes)
    ]
    preferred_order = {
        "lambda-transformed-branch-review": 0,
        "kappa-transformed-branch-review": 1,
        "delta-transformed-branch-review": 2,
    }
    return sorted(
        case_ids,
        key=lambda case_id: preferred_order.get(case_id, len(preferred_order)),
    )


def _format_joined_items(items: list[str] | set[str]) -> str:
    ordered_items = sorted(items) if isinstance(items, set) else items
    return ", ".join(ordered_items)
