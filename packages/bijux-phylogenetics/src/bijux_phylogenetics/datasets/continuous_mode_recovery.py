from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.comparative.continuous_mode_recovery import (
    ContinuousModeRecoveryReport,
    ContinuousModeRecoveryScenario,
    run_continuous_mode_recovery,
    write_continuous_mode_recovery_execution_table,
    write_continuous_mode_recovery_model_choice_table,
    write_continuous_mode_recovery_parameter_comparison_table,
    write_continuous_mode_recovery_parameter_table,
    write_continuous_mode_recovery_summary_table,
    write_continuous_mode_recovery_warning_table,
    write_geiger_fitcontinuous_recovery_reference_payload_table,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import write_continuous_trait_table

_DATASET_ID = "continuous_mode_recovery_panel"
_DATASET_LABEL = "Continuous trait-model recovery panel"
_DEFAULT_TREE_FILE = "trees/reference-tree-twelve-taxa.nwk"


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


def load_continuous_mode_recovery_panel_dataset() -> ContinuousModeRecoveryPanelDataset:
    """Expose the packaged continuous-mode recovery panel as a first-class surface."""
    dataset_root = _resource_root()
    default_tree_path = dataset_root / _DEFAULT_TREE_FILE
    reference_tree_paths = sorted((dataset_root / "trees").glob("*.nwk"))
    simulation_cases_path = dataset_root / "simulation-cases.tsv"
    taxon_count = max(load_tree(path).tip_count for path in reference_tree_paths)
    case_count = len(_load_scenarios(simulation_cases_path, dataset_root))
    return ContinuousModeRecoveryPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        default_tree_path=default_tree_path,
        reference_tree_paths=reference_tree_paths,
        simulation_cases_path=simulation_cases_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=taxon_count,
        tree_count=len(reference_tree_paths),
        case_count=case_count,
        source_summary=(
            "Deterministic continuous-trait recovery panel with one governed "
            "twelve-taxon review tree for BM, OU, EB, and weak-OU identifiability, "
            "plus one governed twenty-four-taxon ultrametric review tree for "
            "Pagel-lambda, Pagel-kappa, and Pagel-delta transformed-branch "
            "recovery comparisons against stored local geiger references."
        ),
    )


def export_continuous_mode_recovery_panel_dataset(
    destination: Path,
) -> ContinuousModeRecoveryPanelExportResult:
    """Copy the packaged continuous-mode recovery panel and reference outputs."""
    dataset = load_continuous_mode_recovery_panel_dataset()
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
    return ContinuousModeRecoveryPanelExportResult(
        output_root=destination,
        readme_path=Path(readme_path),
        default_tree_path=reference_tree_root / Path(_DEFAULT_TREE_FILE).name,
        reference_tree_root=reference_tree_root,
        simulation_cases_path=Path(simulation_cases_path),
        expected_output_root=expected_output_root,
    )


def run_continuous_mode_recovery_panel_workflow() -> (
    ContinuousModeRecoveryPanelWorkflowReport
):
    """Run the governed recovery workflow over the packaged continuous-mode panel."""
    dataset = load_continuous_mode_recovery_panel_dataset()
    scenarios = _load_scenarios(dataset.simulation_cases_path, dataset.dataset_root)
    recovery_report = run_continuous_mode_recovery(
        dataset.default_tree_path,
        scenarios,
    )
    return ContinuousModeRecoveryPanelWorkflowReport(
        dataset=dataset,
        recovery_report=recovery_report,
    )


def write_continuous_mode_recovery_panel_workflow_bundle(
    output_root: Path,
    report: ContinuousModeRecoveryPanelWorkflowReport,
) -> ContinuousModeRecoveryPanelWorkflowBundle:
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
    parameter_comparison_rows = [
        row
        for case in report.recovery_report.case_reports
        for row in case.parameter_comparison_rows
    ]
    parameter_pass_count = sum(1 for row in parameter_rows if row.within_tolerance)
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
        parameter_row_count=len(parameter_rows),
        parameter_comparison_row_count=len(parameter_comparison_rows),
        parameter_closer_to_truth_count_bijux=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "bijux"
        ),
        parameter_closer_to_truth_count_geiger=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "geiger"
        ),
        expected_warning_case_count=expected_warning_case_count,
        expected_warning_present_count=expected_warning_present_count,
    )
    recovery_summary_path = write_continuous_mode_recovery_summary_table(
        output_root / "recovery-summary.tsv",
        report.recovery_report,
    )
    parameter_recovery_path = write_continuous_mode_recovery_parameter_table(
        output_root / "parameter-recovery.tsv",
        report.recovery_report,
    )
    parameter_comparison_path = write_continuous_mode_recovery_parameter_comparison_table(
        output_root / "parameter-comparison.tsv",
        report.recovery_report,
    )
    model_choice_path = write_continuous_mode_recovery_model_choice_table(
        output_root / "model-choice.tsv",
        report.recovery_report,
    )
    execution_review_path = write_continuous_mode_recovery_execution_table(
        output_root / "execution-review.tsv",
        report.recovery_report,
    )
    warning_review_path = write_continuous_mode_recovery_warning_table(
        output_root / "warning-review.tsv",
        report.recovery_report,
    )
    geiger_reference_path = write_geiger_fitcontinuous_recovery_reference_payload_table(
        output_root / "geiger-reference.tsv",
        report.recovery_report,
    )
    simulated_traits_root = output_root / "simulated-traits"
    simulated_traits_root.mkdir(parents=True, exist_ok=True)
    for case in report.recovery_report.case_reports:
        write_continuous_trait_table(
            simulated_traits_root / f"{case.scenario.case_id}.tsv",
            case.simulation,
        )
    return ContinuousModeRecoveryPanelWorkflowBundle(
        output_root=output_root,
        selection_review_case_count=selection_review_case_count,
        selection_match_count=selection_match_count,
        geiger_selection_match_count=geiger_selection_match_count,
        parameter_pass_count=parameter_pass_count,
        parameter_row_count=len(parameter_rows),
        parameter_comparison_row_count=len(parameter_comparison_rows),
        parameter_closer_to_truth_count_bijux=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "bijux"
        ),
        parameter_closer_to_truth_count_geiger=sum(
            1 for row in parameter_comparison_rows if row.closer_engine == "geiger"
        ),
        expected_warning_case_count=expected_warning_case_count,
        expected_warning_present_count=expected_warning_present_count,
        workflow_summary_path=workflow_summary_path,
        recovery_summary_path=recovery_summary_path,
        parameter_recovery_path=parameter_recovery_path,
        parameter_comparison_path=parameter_comparison_path,
        model_choice_path=model_choice_path,
        execution_review_path=execution_review_path,
        warning_review_path=warning_review_path,
        geiger_reference_path=geiger_reference_path,
        simulated_traits_root=simulated_traits_root,
    )


def run_continuous_mode_recovery_panel_demo(
    output_root: Path,
) -> ContinuousModeRecoveryPanelDemoResult:
    """Materialize the packaged panel and rerun the recovery outputs."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    workflow_report = run_continuous_mode_recovery_panel_workflow()
    dataset_export = export_continuous_mode_recovery_panel_dataset(
        output_root / "dataset"
    )
    workflow_bundle = write_continuous_mode_recovery_panel_workflow_bundle(
        output_root / "workflow",
        workflow_report,
    )
    overview_path = _write_overview(
        output_root / "overview.md",
        workflow_report,
        workflow_bundle,
    )
    return ContinuousModeRecoveryPanelDemoResult(
        output_root=output_root,
        dataset=workflow_report.dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _write_workflow_summary_table(
    path: Path,
    *,
    report: ContinuousModeRecoveryPanelWorkflowReport,
    selection_review_case_count: int,
    selection_match_count: int,
    geiger_selection_match_count: int,
    parameter_pass_count: int,
    parameter_row_count: int,
    parameter_comparison_row_count: int,
    parameter_closer_to_truth_count_bijux: int,
    parameter_closer_to_truth_count_geiger: int,
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
                "parameter_row_count",
                "parameter_comparison_row_count",
                "parameter_closer_to_truth_count_bijux",
                "parameter_closer_to_truth_count_geiger",
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
                str(parameter_row_count),
                str(parameter_comparison_row_count),
                str(parameter_closer_to_truth_count_bijux),
                str(parameter_closer_to_truth_count_geiger),
                str(expected_warning_case_count),
                str(expected_warning_present_count),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_overview(
    path: Path,
    report: ContinuousModeRecoveryPanelWorkflowReport,
    bundle: ContinuousModeRecoveryPanelWorkflowBundle,
) -> Path:
    lines = [
        "# Continuous Trait-Model Recovery Demo",
        "",
        f"- dataset id: `{report.dataset.dataset_id}`",
        f"- governed trees: `{report.dataset.tree_count}`",
        f"- largest taxon count: `{report.dataset.taxon_count}`",
        f"- recovery cases: `{report.dataset.case_count}`",
        f"- selection review cases: `{bundle.selection_review_case_count}`",
        f"- Bijux model-selection matches expectation: `{bundle.selection_match_count}`",
        f"- geiger model-selection matches expectation: `{bundle.geiger_selection_match_count}`",
        f"- parameter recoveries within tolerance: `{bundle.parameter_pass_count}/{bundle.parameter_row_count}`",
        f"- paired parameter comparisons: `{bundle.parameter_comparison_row_count}`",
        f"- parameters closer to truth in Bijux: `{bundle.parameter_closer_to_truth_count_bijux}`",
        f"- parameters closer to truth in geiger: `{bundle.parameter_closer_to_truth_count_geiger}`",
        f"- expected warning cases satisfied: `{bundle.expected_warning_present_count}/{bundle.expected_warning_case_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{bundle.workflow_summary_path.name}`",
        f"- recovery summary: `{bundle.recovery_summary_path.name}`",
        f"- parameter recovery ledger: `{bundle.parameter_recovery_path.name}`",
        f"- parameter comparison ledger: `{bundle.parameter_comparison_path.name}`",
        f"- model-choice ledger: `{bundle.model_choice_path.name}`",
        f"- execution review ledger: `{bundle.execution_review_path.name}`",
        f"- warning ledger: `{bundle.warning_review_path.name}`",
        f"- stored geiger reference ledger: `{bundle.geiger_reference_path.name}`",
        f"- simulated traits directory: `{bundle.simulated_traits_root.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _load_scenarios(
    path: Path,
    dataset_root: Path,
) -> list[ContinuousModeRecoveryScenario]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            ContinuousModeRecoveryScenario(
                case_id=row["case_id"],
                label=row["label"],
                generating_model=row["generating_model"],
                expected_selected_model=(
                    None if not row["expected_selected_model"] else row["expected_selected_model"]
                ),
                root_state=float(row["root_state"]),
                sigma=float(row["sigma"]),
                seed=int(row["seed"]),
                tree_path=dataset_root / row["tree_file"],
                alpha=_optional_float(row["alpha"]),
                theta=_optional_float(row["theta"]),
                rate_change=_optional_float(row["rate_change"]),
                lambda_value=_optional_float(row["lambda_value"]),
                kappa=_optional_float(row["kappa"]),
                delta=_optional_float(row["delta"]),
                candidate_modes=tuple(
                    item for item in row["candidate_modes"].split(",") if item
                ),
                lambda_bounds=(float(row["lambda_lower"]), float(row["lambda_upper"])),
                kappa_bounds=(float(row["kappa_lower"]), float(row["kappa_upper"])),
                delta_bounds=(float(row["delta_lower"]), float(row["delta_upper"])),
                ou_bounds=(float(row["ou_lower"]), float(row["ou_upper"])),
                early_burst_bounds=(
                    float(row["early_burst_lower"]),
                    float(row["early_burst_upper"]),
                ),
                parameter_tolerances=_build_parameter_tolerances(row),
                expected_warning_kinds=_split_items(row["expected_warning_kinds"]),
                notes=row["notes"],
            )
            for row in reader
        ]


def _build_parameter_tolerances(row: dict[str, str]) -> dict[str, float]:
    tolerances: dict[str, float] = {}
    for field_name, parameter_name in (
        ("sigma_squared_tolerance", "sigma_squared"),
        ("alpha_tolerance", "alpha"),
        ("theta_tolerance", "theta"),
        ("rate_change_tolerance", "rate_change"),
        ("lambda_tolerance", "lambda"),
        ("kappa_tolerance", "kappa"),
        ("delta_tolerance", "delta"),
    ):
        if row[field_name]:
            tolerances[parameter_name] = float(row[field_name])
    return tolerances


def _optional_float(value: str) -> float | None:
    return None if not value else float(value)


def _split_items(value: str) -> list[str]:
    if not value:
        return []
    return [item for item in value.split(",") if item]


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "simulation"
        / _DATASET_ID
    )
