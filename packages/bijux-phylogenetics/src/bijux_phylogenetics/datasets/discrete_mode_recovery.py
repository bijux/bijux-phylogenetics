from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil

from bijux_phylogenetics.comparative.discrete_mode_recovery import (
    DiscreteModeRecoveryReport,
    DiscreteModeRecoveryScenario,
    run_discrete_mode_recovery,
    write_discrete_mode_recovery_execution_table,
    write_discrete_mode_recovery_model_choice_table,
    write_discrete_mode_recovery_rate_comparison_table,
    write_discrete_mode_recovery_rate_table,
    write_discrete_mode_recovery_summary_table,
    write_discrete_mode_recovery_warning_table,
    write_geiger_fitdiscrete_recovery_reference_payload_table,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.simulation import DiscreteHistoryRateRow, write_discrete_trait_table

_DATASET_ID = "discrete_mode_recovery_panel"
_DATASET_LABEL = "Discrete trait-model recovery panel"
_DEFAULT_TREE_FILE = "trees/reference-tree-twelve-taxa.nwk"


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


def load_discrete_mode_recovery_panel_dataset() -> DiscreteModeRecoveryPanelDataset:
    """Expose the packaged discrete-mode recovery panel as a first-class surface."""
    dataset_root = _resource_root()
    default_tree_path = dataset_root / _DEFAULT_TREE_FILE
    reference_tree_paths = sorted((dataset_root / "trees").glob("*.nwk"))
    simulation_cases_path = dataset_root / "simulation-cases.tsv"
    taxon_count = max(load_tree(path).tip_count for path in reference_tree_paths)
    case_count = len(_load_scenarios(simulation_cases_path, dataset_root))
    return DiscreteModeRecoveryPanelDataset(
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
            "Deterministic discrete-trait recovery panel with one governed "
            "twelve-taxon overparameterized ARD failure tree plus one governed "
            "twenty-four-taxon rooted review tree for stable ER and SYM "
            "selection cases and one weak-identification ARD review surface "
            "compared against stored local geiger references."
        ),
    )


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
        default_tree_path=reference_tree_root / Path(_DEFAULT_TREE_FILE).name,
        reference_tree_root=reference_tree_root,
        simulation_cases_path=Path(simulation_cases_path),
        expected_output_root=expected_output_root,
    )


def run_discrete_mode_recovery_panel_workflow() -> DiscreteModeRecoveryPanelWorkflowReport:
    """Run the governed recovery workflow over the packaged discrete-mode panel."""
    dataset = load_discrete_mode_recovery_panel_dataset()
    scenarios = _load_scenarios(dataset.simulation_cases_path, dataset.dataset_root)
    recovery_report = run_discrete_mode_recovery(
        dataset.default_tree_path,
        scenarios,
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


def _load_scenarios(
    path: Path,
    dataset_root: Path,
) -> list[DiscreteModeRecoveryScenario]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            DiscreteModeRecoveryScenario(
                case_id=row["case_id"],
                label=row["label"],
                generating_model=row["generating_model"],
                expected_selected_model=(
                    None if not row["expected_selected_model"] else row["expected_selected_model"]
                ),
                states=_split_items(row["states"]),
                rate_rows=_parse_rate_rows(row["rate_rows"]),
                root_state=row["root_state"],
                seed=int(row["seed"]),
                tree_path=dataset_root / row["tree_file"],
                transform=_optional_string(row.get("transform", "")),
                transform_parameter_value=_optional_float(
                    row.get("transform_parameter_value", "")
                ),
                rate_tolerance=_optional_float(row["rate_tolerance"]),
                parameter_tolerances=_parse_parameter_tolerances(
                    row.get("parameter_tolerances", "")
                ),
                lambda_bounds=_parse_bounds(
                    row.get("lambda_bounds", ""),
                    default=(0.0, 1.0),
                ),
                kappa_bounds=_parse_bounds(
                    row.get("kappa_bounds", ""),
                    default=(0.0, 1.0),
                ),
                delta_bounds=_parse_bounds(
                    row.get("delta_bounds", ""),
                    default=(0.006737947, 3.0),
                ),
                early_burst_bounds=_parse_bounds(
                    row.get("early_burst_bounds", ""),
                    default=(-10.0, 10.0),
                ),
                expected_overparameterized=(row["expected_overparameterized"] == "true"),
                expected_warning_kinds=_split_items(row["expected_warning_kinds"]),
                notes=row["notes"],
            )
            for row in reader
        ]


def _parse_rate_rows(value: str) -> list[DiscreteHistoryRateRow]:
    rows: list[DiscreteHistoryRateRow] = []
    for item in (entry for entry in value.split(";") if entry):
        pair, rate_text = item.split("=")
        source_state, target_state = pair.split(">")
        rows.append(
            DiscreteHistoryRateRow(
                source_state=source_state,
                target_state=target_state,
                rate=float(rate_text),
            )
        )
    return rows


def _optional_float(value: str) -> float | None:
    return None if not value else float(value)


def _optional_string(value: str) -> str | None:
    return None if not value else value


def _parse_parameter_tolerances(value: str) -> dict[str, float]:
    tolerances: dict[str, float] = {}
    for item in _split_items(value):
        parameter, tolerance = item.split("=")
        tolerances[parameter] = float(tolerance)
    return tolerances


def _parse_bounds(
    value: str,
    *,
    default: tuple[float, float],
) -> tuple[float, float]:
    if not value:
        return default
    left, right = value.split(":")
    return (float(left), float(right))


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
