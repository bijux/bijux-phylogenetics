from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from bijux_phylogenetics.core.alignment import (
    AlignmentRecord,
    FastaInputValidationReport,
    SequenceCompositionOutlier,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.pruning import drop_tree_taxa
from bijux_phylogenetics.core.traits import (
    TraitLinkageReport,
    TraitValidationReport,
    link_tree_to_traits,
    validate_traits_table,
)
from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.diagnostics.validation import (
    TreeInspectionReport,
    TreeValidationReport,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.io.fasta import (
    detect_composition_outlier_sequences,
    load_fasta_alignment,
    validate_fasta_input,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.newick import write_newick

_DATASET_ID = "catarrhine_data_quality_stress_panel"
_DATASET_LABEL = "Catarrhine data quality stress panel"
_RAW_ALIGNMENT_NAME = "alignment.fasta"
_RAW_TREE_NAME = "tree.nwk"
_RAW_TRAITS_NAME = "traits.csv"
_SEQUENCE_TYPE = "dna"
_WORKFLOW_REQUIRED_TRAITS = ("body_mass_g", "gestation_days")
_TREE_BRANCH_FLOOR = 1e-6


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelDataset:
    """Packaged stress panel with explicit sequence, tree, and trait defects."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    raw_alignment_path: Path
    raw_tree_path: Path
    raw_traits_path: Path
    reference_output_root: Path
    taxon_count: int
    raw_trait_row_count: int
    required_traits: tuple[str, ...]
    sequence_type: str
    source_summary: str


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelExportResult:
    """Materialized copy of the packaged stress panel."""

    output_root: Path
    readme_path: Path
    raw_alignment_path: Path
    raw_tree_path: Path
    raw_traits_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class TraitDuplicateResolution:
    """Deterministic duplicate-taxon choice for one raw trait surface."""

    taxon: str
    occurrence_count: int
    selected_row_number: int
    selected_non_missing_field_count: int
    discarded_row_numbers: list[int]
    selected_reason: str


@dataclass(slots=True)
class TraitMissingObservation:
    """One missing trait observation in the raw traits table."""

    taxon: str
    row_number: int
    trait: str
    required_for_analysis: bool
    action: str


@dataclass(slots=True)
class DataQualityRepairAction:
    """One explicit action taken to move from raw inputs to the cleaned subset."""

    action_kind: str
    affected_taxa: list[str]
    affected_nodes: list[str]
    reason: str
    result: str


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelWorkflowReport:
    """Audit and cleanup report for the packaged stress panel."""

    dataset: CatarrhineDataQualityStressPanelDataset
    raw_alignment_validation: FastaInputValidationReport
    sequence_outliers: list[SequenceCompositionOutlier]
    raw_tree_inspection: TreeInspectionReport
    raw_tree_validation: TreeValidationReport
    trait_duplicates: list[TraitDuplicateResolution]
    missing_traits: list[TraitMissingObservation]
    cleaned_trait_validation: TraitValidationReport
    cleaned_tree_validation: TreeValidationReport
    cleaned_linkage: TraitLinkageReport
    cleaned_alignment_validation: FastaInputValidationReport
    cleaned_alignment_records: list[AlignmentRecord]
    cleaned_tree_path: Path
    cleaned_traits_path: Path
    cleaned_alignment_path: Path
    cleaned_taxa: list[str]
    dropped_taxa: list[str]
    repair_actions: list[DataQualityRepairAction]
    repaired_branch_nodes: list[str]


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelWorkflowBundle:
    """Written workflow outputs for the packaged stress panel."""

    output_root: Path
    raw_taxon_count: int
    cleaned_taxon_count: int
    duplicate_trait_taxon_count: int
    missing_trait_value_count: int
    sequence_outlier_count: int
    tree_zero_length_branch_count: int
    tree_long_branch_outlier_count: int
    dropped_taxon_count: int
    repaired_branch_count: int
    workflow_summary_path: Path
    trait_duplicates_path: Path
    trait_missing_values_path: Path
    sequence_outliers_path: Path
    tree_issues_path: Path
    repair_actions_path: Path
    cleaned_traits_path: Path
    cleaned_alignment_path: Path
    cleaned_tree_path: Path
    cleaned_linkage_path: Path
    cleaned_validation_path: Path


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelDemoResult:
    """Dataset export plus written workflow bundle for the public demo."""

    output_root: Path
    dataset: CatarrhineDataQualityStressPanelDataset
    dataset_export: CatarrhineDataQualityStressPanelExportResult
    workflow_bundle: CatarrhineDataQualityStressPanelWorkflowBundle
    overview_path: Path


def load_catarrhine_data_quality_stress_panel_dataset() -> (
    CatarrhineDataQualityStressPanelDataset
):
    """Expose the packaged catarrhine stress panel as a first-class dataset surface."""
    dataset_root = _resource_root()
    raw_root = dataset_root / "raw"
    raw_alignment_path = raw_root / _RAW_ALIGNMENT_NAME
    raw_tree_path = raw_root / _RAW_TREE_NAME
    raw_traits_path = raw_root / _RAW_TRAITS_NAME
    taxon_count = validate_fasta_input(
        raw_alignment_path, sequence_type=_SEQUENCE_TYPE
    ).summary.sequence_count
    raw_trait_row_count = len(_load_permissive_trait_rows(raw_traits_path))
    return CatarrhineDataQualityStressPanelDataset(
        dataset_id=_DATASET_ID,
        label=_DATASET_LABEL,
        dataset_root=dataset_root,
        readme_path=dataset_root / "README.md",
        raw_alignment_path=raw_alignment_path,
        raw_tree_path=raw_tree_path,
        raw_traits_path=raw_traits_path,
        reference_output_root=dataset_root / "expected",
        taxon_count=taxon_count,
        raw_trait_row_count=raw_trait_row_count,
        required_traits=_WORKFLOW_REQUIRED_TRAITS,
        sequence_type=_SEQUENCE_TYPE,
        source_summary=(
            "Real catarrhine mitochondrial sequence and topology material packaged "
            "with deliberate data-quality defects so duplicate traits, missing "
            "values, sequence outliers, and branch-length pathologies remain "
            "visible and reviewable."
        ),
    )


def export_catarrhine_data_quality_stress_panel_dataset(
    destination: Path,
) -> CatarrhineDataQualityStressPanelExportResult:
    """Copy the packaged raw stress inputs and governed expected outputs."""
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    if destination.exists():
        shutil.rmtree(destination)
    destination.mkdir(parents=True, exist_ok=True)
    readme_path = Path(shutil.copy2(dataset.readme_path, destination / "README.md"))
    raw_root = destination / "raw"
    raw_root.mkdir(parents=True, exist_ok=True)
    raw_alignment_path = Path(
        shutil.copy2(dataset.raw_alignment_path, raw_root / _RAW_ALIGNMENT_NAME)
    )
    raw_tree_path = Path(shutil.copy2(dataset.raw_tree_path, raw_root / _RAW_TREE_NAME))
    raw_traits_path = Path(
        shutil.copy2(dataset.raw_traits_path, raw_root / _RAW_TRAITS_NAME)
    )
    expected_output_root = destination / "expected"
    shutil.copytree(dataset.reference_output_root, expected_output_root)
    return CatarrhineDataQualityStressPanelExportResult(
        output_root=destination,
        readme_path=readme_path,
        raw_alignment_path=raw_alignment_path,
        raw_tree_path=raw_tree_path,
        raw_traits_path=raw_traits_path,
        expected_output_root=expected_output_root,
    )


def run_catarrhine_data_quality_stress_panel_workflow(
    out_dir: Path,
) -> CatarrhineDataQualityStressPanelWorkflowReport:
    """Audit the raw stress panel and write one cleaned comparative subset."""
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    out_dir.mkdir(parents=True, exist_ok=True)
    assembled_root = out_dir / "assembled"
    assembled_root.mkdir(parents=True, exist_ok=True)

    raw_alignment_validation = validate_fasta_input(
        dataset.raw_alignment_path,
        sequence_type=dataset.sequence_type,
    )
    sequence_outliers = detect_composition_outlier_sequences(dataset.raw_alignment_path)
    raw_tree_inspection = inspect_tree_path(dataset.raw_tree_path)
    raw_tree_validation = validate_tree_path(dataset.raw_tree_path)

    raw_trait_rows = _load_permissive_trait_rows(dataset.raw_traits_path)
    trait_duplicates = _resolve_duplicate_traits(raw_trait_rows)
    duplicate_lookup = {row.taxon: row for row in trait_duplicates}
    selected_trait_rows = _selected_trait_rows(raw_trait_rows)
    missing_traits = _detect_missing_traits(
        raw_trait_rows,
        required_traits=set(dataset.required_traits),
        duplicate_lookup=duplicate_lookup,
    )

    sequence_outlier_taxa = sorted(
        {row.identifier for row in sequence_outliers if row.identifier}
    )
    tree_outlier_taxa = sorted(raw_tree_inspection.long_branch_taxa)
    missing_required_trait_taxa = sorted(
        {
            row.taxon
            for row in missing_traits
            if row.required_for_analysis
            and row.action == "drop_taxon_from_cleaned_traits"
        }
    )
    excluded_taxa = sorted(
        set(sequence_outlier_taxa)
        | set(tree_outlier_taxa)
        | set(missing_required_trait_taxa)
    )

    cleaned_trait_rows = [
        dict(row)
        for row in selected_trait_rows
        if row["taxon"] not in excluded_taxa
        and all(row[trait] for trait in dataset.required_traits)
    ]
    cleaned_traits_path = write_taxon_rows(
        assembled_root / "cleaned-traits.csv",
        columns=list(cleaned_trait_rows[0].keys()),
        rows=cleaned_trait_rows,
    )

    cleaned_tree, _ = drop_tree_taxa(dataset.raw_tree_path, excluded_taxa)
    repaired_branch_nodes = _apply_branch_length_floor(
        cleaned_tree.root, floor=_TREE_BRANCH_FLOOR
    )
    cleaned_tree_path = write_newick(assembled_root / "cleaned-tree.nwk", cleaned_tree)

    cleaned_alignment_records = [
        record
        for record in load_fasta_alignment(dataset.raw_alignment_path)
        if record.identifier not in excluded_taxa
    ]
    cleaned_alignment_path = write_fasta_alignment(
        assembled_root / "cleaned-alignment.fasta",
        cleaned_alignment_records,
    )

    cleaned_trait_validation = validate_traits_table(cleaned_traits_path)
    cleaned_tree_validation = validate_tree_path(cleaned_tree_path)
    cleaned_alignment_validation = validate_fasta_input(
        cleaned_alignment_path,
        sequence_type=dataset.sequence_type,
    )
    cleaned_linkage = link_tree_to_traits(
        cleaned_tree_path,
        cleaned_traits_path,
        strict=True,
    )
    cleaned_taxa = sorted(cleaned_linkage.usable_taxa)

    repair_actions = _build_repair_actions(
        trait_duplicates=trait_duplicates,
        missing_required_trait_taxa=missing_required_trait_taxa,
        sequence_outlier_taxa=sequence_outlier_taxa,
        tree_outlier_taxa=tree_outlier_taxa,
        repaired_branch_nodes=repaired_branch_nodes,
    )
    return CatarrhineDataQualityStressPanelWorkflowReport(
        dataset=dataset,
        raw_alignment_validation=raw_alignment_validation,
        sequence_outliers=sequence_outliers,
        raw_tree_inspection=raw_tree_inspection,
        raw_tree_validation=raw_tree_validation,
        trait_duplicates=trait_duplicates,
        missing_traits=missing_traits,
        cleaned_trait_validation=cleaned_trait_validation,
        cleaned_tree_validation=cleaned_tree_validation,
        cleaned_linkage=cleaned_linkage,
        cleaned_alignment_validation=cleaned_alignment_validation,
        cleaned_alignment_records=cleaned_alignment_records,
        cleaned_tree_path=cleaned_tree_path,
        cleaned_traits_path=cleaned_traits_path,
        cleaned_alignment_path=cleaned_alignment_path,
        cleaned_taxa=cleaned_taxa,
        dropped_taxa=excluded_taxa,
        repair_actions=repair_actions,
        repaired_branch_nodes=repaired_branch_nodes,
    )


def write_catarrhine_data_quality_stress_panel_workflow_bundle(
    output_root: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> CatarrhineDataQualityStressPanelWorkflowBundle:
    """Write reviewer-facing ledgers for the stress-panel cleanup workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    trait_duplicates_path = _write_trait_duplicates_table(
        output_root / "trait-duplicates.tsv",
        report.trait_duplicates,
    )
    trait_missing_values_path = _write_trait_missing_values_table(
        output_root / "trait-missing-values.tsv",
        report.missing_traits,
    )
    sequence_outliers_path = _write_sequence_outliers_table(
        output_root / "sequence-outliers.tsv",
        report.sequence_outliers,
        dropped_taxa=set(report.dropped_taxa),
    )
    tree_issues_path = _write_tree_issues_table(
        output_root / "tree-issues.tsv",
        report,
    )
    repair_actions_path = _write_repair_actions_table(
        output_root / "repair-actions.tsv",
        report.repair_actions,
    )
    cleaned_traits_path = _copy_output(
        report.cleaned_traits_path,
        output_root / "cleaned-traits.csv",
    )
    cleaned_alignment_path = _copy_output(
        report.cleaned_alignment_path,
        output_root / "cleaned-alignment.fasta",
    )
    cleaned_tree_path = _copy_output(
        report.cleaned_tree_path,
        output_root / "cleaned-tree.nwk",
    )
    cleaned_linkage_path = _write_cleaned_linkage_table(
        output_root / "cleaned-linkage.tsv",
        report,
    )
    cleaned_validation_path = _write_cleaned_validation_table(
        output_root / "cleaned-validation.tsv",
        report,
    )
    return CatarrhineDataQualityStressPanelWorkflowBundle(
        output_root=output_root,
        raw_taxon_count=report.dataset.taxon_count,
        cleaned_taxon_count=len(report.cleaned_taxa),
        duplicate_trait_taxon_count=len(report.trait_duplicates),
        missing_trait_value_count=len(report.missing_traits),
        sequence_outlier_count=len(report.sequence_outliers),
        tree_zero_length_branch_count=report.raw_tree_validation.zero_length_branches,
        tree_long_branch_outlier_count=len(
            report.raw_tree_inspection.long_branch_outliers
        ),
        dropped_taxon_count=len(report.dropped_taxa),
        repaired_branch_count=len(report.repaired_branch_nodes),
        workflow_summary_path=workflow_summary_path,
        trait_duplicates_path=trait_duplicates_path,
        trait_missing_values_path=trait_missing_values_path,
        sequence_outliers_path=sequence_outliers_path,
        tree_issues_path=tree_issues_path,
        repair_actions_path=repair_actions_path,
        cleaned_traits_path=cleaned_traits_path,
        cleaned_alignment_path=cleaned_alignment_path,
        cleaned_tree_path=cleaned_tree_path,
        cleaned_linkage_path=cleaned_linkage_path,
        cleaned_validation_path=cleaned_validation_path,
    )


def run_catarrhine_data_quality_stress_panel_demo(
    output_root: Path,
) -> CatarrhineDataQualityStressPanelDemoResult:
    """Materialize the packaged stress dataset and rerun the governed cleanup workflow."""
    if output_root.exists():
        shutil.rmtree(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    dataset_export = export_catarrhine_data_quality_stress_panel_dataset(
        output_root / "dataset"
    )
    with TemporaryDirectory(prefix="catarrhine-data-quality-stress-") as temporary_root:
        workflow_report = run_catarrhine_data_quality_stress_panel_workflow(
            Path(temporary_root)
        )
        workflow_bundle = write_catarrhine_data_quality_stress_panel_workflow_bundle(
            output_root / "workflow",
            workflow_report,
        )
    overview_path = _write_overview(
        output_root / "overview.md", dataset, workflow_bundle
    )
    return CatarrhineDataQualityStressPanelDemoResult(
        output_root=output_root,
        dataset=dataset,
        dataset_export=dataset_export,
        workflow_bundle=workflow_bundle,
        overview_path=overview_path,
    )


def _copy_output(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    return Path(shutil.copy2(source, destination))


def _load_permissive_trait_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"trait table has no header row: {path}")
        columns = [column.strip() for column in reader.fieldnames]
        return [
            {column: str(row.get(column, "")).strip() for column in columns}
            for row in reader
        ]


def _selected_trait_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[str, list[tuple[int, dict[str, str]]]] = {}
    for index, row in enumerate(rows, start=2):
        grouped.setdefault(row["taxon"], []).append((index, row))
    selected: list[tuple[int, dict[str, str]]] = []
    for entries in grouped.values():
        selected.append(_select_trait_row(entries))
    return [dict(row) for _, row in sorted(selected, key=lambda item: item[0])]


def _resolve_duplicate_traits(
    rows: list[dict[str, str]],
) -> list[TraitDuplicateResolution]:
    grouped: dict[str, list[tuple[int, dict[str, str]]]] = {}
    for index, row in enumerate(rows, start=2):
        grouped.setdefault(row["taxon"], []).append((index, row))
    resolutions: list[TraitDuplicateResolution] = []
    for taxon, entries in sorted(grouped.items()):
        if len(entries) < 2:
            continue
        selected_row_number, selected_row = _select_trait_row(entries)
        resolutions.append(
            TraitDuplicateResolution(
                taxon=taxon,
                occurrence_count=len(entries),
                selected_row_number=selected_row_number,
                selected_non_missing_field_count=_non_missing_field_count(selected_row),
                discarded_row_numbers=[
                    row_number
                    for row_number, _ in entries
                    if row_number != selected_row_number
                ],
                selected_reason="highest_non_missing_field_count_then_first_row",
            )
        )
    return resolutions


def _select_trait_row(
    entries: list[tuple[int, dict[str, str]]],
) -> tuple[int, dict[str, str]]:
    return max(
        entries,
        key=lambda item: (_non_missing_field_count(item[1]), -item[0]),
    )


def _non_missing_field_count(row: dict[str, str]) -> int:
    return sum(1 for key, value in row.items() if key != "taxon" and value)


def _detect_missing_traits(
    rows: list[dict[str, str]],
    *,
    required_traits: set[str],
    duplicate_lookup: dict[str, TraitDuplicateResolution],
) -> list[TraitMissingObservation]:
    observations: list[TraitMissingObservation] = []
    for row_number, row in enumerate(rows, start=2):
        for trait, value in row.items():
            if trait == "taxon" or value:
                continue
            required_for_analysis = trait in required_traits
            if row["taxon"] in duplicate_lookup and (
                row_number != duplicate_lookup[row["taxon"]].selected_row_number
            ):
                action = "dropped_duplicate_row"
            elif required_for_analysis:
                action = "drop_taxon_from_cleaned_traits"
            else:
                action = "preserve_nonrequired_missingness"
            observations.append(
                TraitMissingObservation(
                    taxon=row["taxon"],
                    row_number=row_number,
                    trait=trait,
                    required_for_analysis=required_for_analysis,
                    action=action,
                )
            )
    return observations


def _apply_branch_length_floor(root: TreeNode, *, floor: float) -> list[str]:
    repaired_nodes: list[str] = []

    def visit(node: TreeNode, *, is_root: bool) -> None:
        if not is_root and node.branch_length is not None and node.branch_length <= 0.0:
            node.branch_length = floor
            repaired_nodes.append(node.name or "<internal>")
        for child in node.children:
            visit(child, is_root=False)

    visit(root, is_root=True)
    return repaired_nodes


def _build_repair_actions(
    *,
    trait_duplicates: list[TraitDuplicateResolution],
    missing_required_trait_taxa: list[str],
    sequence_outlier_taxa: list[str],
    tree_outlier_taxa: list[str],
    repaired_branch_nodes: list[str],
) -> list[DataQualityRepairAction]:
    actions: list[DataQualityRepairAction] = []
    for row in trait_duplicates:
        actions.append(
            DataQualityRepairAction(
                action_kind="resolve_duplicate_trait_taxon",
                affected_taxa=[row.taxon],
                affected_nodes=[],
                reason="multiple raw trait rows were present for one taxon",
                result=(
                    f"kept row {row.selected_row_number} and discarded "
                    + ",".join(str(value) for value in row.discarded_row_numbers)
                ),
            )
        )
    if missing_required_trait_taxa:
        actions.append(
            DataQualityRepairAction(
                action_kind="drop_taxa_for_missing_required_traits",
                affected_taxa=missing_required_trait_taxa,
                affected_nodes=[],
                reason="required comparative traits remained missing after duplicate resolution",
                result="excluded taxa from the cleaned comparative subset",
            )
        )
    if sequence_outlier_taxa:
        actions.append(
            DataQualityRepairAction(
                action_kind="drop_sequence_outlier_taxa",
                affected_taxa=sequence_outlier_taxa,
                affected_nodes=[],
                reason="alignment composition outliers exceeded the governed stress threshold",
                result="excluded taxa from the cleaned alignment and downstream subset",
            )
        )
    if tree_outlier_taxa:
        actions.append(
            DataQualityRepairAction(
                action_kind="drop_tree_branch_outlier_taxa",
                affected_taxa=tree_outlier_taxa,
                affected_nodes=[],
                reason="terminal branches were extreme relative to the raw tree baseline",
                result="excluded taxa from the cleaned tree and linked surfaces",
            )
        )
    if repaired_branch_nodes:
        actions.append(
            DataQualityRepairAction(
                action_kind="floor_nonpositive_branch_lengths",
                affected_taxa=[],
                affected_nodes=repaired_branch_nodes,
                reason="nonpositive branch lengths block comparative and tree-distance assumptions",
                result=f"replaced nonpositive branch lengths with {_TREE_BRANCH_FLOOR}",
            )
        )
    return actions


def _write_workflow_summary_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows = [
        "\t".join(
            [
                "dataset_id",
                "raw_taxon_count",
                "raw_trait_row_count",
                "duplicate_trait_taxon_count",
                "missing_trait_value_count",
                "sequence_outlier_count",
                "tree_zero_length_branch_count",
                "tree_long_branch_outlier_count",
                "dropped_taxon_count",
                "cleaned_taxon_count",
                "repaired_branch_count",
            ]
        ),
        "\t".join(
            [
                report.dataset.dataset_id,
                str(report.dataset.taxon_count),
                str(report.dataset.raw_trait_row_count),
                str(len(report.trait_duplicates)),
                str(len(report.missing_traits)),
                str(len(report.sequence_outliers)),
                str(report.raw_tree_validation.zero_length_branches),
                str(len(report.raw_tree_inspection.long_branch_outliers)),
                str(len(report.dropped_taxa)),
                str(len(report.cleaned_taxa)),
                str(len(report.repaired_branch_nodes)),
            ]
        ),
    ]
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_trait_duplicates_table(
    path: Path,
    rows: list[TraitDuplicateResolution],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "occurrence_count",
            "selected_row_number",
            "selected_non_missing_field_count",
            "discarded_row_numbers",
            "selected_reason",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "occurrence_count": str(row.occurrence_count),
                "selected_row_number": str(row.selected_row_number),
                "selected_non_missing_field_count": str(
                    row.selected_non_missing_field_count
                ),
                "discarded_row_numbers": ",".join(
                    str(value) for value in row.discarded_row_numbers
                ),
                "selected_reason": row.selected_reason,
            }
            for row in rows
        ],
    )


def _write_trait_missing_values_table(
    path: Path,
    rows: list[TraitMissingObservation],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "row_number",
            "trait",
            "required_for_analysis",
            "action",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "row_number": str(row.row_number),
                "trait": row.trait,
                "required_for_analysis": str(row.required_for_analysis).lower(),
                "action": row.action,
            }
            for row in rows
        ],
    )


def _write_sequence_outliers_table(
    path: Path,
    rows: list[SequenceCompositionOutlier],
    *,
    dropped_taxa: set[str],
) -> Path:
    return write_taxon_rows(
        path,
        columns=["taxon", "deviation", "robust_z_score", "action"],
        rows=[
            {
                "taxon": row.identifier,
                "deviation": _format_number(row.deviation),
                "robust_z_score": _format_number(row.robust_z_score),
                "action": (
                    "drop_taxon_from_cleaned_alignment"
                    if row.identifier in dropped_taxa
                    else "flag_only"
                ),
            }
            for row in rows
        ],
    )


def _write_tree_issues_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows: list[dict[str, str]] = []
    if report.raw_tree_validation.zero_length_branches:
        rows.append(
            {
                "issue_code": "zero_length_branches",
                "severity": "warning",
                "affected_taxa": ",".join(report.repaired_branch_nodes),
                "affected_nodes": ",".join(report.repaired_branch_nodes),
                "raw_value": str(report.raw_tree_validation.zero_length_branches),
                "action": "apply_branch_length_floor_in_cleaned_tree",
            }
        )
    for outlier in report.raw_tree_inspection.long_branch_outliers:
        rows.append(
            {
                "issue_code": "long_branch_outlier",
                "severity": "warning",
                "affected_taxa": outlier.node
                if outlier.branch_type == "terminal"
                else "",
                "affected_nodes": outlier.node,
                "raw_value": _format_number(outlier.branch_length),
                "action": (
                    "drop_taxon_from_cleaned_tree"
                    if outlier.branch_type == "terminal"
                    else "flag_only"
                ),
            }
        )
    return write_taxon_rows(
        path,
        columns=[
            "issue_code",
            "severity",
            "affected_taxa",
            "affected_nodes",
            "raw_value",
            "action",
        ],
        rows=rows,
    )


def _write_repair_actions_table(
    path: Path,
    rows: list[DataQualityRepairAction],
) -> Path:
    return write_taxon_rows(
        path,
        columns=["action_kind", "affected_taxa", "affected_nodes", "reason", "result"],
        rows=[
            {
                "action_kind": row.action_kind,
                "affected_taxa": ",".join(row.affected_taxa),
                "affected_nodes": ",".join(row.affected_nodes),
                "reason": row.reason,
                "result": row.result,
            }
            for row in rows
        ],
    )


def _write_cleaned_linkage_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    alignment_taxa = {record.identifier for record in report.cleaned_alignment_records}
    tree_taxa = set(report.cleaned_linkage.usable_taxa)
    trait_taxa = set(report.cleaned_linkage.usable_taxa)
    taxa = sorted(alignment_taxa | tree_taxa | trait_taxa)
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "present_in_tree",
            "present_in_alignment",
            "present_in_traits",
        ],
        rows=[
            {
                "taxon": taxon,
                "present_in_tree": str(taxon in tree_taxa).lower(),
                "present_in_alignment": str(taxon in alignment_taxa).lower(),
                "present_in_traits": str(taxon in trait_taxa).lower(),
            }
            for taxon in taxa
        ],
    )


def _write_cleaned_validation_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    alignment_warning_count = len(
        _substantive_alignment_warnings(report.cleaned_alignment_validation.warnings)
    )
    rows = [
        {
            "surface": "alignment",
            "status": "pass" if alignment_warning_count == 0 else "warning",
            "warning_count": str(alignment_warning_count),
            "detail": (
                f"{report.cleaned_alignment_validation.summary.sequence_count} sequences "
                "remain in the cleaned alignment"
            ),
        },
        {
            "surface": "tree",
            "status": "pass"
            if report.cleaned_tree_validation.biologically_safe
            else "warning",
            "warning_count": str(len(report.cleaned_tree_validation.warnings)),
            "detail": report.cleaned_tree_validation.validity_decision,
        },
        {
            "surface": "traits",
            "status": "pass",
            "warning_count": "0",
            "detail": (
                f"{report.cleaned_trait_validation.row_count} cleaned trait rows keep "
                "all required comparative fields populated"
            ),
        },
        {
            "surface": "linkage",
            "status": "pass",
            "warning_count": "0",
            "detail": (
                f"{report.cleaned_linkage.linked_taxa} taxa overlap exactly across the "
                "cleaned tree and trait table"
            ),
        },
    ]
    return write_taxon_rows(
        path,
        columns=["surface", "status", "warning_count", "detail"],
        rows=rows,
    )


def _substantive_alignment_warnings(warnings: list[str]) -> list[str]:
    ignored = {
        "automatic sequence type defaults to dna from nucleotide-like characters that remain protein-compatible by alphabet alone"
    }
    return [warning for warning in warnings if warning not in ignored]


def _write_overview(
    path: Path,
    dataset: CatarrhineDataQualityStressPanelDataset,
    workflow_bundle: CatarrhineDataQualityStressPanelWorkflowBundle,
) -> Path:
    lines = [
        "# Catarrhine Data Quality Stress Demo",
        "",
        f"- dataset id: `{dataset.dataset_id}`",
        f"- raw taxon count: `{workflow_bundle.raw_taxon_count}`",
        f"- cleaned taxon count: `{workflow_bundle.cleaned_taxon_count}`",
        f"- duplicate trait taxa: `{workflow_bundle.duplicate_trait_taxon_count}`",
        f"- sequence outliers: `{workflow_bundle.sequence_outlier_count}`",
        f"- repaired branch count: `{workflow_bundle.repaired_branch_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- trait duplicates: `{workflow_bundle.trait_duplicates_path.name}`",
        f"- trait missing values: `{workflow_bundle.trait_missing_values_path.name}`",
        f"- sequence outliers: `{workflow_bundle.sequence_outliers_path.name}`",
        f"- tree issues: `{workflow_bundle.tree_issues_path.name}`",
        f"- repair actions: `{workflow_bundle.repair_actions_path.name}`",
        f"- cleaned traits: `{workflow_bundle.cleaned_traits_path.name}`",
        f"- cleaned alignment: `{workflow_bundle.cleaned_alignment_path.name}`",
        f"- cleaned tree: `{workflow_bundle.cleaned_tree_path.name}`",
        f"- cleaned validation: `{workflow_bundle.cleaned_validation_path.name}`",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _format_number(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".12g")


def _resource_root() -> Path:
    return (
        Path(__file__).resolve().parent.parent
        / "resources"
        / "datasets"
        / "stress"
        / _DATASET_ID
    )
