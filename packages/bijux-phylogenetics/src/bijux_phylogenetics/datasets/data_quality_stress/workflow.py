from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from .cleanup import build_catarrhine_data_quality_stress_panel_workflow_report
from .export import export_catarrhine_data_quality_stress_panel_dataset
from bijux_phylogenetics.core.alignment import (
    AlignmentRecord,
    CodingSequencePreparationReport,
    FastaInputValidationReport,
    FastaRepairReport,
    SequenceCompositionOutlier,
    SequenceLengthOutlier,
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
from bijux_phylogenetics.runtime.errors import MetadataJoinError
from bijux_phylogenetics.io.fasta import (
    load_fasta_alignment,
    write_fasta_alignment,
)
from bijux_phylogenetics.io.fasta.cleaning import detect_composition_outlier_sequences
from bijux_phylogenetics.io.fasta.coding import prepare_coding_sequences_for_alignment
from bijux_phylogenetics.io.fasta.records import (
    detect_sequence_length_outliers,
    repair_fasta_input,
    validate_fasta_input,
)
from bijux_phylogenetics.io.newick import write_newick
from .panel import (
    load_catarrhine_data_quality_stress_panel_dataset as load_packaged_catarrhine_data_quality_stress_panel_dataset,
)
from .traits import (
    detect_missing_traits,
    load_permissive_trait_rows,
    resolve_duplicate_traits,
    selected_trait_rows,
)

_DATASET_ID = "catarrhine_data_quality_stress_panel"
_DATASET_LABEL = "Catarrhine data quality stress panel"
_RAW_ALIGNMENT_NAME = "alignment.fasta"
_RAW_SEQUENCE_INPUT_NAME = "sequence-input.fasta"
_RAW_CODING_SEQUENCE_NAME = "coding-sequences.fasta"
_RAW_TREE_NAME = "tree.nwk"
_RAW_TRAITS_NAME = "traits.csv"
_RAW_TRAIT_MISMATCH_NAME = "traits-mismatch.csv"
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
    raw_sequence_input_path: Path
    raw_coding_sequences_path: Path
    raw_trait_mismatch_path: Path
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
    raw_sequence_input_path: Path
    raw_coding_sequences_path: Path
    raw_trait_mismatch_path: Path
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
    raw_sequence_input_validation: FastaInputValidationReport
    raw_sequence_input_repair: FastaRepairReport
    raw_sequence_length_outliers: list[SequenceLengthOutlier]
    repaired_sequence_input_validation: FastaInputValidationReport
    coding_sequence_preparation: CodingSequencePreparationReport
    raw_trait_mismatch_linkage: TraitLinkageReport
    raw_trait_mismatch_error: str | None
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
    repaired_sequence_input_path: Path
    prepared_coding_sequences_path: Path
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
    duplicate_sequence_identifier_count: int
    illegal_character_count: int
    empty_sequence_count: int
    raw_sequence_length_outlier_count: int
    duplicate_trait_taxon_count: int
    missing_trait_value_count: int
    sequence_outlier_count: int
    tree_zero_length_branch_count: int
    tree_negative_branch_count: int
    tree_long_branch_outlier_count: int
    coding_frame_error_count: int
    coding_internal_stop_count: int
    raw_trait_missing_from_traits_count: int
    raw_trait_extra_taxon_count: int
    dropped_taxon_count: int
    repaired_branch_count: int
    workflow_summary_path: Path
    raw_sequence_findings_path: Path
    raw_sequence_repair_path: Path
    repaired_sequence_input_path: Path
    repaired_sequence_validation_path: Path
    coding_sequence_exclusions_path: Path
    prepared_coding_sequences_path: Path
    raw_trait_linkage_path: Path
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
    return load_packaged_catarrhine_data_quality_stress_panel_dataset()


def run_catarrhine_data_quality_stress_panel_workflow(
    out_dir: Path,
) -> CatarrhineDataQualityStressPanelWorkflowReport:
    """Audit the raw stress panel and write one cleaned comparative subset."""
    dataset = load_catarrhine_data_quality_stress_panel_dataset()
    return build_catarrhine_data_quality_stress_panel_workflow_report(dataset, out_dir)


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
    raw_sequence_findings_path = _write_raw_sequence_findings_table(
        output_root / "raw-sequence-findings.tsv",
        report,
    )
    raw_sequence_repair_path = _write_raw_sequence_repair_table(
        output_root / "raw-sequence-repair.tsv",
        report,
    )
    repaired_sequence_input_path = _copy_output(
        report.repaired_sequence_input_path,
        output_root / "repaired-sequence-input.fasta",
    )
    repaired_sequence_validation_path = _write_repaired_sequence_validation_table(
        output_root / "repaired-sequence-validation.tsv",
        report,
    )
    coding_sequence_exclusions_path = _write_coding_sequence_exclusions_table(
        output_root / "coding-sequence-exclusions.tsv",
        report.coding_sequence_preparation,
    )
    prepared_coding_sequences_path = _copy_output(
        report.prepared_coding_sequences_path,
        output_root / "prepared-coding-sequences.fasta",
    )
    raw_trait_linkage_path = _write_raw_trait_linkage_table(
        output_root / "raw-trait-linkage.tsv",
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
        duplicate_sequence_identifier_count=len(
            report.raw_sequence_input_validation.duplicate_identifiers
        ),
        illegal_character_count=len(
            report.raw_sequence_input_validation.illegal_characters
        ),
        empty_sequence_count=len(report.raw_sequence_input_validation.empty_sequences),
        raw_sequence_length_outlier_count=len(report.raw_sequence_length_outliers),
        duplicate_trait_taxon_count=len(report.trait_duplicates),
        missing_trait_value_count=len(report.missing_traits),
        sequence_outlier_count=len(report.sequence_outliers),
        tree_zero_length_branch_count=report.raw_tree_validation.zero_length_branches,
        tree_negative_branch_count=report.raw_tree_validation.negative_branch_lengths,
        tree_long_branch_outlier_count=len(
            report.raw_tree_inspection.long_branch_outliers
        ),
        coding_frame_error_count=sum(
            1
            for row in report.coding_sequence_preparation.excluded_sequences
            if row.reason == "frame-error"
        ),
        coding_internal_stop_count=sum(
            1
            for row in report.coding_sequence_preparation.excluded_sequences
            if row.reason == "internal-stop-codon"
        ),
        raw_trait_missing_from_traits_count=len(
            report.raw_trait_mismatch_linkage.missing_from_traits
        ),
        raw_trait_extra_taxon_count=len(
            report.raw_trait_mismatch_linkage.extra_trait_taxa
        ),
        dropped_taxon_count=len(report.dropped_taxa),
        repaired_branch_count=len(report.repaired_branch_nodes),
        workflow_summary_path=workflow_summary_path,
        raw_sequence_findings_path=raw_sequence_findings_path,
        raw_sequence_repair_path=raw_sequence_repair_path,
        repaired_sequence_input_path=repaired_sequence_input_path,
        repaired_sequence_validation_path=repaired_sequence_validation_path,
        coding_sequence_exclusions_path=coding_sequence_exclusions_path,
        prepared_coding_sequences_path=prepared_coding_sequences_path,
        raw_trait_linkage_path=raw_trait_linkage_path,
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


def _write_workflow_summary_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows = [
        "\t".join(
            [
                "dataset_id",
                "raw_taxon_count",
                "duplicate_sequence_identifier_count",
                "illegal_character_count",
                "empty_sequence_count",
                "raw_sequence_length_outlier_count",
                "coding_frame_error_count",
                "coding_internal_stop_count",
                "raw_trait_missing_from_traits_count",
                "raw_trait_extra_taxon_count",
                "raw_trait_row_count",
                "duplicate_trait_taxon_count",
                "missing_trait_value_count",
                "sequence_outlier_count",
                "tree_zero_length_branch_count",
                "tree_negative_branch_count",
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
                str(len(report.raw_sequence_input_validation.duplicate_identifiers)),
                str(len(report.raw_sequence_input_validation.illegal_characters)),
                str(len(report.raw_sequence_input_validation.empty_sequences)),
                str(len(report.raw_sequence_length_outliers)),
                str(
                    sum(
                        1
                        for row in report.coding_sequence_preparation.excluded_sequences
                        if row.reason == "frame-error"
                    )
                ),
                str(
                    sum(
                        1
                        for row in report.coding_sequence_preparation.excluded_sequences
                        if row.reason == "internal-stop-codon"
                    )
                ),
                str(len(report.raw_trait_mismatch_linkage.missing_from_traits)),
                str(len(report.raw_trait_mismatch_linkage.extra_trait_taxa)),
                str(report.dataset.raw_trait_row_count),
                str(len(report.trait_duplicates)),
                str(len(report.missing_traits)),
                str(len(report.sequence_outliers)),
                str(report.raw_tree_validation.zero_length_branches),
                str(report.raw_tree_validation.negative_branch_lengths),
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


def _write_raw_sequence_findings_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows: list[dict[str, str]] = []
    for row in report.raw_sequence_input_validation.duplicate_identifiers:
        rows.append(
            {
                "issue_kind": "duplicate_identifier",
                "identifier": row.identifier,
                "detail": (
                    f"occurrences={row.occurrences};record_indices="
                    + ",".join(str(value) for value in row.record_indices)
                ),
                "action": "normalize_identifier_collision",
            }
        )
    for row in report.raw_sequence_input_validation.illegal_characters:
        rows.append(
            {
                "issue_kind": "illegal_character",
                "identifier": row.identifier,
                "detail": f"record_index={row.record_index};position={row.position};character={row.character}",
                "action": "remove_invalid_record",
            }
        )
    for row in report.raw_sequence_input_validation.empty_sequences:
        rows.append(
            {
                "issue_kind": "empty_sequence",
                "identifier": row.identifier,
                "detail": f"record_index={row.record_index}",
                "action": "remove_invalid_record",
            }
        )
    for row in report.raw_sequence_length_outliers:
        rows.append(
            {
                "issue_kind": "sequence_length_outlier",
                "identifier": row.identifier,
                "detail": (
                    f"raw_length={row.raw_length};median_length="
                    f"{_format_number(row.median_length)};note={row.note}"
                ),
                "action": "drop_length_outlier",
            }
        )
    return write_taxon_rows(
        path,
        columns=["issue_kind", "identifier", "detail", "action"],
        rows=rows,
    )


def _write_raw_sequence_repair_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    rows: list[dict[str, str]] = []
    for row in report.raw_sequence_input_repair.normalized_identifiers:
        rows.append(
            {
                "action_kind": "rename_identifier",
                "identifier": row.original_identifier,
                "result_identifier": row.repaired_identifier,
                "detail": f"record_index={row.record_index};note={row.note}",
            }
        )
    for row in report.raw_sequence_input_repair.removed_records:
        rows.append(
            {
                "action_kind": "remove_record",
                "identifier": row.identifier,
                "result_identifier": "",
                "detail": f"record_index={row.record_index};reason={row.reason}",
            }
        )
    for row in report.raw_sequence_length_outliers:
        rows.append(
            {
                "action_kind": "drop_length_outlier",
                "identifier": row.identifier,
                "result_identifier": "",
                "detail": f"raw_length={row.raw_length};note={row.note}",
            }
        )
    return write_taxon_rows(
        path,
        columns=["action_kind", "identifier", "result_identifier", "detail"],
        rows=rows,
    )


def _write_repaired_sequence_validation_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    substantive_warnings = _substantive_alignment_warnings(
        report.repaired_sequence_input_validation.warnings
    )
    return write_taxon_rows(
        path,
        columns=[
            "surface",
            "sequence_count",
            "duplicate_identifier_count",
            "illegal_character_count",
            "empty_sequence_count",
            "length_outlier_count",
            "warning_count",
            "detail",
        ],
        rows=[
            {
                "surface": "repaired_sequence_input",
                "sequence_count": str(
                    report.repaired_sequence_input_validation.summary.sequence_count
                ),
                "duplicate_identifier_count": str(
                    len(report.repaired_sequence_input_validation.duplicate_identifiers)
                ),
                "illegal_character_count": str(
                    len(report.repaired_sequence_input_validation.illegal_characters)
                ),
                "empty_sequence_count": str(
                    len(report.repaired_sequence_input_validation.empty_sequences)
                ),
                "length_outlier_count": str(
                    len(report.repaired_sequence_input_validation.length_outliers)
                ),
                "warning_count": str(len(substantive_warnings)),
                "detail": (
                    "repaired FASTA input retains only unique identifiers and "
                    "records without illegal characters, empty bodies, or retained "
                    "length outliers"
                ),
            }
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


def _write_coding_sequence_exclusions_table(
    path: Path,
    report: CodingSequencePreparationReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "identifier",
            "reason",
            "comparable_length",
            "invalid_codon_count",
            "premature_stop_count",
            "terminal_stop_count",
            "trailing_bases",
            "note",
        ],
        rows=[
            {
                "identifier": row.identifier,
                "reason": row.reason,
                "comparable_length": str(row.comparable_length),
                "invalid_codon_count": str(row.invalid_codon_count),
                "premature_stop_count": str(row.premature_stop_count),
                "terminal_stop_count": str(row.terminal_stop_count),
                "trailing_bases": str(row.trailing_bases),
                "note": row.note,
            }
            for row in report.excluded_sequences
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
    if report.raw_tree_validation.negative_branch_lengths:
        rows.append(
            {
                "issue_code": "negative_branch_lengths",
                "severity": "warning",
                "affected_taxa": "",
                "affected_nodes": ",".join(
                    _tree_warning_nodes(
                        report.raw_tree_validation,
                        warning_code="negative_branch_lengths",
                    )
                ),
                "raw_value": str(report.raw_tree_validation.negative_branch_lengths),
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


def _write_raw_trait_linkage_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    linkage = report.raw_trait_mismatch_linkage
    return write_taxon_rows(
        path,
        columns=[
            "surface",
            "tree_taxa",
            "trait_taxa",
            "linked_taxa",
            "missing_from_traits",
            "extra_trait_taxa",
            "strict_status",
            "detail",
        ],
        rows=[
            {
                "surface": "raw_trait_mismatch",
                "tree_taxa": str(linkage.tree_taxa),
                "trait_taxa": str(linkage.trait_taxa),
                "linked_taxa": str(linkage.linked_taxa),
                "missing_from_traits": ",".join(linkage.missing_from_traits),
                "extra_trait_taxa": ",".join(linkage.extra_trait_taxa),
                "strict_status": (
                    "failed"
                    if report.raw_trait_mismatch_error is not None
                    else "passed"
                ),
                "detail": report.raw_trait_mismatch_error or "raw linkage passed",
            }
        ],
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


def _tree_warning_nodes(
    report: TreeValidationReport,
    *,
    warning_code: str,
) -> list[str]:
    affected_nodes = [
        node
        for warning in report.warning_details
        if warning.code == warning_code
        for node in warning.affected_nodes
    ]
    return sorted(dict.fromkeys(affected_nodes))


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
        f"- duplicate sequence identifiers: `{workflow_bundle.duplicate_sequence_identifier_count}`",
        f"- illegal FASTA characters: `{workflow_bundle.illegal_character_count}`",
        f"- empty FASTA records: `{workflow_bundle.empty_sequence_count}`",
        f"- raw-sequence length outliers: `{workflow_bundle.raw_sequence_length_outlier_count}`",
        f"- coding frame errors: `{workflow_bundle.coding_frame_error_count}`",
        f"- coding internal stop codons: `{workflow_bundle.coding_internal_stop_count}`",
        f"- duplicate trait taxa: `{workflow_bundle.duplicate_trait_taxon_count}`",
        f"- sequence outliers: `{workflow_bundle.sequence_outlier_count}`",
        f"- repaired branch count: `{workflow_bundle.repaired_branch_count}`",
        "",
        "Generated outputs:",
        "",
        f"- workflow summary: `{workflow_bundle.workflow_summary_path.name}`",
        f"- raw sequence findings: `{workflow_bundle.raw_sequence_findings_path.name}`",
        f"- raw sequence repair ledger: `{workflow_bundle.raw_sequence_repair_path.name}`",
        f"- repaired sequence input: `{workflow_bundle.repaired_sequence_input_path.name}`",
        f"- repaired sequence validation: `{workflow_bundle.repaired_sequence_validation_path.name}`",
        f"- coding sequence exclusions: `{workflow_bundle.coding_sequence_exclusions_path.name}`",
        f"- prepared coding sequences: `{workflow_bundle.prepared_coding_sequences_path.name}`",
        f"- raw trait linkage mismatch: `{workflow_bundle.raw_trait_linkage_path.name}`",
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
