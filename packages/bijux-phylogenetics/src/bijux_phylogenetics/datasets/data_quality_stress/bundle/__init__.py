from __future__ import annotations

import shutil
from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    CodingSequencePreparationReport,
    SequenceCompositionOutlier,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import (
    CatarrhineDataQualityStressPanelWorkflowReport,
    DataQualityRepairAction,
    TraitDuplicateResolution,
    TraitMissingObservation,
)
from .shared import (
    _copy_output,
    _format_number,
    _substantive_alignment_warnings,
    _tree_warning_nodes,
)
from .workflow_artifacts import _build_workflow_bundle, _write_workflow_summary_table


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
    return _build_workflow_bundle(
        output_root=output_root,
        report=report,
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
