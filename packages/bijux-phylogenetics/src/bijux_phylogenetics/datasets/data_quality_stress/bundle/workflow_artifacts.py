from __future__ import annotations

from pathlib import Path

from ..models import (
    CatarrhineDataQualityStressPanelWorkflowBundle,
    CatarrhineDataQualityStressPanelWorkflowReport,
)


def _build_workflow_bundle(
    *,
    output_root: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
    workflow_summary_path: Path,
    raw_sequence_findings_path: Path,
    raw_sequence_repair_path: Path,
    repaired_sequence_input_path: Path,
    repaired_sequence_validation_path: Path,
    coding_sequence_exclusions_path: Path,
    prepared_coding_sequences_path: Path,
    raw_trait_linkage_path: Path,
    trait_duplicates_path: Path,
    trait_missing_values_path: Path,
    sequence_outliers_path: Path,
    tree_issues_path: Path,
    repair_actions_path: Path,
    cleaned_traits_path: Path,
    cleaned_alignment_path: Path,
    cleaned_tree_path: Path,
    cleaned_linkage_path: Path,
    cleaned_validation_path: Path,
) -> CatarrhineDataQualityStressPanelWorkflowBundle:
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
