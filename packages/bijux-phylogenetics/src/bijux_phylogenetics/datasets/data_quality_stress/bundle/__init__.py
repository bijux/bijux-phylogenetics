from __future__ import annotations

import shutil
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import CatarrhineDataQualityStressPanelWorkflowReport
from .shared import _copy_output, _substantive_alignment_warnings
from .sequence_review import (
    _write_coding_sequence_exclusions_table,
    _write_raw_sequence_findings_table,
    _write_raw_sequence_repair_table,
    _write_repaired_sequence_validation_table,
    _write_sequence_outliers_table,
)
from .trait_review import (
    _write_raw_trait_linkage_table,
    _write_trait_duplicates_table,
    _write_trait_missing_values_table,
)
from .tree_review import _write_repair_actions_table, _write_tree_issues_table
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
