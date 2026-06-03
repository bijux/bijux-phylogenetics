from __future__ import annotations

from pathlib import Path
import shutil

from ..models import (
    CatarrhineDataQualityStressPanelWorkflowBundle,
    CatarrhineDataQualityStressPanelWorkflowReport,
)
from .cleaned_review import (
    _write_cleaned_linkage_table,
    _write_cleaned_validation_table,
)
from .sequence_review import (
    _write_coding_sequence_exclusions_table,
    _write_raw_sequence_findings_table,
    _write_raw_sequence_repair_table,
    _write_repaired_sequence_validation_table,
    _write_sequence_outliers_table,
)
from .shared import _copy_output
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
