from __future__ import annotations

from pathlib import Path

from ..models import (
    CatarrhineDataQualityStressPanelDataset,
    CatarrhineDataQualityStressPanelWorkflowReport,
)
from .cleaned_subset import assemble_cleaned_subset
from .repair_actions import build_repair_actions
from .sequence_processing import prepare_sequence_surfaces
from .trait_processing import prepare_trait_surfaces


def build_catarrhine_data_quality_stress_panel_workflow_report(
    dataset: CatarrhineDataQualityStressPanelDataset,
    out_dir: Path,
) -> CatarrhineDataQualityStressPanelWorkflowReport:
    out_dir.mkdir(parents=True, exist_ok=True)
    assembled_root = out_dir / "assembled"
    assembled_root.mkdir(parents=True, exist_ok=True)

    sequence_surfaces = prepare_sequence_surfaces(
        dataset=dataset,
        assembled_root=assembled_root,
    )
    trait_surfaces = prepare_trait_surfaces(
        dataset=dataset,
    )
    excluded_taxa = sorted(
        set(sequence_surfaces.sequence_outlier_taxa)
        | set(sequence_surfaces.tree_outlier_taxa)
        | set(trait_surfaces.missing_required_trait_taxa)
    )

    cleaned_subset = assemble_cleaned_subset(
        dataset=dataset,
        assembled_root=assembled_root,
        excluded_taxa=excluded_taxa,
        cleaned_trait_candidates=trait_surfaces.cleaned_trait_candidates,
        validate_cleaned_tree=sequence_surfaces.validate_cleaned_tree,
        validate_cleaned_alignment=sequence_surfaces.validate_cleaned_alignment,
    )

    repair_actions = build_repair_actions(
        raw_sequence_input_repair=sequence_surfaces.raw_sequence_input_repair,
        raw_sequence_length_outliers=sequence_surfaces.raw_sequence_length_outliers,
        coding_sequence_preparation=sequence_surfaces.coding_sequence_preparation,
        raw_trait_mismatch_error=trait_surfaces.raw_trait_mismatch_error,
        trait_duplicates=trait_surfaces.trait_duplicates,
        missing_required_trait_taxa=trait_surfaces.missing_required_trait_taxa,
        sequence_outlier_taxa=sequence_surfaces.sequence_outlier_taxa,
        tree_outlier_taxa=sequence_surfaces.tree_outlier_taxa,
        repaired_branch_nodes=cleaned_subset.repaired_branch_nodes,
    )
    return CatarrhineDataQualityStressPanelWorkflowReport(
        dataset=dataset,
        raw_sequence_input_validation=sequence_surfaces.raw_sequence_input_validation,
        raw_sequence_input_repair=sequence_surfaces.raw_sequence_input_repair,
        raw_sequence_length_outliers=sequence_surfaces.raw_sequence_length_outliers,
        repaired_sequence_input_validation=(
            sequence_surfaces.repaired_sequence_input_validation
        ),
        coding_sequence_preparation=sequence_surfaces.coding_sequence_preparation,
        raw_trait_mismatch_linkage=trait_surfaces.raw_trait_mismatch_linkage,
        raw_trait_mismatch_error=trait_surfaces.raw_trait_mismatch_error,
        raw_alignment_validation=sequence_surfaces.raw_alignment_validation,
        sequence_outliers=sequence_surfaces.sequence_outliers,
        raw_tree_inspection=sequence_surfaces.raw_tree_inspection,
        raw_tree_validation=sequence_surfaces.raw_tree_validation,
        trait_duplicates=trait_surfaces.trait_duplicates,
        missing_traits=trait_surfaces.missing_traits,
        cleaned_trait_validation=cleaned_subset.cleaned_trait_validation,
        cleaned_tree_validation=cleaned_subset.cleaned_tree_validation,
        cleaned_linkage=cleaned_subset.cleaned_linkage,
        cleaned_alignment_validation=cleaned_subset.cleaned_alignment_validation,
        cleaned_alignment_records=cleaned_subset.cleaned_alignment_records,
        repaired_sequence_input_path=sequence_surfaces.repaired_sequence_input_path,
        prepared_coding_sequences_path=sequence_surfaces.prepared_coding_sequences_path,
        cleaned_tree_path=cleaned_subset.cleaned_tree_path,
        cleaned_traits_path=cleaned_subset.cleaned_traits_path,
        cleaned_alignment_path=cleaned_subset.cleaned_alignment_path,
        cleaned_taxa=cleaned_subset.cleaned_taxa,
        dropped_taxa=excluded_taxa,
        repair_actions=repair_actions,
        repaired_branch_nodes=cleaned_subset.repaired_branch_nodes,
    )
