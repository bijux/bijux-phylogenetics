from __future__ import annotations

from ..models import (
    CatarrhineDataQualityStressPanelDataset,
    CatarrhineDataQualityStressPanelWorkflowReport,
    DataQualityRepairAction,
)
from .cleaned_subset import CleanedSubsetSurfaces
from .sequence_processing import SequenceCleanupSurfaces
from .trait_processing import TraitCleanupSurfaces


def build_workflow_report(
    *,
    dataset: CatarrhineDataQualityStressPanelDataset,
    sequence_surfaces: SequenceCleanupSurfaces,
    trait_surfaces: TraitCleanupSurfaces,
    cleaned_subset: CleanedSubsetSurfaces,
    dropped_taxa: list[str],
    repair_actions: list[DataQualityRepairAction],
) -> CatarrhineDataQualityStressPanelWorkflowReport:
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
        dropped_taxa=dropped_taxa,
        repair_actions=repair_actions,
        repaired_branch_nodes=cleaned_subset.repaired_branch_nodes,
    )
