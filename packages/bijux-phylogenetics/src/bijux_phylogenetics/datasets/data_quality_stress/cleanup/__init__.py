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
from .workflow_report import build_workflow_report


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
    return build_workflow_report(
        dataset=dataset,
        sequence_surfaces=sequence_surfaces,
        trait_surfaces=trait_surfaces,
        cleaned_subset=cleaned_subset,
        dropped_taxa=excluded_taxa,
        repair_actions=repair_actions,
    )
