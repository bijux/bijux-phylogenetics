from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.datasets.study_inputs import (
    TraitLinkageReport,
    link_tree_to_traits,
)
from bijux_phylogenetics.runtime.errors import MetadataJoinError

from ..models import (
    CatarrhineDataQualityStressPanelDataset,
    TraitDuplicateResolution,
    TraitMissingObservation,
)
from ..traits import (
    detect_missing_traits,
    load_permissive_trait_rows,
    resolve_duplicate_traits,
    selected_trait_rows,
)


@dataclass(slots=True)
class TraitCleanupSurfaces:
    raw_trait_mismatch_linkage: TraitLinkageReport
    raw_trait_mismatch_error: str | None
    trait_duplicates: list[TraitDuplicateResolution]
    missing_traits: list[TraitMissingObservation]
    cleaned_trait_candidates: list[dict[str, str]]
    missing_required_trait_taxa: list[str]


def prepare_trait_surfaces(
    *,
    dataset: CatarrhineDataQualityStressPanelDataset,
) -> TraitCleanupSurfaces:
    raw_trait_mismatch_linkage = link_tree_to_traits(
        dataset.raw_tree_path,
        dataset.raw_trait_mismatch_path,
        strict=False,
    )
    raw_trait_mismatch_error: str | None = None
    try:
        link_tree_to_traits(
            dataset.raw_tree_path,
            dataset.raw_trait_mismatch_path,
            strict=True,
        )
    except MetadataJoinError as error:
        raw_trait_mismatch_error = str(error)

    raw_trait_rows = load_permissive_trait_rows(dataset.raw_traits_path)
    trait_duplicates = resolve_duplicate_traits(raw_trait_rows)
    duplicate_lookup = {row.taxon: row for row in trait_duplicates}
    cleaned_trait_candidates = selected_trait_rows(raw_trait_rows)
    missing_traits = detect_missing_traits(
        raw_trait_rows,
        required_traits=set(dataset.required_traits),
        duplicate_lookup=duplicate_lookup,
    )
    missing_required_trait_taxa = sorted(
        {
            row.taxon
            for row in missing_traits
            if row.required_for_analysis
            and row.action == "drop_taxon_from_cleaned_traits"
        }
    )
    return TraitCleanupSurfaces(
        raw_trait_mismatch_linkage=raw_trait_mismatch_linkage,
        raw_trait_mismatch_error=raw_trait_mismatch_error,
        trait_duplicates=trait_duplicates,
        missing_traits=missing_traits,
        cleaned_trait_candidates=cleaned_trait_candidates,
        missing_required_trait_taxa=missing_required_trait_taxa,
    )
