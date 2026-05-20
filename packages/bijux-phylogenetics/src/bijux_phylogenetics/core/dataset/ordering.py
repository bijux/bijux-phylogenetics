from __future__ import annotations

from pathlib import Path

from .context import _load_dataset_context, _ordered_taxa
from .models import DatasetOrderingAudit, DatasetOrderingConflict


def _ordering_conflicts_for_surface(
    *,
    surface: str,
    canonical_order: list[str],
    observed_order: list[str],
) -> list[DatasetOrderingConflict]:
    canonical_shared = [
        taxon for taxon in canonical_order if taxon in set(observed_order)
    ]
    observed_shared = [
        taxon for taxon in observed_order if taxon in set(canonical_order)
    ]
    if canonical_shared == observed_shared:
        return []

    expected_index = {
        taxon: index for index, taxon in enumerate(canonical_shared, start=1)
    }
    return [
        DatasetOrderingConflict(
            surface=surface,
            taxon=taxon,
            expected_index=expected_index[taxon],
            observed_index=index,
        )
        for index, taxon in enumerate(observed_shared, start=1)
        if expected_index[taxon] != index
    ]


def audit_dataset_taxon_ordering(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
) -> DatasetOrderingAudit:
    """Detect silent taxon-order drift across dataset surfaces."""
    context = _load_dataset_context(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
    )
    conflicts: list[DatasetOrderingConflict] = []
    drifted_surfaces: list[str] = []
    surface_orders = {
        "metadata": _ordered_taxa(context.metadata_table),
        "traits": _ordered_taxa(context.traits_table),
    }
    if alignment_path is not None:
        surface_orders["alignment"] = context.alignment_ids
    if tip_dates_path is not None:
        surface_orders["tip_dates"] = context.tip_date_taxa

    for surface, observed_order in surface_orders.items():
        surface_conflicts = _ordering_conflicts_for_surface(
            surface=surface,
            canonical_order=context.tree_taxa,
            observed_order=observed_order,
        )
        if surface_conflicts:
            drifted_surfaces.append(surface)
            conflicts.extend(surface_conflicts)

    return DatasetOrderingAudit(
        canonical_surface="tree",
        consistent=not conflicts,
        drifted_surfaces=sorted(drifted_surfaces),
        conflicts=conflicts,
    )
