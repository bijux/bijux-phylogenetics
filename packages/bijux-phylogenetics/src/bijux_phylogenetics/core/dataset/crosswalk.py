from __future__ import annotations

from pathlib import Path

from .context import _collect_external_ids, _load_dataset_context, _ordered_taxa
from .models import (
    DatasetCompletenessMatrix,
    DatasetCompletenessRow,
    DatasetCrosswalkReport,
    DatasetCrosswalkRow,
    DatasetMismatchReport,
    DatasetMismatchRow,
)


def build_dataset_crosswalk(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> DatasetCrosswalkReport:
    """Generate an explicit taxon crosswalk across the main dataset surfaces."""
    context = _load_dataset_context(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    metadata_taxa = set(_ordered_taxa(context.metadata_table))
    trait_taxa = set(_ordered_taxa(context.traits_table))
    alignment_taxa = set(context.alignment_ids)
    tip_date_taxa = set(context.tip_date_taxa)
    tree_taxa = set(context.tree_taxa)
    union_taxa = sorted(
        tree_taxa | metadata_taxa | trait_taxa | alignment_taxa | tip_date_taxa
    )
    rows = [
        DatasetCrosswalkRow(
            taxon=taxon,
            tree_tip=taxon if taxon in tree_taxa else None,
            alignment_id=taxon if taxon in alignment_taxa else None,
            metadata_id=taxon if taxon in metadata_taxa else None,
            trait_id=taxon if taxon in trait_taxa else None,
            tip_date_id=taxon if taxon in tip_date_taxa else None,
            geography_source="geography" if taxon in context.geography_taxa else None,
            calibration_targets=context.calibration_taxa_to_targets.get(taxon, []),
            external_taxonomy_ids=_collect_external_ids(
                taxon,
                context.metadata_table,
                context.traits_table,
            ),
        )
        for taxon in union_taxa
    ]
    return DatasetCrosswalkReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
        rows=rows,
    )


def build_dataset_completeness_matrix(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> DatasetCompletenessMatrix:
    """Build a taxon-by-surface completeness matrix for one dataset."""
    crosswalk = build_dataset_crosswalk(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    rows = [
        DatasetCompletenessRow(
            taxon=row.taxon,
            in_tree=row.tree_tip is not None,
            in_alignment=row.alignment_id is not None,
            in_metadata=row.metadata_id is not None,
            in_traits=row.trait_id is not None,
            in_tip_dates=row.tip_date_id is not None,
            in_geography=row.geography_source is not None,
            in_calibrations=bool(row.calibration_targets),
        )
        for row in crosswalk.rows
    ]
    surface_counts = {
        "tree": sum(1 for row in rows if row.in_tree),
        "alignment": sum(1 for row in rows if row.in_alignment),
        "metadata": sum(1 for row in rows if row.in_metadata),
        "traits": sum(1 for row in rows if row.in_traits),
        "tip_dates": sum(1 for row in rows if row.in_tip_dates),
        "geography": sum(1 for row in rows if row.in_geography),
        "calibrations": sum(1 for row in rows if row.in_calibrations),
    }
    context = _load_dataset_context(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    return DatasetCompletenessMatrix(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
        geography_columns=context.geography_columns,
        rows=rows,
        surface_counts=surface_counts,
    )


def build_dataset_mismatch_report(
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    *,
    alignment_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> DatasetMismatchReport:
    """Show which taxa are missing from which requested dataset surfaces."""
    matrix = build_dataset_completeness_matrix(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    requested_surfaces = ["tree", "metadata", "traits"]
    if alignment_path is not None:
        requested_surfaces.append("alignment")
    if tip_dates_path is not None:
        requested_surfaces.append("tip_dates")
    if calibration_path is not None:
        requested_surfaces.append("calibrations")

    rows: list[DatasetMismatchRow] = []
    mismatch_counts = dict.fromkeys(requested_surfaces, 0)
    for row in matrix.rows:
        present_surfaces: list[str] = []
        missing_surfaces: list[str] = []
        for surface in requested_surfaces:
            present = {
                "tree": row.in_tree,
                "metadata": row.in_metadata,
                "traits": row.in_traits,
                "alignment": row.in_alignment,
                "tip_dates": row.in_tip_dates,
                "calibrations": row.in_calibrations,
            }[surface]
            if present:
                present_surfaces.append(surface)
            else:
                missing_surfaces.append(surface)
                mismatch_counts[surface] += 1
        if missing_surfaces and present_surfaces:
            rows.append(
                DatasetMismatchRow(
                    taxon=row.taxon,
                    present_surfaces=present_surfaces,
                    missing_surfaces=missing_surfaces,
                    message=f"taxon is missing from {', '.join(missing_surfaces)} while present in {', '.join(present_surfaces)}",
                )
            )
    return DatasetMismatchReport(
        requested_surfaces=requested_surfaces,
        rows=rows,
        mismatch_counts=mismatch_counts,
    )
