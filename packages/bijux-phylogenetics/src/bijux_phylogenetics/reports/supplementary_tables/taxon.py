from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.core.dataset import (
    DatasetAuditReport,
    DatasetCompletenessRow,
    DatasetCrosswalkRow,
    DatasetExclusionRow,
    audit_dataset_inputs,
)
from bijux_phylogenetics.datasets.study_inputs import (
    TaxonTable,
    load_taxon_table,
    write_taxon_rows,
)
from bijux_phylogenetics.phylo.taxa import (
    TaxonWorkflowLossReport,
    TaxonWorkflowLossRow,
    build_taxon_workflow_loss_report,
)

from .columns import taxon_table_columns
from .models import SupplementaryTaxonTableResult, SupplementaryTaxonTableRow
from .shared import row_lookup, stringify_list, stringify_mapping


def _serialize_row(
    row: SupplementaryTaxonTableRow,
    *,
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> dict[str, object]:
    serialized: dict[str, object] = {
        "taxon": row.taxon,
        "tree_tip_id": row.tree_tip_id or "",
        "alignment_id": row.alignment_id or "",
        "metadata_id": row.metadata_id or "",
        "trait_id": row.trait_id or "",
        "tip_date_id": row.tip_date_id or "",
        "geography_source": row.geography_source or "",
        "calibration_targets": stringify_list(row.calibration_targets),
        "external_taxonomy_ids": stringify_mapping(row.external_taxonomy_ids),
        "analysis_status": row.analysis_status,
        "analysis_exclusion_reason": row.analysis_exclusion_reason or "",
        "analysis_exclusion_causes": stringify_list(row.analysis_exclusion_causes),
        "analysis_first_failed_surface": row.analysis_first_failed_surface or "",
        "affected_analyses": stringify_list(row.affected_analyses),
        "reporting_status": row.reporting_status,
        "reporting_loss_reason": row.reporting_loss_reason or "",
        "workflow_first_loss_stage": row.workflow_first_loss_stage or "",
        "workflow_loss_stages": stringify_list(row.workflow_loss_stages),
        "workflow_loss_reasons": stringify_list(row.workflow_loss_reasons),
        "in_tree": row.in_tree,
        "in_alignment": row.in_alignment,
        "in_metadata": row.in_metadata,
        "in_traits": row.in_traits,
        "in_tip_dates": row.in_tip_dates,
        "in_geography": row.in_geography,
        "in_calibrations": row.in_calibrations,
    }
    for column in metadata_table.columns:
        if column == metadata_table.taxon_column:
            continue
        serialized[f"metadata_{column}"] = row.metadata_values.get(column, "")
    for column in traits_table.columns:
        if column == traits_table.taxon_column:
            continue
        serialized[f"trait_{column}"] = row.trait_values.get(column, "")
    return serialized


def _build_supplementary_row(
    *,
    taxon: str,
    analysis_taxa: set[str],
    crosswalk: DatasetCrosswalkRow | None,
    completeness: DatasetCompletenessRow | None,
    exclusion: DatasetExclusionRow | None,
    workflow_loss: TaxonWorkflowLossRow | None,
    metadata_values: dict[str, str],
    trait_values: dict[str, str],
) -> SupplementaryTaxonTableRow:
    workflow_events = [] if workflow_loss is None else workflow_loss.loss_events
    analysis_exclusion_reason = (
        None if exclusion is None or not exclusion.causes else exclusion.causes[0]
    )
    reporting_loss_reason = (
        None
        if not workflow_events
        or workflow_loss is None
        or workflow_loss.retained_for_reporting
        else f"{workflow_events[0].stage}:{workflow_events[0].reason}"
    )
    return SupplementaryTaxonTableRow(
        taxon=taxon,
        tree_tip_id=None if crosswalk is None else crosswalk.tree_tip,
        alignment_id=None if crosswalk is None else crosswalk.alignment_id,
        metadata_id=None if crosswalk is None else crosswalk.metadata_id,
        trait_id=None if crosswalk is None else crosswalk.trait_id,
        tip_date_id=None if crosswalk is None else crosswalk.tip_date_id,
        geography_source=None if crosswalk is None else crosswalk.geography_source,
        calibration_targets=[] if crosswalk is None else crosswalk.calibration_targets,
        external_taxonomy_ids={}
        if crosswalk is None
        else crosswalk.external_taxonomy_ids,
        analysis_status="included" if taxon in analysis_taxa else "excluded",
        analysis_exclusion_reason=analysis_exclusion_reason,
        analysis_exclusion_causes=[] if exclusion is None else exclusion.causes,
        analysis_first_failed_surface=None
        if exclusion is None
        else exclusion.first_failed_surface,
        affected_analyses=[] if exclusion is None else exclusion.affected_analyses,
        reporting_status="retained"
        if workflow_loss is None or workflow_loss.retained_for_reporting
        else "dropped",
        reporting_loss_reason=reporting_loss_reason,
        workflow_first_loss_stage=None
        if workflow_loss is None
        else workflow_loss.first_loss_stage,
        workflow_loss_stages=[event.stage for event in workflow_events],
        workflow_loss_reasons=[event.reason for event in workflow_events],
        in_tree=False if completeness is None else completeness.in_tree,
        in_alignment=False if completeness is None else completeness.in_alignment,
        in_metadata=False if completeness is None else completeness.in_metadata,
        in_traits=False if completeness is None else completeness.in_traits,
        in_tip_dates=False if completeness is None else completeness.in_tip_dates,
        in_geography=False if completeness is None else completeness.in_geography,
        in_calibrations=False if completeness is None else completeness.in_calibrations,
        metadata_values=metadata_values,
        trait_values=trait_values,
    )


def _build_taxon_rows(
    *,
    dataset_audit: DatasetAuditReport,
    workflow_loss: TaxonWorkflowLossReport,
    metadata_table: TaxonTable,
    traits_table: TaxonTable,
) -> list[SupplementaryTaxonTableRow]:
    analysis_taxa = set(dataset_audit.analysis_taxa)
    crosswalk_by_taxon = {row.taxon: row for row in dataset_audit.crosswalk.rows}
    completeness_by_taxon = {
        row.taxon: row for row in dataset_audit.completeness_matrix.rows
    }
    exclusion_by_taxon = {row.taxon: row for row in dataset_audit.exclusion_table.rows}
    workflow_loss_by_taxon = {row.taxon: row for row in workflow_loss.rows}
    metadata_by_taxon = row_lookup(metadata_table)
    traits_by_taxon = row_lookup(traits_table)
    taxa = sorted(
        set(crosswalk_by_taxon)
        | set(completeness_by_taxon)
        | set(exclusion_by_taxon)
        | set(workflow_loss_by_taxon)
        | set(metadata_by_taxon)
        | set(traits_by_taxon)
    )
    return [
        _build_supplementary_row(
            taxon=taxon,
            analysis_taxa=analysis_taxa,
            crosswalk=crosswalk_by_taxon.get(taxon),
            completeness=completeness_by_taxon.get(taxon),
            exclusion=exclusion_by_taxon.get(taxon),
            workflow_loss=workflow_loss_by_taxon.get(taxon),
            metadata_values=metadata_by_taxon.get(taxon, {}),
            trait_values=traits_by_taxon.get(taxon, {}),
        )
        for taxon in taxa
    ]


def write_supplementary_taxon_table(
    path: Path,
    *,
    tree_path: Path,
    metadata_path: Path,
    traits_path: Path,
    alignment_path: Path | None = None,
    filtered_alignment_path: Path | None = None,
    inference_tree_path: Path | None = None,
    reported_taxa_path: Path | None = None,
    tip_dates_path: Path | None = None,
    calibration_path: Path | None = None,
) -> SupplementaryTaxonTableResult:
    """Write one supplementary taxon table with IDs, annotations, and exclusion evidence."""
    dataset_audit = audit_dataset_inputs(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        tip_dates_path=tip_dates_path,
        calibration_path=calibration_path,
    )
    workflow_loss = build_taxon_workflow_loss_report(
        tree_path,
        metadata_path,
        traits_path,
        alignment_path=alignment_path,
        filtered_alignment_path=filtered_alignment_path,
        inference_tree_path=inference_tree_path,
        reported_taxa_path=reported_taxa_path,
    )
    metadata_table = load_taxon_table(metadata_path)
    traits_table = load_taxon_table(traits_path)
    rows = _build_taxon_rows(
        dataset_audit=dataset_audit,
        workflow_loss=workflow_loss,
        metadata_table=metadata_table,
        traits_table=traits_table,
    )
    columns = taxon_table_columns(metadata_table, traits_table)
    write_taxon_rows(
        path,
        columns=columns,
        rows=[
            _serialize_row(
                row,
                metadata_table=metadata_table,
                traits_table=traits_table,
            )
            for row in rows
        ],
    )
    analysis_included_count = sum(
        1 for row in rows if row.analysis_status == "included"
    )
    reporting_retained_count = sum(
        1 for row in rows if row.reporting_status == "retained"
    )
    return SupplementaryTaxonTableResult(
        output_path=path,
        row_count=len(rows),
        analysis_included_count=analysis_included_count,
        analysis_excluded_count=len(rows) - analysis_included_count,
        reporting_retained_count=reporting_retained_count,
        reporting_dropped_count=len(rows) - reporting_retained_count,
        metadata_column_count=len(metadata_table.columns) - 1,
        trait_column_count=len(traits_table.columns) - 1,
        columns=columns,
        rows=rows,
    )
