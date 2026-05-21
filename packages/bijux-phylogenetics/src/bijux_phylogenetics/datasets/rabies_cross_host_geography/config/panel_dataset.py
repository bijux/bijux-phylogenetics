from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.records import validate_fasta_input

from ..models import _SOURCE_ACCESSIONS, RabiesCrossHostGeographyPanelDataset
from .audit import (
    _build_workflow_config_audit_rows,
    _raise_for_failed_config_audit,
    _read_observed_groups,
)
from .workflow_config import _load_workflow_config


def load_rabies_cross_host_geography_panel_dataset(
    config_path: Path | None = None,
) -> RabiesCrossHostGeographyPanelDataset:
    """Expose the packaged rabies host-and-geography panel as one owned surface."""
    resolved_config = _load_workflow_config(config_path)
    _raise_for_failed_config_audit(_build_workflow_config_audit_rows(resolved_config))
    dataset_root = resolved_config.config_path.parent
    validation = validate_fasta_input(
        resolved_config.sequences_path,
        sequence_type=resolved_config.sequence_type,
    )
    observed_host_groups, observed_region_groups = _read_observed_groups(
        resolved_config.metadata_path,
        host_trait=resolved_config.host_trait,
        geography_trait=resolved_config.geography_trait,
    )
    return RabiesCrossHostGeographyPanelDataset(
        dataset_id=resolved_config.dataset_id,
        label=resolved_config.label,
        dataset_root=dataset_root,
        workflow_config_path=resolved_config.config_path,
        sequences_path=resolved_config.sequences_path,
        metadata_path=resolved_config.metadata_path,
        centroids_path=resolved_config.centroids_path,
        reference_output_root=dataset_root / "expected",
        sequence_count=validation.summary.sequence_count,
        sequence_type=resolved_config.sequence_type,
        workflow_prefix=resolved_config.workflow_prefix,
        host_trait=resolved_config.host_trait,
        geography_trait=resolved_config.geography_trait,
        host_model=resolved_config.host_model,
        geography_model=resolved_config.geography_model,
        iqtree_seed=resolved_config.iqtree_seed,
        iqtree_threads=resolved_config.iqtree_threads,
        bootstrap_replicates=resolved_config.bootstrap_replicates,
        timeout_seconds=resolved_config.timeout_seconds,
        max_bootstrap_tree_count=resolved_config.max_bootstrap_tree_count,
        max_report_table_rows=resolved_config.max_report_table_rows,
        memory_warning_threshold_bytes=resolved_config.memory_warning_threshold_bytes,
        outgroup_taxa=resolved_config.outgroup_taxa,
        observed_host_group_count=len(observed_host_groups),
        observed_region_group_count=len(observed_region_groups),
        clade_metadata_columns=resolved_config.clade_metadata_columns,
        comparative_formula=resolved_config.comparative_formula,
        comparative_response=resolved_config.comparative_response,
        comparative_branch_length_floor=resolved_config.comparative_branch_length_floor,
        source_accessions=_SOURCE_ACCESSIONS,
        source_summary=(
            "Real rabies virus nucleoprotein sequences paired with grouped host "
            "and macroregion metadata so one governed workflow can rerun tree "
            "inference, host switching, geography review, bootstrap topology "
            "summary, clade extraction, and one comparative model from raw "
            "sequence inputs."
        ),
    )
