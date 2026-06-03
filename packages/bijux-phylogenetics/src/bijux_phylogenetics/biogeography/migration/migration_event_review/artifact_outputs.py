from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.biogeography.state_models import GeographicExcludedTaxonRow
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import (
    GeographicMigrationEventReport,
    GeographicMigrationTreeSetReport,
)
from .shared import stringify_optional_float


def write_geographic_migration_event_summary_table(
    path: Path,
    report: GeographicMigrationEventReport,
) -> Path:
    """Write one summary ledger for geographic movement events on one tree."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "likelihood_method",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "tree_depth",
            "event_count",
            "strongly_supported_event_count",
            "mean_event_support",
            "earliest_midpoint_depth",
            "latest_midpoint_depth",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "likelihood_method": summary.likelihood_method,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "tree_depth": str(summary.tree_depth),
                "event_count": str(summary.event_count),
                "strongly_supported_event_count": str(
                    summary.strongly_supported_event_count
                ),
                "mean_event_support": str(summary.mean_event_support),
                "earliest_midpoint_depth": stringify_optional_float(
                    summary.earliest_midpoint_depth
                ),
                "latest_midpoint_depth": stringify_optional_float(
                    summary.latest_midpoint_depth
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_geographic_migration_event_table(
    path: Path,
    report: GeographicMigrationEventReport,
) -> Path:
    """Write one branchwise geographic movement-event ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_depth",
            "child_depth",
            "midpoint_depth",
            "source_region",
            "target_region",
            "support",
            "strongly_supported",
            "confidence_class",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": str(row.branch_length),
                "parent_depth": str(row.parent_depth),
                "child_depth": str(row.child_depth),
                "midpoint_depth": str(row.midpoint_depth),
                "source_region": row.source_region,
                "target_region": row.target_region,
                "support": str(row.support),
                "strongly_supported": str(row.strongly_supported).lower(),
                "confidence_class": row.confidence_class,
            }
            for row in report.event_rows
        ],
    )


def write_geographic_migration_exclusion_table(
    path: Path,
    report: GeographicMigrationEventReport,
) -> Path:
    """Write one excluded-taxa ledger for one-tree geographic movement events."""
    return write_geographic_exclusion_rows(path, report.exclusion_rows)


def write_geographic_migration_tree_set_summary_table(
    path: Path,
    report: GeographicMigrationTreeSetReport,
) -> Path:
    """Write one overall summary ledger for geographic movement-event tree-set review."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "total_tree_count",
            "burnin_tree_count",
            "kept_tree_count",
            "shared_tree_taxon_count",
            "analysis_taxon_count",
            "rooted_topology_count",
            "unrooted_topology_count",
            "event_row_count",
            "event_summary_count",
            "topology_sensitive_event_count",
            "low_support_event_count",
            "excluded_taxon_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "total_tree_count": str(summary.total_tree_count),
                "burnin_tree_count": str(summary.burnin_tree_count),
                "kept_tree_count": str(summary.kept_tree_count),
                "shared_tree_taxon_count": str(summary.shared_tree_taxon_count),
                "analysis_taxon_count": str(summary.analysis_taxon_count),
                "rooted_topology_count": str(summary.rooted_topology_count),
                "unrooted_topology_count": str(summary.unrooted_topology_count),
                "event_row_count": str(summary.event_row_count),
                "event_summary_count": str(summary.event_summary_count),
                "topology_sensitive_event_count": str(
                    summary.topology_sensitive_event_count
                ),
                "low_support_event_count": str(summary.low_support_event_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_geographic_migration_tree_set_tree_table(
    path: Path,
    report: GeographicMigrationTreeSetReport,
) -> Path:
    """Write one retained-tree ledger for geographic movement-event tree-set review."""
    return write_taxon_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "event_count",
            "strongly_supported_event_count",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "event_count": str(row.event_count),
                "strongly_supported_event_count": str(
                    row.strongly_supported_event_count
                ),
            }
            for row in report.tree_rows
        ],
    )


def write_geographic_migration_tree_set_event_table(
    path: Path,
    report: GeographicMigrationTreeSetReport,
) -> Path:
    """Write one per-tree event ledger for geographic movement-event tree-set review."""
    return write_taxon_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_depth",
            "child_depth",
            "midpoint_depth",
            "source_region",
            "target_region",
            "support",
            "strongly_supported",
            "confidence_class",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": str(row.branch_length),
                "parent_depth": str(row.parent_depth),
                "child_depth": str(row.child_depth),
                "midpoint_depth": str(row.midpoint_depth),
                "source_region": row.source_region,
                "target_region": row.target_region,
                "support": str(row.support),
                "strongly_supported": str(row.strongly_supported).lower(),
                "confidence_class": row.confidence_class,
            }
            for row in report.event_rows
        ],
    )


def write_geographic_migration_tree_set_event_summary_table(
    path: Path,
    report: GeographicMigrationTreeSetReport,
) -> Path:
    """Write one comparable-event summary ledger across retained trees."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "child_descendant_taxa",
            "source_region",
            "target_region",
            "tree_presence_count",
            "tree_presence_fraction",
            "strongly_supported_tree_count",
            "strongly_supported_tree_fraction",
            "mean_support",
            "lower_95_midpoint_depth",
            "upper_95_midpoint_depth",
            "minimum_parent_depth",
            "maximum_child_depth",
            "stability_class",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "source_region": row.source_region,
                "target_region": row.target_region,
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "strongly_supported_tree_count": str(row.strongly_supported_tree_count),
                "strongly_supported_tree_fraction": str(
                    row.strongly_supported_tree_fraction
                ),
                "mean_support": str(row.mean_support),
                "lower_95_midpoint_depth": str(row.lower_95_midpoint_depth),
                "upper_95_midpoint_depth": str(row.upper_95_midpoint_depth),
                "minimum_parent_depth": str(row.minimum_parent_depth),
                "maximum_child_depth": str(row.maximum_child_depth),
                "stability_class": row.stability_class,
            }
            for row in report.event_summaries
        ],
    )


def write_geographic_migration_tree_set_exclusion_table(
    path: Path,
    report: GeographicMigrationTreeSetReport,
) -> Path:
    """Write one excluded-taxa ledger for geographic movement-event tree-set review."""
    return write_geographic_exclusion_rows(path, report.exclusion_rows)


def write_geographic_exclusion_rows(
    path: Path,
    rows: list[GeographicExcludedTaxonRow],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_region",
            "normalized_region",
            "reason",
            "note",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "raw_region": row.raw_region,
                "normalized_region": row.normalized_region or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in rows
        ],
    )
