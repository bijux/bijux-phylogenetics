from __future__ import annotations

from pathlib import Path
import statistics

from bijux_phylogenetics.ancestral.common import stable_value
from bijux_phylogenetics.biogeography.state_models import (
    summarize_geographic_state_model,
)
from bijux_phylogenetics.io.trees import load_tree

from .contracts import (
    GeographicMigrationEventReport,
    GeographicMigrationEventSummary,
)
from .shared import (
    build_migration_event_rows,
    tree_depth,
)


def summarize_geographic_migration_events(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
) -> GeographicMigrationEventReport:
    """Extract inferred geographic movement events from one rooted tree."""
    base_report = summarize_geographic_state_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_regions=allowed_regions,
    )
    tree = load_tree(tree_path)
    event_rows = build_migration_event_rows(base_report, tree)
    tree_depth_value = tree_depth(tree)
    supports = [row.support for row in event_rows]
    summary = GeographicMigrationEventSummary(
        trait=base_report.trait,
        taxon_column=base_report.taxon_column,
        model=base_report.model,
        internal_model=base_report.internal_model,
        likelihood_method=base_report.likelihood_method,
        analyzed_taxon_count=base_report.summary.analyzed_taxon_count,
        excluded_taxon_count=base_report.summary.excluded_taxon_count,
        tree_depth=tree_depth_value,
        event_count=len(event_rows),
        strongly_supported_event_count=sum(
            row.strongly_supported for row in event_rows
        ),
        mean_event_support=stable_value(statistics.fmean(supports))
        if supports
        else 0.0,
        earliest_midpoint_depth=(
            stable_value(min(row.midpoint_depth for row in event_rows))
            if event_rows
            else None
        ),
        latest_midpoint_depth=(
            stable_value(max(row.midpoint_depth for row in event_rows))
            if event_rows
            else None
        ),
        warning_count=len(base_report.warnings),
    )
    return GeographicMigrationEventReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=base_report.trait,
        taxon_column=base_report.taxon_column,
        model=base_report.model,
        internal_model=base_report.internal_model,
        likelihood_method=base_report.likelihood_method,
        summary=summary,
        event_rows=event_rows,
        exclusion_rows=list(base_report.exclusion_rows),
        warnings=list(base_report.warnings),
    )
