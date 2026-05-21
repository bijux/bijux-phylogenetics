from __future__ import annotations

from pathlib import Path
import statistics
import tempfile

from bijux_phylogenetics.ancestral.common import stable_value
from bijux_phylogenetics.ancestral.tree_set.preparation import (
    load_tree_set_trees,
    prepare_analysis_tree_set,
    shared_taxa,
    validate_burnin_fraction,
)
from bijux_phylogenetics.biogeography.state_models import (
    GeographicExcludedTaxonRow,
    summarize_geographic_state_model,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from .migration_event_review import (
    build_migration_event_rows,
    GeographicMigrationEventReport,
    GeographicMigrationEventSummary,
    GeographicMigrationTreeRow,
    GeographicMigrationTreeSetEventRow,
    GeographicMigrationTreeSetReport,
    GeographicMigrationTreeSetSummary,
    stringify_optional_float,
    summarize_tree_set_events,
    tree_depth,
    tree_set_support_warnings,
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


def summarize_geographic_migration_event_tree_set(
    tree_set_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
    burnin_fraction: float = 0.0,
) -> GeographicMigrationTreeSetReport:
    """Extract and summarize inferred geographic movement events across a tree set."""
    validate_burnin_fraction(burnin_fraction)
    _source_format, trees = load_tree_set_trees(tree_set_path)
    total_tree_count = len(trees)
    burnin_tree_count = int(total_tree_count * burnin_fraction)
    kept_tree_entries = [
        (source_tree_index, tree)
        for source_tree_index, tree in enumerate(trees, start=1)
    ][burnin_tree_count:]
    if not kept_tree_entries:
        raise ValueError(
            "geographic movement-event tree-set analysis retains no trees after burn-in removal"
        )
    kept_trees = [tree for _, tree in kept_tree_entries]
    shared_tree_taxa = sorted(shared_taxa(kept_trees))
    warnings: list[str] = []
    if any(set(tree.tip_names) != set(shared_tree_taxa) for tree in kept_trees):
        warnings.append(
            "retained trees do not share identical tip sets and were reduced to their shared taxa"
        )
    (
        analysis_trees,
        topology_summary,
        analysis_taxa,
        raw_exclusions,
        dataset_warnings,
        resolved_taxon_column,
    ) = prepare_analysis_tree_set(
        traits_path=traits_path,
        taxon_column=taxon_column,
        trait=trait,
        kept_tree_entries=kept_tree_entries,
        shared_tree_taxa=shared_tree_taxa,
        dataset_kind="discrete",
    )
    warnings.extend(dataset_warnings)
    tree_rows: list[GeographicMigrationTreeRow] = []
    event_rows: list[GeographicMigrationTreeSetEventRow] = []
    internal_model: str | None = None
    with tempfile.TemporaryDirectory(
        prefix="bijux-phylogenetics-biogeography-event-tree-set-"
    ) as tmp_dir:
        current_tree_path = Path(tmp_dir) / "biogeography-event-current-tree.nwk"
        for (source_tree_index, analysis_tree), topology_record in zip(
            analysis_trees,
            topology_summary.records,
            strict=True,
        ):
            current_tree_path.write_text(
                dumps_newick(analysis_tree) + "\n",
                encoding="utf-8",
            )
            report = summarize_geographic_migration_events(
                current_tree_path,
                traits_path,
                trait=trait,
                taxon_column=resolved_taxon_column,
                model=model,
                allowed_regions=allowed_regions,
            )
            internal_model = report.internal_model
            tree_rows.append(
                GeographicMigrationTreeRow(
                    source_tree_index=source_tree_index,
                    post_burnin_index=topology_record.index,
                    rooted_topology_id=topology_record.rooted_topology_id,
                    unrooted_topology_id=topology_record.unrooted_topology_id,
                    event_count=report.summary.event_count,
                    strongly_supported_event_count=(
                        report.summary.strongly_supported_event_count
                    ),
                )
            )
            event_rows.extend(
                GeographicMigrationTreeSetEventRow(
                    source_tree_index=source_tree_index,
                    post_burnin_index=topology_record.index,
                    rooted_topology_id=topology_record.rooted_topology_id,
                    unrooted_topology_id=topology_record.unrooted_topology_id,
                    branch_id=row.branch_id,
                    parent_node=row.parent_node,
                    child_node=row.child_node,
                    child_descendant_taxa=list(row.child_descendant_taxa),
                    branch_length=row.branch_length,
                    parent_depth=row.parent_depth,
                    child_depth=row.child_depth,
                    midpoint_depth=row.midpoint_depth,
                    source_region=row.source_region,
                    target_region=row.target_region,
                    support=row.support,
                    strongly_supported=row.strongly_supported,
                    confidence_class=row.confidence_class,
                )
                for row in report.event_rows
            )
    event_summaries = summarize_tree_set_events(
        event_rows,
        kept_tree_count=len(analysis_trees),
    )
    warnings.extend(tree_set_support_warnings(event_summaries))
    exclusion_rows = [
        GeographicExcludedTaxonRow(
            taxon=row.taxon,
            raw_region="",
            normalized_region=None,
            reason=row.reason,
            note="",
        )
        for row in raw_exclusions
    ]
    summary = GeographicMigrationTreeSetSummary(
        trait=trait,
        taxon_column=resolved_taxon_column,
        model=model,
        internal_model=internal_model or model,
        total_tree_count=total_tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(analysis_trees),
        shared_tree_taxon_count=len(shared_tree_taxa),
        analysis_taxon_count=len(analysis_taxa),
        rooted_topology_count=len({row.rooted_topology_id for row in tree_rows}),
        unrooted_topology_count=len({row.unrooted_topology_id for row in tree_rows}),
        event_row_count=len(event_rows),
        event_summary_count=len(event_summaries),
        topology_sensitive_event_count=sum(
            row.stability_class == "topology_sensitive" for row in event_summaries
        ),
        low_support_event_count=sum(
            row.stability_class == "low_support" for row in event_summaries
        ),
        excluded_taxon_count=len(exclusion_rows),
        warning_count=len(warnings),
    )
    return GeographicMigrationTreeSetReport(
        tree_set_path=tree_set_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=resolved_taxon_column,
        model=model,
        internal_model=internal_model or model,
        burnin_fraction=burnin_fraction,
        total_tree_count=total_tree_count,
        burnin_tree_count=burnin_tree_count,
        kept_tree_count=len(analysis_trees),
        shared_tree_taxa=shared_tree_taxa,
        analysis_taxa=analysis_taxa,
        rooted_topology_count=summary.rooted_topology_count,
        unrooted_topology_count=summary.unrooted_topology_count,
        summary=summary,
        tree_rows=tree_rows,
        event_rows=event_rows,
        event_summaries=event_summaries,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )


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
    return _write_geographic_exclusion_rows(path, report.exclusion_rows)


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
    return _write_geographic_exclusion_rows(path, report.exclusion_rows)


def _write_geographic_exclusion_rows(
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
