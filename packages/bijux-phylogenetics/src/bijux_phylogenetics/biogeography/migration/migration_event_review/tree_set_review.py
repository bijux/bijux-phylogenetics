from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.tree_set.preparation import (
    load_tree_set_trees,
    prepare_analysis_tree_set,
    shared_taxa,
    validate_burnin_fraction,
)
from bijux_phylogenetics.biogeography.state_models import GeographicExcludedTaxonRow
from bijux_phylogenetics.io.newick import dumps_newick

from .contracts import (
    GeographicMigrationTreeRow,
    GeographicMigrationTreeSetEventRow,
    GeographicMigrationTreeSetReport,
    GeographicMigrationTreeSetSummary,
)
from .shared import (
    summarize_tree_set_events,
    tree_set_support_warnings,
)
from .single_tree_review import summarize_geographic_migration_events


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
