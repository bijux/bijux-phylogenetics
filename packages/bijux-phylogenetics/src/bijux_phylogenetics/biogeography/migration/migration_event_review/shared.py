from __future__ import annotations

import statistics

from bijux_phylogenetics.ancestral.common import (
    node_descendant_taxa,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .contracts import (
    GeographicMigrationEventRow,
    GeographicMigrationTreeSetEventRow,
    GeographicMigrationTreeSetEventSummaryRow,
)


def build_migration_event_rows(
    base_report: GeographicStateModelReport,
    tree: PhyloTree,
) -> list[GeographicMigrationEventRow]:
    node_by_signature = {node_signature(node): node for node in tree.iter_nodes()}
    depth_by_node = node_depths(tree)
    rows: list[GeographicMigrationEventRow] = []
    for event in base_report.transition_event_rows:
        if not event.changed:
            continue
        child_node = node_by_signature[event.child_node]
        child_taxa = node_descendant_taxa(child_node)
        parent_depth = depth_by_node[event.parent_node]
        child_depth = depth_by_node[event.child_node]
        branch_length = stable_value(child_depth - parent_depth)
        rows.append(
            GeographicMigrationEventRow(
                branch_id=event.child_node,
                parent_node=event.parent_node,
                child_node=event.child_node,
                child_descendant_taxa=child_taxa,
                branch_length=branch_length,
                parent_depth=parent_depth,
                child_depth=child_depth,
                midpoint_depth=stable_value((parent_depth + child_depth) / 2.0),
                source_region=event.source_region,
                target_region=event.target_region,
                support=stable_value(event.support),
                strongly_supported=event.strongly_supported,
                confidence_class=classify_support(event.support),
            )
        )
    return sorted(
        rows,
        key=lambda row: (
            row.midpoint_depth,
            row.source_region,
            row.target_region,
            row.branch_id,
        ),
    )


def node_depths(tree: PhyloTree) -> dict[str, float]:
    depths = {node_signature(tree.root): 0.0}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            child_depth = stable_value(depth + float(child.branch_length or 0.0))
            depths[node_signature(child)] = child_depth
            visit(child, child_depth)

    visit(tree.root, 0.0)
    return depths


def tree_depth(tree: PhyloTree) -> float:
    return stable_value(
        max((distance or 0.0) for distance in tree.root_to_tip_lengths())
    )


def classify_support(support: float) -> str:
    if support >= 0.9:
        return "strong"
    if support >= 0.6:
        return "moderate"
    return "weak"


def summarize_tree_set_events(
    rows: list[GeographicMigrationTreeSetEventRow],
    *,
    kept_tree_count: int,
) -> list[GeographicMigrationTreeSetEventSummaryRow]:
    grouped: dict[tuple[str, str, str], list[GeographicMigrationTreeSetEventRow]] = {}
    for row in rows:
        grouped.setdefault(
            (row.branch_id, row.source_region, row.target_region),
            [],
        ).append(row)
    summaries: list[GeographicMigrationTreeSetEventSummaryRow] = []
    for (branch_id, source_region, target_region), event_rows in sorted(
        grouped.items()
    ):
        presence_fraction = stable_value(len(event_rows) / kept_tree_count)
        strongly_supported_tree_count = sum(
            row.strongly_supported for row in event_rows
        )
        strongly_supported_tree_fraction = stable_value(
            strongly_supported_tree_count / len(event_rows)
        )
        midpoint_depths = [row.midpoint_depth for row in event_rows]
        mean_support = stable_value(statistics.fmean(row.support for row in event_rows))
        if presence_fraction < 1.0:
            stability_class = "topology_sensitive"
        elif strongly_supported_tree_fraction < 0.5 or mean_support < 0.5:
            stability_class = "low_support"
        else:
            stability_class = "stable"
        summaries.append(
            GeographicMigrationTreeSetEventSummaryRow(
                branch_id=branch_id,
                child_descendant_taxa=event_rows[0].child_descendant_taxa,
                source_region=source_region,
                target_region=target_region,
                tree_presence_count=len(event_rows),
                tree_presence_fraction=presence_fraction,
                strongly_supported_tree_count=strongly_supported_tree_count,
                strongly_supported_tree_fraction=strongly_supported_tree_fraction,
                mean_support=mean_support,
                lower_95_midpoint_depth=stable_value(
                    empirical_quantile(midpoint_depths, 0.025)
                ),
                upper_95_midpoint_depth=stable_value(
                    empirical_quantile(midpoint_depths, 0.975)
                ),
                minimum_parent_depth=stable_value(
                    min(row.parent_depth for row in event_rows)
                ),
                maximum_child_depth=stable_value(
                    max(row.child_depth for row in event_rows)
                ),
                stability_class=stability_class,
            )
        )
    return summaries


def tree_set_support_warnings(
    summaries: list[GeographicMigrationTreeSetEventSummaryRow],
) -> list[str]:
    warnings: list[str] = []
    if any(row.stability_class == "topology_sensitive" for row in summaries):
        warnings.append(
            "one or more inferred geographic movement events are topology-sensitive across retained trees"
        )
    if any(row.stability_class == "low_support" for row in summaries):
        warnings.append(
            "one or more inferred geographic movement events remain weakly supported across retained trees"
        )
    return warnings


def stringify_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)


def empirical_quantile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    index = (len(ordered) - 1) * probability
    lower = int(index)
    upper = min(lower + 1, len(ordered) - 1)
    if lower == upper:
        return ordered[lower]
    fraction = index - lower
    return ordered[lower] + ((ordered[upper] - ordered[lower]) * fraction)
