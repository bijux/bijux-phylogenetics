from __future__ import annotations

from bijux_phylogenetics.ancestral.common import node_descendant_taxa, node_signature
from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.io.newick import loads_newick

from .contracts import (
    AncestralContinuousChangeBranchRow,
    AncestralContinuousChangeCountRow,
)


def summarize_continuous_change_branches(
    report: ContinuousAncestralReport,
) -> list[AncestralContinuousChangeBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimates_by_node = {estimate.node: estimate for estimate in report.estimates}
    traversed_pairs: list[tuple[object, object]] = []

    def collect_pairs(parent_node) -> None:
        for child in parent_node.children:
            traversed_pairs.append((parent_node, child))
            collect_pairs(child)

    collect_pairs(tree.root)
    all_deltas = [
        abs(
            estimates_by_node[node_signature(child_node)].estimate
            - estimates_by_node[node_signature(parent_node)].estimate
        )
        for parent_node, child_node in traversed_pairs
    ]
    max_delta = max(all_deltas) if all_deltas else 0.0
    tolerance = max(1e-9, max_delta * 1e-6)
    rows: list[AncestralContinuousChangeBranchRow] = []
    for parent_node, child_node in traversed_pairs:
        child_signature = node_signature(child_node)
        parent_signature = node_signature(parent_node)
        child_estimate = estimates_by_node[child_signature].estimate
        parent_estimate = estimates_by_node[parent_signature].estimate
        delta = child_estimate - parent_estimate
        if abs(delta) <= tolerance:
            direction = "stable"
        elif delta > 0.0:
            direction = "increase"
        else:
            direction = "decrease"
        rows.append(
            AncestralContinuousChangeBranchRow(
                parent_node=parent_signature,
                child_node=child_signature,
                child_descendant_taxa=tuple(node_descendant_taxa(child_node)),
                branch_length=child_node.branch_length,
                parent_estimate=parent_estimate,
                child_estimate=child_estimate,
                delta=delta,
                absolute_delta=abs(delta),
                direction=direction,
            )
        )
    return rows


def summarize_continuous_change_counts(
    rows: list[AncestralContinuousChangeBranchRow],
) -> list[AncestralContinuousChangeCountRow]:
    total = len(rows)
    counts: list[AncestralContinuousChangeCountRow] = []
    for direction in ("increase", "decrease", "stable"):
        direction_rows = [row for row in rows if row.direction == direction]
        if direction_rows:
            deltas = [row.delta for row in direction_rows]
            minimum_delta = min(deltas)
            maximum_delta = max(deltas)
            mean_delta = sum(deltas) / len(deltas)
        else:
            minimum_delta = 0.0
            maximum_delta = 0.0
            mean_delta = 0.0
        counts.append(
            AncestralContinuousChangeCountRow(
                direction=direction,
                branch_count=len(direction_rows),
                branch_fraction=0.0 if total == 0 else len(direction_rows) / total,
                mean_delta=mean_delta,
                minimum_delta=minimum_delta,
                maximum_delta=maximum_delta,
            )
        )
    return counts
