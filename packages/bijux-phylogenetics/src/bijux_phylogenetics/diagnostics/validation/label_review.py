from __future__ import annotations

from bijux_phylogenetics.phylo.taxa import inspect_tree_taxa_safety
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import (
    InternalLabelInterpretation,
    InternalNodeLabelConflict,
    UnsafeExternalLabel,
)
from .structure import _node_signature


def _parse_internal_label_numeric(label: str) -> float | None:
    try:
        return float(label)
    except ValueError:
        return None


def _internal_label_diagnostics(
    tree: PhyloTree,
) -> tuple[
    list[InternalLabelInterpretation],
    list[InternalLabelInterpretation],
    list[str],
    bool,
]:
    likely_support_labels: list[InternalLabelInterpretation] = []
    likely_named_internal_labels: list[InternalLabelInterpretation] = []

    for node in tree.iter_nodes():
        if node.is_leaf() or node.name is None or not node.name.strip():
            continue
        stripped = node.name.strip()
        numeric_value = _parse_internal_label_numeric(stripped)
        if numeric_value is None:
            interpretation = (
                "ambiguous_alphanumeric_label"
                if any(character.isdigit() for character in stripped)
                else "named_internal_label"
            )
        elif numeric_value < 0 or numeric_value > 100:
            interpretation = "out_of_range_support"
        elif 0.0 <= numeric_value <= 1.0:
            interpretation = "fractional_support"
        else:
            interpretation = "percentage_support"
        interpretation_row = InternalLabelInterpretation(
            node=_node_signature(node),
            node_id=_node_signature(node),
            label=node.name,
            interpretation=interpretation,
            numeric_value=None if numeric_value is None else round(numeric_value, 15),
        )
        if numeric_value is not None:
            likely_support_labels.append(interpretation_row)
        else:
            likely_named_internal_labels.append(interpretation_row)

    suspicious_ranges: list[str] = []
    support_values = [
        row.numeric_value
        for row in likely_support_labels
        if row.numeric_value is not None
    ]
    for row in likely_support_labels:
        if row.numeric_value is None:
            continue
        if row.numeric_value < 0:
            suspicious_ranges.append(
                f"support value {row.numeric_value:g} at {row.node} is negative"
            )
        elif row.numeric_value > 100:
            suspicious_ranges.append(
                f"support value {row.numeric_value:g} at {row.node} exceeds 100"
            )

    has_fraction_scale = any(0 <= value <= 1 for value in support_values)
    has_percent_scale = any(1 < value <= 100 for value in support_values)
    mixed_scales = has_fraction_scale and has_percent_scale
    return (
        sorted(likely_support_labels, key=lambda row: (row.node, row.label)),
        sorted(likely_named_internal_labels, key=lambda row: (row.node, row.label)),
        suspicious_ranges,
        mixed_scales,
    )


def _internal_label_conflicts(
    likely_support_labels: list[InternalLabelInterpretation],
    likely_named_internal_labels: list[InternalLabelInterpretation],
    suspicious_support_value_ranges: list[str],
    mixed_support_scales: bool,
) -> list[InternalNodeLabelConflict]:
    conflicts: list[InternalNodeLabelConflict] = []
    for row in likely_support_labels:
        if row.numeric_value is None:
            continue
        if row.numeric_value < 0 or row.numeric_value > 100:
            conflicts.append(
                InternalNodeLabelConflict(
                    node_id=row.node_id,
                    label=row.label,
                    conflict_type="support_out_of_range",
                    detail="numeric internal label falls outside accepted support scales",
                )
            )
        elif 0.0 <= row.numeric_value <= 1.0:
            conflicts.append(
                InternalNodeLabelConflict(
                    node_id=row.node_id,
                    label=row.label,
                    conflict_type="ambiguous_fraction_support",
                    detail="numeric internal label could be posterior support or a named code on a 0-1 scale",
                )
            )
    for row in likely_named_internal_labels:
        stripped = row.label.strip()
        if any(character.isdigit() for character in stripped):
            conflicts.append(
                InternalNodeLabelConflict(
                    node_id=row.node_id,
                    label=row.label,
                    conflict_type="alphanumeric_internal_label",
                    detail="non-numeric internal label contains digits and may conflate clade naming with support annotation",
                )
            )
    if mixed_support_scales:
        conflicts.append(
            InternalNodeLabelConflict(
                node_id="<tree>",
                label="mixed-scales",
                conflict_type="mixed_support_scales",
                detail="tree mixes fraction-like and percentage-like support labels",
            )
        )
    if suspicious_support_value_ranges:
        conflicts.append(
            InternalNodeLabelConflict(
                node_id="<tree>",
                label="support-range",
                conflict_type="support_range_summary",
                detail="one or more internal support-like labels fall outside standard ranges",
            )
        )
    return sorted(
        conflicts, key=lambda row: (row.node_id, row.conflict_type, row.label)
    )


def _unsafe_external_labels(tree: PhyloTree) -> list[UnsafeExternalLabel]:
    report = inspect_tree_taxa_safety(tree, policy="spaces-to-underscores")
    labels: list[UnsafeExternalLabel] = []
    default_engines = ["iqtree", "raxml", "mrbayes", "beast", "r", "shell"]
    for entry in report.unsafe_taxa:
        labels.append(
            UnsafeExternalLabel(
                raw_label=entry.raw_label,
                normalized_label=entry.normalized_label,
                engines=list(default_engines),
                reasons=entry.reasons,
            )
        )
    for collision in report.collisions:
        for raw_label in collision.raw_labels:
            if raw_label in {item.raw_label for item in labels}:
                continue
            labels.append(
                UnsafeExternalLabel(
                    raw_label=raw_label,
                    normalized_label=collision.normalized_label,
                    engines=list(default_engines),
                    reasons=["collides with another label after normalization"],
                )
            )
    return sorted(labels, key=lambda row: row.raw_label)
