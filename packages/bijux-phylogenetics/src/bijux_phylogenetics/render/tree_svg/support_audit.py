from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.assumptions import standardize_support_labels
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path

from .contracts import SupportLabelRenderAudit


def format_branch_value(value: float) -> str:
    """Format one branch-length or support value for stable SVG rendering."""
    return format(round(value, 15), ".15g")


def coerce_support_label(raw: str | None) -> str | None:
    """Normalize one raw node label into a render-safe support label."""
    if raw is None or not raw.strip():
        return None
    try:
        return format_branch_value(float(raw))
    except ValueError:
        return None


def audit_support_label_rendering(tree_path: Path) -> SupportLabelRenderAudit:
    """Audit whether support labels can be rendered safely for reviewer-facing figures."""
    inspection = inspect_tree_path(tree_path)
    if inspection.suspicious_support_value_ranges:
        return SupportLabelRenderAudit(
            validated=False,
            labels_by_node={},
            warnings=[
                "support labels were withheld because one or more values fall outside interpretable support ranges",
                *inspection.suspicious_support_value_ranges,
            ],
        )

    standardized = standardize_support_labels(tree_path)
    labels_by_node = {
        row.node: format_branch_value(row.support_percent) for row in standardized
    }
    warnings: list[str] = []
    if inspection.mixed_support_scales:
        warnings.append(
            "support labels were standardized from mixed input scales before rendering"
        )
    return SupportLabelRenderAudit(
        validated=True,
        labels_by_node=labels_by_node,
        warnings=warnings,
    )
