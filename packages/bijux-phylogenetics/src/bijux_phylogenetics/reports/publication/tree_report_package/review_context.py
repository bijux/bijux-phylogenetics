from __future__ import annotations

from bijux_phylogenetics.diagnostics.validation import TreeInspectionReport
from bijux_phylogenetics.io.iqtree_support import support_fraction
from bijux_phylogenetics.render.tree_svg import SupportLabelRenderAudit
from bijux_phylogenetics.trees import (
    BranchLengthDistributionReport,
    CladeTableReport,
)

from .contracts import TreeBranchStatisticsRow, TreeSupportRow


def support_class(value: float | None) -> str:
    """Classify one support value into reviewer-facing support bands."""
    fraction = support_fraction(value)
    if fraction is None:
        return "missing"
    if fraction >= 0.95:
        return "strong"
    if fraction >= 0.80:
        return "moderate"
    return "weak"


def summarize_tree_support(clades: CladeTableReport) -> list[TreeSupportRow]:
    """Derive reviewer-facing support rows from one clade report."""
    rows: list[TreeSupportRow] = []
    for row in clades.rows:
        if row.node_kind == "tip":
            continue
        rows.append(
            TreeSupportRow(
                node_kind=row.node_kind,
                node=row.clade_id,
                node_label=row.node_label,
                descendant_taxa=tuple(row.taxa),
                support=row.support,
                support_fraction=row.support_fraction,
                support_class=support_class(row.support),
                branch_length=row.branch_length,
                root_depth=row.root_depth,
            )
        )
    return rows


def summarize_tree_branch_statistics(
    branch_lengths: BranchLengthDistributionReport,
) -> TreeBranchStatisticsRow:
    """Derive reviewer-facing branch statistics from one branch report."""
    aggregate = branch_lengths.aggregate
    return TreeBranchStatisticsRow(
        branch_count=aggregate.branch_count,
        defined_branch_count=aggregate.defined_branch_count,
        missing_branch_count=aggregate.missing_branch_count,
        zero_length_branch_count=aggregate.zero_length_branch_count,
        negative_branch_count=aggregate.negative_branch_count,
        positive_branch_count=aggregate.positive_branch_count,
        long_outlier_count=aggregate.long_outlier_count,
        short_outlier_count=aggregate.short_outlier_count,
        minimum_branch_length=aggregate.minimum_branch_length,
        maximum_branch_length=aggregate.maximum_branch_length,
        mean_branch_length=aggregate.mean_branch_length,
        median_branch_length=aggregate.median_branch_length,
        positive_branch_median=aggregate.positive_branch_median,
    )


def build_reviewer_summary(
    *,
    inspection: TreeInspectionReport,
    support_rows: list[TreeSupportRow],
    branch_stats: TreeBranchStatisticsRow,
    support_audit: SupportLabelRenderAudit,
) -> tuple[list[str], list[str]]:
    """Build reviewer summary bullets and explicit tree-report limitations."""
    supported_rows = [row for row in support_rows if row.support is not None]
    strong_rows = [row for row in support_rows if row.support_class == "strong"]
    limitations: list[str] = []
    if not support_audit.validated:
        limitations.append(
            "support labels were withheld from the rendered figure because the input support surface was not safe to standardize"
        )
    if branch_stats.missing_branch_count:
        limitations.append(
            "branch-length summaries include missing lengths, so weighted interpretation is incomplete"
        )
    if branch_stats.negative_branch_count:
        limitations.append(
            "negative branch lengths remain in the source tree and should be corrected before downstream weighted analysis"
        )
    summary = [
        f"tree quality score: {inspection.tree_quality_score}",
        f"tip count: {inspection.tip_count}",
        f"internal clade count: {sum(1 for row in support_rows if row.node_kind == 'internal')}",
        f"supported branch count: {len(supported_rows)}",
        f"strong-support branch count: {len(strong_rows)}",
        f"long-branch outlier count: {branch_stats.long_outlier_count}",
    ]
    if support_audit.validated and support_audit.warnings:
        summary.extend(support_audit.warnings)
    return summary, limitations
