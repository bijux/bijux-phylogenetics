from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import statistics
import tempfile

from bijux_phylogenetics.ancestral.common import (
    node_descendant_taxa,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.ancestral.tree_set.preparation import (
    load_tree_set_trees,
    prepare_analysis_tree_set,
    shared_taxa,
    validate_burnin_fraction,
)
from bijux_phylogenetics.biogeography.geographic_states import (
    GeographicExcludedTaxonRow,
    GeographicStateModelReport,
    summarize_geographic_state_model,
)
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree


@dataclass(frozen=True, slots=True)
class GeographicMigrationEventRow:
    """One inferred geographic movement event on one analyzed rooted tree."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float
    parent_depth: float
    child_depth: float
    midpoint_depth: float
    source_region: str
    target_region: str
    support: float
    strongly_supported: bool
    confidence_class: str


@dataclass(frozen=True, slots=True)
class GeographicMigrationEventSummary:
    """Reviewer-facing summary for geographic movement events on one tree."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    tree_depth: float
    event_count: int
    strongly_supported_event_count: int
    mean_event_support: float
    earliest_midpoint_depth: float | None
    latest_midpoint_depth: float | None
    warning_count: int


@dataclass(slots=True)
class GeographicMigrationEventReport:
    """Owned geographic movement-event review surface for one rooted tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    summary: GeographicMigrationEventSummary
    event_rows: list[GeographicMigrationEventRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeRow:
    """One retained tree summary from geographic movement-event tree-set review."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    event_count: int
    strongly_supported_event_count: int


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeSetEventRow:
    """One inferred geographic movement event from one retained tree."""

    source_tree_index: int
    post_burnin_index: int
    rooted_topology_id: str
    unrooted_topology_id: str
    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float
    parent_depth: float
    child_depth: float
    midpoint_depth: float
    source_region: str
    target_region: str
    support: float
    strongly_supported: bool
    confidence_class: str


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeSetEventSummaryRow:
    """One comparable event summary across retained trees."""

    branch_id: str
    child_descendant_taxa: list[str]
    source_region: str
    target_region: str
    tree_presence_count: int
    tree_presence_fraction: float
    strongly_supported_tree_count: int
    strongly_supported_tree_fraction: float
    mean_support: float
    lower_95_midpoint_depth: float
    upper_95_midpoint_depth: float
    minimum_parent_depth: float
    maximum_child_depth: float
    stability_class: str


@dataclass(frozen=True, slots=True)
class GeographicMigrationTreeSetSummary:
    """Reviewer-facing summary for geographic movement events across a tree set."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxon_count: int
    analysis_taxon_count: int
    rooted_topology_count: int
    unrooted_topology_count: int
    event_row_count: int
    event_summary_count: int
    topology_sensitive_event_count: int
    low_support_event_count: int
    excluded_taxon_count: int
    warning_count: int


@dataclass(slots=True)
class GeographicMigrationTreeSetReport:
    """Geographic movement-event stability across a posterior or bootstrap tree set."""

    tree_set_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    kept_tree_count: int
    shared_tree_taxa: list[str]
    analysis_taxa: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    summary: GeographicMigrationTreeSetSummary
    tree_rows: list[GeographicMigrationTreeRow]
    event_rows: list[GeographicMigrationTreeSetEventRow]
    event_summaries: list[GeographicMigrationTreeSetEventSummaryRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


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
    event_rows = _build_migration_event_rows(base_report, tree)
    tree_depth = _tree_depth(tree)
    supports = [row.support for row in event_rows]
    summary = GeographicMigrationEventSummary(
        trait=base_report.trait,
        taxon_column=base_report.taxon_column,
        model=base_report.model,
        internal_model=base_report.internal_model,
        likelihood_method=base_report.likelihood_method,
        analyzed_taxon_count=base_report.summary.analyzed_taxon_count,
        excluded_taxon_count=base_report.summary.excluded_taxon_count,
        tree_depth=tree_depth,
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
    event_summaries = _summarize_tree_set_events(
        event_rows,
        kept_tree_count=len(analysis_trees),
    )
    warnings.extend(_tree_set_support_warnings(event_summaries))
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
                "earliest_midpoint_depth": _stringify_optional_float(
                    summary.earliest_midpoint_depth
                ),
                "latest_midpoint_depth": _stringify_optional_float(
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


def _build_migration_event_rows(
    base_report: GeographicStateModelReport,
    tree: PhyloTree,
) -> list[GeographicMigrationEventRow]:
    node_by_signature = {node_signature(node): node for node in tree.iter_nodes()}
    depth_by_node = _node_depths(tree)
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
                confidence_class=_classify_support(event.support),
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


def _node_depths(tree: PhyloTree) -> dict[str, float]:
    depths = {node_signature(tree.root): 0.0}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            child_depth = stable_value(depth + float(child.branch_length or 0.0))
            depths[node_signature(child)] = child_depth
            visit(child, child_depth)

    visit(tree.root, 0.0)
    return depths


def _tree_depth(tree: PhyloTree) -> float:
    return stable_value(
        max((distance or 0.0) for distance in tree.root_to_tip_lengths())
    )


def _classify_support(support: float) -> str:
    if support >= 0.9:
        return "strong"
    if support >= 0.6:
        return "moderate"
    return "weak"


def _summarize_tree_set_events(
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
                    _empirical_quantile(midpoint_depths, 0.025)
                ),
                upper_95_midpoint_depth=stable_value(
                    _empirical_quantile(midpoint_depths, 0.975)
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


def _tree_set_support_warnings(
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


def _stringify_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)


def _empirical_quantile(values: list[float], probability: float) -> float:
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
