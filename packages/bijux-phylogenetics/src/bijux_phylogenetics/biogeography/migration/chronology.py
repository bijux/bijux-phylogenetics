from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    node_descendant_taxa,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.biogeography.migration.migration_event_review import (
    GeographicMigrationEventReport,
    summarize_geographic_migration_events,
)
from bijux_phylogenetics.biogeography.state_models import (
    GeographicExcludedTaxonRow,
    GeographicStateModelReport,
    summarize_geographic_state_model,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError


@dataclass(frozen=True, slots=True)
class DatedBiogeographyNodeRow:
    """One node or tip age row in dated biogeography review."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    branch_length: float | None
    depth_from_root: float
    age_before_present: float
    most_likely_region: str
    region_confidence: float
    ambiguous: bool
    is_root: bool


@dataclass(frozen=True, slots=True)
class DatedBiogeographyEventRow:
    """One geographic transition placed into dated-tree age context."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    source_region: str
    target_region: str
    branch_length: float
    parent_depth: float
    child_depth: float
    parent_age_before_present: float
    child_age_before_present: float
    midpoint_age_before_present: float
    time_bin_label: str
    support: float
    strongly_supported: bool
    confidence_class: str


@dataclass(frozen=True, slots=True)
class DatedBiogeographyTimeBinRow:
    """One age-bin summary with explicit support uncertainty."""

    time_bin_label: str
    start_age_before_present: float
    end_age_before_present: float
    event_count: int
    strongly_supported_event_count: int
    low_support_event_count: int
    support_weight_total: float
    mean_support: float | None
    support_uncertainty: float | None
    earliest_event_age_before_present: float | None
    latest_event_age_before_present: float | None
    dominant_transition: str | None
    transition_diversity: int
    uncertainty_class: str


@dataclass(frozen=True, slots=True)
class DatedBiogeographySummary:
    """Reviewer-facing summary for dated biogeography review."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    rooted: bool
    branch_length_status: str
    tree_is_time_scaled: bool
    tip_count: int
    node_age_row_count: int
    root_age: float
    event_count: int
    time_bin_count: int
    empty_time_bin_count: int
    high_uncertainty_bin_count: int
    warning_count: int


@dataclass(slots=True)
class DatedBiogeographyReport:
    """Owned dated-tree biogeography review surface."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    likelihood_method: str
    time_bin_count: int
    summary: DatedBiogeographySummary
    node_rows: list[DatedBiogeographyNodeRow]
    event_rows: list[DatedBiogeographyEventRow]
    time_bin_rows: list[DatedBiogeographyTimeBinRow]
    exclusion_rows: list[GeographicExcludedTaxonRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class _AgeBin:
    label: str
    start_age: float
    end_age: float


def summarize_biogeographic_transition_chronology(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
    time_bin_count: int = 4,
) -> DatedBiogeographyReport:
    """Place inferred geographic transitions into dated-tree age context."""
    if time_bin_count < 1:
        raise ValueError("dated biogeography requires at least one time bin")
    inspection = inspect_tree_path(tree_path)
    if not inspection.rooted:
        raise AncestralReconstructionError(
            "dated biogeography requires a rooted tree with complete branch lengths"
        )
    if inspection.branch_length_status != "complete":
        raise AncestralReconstructionError(
            "dated biogeography requires complete branch lengths"
        )
    if inspection.is_ultrametric is not True:
        raise AncestralReconstructionError(
            "dated biogeography requires an ultrametric time tree"
        )
    tree = load_tree(tree_path)
    root_age = _root_age(tree)
    if root_age <= 0.0:
        raise AncestralReconstructionError(
            "dated biogeography requires positive time depth from root to tips"
        )
    state_report = summarize_geographic_state_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_regions=allowed_regions,
    )
    event_report = summarize_geographic_migration_events(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_regions=allowed_regions,
    )
    depth_by_node = _node_depths(tree)
    observed_tip_regions = _observed_tip_regions(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=state_report.taxon_column,
        trait=trait,
        excluded_taxa={row.taxon for row in state_report.exclusion_rows},
    )
    node_rows = _build_node_rows(
        tree=tree,
        root_age=root_age,
        depth_by_node=depth_by_node,
        state_report=state_report,
        observed_tip_regions=observed_tip_regions,
    )
    age_bins = _build_age_bins(root_age=root_age, time_bin_count=time_bin_count)
    event_rows = _build_event_rows(
        event_report=event_report,
        root_age=root_age,
        age_bins=age_bins,
    )
    time_bin_rows = _build_time_bin_rows(age_bins=age_bins, event_rows=event_rows)
    warnings = list(dict.fromkeys([*state_report.warnings, *event_report.warnings]))
    if any(row.event_count == 0 for row in time_bin_rows):
        warnings.append(
            "one or more dated biogeography age bins contain no inferred geographic transitions"
        )
    if any(
        row.uncertainty_class in {"low_support", "mixed_support"}
        for row in time_bin_rows
    ):
        warnings.append(
            "one or more dated biogeography age bins contain weakly supported or mixed-support transitions"
        )
    summary = DatedBiogeographySummary(
        trait=state_report.trait,
        taxon_column=state_report.taxon_column,
        model=state_report.model,
        internal_model=state_report.internal_model,
        likelihood_method=state_report.likelihood_method,
        analyzed_taxon_count=state_report.summary.analyzed_taxon_count,
        excluded_taxon_count=state_report.summary.excluded_taxon_count,
        rooted=inspection.rooted,
        branch_length_status=inspection.branch_length_status,
        tree_is_time_scaled=inspection.is_ultrametric is True,
        tip_count=inspection.tip_count,
        node_age_row_count=len(node_rows),
        root_age=root_age,
        event_count=len(event_rows),
        time_bin_count=len(time_bin_rows),
        empty_time_bin_count=sum(row.event_count == 0 for row in time_bin_rows),
        high_uncertainty_bin_count=sum(
            row.uncertainty_class in {"low_support", "mixed_support", "no_events"}
            for row in time_bin_rows
        ),
        warning_count=len(warnings),
    )
    return DatedBiogeographyReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=state_report.trait,
        taxon_column=state_report.taxon_column,
        model=state_report.model,
        internal_model=state_report.internal_model,
        likelihood_method=state_report.likelihood_method,
        time_bin_count=time_bin_count,
        summary=summary,
        node_rows=node_rows,
        event_rows=event_rows,
        time_bin_rows=time_bin_rows,
        exclusion_rows=list(state_report.exclusion_rows),
        warnings=warnings,
    )


def write_dated_biogeography_summary_table(
    path: Path,
    report: DatedBiogeographyReport,
) -> Path:
    """Write one overall dated biogeography summary ledger."""
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
            "rooted",
            "branch_length_status",
            "tree_is_time_scaled",
            "tip_count",
            "node_age_row_count",
            "root_age",
            "event_count",
            "time_bin_count",
            "empty_time_bin_count",
            "high_uncertainty_bin_count",
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
                "rooted": str(summary.rooted).lower(),
                "branch_length_status": summary.branch_length_status,
                "tree_is_time_scaled": str(summary.tree_is_time_scaled).lower(),
                "tip_count": str(summary.tip_count),
                "node_age_row_count": str(summary.node_age_row_count),
                "root_age": str(summary.root_age),
                "event_count": str(summary.event_count),
                "time_bin_count": str(summary.time_bin_count),
                "empty_time_bin_count": str(summary.empty_time_bin_count),
                "high_uncertainty_bin_count": str(summary.high_uncertainty_bin_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_dated_biogeography_node_table(
    path: Path,
    report: DatedBiogeographyReport,
) -> Path:
    """Write one node-age ledger for dated biogeography."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "is_tip",
            "descendant_taxa",
            "branch_length",
            "depth_from_root",
            "age_before_present",
            "most_likely_region",
            "region_confidence",
            "ambiguous",
            "is_root",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "is_tip": str(row.is_tip).lower(),
                "descendant_taxa": ",".join(row.descendant_taxa),
                "branch_length": _stringify_optional_float(row.branch_length),
                "depth_from_root": str(row.depth_from_root),
                "age_before_present": str(row.age_before_present),
                "most_likely_region": row.most_likely_region,
                "region_confidence": str(row.region_confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "is_root": str(row.is_root).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_dated_biogeography_event_table(
    path: Path,
    report: DatedBiogeographyReport,
) -> Path:
    """Write one dated-event ledger for geographic transitions."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "source_region",
            "target_region",
            "branch_length",
            "parent_depth",
            "child_depth",
            "parent_age_before_present",
            "child_age_before_present",
            "midpoint_age_before_present",
            "time_bin_label",
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
                "source_region": row.source_region,
                "target_region": row.target_region,
                "branch_length": str(row.branch_length),
                "parent_depth": str(row.parent_depth),
                "child_depth": str(row.child_depth),
                "parent_age_before_present": str(row.parent_age_before_present),
                "child_age_before_present": str(row.child_age_before_present),
                "midpoint_age_before_present": str(row.midpoint_age_before_present),
                "time_bin_label": row.time_bin_label,
                "support": str(row.support),
                "strongly_supported": str(row.strongly_supported).lower(),
                "confidence_class": row.confidence_class,
            }
            for row in report.event_rows
        ],
    )


def write_dated_biogeography_time_bin_table(
    path: Path,
    report: DatedBiogeographyReport,
) -> Path:
    """Write one age-bin uncertainty ledger for dated biogeography."""
    return write_taxon_rows(
        path,
        columns=[
            "time_bin_label",
            "start_age_before_present",
            "end_age_before_present",
            "event_count",
            "strongly_supported_event_count",
            "low_support_event_count",
            "support_weight_total",
            "mean_support",
            "support_uncertainty",
            "earliest_event_age_before_present",
            "latest_event_age_before_present",
            "dominant_transition",
            "transition_diversity",
            "uncertainty_class",
        ],
        rows=[
            {
                "time_bin_label": row.time_bin_label,
                "start_age_before_present": str(row.start_age_before_present),
                "end_age_before_present": str(row.end_age_before_present),
                "event_count": str(row.event_count),
                "strongly_supported_event_count": str(
                    row.strongly_supported_event_count
                ),
                "low_support_event_count": str(row.low_support_event_count),
                "support_weight_total": str(row.support_weight_total),
                "mean_support": _stringify_optional_float(row.mean_support),
                "support_uncertainty": _stringify_optional_float(
                    row.support_uncertainty
                ),
                "earliest_event_age_before_present": _stringify_optional_float(
                    row.earliest_event_age_before_present
                ),
                "latest_event_age_before_present": _stringify_optional_float(
                    row.latest_event_age_before_present
                ),
                "dominant_transition": row.dominant_transition or "",
                "transition_diversity": str(row.transition_diversity),
                "uncertainty_class": row.uncertainty_class,
            }
            for row in report.time_bin_rows
        ],
    )


def write_dated_biogeography_exclusion_table(
    path: Path,
    report: DatedBiogeographyReport,
) -> Path:
    """Write one exclusion ledger for dated biogeography."""
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
            for row in report.exclusion_rows
        ],
    )


def _node_depths(tree: PhyloTree) -> dict[str, float]:
    depths: dict[str, float] = {}

    def visit(node: TreeNode, depth: float) -> None:
        depths[node_signature(node)] = stable_value(depth)
        for child in node.children:
            visit(child, depth + float(child.branch_length or 0.0))

    visit(tree.root, 0.0)
    return depths


def _root_age(tree: PhyloTree) -> float:
    distances = [
        distance for _tip, distance in tree.root_to_tip_pairs() if distance is not None
    ]
    if not distances:
        raise AncestralReconstructionError(
            "dated biogeography requires complete root-to-tip distances"
        )
    return stable_value(max(distances))


def _observed_tip_regions(
    *,
    tree_path: Path,
    traits_path: Path,
    taxon_column: str,
    trait: str,
    excluded_taxa: set[str],
) -> dict[str, str]:
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    tree = load_tree(tree_path)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    observed: dict[str, str] = {}
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None or taxon in excluded_taxa:
            continue
        raw_state = row[trait].strip()
        if raw_state:
            observed[taxon] = raw_state
    return observed


def _build_node_rows(
    *,
    tree: PhyloTree,
    root_age: float,
    depth_by_node: dict[str, float],
    state_report: GeographicStateModelReport,
    observed_tip_regions: dict[str, str],
) -> list[DatedBiogeographyNodeRow]:
    internal_by_node = {row.node: row for row in state_report.node_rows}
    rows: list[DatedBiogeographyNodeRow] = []
    for node in tree.iter_nodes():
        node_id = node_signature(node)
        depth = depth_by_node[node_id]
        age = stable_value(root_age - depth)
        if node.is_leaf():
            region = observed_tip_regions.get(node.name or "", "")
            confidence = 1.0 if region else 0.0
            ambiguous = False
        else:
            region_row = internal_by_node[node_id]
            region = region_row.most_likely_region
            confidence = region_row.confidence
            ambiguous = region_row.ambiguous
        rows.append(
            DatedBiogeographyNodeRow(
                node=node_id,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                branch_length=(
                    None
                    if node is tree.root
                    else stable_value(node.branch_length or 0.0)
                ),
                depth_from_root=depth,
                age_before_present=age,
                most_likely_region=region,
                region_confidence=stable_value(confidence),
                ambiguous=ambiguous,
                is_root=node is tree.root,
            )
        )
    rows.sort(
        key=lambda row: (
            -row.age_before_present,
            row.is_tip,
            row.node,
        )
    )
    return rows


def _build_age_bins(*, root_age: float, time_bin_count: int) -> list[_AgeBin]:
    width = root_age / time_bin_count
    bins: list[_AgeBin] = []
    for index in range(time_bin_count):
        start_age = stable_value(index * width)
        end_age = stable_value(
            root_age if index == time_bin_count - 1 else (index + 1) * width
        )
        bins.append(
            _AgeBin(
                label=f"{start_age:g}-{end_age:g}",
                start_age=start_age,
                end_age=end_age,
            )
        )
    return bins


def _assign_age_bin(age_bins: list[_AgeBin], midpoint_age: float) -> _AgeBin:
    for row in age_bins[:-1]:
        if row.start_age <= midpoint_age < row.end_age:
            return row
    return age_bins[-1]


def _build_event_rows(
    *,
    event_report: GeographicMigrationEventReport,
    root_age: float,
    age_bins: list[_AgeBin],
) -> list[DatedBiogeographyEventRow]:
    rows: list[DatedBiogeographyEventRow] = []
    for row in event_report.event_rows:
        parent_age = stable_value(root_age - row.parent_depth)
        child_age = stable_value(root_age - row.child_depth)
        midpoint_age = stable_value(root_age - row.midpoint_depth)
        age_bin = _assign_age_bin(age_bins, midpoint_age)
        rows.append(
            DatedBiogeographyEventRow(
                branch_id=row.branch_id,
                parent_node=row.parent_node,
                child_node=row.child_node,
                child_descendant_taxa=list(row.child_descendant_taxa),
                source_region=row.source_region,
                target_region=row.target_region,
                branch_length=row.branch_length,
                parent_depth=row.parent_depth,
                child_depth=row.child_depth,
                parent_age_before_present=parent_age,
                child_age_before_present=child_age,
                midpoint_age_before_present=midpoint_age,
                time_bin_label=age_bin.label,
                support=row.support,
                strongly_supported=row.strongly_supported,
                confidence_class=row.confidence_class,
            )
        )
    rows.sort(key=lambda row: (-row.midpoint_age_before_present, row.branch_id))
    return rows


def _build_time_bin_rows(
    *,
    age_bins: list[_AgeBin],
    event_rows: list[DatedBiogeographyEventRow],
) -> list[DatedBiogeographyTimeBinRow]:
    rows_by_bin: dict[str, list[DatedBiogeographyEventRow]] = {
        row.label: [] for row in age_bins
    }
    for row in event_rows:
        rows_by_bin[row.time_bin_label].append(row)
    built: list[DatedBiogeographyTimeBinRow] = []
    for age_bin in age_bins:
        members = rows_by_bin[age_bin.label]
        supports = [row.support for row in members]
        transitions = Counter(
            f"{row.source_region}->{row.target_region}" for row in members
        )
        mean_support = stable_value(sum(supports) / len(supports)) if supports else None
        support_uncertainty = (
            stable_value(1.0 - mean_support) if mean_support is not None else None
        )
        low_support_event_count = sum(row.support < 0.75 for row in members)
        if not members:
            uncertainty_class = "no_events"
        elif (
            mean_support is not None
            and mean_support >= 0.9
            and low_support_event_count == 0
        ):
            uncertainty_class = "stable"
        elif mean_support is not None and mean_support >= 0.75:
            uncertainty_class = "mixed_support"
        else:
            uncertainty_class = "low_support"
        dominant_transition = None
        if transitions:
            dominant_transition = sorted(
                transitions.items(),
                key=lambda item: (-item[1], item[0]),
            )[0][0]
        built.append(
            DatedBiogeographyTimeBinRow(
                time_bin_label=age_bin.label,
                start_age_before_present=age_bin.start_age,
                end_age_before_present=age_bin.end_age,
                event_count=len(members),
                strongly_supported_event_count=sum(
                    row.strongly_supported for row in members
                ),
                low_support_event_count=low_support_event_count,
                support_weight_total=stable_value(sum(supports)) if supports else 0.0,
                mean_support=mean_support,
                support_uncertainty=support_uncertainty,
                earliest_event_age_before_present=(
                    stable_value(
                        max(row.midpoint_age_before_present for row in members)
                    )
                    if members
                    else None
                ),
                latest_event_age_before_present=(
                    stable_value(
                        min(row.midpoint_age_before_present for row in members)
                    )
                    if members
                    else None
                ),
                dominant_transition=dominant_transition,
                transition_diversity=len(transitions),
                uncertainty_class=uncertainty_class,
            )
        )
    return built


def _stringify_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(stable_value(value))
