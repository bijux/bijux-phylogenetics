from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.core.branching_times import compute_tree_branching_times
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.core.ultrametric import assess_tree_ultrametricity
from bijux_phylogenetics.diagnostics.validation import (
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.runtime.errors import (
    DiversificationAnalysisError,
    NonUltrametricTreeError,
    UnrootedTreeError,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.html import write_html_report

_SAMPLING_COLUMNS = (
    "sampling_fraction",
    "sampling_proportion",
    "sampling_probability",
    "sampling_prob",
)


@dataclass(slots=True)
class TimeTreeValidationReport:
    tree_path: Path
    rooted: bool
    ultrametric: bool
    branch_length_status: str
    tip_count: int
    root_age: float
    warnings: list[str]


@dataclass(slots=True)
class LineageThroughTimePoint:
    node: str
    time_before_present: float
    lineage_count: int
    event: str


@dataclass(slots=True)
class LineageThroughTimeReport:
    tree_path: Path
    tip_count: int
    root_age: float
    points: list[LineageThroughTimePoint]


@dataclass(slots=True)
class SamplingFractionIssue:
    taxon: str
    code: str
    raw_value: str
    message: str


@dataclass(slots=True)
class SamplingFractionReport:
    tree_path: Path
    metadata_path: Path
    taxon_column: str
    sampling_column: str | None
    complete: bool
    matched_taxa: list[str]
    missing_taxa: list[str]
    invalid_rows: list[SamplingFractionIssue]
    sampling_fraction: float | None
    heterogeneous_values: bool
    warnings: list[str]


@dataclass(slots=True)
class DiversificationRateReport:
    tree_path: Path
    model: str
    crown_age: float
    observed_tip_count: int
    corrected_tip_count: float
    sampling_fraction: float
    birth_rate: float
    death_rate: float
    net_diversification_rate: float
    relative_extinction: float
    likelihood_kind: str
    log_likelihood: float
    aic: float
    assumptions: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DiversificationGammaStatisticReport:
    tree_path: Path
    tip_count: int
    rooted: bool
    ultrametric: bool
    bifurcating: bool
    root_age: float
    branching_time_count: int
    interval_count: int
    minimum_branching_time: float
    maximum_branching_time: float
    gamma_statistic: float
    sampling_fraction: float | None
    assumptions: list[str]
    warnings: list[str]


@dataclass(slots=True)
class DiversificationModelComparisonRow:
    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    sampling_fraction: float
    net_diversification_rate: float
    relative_extinction: float


@dataclass(slots=True)
class DiversificationModelComparisonReport:
    tree_path: Path
    metadata_path: Path | None
    better_model: str
    rows: list[DiversificationModelComparisonRow]


@dataclass(slots=True)
class CladeDiversificationObservation:
    node: str
    node_name: str | None
    descendant_taxa: list[str]
    tip_count: int
    crown_age: float
    diversification_rate: float
    z_score: float
    classification: str


@dataclass(slots=True)
class CladeDiversificationScanReport:
    tree_path: Path
    model: str
    global_rate: float
    observations: list[CladeDiversificationObservation]
    high_diversification_clades: list[CladeDiversificationObservation]
    low_diversification_clades: list[CladeDiversificationObservation]
    warnings: list[str]


@dataclass(slots=True)
class TraitDependentDiversificationState:
    state: str
    taxon_count: int
    taxa: list[str]
    monophyletic: bool
    crown_age: float | None
    diversification_rate: float | None
    warnings: list[str]


@dataclass(slots=True)
class TraitDependentDiversificationReport:
    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    observed_states: list[str]
    states: list[TraitDependentDiversificationState]
    warnings: list[str]


@dataclass(slots=True)
class DiversificationReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    tree_path: Path
    machine_manifest: dict[str, object]
    methods_summary_text: str
    methods_summary_warning_count: int
    methods_summary_path: Path | None
    report: DiversificationMethodReport


@dataclass(slots=True)
class DiversificationMethodReport:
    tree_path: Path
    metadata_path: Path | None
    traits_path: Path | None
    taxon_column: str | None
    sampling_column: str | None
    estimate_model: str
    clade_model: str
    clade_min_tip_count: int
    validation: TimeTreeValidationReport
    lineage: LineageThroughTimeReport
    gamma_statistic: DiversificationGammaStatisticReport
    primary_estimate: DiversificationRateReport
    model_comparison: DiversificationModelComparisonReport
    clade_scan: CladeDiversificationScanReport
    sampling_report: SamplingFractionReport | None
    trait_report: TraitDependentDiversificationReport | None


@dataclass(slots=True)
class DiversificationMethodsSummaryTextResult:
    output_path: Path
    title: str
    warning_count: int
    better_model: str
    sampling_metadata_complete: bool | None
    clade_observation_count: int
    text: str
    report: DiversificationMethodReport


def _node_depths(tree: PhyloTree) -> dict[str, float]:
    depths: dict[str, float] = {node_signature(tree.root): 0.0}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            next_depth = depth + float(child.branch_length or 0.0)
            depths[node_signature(child)] = next_depth
            if not child.is_leaf():
                visit(child, next_depth)

    visit(tree.root, 0.0)
    return depths


def _root_age(tree: PhyloTree) -> float:
    distances = [
        distance for _tip, distance in tree.root_to_tip_pairs() if distance is not None
    ]
    if not distances:
        raise DiversificationAnalysisError(
            "diversification analysis requires complete root-to-tip distances"
        )
    return float(format(max(distances), ".15g"))


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _sampling_fraction_from_rows(rows: list[float]) -> float:
    return float(format(sum(rows) / max(len(rows), 1), ".15g"))


def _bullet_list(values: list[str]) -> str:
    if not values:
        return "none"
    return ", ".join(f"`{value}`" for value in values)


def _deduplicate_text(values: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        normalized = value.strip()
        if not normalized or normalized in seen:
            continue
        seen.add(normalized)
        ordered.append(normalized)
    return ordered


def _is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(
        len(node.children) == 2 for node in tree.iter_nodes() if not node.is_leaf()
    )


def _resolve_sampling_column(columns: list[str], requested: str | None) -> str | None:
    if requested is not None:
        return requested if requested in columns else None
    return next((column for column in _SAMPLING_COLUMNS if column in columns), None)


def _lineages_at_time(
    report: LineageThroughTimeReport, time_before_present: float
) -> int:
    current = report.points[0].lineage_count
    for point in report.points[1:]:
        if point.time_before_present < time_before_present:
            break
        current = point.lineage_count
    return current


def _interval_log_likelihood(
    tree: PhyloTree, *, birth_rate: float, death_rate: float, relative_extinction: float
) -> float:
    root_age = _root_age(tree)
    depths = _node_depths(tree)
    event_rows = [
        (
            float(format(root_age - depths[node_signature(node)], ".15g")),
            max(len(node.children) - 1, 0),
        )
        for node in tree.iter_nodes()
        if node is not tree.root and not node.is_leaf()
    ]
    event_rows.sort(reverse=True)
    lineage_count = max(len(tree.root.children), 1)
    previous_age = root_age
    log_likelihood = 0.0
    effective_birth = max(
        birth_rate * max(1.0 - (relative_extinction / 2.0), 1e-9), 1e-9
    )
    turnover = max(birth_rate + death_rate, 1e-9)
    for event_age, gained_lineages in event_rows:
        interval = previous_age - event_age
        event_rate = max(lineage_count * effective_birth, 1e-9)
        log_likelihood += (gained_lineages * math.log(event_rate)) - (
            lineage_count * turnover * interval
        )
        lineage_count += gained_lineages
        previous_age = event_age
    log_likelihood -= lineage_count * turnover * previous_age
    return float(format(log_likelihood, ".15g"))


def _node_age(tree: PhyloTree, depths: dict[str, float], node: TreeNode) -> float:
    return float(format(_root_age(tree) - depths[node_signature(node)], ".15g"))


def _find_smallest_covering_node(
    tree: PhyloTree, taxa: set[str]
) -> tuple[TreeNode, list[str]]:
    best_node = tree.root
    best_taxa = _descendant_taxa(tree.root)
    for node in tree.iter_nodes():
        descendant_taxa = _descendant_taxa(node)
        descendant_set = set(descendant_taxa)
        if taxa <= descendant_set and len(descendant_set) <= len(best_taxa):
            best_node = node
            best_taxa = descendant_taxa
    return best_node, best_taxa


def validate_time_tree_for_diversification(tree_path: Path) -> TimeTreeValidationReport:
    """Validate the rooted ultrametric time-tree contract required for diversification analysis."""
    validation = validate_tree_path(tree_path, require_rooted=True)
    if validation.branch_length_status != "complete":
        raise DiversificationAnalysisError(
            "diversification analysis requires complete branch lengths"
        )
    ultrametric_report = assess_tree_ultrametricity(tree_path)
    if ultrametric_report.rooted is not True:
        raise UnrootedTreeError(f"tree is not rooted: {tree_path}")
    if not ultrametric_report.ultrametric:
        raise NonUltrametricTreeError(
            f"tree is not ultrametric within ape-style diversification tolerance: {tree_path}",
            details={
                "tolerance": ultrametric_report.tolerance,
                "criterion_name": ultrametric_report.criterion_name,
                "criterion_value": ultrametric_report.criterion_value,
                "max_tip_depth_deviation": ultrametric_report.max_tip_depth_deviation,
                "offending_taxa": list(ultrametric_report.offending_taxa),
            },
        )
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=validation.rooted,
        ultrametric=True,
        branch_length_status=validation.branch_length_status,
        tip_count=validation.tip_count,
        root_age=float(format(ultrametric_report.root_age, ".15g")),
        warnings=list(validation.warnings),
    )


def compute_lineage_through_time_curve(tree_path: Path) -> LineageThroughTimeReport:
    """Compute a deterministic lineage-through-time summary for one valid time tree."""
    validation = validate_time_tree_for_diversification(tree_path)
    tree = load_tree(tree_path)
    depths = _node_depths(tree)
    root_age = validation.root_age
    events = [
        (
            float(format(root_age - depths[node_signature(node)], ".15g")),
            node_signature(node),
            max(len(node.children) - 1, 0),
        )
        for node in tree.iter_nodes()
        if node is not tree.root and not node.is_leaf()
    ]
    events.sort(key=lambda item: (-item[0], item[1]))

    lineage_count = max(len(tree.root.children), 1)
    points = [
        LineageThroughTimePoint(
            node=node_signature(tree.root),
            time_before_present=root_age,
            lineage_count=lineage_count,
            event="root",
        )
    ]
    for event_age, signature, gained_lineages in events:
        lineage_count += gained_lineages
        points.append(
            LineageThroughTimePoint(
                node=signature,
                time_before_present=event_age,
                lineage_count=lineage_count,
                event="speciation",
            )
        )
    if points[-1].time_before_present != 0.0:
        points.append(
            LineageThroughTimePoint(
                node="present",
                time_before_present=0.0,
                lineage_count=tree.tip_count,
                event="present",
            )
        )
    return LineageThroughTimeReport(
        tree_path=tree_path,
        tip_count=tree.tip_count,
        root_age=root_age,
        points=points,
    )


def write_lineage_through_time_table(
    path: Path, report: LineageThroughTimeReport
) -> Path:
    """Export a lineage-through-time curve as a deterministic table."""
    rows = [
        {
            "node": point.node,
            "time_before_present": format(point.time_before_present, ".15g"),
            "lineage_count": str(point.lineage_count),
            "event": point.event,
        }
        for point in report.points
    ]
    return write_taxon_rows(
        path,
        columns=["node", "time_before_present", "lineage_count", "event"],
        rows=rows,
    )


def inspect_diversification_time_tree(tree_path: Path) -> TimeTreeValidationReport:
    """Inspect time-tree readiness with explicit diversification semantics."""
    inspection = inspect_tree_path(tree_path)
    if inspection.branch_length_status != "complete":
        raise DiversificationAnalysisError(
            "diversification analysis requires complete branch lengths"
        )
    if not inspection.rooted:
        raise DiversificationAnalysisError(
            "diversification analysis requires a rooted tree"
        )
    ultrametric_report = assess_tree_ultrametricity(tree_path)
    if not ultrametric_report.ultrametric:
        raise DiversificationAnalysisError(
            "diversification analysis requires an ultrametric time tree"
        )
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=inspection.rooted,
        ultrametric=True,
        branch_length_status=inspection.branch_length_status,
        tip_count=inspection.tip_count,
        root_age=float(format(ultrametric_report.root_age, ".15g")),
        warnings=list(inspection.warnings),
    )


def detect_incomplete_taxon_sampling_metadata(
    tree_path: Path,
    metadata_path: Path,
    *,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
) -> SamplingFractionReport:
    """Inspect taxon sampling fractions keyed to the tree tips."""
    tree = load_tree(tree_path)
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    resolved_sampling_column = _resolve_sampling_column(table.columns, sampling_column)
    if resolved_sampling_column is None:
        return SamplingFractionReport(
            tree_path=tree_path,
            metadata_path=metadata_path,
            taxon_column=table.taxon_column,
            sampling_column=None,
            complete=False,
            matched_taxa=[],
            missing_taxa=sorted(tree.tip_names),
            invalid_rows=[],
            sampling_fraction=None,
            heterogeneous_values=False,
            warnings=["metadata does not declare a sampling-fraction column"],
        )

    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    matched_taxa: list[str] = []
    missing_taxa: list[str] = []
    invalid_rows: list[SamplingFractionIssue] = []
    fractions: list[float] = []
    for taxon in sorted(tree.tip_names):
        row = rows_by_taxon.get(taxon)
        if row is None:
            missing_taxa.append(taxon)
            continue
        matched_taxa.append(taxon)
        raw_value = row[resolved_sampling_column].strip()
        if not raw_value:
            invalid_rows.append(
                SamplingFractionIssue(
                    taxon=taxon,
                    code="missing-sampling-fraction",
                    raw_value=raw_value,
                    message="sampling fraction is missing for this taxon",
                )
            )
            continue
        try:
            value = float(raw_value)
        except ValueError:
            invalid_rows.append(
                SamplingFractionIssue(
                    taxon=taxon,
                    code="invalid-sampling-fraction",
                    raw_value=raw_value,
                    message="sampling fraction must be numeric",
                )
            )
            continue
        if value <= 0.0 or value > 1.0:
            invalid_rows.append(
                SamplingFractionIssue(
                    taxon=taxon,
                    code="out-of-range-sampling-fraction",
                    raw_value=raw_value,
                    message="sampling fraction must be greater than 0 and at most 1",
                )
            )
            continue
        fractions.append(value)

    heterogeneous_values = len({format(value, ".9g") for value in fractions}) > 1
    warnings: list[str] = []
    if missing_taxa:
        warnings.append("sampling metadata does not cover every tree tip")
    if invalid_rows:
        warnings.append(
            "one or more sampling fractions are missing, invalid, or out of range"
        )
    if heterogeneous_values:
        warnings.append(
            "sampling fractions vary across taxa; mean sampling fraction will be used when correction is applied"
        )
    return SamplingFractionReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=table.taxon_column,
        sampling_column=resolved_sampling_column,
        complete=not missing_taxa and not invalid_rows,
        matched_taxa=matched_taxa,
        missing_taxa=missing_taxa,
        invalid_rows=invalid_rows,
        sampling_fraction=_sampling_fraction_from_rows(fractions)
        if fractions
        else None,
        heterogeneous_values=heterogeneous_values,
        warnings=warnings,
    )


def compute_diversification_gamma_statistic(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
) -> DiversificationGammaStatisticReport:
    """Compute the Pybus-Harvey diversification gamma statistic on one time tree."""
    validation = validate_time_tree_for_diversification(tree_path)
    tree = load_tree(tree_path)
    if tree.tip_count < 3:
        raise DiversificationAnalysisError(
            "diversification gamma statistic requires at least three tips",
            code="diversification_gamma_statistic_requires_three_or_more_tips",
            details={"tip_count": tree.tip_count},
        )
    if not _is_strictly_bifurcating(tree):
        raise DiversificationAnalysisError(
            "diversification gamma statistic requires a fully bifurcating tree",
            code="diversification_gamma_statistic_requires_bifurcating_tree",
            details={
                "tip_count": tree.tip_count,
                "internal_child_counts": [
                    len(node.children)
                    for node in tree.iter_nodes()
                    if not node.is_leaf()
                ],
            },
        )
    branching_time_report = compute_tree_branching_times(tree_path)
    branching_times = sorted(row.branching_time for row in branching_time_report.rows)
    if len(branching_times) != tree.tip_count - 1:
        raise DiversificationAnalysisError(
            "diversification gamma statistic requires one branching time per internal node on a fully bifurcating tree",
            code="diversification_gamma_statistic_requires_complete_branching_times",
            details={
                "tip_count": tree.tip_count,
                "branching_time_count": len(branching_times),
            },
        )
    intervals = [branching_times[0]]
    intervals.extend(
        branching_times[index] - branching_times[index - 1]
        for index in range(1, len(branching_times))
    )
    waiting_times = list(reversed(intervals))
    total_span = sum(
        multiplier * interval
        for multiplier, interval in zip(
            range(2, tree.tip_count + 1), waiting_times, strict=True
        )
    )
    if total_span <= 0.0:
        raise DiversificationAnalysisError(
            "diversification gamma statistic requires positive branching-time span",
            code="diversification_gamma_statistic_requires_positive_branching_time_span",
            details={"total_span": total_span},
        )
    running = 0.0
    cumulative_total = 0.0
    for multiplier, interval in zip(
        range(2, tree.tip_count), waiting_times[:-1], strict=True
    ):
        running += multiplier * interval
        cumulative_total += running
    statistic_mean = total_span / 2.0
    statistic_standard_deviation = total_span * math.sqrt(
        1.0 / (12.0 * (tree.tip_count - 2))
    )
    gamma_statistic = float(
        format(
            ((cumulative_total / (tree.tip_count - 2)) - statistic_mean)
            / statistic_standard_deviation,
            ".15g",
        )
    )

    sampling_fraction: float | None = None
    warnings = list(validation.warnings)
    if metadata_path is not None:
        sampling_report = detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
        sampling_fraction = sampling_report.sampling_fraction
        warnings.extend(sampling_report.warnings)
        if sampling_fraction is None:
            warnings.append(
                "gamma statistic assumes complete taxon sampling and the supplied metadata did not provide one valid complete-sampling estimate"
            )
        elif sampling_fraction < 1.0:
            warnings.append(
                "gamma statistic assumes complete taxon sampling and may be biased when the supplied sampling fractions indicate incomplete sampling"
            )

    assumptions = [
        "gamma statistic is interpreted on a rooted ultrametric crown tree",
        "gamma statistic is only computed for fully bifurcating trees with three or more tips",
    ]
    return DiversificationGammaStatisticReport(
        tree_path=tree_path,
        tip_count=tree.tip_count,
        rooted=True,
        ultrametric=True,
        bifurcating=True,
        root_age=validation.root_age,
        branching_time_count=len(branching_times),
        interval_count=len(waiting_times),
        minimum_branching_time=float(format(min(branching_times), ".15g")),
        maximum_branching_time=float(format(max(branching_times), ".15g")),
        gamma_statistic=gamma_statistic,
        sampling_fraction=sampling_fraction,
        assumptions=assumptions,
        warnings=sorted(set(warnings)),
    )


def estimate_diversification_rate(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    model: str = "birth-death",
) -> DiversificationRateReport:
    """Estimate one simple diversification model from a rooted ultrametric time tree."""
    if model not in {"yule", "birth-death"}:
        raise ValueError(f"unsupported diversification model: {model}")
    validation = validate_time_tree_for_diversification(tree_path)
    tree = load_tree(tree_path)
    ltt = compute_lineage_through_time_curve(tree_path)
    sampling_report = (
        detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
        if metadata_path is not None
        else None
    )
    sampling_fraction = (
        sampling_report.sampling_fraction
        if sampling_report and sampling_report.sampling_fraction
        else 1.0
    )
    corrected_tip_count = float(
        format(validation.tip_count / sampling_fraction, ".15g")
    )
    crown_age = validation.root_age
    net_diversification_rate = float(
        format(math.log(max(corrected_tip_count, 1.000000001)) / crown_age, ".15g")
    )
    midpoint_lineages = _lineages_at_time(ltt, crown_age / 2.0)
    expected_midpoint_lineages = math.sqrt(max(validation.tip_count, 1))
    slowdown = max(
        0.0, 1.0 - (midpoint_lineages / max(expected_midpoint_lineages, 1.0))
    )
    relative_extinction = (
        0.0 if model == "yule" else float(format(min(slowdown, 0.95), ".15g"))
    )
    birth_rate = (
        net_diversification_rate
        if model == "yule"
        else float(
            format(
                net_diversification_rate / max(1.0 - relative_extinction, 1e-9), ".15g"
            )
        )
    )
    death_rate = (
        0.0
        if model == "yule"
        else float(format(max(birth_rate - net_diversification_rate, 0.0), ".15g"))
    )
    log_likelihood = _interval_log_likelihood(
        tree,
        birth_rate=birth_rate,
        death_rate=death_rate,
        relative_extinction=relative_extinction,
    )
    parameter_count = 1 if model == "yule" else 2
    aic = float(format((2.0 * parameter_count) - (2.0 * log_likelihood), ".15g"))
    assumptions = [
        "tree is treated as a rooted ultrametric crown tree",
        "rates are estimated from extant tip counts and crown age only",
    ]
    warnings = list(validation.warnings)
    if sampling_report is not None:
        warnings.extend(sampling_report.warnings)
        if sampling_report.sampling_fraction is None:
            warnings.append(
                "sampling correction could not be applied because no valid sampling fractions were available"
            )
    if model == "birth-death":
        assumptions.append(
            "relative extinction is approximated from lineage-through-time slowdown"
        )
    else:
        assumptions.append("yule model assumes zero extinction")
    return DiversificationRateReport(
        tree_path=tree_path,
        model=model,
        crown_age=crown_age,
        observed_tip_count=validation.tip_count,
        corrected_tip_count=corrected_tip_count,
        sampling_fraction=sampling_fraction,
        birth_rate=birth_rate,
        death_rate=death_rate,
        net_diversification_rate=net_diversification_rate,
        relative_extinction=relative_extinction,
        likelihood_kind="heuristic-interval-log-likelihood",
        log_likelihood=log_likelihood,
        aic=aic,
        assumptions=assumptions,
        warnings=warnings,
    )


def compare_diversification_models(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
) -> DiversificationModelComparisonReport:
    """Compare Yule and simple birth-death diversification fits on one time tree."""
    yule = estimate_diversification_rate(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        model="yule",
    )
    birth_death = estimate_diversification_rate(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        model="birth-death",
    )
    rows = [
        DiversificationModelComparisonRow(
            model="yule",
            parameter_count=1,
            log_likelihood=yule.log_likelihood,
            aic=yule.aic,
            sampling_fraction=yule.sampling_fraction,
            net_diversification_rate=yule.net_diversification_rate,
            relative_extinction=yule.relative_extinction,
        ),
        DiversificationModelComparisonRow(
            model="birth-death",
            parameter_count=2,
            log_likelihood=birth_death.log_likelihood,
            aic=birth_death.aic,
            sampling_fraction=birth_death.sampling_fraction,
            net_diversification_rate=birth_death.net_diversification_rate,
            relative_extinction=birth_death.relative_extinction,
        ),
    ]
    better_model = min(rows, key=lambda row: row.aic).model
    return DiversificationModelComparisonReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        better_model=better_model,
        rows=rows,
    )


def detect_diversification_outlier_clades(
    tree_path: Path,
    *,
    min_tip_count: int = 2,
    model: str = "birth-death",
) -> CladeDiversificationScanReport:
    """Flag clades whose diversification rate is high or low relative to the tree-wide baseline."""
    global_report = estimate_diversification_rate(tree_path, model=model)
    tree = load_tree(tree_path)
    depths = _node_depths(tree)
    observations: list[CladeDiversificationObservation] = []
    raw_rows: list[tuple[TreeNode, list[str], float, float]] = []
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendant_taxa = _descendant_taxa(node)
        if len(descendant_taxa) < min_tip_count:
            continue
        crown_age = _node_age(tree, depths, node)
        if crown_age <= 0.0:
            continue
        diversification_rate = float(
            format(math.log(len(descendant_taxa)) / crown_age, ".15g")
        )
        raw_rows.append((node, descendant_taxa, crown_age, diversification_rate))
    rates = [row[3] for row in raw_rows]
    mean_rate = sum(rates) / max(len(rates), 1)
    variance = sum((rate - mean_rate) ** 2 for rate in rates) / max(len(rates), 1)
    standard_deviation = math.sqrt(variance)
    for node, descendant_taxa, crown_age, diversification_rate in raw_rows:
        z_score = (
            float(
                format((diversification_rate - mean_rate) / standard_deviation, ".15g")
            )
            if standard_deviation > 0.0
            else 0.0
        )
        if z_score >= 1.0:
            classification = "high"
        elif z_score <= -1.0:
            classification = "low"
        else:
            classification = "baseline"
        observations.append(
            CladeDiversificationObservation(
                node=node_signature(node),
                node_name=node.name,
                descendant_taxa=descendant_taxa,
                tip_count=len(descendant_taxa),
                crown_age=crown_age,
                diversification_rate=diversification_rate,
                z_score=z_score,
                classification=classification,
            )
        )
    high = [row for row in observations if row.classification == "high"]
    low = [row for row in observations if row.classification == "low"]
    return CladeDiversificationScanReport(
        tree_path=tree_path,
        model=model,
        global_rate=global_report.net_diversification_rate,
        observations=observations,
        high_diversification_clades=high,
        low_diversification_clades=low,
        warnings=list(global_report.warnings),
    )


def run_trait_dependent_diversification_analysis(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> TraitDependentDiversificationReport:
    """Summarize simple state-linked diversification rates when trait states form interpretable clades."""
    validate_time_tree_for_diversification(tree_path)
    tree = load_tree(tree_path)
    table = load_taxon_table(traits_path, taxon_column=taxon_column)
    if trait not in table.columns:
        raise DiversificationAnalysisError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    rows_by_taxon = {
        row[table.taxon_column]: row
        for row in table.rows
        if row[table.taxon_column] in tree_taxa and row[trait].strip()
    }
    observed_states = sorted({row[trait].strip() for row in rows_by_taxon.values()})
    depths = _node_depths(tree)
    states: list[TraitDependentDiversificationState] = []
    warnings: list[str] = []
    for state in observed_states:
        taxa = sorted(
            taxon for taxon, row in rows_by_taxon.items() if row[trait].strip() == state
        )
        state_warnings: list[str] = []
        if len(taxa) < 2:
            state_warnings.append("state is represented by fewer than two taxa")
            states.append(
                TraitDependentDiversificationState(
                    state=state,
                    taxon_count=len(taxa),
                    taxa=taxa,
                    monophyletic=False,
                    crown_age=None,
                    diversification_rate=None,
                    warnings=state_warnings,
                )
            )
            continue
        covering_node, descendant_taxa = _find_smallest_covering_node(tree, set(taxa))
        monophyletic = descendant_taxa == taxa
        crown_age = _node_age(tree, depths, covering_node)
        diversification_rate = (
            float(format(math.log(len(taxa)) / crown_age, ".15g"))
            if monophyletic and crown_age > 0.0
            else None
        )
        if not monophyletic:
            state_warnings.append("state taxa are not monophyletic in the input tree")
        states.append(
            TraitDependentDiversificationState(
                state=state,
                taxon_count=len(taxa),
                taxa=taxa,
                monophyletic=monophyletic,
                crown_age=crown_age if monophyletic else None,
                diversification_rate=diversification_rate,
                warnings=state_warnings,
            )
        )
        warnings.extend(state_warnings)
    return TraitDependentDiversificationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        observed_states=observed_states,
        states=states,
        warnings=sorted(set(warnings)),
    )


def build_diversification_method_report(
    tree_path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    traits_path: Path | None = None,
    trait: str | None = None,
    estimate_model: str = "birth-death",
    clade_model: str = "birth-death",
    clade_min_tip_count: int = 2,
) -> DiversificationMethodReport:
    """Build one integrated diversification-method report from governed evidence."""
    validation = validate_time_tree_for_diversification(tree_path)
    lineage = compute_lineage_through_time_curve(tree_path)
    gamma_statistic = compute_diversification_gamma_statistic(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
    )
    primary_estimate = estimate_diversification_rate(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        model=estimate_model,
    )
    model_comparison = compare_diversification_models(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
    )
    clade_scan = detect_diversification_outlier_clades(
        tree_path,
        min_tip_count=clade_min_tip_count,
        model=clade_model,
    )
    sampling_report = (
        detect_incomplete_taxon_sampling_metadata(
            tree_path,
            metadata_path,
            taxon_column=taxon_column,
            sampling_column=sampling_column,
        )
        if metadata_path is not None
        else None
    )
    trait_report = (
        run_trait_dependent_diversification_analysis(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
        )
        if traits_path is not None and trait is not None
        else None
    )
    return DiversificationMethodReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        estimate_model=estimate_model,
        clade_model=clade_model,
        clade_min_tip_count=clade_min_tip_count,
        validation=validation,
        lineage=lineage,
        gamma_statistic=gamma_statistic,
        primary_estimate=primary_estimate,
        model_comparison=model_comparison,
        clade_scan=clade_scan,
        sampling_report=sampling_report,
        trait_report=trait_report,
    )


def _diversification_methods_summary_assumptions(
    report: DiversificationMethodReport,
) -> list[str]:
    assumptions = [
        *report.gamma_statistic.assumptions,
        *report.primary_estimate.assumptions,
    ]
    if report.sampling_report is None:
        assumptions.append(
            "no sampling metadata was supplied, so all rate summaries default to the complete-sampling assumption"
        )
    elif report.sampling_report.heterogeneous_values:
        assumptions.append(
            "heterogeneous taxon sampling fractions are collapsed to the mean sampling fraction before correction is applied"
        )
    if report.trait_report is not None:
        assumptions.append(
            "trait-linked diversification summaries remain descriptive and only convert state-specific crown ages into rates when the state taxa are monophyletic"
        )
    return _deduplicate_text(assumptions)


def _diversification_methods_summary_warnings(
    report: DiversificationMethodReport,
) -> list[str]:
    warnings = [
        *report.validation.warnings,
        *report.gamma_statistic.warnings,
        *report.primary_estimate.warnings,
        *report.clade_scan.warnings,
    ]
    if report.sampling_report is not None:
        warnings.extend(report.sampling_report.warnings)
    if report.trait_report is not None:
        warnings.extend(report.trait_report.warnings)
    return _deduplicate_text(warnings)


def build_diversification_methods_summary_text(
    report: DiversificationMethodReport,
) -> str:
    """Build reviewer-facing Markdown methods text for one diversification analysis."""
    comparison_rows = sorted(report.model_comparison.rows, key=lambda row: row.aic)
    better_row = comparison_rows[0]
    runner_up_delta = (
        comparison_rows[1].aic - comparison_rows[0].aic
        if len(comparison_rows) > 1
        else 0.0
    )
    sampling_report = report.sampling_report
    trait_report = report.trait_report
    warnings = _diversification_methods_summary_warnings(report)
    assumptions = _diversification_methods_summary_assumptions(report)
    sampling_fraction_text = (
        "not available"
        if sampling_report is None or sampling_report.sampling_fraction is None
        else format(sampling_report.sampling_fraction, ".15g")
    )
    trait_state_count = 0 if trait_report is None else len(trait_report.states)
    trait_monophyletic_count = (
        0 if trait_report is None else sum(1 for row in trait_report.states if row.monophyletic)
    )
    return (
        "# Diversification Analysis Methods Summary\n\n"
        f"This diversification analysis evaluated rooted ultrametric time tree `{report.tree_path.name}`"
        + (
            f" with sampling metadata `{report.metadata_path.name}`"
            if report.metadata_path is not None
            else " without an external sampling-metadata table"
        )
        + (
            f" and trait table `{report.traits_path.name}` for state column `{trait_report.trait}`."
            if trait_report is not None
            else "."
        )
        + "\n\n## Time Tree And Inputs\n\n"
        f"- rooted time tree required: `{'yes' if report.validation.rooted else 'no'}`\n"
        f"- ultrametric time tree required: `{'yes' if report.validation.ultrametric else 'no'}`\n"
        f"- branch-length completeness: `{report.validation.branch_length_status}`\n"
        f"- analyzed tip count: `{report.validation.tip_count}`\n"
        f"- crown age: `{format(report.validation.root_age, '.15g')}`\n"
        f"- lineage-through-time points retained: `{len(report.lineage.points)}`\n"
        f"- tree validation warnings: {_bullet_list(report.validation.warnings)}\n\n"
        "## Sampling Correction\n\n"
        + (
            "- sampling metadata: not supplied, so diversification rates follow the complete-sampling assumption\n"
            if sampling_report is None
            else (
                f"- taxon column: `{sampling_report.taxon_column}`\n"
                + f"- sampling column: `{sampling_report.sampling_column or 'missing'}`\n"
                + f"- matched taxa with sampling rows: `{len(sampling_report.matched_taxa)}`\n"
                + f"- tree tips missing sampling rows: `{len(sampling_report.missing_taxa)}`\n"
                + f"- invalid sampling rows: `{len(sampling_report.invalid_rows)}`\n"
                + f"- sampling metadata complete: `{'yes' if sampling_report.complete else 'no'}`\n"
                + f"- mean sampling fraction used for correction: `{sampling_fraction_text}`\n"
                + f"- heterogeneous sampling fractions: `{'yes' if sampling_report.heterogeneous_values else 'no'}`\n"
                + f"- sampling warnings: {_bullet_list(sampling_report.warnings)}\n"
            )
        )
        + "\n## Diversification Models And Rates\n\n"
        f"- primary reported rate model: `{report.estimate_model}`\n"
        f"- compared candidate models: {_bullet_list([row.model for row in comparison_rows])}\n"
        f"- better-supported model by AIC: `{report.model_comparison.better_model}`\n"
        f"- better-model AIC: `{format(better_row.aic, '.15g')}`\n"
        f"- runner-up delta AIC: `{format(runner_up_delta, '.15g')}`\n"
        f"- reported net diversification rate: `{format(report.primary_estimate.net_diversification_rate, '.15g')}`\n"
        f"- reported relative extinction: `{format(report.primary_estimate.relative_extinction, '.15g')}`\n"
        f"- corrected tip count under the reported model: `{format(report.primary_estimate.corrected_tip_count, '.15g')}`\n"
        f"- Pybus-Harvey gamma statistic: `{format(report.gamma_statistic.gamma_statistic, '.15g')}`\n\n"
        "## Clade And Trait Review\n\n"
        f"- clade outlier scan model: `{report.clade_model}`\n"
        f"- minimum clade size for outlier review: `{report.clade_min_tip_count}`\n"
        f"- evaluated clades: `{len(report.clade_scan.observations)}`\n"
        f"- high-diversification outliers: `{len(report.clade_scan.high_diversification_clades)}`\n"
        f"- low-diversification outliers: `{len(report.clade_scan.low_diversification_clades)}`\n"
        + (
            "- trait-linked diversification surface: not requested\n"
            if trait_report is None
            else (
                f"- trait-linked diversification trait: `{trait_report.trait}`\n"
                + f"- observed trait states reviewed: `{trait_state_count}`\n"
                + f"- monophyletic states with interpretable crown ages: `{trait_monophyletic_count}`\n"
                + f"- trait-linked warnings: {_bullet_list(trait_report.warnings)}\n"
            )
        )
        + "\n## Assumptions And Caveats\n\n"
        + "\n".join(f"- {item}" for item in assumptions)
        + "\n\n## Reviewer Warnings\n\n"
        + f"- combined warning count: `{len(warnings)}`\n"
        + f"- warning details: {_bullet_list(warnings)}\n"
    )


def write_diversification_methods_summary_text(
    path: Path,
    report: DiversificationMethodReport,
) -> DiversificationMethodsSummaryTextResult:
    """Write reviewer-facing Markdown methods text for one diversification analysis."""
    text = build_diversification_methods_summary_text(report)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return DiversificationMethodsSummaryTextResult(
        output_path=path,
        title="Diversification Analysis Methods Summary",
        warning_count=len(_diversification_methods_summary_warnings(report)),
        better_model=report.model_comparison.better_model,
        sampling_metadata_complete=(
            None if report.sampling_report is None else report.sampling_report.complete
        ),
        clade_observation_count=len(report.clade_scan.observations),
        text=text,
        report=report,
    )


def write_clade_diversification_table(
    path: Path, report: CladeDiversificationScanReport
) -> Path:
    """Export clade diversification summaries as a deterministic TSV."""
    rows = [
        {
            "node": row.node,
            "node_name": row.node_name or "",
            "descendant_taxa": ",".join(row.descendant_taxa),
            "tip_count": str(row.tip_count),
            "crown_age": format(row.crown_age, ".15g"),
            "diversification_rate": format(row.diversification_rate, ".15g"),
            "z_score": format(row.z_score, ".15g"),
            "classification": row.classification,
        }
        for row in report.observations
    ]
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "tip_count",
            "crown_age",
            "diversification_rate",
            "z_score",
            "classification",
        ],
        rows=rows,
    )


def write_diversification_gamma_statistic_table(
    path: Path, report: DiversificationGammaStatisticReport
) -> Path:
    """Export one deterministic diversification gamma-statistic ledger."""
    rows = [
        {
            "tip_count": str(report.tip_count),
            "rooted": str(report.rooted).lower(),
            "ultrametric": str(report.ultrametric).lower(),
            "bifurcating": str(report.bifurcating).lower(),
            "root_age": format(report.root_age, ".15g"),
            "branching_time_count": str(report.branching_time_count),
            "interval_count": str(report.interval_count),
            "minimum_branching_time": format(report.minimum_branching_time, ".15g"),
            "maximum_branching_time": format(report.maximum_branching_time, ".15g"),
            "gamma_statistic": format(report.gamma_statistic, ".15g"),
            "sampling_fraction": ""
            if report.sampling_fraction is None
            else format(report.sampling_fraction, ".15g"),
            "assumptions": "; ".join(report.assumptions),
            "warnings": "; ".join(report.warnings),
        }
    ]
    return write_taxon_rows(
        path,
        columns=[
            "tip_count",
            "rooted",
            "ultrametric",
            "bifurcating",
            "root_age",
            "branching_time_count",
            "interval_count",
            "minimum_branching_time",
            "maximum_branching_time",
            "gamma_statistic",
            "sampling_fraction",
            "assumptions",
            "warnings",
        ],
        rows=rows,
    )


def write_diversification_model_comparison_table(
    path: Path, report: DiversificationModelComparisonReport
) -> Path:
    """Export one deterministic diversification model-comparison ledger."""
    rows = [
        {
            "model": row.model,
            "parameter_count": str(row.parameter_count),
            "log_likelihood": format(row.log_likelihood, ".15g"),
            "aic": format(row.aic, ".15g"),
            "sampling_fraction": format(row.sampling_fraction, ".15g"),
            "net_diversification_rate": format(
                row.net_diversification_rate, ".15g"
            ),
            "relative_extinction": format(row.relative_extinction, ".15g"),
            "better_model": str(row.model == report.better_model).lower(),
        }
        for row in report.rows
    ]
    return write_taxon_rows(
        path,
        columns=[
            "model",
            "parameter_count",
            "log_likelihood",
            "aic",
            "sampling_fraction",
            "net_diversification_rate",
            "relative_extinction",
            "better_model",
        ],
        rows=rows,
    )


def write_trait_dependent_diversification_table(
    path: Path, report: TraitDependentDiversificationReport
) -> Path:
    """Export state-linked diversification summaries as a deterministic TSV."""
    rows = [
        {
            "state": row.state,
            "taxon_count": str(row.taxon_count),
            "taxa": ",".join(row.taxa),
            "monophyletic": str(row.monophyletic).lower(),
            "crown_age": "" if row.crown_age is None else format(row.crown_age, ".15g"),
            "diversification_rate": ""
            if row.diversification_rate is None
            else format(row.diversification_rate, ".15g"),
            "warnings": "; ".join(row.warnings),
        }
        for row in report.states
    ]
    return write_taxon_rows(
        path,
        columns=[
            "state",
            "taxon_count",
            "taxa",
            "monophyletic",
            "crown_age",
            "diversification_rate",
            "warnings",
        ],
        rows=rows,
    )


def render_diversification_report(
    *,
    tree_path: Path,
    out_path: Path,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    traits_path: Path | None = None,
    trait: str | None = None,
    methods_summary_path: Path | None = None,
) -> DiversificationReportBuildResult:
    """Render a deterministic HTML report for diversification and macroevolution summaries."""
    report = build_diversification_method_report(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        traits_path=traits_path,
        trait=trait,
        estimate_model="birth-death",
        clade_model="birth-death",
        clade_min_tip_count=2,
    )
    methods_summary_text = build_diversification_methods_summary_text(report)
    methods_summary = (
        None
        if methods_summary_path is None
        else write_diversification_methods_summary_text(methods_summary_path, report)
    )
    sections = [
        ("methods-summary-text", methods_summary_text),
        ("lineage-through-time", json.dumps(report.lineage, default=str, indent=2)),
        (
            "diversification-gamma-statistic",
            json.dumps(report.gamma_statistic, default=str, indent=2),
        ),
        (
            "diversification-estimate",
            json.dumps(report.primary_estimate, default=str, indent=2),
        ),
        (
            "diversification-model-comparison",
            json.dumps(report.model_comparison, default=str, indent=2),
        ),
        (
            "clade-diversification-scan",
            json.dumps(report.clade_scan, default=str, indent=2),
        ),
    ]
    if report.trait_report is not None:
        sections.append(
            (
                "trait-dependent-diversification",
                json.dumps(report.trait_report, default=str, indent=2),
            )
        )
    title = "Bijux Diversification Report"
    manifest = {
        "report_kind": "diversification",
        "tree_path": str(tree_path),
        "metadata_path": None if metadata_path is None else str(metadata_path),
        "traits_path": None if traits_path is None else str(traits_path),
        "trait": trait,
        "sections": [name for name, _value in sections],
        "outputs": {
            "methods_summary_path": None
            if methods_summary is None
            else str(methods_summary.output_path)
        },
        "metrics": {
            "methods_summary_warning_count": len(
                _diversification_methods_summary_warnings(report)
            ),
            "better_model": report.model_comparison.better_model,
        },
    }
    write_html_report(
        title=title,
        sections=sections,
        out_path=out_path,
        embedded_json=manifest,
    )
    return DiversificationReportBuildResult(
        output_path=out_path,
        report_kind="diversification",
        title=title,
        tree_path=tree_path,
        machine_manifest=manifest,
        methods_summary_text=methods_summary_text,
        methods_summary_warning_count=len(
            _diversification_methods_summary_warnings(report)
        ),
        methods_summary_path=None if methods_summary is None else methods_summary.output_path,
        report=report,
    )
