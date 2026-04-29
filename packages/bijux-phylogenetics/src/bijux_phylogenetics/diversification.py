from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import inspect_tree_path, validate_tree_path
from bijux_phylogenetics.errors import DiversificationAnalysisError
from bijux_phylogenetics.io.trees import load_tree

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
    distances = [distance for _tip, distance in tree.root_to_tip_pairs() if distance is not None]
    if not distances:
        raise DiversificationAnalysisError("diversification analysis requires complete root-to-tip distances")
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


def _resolve_sampling_column(columns: list[str], requested: str | None) -> str | None:
    if requested is not None:
        return requested if requested in columns else None
    return next((column for column in _SAMPLING_COLUMNS if column in columns), None)


def _lineages_at_time(report: LineageThroughTimeReport, time_before_present: float) -> int:
    current = report.points[0].lineage_count
    for point in report.points[1:]:
        if point.time_before_present < time_before_present:
            break
        current = point.lineage_count
    return current


def _interval_log_likelihood(tree: PhyloTree, *, birth_rate: float, death_rate: float, relative_extinction: float) -> float:
    root_age = _root_age(tree)
    depths = _node_depths(tree)
    event_rows = [
        (float(format(root_age - depths[node_signature(node)], ".15g")), max(len(node.children) - 1, 0))
        for node in tree.iter_nodes()
        if node is not tree.root and not node.is_leaf()
    ]
    event_rows.sort(reverse=True)
    lineage_count = max(len(tree.root.children), 1)
    previous_age = root_age
    log_likelihood = 0.0
    effective_birth = max(birth_rate * max(1.0 - (relative_extinction / 2.0), 1e-9), 1e-9)
    turnover = max(birth_rate + death_rate, 1e-9)
    for event_age, gained_lineages in event_rows:
        interval = previous_age - event_age
        event_rate = max(lineage_count * effective_birth, 1e-9)
        log_likelihood += (gained_lineages * math.log(event_rate)) - (lineage_count * turnover * interval)
        lineage_count += gained_lineages
        previous_age = event_age
    log_likelihood -= lineage_count * turnover * previous_age
    return float(format(log_likelihood, ".15g"))


def _node_age(tree: PhyloTree, depths: dict[str, float], node: TreeNode) -> float:
    return float(format(_root_age(tree) - depths[node_signature(node)], ".15g"))


def validate_time_tree_for_diversification(tree_path: Path) -> TimeTreeValidationReport:
    """Validate the rooted ultrametric time-tree contract required for diversification analysis."""
    validation = validate_tree_path(
        tree_path,
        require_rooted=True,
        require_ultrametric=True,
    )
    if validation.branch_length_status != "complete":
        raise DiversificationAnalysisError("diversification analysis requires complete branch lengths")
    tree = load_tree(tree_path)
    root_age = _root_age(tree)
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=validation.rooted,
        ultrametric=validation.ultrametric is True,
        branch_length_status=validation.branch_length_status,
        tip_count=validation.tip_count,
        root_age=root_age,
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


def write_lineage_through_time_table(path: Path, report: LineageThroughTimeReport) -> Path:
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
        raise DiversificationAnalysisError("diversification analysis requires complete branch lengths")
    if not inspection.rooted:
        raise DiversificationAnalysisError("diversification analysis requires a rooted tree")
    if inspection.is_ultrametric is not True:
        raise DiversificationAnalysisError("diversification analysis requires an ultrametric time tree")
    tree = load_tree(tree_path)
    return TimeTreeValidationReport(
        tree_path=tree_path,
        rooted=inspection.rooted,
        ultrametric=inspection.is_ultrametric is True,
        branch_length_status=inspection.branch_length_status,
        tip_count=inspection.tip_count,
        root_age=_root_age(tree),
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
        warnings.append("one or more sampling fractions are missing, invalid, or out of range")
    if heterogeneous_values:
        warnings.append("sampling fractions vary across taxa; mean sampling fraction will be used when correction is applied")
    return SamplingFractionReport(
        tree_path=tree_path,
        metadata_path=metadata_path,
        taxon_column=table.taxon_column,
        sampling_column=resolved_sampling_column,
        complete=not missing_taxa and not invalid_rows,
        matched_taxa=matched_taxa,
        missing_taxa=missing_taxa,
        invalid_rows=invalid_rows,
        sampling_fraction=_sampling_fraction_from_rows(fractions) if fractions else None,
        heterogeneous_values=heterogeneous_values,
        warnings=warnings,
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
    sampling_fraction = sampling_report.sampling_fraction if sampling_report and sampling_report.sampling_fraction else 1.0
    corrected_tip_count = float(format(validation.tip_count / sampling_fraction, ".15g"))
    crown_age = validation.root_age
    net_diversification_rate = float(format(math.log(max(corrected_tip_count, 1.000000001)) / crown_age, ".15g"))
    midpoint_lineages = _lineages_at_time(ltt, crown_age / 2.0)
    expected_midpoint_lineages = math.sqrt(max(validation.tip_count, 1))
    slowdown = max(0.0, 1.0 - (midpoint_lineages / max(expected_midpoint_lineages, 1.0)))
    relative_extinction = 0.0 if model == "yule" else float(format(min(slowdown, 0.95), ".15g"))
    birth_rate = (
        net_diversification_rate
        if model == "yule"
        else float(format(net_diversification_rate / max(1.0 - relative_extinction, 1e-9), ".15g"))
    )
    death_rate = 0.0 if model == "yule" else float(format(max(birth_rate - net_diversification_rate, 0.0), ".15g"))
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
            warnings.append("sampling correction could not be applied because no valid sampling fractions were available")
    if model == "birth-death":
        assumptions.append("relative extinction is approximated from lineage-through-time slowdown")
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
        diversification_rate = float(format(math.log(len(descendant_taxa)) / crown_age, ".15g"))
        raw_rows.append((node, descendant_taxa, crown_age, diversification_rate))
    rates = [row[3] for row in raw_rows]
    mean_rate = sum(rates) / max(len(rates), 1)
    variance = sum((rate - mean_rate) ** 2 for rate in rates) / max(len(rates), 1)
    standard_deviation = math.sqrt(variance)
    for node, descendant_taxa, crown_age, diversification_rate in raw_rows:
        z_score = (
            float(format((diversification_rate - mean_rate) / standard_deviation, ".15g"))
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
