from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.branch_lengths.branching_times import (
    compute_tree_branching_times,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import DiversificationAnalysisError

from .lineage import compute_lineage_through_time_curve
from .models import (
    DiversificationGammaStatisticReport,
    DiversificationModelComparisonReport,
    DiversificationModelComparisonRow,
    DiversificationRateReport,
    LineageThroughTimeReport,
)
from .sampling import detect_incomplete_taxon_sampling_metadata
from .trees import (
    is_strictly_bifurcating,
    node_depths,
    root_age,
    validate_time_tree_for_diversification,
)


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
    tree_root_age = root_age(tree)
    depths = node_depths(tree)
    event_rows = [
        (
            float(format(tree_root_age - depths[node_signature(node)], ".15g")),
            max(len(node.children) - 1, 0),
        )
        for node in tree.iter_nodes()
        if node is not tree.root and not node.is_leaf()
    ]
    event_rows.sort(reverse=True)
    lineage_count = max(len(tree.root.children), 1)
    previous_age = tree_root_age
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


def geiger_birth_death_exclusion_reason() -> str:
    return (
        "geiger::bd.ms birth-death parity is explicitly excluded in this round "
        "because local geiger exposes bd.ms as a simulation-oriented surface with "
        "arguments `phy`, `time`, `n`, `missing`, `crown`, and `epsilon`, while "
        "bijux currently reports heuristic Yule and birth-death diversification "
        "summaries rather than one optimized geiger-matched diversification "
        "likelihood contract"
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
    if not is_strictly_bifurcating(tree):
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
    lineage_report = compute_lineage_through_time_curve(tree_path)
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
    midpoint_lineages = _lineages_at_time(lineage_report, crown_age / 2.0)
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
        geiger_birth_death_exclusion_reason(),
    ]
    warnings = list(validation.warnings)
    warnings.append(geiger_birth_death_exclusion_reason())
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
