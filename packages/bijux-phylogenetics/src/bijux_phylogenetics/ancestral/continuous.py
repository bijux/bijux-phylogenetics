from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralContinuousDataset,
    dump_pruned_tree,
    load_continuous_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
    write_ancestral_rows,
)
from bijux_phylogenetics.core.tree import PhyloTree

_NORMAL_95_CRITICAL = 1.959963984540054


@dataclass(slots=True)
class ContinuousAncestralEstimate:
    """One continuous ancestral-state estimate for a tree node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    estimate: float
    standard_error: float
    lower_95_interval: float
    upper_95_interval: float
    uncertainty_width: float
    confidence: float
    interpretation: str
    unstable: bool
    downstream_risks: list[str]


@dataclass(slots=True)
class ContinuousAncestralReport:
    """Continuous ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    alpha: float
    taxon_count: int
    analysis_tree_newick: str
    dropped_missing_taxa: list[str]
    dropped_non_numeric_taxa: list[str]
    warnings: list[str]
    unstable_nodes: list[str]
    weak_support_nodes: list[str]
    estimates: list[ContinuousAncestralEstimate]


@dataclass(slots=True)
class ContinuousAncestralSummary:
    """Reviewer-facing summary for one continuous ancestral reconstruction."""

    trait: str
    taxon_column: str
    model: str
    alpha: float
    analyzed_taxon_count: int
    excluded_taxon_count: int
    missing_tip_taxon_count: int
    non_numeric_tip_taxon_count: int
    internal_node_count: int
    unstable_node_count: int
    weak_support_node_count: int
    root_node: str
    root_estimate: float
    root_standard_error: float
    root_lower_95_interval: float
    root_upper_95_interval: float
    warning_count: int


@dataclass(slots=True)
class ContinuousAncestralExclusion:
    """One excluded tip from a continuous ancestral reconstruction."""

    taxon: str
    reason: str


def reconstruct_continuous_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
) -> ContinuousAncestralReport:
    """Reconstruct continuous ancestral states under a Brownian or OU-style model."""
    if model not in {"brownian", "ou"}:
        raise ValueError(f"unsupported continuous ancestral model: {model}")
    if alpha <= 0:
        raise ValueError(
            f"alpha must be positive for continuous ancestral reconstruction, got {alpha}"
        )
    dataset = load_continuous_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return _reconstruct_continuous_from_dataset(
        dataset,
        working_tree=dataset.tree,
        model=model,
        alpha=alpha,
    )


def _reconstruct_continuous_from_dataset(
    dataset: AncestralContinuousDataset,
    *,
    working_tree: PhyloTree,
    model: str,
    alpha: float,
) -> ContinuousAncestralReport:
    global_mean = sum(dataset.values_by_taxon[taxon] for taxon in dataset.taxa) / len(
        dataset.taxa
    )
    sigma = _sample_standard_deviation(
        [dataset.values_by_taxon[taxon] for taxon in dataset.taxa]
    )
    trait_range = (
        max(dataset.values_by_taxon.values()) - min(dataset.values_by_taxon.values())
        if dataset.values_by_taxon
        else 0.0
    )
    estimates: list[ContinuousAncestralEstimate] = []

    def visit(node) -> tuple[float, float]:
        if node.is_leaf():
            estimate = dataset.values_by_taxon[node.name]
            estimates.append(
                ContinuousAncestralEstimate(
                    node=node_signature(node),
                    node_name=node.name,
                    is_tip=True,
                    descendant_taxa=node_descendant_taxa(node),
                    estimate=stable_value(estimate),
                    standard_error=0.0,
                    lower_95_interval=stable_value(estimate),
                    upper_95_interval=stable_value(estimate),
                    uncertainty_width=0.0,
                    confidence=1.0,
                    interpretation="observed tip value",
                    unstable=False,
                    downstream_risks=[],
                )
            )
            return estimate, float(node.branch_length or 0.0)

        if len(node.children) != 2:
            raise ValueError(
                "continuous ancestral reconstruction requires a fully dichotomous rooted tree"
            )

        left_child, right_child = node.children
        left_estimate, left_working_length = visit(left_child)
        right_estimate, right_working_length = visit(right_child)

        if model == "brownian":
            sum_working_lengths = max(
                left_working_length + right_working_length,
                1e-12,
            )
            estimate = (
                left_estimate * right_working_length
                + right_estimate * left_working_length
            ) / sum_working_lengths
            variance = sum_working_lengths
            returned_length = (
                float(node.branch_length or 0.0)
                + (left_working_length * right_working_length) / sum_working_lengths
            )
        else:
            child_payloads: list[tuple[float, float]] = []
            for child, child_estimate, child_variance in (
                (left_child, left_estimate, left_working_length),
                (right_child, right_estimate, right_working_length),
            ):
                branch_length = float(child.branch_length or 0.0)
                shrink = math.exp(-alpha * branch_length)
                transformed_estimate = (
                    shrink * child_estimate + (1.0 - shrink) * global_mean
                )
                stationary_variance = ((sigma**2) / (2.0 * alpha)) * (
                    1.0 - math.exp(-2.0 * alpha * branch_length)
                )
                propagated_variance = (
                    child_variance * math.exp(-2.0 * alpha * branch_length)
                ) + stationary_variance
                child_payloads.append(
                    (transformed_estimate, max(propagated_variance, 1e-12))
                )
            weight_sum = sum(
                1.0 / child_variance for _, child_variance in child_payloads
            )
            estimate = (
                sum(
                    (value / child_variance) for value, child_variance in child_payloads
                )
                / weight_sum
            )
            variance = 1.0 / weight_sum
            returned_length = variance

        standard_error = math.sqrt(max(variance, 0.0))
        lower = estimate - _NORMAL_95_CRITICAL * standard_error
        upper = estimate + _NORMAL_95_CRITICAL * standard_error
        uncertainty_width = max(0.0, upper - lower)
        confidence, unstable = _continuous_confidence(uncertainty_width, trait_range)
        interpretation = _continuous_interpretation(
            uncertainty_width, trait_range, unstable=unstable
        )
        estimates.append(
            ContinuousAncestralEstimate(
                node=node_signature(node),
                node_name=node.name,
                is_tip=False,
                descendant_taxa=node_descendant_taxa(node),
                estimate=stable_value(estimate),
                standard_error=stable_value(standard_error),
                lower_95_interval=stable_value(lower),
                upper_95_interval=stable_value(upper),
                uncertainty_width=stable_value(uncertainty_width),
                confidence=stable_value(confidence),
                interpretation=interpretation,
                unstable=unstable,
                downstream_risks=_continuous_downstream_risks(unstable),
            )
        )
        return estimate, returned_length

    visit(working_tree.root)
    ordered_estimates = _ordered_estimates(dataset, estimates)
    unstable_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.unstable
    ]
    weak_support_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.confidence < 0.75
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append(
            "one or more continuous ancestral estimates have broad uncertainty intervals"
        )
    if weak_support_nodes:
        warnings.append(
            "low-confidence ancestral estimates should not be overinterpreted for evolutionary timing or trait polarity"
        )
    return ContinuousAncestralReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=model,
        alpha=stable_value(alpha),
        taxon_count=len(dataset.taxa),
        analysis_tree_newick=dump_pruned_tree(working_tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        dropped_non_numeric_taxa=dataset.dropped_non_numeric_taxa,
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        estimates=ordered_estimates,
    )
def summarize_continuous_ancestral_report(
    report: ContinuousAncestralReport,
) -> ContinuousAncestralSummary:
    """Summarize the main review facts for one continuous ancestral report."""
    internal_estimates = [estimate for estimate in report.estimates if not estimate.is_tip]
    if not internal_estimates:
        raise ValueError(
            "continuous ancestral summary requires at least one internal-node estimate"
        )
    root_estimate = max(
        internal_estimates,
        key=lambda estimate: (
            len(estimate.descendant_taxa),
            estimate.node,
        ),
    )
    return ContinuousAncestralSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        alpha=report.alpha,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.dropped_missing_taxa)
        + len(report.dropped_non_numeric_taxa),
        missing_tip_taxon_count=len(report.dropped_missing_taxa),
        non_numeric_tip_taxon_count=len(report.dropped_non_numeric_taxa),
        internal_node_count=len(internal_estimates),
        unstable_node_count=len(report.unstable_nodes),
        weak_support_node_count=len(report.weak_support_nodes),
        root_node=root_estimate.node,
        root_estimate=root_estimate.estimate,
        root_standard_error=root_estimate.standard_error,
        root_lower_95_interval=root_estimate.lower_95_interval,
        root_upper_95_interval=root_estimate.upper_95_interval,
        warning_count=len(report.warnings),
    )


def continuous_ancestral_exclusions(
    report: ContinuousAncestralReport,
) -> list[ContinuousAncestralExclusion]:
    """Return one explicit exclusion row per dropped tip taxon."""
    rows = [
        ContinuousAncestralExclusion(
            taxon=taxon,
            reason="missing_trait_value",
        )
        for taxon in report.dropped_missing_taxa
    ]
    rows.extend(
        ContinuousAncestralExclusion(
            taxon=taxon,
            reason="non_numeric_trait_value",
        )
        for taxon in report.dropped_non_numeric_taxa
    )
    return rows


def write_continuous_ancestral_summary_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one summary ledger for a continuous ancestral reconstruction."""
    summary = summarize_continuous_ancestral_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "alpha",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "missing_tip_taxon_count",
            "non_numeric_tip_taxon_count",
            "internal_node_count",
            "unstable_node_count",
            "weak_support_node_count",
            "root_node",
            "root_estimate",
            "root_standard_error",
            "root_lower_95_interval",
            "root_upper_95_interval",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "alpha": str(summary.alpha),
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "missing_tip_taxon_count": str(summary.missing_tip_taxon_count),
                "non_numeric_tip_taxon_count": str(
                    summary.non_numeric_tip_taxon_count
                ),
                "internal_node_count": str(summary.internal_node_count),
                "unstable_node_count": str(summary.unstable_node_count),
                "weak_support_node_count": str(summary.weak_support_node_count),
                "root_node": summary.root_node,
                "root_estimate": str(summary.root_estimate),
                "root_standard_error": str(summary.root_standard_error),
                "root_lower_95_interval": str(summary.root_lower_95_interval),
                "root_upper_95_interval": str(summary.root_upper_95_interval),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_continuous_ancestral_uncertainty_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one internal-node uncertainty ledger for a continuous reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "estimate",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "uncertainty_width",
            "confidence",
            "interpretation",
            "unstable",
        ],
        rows=[
            {
                "node": estimate.node,
                "node_name": estimate.node_name or "",
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "estimate": str(estimate.estimate),
                "standard_error": str(estimate.standard_error),
                "lower_95_interval": str(estimate.lower_95_interval),
                "upper_95_interval": str(estimate.upper_95_interval),
                "uncertainty_width": str(estimate.uncertainty_width),
                "confidence": str(estimate.confidence),
                "interpretation": estimate.interpretation,
                "unstable": str(estimate.unstable).lower(),
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def write_continuous_ancestral_exclusion_table(
    path: Path, report: ContinuousAncestralReport
) -> Path:
    """Write one explicit excluded-tip ledger for a continuous reconstruction."""
    exclusions = continuous_ancestral_exclusions(report)
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in exclusions
        ],
    )


def _sample_standard_deviation(values: list[float]) -> float:
    if len(values) < 2:
        return 1.0
    mean_value = sum(values) / len(values)
    variance = sum((value - mean_value) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(max(variance, 1e-12))


def _ordered_estimates(
    dataset: AncestralContinuousDataset,
    estimates: list[ContinuousAncestralEstimate],
) -> list[ContinuousAncestralEstimate]:
    node_order = {
        node_signature(node): index
        for index, node in enumerate(dataset.tree.iter_nodes())
    }
    return sorted(estimates, key=lambda estimate: node_order[estimate.node])


def _continuous_confidence(
    uncertainty_width: float, trait_range: float
) -> tuple[float, bool]:
    if uncertainty_width == 0.0:
        return 1.0, False
    scale = max(trait_range, 1e-12)
    relative_width = uncertainty_width / scale
    confidence = max(0.0, min(1.0, 1.0 - min(relative_width, 1.0)))
    return stable_value(confidence), relative_width > 0.6 or confidence < 0.55


def _continuous_interpretation(
    uncertainty_width: float, trait_range: float, *, unstable: bool
) -> str:
    if uncertainty_width == 0.0:
        return "observed tip value"
    scale = max(trait_range, 1e-12)
    relative_width = uncertainty_width / scale
    if unstable:
        return "unstable node estimate"
    if relative_width <= 0.25:
        return "narrow uncertainty"
    if relative_width <= 0.6:
        return "moderate uncertainty"
    return "broad uncertainty"


def _continuous_downstream_risks(unstable: bool) -> list[str]:
    if not unstable:
        return []
    return [
        "node ordering and trait-polarity interpretations may change across alternative trees or models",
        "publication claims about deep ancestral values should be treated as provisional",
    ]
