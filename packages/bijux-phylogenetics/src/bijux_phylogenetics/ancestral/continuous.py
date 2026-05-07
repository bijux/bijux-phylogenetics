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
)

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
            weight_sum = sum(1.0 / child_variance for _, child_variance in child_payloads)
            estimate = (
                sum((value / child_variance) for value, child_variance in child_payloads)
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

    visit(dataset.tree.root)
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
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        alpha=stable_value(alpha),
        taxon_count=len(dataset.taxa),
        analysis_tree_newick=dump_pruned_tree(dataset.tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        dropped_non_numeric_taxa=dataset.dropped_non_numeric_taxa,
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        estimates=ordered_estimates,
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
