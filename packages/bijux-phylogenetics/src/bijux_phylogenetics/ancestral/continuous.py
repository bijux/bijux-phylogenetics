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
    interpretation: str


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
        raise ValueError(f"alpha must be positive for continuous ancestral reconstruction, got {alpha}")
    dataset = load_continuous_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    global_mean = sum(dataset.values_by_taxon[taxon] for taxon in dataset.taxa) / len(dataset.taxa)
    sigma = _sample_standard_deviation([dataset.values_by_taxon[taxon] for taxon in dataset.taxa])
    trait_range = max(dataset.values_by_taxon.values()) - min(dataset.values_by_taxon.values()) if dataset.values_by_taxon else 0.0
    estimates: list[ContinuousAncestralEstimate] = []

    def visit(node) -> tuple[float, float]:
        if node.is_leaf():
            estimate = dataset.values_by_taxon[node.name]
            standard_error = 0.0
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
                    interpretation="observed tip value",
                )
            )
            return estimate, standard_error**2

        child_payloads: list[tuple[float, float]] = []
        for child in node.children:
            child_estimate, child_variance = visit(child)
            branch_length = float(child.branch_length or 0.0)
            if model == "brownian":
                transformed_estimate = child_estimate
                propagated_variance = child_variance + branch_length
            else:
                shrink = math.exp(-alpha * branch_length)
                transformed_estimate = shrink * child_estimate + (1.0 - shrink) * global_mean
                stationary_variance = ((sigma**2) / (2.0 * alpha)) * (1.0 - math.exp(-2.0 * alpha * branch_length))
                propagated_variance = (child_variance * math.exp(-2.0 * alpha * branch_length)) + stationary_variance
            child_payloads.append((transformed_estimate, max(propagated_variance, 1e-12)))

        weight_sum = sum(1.0 / variance for _, variance in child_payloads)
        estimate = sum((value / variance) for value, variance in child_payloads) / weight_sum
        variance = 1.0 / weight_sum
        standard_error = math.sqrt(variance)
        lower = estimate - 1.96 * standard_error
        upper = estimate + 1.96 * standard_error
        uncertainty_width = max(0.0, upper - lower)
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
                interpretation=_continuous_interpretation(uncertainty_width, trait_range),
            )
        )
        return estimate, variance

    visit(dataset.tree.root)
    ordered_estimates = _ordered_estimates(dataset, estimates)
    unstable_nodes = [
        estimate.node
        for estimate in ordered_estimates
        if not estimate.is_tip and estimate.interpretation == "broad uncertainty"
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append("one or more continuous ancestral estimates have broad uncertainty intervals")
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


def _continuous_interpretation(uncertainty_width: float, trait_range: float) -> str:
    if uncertainty_width == 0.0:
        return "observed tip value"
    scale = max(trait_range, 1e-12)
    relative_width = uncertainty_width / scale
    if relative_width <= 0.25:
        return "narrow uncertainty"
    if relative_width <= 0.6:
        return "moderate uncertainty"
    return "broad uncertainty"
