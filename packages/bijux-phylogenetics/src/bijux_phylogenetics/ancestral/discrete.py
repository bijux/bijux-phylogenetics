from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    dump_pruned_tree,
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
)
from bijux_phylogenetics.discrete_evolution import run_discrete_state_transition_model


@dataclass(slots=True)
class DiscreteAncestralEstimate:
    """One discrete ancestral-state estimate for a tree node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    state_set: list[str]
    most_likely_state: str
    state_probabilities: dict[str, float]
    ambiguous: bool
    confidence: float
    interpretation: str
    unstable: bool
    downstream_risks: list[str]


@dataclass(slots=True)
class DiscreteAncestralReport:
    """Discrete ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    analysis_tree_newick: str
    dropped_missing_taxa: list[str]
    warnings: list[str]
    unstable_nodes: list[str]
    weak_support_nodes: list[str]
    estimates: list[DiscreteAncestralEstimate]


def reconstruct_discrete_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteAncestralReport:
    """Reconstruct discrete ancestral states under Fitch parsimony."""
    resolved_model = _resolve_discrete_model_name(model)
    if resolved_model == "fitch" and state_ordering != "unordered":
        raise ValueError(
            "ordered discrete ancestral reconstruction requires a likelihood model"
        )
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if resolved_model != "fitch":
        likelihood_report = run_discrete_state_transition_model(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=resolved_model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
        )
        estimates = [
            _build_discrete_estimate(
                node=estimate.node,
                node_name=estimate.node_name,
                is_tip=estimate.is_tip,
                descendant_taxa=estimate.descendant_taxa,
                most_likely_state=estimate.most_likely_state,
                state_probabilities=estimate.state_probabilities,
            )
            for estimate in likelihood_report.estimates
        ]
        unstable_nodes = [
            estimate.node
            for estimate in estimates
            if estimate.unstable and not estimate.is_tip
        ]
        weak_support_nodes = [
            estimate.node
            for estimate in estimates
            if not estimate.is_tip and estimate.confidence < 0.75
        ]
        warnings = list(dataset.warnings)
        if unstable_nodes:
            warnings.append(
                "one or more discrete ancestral nodes remain unstable across candidate states"
            )
        if weak_support_nodes:
            warnings.append(
                "low-confidence ancestral state assignments should not be overinterpreted as definitive transitions"
            )
        return DiscreteAncestralReport(
            tree_path=tree_path,
            traits_path=traits_path,
            taxon_column=dataset.taxon_column,
            trait=trait,
            model=resolved_model,
            state_ordering=state_ordering,
            ordered_states=likelihood_report.ordered_states,
            taxon_count=len(dataset.taxa),
            observed_states=dataset.observed_states,
            state_counts=dataset.state_counts,
            sparse_states=dataset.sparse_states,
            analysis_tree_newick=dump_pruned_tree(dataset.tree),
            dropped_missing_taxa=dataset.dropped_missing_taxa,
            warnings=warnings,
            unstable_nodes=unstable_nodes,
            weak_support_nodes=weak_support_nodes,
            estimates=estimates,
        )
    estimates: list[DiscreteAncestralEstimate] = []

    def downpass(node) -> set[str]:
        if node.is_leaf():
            state = dataset.states_by_taxon[node.name]
            return {state}
        child_sets = [downpass(child) for child in node.children]
        candidate = set(child_sets[0])
        for child_set in child_sets[1:]:
            intersection = candidate & child_set
            if intersection:
                candidate = intersection
            else:
                candidate |= child_set
        return candidate

    candidate_sets = {
        node_signature(node): downpass(node) for node in dataset.tree.iter_nodes()
    }

    for node in dataset.tree.iter_nodes():
        signature = node_signature(node)
        if node.is_leaf():
            resolved_state = dataset.states_by_taxon[node.name]
            probabilities = {resolved_state: 1.0}
            state_set = [resolved_state]
        else:
            state_set = sorted(candidate_sets[signature])
            probabilities = {state: 1.0 / len(state_set) for state in state_set}
            resolved_state = state_set[0]
        estimates.append(
            _build_discrete_estimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                most_likely_state=resolved_state,
                state_probabilities=probabilities,
            )
        )
    unstable_nodes = [
        estimate.node
        for estimate in estimates
        if estimate.unstable and not estimate.is_tip
    ]
    weak_support_nodes = [
        estimate.node
        for estimate in estimates
        if not estimate.is_tip and estimate.confidence < 0.75
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append(
            "one or more discrete ancestral nodes remain unstable across candidate states"
        )
    if weak_support_nodes:
        warnings.append(
            "low-confidence ancestral state assignments should not be overinterpreted as definitive transitions"
        )

    return DiscreteAncestralReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=resolved_model,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        taxon_count=len(dataset.taxa),
        observed_states=dataset.observed_states,
        state_counts=dataset.state_counts,
        sparse_states=dataset.sparse_states,
        analysis_tree_newick=dump_pruned_tree(dataset.tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        estimates=estimates,
    )


def _resolve_discrete_model_name(model: str) -> str:
    aliases = {
        "fitch": "fitch",
        "equal-rates": "equal-rates",
        "er": "equal-rates",
        "symmetric": "symmetric",
        "sym": "symmetric",
        "all-rates-different": "all-rates-different",
        "ard": "all-rates-different",
    }
    resolved = aliases.get(model)
    if resolved is None:
        raise ValueError(f"unsupported discrete ancestral model: {model}")
    return resolved


def _build_discrete_estimate(
    *,
    node: str,
    node_name: str | None,
    is_tip: bool,
    descendant_taxa: list[str],
    most_likely_state: str,
    state_probabilities: dict[str, float],
) -> DiscreteAncestralEstimate:
    state_set = sorted(state_probabilities)
    ordered_probabilities = sorted(state_probabilities.values(), reverse=True)
    confidence = ordered_probabilities[0] if ordered_probabilities else 0.0
    runner_up = ordered_probabilities[1] if len(ordered_probabilities) > 1 else 0.0
    unstable = not is_tip and ((confidence - runner_up) < 0.15 or confidence < 0.7)
    if is_tip:
        interpretation = "observed tip state"
    elif unstable:
        interpretation = "unstable node state"
    elif confidence >= 0.9:
        interpretation = "strongly supported node state"
    else:
        interpretation = "moderately supported node state"
    return DiscreteAncestralEstimate(
        node=node,
        node_name=node_name,
        is_tip=is_tip,
        descendant_taxa=descendant_taxa,
        state_set=state_set,
        most_likely_state=most_likely_state,
        state_probabilities=state_probabilities,
        ambiguous=len(state_set) > 1,
        confidence=confidence,
        interpretation=interpretation,
        unstable=unstable,
        downstream_risks=_discrete_downstream_risks(unstable),
    )


def _discrete_downstream_risks(unstable: bool) -> list[str]:
    if not unstable:
        return []
    return [
        "transition counts and inferred ancestral geography may change under alternative state models",
        "biological narratives about ancestral states should be treated as provisional for this node",
    ]
