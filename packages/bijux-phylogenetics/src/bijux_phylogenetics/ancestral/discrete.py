from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    dump_pruned_tree,
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
)


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


@dataclass(slots=True)
class DiscreteAncestralReport:
    """Discrete ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    analysis_tree_newick: str
    dropped_missing_taxa: list[str]
    warnings: list[str]
    estimates: list[DiscreteAncestralEstimate]


def reconstruct_discrete_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
) -> DiscreteAncestralReport:
    """Reconstruct discrete ancestral states under Fitch parsimony."""
    if model != "fitch":
        raise ValueError(f"unsupported discrete ancestral model: {model}")
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
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
        node_signature(node): downpass(node)
        for node in dataset.tree.iter_nodes()
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
            DiscreteAncestralEstimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                state_set=state_set,
                most_likely_state=resolved_state,
                state_probabilities=probabilities,
                ambiguous=len(state_set) > 1,
            )
        )

    return DiscreteAncestralReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        taxon_count=len(dataset.taxa),
        observed_states=dataset.observed_states,
        state_counts=dataset.state_counts,
        sparse_states=dataset.sparse_states,
        analysis_tree_newick=dump_pruned_tree(dataset.tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        warnings=list(dataset.warnings),
        estimates=estimates,
    )
