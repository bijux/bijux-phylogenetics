from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field

import numpy

from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

FiniteStateLeafVector = Callable[[str, str], numpy.ndarray]
FiniteStateTransitionMatrix = Callable[[float], numpy.ndarray]


@dataclass(frozen=True, slots=True)
class FiniteStateLikelihoodSurface:
    """One reusable finite-state CTMC model surface over rooted-tree tips."""

    model_name: str
    state_labels: tuple[str, ...]
    root_prior: numpy.ndarray
    root_prior_source: str
    observation_policy: str
    transition_matrix_for_branch_length: FiniteStateTransitionMatrix
    leaf_likelihood_vector_for_state: FiniteStateLeafVector
    parameter_values: Mapping[str, float] = field(default_factory=dict)

    @property
    def state_count(self) -> int:
        return len(self.state_labels)


@dataclass(slots=True)
class CompiledFiniteStateLikelihoodContext:
    """One tree-bound finite-state likelihood context with cached branches."""

    tree: PhyloTree
    surface: FiniteStateLikelihoodSurface
    transition_matrix_by_node_id: dict[str, numpy.ndarray]


def compile_finite_state_likelihood_context(
    tree: PhyloTree,
    surface: FiniteStateLikelihoodSurface,
) -> CompiledFiniteStateLikelihoodContext:
    """Compile one rooted tree plus one finite-state model into cached branches."""
    transition_matrix_by_node_id: dict[str, numpy.ndarray] = {}
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        transition_matrix_by_node_id[child.node_id] = (
            surface.transition_matrix_for_branch_length(
                max(float(child.branch_length or 0.0), 0.0)
            )
        )
    return CompiledFiniteStateLikelihoodContext(
        tree=tree,
        surface=surface,
        transition_matrix_by_node_id=transition_matrix_by_node_id,
    )


def evaluate_finite_state_site_log_likelihood(
    context: CompiledFiniteStateLikelihoodContext,
    states: tuple[str, ...],
    *,
    taxon_order: list[str],
) -> float:
    """Evaluate one aligned site against one compiled finite-state context."""
    states_by_taxon = dict(zip(taxon_order, states, strict=True))
    pruning_pass = postorder_conditional_likelihoods(
        context.tree,
        state_count=context.surface.state_count,
        leaf_likelihood=lambda node: _leaf_likelihood_vector(
            context.surface,
            states_by_taxon=states_by_taxon,
            node_name=node.name,
        ),
        transition_matrix_for_child=lambda child: _transition_matrix_for_child(
            context,
            child_node_id=child.node_id,
        ),
    )
    return log_likelihood_from_root_prior(
        context.tree,
        pruning_pass,
        root_prior=context.surface.root_prior,
    )


def evaluate_finite_state_tree_log_likelihood_from_patterns(
    context: CompiledFiniteStateLikelihoodContext,
    compressed_patterns: CompressedAlignmentSitePatterns,
) -> float:
    """Evaluate one compressed alignment against one compiled finite-state model."""
    total = 0.0
    for pattern in compressed_patterns.patterns:
        total += pattern.weight * evaluate_finite_state_site_log_likelihood(
            context,
            pattern.states,
            taxon_order=compressed_patterns.taxon_order,
        )
    return total


def evaluate_finite_state_tree_log_likelihood(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    surface: FiniteStateLikelihoodSurface,
) -> float:
    """Evaluate one compressed alignment on one shared finite-state surface."""
    validate_explicit_branch_lengths(tree, model_name=surface.model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=surface.model_name,
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    return evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )


def _leaf_likelihood_vector(
    surface: FiniteStateLikelihoodSurface,
    *,
    states_by_taxon: dict[str, str],
    node_name: str | None,
) -> numpy.ndarray:
    if node_name is None:
        raise ValueError(
            f"{surface.model_name} likelihood requires named tree tips for alignment lookup"
        )
    return surface.leaf_likelihood_vector_for_state(
        node_name,
        states_by_taxon[node_name],
    )


def _transition_matrix_for_child(
    context: CompiledFiniteStateLikelihoodContext,
    *,
    child_node_id: str | None,
) -> numpy.ndarray:
    if child_node_id is None:
        raise ValueError("tree node is missing a stable node_id")
    return context.transition_matrix_by_node_id[child_node_id]
