from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    DuplicateTaxonError,
    PhylogeneticsError,
    UnrootedTreeError,
)
from bijux_phylogenetics.trees.tree_sets.topology import _rooted_topology_id

TREE_TOPOLOGY_PRIOR_FAMILIES = ("uniform-rooted-labeled-bifurcating",)


@dataclass(frozen=True, slots=True)
class TreeTopologyPriorModel:
    """One validated topology prior over a fixed labeled taxon set."""

    family: str
    taxa: list[str]
    topology_count: int
    log_topology_probability: float


@dataclass(frozen=True, slots=True)
class TreeTopologyPriorEvaluationReport:
    """One fixed-tree topology-prior evaluation report."""

    family: str
    taxa: list[str]
    topology_id: str
    topology_count: int
    tree_newick: str
    log_prior: float


def build_uniform_rooted_tree_topology_prior(
    taxa: list[str],
) -> TreeTopologyPriorModel:
    """Build the uniform prior over rooted labeled bifurcating trees."""
    validated_taxa = validate_tree_topology_prior_taxa(taxa)
    topology_count = count_rooted_labeled_bifurcating_topologies(len(validated_taxa))
    return TreeTopologyPriorModel(
        family="uniform-rooted-labeled-bifurcating",
        taxa=validated_taxa,
        topology_count=topology_count,
        log_topology_probability=-_log_rooted_labeled_bifurcating_topology_count(
            len(validated_taxa)
        ),
    )


def validate_tree_topology_prior_taxa(taxa: list[str]) -> list[str]:
    """Require an explicit distinct non-blank taxon set for topology priors."""
    if len(taxa) < 2:
        raise PhylogeneticsError(
            "tree topology prior requires at least two taxa",
            code="tree_topology_prior_requires_two_or_more_taxa",
        )
    blank_taxa = sorted({taxon for taxon in taxa if not taxon.strip()})
    if blank_taxa:
        raise PhylogeneticsError(
            "tree topology prior does not allow blank taxon labels",
            code="tree_topology_prior_blank_taxon_label",
            details={"blank_taxa": blank_taxa},
        )
    duplicate_taxa = sorted({taxon for taxon in taxa if taxa.count(taxon) > 1})
    if duplicate_taxa:
        raise DuplicateTaxonError(
            "tree topology prior requires distinct taxa; duplicates: "
            + ", ".join(duplicate_taxa),
            code="tree_topology_prior_duplicate_taxa",
            details={"duplicate_taxa": duplicate_taxa},
        )
    return sorted(taxa)


def count_rooted_labeled_bifurcating_topologies(taxon_count: int) -> int:
    """Count rooted labeled bifurcating trees on a fixed taxon set."""
    if taxon_count < 2:
        raise PhylogeneticsError(
            "rooted labeled bifurcating tree counts require at least two taxa",
            code="rooted_labeled_bifurcating_topology_count_requires_two_or_more_taxa",
        )
    count = 1
    for odd_value in range(1, (2 * taxon_count) - 2, 2):
        count *= odd_value
    return count


def evaluate_tree_topology_log_prior(
    tree: PhyloTree,
    prior_model: TreeTopologyPriorModel,
) -> TreeTopologyPriorEvaluationReport:
    """Evaluate one rooted labeled bifurcating tree under a fixed-set topology prior."""
    _validate_tree_topology_prior_tree(tree, prior_model)
    topology_id = _rooted_topology_id(tree, set(prior_model.taxa))
    return TreeTopologyPriorEvaluationReport(
        family=prior_model.family,
        taxa=list(prior_model.taxa),
        topology_id=topology_id,
        topology_count=prior_model.topology_count,
        tree_newick=tree.to_newick(),
        log_prior=prior_model.log_topology_probability,
    )


def _validate_tree_topology_prior_tree(
    tree: PhyloTree,
    prior_model: TreeTopologyPriorModel,
) -> None:
    if tree.rooted is not True:
        raise UnrootedTreeError("tree topology prior requires a rooted tree")
    if sorted(tree.tip_names) != prior_model.taxa:
        raise PhylogeneticsError(
            "tree topology prior requires the exact prior taxon set",
            code="tree_topology_prior_taxa_mismatch",
            details={
                "expected_taxa": list(prior_model.taxa),
                "observed_taxa": sorted(tree.tip_names),
            },
        )
    if not _tree_is_strictly_bifurcating(tree):
        raise PhylogeneticsError(
            "tree topology prior requires a strictly bifurcating tree",
            code="tree_topology_prior_requires_strictly_bifurcating_tree",
        )


def _tree_is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(
        len(node.children) == 2 for node in tree.iter_nodes() if not node.is_leaf()
    )


def _log_rooted_labeled_bifurcating_topology_count(taxon_count: int) -> float:
    return math.fsum(
        math.log(odd_value) for odd_value in range(1, (2 * taxon_count) - 2, 2)
    )
