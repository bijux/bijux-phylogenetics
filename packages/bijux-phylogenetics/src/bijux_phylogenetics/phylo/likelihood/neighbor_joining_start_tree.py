from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    AmbiguityPolicy,
    DistanceModel,
    GapHandlingMode,
    build_distance_tree_from_genetic_distance_matrix,
    compute_pairwise_genetic_distance_matrix_from_records,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology import (
    root_unrooted_tree_deterministically,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import NucleotideLikelihoodNeighborJoiningStartTreeReport
from .topology_search import resolve_nucleotide_topology_search_records


def build_likelihood_neighbor_joining_tree(
    records: list[AlignmentRecord] | Path,
    *,
    distance_model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> tuple[PhyloTree, NucleotideLikelihoodNeighborJoiningStartTreeReport]:
    """Build one rooted likelihood starting tree from native alignment-derived NJ distances."""
    resolved_records, resolved_alignment_path = (
        resolve_nucleotide_topology_search_records(records)
    )
    distance_report = compute_pairwise_genetic_distance_matrix_from_records(
        resolved_records,
        model=distance_model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    unrooted_tree, _distance_tree_report = build_distance_tree_from_genetic_distance_matrix(
        distance_report,
        method="neighbor-joining",
    )
    rooted_tree = _root_neighbor_joining_start_tree(unrooted_tree)
    report = NucleotideLikelihoodNeighborJoiningStartTreeReport(
        algorithm="nucleotide-likelihood-neighbor-joining-start-tree",
        distance_method="neighbor-joining",
        distance_model=distance_report.model,
        gap_handling=distance_report.gap_handling,
        ambiguity_policy=distance_report.ambiguity_policy,
        alignment_path=None
        if resolved_alignment_path is None
        else str(resolved_alignment_path),
        taxon_count=len(distance_report.identifiers),
        site_count=distance_report.alignment_length,
        pair_count=sum(
            1
            for pair in distance_report.pairs
            if pair.left_identifier != pair.right_identifier
        ),
        distance_warnings=list(distance_report.warnings),
        rooting_policy="largest-root-clade",
        unrooted_tree_newick=dumps_newick(unrooted_tree),
        rooted_tree_newick=dumps_newick(rooted_tree),
    )
    return rooted_tree, report


def build_likelihood_neighbor_joining_tree_from_alignment(
    alignment_path: Path,
    *,
    distance_model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> tuple[PhyloTree, NucleotideLikelihoodNeighborJoiningStartTreeReport]:
    """Build one rooted likelihood NJ starting tree from one FASTA alignment path."""
    return build_likelihood_neighbor_joining_tree(
        alignment_path,
        distance_model=distance_model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )


def _root_neighbor_joining_start_tree(tree: PhyloTree) -> PhyloTree:
    working_tree = tree.copy().refresh()
    if len(working_tree.root.children) == 2:
        working_tree.rooted = True
        return working_tree
    return root_unrooted_tree_deterministically(working_tree)
