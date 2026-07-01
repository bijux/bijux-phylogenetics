from __future__ import annotations

import json
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


def write_nucleotide_likelihood_neighbor_joining_start_tree_run_json(
    path: Path,
    report: NucleotideLikelihoodNeighborJoiningStartTreeReport,
) -> Path:
    """Write one machine-readable native NJ starting-tree payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "algorithm": report.algorithm,
        "distance_method": report.distance_method,
        "distance_model": report.distance_model,
        "gap_handling": report.gap_handling,
        "ambiguity_policy": report.ambiguity_policy,
        "alignment_path": report.alignment_path,
        "taxon_count": report.taxon_count,
        "site_count": report.site_count,
        "pair_count": report.pair_count,
        "distance_warnings": report.distance_warnings,
        "rooting_policy": report.rooting_policy,
        "unrooted_tree_newick": report.unrooted_tree_newick,
        "rooted_tree_newick": report.rooted_tree_newick,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path


def write_nucleotide_likelihood_neighbor_joining_start_tree_artifacts(
    out_dir: Path,
    report: NucleotideLikelihoodNeighborJoiningStartTreeReport,
) -> dict[str, Path]:
    """Write the governed artifact family for one native NJ starting tree."""
    out_dir.mkdir(parents=True, exist_ok=True)
    unrooted_tree_path = out_dir / "unrooted_tree.nwk"
    unrooted_tree_path.write_text(report.unrooted_tree_newick + "\n", encoding="utf-8")
    rooted_tree_path = out_dir / "rooted_tree.nwk"
    rooted_tree_path.write_text(report.rooted_tree_newick + "\n", encoding="utf-8")
    run_json_path = write_nucleotide_likelihood_neighbor_joining_start_tree_run_json(
        out_dir / "run.json",
        report,
    )
    return {
        "unrooted_tree_path": unrooted_tree_path,
        "rooted_tree_path": rooted_tree_path,
        "run_json_path": run_json_path,
    }


def _root_neighbor_joining_start_tree(tree: PhyloTree) -> PhyloTree:
    working_tree = tree.copy().refresh()
    if len(working_tree.root.children) == 2:
        working_tree.rooted = True
        return working_tree
    return root_unrooted_tree_deterministically(working_tree)
