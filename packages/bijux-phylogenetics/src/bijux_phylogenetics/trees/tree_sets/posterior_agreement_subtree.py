from __future__ import annotations

import csv
from itertools import combinations
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_object_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .contracts import (
    PosteriorAgreementSubtreeCandidateRow,
    PosteriorAgreementSubtreeReport,
)
from .inventory import _analyze_tree_set, _require_exact_taxa
from .topology import _rooted_topology_id

POSTERIOR_AGREEMENT_SUBTREE_SEARCH_STRATEGY = "exact-descending-retained-subsets"


def _candidate_subset_count(shared_taxon_count: int) -> int:
    return sum(
        math.comb(shared_taxon_count, retained_taxon_count)
        for retained_taxon_count in range(2, shared_taxon_count + 1)
    )


def _iter_retained_taxon_subsets(shared_taxa: list[str]):
    ordered_taxa = sorted(shared_taxa)
    for retained_taxon_count in range(len(ordered_taxa), 1, -1):
        for retained_taxa in combinations(ordered_taxa, retained_taxon_count):
            yield list(retained_taxa)


def summarize_posterior_agreement_subtree(
    path: Path,
) -> tuple[PhyloTree, PosteriorAgreementSubtreeReport]:
    """Find the largest retained taxon subset whose pruned posterior trees share one rooted topology."""
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    candidate_rows: list[PosteriorAgreementSubtreeCandidateRow] = []
    agreement_tree: PhyloTree | None = None
    retained_taxa: list[str] | None = None
    stable_rooted_topology_id: str | None = None
    for candidate_index, retained_subset in enumerate(
        _iter_retained_taxon_subsets(exact_taxa),
        start=1,
    ):
        pruned_trees = [
            prune_tree_object_to_requested_taxa(tree, retained_subset)
            for tree in analysis.trees
        ]
        retained_taxa_set = set(retained_subset)
        rooted_topology_counts: dict[str, int] = {}
        rooted_representatives: dict[str, PhyloTree] = {}
        for tree in pruned_trees:
            topology_id = _rooted_topology_id(tree, retained_taxa_set)
            rooted_topology_counts[topology_id] = (
                rooted_topology_counts.get(topology_id, 0) + 1
            )
            rooted_representatives.setdefault(topology_id, tree)
        dominant_rooted_topology_id, dominant_rooted_topology_count = max(
            rooted_topology_counts.items(),
            key=lambda item: (item[1], item[0]),
        )
        stable_topology_reached = len(rooted_topology_counts) == 1
        candidate_rows.append(
            PosteriorAgreementSubtreeCandidateRow(
                candidate_index=candidate_index,
                retained_taxon_count=len(retained_subset),
                retained_taxa=retained_subset,
                removed_taxa=sorted(set(exact_taxa) - set(retained_subset)),
                rooted_topology_count=len(rooted_topology_counts),
                dominant_rooted_topology_frequency=round(
                    dominant_rooted_topology_count / len(pruned_trees),
                    15,
                ),
                stable_topology_reached=stable_topology_reached,
            )
        )
        if stable_topology_reached:
            retained_taxa = retained_subset
            stable_rooted_topology_id = dominant_rooted_topology_id
            agreement_tree = rooted_representatives[dominant_rooted_topology_id]
            break
    if (
        agreement_tree is None
        or retained_taxa is None
        or stable_rooted_topology_id is None
    ):
        raise AssertionError(
            "posterior agreement subtree summary must find a stable retained taxon subset"
        )
    report = PosteriorAgreementSubtreeReport(
        path=analysis.path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        search_strategy=POSTERIOR_AGREEMENT_SUBTREE_SEARCH_STRATEGY,
        possible_retained_subset_count=_candidate_subset_count(len(exact_taxa)),
        evaluated_candidate_count=len(candidate_rows),
        retained_taxa=retained_taxa,
        agreement_removed_taxa=sorted(set(exact_taxa) - set(retained_taxa)),
        stable_rooted_topology_id=stable_rooted_topology_id,
        agreement_subtree_newick=dumps_newick(agreement_tree),
        candidate_rows=candidate_rows,
    )
    return agreement_tree, report


def write_posterior_agreement_subtree_summary_table(
    path: Path,
    report: PosteriorAgreementSubtreeReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_count",
                "search_strategy",
                "possible_retained_subset_count",
                "evaluated_candidate_count",
                "retained_taxa",
                "agreement_removed_taxa",
                "stable_rooted_topology_id",
                "agreement_subtree_newick",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "tree_count": report.tree_count,
                "search_strategy": report.search_strategy,
                "possible_retained_subset_count": report.possible_retained_subset_count,
                "evaluated_candidate_count": report.evaluated_candidate_count,
                "retained_taxa": "|".join(report.retained_taxa),
                "agreement_removed_taxa": "|".join(report.agreement_removed_taxa),
                "stable_rooted_topology_id": report.stable_rooted_topology_id,
                "agreement_subtree_newick": report.agreement_subtree_newick,
            }
        )
    return path


def write_posterior_agreement_subtree_removed_taxa_table(
    path: Path,
    report: PosteriorAgreementSubtreeReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["taxon", "removed_for_agreement_subtree"],
            delimiter="\t",
        )
        writer.writeheader()
        for taxon in report.shared_taxa:
            writer.writerow(
                {
                    "taxon": taxon,
                    "removed_for_agreement_subtree": str(
                        taxon in set(report.agreement_removed_taxa)
                    ).lower(),
                }
            )
    return path


def write_posterior_agreement_subtree_search_table(
    path: Path,
    report: PosteriorAgreementSubtreeReport,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "candidate_index",
                "retained_taxon_count",
                "retained_taxa",
                "removed_taxa",
                "rooted_topology_count",
                "dominant_rooted_topology_frequency",
                "stable_topology_reached",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.candidate_rows:
            writer.writerow(
                {
                    "candidate_index": row.candidate_index,
                    "retained_taxon_count": row.retained_taxon_count,
                    "retained_taxa": "|".join(row.retained_taxa),
                    "removed_taxa": "|".join(row.removed_taxa),
                    "rooted_topology_count": row.rooted_topology_count,
                    "dominant_rooted_topology_frequency": format(
                        row.dominant_rooted_topology_frequency,
                        ".15g",
                    ),
                    "stable_topology_reached": str(row.stable_topology_reached).lower(),
                }
            )
    return path


def write_posterior_agreement_subtree_artifacts(
    out_dir: Path,
    tree: PhyloTree,
    report: PosteriorAgreementSubtreeReport,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    tree_path = write_newick(out_dir / "posterior-agreement-subtree.nwk", tree)
    summary_path = write_posterior_agreement_subtree_summary_table(
        out_dir / "posterior-agreement-subtree-summary.tsv",
        report,
    )
    removed_taxa_path = write_posterior_agreement_subtree_removed_taxa_table(
        out_dir / "posterior-agreement-subtree-removed-taxa.tsv",
        report,
    )
    search_path = write_posterior_agreement_subtree_search_table(
        out_dir / "posterior-agreement-subtree-search.tsv",
        report,
    )
    return {
        "posterior_agreement_subtree_path": tree_path,
        "posterior_agreement_subtree_summary_path": summary_path,
        "posterior_agreement_subtree_removed_taxa_path": removed_taxa_path,
        "posterior_agreement_subtree_search_path": search_path,
    }
