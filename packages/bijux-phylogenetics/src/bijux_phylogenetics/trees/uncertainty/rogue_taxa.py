from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_object_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from ..tree_sets.consensus import _build_consensus_tree_with_threshold_from_trees
from ..tree_sets.inventory import _analyze_tree_set, _require_exact_taxa
from ..tree_sets.topology import _clade_counts, _rooted_topology_id, _tree_distance
from .models import RogueTaxonDetectionReport, RogueTaxonScoreRow

ROGUE_TAXON_RANKING_OBJECTIVE = (
    "maximize-consensus-resolution-then-support-then-rf-stability"
)


def _mean(values: list[float]) -> float:
    if not values:
        return 0.0
    return round(sum(values) / len(values), 15)


def _consensus_resolution(*, included_clade_count: int, taxon_count: int) -> float:
    # Normalize by the maximum informative rooted clades for the current taxon set.
    return round(
        included_clade_count / max(taxon_count - 2, 1),
        15,
    )


def _consensus_metrics_from_trees(
    trees: list[PhyloTree],
    *,
    threshold: float,
) -> tuple[PhyloTree, float, float]:
    exact_taxa = sorted(trees[0].tip_names)
    exact_taxa_set = set(exact_taxa)
    clade_counts = _clade_counts(trees, exact_taxa_set)
    consensus_tree, included_clade_count = (
        _build_consensus_tree_with_threshold_from_trees(
            trees,
            threshold=threshold,
        )
    )
    included_support_percents = [
        round((count / len(trees)) * 100.0, 15)
        for count in clade_counts.values()
        if count / len(trees) >= threshold
    ]
    return (
        consensus_tree,
        _consensus_resolution(
            included_clade_count=included_clade_count,
            taxon_count=len(exact_taxa),
        ),
        _mean(included_support_percents),
    )


def _mean_pairwise_normalized_robinson_foulds(trees: list[PhyloTree]) -> float:
    exact_taxa_set = set(trees[0].tip_names)
    normalized_distances: list[float] = []
    for left_index, left in enumerate(trees):
        for right in trees[left_index + 1 :]:
            _distance, normalized = _tree_distance(left, right, exact_taxa_set)
            normalized_distances.append(normalized)
    return _mean(normalized_distances)


def _rooted_topology_summary(
    trees: list[PhyloTree],
) -> tuple[int, float]:
    exact_taxa_set = set(trees[0].tip_names)
    counts: dict[str, int] = {}
    for tree in trees:
        topology_id = _rooted_topology_id(tree, exact_taxa_set)
        counts[topology_id] = counts.get(topology_id, 0) + 1
    dominant_count = max(counts.values(), default=0)
    dominant_frequency = 0.0 if not trees else round(dominant_count / len(trees), 15)
    return len(counts), dominant_frequency


def detect_rogue_taxa(
    path: Path,
    *,
    consensus_threshold: float = 0.5,
) -> RogueTaxonDetectionReport:
    """Rank taxa by the consensus and RF-stability improvement gained by removing them."""
    if not 0.0 < consensus_threshold <= 1.0:
        raise ValueError(
            "consensus_threshold must be greater than 0 and at most 1, "
            f"got {consensus_threshold}"
        )
    analysis = _analyze_tree_set(path)
    exact_taxa = _require_exact_taxa(analysis)
    if len(exact_taxa) < 4:
        raise InvalidAlignmentError(
            "rogue taxon detection requires at least four taxa shared across the tree set"
        )

    baseline_consensus_tree, baseline_consensus_resolution, baseline_mean_support = (
        _consensus_metrics_from_trees(
            analysis.trees,
            threshold=consensus_threshold,
        )
    )
    baseline_mean_normalized_rf = _mean_pairwise_normalized_robinson_foulds(
        analysis.trees
    )
    baseline_rooted_topology_count, baseline_dominant_topology_frequency = (
        _rooted_topology_summary(analysis.trees)
    )

    scored_rows: list[
        tuple[tuple[float, float, float, float, str], dict[str, object]]
    ] = []
    for taxon in exact_taxa:
        retained_taxa = [name for name in exact_taxa if name != taxon]
        pruned_trees = [
            prune_tree_object_to_requested_taxa(tree, retained_taxa)
            for tree in analysis.trees
        ]
        (
            pruned_consensus_tree,
            pruned_consensus_resolution,
            pruned_mean_support,
        ) = _consensus_metrics_from_trees(
            pruned_trees,
            threshold=consensus_threshold,
        )
        pruned_mean_normalized_rf = _mean_pairwise_normalized_robinson_foulds(
            pruned_trees
        )
        pruned_rooted_topology_count, pruned_dominant_topology_frequency = (
            _rooted_topology_summary(pruned_trees)
        )
        row_values = {
            "taxon": taxon,
            "mean_terminal_branch_length": (
                None
                if taxon not in analysis.terminal_lengths
                else _mean(analysis.terminal_lengths[taxon])
            ),
            "baseline_consensus_resolution": baseline_consensus_resolution,
            "pruned_consensus_resolution": pruned_consensus_resolution,
            "consensus_resolution_delta": round(
                pruned_consensus_resolution - baseline_consensus_resolution,
                15,
            ),
            "baseline_mean_support_percent": baseline_mean_support,
            "pruned_mean_support_percent": pruned_mean_support,
            "mean_support_percent_delta": round(
                pruned_mean_support - baseline_mean_support,
                15,
            ),
            "baseline_mean_normalized_robinson_foulds": baseline_mean_normalized_rf,
            "pruned_mean_normalized_robinson_foulds": pruned_mean_normalized_rf,
            "normalized_robinson_foulds_stability_delta": round(
                baseline_mean_normalized_rf - pruned_mean_normalized_rf,
                15,
            ),
            "baseline_rooted_topology_count": baseline_rooted_topology_count,
            "pruned_rooted_topology_count": pruned_rooted_topology_count,
            "rooted_topology_count_delta": (
                baseline_rooted_topology_count - pruned_rooted_topology_count
            ),
            "baseline_dominant_topology_frequency": (
                baseline_dominant_topology_frequency
            ),
            "pruned_dominant_topology_frequency": pruned_dominant_topology_frequency,
            "dominant_topology_frequency_delta": round(
                pruned_dominant_topology_frequency
                - baseline_dominant_topology_frequency,
                15,
            ),
            "pruned_consensus_newick": dumps_newick(pruned_consensus_tree),
        }
        ranking_key = (
            -float(row_values["consensus_resolution_delta"]),
            -float(row_values["mean_support_percent_delta"]),
            -float(row_values["normalized_robinson_foulds_stability_delta"]),
            float(row_values["pruned_mean_normalized_robinson_foulds"]),
            str(row_values["taxon"]),
        )
        scored_rows.append((ranking_key, row_values))

    scored_rows.sort(key=lambda item: item[0])
    rows = [
        RogueTaxonScoreRow(
            taxon=str(row_values["taxon"]),
            rank=rank,
            mean_terminal_branch_length=float(
                row_values["mean_terminal_branch_length"]
            ),
            baseline_consensus_resolution=float(
                row_values["baseline_consensus_resolution"]
            ),
            pruned_consensus_resolution=float(
                row_values["pruned_consensus_resolution"]
            ),
            consensus_resolution_delta=float(row_values["consensus_resolution_delta"]),
            baseline_mean_support_percent=float(
                row_values["baseline_mean_support_percent"]
            ),
            pruned_mean_support_percent=float(
                row_values["pruned_mean_support_percent"]
            ),
            mean_support_percent_delta=float(row_values["mean_support_percent_delta"]),
            baseline_mean_normalized_robinson_foulds=float(
                row_values["baseline_mean_normalized_robinson_foulds"]
            ),
            pruned_mean_normalized_robinson_foulds=float(
                row_values["pruned_mean_normalized_robinson_foulds"]
            ),
            normalized_robinson_foulds_stability_delta=float(
                row_values["normalized_robinson_foulds_stability_delta"]
            ),
            baseline_rooted_topology_count=int(
                row_values["baseline_rooted_topology_count"]
            ),
            pruned_rooted_topology_count=int(
                row_values["pruned_rooted_topology_count"]
            ),
            rooted_topology_count_delta=int(row_values["rooted_topology_count_delta"]),
            baseline_dominant_topology_frequency=float(
                row_values["baseline_dominant_topology_frequency"]
            ),
            pruned_dominant_topology_frequency=float(
                row_values["pruned_dominant_topology_frequency"]
            ),
            dominant_topology_frequency_delta=float(
                row_values["dominant_topology_frequency_delta"]
            ),
            pruned_consensus_newick=str(row_values["pruned_consensus_newick"]),
        )
        for rank, (_ranking_key, row_values) in enumerate(scored_rows, start=1)
    ]
    return RogueTaxonDetectionReport(
        path=path,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=exact_taxa,
        consensus_threshold=consensus_threshold,
        ranking_objective=ROGUE_TAXON_RANKING_OBJECTIVE,
        baseline_consensus_newick=dumps_newick(baseline_consensus_tree),
        baseline_consensus_resolution=baseline_consensus_resolution,
        baseline_mean_support_percent=baseline_mean_support,
        baseline_mean_normalized_robinson_foulds=baseline_mean_normalized_rf,
        baseline_rooted_topology_count=baseline_rooted_topology_count,
        baseline_dominant_topology_frequency=baseline_dominant_topology_frequency,
        rows=rows,
    )


def write_rogue_taxon_table(path: Path, report: RogueTaxonDetectionReport) -> Path:
    """Write one ranked rogue-taxon table with consensus and RF-stability deltas."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "taxon",
                "mean_terminal_branch_length",
                "baseline_consensus_resolution",
                "pruned_consensus_resolution",
                "consensus_resolution_delta",
                "baseline_mean_support_percent",
                "pruned_mean_support_percent",
                "mean_support_percent_delta",
                "baseline_mean_normalized_robinson_foulds",
                "pruned_mean_normalized_robinson_foulds",
                "normalized_robinson_foulds_stability_delta",
                "baseline_rooted_topology_count",
                "pruned_rooted_topology_count",
                "rooted_topology_count_delta",
                "baseline_dominant_topology_frequency",
                "pruned_dominant_topology_frequency",
                "dominant_topology_frequency_delta",
                "pruned_consensus_newick",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rows:
            writer.writerow(
                {
                    "rank": row.rank,
                    "taxon": row.taxon,
                    "mean_terminal_branch_length": (
                        ""
                        if row.mean_terminal_branch_length is None
                        else format(row.mean_terminal_branch_length, ".15g")
                    ),
                    "baseline_consensus_resolution": format(
                        row.baseline_consensus_resolution, ".15g"
                    ),
                    "pruned_consensus_resolution": format(
                        row.pruned_consensus_resolution, ".15g"
                    ),
                    "consensus_resolution_delta": format(
                        row.consensus_resolution_delta, ".15g"
                    ),
                    "baseline_mean_support_percent": format(
                        row.baseline_mean_support_percent, ".15g"
                    ),
                    "pruned_mean_support_percent": format(
                        row.pruned_mean_support_percent, ".15g"
                    ),
                    "mean_support_percent_delta": format(
                        row.mean_support_percent_delta, ".15g"
                    ),
                    "baseline_mean_normalized_robinson_foulds": format(
                        row.baseline_mean_normalized_robinson_foulds, ".15g"
                    ),
                    "pruned_mean_normalized_robinson_foulds": format(
                        row.pruned_mean_normalized_robinson_foulds, ".15g"
                    ),
                    "normalized_robinson_foulds_stability_delta": format(
                        row.normalized_robinson_foulds_stability_delta, ".15g"
                    ),
                    "baseline_rooted_topology_count": row.baseline_rooted_topology_count,
                    "pruned_rooted_topology_count": row.pruned_rooted_topology_count,
                    "rooted_topology_count_delta": row.rooted_topology_count_delta,
                    "baseline_dominant_topology_frequency": format(
                        row.baseline_dominant_topology_frequency, ".15g"
                    ),
                    "pruned_dominant_topology_frequency": format(
                        row.pruned_dominant_topology_frequency, ".15g"
                    ),
                    "dominant_topology_frequency_delta": format(
                        row.dominant_topology_frequency_delta, ".15g"
                    ),
                    "pruned_consensus_newick": row.pruned_consensus_newick,
                }
            )
    return path
