from __future__ import annotations

import csv
import math
from pathlib import Path

from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from ..tree_sets.contracts import TreeSetReport
from ..tree_sets.inventory import (
    _analyze_tree_set,
    _require_exact_taxa,
    _require_tree_set,
    _TreeSetAnalysis,
    _validate_same_taxa,
)
from ..tree_sets.topology import _tree_distance
from .instability import _build_unstable_clade_report
from .models import (
    PosteriorTopologicalDiversityComparisonReport,
    PosteriorTopologicalDiversitySummary,
    PosteriorTopologyDiversityReport,
    PosteriorTopologyMode,
    PosteriorTopologyMultimodalityReport,
    TreeDistanceDistributionRow,
    TreeTopologyCluster,
    TreeTopologyClusterReport,
)


def _shannon_effective_count(frequencies: list[float]) -> float:
    if not frequencies:
        return 0.0
    return round(
        math.exp(
            -sum(
                frequency * math.log(frequency)
                for frequency in frequencies
                if frequency > 0.0
            )
        ),
        15,
    )


def _mean_pairwise_distance(
    trees: list[PhyloTree],
    shared_taxa: set[str],
) -> tuple[float, float]:
    comparisons: list[tuple[int, float]] = []
    for left_index, left in enumerate(trees):
        for right in trees[left_index + 1 :]:
            comparisons.append(_tree_distance(left, right, shared_taxa))
    if not comparisons:
        return 0.0, 0.0
    return (
        round(sum(distance for distance, _ in comparisons) / len(comparisons), 15),
        round(sum(normalized for _, normalized in comparisons) / len(comparisons), 15),
    )


def _topology_modes_from_clusters(
    trees: list[PhyloTree],
    clusters: list[TreeTopologyCluster],
    *,
    min_mode_frequency: float,
) -> list[PosteriorTopologyMode]:
    return [
        PosteriorTopologyMode(
            rooted_topology_id=cluster.rooted_topology_id,
            representative_index=cluster.representative_index,
            representative_newick=cluster.representative_newick,
            tree_indices=cluster.tree_indices,
            tree_count=cluster.tree_count,
            frequency=cluster.frequency,
        )
        for cluster in clusters
        if cluster.frequency >= min_mode_frequency
    ]


def _build_topology_cluster_report(
    analysis: _TreeSetAnalysis,
) -> TreeTopologyClusterReport:
    tree_count = len(analysis.trees)
    clusters: list[TreeTopologyCluster] = []
    indices_by_topology: dict[str, list[int]] = {}
    for record in analysis.records:
        indices_by_topology.setdefault(record.rooted_topology_id, []).append(
            record.index
        )
    for topology_id, indices in sorted(
        indices_by_topology.items(),
        key=lambda item: (-len(item[1]), item[1][0]),
    ):
        representative_index, representative_newick, _tree = (
            analysis.rooted_representatives[topology_id]
        )
        clusters.append(
            TreeTopologyCluster(
                rooted_topology_id=topology_id,
                tree_indices=indices,
                tree_count=len(indices),
                frequency=round(len(indices) / tree_count, 15),
                representative_index=representative_index,
                representative_newick=representative_newick,
            )
        )
    return TreeTopologyClusterReport(
        path=analysis.path,
        tree_count=tree_count,
        processing=analysis.processing,
        rooted_topology_count=len(analysis.rooted_topology_counts),
        clusters=clusters,
    )


def _rf_distribution_from_analysis(
    analysis: _TreeSetAnalysis,
) -> tuple[list[TreeDistanceDistributionRow], int, float, float, int, float]:
    exact_taxa_set = set(_require_exact_taxa(analysis))
    representatives = [
        (
            topology_id,
            analysis.rooted_topology_counts[topology_id],
            analysis.rooted_representatives[topology_id][2],
        )
        for topology_id in sorted(
            analysis.rooted_representatives,
            key=lambda topology_id: analysis.rooted_representatives[topology_id][0],
        )
    ]
    pair_counts: dict[tuple[int, float], int] = {}
    total_pairs = 0
    for index, (_left_id, left_count, left_tree) in enumerate(representatives):
        for right_index, (_right_id, right_count, right_tree) in enumerate(
            representatives[index:], start=index
        ):
            if right_index == index:
                pair_count = (left_count * (left_count - 1)) // 2
            else:
                pair_count = left_count * right_count
            if pair_count == 0:
                continue
            distance, normalized = _tree_distance(left_tree, right_tree, exact_taxa_set)
            key = (distance, normalized)
            pair_counts[key] = pair_counts.get(key, 0) + pair_count
            total_pairs += pair_count
    if total_pairs == 0:
        return [], 0, 0.0, 0.0, 0, 0.0
    rows = [
        TreeDistanceDistributionRow(
            robinson_foulds_distance=distance,
            normalized_robinson_foulds=normalized,
            pair_count=count,
            frequency=round(count / total_pairs, 15),
        )
        for (distance, normalized), count in sorted(pair_counts.items())
    ]
    mean_rf = round(
        sum(row.robinson_foulds_distance * row.pair_count for row in rows)
        / total_pairs,
        15,
    )
    mean_normalized_rf = round(
        sum(row.normalized_robinson_foulds * row.pair_count for row in rows)
        / total_pairs,
        15,
    )
    maximum_rf = max(row.robinson_foulds_distance for row in rows)
    maximum_normalized_rf = round(
        max(row.normalized_robinson_foulds for row in rows),
        15,
    )
    return (
        rows,
        total_pairs,
        mean_rf,
        mean_normalized_rf,
        maximum_rf,
        maximum_normalized_rf,
    )


def _build_posterior_topology_diversity_report(
    analysis: _TreeSetAnalysis,
) -> PosteriorTopologyDiversityReport:
    summary = TreeSetReport(
        path=analysis.path,
        source_format=analysis.source_format,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=analysis.shared_taxa,
        taxa_union=analysis.taxa_union,
        rooted_topology_count=len(analysis.rooted_topology_counts),
        unrooted_topology_count=len(analysis.unrooted_topology_counts),
        records=analysis.records,
    )
    clusters = _build_topology_cluster_report(analysis)
    unstable_clades = _build_unstable_clade_report(analysis)
    (
        distribution,
        pair_count,
        mean_rf,
        mean_normalized_rf,
        maximum_rf,
        maximum_normalized_rf,
    ) = _rf_distribution_from_analysis(analysis)
    dominant_topology_frequency = (
        0.0 if not clusters.clusters else clusters.clusters[0].frequency
    )
    return PosteriorTopologyDiversityReport(
        path=analysis.path,
        tree_count=summary.tree_count,
        processing=analysis.processing,
        rooted_topology_count=summary.rooted_topology_count,
        dominant_topology_frequency=dominant_topology_frequency,
        effective_topology_count=_shannon_effective_count(
            [cluster.frequency for cluster in clusters.clusters]
        ),
        pair_count=pair_count,
        mean_robinson_foulds_distance=mean_rf,
        mean_normalized_robinson_foulds_distance=mean_normalized_rf,
        maximum_robinson_foulds_distance=maximum_rf,
        maximum_normalized_robinson_foulds_distance=maximum_normalized_rf,
        unstable_clade_count=len(unstable_clades.clades),
        rf_distribution=distribution,
    )


def write_tree_distance_distribution_table(
    path: Path,
    report: PosteriorTopologyDiversityReport,
) -> Path:
    """Write the pairwise RF-distance distribution as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
                "pair_count",
                "frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.rf_distribution:
            writer.writerow(
                {
                    "robinson_foulds_distance": row.robinson_foulds_distance,
                    "normalized_robinson_foulds": format(
                        row.normalized_robinson_foulds,
                        ".15g",
                    ),
                    "pair_count": row.pair_count,
                    "frequency": format(row.frequency, ".15g"),
                }
            )
    return path


def summarize_posterior_topology_diversity(
    path: Path,
) -> PosteriorTopologyDiversityReport:
    """Summarize topology dispersion and instability across one posterior tree set."""
    return _build_posterior_topology_diversity_report(_analyze_tree_set(path))


def write_topology_cluster_table(path: Path, report: TreeTopologyClusterReport) -> Path:
    """Write rooted topology clusters as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rooted_topology_id",
                "tree_indices",
                "tree_count",
                "frequency",
                "representative_index",
                "representative_newick",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.clusters:
            writer.writerow(
                {
                    "rooted_topology_id": row.rooted_topology_id,
                    "tree_indices": ",".join(str(index) for index in row.tree_indices),
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "representative_index": row.representative_index,
                    "representative_newick": row.representative_newick,
                }
            )
    return path


def cluster_trees_by_topology(path: Path) -> TreeTopologyClusterReport:
    """Cluster trees by identical rooted topology signatures."""
    return _build_topology_cluster_report(_analyze_tree_set(path))


def compare_posterior_topological_diversity(
    left_path: Path,
    right_path: Path,
) -> PosteriorTopologicalDiversityComparisonReport:
    """Compare topology diversity and dispersion across two posterior tree sets."""
    left_clusters = cluster_trees_by_topology(left_path)
    right_clusters = cluster_trees_by_topology(right_path)
    _, left_trees = _require_tree_set(left_path)
    _, right_trees = _require_tree_set(right_path)
    left_taxa = set(_validate_same_taxa(left_trees))
    right_taxa = set(_validate_same_taxa(right_trees))
    if left_taxa != right_taxa:
        raise InvalidAlignmentError(
            "posterior diversity comparison requires identical taxon sets across both inputs"
        )
    left_mean_rf, left_mean_normalized = _mean_pairwise_distance(left_trees, left_taxa)
    right_mean_rf, right_mean_normalized = _mean_pairwise_distance(
        right_trees,
        right_taxa,
    )
    left_summary = PosteriorTopologicalDiversitySummary(
        tree_count=len(left_trees),
        rooted_topology_count=left_clusters.rooted_topology_count,
        dominant_topology_frequency=left_clusters.clusters[0].frequency,
        effective_topology_count=_shannon_effective_count(
            [cluster.frequency for cluster in left_clusters.clusters]
        ),
        mean_within_set_robinson_foulds=left_mean_rf,
        mean_within_set_normalized_robinson_foulds=left_mean_normalized,
    )
    right_summary = PosteriorTopologicalDiversitySummary(
        tree_count=len(right_trees),
        rooted_topology_count=right_clusters.rooted_topology_count,
        dominant_topology_frequency=right_clusters.clusters[0].frequency,
        effective_topology_count=_shannon_effective_count(
            [cluster.frequency for cluster in right_clusters.clusters]
        ),
        mean_within_set_robinson_foulds=right_mean_rf,
        mean_within_set_normalized_robinson_foulds=right_mean_normalized,
    )
    warnings: list[str] = []
    if left_summary.effective_topology_count != right_summary.effective_topology_count:
        richer = (
            "left"
            if left_summary.effective_topology_count
            > right_summary.effective_topology_count
            else "right"
        )
        warnings.append(
            f"{richer} analysis spans a broader effective topology spectrum"
        )
    if (
        abs(
            left_summary.mean_within_set_normalized_robinson_foulds
            - right_summary.mean_within_set_normalized_robinson_foulds
        )
        >= 0.15
    ):
        warnings.append(
            "posterior analyses differ materially in within-set topological dispersion"
        )
    return PosteriorTopologicalDiversityComparisonReport(
        left_path=left_path,
        right_path=right_path,
        left_summary=left_summary,
        right_summary=right_summary,
        warnings=warnings,
    )


def detect_posterior_topology_multimodality(
    path: Path,
    *,
    min_mode_frequency: float = 0.2,
    min_mode_count: int = 2,
) -> PosteriorTopologyMultimodalityReport:
    """Report whether a posterior tree set contains multiple high-frequency topology modes."""
    if not 0.0 < min_mode_frequency <= 1.0:
        raise ValueError(
            f"min_mode_frequency must be between 0 and 1, got {min_mode_frequency}"
        )
    if min_mode_count < 2:
        raise ValueError(f"min_mode_count must be at least 2, got {min_mode_count}")
    clusters = cluster_trees_by_topology(path)
    _, trees = _require_tree_set(path)
    modes = _topology_modes_from_clusters(
        trees,
        clusters.clusters,
        min_mode_frequency=min_mode_frequency,
    )
    multimodal = len(modes) >= min_mode_count
    warnings: list[str] = []
    if multimodal:
        warnings.append(
            "posterior topology distribution contains multiple high-frequency modes"
        )
    if clusters.clusters and clusters.clusters[0].frequency < 0.75:
        warnings.append("no single topology dominates the posterior tree set")
    return PosteriorTopologyMultimodalityReport(
        path=path,
        tree_count=clusters.tree_count,
        rooted_topology_count=clusters.rooted_topology_count,
        dominant_mode_frequency=0.0
        if not clusters.clusters
        else clusters.clusters[0].frequency,
        mode_count=len(modes),
        multimodal=multimodal,
        modes=modes,
        warnings=warnings,
    )
