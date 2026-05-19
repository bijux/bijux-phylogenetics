from __future__ import annotations

import csv
import math
from pathlib import Path
import tempfile
from time import perf_counter

from bijux_phylogenetics.core.clade_sets import (
    informative_rooted_clade_nodes,
    informative_rooted_clades,
)
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.io.iqtree_support import (
    parse_iqtree_branch_support_label,
)
from bijux_phylogenetics.io.iqtree_support import (
    support_fraction as normalize_support_fraction,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

from .models import (
    BootstrapPosteriorCladeComparison,
    BootstrapPosteriorSupportComparisonReport,
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    BootstrapUnstableBranch,
    CladeCredibilityConflict,
    CladeCredibilityConflictReport,
    CladeFrequencyDelta,
    ConsensusThresholdSensitivityReport,
    ConsensusThresholdSensitivityRow,
    PosteriorTopologicalDiversityComparisonReport,
    PosteriorTopologicalDiversitySummary,
    PosteriorTopologyDiversityReport,
    PosteriorTopologyMode,
    PosteriorTopologyMultimodalityReport,
    PosteriorTreeSetComparisonReport,
    TaxonPlacementSignature,
    TreeDistanceDistributionRow,
    TreeSetBenchmarkRow,
    TreeSetMaturityGateCheck,
    TreeSetMaturityGateReport,
    TreeSetScalingBenchmarkReport,
    TreeSetStorageRiskReport,
    TreeSetThinningSensitivityReport,
    TreeSetThinningSensitivityRow,
    TreeTopologyCluster,
    TreeTopologyClusterReport,
    UncertaintyAwareCladeConclusion,
    UncertaintyAwareConclusionSummaryReport,
    UnstableClade,
    UnstableCladeReport,
    UnstableTaxaReport,
    UnstableTaxon,
)
from ..tree_set import (
    CladeFrequencyReport,
    ConsensusTreeReport,
    TreeSetReport,
    _TreeSetAnalysis,
    _analyze_tree_set,
    _build_clade_frequency_report,
    _build_consensus_tree_with_threshold,
    _build_tree_distance_matrix_report,
    _clade_counts,
    _clade_signature,
    _clades_conflict,
    _format_clade,
    _require_exact_taxa,
    _require_tree_set,
    _rooted_topology_id,
    _support_classification,
    _tree_distance,
    _validate_same_taxa,
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
    compute_clade_frequency_table,
    compute_consensus_tree_with_threshold,
    enforce_tree_set_tree_budget,
    load_tree_set,
    write_clade_frequency_table,
    write_consensus_tree,
    write_tree_distance_matrix,
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
    trees: list[PhyloTree], shared_taxa: set[str]
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


def _build_unstable_clade_report(analysis: _TreeSetAnalysis) -> UnstableCladeReport:
    _require_exact_taxa(analysis)
    counts = analysis.clade_counts or {}
    all_clades = set(counts)
    tree_count = len(analysis.trees)
    unstable_clades = [
        UnstableClade(
            clade=_format_clade(clade),
            tree_count=count,
            frequency=round(count / tree_count, 15),
            conflict_count=len(
                conflicts := sorted(
                    _format_clade(other)
                    for other in all_clades
                    if _clades_conflict(clade, other)
                )
            ),
            instability_score=round(
                min(count / tree_count, 1.0 - (count / tree_count)), 15
            ),
            support_classification=_support_classification(
                round(count / tree_count, 15),
                len(conflicts),
            ),
            conflicting_clades=conflicts,
        )
        for clade, count in sorted(
            counts.items(), key=lambda item: _format_clade(item[0])
        )
        if count < tree_count
    ]
    unstable_clades.sort(
        key=lambda row: (-row.instability_score, -row.conflict_count, row.clade)
    )
    return UnstableCladeReport(
        path=analysis.path,
        tree_count=tree_count,
        processing=analysis.processing,
        clades=unstable_clades,
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


def detect_unstable_taxa(path: Path) -> UnstableTaxaReport:
    """Report taxa whose placement signatures vary across trees in a set."""
    analysis = _analyze_tree_set(path)
    trees = analysis.trees
    shared_taxa = set(_require_exact_taxa(analysis))
    taxa: list[UnstableTaxon] = []
    for taxon in sorted(shared_taxa):
        signature_counts: dict[str, int] = {}
        for tree in trees:
            signature = _clade_signature(tree, shared_taxa, taxon)
            signature_counts[signature] = signature_counts.get(signature, 0) + 1
        if len(signature_counts) < 2:
            continue
        placements = [
            TaxonPlacementSignature(
                signature=signature,
                tree_count=count,
                frequency=round(count / len(trees), 15),
            )
            for signature, count in sorted(
                signature_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
        ]
        taxa.append(
            UnstableTaxon(
                taxon=taxon,
                unique_placements=len(signature_counts),
                dominant_frequency=placements[0].frequency,
                instability_score=round(1.0 - placements[0].frequency, 15),
                placements=placements,
            )
        )
    taxa.sort(
        key=lambda row: (-row.instability_score, -row.unique_placements, row.taxon)
    )
    return UnstableTaxaReport(
        path=path,
        tree_count=len(trees),
        processing=analysis.processing,
        taxa=taxa,
    )


def detect_unstable_clades(path: Path) -> UnstableCladeReport:
    """Report non-unanimous clades and their conflicting alternatives."""
    return _build_unstable_clade_report(_analyze_tree_set(path))


def summarize_bootstrap_tree_set(
    path: Path,
    *,
    consensus_threshold: float = 0.5,
    robust_support_threshold: float = 0.9,
) -> BootstrapTreeSetSummaryReport:
    """Summarize bootstrap replicate trees through one review-oriented report."""
    return _build_bootstrap_tree_set_summary_report(
        _analyze_tree_set(path),
        consensus_threshold=consensus_threshold,
        robust_support_threshold=robust_support_threshold,
    )


def _build_bootstrap_tree_set_summary_report(
    analysis: _TreeSetAnalysis,
    *,
    consensus_threshold: float = 0.5,
    robust_support_threshold: float = 0.9,
) -> BootstrapTreeSetSummaryReport:
    if not 0.0 < consensus_threshold < 1.0:
        raise ValueError(
            f"consensus_threshold must be between 0 and 1, got {consensus_threshold}"
        )
    if not 0.0 < robust_support_threshold <= 1.0:
        raise ValueError(
            "robust_support_threshold must be between 0 and 1, "
            f"got {robust_support_threshold}"
        )
    path = analysis.path
    summary = TreeSetReport(
        path=path,
        source_format=analysis.source_format,
        tree_count=len(analysis.trees),
        processing=analysis.processing,
        shared_taxa=analysis.shared_taxa,
        taxa_union=analysis.taxa_union,
        rooted_topology_count=len(analysis.rooted_topology_counts),
        unrooted_topology_count=len(analysis.unrooted_topology_counts),
        records=analysis.records,
    )
    clade_frequencies = _build_clade_frequency_report(analysis)
    consensus_tree, consensus = _build_consensus_tree_with_threshold(
        analysis,
        threshold=consensus_threshold,
    )
    diversity = _build_posterior_topology_diversity_report(analysis)
    unstable_clades = _build_unstable_clade_report(analysis)
    shared_taxa = set(summary.shared_taxa)
    consensus_clades = informative_rooted_clade_nodes(consensus_tree, shared_taxa)
    frequencies_by_clade = {
        row.clade: row for row in clade_frequencies.clade_frequencies
    }
    unstable_by_clade = {row.clade: row for row in unstable_clades.clades}
    unstable_branches: list[BootstrapUnstableBranch] = []
    for clade in sorted(consensus_clades, key=_format_clade):
        clade_id = _format_clade(clade)
        frequency = frequencies_by_clade[clade_id]
        unstable_row = unstable_by_clade.get(clade_id)
        conflict_count = 0 if unstable_row is None else unstable_row.conflict_count
        instability_score = (
            0.0 if unstable_row is None else unstable_row.instability_score
        )
        support_classification = _support_classification(
            frequency.frequency, conflict_count
        )
        if frequency.frequency >= robust_support_threshold and conflict_count == 0:
            continue
        unstable_branches.append(
            BootstrapUnstableBranch(
                clade=clade_id,
                bootstrap_tree_count=frequency.tree_count,
                bootstrap_frequency=frequency.frequency,
                bootstrap_support_percent=round(frequency.frequency * 100.0, 15),
                conflict_count=conflict_count,
                instability_score=instability_score,
                support_classification=support_classification,
                conflicting_clades=(
                    [] if unstable_row is None else unstable_row.conflicting_clades
                ),
            )
        )
    unstable_branches.sort(
        key=lambda row: (
            -row.instability_score,
            row.bootstrap_frequency,
            -row.conflict_count,
            row.clade,
        )
    )
    warnings: list[str] = []
    if diversity.rooted_topology_count > 1:
        warnings.append("bootstrap replicate trees contain multiple rooted topologies")
    if unstable_branches:
        warnings.append(
            "consensus tree contains branches below the robust bootstrap threshold or with conflicting alternatives"
        )
    return BootstrapTreeSetSummaryReport(
        path=path,
        consensus_threshold=consensus_threshold,
        robust_support_threshold=robust_support_threshold,
        tree_count=summary.tree_count,
        processing=analysis.processing,
        shared_taxa=summary.shared_taxa,
        summary=summary,
        clade_frequencies=clade_frequencies,
        consensus=consensus,
        diversity=diversity,
        unstable_clades=unstable_clades,
        unstable_branch_count=len(unstable_branches),
        unstable_branches=unstable_branches,
        warnings=warnings,
    )


def write_bootstrap_tree_set_summary_table(
    path: Path, report: BootstrapTreeSetSummaryReport
) -> Path:
    """Write a one-row TSV summary for one bootstrap tree set."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tree_count",
                "runtime_seconds",
                "peak_memory_bytes",
                "skipped_malformed_tree_count",
                "shared_taxon_count",
                "rooted_topology_count",
                "dominant_topology_frequency",
                "effective_topology_count",
                "mean_robinson_foulds_distance",
                "mean_normalized_robinson_foulds_distance",
                "consensus_threshold",
                "robust_support_threshold",
                "unstable_branch_count",
                "warning_count",
                "consensus_newick",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "tree_count": report.tree_count,
                "runtime_seconds": format(report.processing.runtime_seconds, ".15g"),
                "peak_memory_bytes": report.processing.peak_memory_bytes,
                "skipped_malformed_tree_count": (
                    report.processing.skipped_malformed_tree_count
                ),
                "shared_taxon_count": len(report.shared_taxa),
                "rooted_topology_count": report.diversity.rooted_topology_count,
                "dominant_topology_frequency": format(
                    report.diversity.dominant_topology_frequency, ".15g"
                ),
                "effective_topology_count": format(
                    report.diversity.effective_topology_count, ".15g"
                ),
                "mean_robinson_foulds_distance": format(
                    report.diversity.mean_robinson_foulds_distance, ".15g"
                ),
                "mean_normalized_robinson_foulds_distance": format(
                    report.diversity.mean_normalized_robinson_foulds_distance, ".15g"
                ),
                "consensus_threshold": format(report.consensus_threshold, ".15g"),
                "robust_support_threshold": format(
                    report.robust_support_threshold, ".15g"
                ),
                "unstable_branch_count": report.unstable_branch_count,
                "warning_count": len(report.warnings),
                "consensus_newick": report.consensus.consensus_newick,
            }
        )
    return path


def write_bootstrap_unstable_branch_table(
    path: Path, report: BootstrapTreeSetSummaryReport
) -> Path:
    """Write consensus-branch instability evidence for one bootstrap tree set."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade",
                "bootstrap_tree_count",
                "bootstrap_frequency",
                "bootstrap_support_percent",
                "conflict_count",
                "instability_score",
                "support_classification",
                "conflicting_clades",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.unstable_branches:
            writer.writerow(
                {
                    "clade": row.clade,
                    "bootstrap_tree_count": row.bootstrap_tree_count,
                    "bootstrap_frequency": format(row.bootstrap_frequency, ".15g"),
                    "bootstrap_support_percent": format(
                        row.bootstrap_support_percent, ".15g"
                    ),
                    "conflict_count": row.conflict_count,
                    "instability_score": format(row.instability_score, ".15g"),
                    "support_classification": row.support_classification,
                    "conflicting_clades": ",".join(row.conflicting_clades),
                }
            )
    return path


def write_bootstrap_tree_set_artifacts(
    tree_set_path: Path,
    *,
    out_dir: Path,
    prefix: str = "bootstrap-tree-set",
    consensus_threshold: float = 0.5,
    robust_support_threshold: float = 0.9,
    max_tree_count: int | None = None,
    memory_warning_threshold_bytes: int | None = None,
) -> BootstrapTreeSetArtifactReport:
    """Write a governed artifact set for one bootstrap replicate tree file."""
    budget = build_tree_set_workflow_budget(
        max_tree_count=max_tree_count,
        memory_warning_threshold_bytes=memory_warning_threshold_bytes,
    )
    analysis = _analyze_tree_set(tree_set_path)
    enforce_tree_set_tree_budget(
        tree_count=len(analysis.trees),
        budget=budget,
        workflow_name="bootstrap tree-set artifact workflow",
        source_path=tree_set_path,
    )
    summary_report = _build_bootstrap_tree_set_summary_report(
        analysis,
        consensus_threshold=consensus_threshold,
        robust_support_threshold=robust_support_threshold,
    )
    consensus_tree, _ = _build_consensus_tree_with_threshold(
        analysis,
        threshold=consensus_threshold,
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    base_path = out_dir / prefix
    output_paths = {
        "summary_table": write_bootstrap_tree_set_summary_table(
            base_path.with_suffix(".summary.tsv"), summary_report
        ),
        "consensus_tree": write_consensus_tree(
            base_path.with_suffix(".consensus.nwk"), consensus_tree
        ),
        "clade_frequencies": write_clade_frequency_table(
            base_path.with_suffix(".clade-frequencies.tsv"),
            summary_report.clade_frequencies,
        ),
        "unstable_branches": write_bootstrap_unstable_branch_table(
            base_path.with_suffix(".unstable-branches.tsv"), summary_report
        ),
        "unstable_clades": write_unstable_clade_table(
            base_path.with_suffix(".unstable-clades.tsv"),
            summary_report.unstable_clades,
        ),
        "distance_matrix": write_tree_distance_matrix(
            base_path.with_suffix(".distance-matrix.tsv"),
            _build_tree_distance_matrix_report(analysis),
        ),
        "rf_distribution": write_tree_distance_distribution_table(
            base_path.with_suffix(".rf-distribution.tsv"),
            summary_report.diversity,
        ),
        "topology_clusters": write_topology_cluster_table(
            base_path.with_suffix(".topology-clusters.tsv"),
            _build_topology_cluster_report(analysis),
        ),
    }
    budget_report = build_tree_set_budget_report(
        budget=budget,
        peak_memory_bytes=summary_report.processing.peak_memory_bytes,
    )
    return BootstrapTreeSetArtifactReport(
        input_path=tree_set_path,
        out_dir=out_dir,
        prefix=prefix,
        summary_report=summary_report,
        budget_report=budget_report,
        output_paths=output_paths,
    )


def write_unstable_clade_table(path: Path, report: UnstableCladeReport) -> Path:
    """Write unstable clades and their conflicting alternatives as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade",
                "tree_count",
                "frequency",
                "conflict_count",
                "instability_score",
                "support_classification",
                "conflicting_clades",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.clades:
            writer.writerow(
                {
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "conflict_count": row.conflict_count,
                    "instability_score": format(row.instability_score, ".15g"),
                    "support_classification": row.support_classification,
                    "conflicting_clades": ",".join(row.conflicting_clades),
                }
            )
    return path


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
        right_trees, right_taxa
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
        trees, clusters.clusters, min_mode_frequency=min_mode_frequency
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


def summarize_clade_credibility_conflicts(
    path: Path,
    *,
    credibility_threshold: float = 0.5,
) -> CladeCredibilityConflictReport:
    """Identify mutually incompatible clades that both achieve high posterior credibility."""
    if not 0.0 < credibility_threshold < 1.0:
        raise ValueError(
            f"credibility_threshold must be between 0 and 1, got {credibility_threshold}"
        )
    _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    counts = _clade_counts(trees, shared_taxa)
    frequencies = {
        clade: round(count / len(trees), 15) for clade, count in counts.items()
    }
    high_credibility = [
        clade
        for clade, frequency in frequencies.items()
        if frequency >= credibility_threshold
    ]
    conflicts: list[CladeCredibilityConflict] = []
    for index, left_clade in enumerate(sorted(high_credibility, key=_format_clade)):
        for right_clade in sorted(high_credibility[index + 1 :], key=_format_clade):
            if not _clades_conflict(left_clade, right_clade):
                continue
            conflicts.append(
                CladeCredibilityConflict(
                    left_clade=_format_clade(left_clade),
                    left_frequency=frequencies[left_clade],
                    right_clade=_format_clade(right_clade),
                    right_frequency=frequencies[right_clade],
                    combined_frequency=round(
                        frequencies[left_clade] + frequencies[right_clade], 15
                    ),
                )
            )
    conflicts.sort(
        key=lambda row: (-row.combined_frequency, row.left_clade, row.right_clade)
    )
    return CladeCredibilityConflictReport(
        path=path,
        tree_count=len(trees),
        credibility_threshold=credibility_threshold,
        high_credibility_clade_count=len(high_credibility),
        conflict_count=len(conflicts),
        conflicts=conflicts,
    )


def write_clade_credibility_conflict_table(
    path: Path, report: CladeCredibilityConflictReport
) -> Path:
    """Write high-credibility clade conflicts as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_clade",
                "left_frequency",
                "right_clade",
                "right_frequency",
                "combined_frequency",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.conflicts:
            writer.writerow(
                {
                    "left_clade": row.left_clade,
                    "left_frequency": format(row.left_frequency, ".15g"),
                    "right_clade": row.right_clade,
                    "right_frequency": format(row.right_frequency, ".15g"),
                    "combined_frequency": format(row.combined_frequency, ".15g"),
                }
            )
    return path


def summarize_uncertainty_aware_conclusions(
    path: Path,
    *,
    robust_threshold: float = 0.9,
    uncertain_min_frequency: float = 0.3,
    uncertain_max_frequency: float = 0.7,
    credibility_threshold: float = 0.5,
) -> UncertaintyAwareConclusionSummaryReport:
    """Classify clade-level conclusions as robust, uncertain, or conflict-prone."""
    if not 0.0 < robust_threshold <= 1.0:
        raise ValueError(
            f"robust_threshold must be between 0 and 1, got {robust_threshold}"
        )
    detect_unstable_clades(path)
    conflict_report = summarize_clade_credibility_conflicts(
        path, credibility_threshold=credibility_threshold
    )
    conflict_clades = {row.left_clade for row in conflict_report.conflicts} | {
        row.right_clade for row in conflict_report.conflicts
    }
    frequency_report = compute_clade_frequency_table(path)
    robust_clades: list[UncertaintyAwareCladeConclusion] = []
    uncertain_clades: list[UncertaintyAwareCladeConclusion] = []
    conflicting_clades: list[UncertaintyAwareCladeConclusion] = []
    for row in frequency_report.clade_frequencies:
        if row.clade in conflict_clades:
            conflicting_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="conflict-prone",
                    rationale="clade reaches high posterior frequency but is incompatible with another high-credibility clade",
                )
            )
            continue
        if row.frequency >= robust_threshold:
            robust_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="robust",
                    rationale="clade remains near-fixed across the posterior tree set",
                )
            )
            continue
        if uncertain_min_frequency <= row.frequency <= uncertain_max_frequency:
            uncertain_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="uncertain",
                    rationale="clade holds intermediate support and should not anchor strong biological interpretation",
                )
            )
    robust_clades.sort(key=lambda row: (-row.frequency, row.clade))
    uncertain_clades.sort(key=lambda row: (-row.frequency, row.clade))
    conflicting_clades.sort(key=lambda row: (-row.frequency, row.clade))
    return UncertaintyAwareConclusionSummaryReport(
        path=path,
        tree_count=frequency_report.tree_count,
        robust_clade_count=len(robust_clades),
        uncertain_clade_count=len(uncertain_clades),
        conflicting_clade_count=len(conflicting_clades),
        robust_clades=robust_clades,
        uncertain_clades=uncertain_clades,
        conflicting_clades=conflicting_clades,
    )


def write_uncertainty_conclusion_table(
    path: Path, report: UncertaintyAwareConclusionSummaryReport
) -> Path:
    """Write robust, uncertain, and conflict-prone clade conclusions as a TSV table."""
    rows = [
        *report.robust_clades,
        *report.uncertain_clades,
        *report.conflicting_clades,
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["clade", "frequency", "conclusion", "rationale"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "clade": row.clade,
                    "frequency": format(row.frequency, ".15g"),
                    "conclusion": row.conclusion,
                    "rationale": row.rationale,
                }
            )
    return path


def compare_posterior_tree_sets(
    left_path: Path, right_path: Path
) -> PosteriorTreeSetComparisonReport:
    """Compare two tree sets over shared taxa, clade support, and cross-set topology distance."""
    _, left_trees = _require_tree_set(left_path)
    _, right_trees = _require_tree_set(right_path)
    left_taxa = _validate_same_taxa(left_trees)
    right_taxa = _validate_same_taxa(right_trees)
    if left_taxa != right_taxa:
        raise InvalidAlignmentError(
            "posterior tree-set comparison requires identical taxon sets across both inputs"
        )

    shared_taxa = set(left_taxa)
    left_counts = _clade_counts(left_trees, shared_taxa)
    right_counts = _clade_counts(right_trees, shared_taxa)
    all_clades = left_counts.keys() | right_counts.keys()
    deltas = [
        CladeFrequencyDelta(
            clade=_format_clade(clade),
            left_frequency=round(left_counts.get(clade, 0) / len(left_trees), 15),
            right_frequency=round(right_counts.get(clade, 0) / len(right_trees), 15),
            delta=round(
                (right_counts.get(clade, 0) / len(right_trees))
                - (left_counts.get(clade, 0) / len(left_trees)),
                15,
            ),
        )
        for clade in sorted(all_clades, key=_format_clade)
    ]
    comparisons = [
        _tree_distance(left, right, shared_taxa)
        for left in left_trees
        for right in right_trees
    ]
    left_topologies = {_rooted_topology_id(tree, shared_taxa) for tree in left_trees}
    right_topologies = {_rooted_topology_id(tree, shared_taxa) for tree in right_trees}
    return PosteriorTreeSetComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=left_taxa,
        left_tree_count=len(left_trees),
        right_tree_count=len(right_trees),
        left_rooted_topology_count=len(left_topologies),
        right_rooted_topology_count=len(right_topologies),
        shared_rooted_topology_count=len(left_topologies & right_topologies),
        mean_between_set_robinson_foulds=round(
            sum(distance for distance, _ in comparisons) / len(comparisons), 15
        ),
        mean_between_set_normalized_robinson_foulds=round(
            sum(normalized for _, normalized in comparisons) / len(comparisons),
            15,
        ),
        clade_frequency_deltas=deltas,
    )


def _require_tree(path: Path) -> PhyloTree:
    if not path.exists():
        raise FileNotFoundError(f"tree file not found: {path}")
    return load_tree(path)


def _parse_support_label(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    parsed_iqtree_label = parse_iqtree_branch_support_label(text)
    if parsed_iqtree_label is not None:
        support_value = (
            parsed_iqtree_label.ufboot_support
            if parsed_iqtree_label.ufboot_support is not None
            else parsed_iqtree_label.sh_alrt_support
        )
        if support_value is None:
            return None
        normalized = normalize_support_fraction(support_value)
        return None if normalized is None else round(normalized, 15)
    try:
        parsed = float(text)
    except ValueError:
        return None
    return round(parsed / 100.0, 15) if parsed > 1.0 else round(parsed, 15)


def _support_agreement_label(
    bootstrap_support: float | None,
    posterior_frequency: float | None,
    absolute_delta: float | None,
) -> str:
    if bootstrap_support is None and posterior_frequency is None:
        return "not_observed"
    if bootstrap_support is None or posterior_frequency is None:
        return "method_specific"
    if absolute_delta is None:
        return "not_comparable"
    if absolute_delta >= 0.35:
        return "strong_conflict"
    if absolute_delta >= 0.15:
        return "moderate_difference"
    return "broad_agreement"


def compare_bootstrap_and_posterior_uncertainty(
    bootstrap_tree_path: Path,
    posterior_tree_set_path: Path,
) -> BootstrapPosteriorSupportComparisonReport:
    """Compare bootstrap support on one summary tree against posterior clade frequencies from a tree set."""
    bootstrap_tree = _require_tree(bootstrap_tree_path)
    posterior_report = compute_clade_frequency_table(posterior_tree_set_path)
    shared_taxa = set(bootstrap_tree.tip_names)
    posterior_taxa = set(posterior_report.shared_taxa)
    if shared_taxa != posterior_taxa:
        raise InvalidAlignmentError(
            "bootstrap versus posterior comparison requires identical taxon sets"
        )
    bootstrap_nodes = informative_rooted_clade_nodes(bootstrap_tree, shared_taxa)
    bootstrap_support_by_clade = {
        _format_clade(clade): _parse_support_label(node.name)
        for clade, node in bootstrap_nodes.items()
    }
    posterior_frequency_by_clade = {
        row.clade: row.frequency for row in posterior_report.clade_frequencies
    }
    all_clades = sorted(
        set(bootstrap_support_by_clade) | set(posterior_frequency_by_clade)
    )
    rows: list[BootstrapPosteriorCladeComparison] = []
    for clade in all_clades:
        bootstrap_support = bootstrap_support_by_clade.get(clade)
        posterior_frequency = posterior_frequency_by_clade.get(clade)
        absolute_delta = None
        if bootstrap_support is not None and posterior_frequency is not None:
            absolute_delta = abs(bootstrap_support - posterior_frequency)
        agreement = _support_agreement_label(
            bootstrap_support, posterior_frequency, absolute_delta
        )
        rows.append(
            BootstrapPosteriorCladeComparison(
                clade=clade,
                bootstrap_support=bootstrap_support,
                posterior_frequency=posterior_frequency,
                absolute_delta=absolute_delta,
                agreement=agreement,
            )
        )
    topology_mismatch_clade_count = sum(
        1 for row in rows if row.agreement == "method_specific"
    )
    return BootstrapPosteriorSupportComparisonReport(
        bootstrap_tree_path=bootstrap_tree_path,
        posterior_tree_set_path=posterior_tree_set_path,
        posterior_tree_count=posterior_report.tree_count,
        shared_taxa=posterior_report.shared_taxa,
        high_conflict_clade_count=sum(
            1 for row in rows if row.agreement == "strong_conflict"
        ),
        topology_mismatch_detected=topology_mismatch_clade_count > 0,
        topology_mismatch_clade_count=topology_mismatch_clade_count,
        rows=rows,
    )


def benchmark_tree_set_uncertainty(
    *,
    tree_counts: list[int] | None = None,
    taxon_counts: list[int] | None = None,
    replicates: int = 1,
    seed: int = 1,
) -> TreeSetScalingBenchmarkReport:
    """Benchmark core posterior-uncertainty summaries across tree-count and taxon-count scaling."""
    from bijux_phylogenetics.simulation import (
        simulate_birth_death_trees,
        write_tree_set,
    )

    counts = tree_counts or [8, 32, 128]
    taxa = taxon_counts or [8, 32, 64]
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    rows: list[TreeSetBenchmarkRow] = []
    temp_root = Path(tempfile.mkdtemp(prefix="bijux-tree-set-benchmark-"))
    try:
        for taxon_count in taxa:
            for tree_count in counts:
                for replicate in range(1, replicates + 1):
                    trees, _ = simulate_birth_death_trees(
                        tree_count=tree_count,
                        tip_count=taxon_count,
                        seed=seed + len(rows),
                    )
                    tree_set_path = write_tree_set(
                        temp_root
                        / f"trees-{tree_count}-taxa-{taxon_count}-replicate-{replicate}.nwk",
                        trees,
                    )
                    started = perf_counter()
                    summary = load_tree_set(tree_set_path)
                    unstable_taxa = detect_unstable_taxa(tree_set_path)
                    unstable_clades = detect_unstable_clades(tree_set_path)
                    conclusions = summarize_uncertainty_aware_conclusions(tree_set_path)
                    elapsed = perf_counter() - started
                    rows.append(
                        TreeSetBenchmarkRow(
                            tree_count=tree_count,
                            taxon_count=taxon_count,
                            replicate=replicate,
                            elapsed_seconds=round(elapsed, 6),
                            peak_memory_bytes=max(
                                summary.processing.peak_memory_bytes,
                                unstable_taxa.processing.peak_memory_bytes,
                                unstable_clades.processing.peak_memory_bytes,
                            ),
                            rooted_topology_count=summary.rooted_topology_count,
                            unstable_taxon_count=len(unstable_taxa.taxa),
                            unstable_clade_count=len(unstable_clades.clades),
                            robust_clade_count=conclusions.robust_clade_count,
                        )
                    )
    finally:
        for path in sorted(temp_root.rglob("*"), reverse=True):
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                path.rmdir()
        temp_root.rmdir()
    return TreeSetScalingBenchmarkReport(
        tree_counts=sorted(counts),
        taxon_counts=sorted(taxa),
        rows=rows,
    )


def assess_tree_set_storage_risk(path: Path) -> TreeSetStorageRiskReport:
    """Flag large posterior/bootstrap outputs that may be awkward to store or review."""
    summary = load_tree_set(path)
    file_size_bytes = path.stat().st_size
    file_size_megabytes = round(file_size_bytes / (1024 * 1024), 6)
    mean_bytes_per_tree = round(file_size_bytes / max(summary.tree_count, 1), 6)
    warnings: list[str] = []
    risk_level = "low"
    if summary.tree_count >= 1000 or file_size_megabytes >= 10.0:
        risk_level = "high"
        warnings.append(
            "tree set is large enough to merit explicit storage and review planning"
        )
    elif summary.tree_count >= 250 or file_size_megabytes >= 1.0:
        risk_level = "moderate"
        warnings.append(
            "tree set is large enough that reviewer-facing summaries should be preferred over raw browsing"
        )
    if summary.rooted_topology_count >= 100:
        warnings.append(
            "high topology diversity increases the cost of manually inspecting individual trees"
        )
    return TreeSetStorageRiskReport(
        path=path,
        file_size_bytes=file_size_bytes,
        file_size_megabytes=file_size_megabytes,
        tree_count=summary.tree_count,
        rooted_topology_count=summary.rooted_topology_count,
        shared_taxon_count=len(summary.shared_taxa),
        mean_bytes_per_tree=mean_bytes_per_tree,
        risk_level=risk_level,
        warnings=warnings,
    )


def assess_tree_set_thinning_sensitivity(
    path: Path,
    *,
    thinning_intervals: list[int] | None = None,
) -> TreeSetThinningSensitivityReport:
    """Compare core posterior conclusions before and after deterministic thinning intervals."""
    from bijux_phylogenetics.simulation import write_tree_set

    intervals = thinning_intervals or [2, 5, 10]
    if any(interval < 1 for interval in intervals):
        raise ValueError("all thinning intervals must be at least 1")
    baseline_summary = load_tree_set(path)
    baseline_clusters = cluster_trees_by_topology(path)
    summarize_uncertainty_aware_conclusions(path)
    baseline_dominant = (
        0.0
        if not baseline_clusters.clusters
        else baseline_clusters.clusters[0].frequency
    )
    rows: list[TreeSetThinningSensitivityRow] = []
    warnings: list[str] = []
    temp_root = Path(tempfile.mkdtemp(prefix="bijux-tree-set-thinning-"))
    try:
        _, trees = _require_tree_set(path)
        for interval in sorted(set(intervals)):
            retained = trees[::interval]
            retained_path = write_tree_set(
                temp_root / f"thinned-{interval}.nwk", retained
            )
            summary = load_tree_set(retained_path)
            clusters = cluster_trees_by_topology(retained_path)
            conclusions = summarize_uncertainty_aware_conclusions(retained_path)
            comparison = compare_posterior_tree_sets(path, retained_path)
            dominant = 0.0 if not clusters.clusters else clusters.clusters[0].frequency
            row_warnings: list[str] = []
            if (
                comparison.shared_rooted_topology_count
                < baseline_summary.rooted_topology_count
            ):
                row_warnings.append(
                    "thinning drops one or more rooted topology modes observed in the full tree set"
                )
            if abs(dominant - baseline_dominant) >= 0.2:
                row_warnings.append(
                    "thinning changes the dominant topology frequency materially"
                )
            row = TreeSetThinningSensitivityRow(
                thinning_interval=interval,
                retained_tree_count=summary.tree_count,
                retained_fraction=round(
                    summary.tree_count / baseline_summary.tree_count, 15
                ),
                rooted_topology_count=summary.rooted_topology_count,
                shared_rooted_topology_count=comparison.shared_rooted_topology_count,
                dominant_topology_frequency=dominant,
                dominant_topology_delta=round(dominant - baseline_dominant, 15),
                robust_clade_count=conclusions.robust_clade_count,
                uncertain_clade_count=conclusions.uncertain_clade_count,
                conflicting_clade_count=conclusions.conflicting_clade_count,
                warnings=row_warnings,
            )
            rows.append(row)
            warnings.extend(row_warnings)
    finally:
        for file_path in sorted(temp_root.glob("*"), reverse=True):
            if file_path.is_file():
                file_path.unlink()
        temp_root.rmdir()
    return TreeSetThinningSensitivityReport(
        path=path,
        original_tree_count=baseline_summary.tree_count,
        original_rooted_topology_count=baseline_summary.rooted_topology_count,
        original_dominant_topology_frequency=baseline_dominant,
        rows=rows,
        warnings=sorted(dict.fromkeys(warnings)),
    )


def compare_consensus_thresholds(
    path: Path,
    *,
    thresholds: list[float] | None = None,
) -> ConsensusThresholdSensitivityReport:
    """Compare consensus trees across multiple posterior clade-frequency thresholds."""
    threshold_values = thresholds or [0.5, 0.75, 0.9]
    if any(not 0.0 < threshold < 1.0 for threshold in threshold_values):
        raise ValueError("all consensus thresholds must be between 0 and 1")
    summary = load_tree_set(path)
    rows: list[ConsensusThresholdSensitivityRow] = []
    warnings: list[str] = []
    for threshold in sorted(set(threshold_values)):
        tree, report = compute_consensus_tree_with_threshold(path, threshold=threshold)
        topology_id = _rooted_topology_id(tree, set(summary.shared_taxa))
        informative_clade_count = len(
            informative_rooted_clades(tree, set(summary.shared_taxa))
        )
        row_warnings: list[str] = []
        if informative_clade_count == 0:
            row_warnings.append(
                "threshold collapses all informative internal clades from the consensus summary"
            )
        rows.append(
            ConsensusThresholdSensitivityRow(
                threshold=threshold,
                informative_clade_count=informative_clade_count,
                rooted_topology_id=topology_id,
                consensus_newick=report.consensus_newick,
                warnings=row_warnings,
            )
        )
        warnings.extend(row_warnings)
    if len({row.rooted_topology_id for row in rows}) > 1:
        warnings.append("consensus topology changes across tested frequency thresholds")
    return ConsensusThresholdSensitivityReport(
        path=path,
        tree_count=summary.tree_count,
        rows=rows,
        warnings=sorted(dict.fromkeys(warnings)),
    )


def assess_tree_set_maturity(
    path: Path,
    *,
    thinning_intervals: list[int] | None = None,
    consensus_thresholds: list[float] | None = None,
) -> TreeSetMaturityGateReport:
    """Classify whether a tree-set uncertainty workflow is merely exploratory or reviewer-capable."""
    summary = load_tree_set(path)
    storage = assess_tree_set_storage_risk(path)
    thinning = assess_tree_set_thinning_sensitivity(
        path, thinning_intervals=thinning_intervals
    )
    consensus = compare_consensus_thresholds(path, thresholds=consensus_thresholds)
    conclusions = summarize_uncertainty_aware_conclusions(path)
    checks = [
        TreeSetMaturityGateCheck(
            name="shared_taxa",
            satisfied=all(
                record.taxa == summary.records[0].taxa for record in summary.records
            ),
            details="all trees in the set share the exact same taxon inventory",
        ),
        TreeSetMaturityGateCheck(
            name="storage_risk",
            satisfied=storage.risk_level != "high",
            details=f"storage risk classified as {storage.risk_level}",
        ),
        TreeSetMaturityGateCheck(
            name="thinning_stability",
            satisfied=not thinning.warnings,
            details="tested thinning intervals preserve the dominant conclusions"
            if not thinning.warnings
            else "; ".join(thinning.warnings),
        ),
        TreeSetMaturityGateCheck(
            name="consensus_stability",
            satisfied=not consensus.warnings,
            details="tested consensus thresholds preserve one topology summary"
            if not consensus.warnings
            else "; ".join(consensus.warnings),
        ),
        TreeSetMaturityGateCheck(
            name="uncertainty_summary",
            satisfied=summary.tree_count >= 2
            and (
                conclusions.robust_clade_count
                + conclusions.uncertain_clade_count
                + conclusions.conflicting_clade_count
            )
            >= 1,
            details="uncertainty-aware clade conclusions are available for reviewer-facing interpretation",
        ),
    ]
    failed = sum(1 for check in checks if not check.satisfied)
    if failed == 0:
        decision = "production_capable"
    elif failed <= 2:
        decision = "usable"
    else:
        decision = "experimental"
    warnings = sorted(
        dict.fromkeys(storage.warnings + thinning.warnings + consensus.warnings)
    )
    return TreeSetMaturityGateReport(
        path=path,
        decision=decision,
        checks=checks,
        warnings=warnings,
    )
