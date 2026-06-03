from __future__ import annotations

import csv
from pathlib import Path

from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clade_nodes

from ..tree_sets.budgets import (
    build_tree_set_budget_report,
    build_tree_set_workflow_budget,
    enforce_tree_set_tree_budget,
)
from ..tree_sets.clade_support import (
    _build_clade_frequency_report,
    _support_classification,
    write_clade_frequency_table,
)
from ..tree_sets.consensus import (
    _build_consensus_tree_with_threshold,
    write_consensus_tree,
)
from ..tree_sets.contracts import TreeSetReport
from ..tree_sets.distances import (
    _build_tree_distance_matrix_report,
    write_tree_distance_matrix,
)
from ..tree_sets.inventory import (
    _analyze_tree_set,
    _TreeSetAnalysis,
)
from ..tree_sets.topology import _format_clade
from .instability import _build_unstable_clade_report, write_unstable_clade_table
from .models import (
    BootstrapTreeSetArtifactReport,
    BootstrapTreeSetSummaryReport,
    BootstrapUnstableBranch,
)
from .topology_diversity import (
    _build_posterior_topology_diversity_report,
    _build_topology_cluster_report,
    write_topology_cluster_table,
    write_tree_distance_distribution_table,
)


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
            frequency.frequency,
            conflict_count,
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
    path: Path,
    report: BootstrapTreeSetSummaryReport,
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
                    report.diversity.dominant_topology_frequency,
                    ".15g",
                ),
                "effective_topology_count": format(
                    report.diversity.effective_topology_count,
                    ".15g",
                ),
                "mean_robinson_foulds_distance": format(
                    report.diversity.mean_robinson_foulds_distance,
                    ".15g",
                ),
                "mean_normalized_robinson_foulds_distance": format(
                    report.diversity.mean_normalized_robinson_foulds_distance,
                    ".15g",
                ),
                "consensus_threshold": format(report.consensus_threshold, ".15g"),
                "robust_support_threshold": format(
                    report.robust_support_threshold,
                    ".15g",
                ),
                "unstable_branch_count": report.unstable_branch_count,
                "warning_count": len(report.warnings),
                "consensus_newick": report.consensus.consensus_newick,
            }
        )
    return path


def write_bootstrap_unstable_branch_table(
    path: Path,
    report: BootstrapTreeSetSummaryReport,
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
                        row.bootstrap_support_percent,
                        ".15g",
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
            base_path.with_suffix(".summary.tsv"),
            summary_report,
        ),
        "consensus_tree": write_consensus_tree(
            base_path.with_suffix(".consensus.nwk"),
            consensus_tree,
        ),
        "clade_frequencies": write_clade_frequency_table(
            base_path.with_suffix(".clade-frequencies.tsv"),
            summary_report.clade_frequencies,
        ),
        "unstable_branches": write_bootstrap_unstable_branch_table(
            base_path.with_suffix(".unstable-branches.tsv"),
            summary_report,
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
