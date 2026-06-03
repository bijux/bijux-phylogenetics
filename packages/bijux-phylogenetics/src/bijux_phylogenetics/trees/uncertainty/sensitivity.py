from __future__ import annotations

from pathlib import Path
import tempfile
from time import perf_counter

from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades

from ..tree_sets import (
    compute_consensus_tree_with_threshold,
    load_tree_set,
)
from ..tree_sets.inventory import (
    _require_tree_set,
)
from ..tree_sets.topology import _rooted_topology_id
from .comparisons import compare_posterior_tree_sets
from .instability import (
    detect_unstable_clades,
    detect_unstable_taxa,
    summarize_uncertainty_aware_conclusions,
)
from .models import (
    ConsensusThresholdSensitivityReport,
    ConsensusThresholdSensitivityRow,
    TreeSetBenchmarkRow,
    TreeSetMaturityGateCheck,
    TreeSetMaturityGateReport,
    TreeSetScalingBenchmarkReport,
    TreeSetStorageRiskReport,
    TreeSetThinningSensitivityReport,
    TreeSetThinningSensitivityRow,
)
from .topology_diversity import cluster_trees_by_topology


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
                temp_root / f"thinned-{interval}.nwk",
                retained,
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
                    summary.tree_count / baseline_summary.tree_count,
                    15,
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
        path,
        thinning_intervals=thinning_intervals,
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
