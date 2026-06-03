from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.distance import build_distance_method_report
from bijux_phylogenetics.engines.workflows.alignment import run_alignment_trimming
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    summarize_alignment_readiness,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.render.tree_svg import render_tree_svg
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    write_tree_set,
)
from bijux_phylogenetics.trees import (
    cluster_trees_by_topology,
    compute_consensus_tree,
    detect_unstable_clades,
    detect_unstable_taxa,
    load_tree_set,
    summarize_posterior_topology_diversity,
    summarize_uncertainty_aware_conclusions,
)

from .._fixtures import (
    LARGE_ALIGNMENT_SCALING_CLASSES,
    LARGE_TREE_SCALING_TIP_COUNTS,
    LARGE_TREE_SET_SCALING_CLASSES,
    build_balanced_tree,
    interleaved_taxa,
    write_large_alignment,
    write_named_balanced_tree,
    write_trimal_benchmark_fixture,
)
from .._measurement import (
    measure,
    measure_large_alignment_observation,
    measure_large_tree_set_observation,
)
from ..contracts import (
    BenchmarkObservation,
    LargeAlignmentScalingBenchmarkReport,
    LargeAlignmentScalingObservation,
    LargeAlignmentScalingWorkflowBenchmark,
    LargeTreeScalingBenchmarkReport,
    LargeTreeScalingWorkflowBenchmark,
    LargeTreeSetScalingBenchmarkReport,
    LargeTreeSetScalingObservation,
    LargeTreeSetScalingWorkflowBenchmark,
)


def benchmark_large_tree_scaling(
    *,
    replicates: int = 1,
    tip_counts: list[int] | None = None,
) -> LargeTreeScalingBenchmarkReport:
    """Benchmark large-tree validation, comparison, rendering, and reporting."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    counts = list(tip_counts or LARGE_TREE_SCALING_TIP_COUNTS)
    if not counts:
        raise ValueError("tip_counts must contain at least one taxon count")
    if any(count < 2 for count in counts):
        raise ValueError("tip_counts must all be at least 2")

    validation_observations: list[BenchmarkObservation] = []
    comparison_observations: list[BenchmarkObservation] = []
    rendering_observations: list[BenchmarkObservation] = []
    reporting_observations: list[BenchmarkObservation] = []

    from bijux_phylogenetics.reports.service import render_tree_report

    with tempfile.TemporaryDirectory(prefix="bijux-large-tree-scaling-") as tmpdir:
        tmp_path = Path(tmpdir)
        for tip_count in counts:
            balanced_tree_path = write_newick(
                tmp_path / f"large-tree-balanced-{tip_count}.nwk",
                build_balanced_tree(tip_count, prefix="LargeTaxon"),
            )
            comparison_tree_path = write_named_balanced_tree(
                tmp_path / f"large-tree-permuted-balanced-{tip_count}.nwk",
                interleaved_taxa(tip_count, prefix="LargeTaxon"),
            )
            render_output_path = tmp_path / f"large-tree-render-{tip_count}.svg"
            report_output_path = tmp_path / f"large-tree-report-{tip_count}.html"

            validation_observations.append(
                measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda path=balanced_tree_path: validate_tree_path(path),
                )
            )
            comparison_observations.append(
                measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda left=balanced_tree_path, right=comparison_tree_path: (
                        compare_tree_paths(left, right)
                    ),
                )
            )
            rendering_observations.append(
                measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda path=balanced_tree_path, out_path=render_output_path: (
                        render_tree_svg(path, out_path=out_path)
                    ),
                )
            )
            reporting_observations.append(
                measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda path=balanced_tree_path, out_path=report_output_path: (
                        render_tree_report(tree_path=path, out_path=out_path)
                    ),
                )
            )

    workflows = [
        LargeTreeScalingWorkflowBenchmark(
            workflow="tree-validation",
            scaling_axis="taxa",
            observations=validation_observations,
            notes=[
                "measures full structural validation on deterministic balanced trees",
            ],
        ),
        LargeTreeScalingWorkflowBenchmark(
            workflow="tree-comparison",
            scaling_axis="taxa",
            observations=comparison_observations,
            notes=[
                "compares one balanced tree against one deterministically permuted balanced tree across the same shared taxa",
            ],
        ),
        LargeTreeScalingWorkflowBenchmark(
            workflow="tree-rendering",
            scaling_axis="taxa",
            observations=rendering_observations,
            notes=[
                "renders reviewer-facing SVG output with tip labels for each governed tree size",
            ],
        ),
        LargeTreeScalingWorkflowBenchmark(
            workflow="tree-reporting",
            scaling_axis="taxa",
            observations=reporting_observations,
            notes=[
                "builds the full HTML tree report, including validation, inspection, forensic review, and machine manifest output",
            ],
        ),
    ]
    return LargeTreeScalingBenchmarkReport(
        replicates=replicates,
        tip_counts=counts,
        workflows=workflows,
        limitations=[
            "large-tree scaling numbers are local benchmark observations and should be re-run on target hardware before operational promises are made",
            "benchmarks use deterministic synthetic trees so they measure owned workflow cost without conflating external dataset quirks",
        ],
    )


def benchmark_large_alignment_scaling(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
) -> LargeAlignmentScalingBenchmarkReport:
    """Benchmark large-alignment diagnostics, trimming, distance, and readiness."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    classes = list(size_classes or LARGE_ALIGNMENT_SCALING_CLASSES)
    if not classes:
        raise ValueError("size_classes must contain at least one size class")
    if any(
        sequence_count < 2 or alignment_length < 2
        for _, sequence_count, alignment_length in classes
    ):
        raise ValueError(
            "alignment size classes must use at least two sequences and two sites"
        )

    diagnostics_observations: list[LargeAlignmentScalingObservation] = []
    trimming_observations: list[LargeAlignmentScalingObservation] = []
    distance_observations: list[LargeAlignmentScalingObservation] = []
    readiness_observations: list[LargeAlignmentScalingObservation] = []

    with tempfile.TemporaryDirectory(prefix="bijux-large-alignment-scaling-") as tmpdir:
        tmp_path = Path(tmpdir)
        trimal_executable = write_trimal_benchmark_fixture(
            tmp_path / "trimal-benchmark-fixture"
        )
        for label, sequence_count, alignment_length in classes:
            alignment_path = write_large_alignment(
                tmp_path / f"{label}.fasta",
                sequence_count=sequence_count,
                sequence_length=alignment_length,
            )
            trimmed_path = tmp_path / f"{label}.trimmed.fasta"

            diagnostics_observations.append(
                measure_large_alignment_observation(
                    label,
                    sequence_count=sequence_count,
                    alignment_length=alignment_length,
                    replicates=replicates,
                    callback=lambda path=alignment_path: build_alignment_quality_report(
                        path
                    ),
                )
            )
            trimming_observations.append(
                measure_large_alignment_observation(
                    label,
                    sequence_count=sequence_count,
                    alignment_length=alignment_length,
                    replicates=replicates,
                    callback=lambda path=alignment_path, out_path=trimmed_path: (
                        run_alignment_trimming(
                            path,
                            out_path,
                            executable=trimal_executable,
                            mode="gap-threshold",
                        )
                    ),
                )
            )
            distance_observations.append(
                measure_large_alignment_observation(
                    label,
                    sequence_count=sequence_count,
                    alignment_length=alignment_length,
                    replicates=replicates,
                    callback=lambda path=alignment_path: build_distance_method_report(
                        path,
                        model="amino-acid-p-distance",
                        bootstrap_replicates=5,
                    ),
                )
            )
            readiness_observations.append(
                measure_large_alignment_observation(
                    label,
                    sequence_count=sequence_count,
                    alignment_length=alignment_length,
                    replicates=replicates,
                    callback=lambda path=alignment_path: summarize_alignment_readiness(
                        path
                    ),
                )
            )

    workflows = [
        LargeAlignmentScalingWorkflowBenchmark(
            workflow="alignment-diagnostics",
            scaling_axis="aligned_sites",
            observations=diagnostics_observations,
            notes=[
                "measures the full owned alignment-quality report on aligned protein FASTA inputs",
            ],
        ),
        LargeAlignmentScalingWorkflowBenchmark(
            workflow="alignment-trimming",
            scaling_axis="aligned_sites",
            observations=trimming_observations,
            notes=[
                "runs the governed trimming workflow through a deterministic trimAl fixture so manifest and output validation costs are included",
            ],
        ),
        LargeAlignmentScalingWorkflowBenchmark(
            workflow="distance-analysis",
            scaling_axis="aligned_sites",
            observations=distance_observations,
            notes=[
                "builds the owned distance-method report with reduced bootstrap replicates so large-alignment scaling stays practical while still exercising matrix, tree, and maturity surfaces",
            ],
        ),
        LargeAlignmentScalingWorkflowBenchmark(
            workflow="alignment-readiness",
            scaling_axis="aligned_sites",
            observations=readiness_observations,
            notes=[
                "measures the reviewer-facing readiness summary used to decide whether large aligned inputs are suitable for downstream inference families",
            ],
        ),
    ]
    return LargeAlignmentScalingBenchmarkReport(
        replicates=replicates,
        sequence_counts=[sequence_count for _, sequence_count, _ in classes],
        alignment_lengths=[alignment_length for _, _, alignment_length in classes],
        workflows=workflows,
        limitations=[
            "large-alignment scaling numbers are local benchmark observations and should be re-run on target hardware before operational promises are made",
            "distance-analysis uses five bootstrap replicates inside the benchmark so the report path is exercised without letting bootstrap resampling dominate the scaling suite",
        ],
    )


def _benchmark_tree_set_uncertainty_summary_workflow(path: Path) -> None:
    load_tree_set(path)
    detect_unstable_taxa(path)
    detect_unstable_clades(path)
    summarize_uncertainty_aware_conclusions(path)


def benchmark_large_tree_set_scaling(
    *,
    replicates: int = 1,
    size_classes: list[tuple[str, int, int]] | None = None,
) -> LargeTreeSetScalingBenchmarkReport:
    """Benchmark large-tree-set consensus, RF diversity, clustering, and uncertainty summaries."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    classes = list(size_classes or LARGE_TREE_SET_SCALING_CLASSES)
    if not classes:
        raise ValueError("size_classes must contain at least one tree-set size class")
    if any(tree_count < 2 or tip_count < 2 for _, tree_count, tip_count in classes):
        raise ValueError(
            "tree-set size classes must use at least two trees and two taxa"
        )

    consensus_observations: list[LargeTreeSetScalingObservation] = []
    rf_observations: list[LargeTreeSetScalingObservation] = []
    clustering_observations: list[LargeTreeSetScalingObservation] = []
    uncertainty_observations: list[LargeTreeSetScalingObservation] = []

    with tempfile.TemporaryDirectory(prefix="bijux-large-tree-set-scaling-") as tmpdir:
        tmp_path = Path(tmpdir)
        for index, (label, tree_count, tip_count) in enumerate(classes):
            trees, _ = simulate_birth_death_trees(
                tree_count=tree_count,
                tip_count=tip_count,
                seed=1_000 + index,
            )
            tree_set_path = write_tree_set(tmp_path / f"{label}.trees", trees)

            consensus_observations.append(
                measure_large_tree_set_observation(
                    label,
                    tree_count=tree_count,
                    tip_count=tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_set_path: compute_consensus_tree(path),
                )
            )
            rf_observations.append(
                measure_large_tree_set_observation(
                    label,
                    tree_count=tree_count,
                    tip_count=tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_set_path: (
                        summarize_posterior_topology_diversity(path)
                    ),
                )
            )
            clustering_observations.append(
                measure_large_tree_set_observation(
                    label,
                    tree_count=tree_count,
                    tip_count=tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_set_path: cluster_trees_by_topology(path),
                )
            )
            uncertainty_observations.append(
                measure_large_tree_set_observation(
                    label,
                    tree_count=tree_count,
                    tip_count=tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_set_path: (
                        _benchmark_tree_set_uncertainty_summary_workflow(path)
                    ),
                )
            )

    workflows = [
        LargeTreeSetScalingWorkflowBenchmark(
            workflow="tree-set-consensus",
            scaling_axis="posterior_samples",
            observations=consensus_observations,
            notes=[
                "computes the owned consensus-tree summary from one simulated posterior tree set at each governed sample and taxon class",
            ],
        ),
        LargeTreeSetScalingWorkflowBenchmark(
            workflow="pairwise-rf-diversity",
            scaling_axis="posterior_samples",
            observations=rf_observations,
            notes=[
                "measures the posterior topology diversity workflow, including pairwise RF-distance aggregation across every retained tree pair",
            ],
        ),
        LargeTreeSetScalingWorkflowBenchmark(
            workflow="topology-clustering",
            scaling_axis="posterior_samples",
            observations=clustering_observations,
            notes=[
                "clusters identical rooted topologies so reviewers can see whether large posterior sets collapse into a few dominant modes",
            ],
        ),
        LargeTreeSetScalingWorkflowBenchmark(
            workflow="uncertainty-summaries",
            scaling_axis="posterior_samples",
            observations=uncertainty_observations,
            notes=[
                "runs the full uncertainty-summary path, including unstable taxa, unstable clades, and reviewer-facing conclusion summaries",
            ],
        ),
    ]
    return LargeTreeSetScalingBenchmarkReport(
        replicates=replicates,
        tree_counts=[tree_count for _, tree_count, _ in classes],
        tip_counts=[tip_count for _, _, tip_count in classes],
        workflows=workflows,
        limitations=[
            "large-tree-set scaling numbers are local benchmark observations and should be re-run on target hardware before operational promises are made",
            "tree-set classes increase posterior sample count and taxon count together so consensus, clustering, RF diversity, and uncertainty summaries are measured across reviewer-relevant large posterior workloads",
        ],
    )
