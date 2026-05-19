from __future__ import annotations

from pathlib import Path
import tempfile
import time
import tracemalloc

from .fixtures import (
    LARGE_ALIGNMENT_SCALING_CLASSES as LARGE_ALIGNMENT_SCALING_CLASSES,
    LARGE_TREE_SCALING_TIP_COUNTS as LARGE_TREE_SCALING_TIP_COUNTS,
    LARGE_TREE_SET_SCALING_CLASSES as LARGE_TREE_SET_SCALING_CLASSES,
    build_balanced_tree as _build_balanced_tree,
    build_caterpillar_tree as _build_caterpillar_tree,
    comparative_stress_payload as _comparative_stress_payload,
    interleaved_taxa as _interleaved_taxa,
    large_alignment_stress_payload as _large_alignment_stress_payload,
    resolve_stress_tier_config as _resolve_stress_tier_config,
    supermatrix_stress_payload as _supermatrix_stress_payload,
    table_generation_stress_payload as _table_generation_stress_payload,
    tree_set_stress_payload as _tree_set_stress_payload,
    write_large_alignment as _write_large_alignment,
    write_named_balanced_tree as _write_named_balanced_tree,
    write_trimal_benchmark_fixture as _write_trimal_benchmark_fixture,
)
from .models import (
    AlignmentDiagnosticsBenchmarkReport as AlignmentDiagnosticsBenchmarkReport,
    AlignmentSiteBenchmarkReport as AlignmentSiteBenchmarkReport,
    BenchmarkObservation as BenchmarkObservation,
    LargeAlignmentScalingBenchmarkReport as LargeAlignmentScalingBenchmarkReport,
    LargeAlignmentScalingObservation as LargeAlignmentScalingObservation,
    LargeAlignmentScalingWorkflowBenchmark as LargeAlignmentScalingWorkflowBenchmark,
    LargeDatasetStressObservation as LargeDatasetStressObservation,
    LargeDatasetStressSuiteReport as LargeDatasetStressSuiteReport,
    LargeTreeScalingBenchmarkReport as LargeTreeScalingBenchmarkReport,
    LargeTreeScalingWorkflowBenchmark as LargeTreeScalingWorkflowBenchmark,
    LargeTreeSetScalingBenchmarkReport as LargeTreeSetScalingBenchmarkReport,
    LargeTreeSetScalingObservation as LargeTreeSetScalingObservation,
    LargeTreeSetScalingWorkflowBenchmark as LargeTreeSetScalingWorkflowBenchmark,
    TreeComparisonBenchmarkReport as TreeComparisonBenchmarkReport,
    TreeSetConsensusBenchmarkReport as TreeSetConsensusBenchmarkReport,
    TreeValidationBenchmarkReport as TreeValidationBenchmarkReport,
    WorkflowPracticalLimitEntry as WorkflowPracticalLimitEntry,
    WorkflowPracticalLimitReport as WorkflowPracticalLimitReport,
    _StressObservationPayload as _StressObservationPayload,
    _StressTierConfig as _StressTierConfig,
)
from .model_fitting import (
    LargeTreeModelFittingBenchmarkBundle as LargeTreeModelFittingBenchmarkBundle,
    LargeTreeModelFittingBenchmarkReport as LargeTreeModelFittingBenchmarkReport,
    LargeTreeModelFittingObservation as LargeTreeModelFittingObservation,
    LargeTreeModelFittingThreshold as LargeTreeModelFittingThreshold,
    benchmark_large_tree_model_fitting as benchmark_large_tree_model_fitting,
    write_large_tree_model_fitting_bundle as write_large_tree_model_fitting_bundle,
    write_large_tree_model_fitting_observation_table as write_large_tree_model_fitting_observation_table,
    write_large_tree_model_fitting_summary_table as write_large_tree_model_fitting_summary_table,
)
from .real_dataset_macroevolution import (
    RealDatasetMacroevolutionAlignmentReviewRow as RealDatasetMacroevolutionAlignmentReviewRow,
    RealDatasetMacroevolutionBenchmarkBundle as RealDatasetMacroevolutionBenchmarkBundle,
    RealDatasetMacroevolutionBenchmarkDemoResult as RealDatasetMacroevolutionBenchmarkDemoResult,
    RealDatasetMacroevolutionBenchmarkReport as RealDatasetMacroevolutionBenchmarkReport,
    RealDatasetMacroevolutionModelRow as RealDatasetMacroevolutionModelRow,
    RealDatasetMacroevolutionParityRow as RealDatasetMacroevolutionParityRow,
    RealDatasetMacroevolutionSummaryRow as RealDatasetMacroevolutionSummaryRow,
    benchmark_real_dataset_macroevolution as benchmark_real_dataset_macroevolution,
    run_real_dataset_macroevolution_benchmark_demo as run_real_dataset_macroevolution_benchmark_demo,
    write_geiger_real_dataset_reference_payload_table as write_geiger_real_dataset_reference_payload_table,
    write_real_dataset_macroevolution_alignment_review_table as write_real_dataset_macroevolution_alignment_review_table,
    write_real_dataset_macroevolution_bundle as write_real_dataset_macroevolution_bundle,
    write_real_dataset_macroevolution_model_table as write_real_dataset_macroevolution_model_table,
    write_real_dataset_macroevolution_parity_table as write_real_dataset_macroevolution_parity_table,
    write_real_dataset_macroevolution_summary_table as write_real_dataset_macroevolution_summary_table,
)
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.distance import build_distance_method_report
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_quality_report,
    summarize_alignment_readiness,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.render.svg import render_tree_svg
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_dna_alignment,
    write_simulated_alignment,
    write_tree_set,
)
from bijux_phylogenetics.trees import cluster_trees_by_topology
from bijux_phylogenetics.trees import compute_consensus_tree
from bijux_phylogenetics.trees import detect_unstable_clades
from bijux_phylogenetics.trees import detect_unstable_taxa
from bijux_phylogenetics.trees import load_tree_set
from bijux_phylogenetics.trees import summarize_posterior_topology_diversity
from bijux_phylogenetics.trees import summarize_uncertainty_aware_conclusions
from bijux_phylogenetics.engines.workflows.alignment import run_alignment_trimming



def _max_runtime_seconds(observations) -> float:
    return round(max(row.runtime_seconds for row in observations), 15)


def _max_peak_memory_bytes(observations) -> int:
    return max(row.peak_memory_bytes for row in observations)


def _measure(
    label: str, item_count: int, *, replicates: int, callback
) -> BenchmarkObservation:
    runtimes: list[float] = []
    peak_memory = 0
    for _ in range(replicates):
        tracemalloc.start()
        started = time.perf_counter()
        callback()
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        runtimes.append(elapsed)
        peak_memory = max(peak_memory, peak)
    return BenchmarkObservation(
        label=label,
        item_count=item_count,
        runtime_seconds=round(sum(runtimes) / len(runtimes), 15),
        peak_memory_bytes=peak_memory,
    )


def _measure_large_alignment_observation(
    label: str,
    *,
    sequence_count: int,
    alignment_length: int,
    replicates: int,
    callback,
) -> LargeAlignmentScalingObservation:
    runtimes: list[float] = []
    peak_memory = 0
    aligned_site_count = sequence_count * alignment_length
    for _ in range(replicates):
        tracemalloc.start()
        started = time.perf_counter()
        callback()
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        runtimes.append(elapsed)
        peak_memory = max(peak_memory, peak)
    return LargeAlignmentScalingObservation(
        label=label,
        sequence_count=sequence_count,
        alignment_length=alignment_length,
        aligned_site_count=aligned_site_count,
        runtime_seconds=round(sum(runtimes) / len(runtimes), 15),
        peak_memory_bytes=peak_memory,
    )


def _measure_large_tree_set_observation(
    label: str,
    *,
    tree_count: int,
    tip_count: int,
    replicates: int,
    callback,
) -> LargeTreeSetScalingObservation:
    runtimes: list[float] = []
    peak_memory = 0
    pair_count = tree_count * max(tree_count - 1, 0) // 2
    for _ in range(replicates):
        tracemalloc.start()
        started = time.perf_counter()
        callback()
        elapsed = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        runtimes.append(elapsed)
        peak_memory = max(peak_memory, peak)
    return LargeTreeSetScalingObservation(
        label=label,
        tree_count=tree_count,
        tip_count=tip_count,
        pair_count=pair_count,
        runtime_seconds=round(sum(runtimes) / len(runtimes), 15),
        peak_memory_bytes=peak_memory,
    )



def _summarize_memory_observation_kinds(kinds: list[str]) -> str:
    distinct = list(dict.fromkeys(kind for kind in kinds if kind))
    if not distinct:
        return "python-tracemalloc"
    if len(distinct) == 1:
        return distinct[0]
    return "mixed"


def _measure_stress_workload(
    callback,
) -> tuple[_StressObservationPayload, float, int, str]:
    tracemalloc.start()
    started = time.perf_counter()
    payload, observed_peak_memory_bytes, observed_memory_kind = callback()
    elapsed_seconds = time.perf_counter() - started
    _, tracemalloc_peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_memory_bytes = max(
        tracemalloc_peak_bytes,
        0 if observed_peak_memory_bytes is None else observed_peak_memory_bytes,
    )
    return (
        payload,
        elapsed_seconds,
        peak_memory_bytes,
        observed_memory_kind or "python-tracemalloc",
    )


def benchmark_large_dataset_stress_suite(
    *,
    tier: str = "small",
) -> LargeDatasetStressSuiteReport:
    """Benchmark large owned workloads across one governed stress tier."""
    config = _resolve_stress_tier_config(tier)
    observations: list[LargeDatasetStressObservation] = []
    limitations = [
        "resource peaks are measured with python tracemalloc where possible and reuse stage-level engine memory observations when an owned workflow already records them",
        "timeout_seconds is a workload budget recorded for review; only engine-backed workflows enforce it internally during execution",
    ]
    with tempfile.TemporaryDirectory(prefix=f"bijux-stress-{config.tier}-") as tmpdir:
        root = Path(tmpdir)
        workloads = [
            lambda: _large_alignment_stress_payload(
                root=root / "alignment", config=config
            ),
            lambda: _supermatrix_stress_payload(
                root=root / "supermatrix", config=config
            ),
            lambda: _tree_set_stress_payload(root=root / "tree-set", config=config),
            lambda: _comparative_stress_payload(
                root=root / "comparative", config=config
            ),
            lambda: _table_generation_stress_payload(
                root=root / "tables", config=config
            ),
        ]
        for workload in workloads:
            payload, runtime_seconds, peak_memory_bytes, memory_observation_kind = (
                _measure_stress_workload(workload)
            )
            observations.append(
                LargeDatasetStressObservation(
                    workload=payload.workload,
                    tier=config.tier,
                    timeout_seconds=config.timeout_seconds,
                    input_size_bytes=payload.input_size_bytes,
                    sequence_count=payload.sequence_count,
                    alignment_length=payload.alignment_length,
                    tree_count=payload.tree_count,
                    taxon_count=payload.taxon_count,
                    locus_count=payload.locus_count,
                    runtime_seconds=round(runtime_seconds, 15),
                    peak_memory_bytes=peak_memory_bytes,
                    memory_observation_kind=memory_observation_kind,
                    output_row_count=payload.output_row_count,
                    notes=payload.notes,
                )
            )
    return LargeDatasetStressSuiteReport(
        tier=config.tier,
        observations=observations,
        limitations=limitations,
    )


def benchmark_tree_validation(
    *,
    replicates: int = 3,
    size_classes: list[tuple[str, int]] | None = None,
) -> TreeValidationBenchmarkReport:
    """Benchmark tree validation across named size classes."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    classes = size_classes or [("small", 16), ("medium", 64), ("large", 256)]
    observations: list[BenchmarkObservation] = []
    with tempfile.TemporaryDirectory(prefix="bijux-tree-validation-") as tmpdir:
        tmp_path = Path(tmpdir)
        for label, tip_count in classes:
            tree_path = write_newick(
                tmp_path / f"{label}.nwk", _build_balanced_tree(tip_count)
            )
            observations.append(
                _measure(
                    label,
                    tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_path: validate_tree_path(path),
                )
            )
    return TreeValidationBenchmarkReport(
        replicates=replicates, observations=observations
    )


def benchmark_tree_comparison(
    *,
    replicates: int = 3,
    taxon_counts: list[int] | None = None,
) -> TreeComparisonBenchmarkReport:
    """Benchmark shared-taxon tree comparison across increasing taxon counts."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    counts = taxon_counts or [8, 16, 32, 64, 128]
    observations: list[BenchmarkObservation] = []
    with tempfile.TemporaryDirectory(prefix="bijux-tree-comparison-") as tmpdir:
        tmp_path = Path(tmpdir)
        for tip_count in counts:
            left_path = write_newick(
                tmp_path / f"compare-left-{tip_count}.nwk",
                _build_balanced_tree(tip_count),
            )
            right_path = write_newick(
                tmp_path / f"compare-right-{tip_count}.nwk",
                _build_caterpillar_tree(tip_count),
            )
            observations.append(
                _measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda left=left_path, right=right_path: (
                        compare_tree_paths(left, right)
                    ),
                )
            )
    return TreeComparisonBenchmarkReport(
        replicates=replicates, observations=observations
    )


def benchmark_alignment_diagnostics(
    *,
    replicates: int = 3,
    sequence_counts: list[int] | None = None,
    sequence_length: int = 128,
) -> AlignmentDiagnosticsBenchmarkReport:
    """Benchmark alignment-quality diagnostics across increasing sequence counts."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    counts = sequence_counts or [8, 16, 32, 64, 128]
    observations: list[BenchmarkObservation] = []
    with tempfile.TemporaryDirectory(prefix="bijux-alignment-diagnostics-") as tmpdir:
        tmp_path = Path(tmpdir)
        for sequence_count in counts:
            tree_path = write_newick(
                tmp_path / f"alignment-tree-{sequence_count}.nwk",
                _build_balanced_tree(sequence_count),
            )
            alignment_report = simulate_dna_alignment(
                tree_path,
                sequence_length=sequence_length,
                substitution_rate=1.0,
                seed=sequence_count,
            )
            alignment_path = write_simulated_alignment(
                tmp_path / f"alignment-{sequence_count}.fasta",
                alignment_report,
            )
            observations.append(
                _measure(
                    f"sequences-{sequence_count}",
                    sequence_count,
                    replicates=replicates,
                    callback=lambda path=alignment_path: build_alignment_quality_report(
                        path
                    ),
                )
            )
    return AlignmentDiagnosticsBenchmarkReport(
        replicates=replicates, observations=observations
    )


def benchmark_alignment_site_scaling(
    *,
    replicates: int = 3,
    site_counts: list[int] | None = None,
    sequence_count: int = 16,
) -> AlignmentSiteBenchmarkReport:
    """Benchmark alignment diagnostics as alignment length increases."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    counts = site_counts or [64, 128, 256, 512]
    observations: list[BenchmarkObservation] = []
    with tempfile.TemporaryDirectory(prefix="bijux-alignment-sites-") as tmpdir:
        tmp_path = Path(tmpdir)
        tree_path = write_newick(
            tmp_path / "alignment-sites-tree.nwk",
            _build_balanced_tree(sequence_count),
        )
        for site_count in counts:
            alignment_report = simulate_dna_alignment(
                tree_path,
                sequence_length=site_count,
                substitution_rate=1.0,
                seed=site_count,
            )
            alignment_path = write_simulated_alignment(
                tmp_path / f"alignment-sites-{site_count}.fasta",
                alignment_report,
            )
            observations.append(
                _measure(
                    f"sites-{site_count}",
                    site_count,
                    replicates=replicates,
                    callback=lambda path=alignment_path: build_alignment_quality_report(
                        path
                    ),
                )
            )
    return AlignmentSiteBenchmarkReport(
        replicates=replicates,
        sequence_count=sequence_count,
        observations=observations,
    )


def benchmark_tree_set_consensus(
    *,
    replicates: int = 3,
    tree_counts: list[int] | None = None,
    tip_count: int = 16,
) -> TreeSetConsensusBenchmarkReport:
    """Benchmark consensus-tree computation as posterior/bootstrap sample counts grow."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    counts = tree_counts or [8, 32, 128, 256]
    observations: list[BenchmarkObservation] = []
    with tempfile.TemporaryDirectory(prefix="bijux-tree-set-consensus-") as tmpdir:
        tmp_path = Path(tmpdir)
        for tree_count in counts:
            trees, _ = simulate_birth_death_trees(
                tree_count=tree_count,
                tip_count=tip_count,
                seed=tree_count,
            )
            tree_set_path = write_tree_set(
                tmp_path / f"tree-set-{tree_count}.trees", trees
            )
            observations.append(
                _measure(
                    f"trees-{tree_count}",
                    tree_count,
                    replicates=replicates,
                    callback=lambda path=tree_set_path: compute_consensus_tree(path),
                )
            )
    return TreeSetConsensusBenchmarkReport(
        replicates=replicates,
        tip_count=tip_count,
        observations=observations,
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
                _build_balanced_tree(tip_count, prefix="LargeTaxon"),
            )
            comparison_tree_path = _write_named_balanced_tree(
                tmp_path / f"large-tree-permuted-balanced-{tip_count}.nwk",
                _interleaved_taxa(tip_count, prefix="LargeTaxon"),
            )
            render_output_path = tmp_path / f"large-tree-render-{tip_count}.svg"
            report_output_path = tmp_path / f"large-tree-report-{tip_count}.html"

            validation_observations.append(
                _measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda path=balanced_tree_path: validate_tree_path(path),
                )
            )
            comparison_observations.append(
                _measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda left=balanced_tree_path, right=comparison_tree_path: (
                        compare_tree_paths(left, right)
                    ),
                )
            )
            rendering_observations.append(
                _measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda path=balanced_tree_path, out_path=render_output_path: (
                        render_tree_svg(path, out_path=out_path)
                    ),
                )
            )
            reporting_observations.append(
                _measure(
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
    if any(sequence_count < 2 or alignment_length < 2 for _, sequence_count, alignment_length in classes):
        raise ValueError("alignment size classes must use at least two sequences and two sites")

    diagnostics_observations: list[LargeAlignmentScalingObservation] = []
    trimming_observations: list[LargeAlignmentScalingObservation] = []
    distance_observations: list[LargeAlignmentScalingObservation] = []
    readiness_observations: list[LargeAlignmentScalingObservation] = []

    with tempfile.TemporaryDirectory(prefix="bijux-large-alignment-scaling-") as tmpdir:
        tmp_path = Path(tmpdir)
        trimal_executable = _write_trimal_benchmark_fixture(
            tmp_path / "trimal-benchmark-fixture"
        )
        for label, sequence_count, alignment_length in classes:
            alignment_path = _write_large_alignment(
                tmp_path / f"{label}.fasta",
                sequence_count=sequence_count,
                sequence_length=alignment_length,
            )
            trimmed_path = tmp_path / f"{label}.trimmed.fasta"

            diagnostics_observations.append(
                _measure_large_alignment_observation(
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
                _measure_large_alignment_observation(
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
                _measure_large_alignment_observation(
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
                _measure_large_alignment_observation(
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
    if any(
        tree_count < 2 or tip_count < 2
        for _, tree_count, tip_count in classes
    ):
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
                _measure_large_tree_set_observation(
                    label,
                    tree_count=tree_count,
                    tip_count=tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_set_path: compute_consensus_tree(path),
                )
            )
            rf_observations.append(
                _measure_large_tree_set_observation(
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
                _measure_large_tree_set_observation(
                    label,
                    tree_count=tree_count,
                    tip_count=tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_set_path: cluster_trees_by_topology(path),
                )
            )
            uncertainty_observations.append(
                _measure_large_tree_set_observation(
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


def _large_tree_limit_entries(
    report: LargeTreeScalingBenchmarkReport,
) -> list[WorkflowPracticalLimitEntry]:
    entries: list[WorkflowPracticalLimitEntry] = []
    for workflow in report.workflows:
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workflow.workflow,
                evidence_source="large-tree-scaling",
                tested_taxon_limit=max(report.tip_counts),
                tested_site_limit=None,
                tested_tree_limit=2 if workflow.workflow == "tree-comparison" else 1,
                tested_posterior_size=None,
                max_runtime_seconds=_max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=_max_peak_memory_bytes(workflow.observations),
                memory_observation_kind="python-tracemalloc",
                notes=workflow.notes,
            )
        )
    return entries


def _large_alignment_limit_entries(
    report: LargeAlignmentScalingBenchmarkReport,
) -> list[WorkflowPracticalLimitEntry]:
    entries: list[WorkflowPracticalLimitEntry] = []
    for workflow in report.workflows:
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workflow.workflow,
                evidence_source="large-alignment-scaling",
                tested_taxon_limit=max(report.sequence_counts),
                tested_site_limit=max(report.alignment_lengths),
                tested_tree_limit=None,
                tested_posterior_size=None,
                max_runtime_seconds=_max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=_max_peak_memory_bytes(workflow.observations),
                memory_observation_kind="python-tracemalloc",
                notes=workflow.notes,
            )
        )
    return entries


def _large_tree_set_limit_entries(
    report: LargeTreeSetScalingBenchmarkReport,
) -> list[WorkflowPracticalLimitEntry]:
    entries: list[WorkflowPracticalLimitEntry] = []
    for workflow in report.workflows:
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workflow.workflow,
                evidence_source="large-tree-set-scaling",
                tested_taxon_limit=max(report.tip_counts),
                tested_site_limit=None,
                tested_tree_limit=max(report.tree_counts),
                tested_posterior_size=max(report.tree_counts),
                max_runtime_seconds=_max_runtime_seconds(workflow.observations),
                max_peak_memory_bytes=_max_peak_memory_bytes(workflow.observations),
                memory_observation_kind="python-tracemalloc",
                notes=workflow.notes,
            )
        )
    return entries


def _stress_limit_entries(
    reports: list[LargeDatasetStressSuiteReport],
) -> list[WorkflowPracticalLimitEntry]:
    rows_by_workload: dict[str, list[LargeDatasetStressObservation]] = {}
    for report in reports:
        for row in report.observations:
            rows_by_workload.setdefault(row.workload, []).append(row)

    entries: list[WorkflowPracticalLimitEntry] = []
    for workload in sorted(rows_by_workload):
        rows = rows_by_workload[workload]
        tiers = sorted({row.tier for row in rows})
        max_runtime = round(max(row.runtime_seconds for row in rows), 15)
        max_peak_memory = max(row.peak_memory_bytes for row in rows)
        memory_kinds = sorted(
            {row.memory_observation_kind for row in rows if row.memory_observation_kind}
        )
        notes: list[str] = []
        for row in rows:
            for note in row.notes:
                if note not in notes:
                    notes.append(note)
        notes.append(
            "tested through governed stress tiers: " + ", ".join(tiers)
        )
        entries.append(
            WorkflowPracticalLimitEntry(
                workflow=workload,
                evidence_source="stress-suite",
                tested_taxon_limit=max(
                    row.taxon_count for row in rows if row.taxon_count is not None
                )
                if any(row.taxon_count is not None for row in rows)
                else None,
                tested_site_limit=max(
                    row.alignment_length
                    for row in rows
                    if row.alignment_length is not None
                )
                if any(row.alignment_length is not None for row in rows)
                else None,
                tested_tree_limit=max(
                    row.tree_count for row in rows if row.tree_count is not None
                )
                if any(row.tree_count is not None for row in rows)
                else None,
                tested_posterior_size=max(
                    row.tree_count
                    for row in rows
                    if row.workload == "posterior-tree-set-consensus"
                    and row.tree_count is not None
                )
                if any(
                    row.workload == "posterior-tree-set-consensus"
                    and row.tree_count is not None
                    for row in rows
                )
                else None,
                max_runtime_seconds=max_runtime,
                max_peak_memory_bytes=max_peak_memory,
                memory_observation_kind=(
                    None if not memory_kinds else ",".join(memory_kinds)
                ),
                notes=notes,
            )
        )
    return entries


def benchmark_workflow_practical_limits(
    *,
    replicates: int = 1,
    tree_tip_counts: list[int] | None = None,
    alignment_size_classes: list[tuple[str, int, int]] | None = None,
    tree_set_size_classes: list[tuple[str, int, int]] | None = None,
    stress_tiers: list[str] | None = None,
) -> WorkflowPracticalLimitReport:
    """Report the largest governed workflow classes currently exercised in benchmark and stress lanes."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    tiers = ["heavy"] if stress_tiers is None else list(stress_tiers)
    if not tiers:
        raise ValueError("stress_tiers must contain at least one governed tier")

    tree_report = benchmark_large_tree_scaling(
        replicates=replicates,
        tip_counts=tree_tip_counts,
    )
    alignment_report = benchmark_large_alignment_scaling(
        replicates=replicates,
        size_classes=alignment_size_classes,
    )
    tree_set_report = benchmark_large_tree_set_scaling(
        replicates=replicates,
        size_classes=tree_set_size_classes,
    )
    stress_reports = [
        benchmark_large_dataset_stress_suite(tier=tier) for tier in tiers
    ]

    entries = [
        *_large_tree_limit_entries(tree_report),
        *_large_alignment_limit_entries(alignment_report),
        *_large_tree_set_limit_entries(tree_set_report),
        *_stress_limit_entries(stress_reports),
    ]
    limitations: list[str] = []
    for report in (tree_report, alignment_report, tree_set_report):
        for item in report.limitations:
            if item not in limitations:
                limitations.append(item)
    for report in stress_reports:
        for item in report.limitations:
            if item not in limitations:
                limitations.append(item)
    limitations.append(
        "practical limits report tested maxima from governed benchmark and stress lanes; it does not claim untested workflows or hardware-specific guarantees"
    )
    return WorkflowPracticalLimitReport(
        replicates=replicates,
        stress_tiers=tiers,
        entries=entries,
        limitations=limitations,
    )
