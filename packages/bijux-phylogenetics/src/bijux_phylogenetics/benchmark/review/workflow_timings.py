from __future__ import annotations

from pathlib import Path
import tempfile

from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta.quality import build_alignment_quality_report
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_dna_alignment,
    write_simulated_alignment,
    write_tree_set,
)
from bijux_phylogenetics.trees import compute_consensus_tree

from .._fixtures import build_balanced_tree, build_caterpillar_tree
from .._measurement import measure
from ..contracts import (
    AlignmentDiagnosticsBenchmarkReport,
    AlignmentSiteBenchmarkReport,
    TreeComparisonBenchmarkReport,
    TreeSetConsensusBenchmarkReport,
    TreeValidationBenchmarkReport,
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
    observations = []
    with tempfile.TemporaryDirectory(prefix="bijux-tree-validation-") as tmpdir:
        tmp_path = Path(tmpdir)
        for label, tip_count in classes:
            tree_path = write_newick(
                tmp_path / f"{label}.nwk", build_balanced_tree(tip_count)
            )
            observations.append(
                measure(
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
    observations = []
    with tempfile.TemporaryDirectory(prefix="bijux-tree-comparison-") as tmpdir:
        tmp_path = Path(tmpdir)
        for tip_count in counts:
            left_path = write_newick(
                tmp_path / f"compare-left-{tip_count}.nwk",
                build_balanced_tree(tip_count),
            )
            right_path = write_newick(
                tmp_path / f"compare-right-{tip_count}.nwk",
                build_caterpillar_tree(tip_count),
            )
            observations.append(
                measure(
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
    observations = []
    with tempfile.TemporaryDirectory(prefix="bijux-alignment-diagnostics-") as tmpdir:
        tmp_path = Path(tmpdir)
        for sequence_count in counts:
            tree_path = write_newick(
                tmp_path / f"alignment-tree-{sequence_count}.nwk",
                build_balanced_tree(sequence_count),
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
                measure(
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
    observations = []
    with tempfile.TemporaryDirectory(prefix="bijux-alignment-sites-") as tmpdir:
        tmp_path = Path(tmpdir)
        tree_path = write_newick(
            tmp_path / "alignment-sites-tree.nwk",
            build_balanced_tree(sequence_count),
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
                measure(
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
    observations = []
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
                measure(
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
