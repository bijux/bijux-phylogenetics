from __future__ import annotations

import tempfile
import time
import tracemalloc
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.fasta import build_alignment_quality_report
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.simulation import simulate_dna_alignment, write_simulated_alignment


@dataclass(frozen=True, slots=True)
class BenchmarkObservation:
    label: str
    item_count: int
    runtime_seconds: float
    peak_memory_bytes: int


@dataclass(slots=True)
class TreeValidationBenchmarkReport:
    replicates: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class TreeComparisonBenchmarkReport:
    replicates: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class AlignmentDiagnosticsBenchmarkReport:
    replicates: int
    observations: list[BenchmarkObservation]


def _measure(label: str, item_count: int, *, replicates: int, callback) -> BenchmarkObservation:
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


def _build_balanced_tree(tip_count: int, *, branch_length: float = 0.1, prefix: str = "Taxon") -> PhyloTree:
    if tip_count < 2:
        raise ValueError(f"tip_count must be at least 2, got {tip_count}")
    leaves = [TreeNode(name=f"{prefix}{index}", branch_length=branch_length) for index in range(1, tip_count + 1)]
    while len(leaves) > 1:
        next_level: list[TreeNode] = []
        for index in range(0, len(leaves), 2):
            left = leaves[index]
            right = leaves[index + 1] if index + 1 < len(leaves) else None
            if right is None:
                left.branch_length = round((left.branch_length or 0.0) + branch_length, 15)
                next_level.append(left)
                continue
            next_level.append(TreeNode(children=[left, right], branch_length=branch_length))
        leaves = next_level
    root = leaves[0]
    root.branch_length = None
    return PhyloTree(root=root, source_format="newick")


def _build_caterpillar_tree(tip_count: int, *, branch_length: float = 0.1, prefix: str = "Taxon") -> PhyloTree:
    if tip_count < 2:
        raise ValueError(f"tip_count must be at least 2, got {tip_count}")
    root = TreeNode(children=[TreeNode(name=f"{prefix}1", branch_length=branch_length), TreeNode(name=f"{prefix}2", branch_length=branch_length)])
    current = root
    for index in range(3, tip_count + 1):
        new_internal = TreeNode(branch_length=branch_length, children=[current.children.pop(), TreeNode(name=f"{prefix}{index}", branch_length=branch_length)])
        current.children.append(new_internal)
        current = new_internal
    root.branch_length = None
    return PhyloTree(root=root, source_format="newick")


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
            tree_path = write_newick(tmp_path / f"{label}.nwk", _build_balanced_tree(tip_count))
            observations.append(
                _measure(
                    label,
                    tip_count,
                    replicates=replicates,
                    callback=lambda path=tree_path: validate_tree_path(path),
                )
            )
    return TreeValidationBenchmarkReport(replicates=replicates, observations=observations)


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
            left_path = write_newick(tmp_path / f"compare-left-{tip_count}.nwk", _build_balanced_tree(tip_count))
            right_path = write_newick(tmp_path / f"compare-right-{tip_count}.nwk", _build_caterpillar_tree(tip_count))
            observations.append(
                _measure(
                    f"taxa-{tip_count}",
                    tip_count,
                    replicates=replicates,
                    callback=lambda left=left_path, right=right_path: compare_tree_paths(left, right),
                )
            )
    return TreeComparisonBenchmarkReport(replicates=replicates, observations=observations)


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
                    callback=lambda path=alignment_path: build_alignment_quality_report(path),
                )
            )
    return AlignmentDiagnosticsBenchmarkReport(replicates=replicates, observations=observations)
