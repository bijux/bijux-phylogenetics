from __future__ import annotations

import tempfile
import time
import tracemalloc
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.newick import write_newick


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
