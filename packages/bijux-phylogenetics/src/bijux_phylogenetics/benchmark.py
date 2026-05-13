from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile
import time
import tracemalloc

from bijux_phylogenetics.branch_lengths import analyze_branch_length_distribution
from bijux_phylogenetics.clades import extract_tree_clades
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.core.concatenation import concatenate_locus_alignments
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.engines.large_alignment_inference import (
    run_large_alignment_inference,
)
from bijux_phylogenetics.io.fasta import build_alignment_quality_report
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_dna_alignment,
    write_simulated_alignment,
    write_tree_set,
)
from bijux_phylogenetics.tree_set import compute_consensus_tree


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


@dataclass(slots=True)
class AlignmentSiteBenchmarkReport:
    replicates: int
    sequence_count: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class TreeSetConsensusBenchmarkReport:
    replicates: int
    tip_count: int
    observations: list[BenchmarkObservation]


@dataclass(frozen=True, slots=True)
class LargeDatasetStressObservation:
    workload: str
    tier: str
    timeout_seconds: float
    input_size_bytes: int
    sequence_count: int | None
    alignment_length: int | None
    tree_count: int | None
    taxon_count: int | None
    locus_count: int | None
    runtime_seconds: float
    peak_memory_bytes: int
    memory_observation_kind: str
    output_row_count: int
    notes: list[str]


@dataclass(slots=True)
class LargeDatasetStressSuiteReport:
    tier: str
    observations: list[LargeDatasetStressObservation]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class _StressObservationPayload:
    workload: str
    input_size_bytes: int
    sequence_count: int | None
    alignment_length: int | None
    tree_count: int | None
    taxon_count: int | None
    locus_count: int | None
    output_row_count: int
    notes: list[str]


@dataclass(frozen=True, slots=True)
class _StressTierConfig:
    tier: str
    timeout_seconds: float
    alignment_sequence_count: int
    alignment_length: int
    supermatrix_taxon_count: int
    supermatrix_locus_lengths: tuple[int, ...]
    tree_set_tree_count: int
    tree_set_tip_count: int
    comparative_taxon_count: int
    table_tip_count: int


_STRESS_TIER_CONFIGS: dict[str, _StressTierConfig] = {
    "small": _StressTierConfig(
        tier="small",
        timeout_seconds=30.0,
        alignment_sequence_count=256,
        alignment_length=512,
        supermatrix_taxon_count=256,
        supermatrix_locus_lengths=(160, 224, 320),
        tree_set_tree_count=256,
        tree_set_tip_count=64,
        comparative_taxon_count=256,
        table_tip_count=256,
    ),
    "heavy": _StressTierConfig(
        tier="heavy",
        timeout_seconds=180.0,
        alignment_sequence_count=1024,
        alignment_length=2048,
        supermatrix_taxon_count=512,
        supermatrix_locus_lengths=(320, 448, 512, 640),
        tree_set_tree_count=1024,
        tree_set_tip_count=96,
        comparative_taxon_count=512,
        table_tip_count=1024,
    ),
}


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


def _resolve_stress_tier_config(tier: str) -> _StressTierConfig:
    try:
        return _STRESS_TIER_CONFIGS[tier]
    except KeyError as error:
        supported = ", ".join(sorted(_STRESS_TIER_CONFIGS))
        raise ValueError(
            f"unsupported stress tier '{tier}'; expected one of: {supported}"
        ) from error


def _build_balanced_tree(
    tip_count: int, *, branch_length: float = 0.1, prefix: str = "Taxon"
) -> PhyloTree:
    if tip_count < 2:
        raise ValueError(f"tip_count must be at least 2, got {tip_count}")
    leaves = [
        TreeNode(name=f"{prefix}{index}", branch_length=branch_length)
        for index in range(1, tip_count + 1)
    ]
    while len(leaves) > 1:
        next_level: list[TreeNode] = []
        for index in range(0, len(leaves), 2):
            left = leaves[index]
            right = leaves[index + 1] if index + 1 < len(leaves) else None
            if right is None:
                left.branch_length = round(
                    (left.branch_length or 0.0) + branch_length, 15
                )
                next_level.append(left)
                continue
            next_level.append(
                TreeNode(children=[left, right], branch_length=branch_length)
            )
        leaves = next_level
    root = leaves[0]
    root.branch_length = None
    return PhyloTree(root=root, source_format="newick")


def _stress_taxa(count: int) -> list[str]:
    return [f"taxon_{index:04d}" for index in range(1, count + 1)]


def _write_named_balanced_tree(
    path: Path,
    taxa: list[str],
    *,
    branch_length: float = 0.1,
) -> Path:
    if len(taxa) < 2:
        raise ValueError(f"expected at least two taxa, got {len(taxa)}")
    leaves = [TreeNode(name=taxon, branch_length=branch_length) for taxon in taxa]
    while len(leaves) > 1:
        next_level: list[TreeNode] = []
        for index in range(0, len(leaves), 2):
            left = leaves[index]
            right = leaves[index + 1] if index + 1 < len(leaves) else None
            if right is None:
                left.branch_length = round(
                    (left.branch_length or 0.0) + branch_length,
                    15,
                )
                next_level.append(left)
                continue
            next_level.append(
                TreeNode(children=[left, right], branch_length=branch_length)
            )
        leaves = next_level
    root = leaves[0]
    root.branch_length = None
    return write_newick(path, PhyloTree(root=root, source_format="newick"))


def _build_caterpillar_tree(
    tip_count: int, *, branch_length: float = 0.1, prefix: str = "Taxon"
) -> PhyloTree:
    if tip_count < 2:
        raise ValueError(f"tip_count must be at least 2, got {tip_count}")
    root = TreeNode(
        children=[
            TreeNode(name=f"{prefix}1", branch_length=branch_length),
            TreeNode(name=f"{prefix}2", branch_length=branch_length),
        ]
    )
    current = root
    for index in range(3, tip_count + 1):
        new_internal = TreeNode(
            branch_length=branch_length,
            children=[
                current.children.pop(),
                TreeNode(name=f"{prefix}{index}", branch_length=branch_length),
            ],
        )
        current.children.append(new_internal)
        current = new_internal
    root.branch_length = None
    return PhyloTree(root=root, source_format="newick")


def _write_fasttree_streaming_fixture(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if not args or "-help" in args:
    print("FastTree Version 2.2 benchmark fixture")
    raise SystemExit(0)

input_path = Path(args[-1])
identifiers = []
with input_path.open(encoding="utf-8") as handle:
    for raw_line in handle:
        line = raw_line.strip()
        if line.startswith(">"):
            identifiers.append(line[1:])

tips = [f"{identifier}:0.1" for identifier in identifiers]
while len(tips) > 1:
    left = tips.pop(0)
    right = tips.pop(0)
    tips.append(f"({left},{right})0.95:0.1")
if not tips:
    raise SystemExit(2)
print(tips[0] + ";")
print(
    "warning: benchmark fixture reconstructs one deterministic balanced tree",
    file=sys.stderr,
)
""",
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path


def _write_large_alignment(
    path: Path,
    *,
    sequence_count: int,
    sequence_length: int,
) -> Path:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for index in range(sequence_count):
            sequence = "".join(
                alphabet[(index + offset) % len(alphabet)]
                for offset in range(sequence_length)
            )
            handle.write(f">taxon_{index + 1:04d}\n{sequence}\n")
    return path


def _write_supermatrix_locus_alignments(
    root: Path,
    *,
    taxon_count: int,
    locus_lengths: tuple[int, ...],
) -> list[Path]:
    taxa = _stress_taxa(taxon_count)
    alphabet = "ACGT"
    root.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for locus_index, locus_length in enumerate(locus_lengths, start=1):
        path = root / f"locus_{locus_index:02d}.fasta"
        with path.open("w", encoding="utf-8") as handle:
            for taxon_index, taxon in enumerate(taxa):
                sequence = "".join(
                    alphabet[(taxon_index + offset + locus_index) % len(alphabet)]
                    for offset in range(locus_length)
                )
                handle.write(f">{taxon}\n{sequence}\n")
        paths.append(path)
    return paths


def _write_comparative_trait_table(path: Path, *, taxon_count: int) -> Path:
    rows: list[dict[str, str]] = []
    for index, taxon in enumerate(_stress_taxa(taxon_count), start=1):
        predictor_size = 10.0 + (index * 0.5)
        predictor_temperature = 18.0 + ((index % 17) * 0.25)
        response_mass = 3.0 + (predictor_size * 1.7) - (predictor_temperature * 0.12)
        response_rate = 1.5 + (predictor_size * 0.35) + (predictor_temperature * 0.2)
        rows.append(
            {
                "taxon": taxon,
                "response_mass": format(response_mass, ".12g"),
                "response_rate": format(response_rate, ".12g"),
                "predictor_size": format(predictor_size, ".12g"),
                "predictor_temperature": format(predictor_temperature, ".12g"),
            }
        )
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "response_mass",
            "response_rate",
            "predictor_size",
            "predictor_temperature",
        ],
        rows=rows,
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


def _large_alignment_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    executable = _write_fasttree_streaming_fixture(root / "FastTree-benchmark-fixture")
    alignment_path = _write_large_alignment(
        root / "large-alignment.fasta",
        sequence_count=config.alignment_sequence_count,
        sequence_length=config.alignment_length,
    )
    report = run_large_alignment_inference(
        alignment_path,
        out_dir=root / "large-alignment-out",
        prefix=f"{config.tier}-alignment",
        sequence_type="protein",
        executable=executable,
        timeout_seconds=config.timeout_seconds,
    )
    stage_peak_memory = max(
        (row.peak_memory_bytes or 0) for row in report.resource_rows
    )
    payload = _StressObservationPayload(
        workload="large-alignment-inference",
        input_size_bytes=report.input_summary.input_bytes,
        sequence_count=report.input_summary.sequence_count,
        alignment_length=report.input_summary.alignment_length,
        tree_count=None,
        taxon_count=report.input_summary.sequence_count,
        locus_count=None,
        output_row_count=len(report.resource_rows) + len(report.support_summary.nodes),
        notes=[
            "uses governed large-alignment inference with a deterministic FastTree fixture",
            "checks streamed rectangularity validation before tree inference",
            "records support-node and resource-ledger row counts from the large-alignment workflow",
        ],
    )
    return (
        payload,
        stage_peak_memory,
        _summarize_memory_observation_kinds(
            [row.memory_observation_kind for row in report.resource_rows]
        ),
    )


def _supermatrix_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    alignment_paths = _write_supermatrix_locus_alignments(
        root / "loci",
        taxon_count=config.supermatrix_taxon_count,
        locus_lengths=config.supermatrix_locus_lengths,
    )
    _, _, report = concatenate_locus_alignments(alignment_paths)
    payload = _StressObservationPayload(
        workload="multi-locus-supermatrix",
        input_size_bytes=sum(path.stat().st_size for path in alignment_paths),
        sequence_count=config.supermatrix_taxon_count,
        alignment_length=report.alignment_length,
        tree_count=None,
        taxon_count=report.taxon_count,
        locus_count=report.locus_count,
        output_row_count=(
            len(report.loci)
            + len(report.occupancy_report.taxa)
            + len(report.occupancy_report.loci)
            + len(report.occupancy_report.cells)
        ),
        notes=[
            "assembles one partitioned supermatrix and occupancy ledger",
            "counts taxon, locus, and occupancy review rows for the concatenated matrix",
        ],
    )
    return payload, None, None


def _tree_set_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    trees, _ = simulate_birth_death_trees(
        tree_count=config.tree_set_tree_count,
        tip_count=config.tree_set_tip_count,
        seed=config.tree_set_tree_count,
        taxon_prefix="tree_taxon_",
    )
    tree_set_path = write_tree_set(root / "tree-set.nwk", trees)
    consensus_tree, consensus_report = compute_consensus_tree(tree_set_path)
    consensus_path = write_newick(root / "consensus.nwk", consensus_tree)
    clade_report = extract_tree_clades(consensus_path)
    payload = _StressObservationPayload(
        workload="posterior-tree-set-consensus",
        input_size_bytes=tree_set_path.stat().st_size,
        sequence_count=None,
        alignment_length=None,
        tree_count=consensus_report.tree_count,
        taxon_count=len(consensus_report.shared_taxa),
        locus_count=None,
        output_row_count=len(clade_report.rows),
        notes=[
            "aggregates one governed consensus tree from a large tree set",
            "counts clade rows on the consensus output rather than materializing one table per sampled tree",
        ],
    )
    return payload, None, None


def _comparative_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    tree_path = _write_named_balanced_tree(
        root / "comparative-tree.nwk",
        _stress_taxa(config.comparative_taxon_count),
    )
    traits_path = _write_comparative_trait_table(
        root / "comparative-traits.tsv",
        taxon_count=config.comparative_taxon_count,
    )
    report = compute_phylogenetic_independent_contrasts(
        tree_path,
        traits_path,
        trait="response_mass",
    )
    payload = _StressObservationPayload(
        workload="comparative-trait-contrasts",
        input_size_bytes=tree_path.stat().st_size + traits_path.stat().st_size,
        sequence_count=None,
        alignment_length=None,
        tree_count=None,
        taxon_count=report.taxon_count,
        locus_count=None,
        output_row_count=len(report.contrasts),
        notes=[
            "runs the governed phylogenetic independent-contrasts surface on a large comparative trait table",
            "counts one reviewer-facing contrast row per internal node",
        ],
    )
    return payload, None, None


def _table_generation_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    tree_path = write_newick(
        root / "table-tree.nwk",
        _build_balanced_tree(config.table_tip_count, prefix="TableTaxon"),
    )
    clade_report = extract_tree_clades(tree_path)
    branch_report = analyze_branch_length_distribution(tree_path)
    payload = _StressObservationPayload(
        workload="tree-annotation-tables",
        input_size_bytes=tree_path.stat().st_size,
        sequence_count=None,
        alignment_length=None,
        tree_count=1,
        taxon_count=config.table_tip_count,
        locus_count=None,
        output_row_count=len(clade_report.rows) + len(branch_report.rows),
        notes=[
            "materializes clade and branch review tables from one large tree",
            "measures row generation directly on the reviewer-facing table surfaces",
        ],
    )
    return payload, None, None


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
