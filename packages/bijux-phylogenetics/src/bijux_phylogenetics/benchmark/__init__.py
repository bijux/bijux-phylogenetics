from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import tempfile
import time
import tracemalloc

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

from bijux_phylogenetics.trees import analyze_branch_length_distribution
from bijux_phylogenetics.trees import extract_tree_clades
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.core.concatenation import concatenate_locus_alignments
from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.distance import build_distance_method_report
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.engines.large_alignment_inference import (
    run_large_alignment_inference,
)
from bijux_phylogenetics.io.fasta import build_alignment_quality_report
from bijux_phylogenetics.io.fasta import summarize_alignment_readiness
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
from bijux_phylogenetics.engines.workflows import run_alignment_trimming


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


@dataclass(slots=True)
class LargeTreeScalingWorkflowBenchmark:
    workflow: str
    scaling_axis: str
    observations: list[BenchmarkObservation]
    notes: list[str]


@dataclass(slots=True)
class LargeTreeScalingBenchmarkReport:
    replicates: int
    tip_counts: list[int]
    workflows: list[LargeTreeScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class LargeAlignmentScalingObservation:
    label: str
    sequence_count: int
    alignment_length: int
    aligned_site_count: int
    runtime_seconds: float
    peak_memory_bytes: int


@dataclass(slots=True)
class LargeAlignmentScalingWorkflowBenchmark:
    workflow: str
    scaling_axis: str
    observations: list[LargeAlignmentScalingObservation]
    notes: list[str]


@dataclass(slots=True)
class LargeAlignmentScalingBenchmarkReport:
    replicates: int
    sequence_counts: list[int]
    alignment_lengths: list[int]
    workflows: list[LargeAlignmentScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class LargeTreeSetScalingObservation:
    label: str
    tree_count: int
    tip_count: int
    pair_count: int
    runtime_seconds: float
    peak_memory_bytes: int


@dataclass(slots=True)
class LargeTreeSetScalingWorkflowBenchmark:
    workflow: str
    scaling_axis: str
    observations: list[LargeTreeSetScalingObservation]
    notes: list[str]


@dataclass(slots=True)
class LargeTreeSetScalingBenchmarkReport:
    replicates: int
    tree_counts: list[int]
    tip_counts: list[int]
    workflows: list[LargeTreeSetScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class WorkflowPracticalLimitEntry:
    workflow: str
    evidence_source: str
    tested_taxon_limit: int | None
    tested_site_limit: int | None
    tested_tree_limit: int | None
    tested_posterior_size: int | None
    max_runtime_seconds: float
    max_peak_memory_bytes: int
    memory_observation_kind: str | None
    notes: list[str]


@dataclass(slots=True)
class WorkflowPracticalLimitReport:
    replicates: int
    stress_tiers: list[str]
    entries: list[WorkflowPracticalLimitEntry]
    limitations: list[str]


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

_LARGE_TREE_SCALING_TIP_COUNTS: tuple[int, ...] = (256, 1024, 2048)
_LARGE_ALIGNMENT_SCALING_CLASSES: tuple[tuple[str, int, int], ...] = (
    ("sequences-256-sites-512", 256, 512),
    ("sequences-512-sites-1024", 512, 1024),
    ("sequences-1024-sites-2048", 1024, 2048),
)
_LARGE_TREE_SET_SCALING_CLASSES: tuple[tuple[str, int, int], ...] = (
    ("trees-128-taxa-48", 128, 48),
    ("trees-256-taxa-64", 256, 64),
    ("trees-384-taxa-96", 384, 96),
)


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
                current.pop_child(),
                TreeNode(name=f"{prefix}{index}", branch_length=branch_length),
            ],
        )
        current.append_child(new_internal)
        current = new_internal
    root.branch_length = None
    return PhyloTree(root=root, source_format="newick")


def _interleaved_taxa(tip_count: int, *, prefix: str = "Taxon") -> list[str]:
    taxa = [f"{prefix}{index}" for index in range(1, tip_count + 1)]
    return [*taxa[::2], *taxa[1::2]]


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


def _write_trimal_benchmark_fixture(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        """#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    print("trimAl v2.0 benchmark fixture")
    raise SystemExit(0)

input_path = Path(args[args.index("-in") + 1])
output_path = Path(args[args.index("-out") + 1])
records = []
identifier = None
sequence = []
for raw_line in input_path.read_text(encoding="utf-8").splitlines():
    line = raw_line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if identifier is not None:
            records.append((identifier, "".join(sequence)))
        identifier = line[1:]
        sequence = []
    else:
        sequence.append(line)
if identifier is not None:
    records.append((identifier, "".join(sequence)))

trim_width = 4 if "-gappyout" in args else 2
with output_path.open("w", encoding="utf-8") as handle:
    for identifier, sequence in records:
        handle.write(f">{identifier}\\n{sequence[:-trim_width]}\\n")
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


def benchmark_large_tree_scaling(
    *,
    replicates: int = 1,
    tip_counts: list[int] | None = None,
) -> LargeTreeScalingBenchmarkReport:
    """Benchmark large-tree validation, comparison, rendering, and reporting."""
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    counts = list(tip_counts or _LARGE_TREE_SCALING_TIP_COUNTS)
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
    classes = list(size_classes or _LARGE_ALIGNMENT_SCALING_CLASSES)
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
    classes = list(size_classes or _LARGE_TREE_SET_SCALING_CLASSES)
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
