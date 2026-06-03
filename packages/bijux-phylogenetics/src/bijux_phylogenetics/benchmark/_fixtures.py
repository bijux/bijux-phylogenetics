from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.engines.inference import (
    run_large_alignment_inference,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.phylo.alignment.concatenation import (
    concatenate_locus_alignments,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    write_tree_set,
)
from bijux_phylogenetics.trees import (
    analyze_branch_length_distribution,
    compute_consensus_tree,
    extract_tree_clades,
)

from .contracts import _StressObservationPayload, _StressTierConfig

STRESS_TIER_CONFIGS: dict[str, _StressTierConfig] = {
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

LARGE_TREE_SCALING_TIP_COUNTS: tuple[int, ...] = (256, 1024, 2048)
LARGE_ALIGNMENT_SCALING_CLASSES: tuple[tuple[str, int, int], ...] = (
    ("sequences-256-sites-512", 256, 512),
    ("sequences-512-sites-1024", 512, 1024),
    ("sequences-1024-sites-2048", 1024, 2048),
)
LARGE_TREE_SET_SCALING_CLASSES: tuple[tuple[str, int, int], ...] = (
    ("trees-128-taxa-48", 128, 48),
    ("trees-256-taxa-64", 256, 64),
    ("trees-384-taxa-96", 384, 96),
)


def resolve_stress_tier_config(tier: str) -> _StressTierConfig:
    try:
        return STRESS_TIER_CONFIGS[tier]
    except KeyError as error:
        supported = ", ".join(sorted(STRESS_TIER_CONFIGS))
        raise ValueError(
            f"unsupported stress tier '{tier}'; expected one of: {supported}"
        ) from error


def build_balanced_tree(
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


def build_caterpillar_tree(
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


def stress_taxa(count: int) -> list[str]:
    return [f"taxon_{index:04d}" for index in range(1, count + 1)]


def interleaved_taxa(tip_count: int, *, prefix: str = "Taxon") -> list[str]:
    taxa = [f"{prefix}{index}" for index in range(1, tip_count + 1)]
    return [*taxa[::2], *taxa[1::2]]


def write_named_balanced_tree(
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


def write_fasttree_streaming_fixture(path: Path) -> Path:
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


def write_trimal_benchmark_fixture(path: Path) -> Path:
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


def write_large_alignment(
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


def write_supermatrix_locus_alignments(
    root: Path,
    *,
    taxon_count: int,
    locus_lengths: tuple[int, ...],
) -> list[Path]:
    taxa = stress_taxa(taxon_count)
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


def write_comparative_trait_table(path: Path, *, taxon_count: int) -> Path:
    rows: list[dict[str, str]] = []
    for index, taxon in enumerate(stress_taxa(taxon_count), start=1):
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


def large_alignment_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    executable = write_fasttree_streaming_fixture(root / "FastTree-benchmark-fixture")
    alignment_path = write_large_alignment(
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
    return payload, stage_peak_memory, None


def supermatrix_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    alignment_paths = write_supermatrix_locus_alignments(
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


def tree_set_stress_payload(
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


def comparative_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    tree_path = write_named_balanced_tree(
        root / "comparative-tree.nwk",
        stress_taxa(config.comparative_taxon_count),
    )
    traits_path = write_comparative_trait_table(
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


def table_generation_stress_payload(
    *,
    root: Path,
    config: _StressTierConfig,
) -> tuple[_StressObservationPayload, int | None, str | None]:
    tree_path = write_newick(
        root / "table-tree.nwk",
        build_balanced_tree(config.table_tip_count, prefix="TableTaxon"),
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
