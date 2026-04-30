from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
from time import perf_counter
import tempfile

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BioTree

from bijux_phylogenetics.compare.topology import _informative_clade_nodes, _informative_clades, _unrooted_splits
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import InvalidAlignmentError
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.simulation import simulate_birth_death_trees, write_tree_set


@dataclass(frozen=True, slots=True)
class TreeSetRecord:
    index: int
    tip_count: int
    taxa: list[str]
    rooted_topology_id: str
    unrooted_topology_id: str


@dataclass(slots=True)
class TreeSetReport:
    path: Path
    source_format: str
    tree_count: int
    shared_taxa: list[str]
    taxa_union: list[str]
    rooted_topology_count: int
    unrooted_topology_count: int
    records: list[TreeSetRecord]


@dataclass(frozen=True, slots=True)
class CladeFrequency:
    clade: str
    tree_count: int
    frequency: float


@dataclass(slots=True)
class CladeFrequencyReport:
    path: Path
    tree_count: int
    shared_taxa: list[str]
    clade_frequencies: list[CladeFrequency]


@dataclass(slots=True)
class ConsensusTreeReport:
    path: Path
    tree_count: int
    shared_taxa: list[str]
    consensus_newick: str


@dataclass(frozen=True, slots=True)
class TreeDistancePair:
    left_index: int
    right_index: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float


@dataclass(slots=True)
class TreeDistanceMatrixReport:
    path: Path
    tree_count: int
    shared_taxa: list[str]
    pairs: list[TreeDistancePair]


@dataclass(frozen=True, slots=True)
class TreeTopologyCluster:
    rooted_topology_id: str
    tree_indices: list[int]
    tree_count: int
    frequency: float
    representative_index: int
    representative_newick: str


@dataclass(slots=True)
class TreeTopologyClusterReport:
    path: Path
    tree_count: int
    rooted_topology_count: int
    clusters: list[TreeTopologyCluster]


@dataclass(frozen=True, slots=True)
class TaxonPlacementSignature:
    signature: str
    tree_count: int
    frequency: float


@dataclass(frozen=True, slots=True)
class UnstableTaxon:
    taxon: str
    unique_placements: int
    dominant_frequency: float
    instability_score: float
    placements: list[TaxonPlacementSignature]


@dataclass(slots=True)
class UnstableTaxaReport:
    path: Path
    tree_count: int
    taxa: list[UnstableTaxon]


@dataclass(frozen=True, slots=True)
class UnstableClade:
    clade: str
    tree_count: int
    frequency: float
    conflict_count: int
    instability_score: float
    support_classification: str
    conflicting_clades: list[str]


@dataclass(slots=True)
class UnstableCladeReport:
    path: Path
    tree_count: int
    clades: list[UnstableClade]


@dataclass(frozen=True, slots=True)
class CladeFrequencyDelta:
    clade: str
    left_frequency: float
    right_frequency: float
    delta: float


@dataclass(slots=True)
class PosteriorTreeSetComparisonReport:
    left_path: Path
    right_path: Path
    shared_taxa: list[str]
    left_tree_count: int
    right_tree_count: int
    left_rooted_topology_count: int
    right_rooted_topology_count: int
    shared_rooted_topology_count: int
    mean_between_set_robinson_foulds: float
    mean_between_set_normalized_robinson_foulds: float
    clade_frequency_deltas: list[CladeFrequencyDelta]


@dataclass(frozen=True, slots=True)
class PosteriorTopologicalDiversitySummary:
    tree_count: int
    rooted_topology_count: int
    dominant_topology_frequency: float
    effective_topology_count: float
    mean_within_set_robinson_foulds: float
    mean_within_set_normalized_robinson_foulds: float


@dataclass(slots=True)
class PosteriorTopologicalDiversityComparisonReport:
    left_path: Path
    right_path: Path
    left_summary: PosteriorTopologicalDiversitySummary
    right_summary: PosteriorTopologicalDiversitySummary
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorTopologyMode:
    rooted_topology_id: str
    representative_index: int
    representative_newick: str
    tree_indices: list[int]
    tree_count: int
    frequency: float


@dataclass(slots=True)
class PosteriorTopologyMultimodalityReport:
    path: Path
    tree_count: int
    rooted_topology_count: int
    dominant_mode_frequency: float
    mode_count: int
    multimodal: bool
    modes: list[PosteriorTopologyMode]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class CladeCredibilityConflict:
    left_clade: str
    left_frequency: float
    right_clade: str
    right_frequency: float
    combined_frequency: float


@dataclass(slots=True)
class CladeCredibilityConflictReport:
    path: Path
    tree_count: int
    credibility_threshold: float
    high_credibility_clade_count: int
    conflict_count: int
    conflicts: list[CladeCredibilityConflict]


@dataclass(frozen=True, slots=True)
class UncertaintyAwareCladeConclusion:
    clade: str
    frequency: float
    conclusion: str
    rationale: str


@dataclass(slots=True)
class UncertaintyAwareConclusionSummaryReport:
    path: Path
    tree_count: int
    robust_clade_count: int
    uncertain_clade_count: int
    conflicting_clade_count: int
    robust_clades: list[UncertaintyAwareCladeConclusion]
    uncertain_clades: list[UncertaintyAwareCladeConclusion]
    conflicting_clades: list[UncertaintyAwareCladeConclusion]


@dataclass(frozen=True, slots=True)
class BootstrapPosteriorCladeComparison:
    clade: str
    bootstrap_support: float | None
    posterior_frequency: float | None
    absolute_delta: float | None
    agreement: str


@dataclass(slots=True)
class BootstrapPosteriorSupportComparisonReport:
    bootstrap_tree_path: Path
    posterior_tree_set_path: Path
    posterior_tree_count: int
    shared_taxa: list[str]
    high_conflict_clade_count: int
    rows: list[BootstrapPosteriorCladeComparison]


@dataclass(frozen=True, slots=True)
class TreeSetBenchmarkRow:
    tree_count: int
    taxon_count: int
    replicate: int
    elapsed_seconds: float
    rooted_topology_count: int
    unstable_taxon_count: int
    unstable_clade_count: int
    robust_clade_count: int


@dataclass(slots=True)
class TreeSetScalingBenchmarkReport:
    tree_counts: list[int]
    taxon_counts: list[int]
    rows: list[TreeSetBenchmarkRow]


@dataclass(slots=True)
class TreeSetStorageRiskReport:
    path: Path
    file_size_bytes: int
    file_size_megabytes: float
    tree_count: int
    rooted_topology_count: int
    shared_taxon_count: int
    mean_bytes_per_tree: float
    risk_level: str
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TreeSetThinningSensitivityRow:
    thinning_interval: int
    retained_tree_count: int
    retained_fraction: float
    rooted_topology_count: int
    shared_rooted_topology_count: int
    dominant_topology_frequency: float
    dominant_topology_delta: float
    robust_clade_count: int
    uncertain_clade_count: int
    conflicting_clade_count: int
    warnings: list[str]


@dataclass(slots=True)
class TreeSetThinningSensitivityReport:
    path: Path
    original_tree_count: int
    original_rooted_topology_count: int
    original_dominant_topology_frequency: float
    rows: list[TreeSetThinningSensitivityRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class ConsensusThresholdSensitivityRow:
    threshold: float
    informative_clade_count: int
    rooted_topology_id: str
    consensus_newick: str
    warnings: list[str]


@dataclass(slots=True)
class ConsensusThresholdSensitivityReport:
    path: Path
    tree_count: int
    rows: list[ConsensusThresholdSensitivityRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class TreeSetMaturityGateCheck:
    name: str
    satisfied: bool
    details: str


@dataclass(slots=True)
class TreeSetMaturityGateReport:
    path: Path
    decision: str
    checks: list[TreeSetMaturityGateCheck]
    warnings: list[str]


def _shared_taxa(trees: list[PhyloTree]) -> set[str]:
    shared = set(trees[0].tip_names)
    for tree in trees[1:]:
        shared &= set(tree.tip_names)
    return shared


def _taxa_union(trees: list[PhyloTree]) -> set[str]:
    taxa: set[str] = set()
    for tree in trees:
        taxa.update(tree.tip_names)
    return taxa


def _require_tree_set(path: Path) -> tuple[str, list[BioTree], list[PhyloTree]]:
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    bio_trees = list(Phylo.parse(path, source_format))
    if not bio_trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    trees = [tree_from_biophylo(tree, source_format=source_format) for tree in bio_trees]
    return source_format, bio_trees, trees


def _format_clade(clade: frozenset[str]) -> str:
    return "|".join(sorted(clade))


def _rooted_topology_id(tree: PhyloTree, shared_taxa: set[str]) -> str:
    return "||".join(sorted(_format_clade(clade) for clade in _informative_clades(tree, shared_taxa)))


def _unrooted_topology_id(tree: PhyloTree, shared_taxa: set[str]) -> str:
    return "||".join(sorted(_format_clade(clade) for clade in _unrooted_splits(tree, shared_taxa)))


def _validate_same_taxa(trees: list[PhyloTree]) -> list[str]:
    first = sorted(trees[0].tip_names)
    for tree in trees[1:]:
        if sorted(tree.tip_names) != first:
            raise InvalidAlignmentError("tree-set analysis requires all trees to share the exact same taxon set")
    return first


def _clade_signature(tree: PhyloTree, shared_taxa: set[str], taxon: str) -> str:
    containing_clades = sorted(
        _format_clade(clade)
        for clade in _informative_clades(tree, shared_taxa)
        if taxon in clade
    )
    if not containing_clades:
        return "(singleton-placement)"
    return "||".join(containing_clades)


def _clade_counts(trees: list[PhyloTree], shared_taxa: set[str]) -> dict[frozenset[str], int]:
    counts: dict[frozenset[str], int] = {}
    for tree in trees:
        for clade in _informative_clades(tree, shared_taxa):
            counts[clade] = counts.get(clade, 0) + 1
    return counts


def _clades_conflict(left: frozenset[str], right: frozenset[str]) -> bool:
    return bool(left & right) and not (left <= right or right <= left)


def _clade_branch_lengths(tree: PhyloTree, shared_taxa: set[str]) -> dict[frozenset[str], float | None]:
    lengths: dict[frozenset[str], float | None] = {}

    def visit(node: TreeNode) -> set[str]:
        if node.is_leaf():
            return {node.name} if node.name is not None and node.name in shared_taxa else set()
        taxa: set[str] = set()
        for child in node.children:
            taxa.update(visit(child))
        if 1 < len(taxa) < len(shared_taxa):
            lengths[frozenset(taxa)] = node.branch_length
        return taxa

    visit(tree.root)
    return lengths


def _terminal_branch_lengths(tree: PhyloTree) -> dict[str, float | None]:
    return {name: length for name, length in tree.terminal_branch_lengths() if name is not None}


def _maximal_nested_clades(parent: frozenset[str], clades: set[frozenset[str]]) -> list[frozenset[str]]:
    nested = [clade for clade in clades if clade < parent]
    return sorted(
        [
            clade
            for clade in nested
            if not any(clade < other < parent for other in nested)
        ],
        key=lambda clade: (len(clade), sorted(clade)),
    )


def _mean(values: list[float]) -> float:
    return round(sum(values) / len(values), 15)


def _shannon_effective_count(frequencies: list[float]) -> float:
    if not frequencies:
        return 0.0
    return round(math.exp(-sum(frequency * math.log(frequency) for frequency in frequencies if frequency > 0.0)), 15)


def _mean_pairwise_distance(trees: list[PhyloTree], shared_taxa: set[str]) -> tuple[float, float]:
    comparisons: list[tuple[int, float]] = []
    for left_index, left in enumerate(trees):
        for right in trees[left_index + 1 :]:
            comparisons.append(_tree_distance(left, right, shared_taxa))
    if not comparisons:
        return 0.0, 0.0
    return (
        round(sum(distance for distance, _ in comparisons) / len(comparisons), 15),
        round(sum(normalized for _, normalized in comparisons) / len(comparisons), 15),
    )


def _support_classification(frequency: float, conflict_count: int) -> str:
    if frequency >= 0.9 and conflict_count == 0:
        return "robust"
    if 0.3 <= frequency <= 0.7:
        return "intermediate-support"
    if conflict_count > 0:
        return "credibility-conflicted"
    return "weak-support"


def _representative_tree_by_indices(trees: list[PhyloTree], indices: list[int]) -> tuple[int, str]:
    representative_index = indices[0]
    representative_tree = trees[representative_index - 1]
    return representative_index, dumps_newick(representative_tree)


def _topology_modes_from_clusters(
    trees: list[PhyloTree],
    clusters: list[TreeTopologyCluster],
    *,
    min_mode_frequency: float,
) -> list[PosteriorTopologyMode]:
    return [
        PosteriorTopologyMode(
            rooted_topology_id=cluster.rooted_topology_id,
            representative_index=cluster.representative_index,
            representative_newick=cluster.representative_newick,
            tree_indices=cluster.tree_indices,
            tree_count=cluster.tree_count,
            frequency=cluster.frequency,
        )
        for cluster in clusters
        if cluster.frequency >= min_mode_frequency
    ]


def _build_consensus_node(
    taxa: frozenset[str],
    *,
    majority_clades: set[frozenset[str]],
    clade_support: dict[frozenset[str], float],
    clade_lengths: dict[frozenset[str], float],
    terminal_lengths: dict[str, float],
    is_root: bool = False,
) -> TreeNode:
    child_clades = _maximal_nested_clades(taxa, majority_clades)
    covered: set[str] = set()
    children: list[TreeNode] = []
    for child_clade in child_clades:
        covered.update(child_clade)
        children.append(
            _build_consensus_node(
                child_clade,
                majority_clades=majority_clades,
                clade_support=clade_support,
                clade_lengths=clade_lengths,
                terminal_lengths=terminal_lengths,
            )
        )
    for taxon in sorted(taxa - covered):
        children.append(TreeNode(name=taxon, branch_length=terminal_lengths.get(taxon)))
    if len(children) == 1:
        return children[0]
    label = None if is_root else format(clade_support[taxa], ".15g")
    return TreeNode(
        name=label,
        branch_length=None if is_root else clade_lengths.get(taxa),
        children=children,
    )


def _tree_distance(left: PhyloTree, right: PhyloTree, shared_taxa: set[str]) -> tuple[int, float]:
    left_clades = _informative_clades(left, shared_taxa)
    right_clades = _informative_clades(right, shared_taxa)
    symmetric_difference = left_clades.symmetric_difference(right_clades)
    denominator = len(left_clades) + len(right_clades)
    normalized = 0.0 if denominator == 0 else len(symmetric_difference) / denominator
    return len(symmetric_difference), normalized


def load_tree_set(path: Path) -> TreeSetReport:
    """Read a set of trees and summarize their topology diversity over shared taxa."""
    source_format, _, trees = _require_tree_set(path)
    shared_taxa = sorted(_shared_taxa(trees))
    taxa_union = sorted(_taxa_union(trees))
    records = [
        TreeSetRecord(
            index=index,
            tip_count=tree.tip_count,
            taxa=sorted(tree.tip_names),
            rooted_topology_id=_rooted_topology_id(tree, set(shared_taxa)),
            unrooted_topology_id=_unrooted_topology_id(tree, set(shared_taxa)),
        )
        for index, tree in enumerate(trees, start=1)
    ]
    return TreeSetReport(
        path=path,
        source_format=source_format,
        tree_count=len(trees),
        shared_taxa=shared_taxa,
        taxa_union=taxa_union,
        rooted_topology_count=len({record.rooted_topology_id for record in records}),
        unrooted_topology_count=len({record.unrooted_topology_id for record in records}),
        records=records,
    )


def compute_clade_frequency_table(path: Path) -> CladeFrequencyReport:
    """Compute informative clade frequencies across a tree set with a shared taxon set."""
    _, _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    counts: dict[str, int] = {}
    for tree in trees:
        for clade in _informative_clades(tree, shared_taxa):
            counts[_format_clade(clade)] = counts.get(_format_clade(clade), 0) + 1
    total = len(trees)
    return CladeFrequencyReport(
        path=path,
        tree_count=total,
        shared_taxa=sorted(shared_taxa),
        clade_frequencies=[
            CladeFrequency(
                clade=clade,
                tree_count=count,
                frequency=round(count / total, 15),
            )
            for clade, count in sorted(counts.items())
        ],
    )


def write_clade_frequency_table(path: Path, report: CladeFrequencyReport) -> Path:
    """Write a clade-frequency table as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["clade", "tree_count", "frequency"], delimiter="\t")
        writer.writeheader()
        for row in report.clade_frequencies:
            writer.writerow(
                {
                    "clade": row.clade,
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                }
            )
    return path


def compute_consensus_tree(path: Path) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a majority-rule consensus tree from a tree set."""
    return compute_consensus_tree_with_threshold(path, threshold=0.5)


def compute_consensus_tree_with_threshold(
    path: Path,
    *,
    threshold: float,
) -> tuple[PhyloTree, ConsensusTreeReport]:
    """Compute a deterministic consensus tree at a caller-supplied clade frequency threshold."""
    if not 0.0 < threshold < 1.0:
        raise ValueError(f"consensus threshold must be between 0 and 1, got {threshold}")
    source_format, _, trees = _require_tree_set(path)
    shared_taxa = _validate_same_taxa(trees)
    universe = frozenset(shared_taxa)
    counts: dict[frozenset[str], int] = {}
    branch_lengths_by_clade: dict[frozenset[str], list[float]] = {}
    terminal_lengths: dict[str, list[float]] = {}
    for tree in trees:
        for clade, length in _clade_branch_lengths(tree, set(shared_taxa)).items():
            counts[clade] = counts.get(clade, 0) + 1
            if length is not None:
                branch_lengths_by_clade.setdefault(clade, []).append(float(length))
        for taxon, length in _terminal_branch_lengths(tree).items():
            if length is not None:
                terminal_lengths.setdefault(taxon, []).append(float(length))
    majority_clades = {
        clade
        for clade, count in counts.items()
        if count / len(trees) >= threshold
    }
    clade_support = {
        clade: round((counts[clade] / len(trees)) * 100.0, 15)
        for clade in majority_clades
    }
    clade_lengths = {
        clade: _mean(lengths)
        for clade, lengths in branch_lengths_by_clade.items()
        if clade in majority_clades and lengths
    }
    terminal_length_means = {
        taxon: _mean(lengths)
        for taxon, lengths in terminal_lengths.items()
        if lengths
    }
    tree = PhyloTree(
        root=_build_consensus_node(
            universe,
            majority_clades=majority_clades,
            clade_support=clade_support,
            clade_lengths=clade_lengths,
            terminal_lengths=terminal_length_means,
            is_root=True,
        ),
        source_format=source_format,
    )
    return tree, ConsensusTreeReport(
        path=path,
        tree_count=len(trees),
        shared_taxa=shared_taxa,
        consensus_newick=dumps_newick(tree),
    )


def write_consensus_tree(path: Path, tree: PhyloTree) -> Path:
    """Write a consensus tree as canonical Newick."""
    return write_newick(path, tree)


def compute_tree_distance_matrix(path: Path) -> TreeDistanceMatrixReport:
    """Compute a pairwise RF-distance matrix across a tree set."""
    _, _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    pairs: list[TreeDistancePair] = []
    for left_index, left in enumerate(trees, start=1):
        for right_index, right in enumerate(trees[left_index - 1 :], start=left_index):
            distance, normalized = _tree_distance(left, right, shared_taxa)
            pairs.append(
                TreeDistancePair(
                    left_index=left_index,
                    right_index=right_index,
                    robinson_foulds_distance=distance,
                    normalized_robinson_foulds=normalized,
                )
            )
    return TreeDistanceMatrixReport(
        path=path,
        tree_count=len(trees),
        shared_taxa=sorted(shared_taxa),
        pairs=pairs,
    )


def write_tree_distance_matrix(path: Path, report: TreeDistanceMatrixReport) -> Path:
    """Write a pairwise tree-distance matrix as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "left_index",
                "right_index",
                "robinson_foulds_distance",
                "normalized_robinson_foulds",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.pairs:
            writer.writerow(
                {
                    "left_index": row.left_index,
                    "right_index": row.right_index,
                    "robinson_foulds_distance": row.robinson_foulds_distance,
                    "normalized_robinson_foulds": format(row.normalized_robinson_foulds, ".15g"),
                }
            )
    return path


def write_topology_cluster_table(path: Path, report: TreeTopologyClusterReport) -> Path:
    """Write rooted topology clusters as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rooted_topology_id",
                "tree_indices",
                "tree_count",
                "frequency",
                "representative_index",
                "representative_newick",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.clusters:
            writer.writerow(
                {
                    "rooted_topology_id": row.rooted_topology_id,
                    "tree_indices": ",".join(str(index) for index in row.tree_indices),
                    "tree_count": row.tree_count,
                    "frequency": format(row.frequency, ".15g"),
                    "representative_index": row.representative_index,
                    "representative_newick": row.representative_newick,
                }
            )
    return path


def cluster_trees_by_topology(path: Path) -> TreeTopologyClusterReport:
    """Cluster trees by identical rooted topology signatures."""
    report = load_tree_set(path)
    clusters_by_id: dict[str, list[int]] = {}
    for record in report.records:
        clusters_by_id.setdefault(record.rooted_topology_id, []).append(record.index)
    _, _, trees = _require_tree_set(path)
    clusters: list[TreeTopologyCluster] = []
    for topology_id, indices in sorted(
        clusters_by_id.items(),
        key=lambda item: (-len(item[1]), item[1][0]),
    ):
        representative_index, representative_newick = _representative_tree_by_indices(trees, indices)
        clusters.append(
            TreeTopologyCluster(
                rooted_topology_id=topology_id,
                tree_indices=indices,
                tree_count=len(indices),
                frequency=round(len(indices) / report.tree_count, 15),
                representative_index=representative_index,
                representative_newick=representative_newick,
            )
        )
    return TreeTopologyClusterReport(
        path=path,
        tree_count=report.tree_count,
        rooted_topology_count=report.rooted_topology_count,
        clusters=clusters,
    )


def detect_unstable_taxa(path: Path) -> UnstableTaxaReport:
    """Report taxa whose placement signatures vary across trees in a set."""
    _, _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    taxa: list[UnstableTaxon] = []
    for taxon in sorted(shared_taxa):
        signature_counts: dict[str, int] = {}
        for tree in trees:
            signature = _clade_signature(tree, shared_taxa, taxon)
            signature_counts[signature] = signature_counts.get(signature, 0) + 1
        if len(signature_counts) < 2:
            continue
        placements = [
            TaxonPlacementSignature(
                signature=signature,
                tree_count=count,
                frequency=round(count / len(trees), 15),
            )
            for signature, count in sorted(
                signature_counts.items(),
                key=lambda item: (-item[1], item[0]),
            )
        ]
        taxa.append(
            UnstableTaxon(
                taxon=taxon,
                unique_placements=len(signature_counts),
                dominant_frequency=placements[0].frequency,
                instability_score=round(1.0 - placements[0].frequency, 15),
                placements=placements,
            )
        )
    taxa.sort(key=lambda row: (-row.instability_score, -row.unique_placements, row.taxon))
    return UnstableTaxaReport(path=path, tree_count=len(trees), taxa=taxa)


def detect_unstable_clades(path: Path) -> UnstableCladeReport:
    """Report non-unanimous clades and their conflicting alternatives."""
    _, _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    counts = _clade_counts(trees, shared_taxa)
    all_clades = set(counts)
    unstable_clades = [
        UnstableClade(
            clade=_format_clade(clade),
            tree_count=count,
            frequency=round(count / len(trees), 15),
            conflict_count=len(conflicts := sorted(_format_clade(other) for other in all_clades if _clades_conflict(clade, other))),
            instability_score=round(min(count / len(trees), 1.0 - (count / len(trees))), 15),
            support_classification=_support_classification(round(count / len(trees), 15), len(conflicts)),
            conflicting_clades=conflicts,
        )
        for clade, count in sorted(
            counts.items(),
            key=lambda item: (_format_clade(item[0])),
        )
        if count < len(trees)
    ]
    unstable_clades.sort(key=lambda row: (-row.instability_score, -row.conflict_count, row.clade))
    return UnstableCladeReport(path=path, tree_count=len(trees), clades=unstable_clades)


def compare_posterior_topological_diversity(
    left_path: Path,
    right_path: Path,
) -> PosteriorTopologicalDiversityComparisonReport:
    """Compare topology diversity and dispersion across two posterior tree sets."""
    left_clusters = cluster_trees_by_topology(left_path)
    right_clusters = cluster_trees_by_topology(right_path)
    _, _, left_trees = _require_tree_set(left_path)
    _, _, right_trees = _require_tree_set(right_path)
    left_taxa = set(_validate_same_taxa(left_trees))
    right_taxa = set(_validate_same_taxa(right_trees))
    if left_taxa != right_taxa:
        raise InvalidAlignmentError("posterior diversity comparison requires identical taxon sets across both inputs")
    left_mean_rf, left_mean_normalized = _mean_pairwise_distance(left_trees, left_taxa)
    right_mean_rf, right_mean_normalized = _mean_pairwise_distance(right_trees, right_taxa)
    left_summary = PosteriorTopologicalDiversitySummary(
        tree_count=len(left_trees),
        rooted_topology_count=left_clusters.rooted_topology_count,
        dominant_topology_frequency=left_clusters.clusters[0].frequency,
        effective_topology_count=_shannon_effective_count([cluster.frequency for cluster in left_clusters.clusters]),
        mean_within_set_robinson_foulds=left_mean_rf,
        mean_within_set_normalized_robinson_foulds=left_mean_normalized,
    )
    right_summary = PosteriorTopologicalDiversitySummary(
        tree_count=len(right_trees),
        rooted_topology_count=right_clusters.rooted_topology_count,
        dominant_topology_frequency=right_clusters.clusters[0].frequency,
        effective_topology_count=_shannon_effective_count([cluster.frequency for cluster in right_clusters.clusters]),
        mean_within_set_robinson_foulds=right_mean_rf,
        mean_within_set_normalized_robinson_foulds=right_mean_normalized,
    )
    warnings: list[str] = []
    if left_summary.effective_topology_count != right_summary.effective_topology_count:
        richer = "left" if left_summary.effective_topology_count > right_summary.effective_topology_count else "right"
        warnings.append(f"{richer} analysis spans a broader effective topology spectrum")
    if abs(left_summary.mean_within_set_normalized_robinson_foulds - right_summary.mean_within_set_normalized_robinson_foulds) >= 0.15:
        warnings.append("posterior analyses differ materially in within-set topological dispersion")
    return PosteriorTopologicalDiversityComparisonReport(
        left_path=left_path,
        right_path=right_path,
        left_summary=left_summary,
        right_summary=right_summary,
        warnings=warnings,
    )


def detect_posterior_topology_multimodality(
    path: Path,
    *,
    min_mode_frequency: float = 0.2,
    min_mode_count: int = 2,
) -> PosteriorTopologyMultimodalityReport:
    """Report whether a posterior tree set contains multiple high-frequency topology modes."""
    if not 0.0 < min_mode_frequency <= 1.0:
        raise ValueError(f"min_mode_frequency must be between 0 and 1, got {min_mode_frequency}")
    if min_mode_count < 2:
        raise ValueError(f"min_mode_count must be at least 2, got {min_mode_count}")
    clusters = cluster_trees_by_topology(path)
    _, _, trees = _require_tree_set(path)
    modes = _topology_modes_from_clusters(trees, clusters.clusters, min_mode_frequency=min_mode_frequency)
    multimodal = len(modes) >= min_mode_count
    warnings: list[str] = []
    if multimodal:
        warnings.append("posterior topology distribution contains multiple high-frequency modes")
    if clusters.clusters and clusters.clusters[0].frequency < 0.75:
        warnings.append("no single topology dominates the posterior tree set")
    return PosteriorTopologyMultimodalityReport(
        path=path,
        tree_count=clusters.tree_count,
        rooted_topology_count=clusters.rooted_topology_count,
        dominant_mode_frequency=0.0 if not clusters.clusters else clusters.clusters[0].frequency,
        mode_count=len(modes),
        multimodal=multimodal,
        modes=modes,
        warnings=warnings,
    )


def summarize_clade_credibility_conflicts(
    path: Path,
    *,
    credibility_threshold: float = 0.5,
) -> CladeCredibilityConflictReport:
    """Identify mutually incompatible clades that both achieve high posterior credibility."""
    if not 0.0 < credibility_threshold < 1.0:
        raise ValueError(f"credibility_threshold must be between 0 and 1, got {credibility_threshold}")
    _, _, trees = _require_tree_set(path)
    shared_taxa = set(_validate_same_taxa(trees))
    counts = _clade_counts(trees, shared_taxa)
    frequencies = {
        clade: round(count / len(trees), 15)
        for clade, count in counts.items()
    }
    high_credibility = [
        clade
        for clade, frequency in frequencies.items()
        if frequency >= credibility_threshold
    ]
    conflicts: list[CladeCredibilityConflict] = []
    for index, left_clade in enumerate(sorted(high_credibility, key=_format_clade)):
        for right_clade in sorted(high_credibility[index + 1 :], key=_format_clade):
            if not _clades_conflict(left_clade, right_clade):
                continue
            conflicts.append(
                CladeCredibilityConflict(
                    left_clade=_format_clade(left_clade),
                    left_frequency=frequencies[left_clade],
                    right_clade=_format_clade(right_clade),
                    right_frequency=frequencies[right_clade],
                    combined_frequency=round(frequencies[left_clade] + frequencies[right_clade], 15),
                )
            )
    conflicts.sort(key=lambda row: (-row.combined_frequency, row.left_clade, row.right_clade))
    return CladeCredibilityConflictReport(
        path=path,
        tree_count=len(trees),
        credibility_threshold=credibility_threshold,
        high_credibility_clade_count=len(high_credibility),
        conflict_count=len(conflicts),
        conflicts=conflicts,
    )


def write_clade_credibility_conflict_table(path: Path, report: CladeCredibilityConflictReport) -> Path:
    """Write high-credibility clade conflicts as a TSV table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["left_clade", "left_frequency", "right_clade", "right_frequency", "combined_frequency"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.conflicts:
            writer.writerow(
                {
                    "left_clade": row.left_clade,
                    "left_frequency": format(row.left_frequency, ".15g"),
                    "right_clade": row.right_clade,
                    "right_frequency": format(row.right_frequency, ".15g"),
                    "combined_frequency": format(row.combined_frequency, ".15g"),
                }
            )
    return path


def summarize_uncertainty_aware_conclusions(
    path: Path,
    *,
    robust_threshold: float = 0.9,
    uncertain_min_frequency: float = 0.3,
    uncertain_max_frequency: float = 0.7,
    credibility_threshold: float = 0.5,
) -> UncertaintyAwareConclusionSummaryReport:
    """Classify clade-level conclusions as robust, uncertain, or conflict-prone."""
    if not 0.0 < robust_threshold <= 1.0:
        raise ValueError(f"robust_threshold must be between 0 and 1, got {robust_threshold}")
    unstable = detect_unstable_clades(path)
    conflict_report = summarize_clade_credibility_conflicts(path, credibility_threshold=credibility_threshold)
    conflict_clades = {
        row.left_clade
        for row in conflict_report.conflicts
    } | {
        row.right_clade
        for row in conflict_report.conflicts
    }
    frequency_report = compute_clade_frequency_table(path)
    robust_clades: list[UncertaintyAwareCladeConclusion] = []
    uncertain_clades: list[UncertaintyAwareCladeConclusion] = []
    conflicting_clades: list[UncertaintyAwareCladeConclusion] = []
    for row in frequency_report.clade_frequencies:
        if row.clade in conflict_clades:
            conflicting_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="conflict-prone",
                    rationale="clade reaches high posterior frequency but is incompatible with another high-credibility clade",
                )
            )
            continue
        if row.frequency >= robust_threshold:
            robust_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="robust",
                    rationale="clade remains near-fixed across the posterior tree set",
                )
            )
            continue
        if uncertain_min_frequency <= row.frequency <= uncertain_max_frequency:
            uncertain_clades.append(
                UncertaintyAwareCladeConclusion(
                    clade=row.clade,
                    frequency=row.frequency,
                    conclusion="uncertain",
                    rationale="clade holds intermediate support and should not anchor strong biological interpretation",
                )
            )
    robust_clades.sort(key=lambda row: (-row.frequency, row.clade))
    uncertain_clades.sort(key=lambda row: (-row.frequency, row.clade))
    conflicting_clades.sort(key=lambda row: (-row.frequency, row.clade))
    return UncertaintyAwareConclusionSummaryReport(
        path=path,
        tree_count=frequency_report.tree_count,
        robust_clade_count=len(robust_clades),
        uncertain_clade_count=len(uncertain_clades),
        conflicting_clade_count=len(conflicting_clades),
        robust_clades=robust_clades,
        uncertain_clades=uncertain_clades,
        conflicting_clades=conflicting_clades,
    )


def write_uncertainty_conclusion_table(path: Path, report: UncertaintyAwareConclusionSummaryReport) -> Path:
    """Write robust, uncertain, and conflict-prone clade conclusions as a TSV table."""
    rows = [
        *report.robust_clades,
        *report.uncertain_clades,
        *report.conflicting_clades,
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["clade", "frequency", "conclusion", "rationale"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "clade": row.clade,
                    "frequency": format(row.frequency, ".15g"),
                    "conclusion": row.conclusion,
                    "rationale": row.rationale,
                }
            )
    return path


def compare_posterior_tree_sets(left_path: Path, right_path: Path) -> PosteriorTreeSetComparisonReport:
    """Compare two tree sets over shared taxa, clade support, and cross-set topology distance."""
    _, _, left_trees = _require_tree_set(left_path)
    _, _, right_trees = _require_tree_set(right_path)
    left_taxa = _validate_same_taxa(left_trees)
    right_taxa = _validate_same_taxa(right_trees)
    if left_taxa != right_taxa:
        raise InvalidAlignmentError("posterior tree-set comparison requires identical taxon sets across both inputs")

    shared_taxa = set(left_taxa)
    left_counts = _clade_counts(left_trees, shared_taxa)
    right_counts = _clade_counts(right_trees, shared_taxa)
    all_clades = left_counts.keys() | right_counts.keys()
    deltas = [
        CladeFrequencyDelta(
            clade=_format_clade(clade),
            left_frequency=round(left_counts.get(clade, 0) / len(left_trees), 15),
            right_frequency=round(right_counts.get(clade, 0) / len(right_trees), 15),
            delta=round((right_counts.get(clade, 0) / len(right_trees)) - (left_counts.get(clade, 0) / len(left_trees)), 15),
        )
        for clade in sorted(all_clades, key=_format_clade)
    ]
    comparisons = [_tree_distance(left, right, shared_taxa) for left in left_trees for right in right_trees]
    left_topologies = {_rooted_topology_id(tree, shared_taxa) for tree in left_trees}
    right_topologies = {_rooted_topology_id(tree, shared_taxa) for tree in right_trees}
    return PosteriorTreeSetComparisonReport(
        left_path=left_path,
        right_path=right_path,
        shared_taxa=left_taxa,
        left_tree_count=len(left_trees),
        right_tree_count=len(right_trees),
        left_rooted_topology_count=len(left_topologies),
        right_rooted_topology_count=len(right_topologies),
        shared_rooted_topology_count=len(left_topologies & right_topologies),
        mean_between_set_robinson_foulds=round(sum(distance for distance, _ in comparisons) / len(comparisons), 15),
        mean_between_set_normalized_robinson_foulds=round(
            sum(normalized for _, normalized in comparisons) / len(comparisons),
            15,
        ),
        clade_frequency_deltas=deltas,
    )


def compare_bootstrap_and_posterior_uncertainty(
    bootstrap_tree_path: Path,
    posterior_tree_set_path: Path,
) -> BootstrapPosteriorSupportComparisonReport:
    """Compare bootstrap support on one summary tree against posterior clade frequencies from a tree set."""
    bootstrap_tree = _require_tree(bootstrap_tree_path)
    posterior_report = compute_clade_frequency_table(posterior_tree_set_path)
    shared_taxa = set(bootstrap_tree.tip_names)
    posterior_taxa = set(posterior_report.shared_taxa)
    if shared_taxa != posterior_taxa:
        raise InvalidAlignmentError("bootstrap versus posterior comparison requires identical taxon sets")
    bootstrap_nodes = _informative_clade_nodes(bootstrap_tree, shared_taxa)
    bootstrap_support_by_clade = {
        _format_clade(clade): _parse_support_label(node.name)
        for clade, node in bootstrap_nodes.items()
    }
    posterior_frequency_by_clade = {
        row.clade: row.frequency
        for row in posterior_report.clade_frequencies
    }
    all_clades = sorted(set(bootstrap_support_by_clade) | set(posterior_frequency_by_clade))
    rows: list[BootstrapPosteriorCladeComparison] = []
    for clade in all_clades:
        bootstrap_support = bootstrap_support_by_clade.get(clade)
        posterior_frequency = posterior_frequency_by_clade.get(clade)
        absolute_delta = None
        if bootstrap_support is not None and posterior_frequency is not None:
            absolute_delta = abs(bootstrap_support - posterior_frequency)
        agreement = _support_agreement_label(bootstrap_support, posterior_frequency, absolute_delta)
        rows.append(
            BootstrapPosteriorCladeComparison(
                clade=clade,
                bootstrap_support=bootstrap_support,
                posterior_frequency=posterior_frequency,
                absolute_delta=absolute_delta,
                agreement=agreement,
            )
        )
    return BootstrapPosteriorSupportComparisonReport(
        bootstrap_tree_path=bootstrap_tree_path,
        posterior_tree_set_path=posterior_tree_set_path,
        posterior_tree_count=posterior_report.tree_count,
        shared_taxa=posterior_report.shared_taxa,
        high_conflict_clade_count=sum(1 for row in rows if row.agreement == "strong_conflict"),
        rows=rows,
    )


def benchmark_tree_set_uncertainty(
    *,
    tree_counts: list[int] | None = None,
    taxon_counts: list[int] | None = None,
    replicates: int = 1,
    seed: int = 1,
) -> TreeSetScalingBenchmarkReport:
    """Benchmark core posterior-uncertainty summaries across tree-count and taxon-count scaling."""
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
                        temp_root / f"trees-{tree_count}-taxa-{taxon_count}-replicate-{replicate}.nwk",
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
        warnings.append("tree set is large enough to merit explicit storage and review planning")
    elif summary.tree_count >= 250 or file_size_megabytes >= 1.0:
        risk_level = "moderate"
        warnings.append("tree set is large enough that reviewer-facing summaries should be preferred over raw browsing")
    if summary.rooted_topology_count >= 100:
        warnings.append("high topology diversity increases the cost of manually inspecting individual trees")
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
    intervals = thinning_intervals or [2, 5, 10]
    if any(interval < 1 for interval in intervals):
        raise ValueError("all thinning intervals must be at least 1")
    baseline_summary = load_tree_set(path)
    baseline_clusters = cluster_trees_by_topology(path)
    baseline_conclusions = summarize_uncertainty_aware_conclusions(path)
    baseline_dominant = 0.0 if not baseline_clusters.clusters else baseline_clusters.clusters[0].frequency
    rows: list[TreeSetThinningSensitivityRow] = []
    warnings: list[str] = []
    temp_root = Path(tempfile.mkdtemp(prefix="bijux-tree-set-thinning-"))
    try:
        _, _, trees = _require_tree_set(path)
        for interval in sorted(set(intervals)):
            retained = trees[::interval]
            retained_path = write_tree_set(temp_root / f"thinned-{interval}.nwk", retained)
            summary = load_tree_set(retained_path)
            clusters = cluster_trees_by_topology(retained_path)
            conclusions = summarize_uncertainty_aware_conclusions(retained_path)
            comparison = compare_posterior_tree_sets(path, retained_path)
            dominant = 0.0 if not clusters.clusters else clusters.clusters[0].frequency
            row_warnings: list[str] = []
            if comparison.shared_rooted_topology_count < baseline_summary.rooted_topology_count:
                row_warnings.append("thinning drops one or more rooted topology modes observed in the full tree set")
            if abs(dominant - baseline_dominant) >= 0.2:
                row_warnings.append("thinning changes the dominant topology frequency materially")
            row = TreeSetThinningSensitivityRow(
                thinning_interval=interval,
                retained_tree_count=summary.tree_count,
                retained_fraction=round(summary.tree_count / baseline_summary.tree_count, 15),
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
        informative_clade_count = len(_informative_clades(tree, set(summary.shared_taxa)))
        row_warnings: list[str] = []
        if informative_clade_count == 0:
            row_warnings.append("threshold collapses all informative internal clades from the consensus summary")
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
    thinning = assess_tree_set_thinning_sensitivity(path, thinning_intervals=thinning_intervals)
    consensus = compare_consensus_thresholds(path, thresholds=consensus_thresholds)
    conclusions = summarize_uncertainty_aware_conclusions(path)
    checks = [
        TreeSetMaturityGateCheck(
            name="shared_taxa",
            satisfied=all(record.taxa == summary.records[0].taxa for record in summary.records),
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
            details="tested thinning intervals preserve the dominant conclusions" if not thinning.warnings else "; ".join(thinning.warnings),
        ),
        TreeSetMaturityGateCheck(
            name="consensus_stability",
            satisfied=not consensus.warnings,
            details="tested consensus thresholds preserve one topology summary" if not consensus.warnings else "; ".join(consensus.warnings),
        ),
        TreeSetMaturityGateCheck(
            name="uncertainty_summary",
            satisfied=summary.tree_count >= 2 and (conclusions.robust_clade_count + conclusions.uncertain_clade_count + conclusions.conflicting_clade_count) >= 1,
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
    warnings = sorted(dict.fromkeys(storage.warnings + thinning.warnings + consensus.warnings))
    return TreeSetMaturityGateReport(
        path=path,
        decision=decision,
        checks=checks,
        warnings=warnings,
    )


def _require_tree(path: Path) -> PhyloTree:
    if not path.exists():
        raise FileNotFoundError(f"tree file not found: {path}")
    source_format = detect_tree_format(path)
    bio_tree = Phylo.read(path, source_format)
    return tree_from_biophylo(bio_tree, source_format=source_format)


def _parse_support_label(value: str | None) -> float | None:
    if value is None:
        return None
    text = value.strip()
    if not text:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    return round(parsed / 100.0, 15) if parsed > 1.0 else round(parsed, 15)


def _support_agreement_label(
    bootstrap_support: float | None,
    posterior_frequency: float | None,
    absolute_delta: float | None,
) -> str:
    if bootstrap_support is None and posterior_frequency is None:
        return "not_observed"
    if bootstrap_support is None or posterior_frequency is None:
        return "method_specific"
    if absolute_delta is None:
        return "not_comparable"
    if absolute_delta >= 0.35:
        return "strong_conflict"
    if absolute_delta >= 0.15:
        return "moderate_difference"
    return "broad_agreement"
