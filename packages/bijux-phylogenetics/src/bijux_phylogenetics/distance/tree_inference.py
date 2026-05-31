from __future__ import annotations

from pathlib import Path
import random
from statistics import median
import tempfile

from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_tree_paths,
)
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.phylo.alignment import AlignmentRecord
from bijux_phylogenetics.phylo.topology.clades import (
    informative_unrooted_splits,
    robinson_foulds_metrics,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    InvalidAlignmentError,
    InvalidDistanceMatrixError,
)
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
)

from .genetic_distance_matrix import (
    _build_alignment_distance_lookup,
    _load_alignment_for_model,
    compute_pairwise_genetic_distance_matrix,
)
from .missing_distance_policy import apply_missing_distance_policy
from .models import (
    AmbiguityPolicy,
    DistanceBootstrapReplicateRow,
    DistanceBootstrapReport,
    DistanceBootstrapSupportRow,
    DistanceBootstrapSupportSummary,
    DistanceMethodAssumptionReport,
    DistanceModel,
    DistanceTreeBuildReport,
    DistanceTreeReferenceComparisonReport,
    DistanceTreeTopologyComparison,
    GapHandlingMode,
    GeneticDistanceMatrix,
    MissingDistancePolicy,
)
from .quality import (
    assess_distance_method_assumptions_from_genetic_distance_matrix,
    inspect_distance_matrix_quality,
)
from .saturation import diagnose_distance_saturation_from_genetic_distance_matrix
from .shared import (
    _build_distance_tree_from_lookup,
    _require_supported_distance_tree_method,
)


def _block_tree_inference_on_saturated_distances(
    report: GeneticDistanceMatrix,
) -> None:
    """Reject corrected distances that are unusable before topology inference."""
    saturation_report = diagnose_distance_saturation_from_genetic_distance_matrix(
        report
    )
    if not saturation_report.blocks_tree_inference:
        return
    blocking_pairs = ", ".join(
        f"{row.left_identifier}/{row.right_identifier} ({row.warning_kind})"
        for row in saturation_report.warning_rows
        if row.blocks_tree_inference
    )
    raise InvalidAlignmentError(
        "distance tree building is blocked before tree inference because one or more corrected distances are undefined or infinite: "
        + blocking_pairs
    )


def build_distance_tree(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    """Build a distance-based tree from an aligned dataset."""
    method_policy = _require_supported_distance_tree_method(method)
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    if quality.method_assessment.decision == "blocked":
        raise InvalidAlignmentError(
            "distance tree building is blocked: "
            + "; ".join(quality.method_assessment.reasons)
        )
    report = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    _block_tree_inference_on_saturated_distances(report)
    return _build_distance_tree_from_genetic_distance_matrix(
        report,
        method=method_policy.method,
        assumptions=quality.assumptions,
        missing_distance_policy=missing_distance_policy,
    )


def _build_distance_tree_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    method: str,
    assumptions: DistanceMethodAssumptionReport,
    missing_distance_policy: MissingDistancePolicy,
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    _block_tree_inference_on_saturated_distances(report)
    if len(report.identifiers) < 2:
        raise InvalidAlignmentError("distance tree building requires at least two taxa")
    try:
        distance_lookup, missing_distance_policy_report = apply_missing_distance_policy(
            report.identifiers,
            _build_alignment_distance_lookup(report),
            policy=missing_distance_policy,
        )
    except InvalidDistanceMatrixError as error:
        raise InvalidAlignmentError(
            error.message,
            code=getattr(error, "code", None),
            details=getattr(error, "details", None),
        ) from error
    tree = _build_distance_tree_from_lookup(
        report.identifiers,
        distance_lookup,
        method=method,
    )
    return tree, DistanceTreeBuildReport(
        alignment_path=report.path,
        model=report.model,
        gap_handling=report.gap_handling,
        ambiguity_policy=report.ambiguity_policy,
        method=method,
        method_policy=_require_supported_distance_tree_method(method),
        taxon_count=len(report.identifiers),
        pair_count=len(report.pairs),
        assumptions=assumptions,
        missing_distance_policy_report=missing_distance_policy_report,
    )


def build_distance_tree_from_genetic_distance_matrix(
    report: GeneticDistanceMatrix,
    *,
    method: str,
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> tuple[PhyloTree, DistanceTreeBuildReport]:
    """Build a distance tree directly from one in-memory genetic distance matrix."""
    method_policy = _require_supported_distance_tree_method(method)
    assumptions = assess_distance_method_assumptions_from_genetic_distance_matrix(
        report
    )
    return _build_distance_tree_from_genetic_distance_matrix(
        report,
        method=method_policy.method,
        assumptions=assumptions,
        missing_distance_policy=missing_distance_policy,
    )


def compare_distance_tree_topologies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> DistanceTreeTopologyComparison:
    """Compare NJ and UPGMA topologies built from the same alignment."""
    nj_tree, _ = build_distance_tree(
        path,
        method="neighbor-joining",
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        missing_distance_policy=missing_distance_policy,
    )
    upgma_tree, _ = build_distance_tree(
        path,
        method="upgma",
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        missing_distance_policy=missing_distance_policy,
    )
    shared_taxa = set(nj_tree.tip_names) & set(upgma_tree.tip_names)
    rooted_metrics = robinson_foulds_metrics(
        nj_tree,
        upgma_tree,
        shared_taxa,
        rf_mode="rooted",
    )
    topology_equal = rooted_metrics.distance == 0
    same_unrooted_topology = informative_unrooted_splits(
        nj_tree,
        shared_taxa,
    ) == informative_unrooted_splits(
        upgma_tree,
        shared_taxa,
    )
    return DistanceTreeTopologyComparison(
        alignment_path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        shared_taxa=sorted(shared_taxa),
        nj_informative_clades=rooted_metrics.left_count,
        upgma_informative_clades=rooted_metrics.right_count,
        robinson_foulds_distance=rooted_metrics.distance,
        normalized_robinson_foulds=rooted_metrics.normalized_distance,
        topology_equal=topology_equal,
        same_unrooted_topology=same_unrooted_topology,
        same_taxa_different_rooting=topology_equal is False and same_unrooted_topology,
    )


def _bootstrap_resample_alignment_columns(
    records: list[AlignmentRecord], *, rng: random.Random
) -> tuple[list[int], list[AlignmentRecord]]:
    sampled_site_indices = [
        rng.randrange(len(records[0].sequence)) for _ in range(len(records[0].sequence))
    ]
    return sampled_site_indices, [
        AlignmentRecord(
            identifier=record.identifier,
            sequence="".join(
                record.sequence[position] for position in sampled_site_indices
            ),
        )
        for record in records
    ]


def _write_bootstrap_alignment(path: Path, records: list[AlignmentRecord]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.identifier}")
        lines.append(record.sequence)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def bootstrap_distance_trees(
    path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    replicates: int = 100,
    seed: int = 1,
) -> tuple[list[PhyloTree], DistanceBootstrapReport]:
    """Bootstrap a distance tree by resampling alignment sites with replacement."""
    method_policy = _require_supported_distance_tree_method(method)
    if replicates < 1:
        raise ValueError(f"replicates must be at least 1, got {replicates}")
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    if quality.method_assessment.decision == "blocked":
        raise InvalidAlignmentError(
            "distance bootstrap is blocked: "
            + "; ".join(quality.method_assessment.reasons)
        )
    records, _ = _load_alignment_for_model(path, model=model)
    rng = random.Random(seed)  # nosec B311
    temp_dir = Path(tempfile.mkdtemp(prefix="bijux-distance-bootstrap-"))
    trees: list[PhyloTree] = []
    replicate_rows: list[DistanceBootstrapReplicateRow] = []
    for index in range(replicates):
        sampled_site_indices, replicate_records = _bootstrap_resample_alignment_columns(
            records,
            rng=rng,
        )
        replicate_path = temp_dir / f"replicate-{index + 1}.fasta"
        _write_bootstrap_alignment(replicate_path, replicate_records)
        tree, _ = build_distance_tree(
            replicate_path,
            method=method_policy.method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )
        trees.append(tree)
        replicate_rows.append(
            DistanceBootstrapReplicateRow(
                replicate_index=index + 1,
                sampled_site_indices=sampled_site_indices,
                tree_newick=dumps_newick(tree),
            )
        )
    tree_set_path = temp_dir / "bootstrap.trees"
    from bijux_phylogenetics.simulation import write_tree_set

    write_tree_set(tree_set_path, trees)
    consensus_tree, _consensus = compute_consensus_tree(tree_set_path)
    support = compute_clade_frequency_table(tree_set_path)
    return trees, DistanceBootstrapReport(
        alignment_path=path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        method=method_policy.method,
        replicates=replicates,
        seed=seed,
        tree_count=len(trees),
        consensus_newick=dumps_newick(consensus_tree),
        replicate_rows=replicate_rows,
        support=[
            DistanceBootstrapSupportRow(
                clade=row.clade,
                tree_count=row.tree_count,
                frequency=row.frequency,
            )
            for row in support.clade_frequencies
        ],
    )


def write_distance_bootstrap_support(
    path: Path, report: DistanceBootstrapReport
) -> Path:
    """Write bootstrap clade support as TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["clade\ttree_count\tfrequency"]
    lines.extend(
        f"{row.clade}\t{row.tree_count}\t{format(row.frequency, '.15g')}"
        for row in report.support
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_distance_bootstrap_draws(path: Path, report: DistanceBootstrapReport) -> Path:
    """Write one deterministic bootstrap draw ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["replicate_index\tsampled_site_indices\ttree_newick"]
    lines.extend(
        (
            f"{row.replicate_index}\t"
            + ",".join(str(index) for index in row.sampled_site_indices)
            + f"\t{row.tree_newick}"
        )
        for row in report.replicate_rows
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def summarize_distance_bootstrap_support(
    report: DistanceBootstrapReport,
    *,
    weak_frequency_threshold: float = 0.5,
) -> DistanceBootstrapSupportSummary:
    """Summarize bootstrap clade frequencies for reviewer-facing reporting."""
    frequencies = sorted(row.frequency for row in report.support)
    weak_clade_count = sum(
        1 for row in report.support if row.frequency < weak_frequency_threshold
    )
    warnings: list[str] = []
    if weak_clade_count:
        warnings.append(
            "one or more consensus clades remain weakly supported across bootstrap replicates"
        )
    if not frequencies:
        warnings.append(
            "bootstrap replicates did not yield any informative internal clades"
        )
    return DistanceBootstrapSupportSummary(
        alignment_path=report.alignment_path,
        method=report.method,
        model=report.model,
        gap_handling=report.gap_handling,
        ambiguity_policy=report.ambiguity_policy,
        replicates=report.replicates,
        clade_count=len(report.support),
        minimum_frequency=None if not frequencies else min(frequencies),
        maximum_frequency=None if not frequencies else max(frequencies),
        median_frequency=None if not frequencies else median(frequencies),
        weak_clade_count=weak_clade_count,
        warnings=warnings,
    )


def compare_distance_tree_to_reference_tree(
    path: Path,
    reference_tree_path: Path,
    *,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    missing_distance_policy: MissingDistancePolicy = "reject",
) -> DistanceTreeReferenceComparisonReport:
    """Compare one built distance tree to an external ML or reviewer-supplied reference tree."""
    tree, _ = build_distance_tree(
        path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        missing_distance_policy=missing_distance_policy,
    )
    temp_dir = Path(tempfile.mkdtemp(prefix="bijux-distance-reference-"))
    built_tree_path = write_newick(temp_dir / "distance-tree.nwk", tree)
    topology = compare_tree_paths(built_tree_path, reference_tree_path)
    branch_lengths = compare_branch_lengths(built_tree_path, reference_tree_path)
    warnings: list[str] = []
    if not topology.topology_equal:
        warnings.append(
            "distance tree disagrees topologically with the supplied reference tree"
        )
    if topology.same_unrooted_topology and not topology.topology_equal:
        warnings.append(
            "distance tree matches the reference on unrooted splits but differs in rooting"
        )
    if topology.same_topology_different_branch_lengths:
        warnings.append(
            "distance tree preserves topology but shifts branch-length interpretation relative to the reference"
        )
    return DistanceTreeReferenceComparisonReport(
        alignment_path=path,
        reference_tree_path=reference_tree_path,
        method=method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        topology=topology,
        branch_lengths=branch_lengths,
        warnings=warnings,
    )
