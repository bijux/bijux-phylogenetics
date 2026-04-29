from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import tempfile

from Bio import Phylo

from bijux_phylogenetics.compare.topology import _informative_clades, compare_tree_paths
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import EngineWorkflowError, InvalidAlignmentError
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.tree_set import compare_posterior_tree_sets, load_tree_set


@dataclass(slots=True)
class MaximumCladeCredibilityTreeReport:
    source_path: Path
    filtered_tree_set_path: Path
    total_tree_count: int
    burnin_fraction: float
    burnin_tree_count: int
    kept_tree_count: int
    rooted_topology_count: int
    shared_taxa: list[str]
    selected_tree_index: int
    clade_credibility_score: float
    mcc_newick: str


@dataclass(slots=True)
class PosteriorTreeSetThinningReport:
    source_path: Path
    output_path: Path
    total_tree_count: int
    burnin_fraction: float
    burnin_tree_count: int
    pre_thinning_tree_count: int
    thinning_interval: int
    retained_tree_count: int
    retained_indices: list[int]


def summarize_maximum_clade_credibility_tree(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> tuple[PhyloTree, MaximumCladeCredibilityTreeReport]:
    """Summarize one posterior tree set by selecting the maximum clade credibility tree."""
    filtered = _filter_tree_set(tree_set_path, burnin_fraction=burnin_fraction)
    clade_frequencies = _clade_frequency_map(filtered.trees)
    selected_tree, selected_index, score = _select_mcc_tree(filtered.trees, clade_frequencies)
    summary = load_tree_set(filtered.filtered_tree_set_path)
    report = MaximumCladeCredibilityTreeReport(
        source_path=tree_set_path,
        filtered_tree_set_path=filtered.filtered_tree_set_path,
        total_tree_count=filtered.total_tree_count,
        burnin_fraction=burnin_fraction,
        burnin_tree_count=filtered.burnin_tree_count,
        kept_tree_count=len(filtered.trees),
        rooted_topology_count=summary.rooted_topology_count,
        shared_taxa=summary.shared_taxa,
        selected_tree_index=selected_index,
        clade_credibility_score=round(score, 15),
        mcc_newick=dumps_newick(selected_tree),
    )
    return selected_tree, report


def thin_posterior_tree_set(
    tree_set_path: Path,
    output_path: Path,
    *,
    thinning_interval: int,
    burnin_fraction: float = 0.0,
) -> PosteriorTreeSetThinningReport:
    """Thin a posterior tree set after optional burn-in removal."""
    if thinning_interval < 1:
        raise ValueError(f"thinning_interval must be at least 1, got {thinning_interval}")
    filtered = _filter_tree_set(tree_set_path, burnin_fraction=burnin_fraction)
    retained = filtered.trees[::thinning_interval]
    if not retained:
        raise EngineWorkflowError("posterior thinning removed every tree")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "".join(dumps_newick(tree) + "\n" for tree in retained),
        encoding="utf-8",
    )
    retained_indices = list(range(filtered.burnin_tree_count + 1, filtered.total_tree_count + 1, thinning_interval))
    return PosteriorTreeSetThinningReport(
        source_path=tree_set_path,
        output_path=output_path,
        total_tree_count=filtered.total_tree_count,
        burnin_fraction=burnin_fraction,
        burnin_tree_count=filtered.burnin_tree_count,
        pre_thinning_tree_count=len(filtered.trees),
        thinning_interval=thinning_interval,
        retained_tree_count=len(retained),
        retained_indices=retained_indices,
    )


@dataclass(slots=True)
class BayesianRunTreeComparison:
    left_path: Path
    right_path: Path
    burnin_fraction: float
    left_mcc: MaximumCladeCredibilityTreeReport
    right_mcc: MaximumCladeCredibilityTreeReport
    tree_set_comparison: object
    mcc_topology: object


def compare_bayesian_tree_sets(
    left_tree_set_path: Path,
    right_tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> BayesianRunTreeComparison:
    """Compare two posterior tree sets after the same burn-in fraction."""
    left_tree, left_mcc = summarize_maximum_clade_credibility_tree(left_tree_set_path, burnin_fraction=burnin_fraction)
    right_tree, right_mcc = summarize_maximum_clade_credibility_tree(right_tree_set_path, burnin_fraction=burnin_fraction)
    tree_set_comparison = compare_posterior_tree_sets(
        left_mcc.filtered_tree_set_path,
        right_mcc.filtered_tree_set_path,
    )
    mcc_topology = _compare_in_memory_trees(left_tree, right_tree)
    return BayesianRunTreeComparison(
        left_path=left_tree_set_path,
        right_path=right_tree_set_path,
        burnin_fraction=burnin_fraction,
        left_mcc=left_mcc,
        right_mcc=right_mcc,
        tree_set_comparison=tree_set_comparison,
        mcc_topology=mcc_topology,
    )


@dataclass(slots=True)
class _FilteredPosteriorTreeSet:
    filtered_tree_set_path: Path
    total_tree_count: int
    burnin_tree_count: int
    trees: list[PhyloTree]


def _filter_tree_set(tree_set_path: Path, *, burnin_fraction: float) -> _FilteredPosteriorTreeSet:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(f"burnin_fraction must be between 0 and 1, got {burnin_fraction}")
    tree_format = detect_tree_format(tree_set_path)
    bio_trees = list(Phylo.parse(tree_set_path, tree_format))
    if not bio_trees:
        raise EngineWorkflowError(f"posterior tree set contains no trees: {tree_set_path}")
    burnin_tree_count = int(len(bio_trees) * burnin_fraction)
    kept_bio_trees = bio_trees[burnin_tree_count:]
    if not kept_bio_trees:
        raise EngineWorkflowError(f"posterior tree set is empty after burn-in filtering: {tree_set_path}")
    trees = [tree_from_biophylo(tree, source_format=tree_format) for tree in kept_bio_trees]
    filtered_tree_set_path = tree_set_path.with_suffix(f".burnin-{_fraction_token(burnin_fraction)}.nwk")
    filtered_tree_set_path.write_text(
        "".join(dumps_newick(tree) + "\n" for tree in trees),
        encoding="utf-8",
    )
    return _FilteredPosteriorTreeSet(
        filtered_tree_set_path=filtered_tree_set_path,
        total_tree_count=len(bio_trees),
        burnin_tree_count=burnin_tree_count,
        trees=trees,
    )


def _clade_frequency_map(trees: list[PhyloTree]) -> dict[frozenset[str], float]:
    taxa_sets = {frozenset(tree.tip_names) for tree in trees}
    if len(taxa_sets) != 1:
        raise InvalidAlignmentError("posterior tree summaries require all trees to share the exact same taxon set")
    shared_taxa = set(next(iter(taxa_sets)))
    counts: dict[frozenset[str], int] = {}
    for tree in trees:
        for clade in _informative_clades(tree, shared_taxa):
            counts[clade] = counts.get(clade, 0) + 1
    return {
        clade: count / len(trees)
        for clade, count in counts.items()
    }


def _select_mcc_tree(
    trees: list[PhyloTree],
    clade_frequencies: dict[frozenset[str], float],
) -> tuple[PhyloTree, int, float]:
    shared_taxa = set(trees[0].tip_names)
    scored = []
    for index, tree in enumerate(trees, start=1):
        score = 0.0
        for clade in _informative_clades(tree, shared_taxa):
            frequency = clade_frequencies.get(clade, 1e-12)
            score += math.log(max(frequency, 1e-12))
        scored.append((score, index, tree))
    best_score, best_index, best_tree = max(scored, key=lambda item: (item[0], -item[1]))
    return best_tree, best_index, best_score


def _fraction_token(value: float) -> str:
    return format(value, ".6f").rstrip("0").rstrip(".").replace(".", "p") or "0"


def _compare_in_memory_trees(left: PhyloTree, right: PhyloTree) -> object:
    left_path = _write_temp_tree(left, suffix="left")
    right_path = _write_temp_tree(right, suffix="right")
    try:
        return compare_tree_paths(left_path, right_path)
    finally:
        left_path.unlink(missing_ok=True)
        right_path.unlink(missing_ok=True)


def _write_temp_tree(tree: PhyloTree, *, suffix: str) -> Path:
    path = Path(tempfile.mkstemp(prefix=f"bijux-bayesian-{suffix}-", suffix=".nwk")[1])
    path.write_text(dumps_newick(tree) + "\n", encoding="utf-8")
    return path
