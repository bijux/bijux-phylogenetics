from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import random
import tempfile

from Bio import Phylo

from bijux_phylogenetics.compare.topology import compare_tree_paths
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.phylo.topology.clades import informative_rooted_clades
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    EngineWorkflowError,
    InvalidAlignmentError,
)
from bijux_phylogenetics.trees import compare_posterior_tree_sets, load_tree_set


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
class PosteriorCladeAgeSummary:
    clade: str
    mean_height: float
    median_height: float
    minimum_height: float
    maximum_height: float
    lower_95_credible_interval: float
    upper_95_credible_interval: float
    tree_count: int


@dataclass(slots=True)
class PosteriorNodeAgeSummaryReport:
    source_path: Path
    filtered_tree_set_path: Path
    burnin_fraction: float
    total_tree_count: int
    kept_tree_count: int
    rows: list[PosteriorCladeAgeSummary]


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


@dataclass(slots=True)
class PosteriorTreeSubsampleEntry:
    retained_order: int
    source_tree_index: int
    post_burnin_index: int
    tree_name: str | None
    state: int | None
    generation: int | None
    rooted: bool | None
    tip_names: list[str]
    newick: str


@dataclass(slots=True)
class PosteriorTreeSubsamplingReport:
    source_path: Path
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    pre_subsampling_tree_count: int
    selection_method: str
    thinning_interval: int | None
    requested_tree_count: int | None
    random_seed: int | None
    retained_tree_count: int
    retained_source_indices: list[int]
    trees: list[PosteriorTreeSubsampleEntry]


def summarize_maximum_clade_credibility_tree(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> tuple[PhyloTree, MaximumCladeCredibilityTreeReport]:
    """Summarize one posterior tree set by selecting the maximum clade credibility tree."""
    filtered = _filter_tree_set(tree_set_path, burnin_fraction=burnin_fraction)
    clade_frequencies = _clade_frequency_map(filtered.trees)
    selected_tree, selected_index, score = _select_mcc_tree(
        filtered.trees, clade_frequencies
    )
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
    report = subsample_posterior_tree_set(
        tree_set_path,
        method="evenly-spaced",
        thinning_interval=thinning_interval,
        burnin_fraction=burnin_fraction,
    )
    write_posterior_tree_subsample(output_path, report)
    return PosteriorTreeSetThinningReport(
        source_path=tree_set_path,
        output_path=output_path,
        total_tree_count=report.total_tree_count,
        burnin_fraction=burnin_fraction,
        burnin_tree_count=report.burnin_tree_count,
        pre_thinning_tree_count=report.pre_subsampling_tree_count,
        thinning_interval=thinning_interval,
        retained_tree_count=report.retained_tree_count,
        retained_indices=report.retained_source_indices,
    )


def summarize_posterior_node_ages(
    tree_set_path: Path,
    *,
    burnin_fraction: float = 0.25,
) -> PosteriorNodeAgeSummaryReport:
    """Summarize clade heights across a posterior tree set after burn-in removal."""
    filtered = _filter_tree_set(tree_set_path, burnin_fraction=burnin_fraction)
    age_rows = _summarize_clade_heights(filtered.trees)
    return PosteriorNodeAgeSummaryReport(
        source_path=tree_set_path,
        filtered_tree_set_path=filtered.filtered_tree_set_path,
        burnin_fraction=burnin_fraction,
        total_tree_count=filtered.total_tree_count,
        kept_tree_count=len(filtered.trees),
        rows=age_rows,
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
    left_tree, left_mcc = summarize_maximum_clade_credibility_tree(
        left_tree_set_path, burnin_fraction=burnin_fraction
    )
    right_tree, right_mcc = summarize_maximum_clade_credibility_tree(
        right_tree_set_path, burnin_fraction=burnin_fraction
    )
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


@dataclass(slots=True)
class _PosteriorTreeSelectionInput:
    source_tree_index: int
    post_burnin_index: int
    tree_name: str | None
    state: int | None
    generation: int | None
    rooted: bool | None
    tip_names: list[str]
    newick: str


@dataclass(slots=True)
class _PosteriorTreeSelectionSource:
    source_path: Path
    burnin_fraction: float
    total_tree_count: int
    burnin_tree_count: int
    trees: list[_PosteriorTreeSelectionInput]


def subsample_posterior_tree_set(
    tree_set_path: Path,
    *,
    method: str,
    thinning_interval: int | None = None,
    sample_count: int | None = None,
    burnin_fraction: float = 0.0,
    random_seed: int | None = None,
) -> PosteriorTreeSubsamplingReport:
    """Subsample a generic posterior tree set after optional burn-in removal."""
    filtered = _filter_tree_set(tree_set_path, burnin_fraction=burnin_fraction)
    source = _PosteriorTreeSelectionSource(
        source_path=tree_set_path,
        burnin_fraction=burnin_fraction,
        total_tree_count=filtered.total_tree_count,
        burnin_tree_count=filtered.burnin_tree_count,
        trees=[
            _PosteriorTreeSelectionInput(
                source_tree_index=filtered.burnin_tree_count + index,
                post_burnin_index=index,
                tree_name=None,
                state=None,
                generation=None,
                rooted=tree.rooted,
                tip_names=tree.tip_names,
                newick=dumps_newick(tree),
            )
            for index, tree in enumerate(filtered.trees, start=1)
        ],
    )
    return _subsample_posterior_tree_source(
        source,
        method=method,
        thinning_interval=thinning_interval,
        sample_count=sample_count,
        random_seed=random_seed,
    )


def subsample_beast_posterior_tree_set(
    tree_set_path: Path,
    *,
    method: str,
    thinning_interval: int | None = None,
    sample_count: int | None = None,
    burnin_fraction: float = 0.0,
    random_seed: int | None = None,
) -> PosteriorTreeSubsamplingReport:
    """Subsample native BEAST posterior trees while preserving state metadata."""
    from bijux_phylogenetics.bayesian.beast.posterior_trees import (
        parse_beast_posterior_tree_samples,
    )

    report = parse_beast_posterior_tree_samples(
        tree_set_path,
        burnin_fraction=burnin_fraction,
    )
    source = _PosteriorTreeSelectionSource(
        source_path=tree_set_path,
        burnin_fraction=report.burnin_fraction,
        total_tree_count=report.total_tree_count,
        burnin_tree_count=report.burnin_tree_count,
        trees=[
            _PosteriorTreeSelectionInput(
                source_tree_index=report.burnin_tree_count + index,
                post_burnin_index=index,
                tree_name=sample.tree_name,
                state=sample.state,
                generation=None,
                rooted=sample.rooted,
                tip_names=sample.tip_names,
                newick=sample.newick,
            )
            for index, sample in enumerate(report.trees, start=1)
        ],
    )
    return _subsample_posterior_tree_source(
        source,
        method=method,
        thinning_interval=thinning_interval,
        sample_count=sample_count,
        random_seed=random_seed,
    )


def subsample_mrbayes_posterior_tree_set(
    tree_set_path: Path,
    *,
    method: str,
    thinning_interval: int | None = None,
    sample_count: int | None = None,
    burnin_fraction: float = 0.0,
    random_seed: int | None = None,
) -> PosteriorTreeSubsamplingReport:
    """Subsample native MrBayes posterior trees while preserving generation metadata."""
    from bijux_phylogenetics.bayesian.mrbayes import (
        parse_mrbayes_posterior_tree_samples,
    )

    _validate_burnin_fraction(burnin_fraction)
    parsed = parse_mrbayes_posterior_tree_samples(tree_set_path)
    burnin_tree_count = int(parsed.tree_count * burnin_fraction)
    kept_trees = parsed.trees[burnin_tree_count:]
    if not kept_trees:
        raise EngineWorkflowError(
            f"MrBayes posterior tree file is empty after burn-in filtering: {tree_set_path}"
        )
    source = _PosteriorTreeSelectionSource(
        source_path=tree_set_path,
        burnin_fraction=burnin_fraction,
        total_tree_count=parsed.tree_count,
        burnin_tree_count=burnin_tree_count,
        trees=[
            _PosteriorTreeSelectionInput(
                source_tree_index=burnin_tree_count + index,
                post_burnin_index=index,
                tree_name=sample.tree_name,
                state=None,
                generation=sample.generation,
                rooted=sample.rooted,
                tip_names=sample.tip_names,
                newick=sample.newick,
            )
            for index, sample in enumerate(kept_trees, start=1)
        ],
    )
    return _subsample_posterior_tree_source(
        source,
        method=method,
        thinning_interval=thinning_interval,
        sample_count=sample_count,
        random_seed=random_seed,
    )


def write_posterior_tree_subsample(
    path: Path, report: PosteriorTreeSubsamplingReport
) -> Path:
    """Write retained posterior trees as a normalized Newick tree set."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(f"{sample.newick}\n" for sample in report.trees),
        encoding="utf-8",
    )
    return path


def write_posterior_tree_subsample_table(
    path: Path, report: PosteriorTreeSubsamplingReport
) -> Path:
    """Write a reviewer-facing ledger for retained posterior-tree samples."""
    return write_taxon_rows(
        path,
        columns=[
            "retained_order",
            "source_tree_index",
            "post_burnin_index",
            "tree_name",
            "state",
            "generation",
            "rooted",
            "tip_count",
            "selection_method",
            "thinning_interval",
            "requested_tree_count",
            "random_seed",
            "burnin_fraction",
            "burnin_tree_count",
            "pre_subsampling_tree_count",
            "retained_tree_count",
        ],
        rows=[
            {
                "retained_order": str(sample.retained_order),
                "source_tree_index": str(sample.source_tree_index),
                "post_burnin_index": str(sample.post_burnin_index),
                "tree_name": "" if sample.tree_name is None else sample.tree_name,
                "state": "" if sample.state is None else str(sample.state),
                "generation": (
                    "" if sample.generation is None else str(sample.generation)
                ),
                "rooted": "" if sample.rooted is None else str(sample.rooted).lower(),
                "tip_count": str(len(sample.tip_names)),
                "selection_method": report.selection_method,
                "thinning_interval": (
                    ""
                    if report.thinning_interval is None
                    else str(report.thinning_interval)
                ),
                "requested_tree_count": (
                    ""
                    if report.requested_tree_count is None
                    else str(report.requested_tree_count)
                ),
                "random_seed": (
                    "" if report.random_seed is None else str(report.random_seed)
                ),
                "burnin_fraction": format(report.burnin_fraction, ".15g"),
                "burnin_tree_count": str(report.burnin_tree_count),
                "pre_subsampling_tree_count": str(report.pre_subsampling_tree_count),
                "retained_tree_count": str(report.retained_tree_count),
            }
            for sample in report.trees
        ],
    )


def _subsample_posterior_tree_source(
    source: _PosteriorTreeSelectionSource,
    *,
    method: str,
    thinning_interval: int | None,
    sample_count: int | None,
    random_seed: int | None,
) -> PosteriorTreeSubsamplingReport:
    retained_inputs, normalized_method = _select_posterior_tree_inputs(
        source.trees,
        method=method,
        thinning_interval=thinning_interval,
        sample_count=sample_count,
        random_seed=random_seed,
    )
    retained_trees = [
        PosteriorTreeSubsampleEntry(
            retained_order=index,
            source_tree_index=tree.source_tree_index,
            post_burnin_index=tree.post_burnin_index,
            tree_name=tree.tree_name,
            state=tree.state,
            generation=tree.generation,
            rooted=tree.rooted,
            tip_names=tree.tip_names,
            newick=tree.newick,
        )
        for index, tree in enumerate(retained_inputs, start=1)
    ]
    return PosteriorTreeSubsamplingReport(
        source_path=source.source_path,
        burnin_fraction=source.burnin_fraction,
        total_tree_count=source.total_tree_count,
        burnin_tree_count=source.burnin_tree_count,
        pre_subsampling_tree_count=len(source.trees),
        selection_method=normalized_method,
        thinning_interval=thinning_interval
        if normalized_method == "evenly-spaced"
        else None,
        requested_tree_count=sample_count if normalized_method == "random" else None,
        random_seed=random_seed if normalized_method == "random" else None,
        retained_tree_count=len(retained_trees),
        retained_source_indices=[tree.source_tree_index for tree in retained_trees],
        trees=retained_trees,
    )


def _select_posterior_tree_inputs(
    trees: list[_PosteriorTreeSelectionInput],
    *,
    method: str,
    thinning_interval: int | None,
    sample_count: int | None,
    random_seed: int | None,
) -> tuple[list[_PosteriorTreeSelectionInput], str]:
    if method not in {"evenly-spaced", "random"}:
        raise ValueError(
            f"method must be one of {{'evenly-spaced', 'random'}}, got {method!r}"
        )
    if method == "evenly-spaced":
        if thinning_interval is None:
            raise ValueError(
                "thinning_interval is required for evenly-spaced posterior subsampling"
            )
        if thinning_interval < 1:
            raise ValueError(
                f"thinning_interval must be at least 1, got {thinning_interval}"
            )
        if sample_count is not None:
            raise ValueError(
                "sample_count is not supported for evenly-spaced posterior subsampling"
            )
        if random_seed is not None:
            raise ValueError(
                "random_seed is only valid for random posterior subsampling"
            )
        retained = trees[::thinning_interval]
        if not retained:
            raise EngineWorkflowError("posterior thinning removed every tree")
        return retained, method
    if thinning_interval is not None:
        raise ValueError(
            "thinning_interval is only valid for evenly-spaced subsampling"
        )
    if sample_count is None:
        raise ValueError("sample_count is required for random posterior subsampling")
    if sample_count < 1:
        raise ValueError(f"sample_count must be at least 1, got {sample_count}")
    if sample_count > len(trees):
        raise ValueError(
            "sample_count cannot exceed the number of post-burn-in trees: "
            f"{sample_count} > {len(trees)}"
        )
    retained_positions = sorted(
        # Deterministic posterior subsampling is required for reproducible review.
        random.Random(random_seed).sample(  # nosec B311
            range(len(trees)),
            sample_count,
        )
    )
    return [trees[position] for position in retained_positions], method


def _filter_tree_set(
    tree_set_path: Path, *, burnin_fraction: float
) -> _FilteredPosteriorTreeSet:
    _validate_burnin_fraction(burnin_fraction)
    tree_format = detect_tree_format(tree_set_path)
    bio_trees = list(Phylo.parse(tree_set_path, tree_format))
    if not bio_trees:
        raise EngineWorkflowError(
            f"posterior tree set contains no trees: {tree_set_path}"
        )
    burnin_tree_count = int(len(bio_trees) * burnin_fraction)
    kept_bio_trees = bio_trees[burnin_tree_count:]
    if not kept_bio_trees:
        raise EngineWorkflowError(
            f"posterior tree set is empty after burn-in filtering: {tree_set_path}"
        )
    trees = [
        tree_from_biophylo(tree, source_format=tree_format) for tree in kept_bio_trees
    ]
    filtered_tree_set_path = Path(
        tempfile.mkstemp(
            prefix=f"{tree_set_path.stem}.burnin-{_fraction_token(burnin_fraction)}-",
            suffix=".nwk",
        )[1]
    )
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


def _validate_burnin_fraction(burnin_fraction: float) -> None:
    if not 0.0 <= burnin_fraction < 1.0:
        raise ValueError(
            f"burnin_fraction must be between 0 and 1, got {burnin_fraction}"
        )


def _clade_frequency_map(trees: list[PhyloTree]) -> dict[frozenset[str], float]:
    taxa_sets = {frozenset(tree.tip_names) for tree in trees}
    if len(taxa_sets) != 1:
        raise InvalidAlignmentError(
            "posterior tree summaries require all trees to share the exact same taxon set"
        )
    shared_taxa = set(next(iter(taxa_sets)))
    counts: dict[frozenset[str], int] = {}
    for tree in trees:
        for clade in informative_rooted_clades(tree, shared_taxa):
            counts[clade] = counts.get(clade, 0) + 1
    return {clade: count / len(trees) for clade, count in counts.items()}


def _select_mcc_tree(
    trees: list[PhyloTree],
    clade_frequencies: dict[frozenset[str], float],
) -> tuple[PhyloTree, int, float]:
    shared_taxa = set(trees[0].tip_names)
    scored = []
    for index, tree in enumerate(trees, start=1):
        score = 0.0
        for clade in informative_rooted_clades(tree, shared_taxa):
            frequency = clade_frequencies.get(clade, 1e-12)
            score += math.log(max(frequency, 1e-12))
        scored.append((score, index, tree))
    best_score, best_index, best_tree = max(
        scored, key=lambda item: (item[0], -item[1])
    )
    return best_tree, best_index, best_score


def _quantile(sorted_values: list[float], fraction: float) -> float:
    if len(sorted_values) == 1:
        return round(sorted_values[0], 15)
    position = max(
        0, min(len(sorted_values) - 1, int(round(fraction * (len(sorted_values) - 1))))
    )
    return round(sorted_values[position], 15)


def _summarize_clade_heights(trees: list[PhyloTree]) -> list[PosteriorCladeAgeSummary]:
    taxa_sets = {frozenset(tree.tip_names) for tree in trees}
    if len(taxa_sets) != 1:
        raise InvalidAlignmentError(
            "posterior age summaries require all trees to share the exact same taxon set"
        )
    shared_taxa = set(next(iter(taxa_sets)))
    clade_heights: dict[frozenset[str], list[float]] = {}
    for tree in trees:
        _collect_clade_heights(
            tree.root,
            shared_taxa=shared_taxa,
            current_height=0.0,
            clade_heights=clade_heights,
        )
    return [
        PosteriorCladeAgeSummary(
            clade="|".join(sorted(clade)),
            mean_height=round(sum(ordered := sorted(heights)) / len(ordered), 15),
            median_height=_quantile(ordered, 0.5),
            minimum_height=round(min(ordered), 15),
            maximum_height=round(max(ordered), 15),
            lower_95_credible_interval=_quantile(ordered, 0.025),
            upper_95_credible_interval=_quantile(ordered, 0.975),
            tree_count=len(ordered),
        )
        for clade, heights in sorted(
            clade_heights.items(), key=lambda item: (len(item[0]), sorted(item[0]))
        )
    ]


def _collect_clade_heights(
    node: TreeNode,
    *,
    shared_taxa: set[str],
    current_height: float,
    clade_heights: dict[frozenset[str], list[float]],
) -> set[str]:
    if node.is_leaf():
        return {node.name} if node.name in shared_taxa else set()
    taxa: set[str] = set()
    for child in node.children:
        branch_length = float(child.branch_length or 0.0)
        taxa.update(
            _collect_clade_heights(
                child,
                shared_taxa=shared_taxa,
                current_height=current_height + branch_length,
                clade_heights=clade_heights,
            )
        )
    if 1 < len(taxa) < len(shared_taxa):
        clade_heights.setdefault(frozenset(taxa), []).append(current_height)
    return taxa


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
