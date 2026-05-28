from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math

from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.node_identity import build_ape_internal_node_map
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    PhylogeneticsError,
    UnrootedTreeError,
)

YULE_TREE_PRIOR_FAMILIES = ("crown-conditioned-yule",)


@dataclass(frozen=True, slots=True)
class YuleTreePriorModel:
    """One validated Yule prior parameterization for rooted crown trees."""

    family: str
    speciation_rate: float


@dataclass(frozen=True, slots=True)
class YuleTreePriorIntervalRow:
    """One deterministic crown-tree interval contribution row."""

    interval_index: int
    older_boundary_age: float
    younger_boundary_age: float
    duration: float
    lineage_count: int
    event_count: int
    interval_log_contribution: float


@dataclass(frozen=True, slots=True)
class YuleTreePriorEvaluationReport:
    """One rooted ultrametric Yule prior evaluation report."""

    family: str
    speciation_rate: float
    tree_newick: str
    tip_count: int
    internal_node_count: int
    post_root_speciation_count: int
    root_age: float
    total_branch_length: float
    ultrametric_tolerance: float
    log_prior: float
    interval_rows: list[YuleTreePriorIntervalRow]


def build_crown_conditioned_yule_tree_prior(
    speciation_rate: float,
) -> YuleTreePriorModel:
    """Build one crown-conditioned pure-birth prior."""
    if not math.isfinite(speciation_rate) or speciation_rate <= 0.0:
        raise PhylogeneticsError(
            "Yule tree prior requires a strictly positive finite speciation rate",
            code="yule_tree_prior_invalid_speciation_rate",
            details={"speciation_rate": speciation_rate},
        )
    return YuleTreePriorModel(
        family="crown-conditioned-yule",
        speciation_rate=speciation_rate,
    )


def evaluate_yule_tree_log_prior(
    tree: PhyloTree,
    prior_model: YuleTreePriorModel,
    *,
    ultrametric_tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> YuleTreePriorEvaluationReport:
    """Evaluate one rooted ultrametric crown tree under a pure-birth Yule prior."""
    _validate_yule_tree(tree)
    tip_depth_by_label = _tip_depth_by_label(tree)
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_depth_by_label,
        tolerance=ultrametric_tolerance,
    )
    if not ultrametric_summary.ultrametric:
        raise NonUltrametricTreeError(
            "Yule tree prior requires an ultrametric tree",
            code="yule_tree_prior_requires_ultrametric_tree",
            details={
                "minimum_tip_depth": ultrametric_summary.minimum_tip_depth,
                "maximum_tip_depth": ultrametric_summary.maximum_tip_depth,
                "max_tip_depth_deviation": ultrametric_summary.max_tip_depth_deviation,
                "offending_taxa": list(ultrametric_summary.offending_taxa),
                "tolerance": ultrametric_summary.tolerance,
            },
        )

    root_age = ultrametric_summary.root_age
    branch_rows = _build_yule_interval_rows(
        tree,
        speciation_rate=prior_model.speciation_rate,
        root_age=root_age,
    )
    log_prior = sum(row.interval_log_contribution for row in branch_rows)
    return YuleTreePriorEvaluationReport(
        family=prior_model.family,
        speciation_rate=prior_model.speciation_rate,
        tree_newick=tree.to_newick(),
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        post_root_speciation_count=max(tree.internal_node_count - 1, 0),
        root_age=float(format(root_age, ".15g")),
        total_branch_length=float(format(tree.total_branch_length(), ".15g")),
        ultrametric_tolerance=ultrametric_tolerance,
        log_prior=float(format(log_prior, ".15g")),
        interval_rows=branch_rows,
    )


def _validate_yule_tree(tree: PhyloTree) -> None:
    if not _tree_is_rooted(tree):
        raise UnrootedTreeError(
            "Yule tree prior requires a rooted tree",
            code="yule_tree_prior_requires_rooted_tree",
        )
    if tree.tip_count < 2:
        raise PhylogeneticsError(
            "Yule tree prior requires at least two tips",
            code="yule_tree_prior_requires_two_or_more_tips",
        )
    for _parent, child in tree.iter_edges():
        if child.branch_length is None:
            raise InvalidBranchLengthError(
                "Yule tree prior requires complete branch lengths"
            )
        if child.branch_length < 0.0:
            raise InvalidBranchLengthError(
                "Yule tree prior requires non-negative branch lengths"
            )
    if not _tree_is_strictly_bifurcating(tree):
        raise PhylogeneticsError(
            "Yule tree prior requires a strictly bifurcating tree",
            code="yule_tree_prior_requires_strictly_bifurcating_tree",
        )


def _tree_is_strictly_bifurcating(tree: PhyloTree) -> bool:
    return all(
        len(node.children) == 2 for node in tree.iter_nodes() if not node.is_leaf()
    )


def _tree_is_rooted(tree: PhyloTree) -> bool:
    if tree.rooted is True:
        return True
    return len(tree.root.children) == 2


def _tip_depth_by_label(tree: PhyloTree) -> dict[str, float]:
    depths = _node_depth_lookup(tree)
    return {
        node.name: depths[node.node_id or ""]
        for node in tree.iter_leaves()
        if node.name is not None
    }


def _node_depth_lookup(tree: PhyloTree) -> dict[str, float]:
    root_id = tree.root.node_id or ""
    depths: dict[str, float] = {root_id: 0.0}

    def visit(node: TreeNode) -> None:
        node_id = node.node_id or ""
        base_depth = depths[node_id]
        for child in node.children:
            depths[child.node_id or ""] = base_depth + float(child.branch_length or 0.0)
            visit(child)

    visit(tree.root)
    return depths


def _build_yule_interval_rows(
    tree: PhyloTree,
    *,
    speciation_rate: float,
    root_age: float,
) -> list[YuleTreePriorIntervalRow]:
    internal_nodes = build_ape_internal_node_map(tree)
    node_depths = _node_depth_lookup(tree)
    branching_ages = [
        float(format(root_age - node_depths[node.node_id or ""], ".15g"))
        for node in internal_nodes.values()
        if node is not tree.root
    ]
    events_by_age = Counter(branching_ages)
    younger_boundaries = sorted(events_by_age, reverse=True) + [0.0]
    older_boundary = float(format(root_age, ".15g"))
    lineage_count = 2
    rows: list[YuleTreePriorIntervalRow] = []
    for interval_index, younger_boundary in enumerate(younger_boundaries, start=1):
        duration = older_boundary - younger_boundary
        event_count = events_by_age.get(younger_boundary, 0)
        interval_log_contribution = (
            (event_count * math.log(speciation_rate))
            - (lineage_count * speciation_rate * duration)
        )
        rows.append(
            YuleTreePriorIntervalRow(
                interval_index=interval_index,
                older_boundary_age=float(format(older_boundary, ".15g")),
                younger_boundary_age=float(format(younger_boundary, ".15g")),
                duration=float(format(duration, ".15g")),
                lineage_count=lineage_count,
                event_count=event_count,
                interval_log_contribution=float(
                    format(interval_log_contribution, ".15g")
                ),
            )
        )
        lineage_count += event_count
        older_boundary = younger_boundary
    return rows
