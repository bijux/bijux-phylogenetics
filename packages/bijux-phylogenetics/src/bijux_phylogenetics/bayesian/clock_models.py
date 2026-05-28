from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import math

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
    UnrootedTreeError,
)

CLOCK_RATE_MODEL_FAMILIES = ("strict-clock",)


@dataclass(frozen=True, slots=True)
class StrictClockRateModel:
    """One validated strict-clock rate model over dated and substitution trees."""

    family: str
    global_clock_rate: float
    branch_length_tolerance: float


@dataclass(frozen=True, slots=True)
class StrictClockRateBranchRow:
    """One branch-level strict-clock rate comparison row."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    dated_time_duration: float
    expected_substitution_branch_length: float
    observed_substitution_branch_length: float
    branch_length_residual: float
    exact_match: bool
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class StrictClockRateModelEvaluationReport:
    """One strict-clock rate evaluation report over a matched rooted tree pair."""

    family: str
    global_clock_rate: float
    branch_length_tolerance: float
    dated_tree_newick: str
    substitution_tree_newick: str
    expected_substitution_tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    branch_count: int
    exact_match_count: int
    mismatch_count: int
    total_log_prior: float
    branch_rows: list[StrictClockRateBranchRow]


def build_strict_clock_rate_model(
    *,
    global_clock_rate: float,
    branch_length_tolerance: float = 1e-12,
) -> StrictClockRateModel:
    """Build one strict-clock rate model with one global rate."""
    if not math.isfinite(global_clock_rate) or global_clock_rate <= 0.0:
        raise PhylogeneticsError(
            "strict-clock rate model requires a strictly positive finite global clock rate",
            code="strict_clock_rate_model_invalid_global_clock_rate",
            details={"global_clock_rate": global_clock_rate},
        )
    if not math.isfinite(branch_length_tolerance) or branch_length_tolerance < 0.0:
        raise PhylogeneticsError(
            "strict-clock rate model requires a non-negative finite branch-length tolerance",
            code="strict_clock_rate_model_invalid_branch_length_tolerance",
            details={"branch_length_tolerance": branch_length_tolerance},
        )
    return StrictClockRateModel(
        family="strict-clock",
        global_clock_rate=float(format(global_clock_rate, ".15g")),
        branch_length_tolerance=float(format(branch_length_tolerance, ".15g")),
    )


def evaluate_strict_clock_tree_log_prior(
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    rate_model: StrictClockRateModel,
) -> StrictClockRateModelEvaluationReport:
    """Evaluate one rooted substitution tree against one rooted dated tree under one global clock rate."""
    if rate_model.family != "strict-clock":
        raise PhylogeneticsError(
            "strict-clock rate model family is unsupported",
            code="strict_clock_rate_model_family_invalid",
            details={
                "family": rate_model.family,
                "allowed_families": list(CLOCK_RATE_MODEL_FAMILIES),
            },
        )
    _validate_rooted_tree(
        substitution_tree,
        workflow_name="strict-clock rate model",
        tree_role="substitution",
    )
    _validate_rooted_tree(
        dated_tree,
        workflow_name="strict-clock rate model",
        tree_role="dated",
    )
    _require_complete_branch_lengths(
        substitution_tree,
        message=(
            "strict-clock rate model requires complete substitution branch lengths"
        ),
    )
    _require_complete_branch_lengths(
        dated_tree,
        message="strict-clock rate model requires complete dated branch durations",
    )
    _require_nonnegative_branch_lengths(
        substitution_tree,
        message=(
            "strict-clock rate model requires non-negative substitution branch lengths"
        ),
    )
    _require_nonnegative_branch_lengths(
        dated_tree,
        message="strict-clock rate model requires non-negative dated branch durations",
    )
    _require_matching_topology(substitution_tree, dated_tree)

    expected_substitution_tree = _copy_with_scaled_branch_lengths(
        dated_tree,
        clock_rate=rate_model.global_clock_rate,
    )
    dated_branch_by_signature = _branch_lookup_by_descendant_taxa(dated_tree)
    branch_rows: list[StrictClockRateBranchRow] = []
    exact_match_count = 0
    for _parent, child in substitution_tree.iter_edges():
        if child.node_id is None:
            raise PhylogeneticsError(
                "strict-clock rate model requires stable branch ids on substitution trees",
                code="strict_clock_rate_model_requires_stable_branch_ids",
            )
        dated_child = dated_branch_by_signature[tuple(child.descendant_taxa)]
        dated_time_duration = float(dated_child.branch_length or 0.0)
        observed_substitution_branch_length = float(child.branch_length or 0.0)
        expected_substitution_branch_length = (
            dated_time_duration * rate_model.global_clock_rate
        )
        branch_length_residual = (
            observed_substitution_branch_length - expected_substitution_branch_length
        )
        exact_match = math.isclose(
            observed_substitution_branch_length,
            expected_substitution_branch_length,
            rel_tol=0.0,
            abs_tol=rate_model.branch_length_tolerance,
        )
        if exact_match:
            exact_match_count += 1
        branch_rows.append(
            StrictClockRateBranchRow(
                branch_id=child.node_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                dated_time_duration=float(format(dated_time_duration, ".15g")),
                expected_substitution_branch_length=float(
                    format(expected_substitution_branch_length, ".15g")
                ),
                observed_substitution_branch_length=float(
                    format(observed_substitution_branch_length, ".15g")
                ),
                branch_length_residual=float(format(branch_length_residual, ".15g")),
                exact_match=exact_match,
                log_prior_contribution=0.0 if exact_match else -math.inf,
            )
        )

    total_log_prior = 0.0 if exact_match_count == len(branch_rows) else -math.inf
    return StrictClockRateModelEvaluationReport(
        family=rate_model.family,
        global_clock_rate=rate_model.global_clock_rate,
        branch_length_tolerance=rate_model.branch_length_tolerance,
        dated_tree_newick=dumps_newick(dated_tree),
        substitution_tree_newick=dumps_newick(substitution_tree),
        expected_substitution_tree_newick=dumps_newick(expected_substitution_tree),
        taxa=sorted(substitution_tree.tip_names),
        tip_count=substitution_tree.tip_count,
        internal_node_count=substitution_tree.internal_node_count,
        branch_count=len(branch_rows),
        exact_match_count=exact_match_count,
        mismatch_count=len(branch_rows) - exact_match_count,
        total_log_prior=total_log_prior,
        branch_rows=branch_rows,
    )


def _validate_rooted_tree(
    tree: PhyloTree,
    *,
    workflow_name: str,
    tree_role: str,
) -> None:
    if tree.rooted is not True:
        raise UnrootedTreeError(
            f"{workflow_name} requires one rooted {tree_role} tree",
            code="strict_clock_rate_model_requires_rooted_tree",
        )


def _require_complete_branch_lengths(tree: PhyloTree, *, message: str) -> None:
    if any(branch_length is None for branch_length in tree.branch_lengths()):
        raise InvalidBranchLengthError(message)


def _require_nonnegative_branch_lengths(tree: PhyloTree, *, message: str) -> None:
    if any(
        float(branch_length or 0.0) < 0.0 for branch_length in tree.branch_lengths()
    ):
        raise InvalidBranchLengthError(message)


def _require_matching_topology(
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
) -> None:
    substitution_taxa = sorted(substitution_tree.tip_names)
    dated_taxa = sorted(dated_tree.tip_names)
    if substitution_taxa != dated_taxa:
        raise PhylogeneticsError(
            "strict-clock rate model requires matching tree taxa",
            code="strict_clock_rate_model_requires_matching_tree_taxa",
        )
    substitution_signatures = set(_branch_lookup_by_descendant_taxa(substitution_tree))
    dated_signatures = set(_branch_lookup_by_descendant_taxa(dated_tree))
    if substitution_signatures != dated_signatures:
        raise PhylogeneticsError(
            "strict-clock rate model requires identical rooted topology between substitution and dated trees",
            code="strict_clock_rate_model_requires_identical_rooted_topology",
        )


def _branch_lookup_by_descendant_taxa(
    tree: PhyloTree,
) -> Mapping[tuple[str, ...], TreeNode]:
    return {
        tuple(child.descendant_taxa): child for _parent, child in tree.iter_edges()
    }


def _copy_with_scaled_branch_lengths(
    dated_tree: PhyloTree,
    *,
    clock_rate: float,
) -> PhyloTree:
    scaled_tree = dated_tree.copy()
    for _parent, child in scaled_tree.iter_edges():
        child.branch_length = float(child.branch_length or 0.0) * clock_rate
    return scaled_tree
