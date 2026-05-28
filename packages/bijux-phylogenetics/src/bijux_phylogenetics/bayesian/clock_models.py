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

CLOCK_RATE_MODEL_FAMILIES = ("strict-clock", "relaxed-lognormal")
RELAXED_CLOCK_RATE_POLICIES = ("independent", "autocorrelated")


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


@dataclass(frozen=True, slots=True)
class RelaxedLognormalClockModel:
    """One validated relaxed lognormal clock model over dated and substitution trees."""

    family: str
    rate_policy: str
    mean_clock_rate: float
    log_standard_deviation: float


@dataclass(frozen=True, slots=True)
class RelaxedLognormalClockBranchRow:
    """One branch-level relaxed clock rate prior contribution row."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    anchor_branch_id: str | None
    dated_time_duration: float
    observed_substitution_branch_length: float
    anchor_rate: float
    branch_rate: float
    expected_substitution_branch_length: float
    branch_rate_deviation: float
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class RelaxedLognormalClockEvaluationReport:
    """One relaxed lognormal clock prior evaluation report."""

    family: str
    rate_policy: str
    mean_clock_rate: float
    log_standard_deviation: float
    dated_tree_newick: str
    substitution_tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    branch_count: int
    minimum_branch_rate: float
    maximum_branch_rate: float
    total_log_prior: float
    branch_rows: list[RelaxedLognormalClockBranchRow]


@dataclass(frozen=True, slots=True)
class _MatchedClockBranch:
    """One matched rooted branch pair between dated and substitution trees."""

    branch_id: str
    parent_branch_id: str | None
    child_name: str | None
    descendant_taxa: list[str]
    dated_time_duration: float
    observed_substitution_branch_length: float


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


def build_relaxed_lognormal_clock_model(
    *,
    rate_policy: str,
    mean_clock_rate: float,
    log_standard_deviation: float,
) -> RelaxedLognormalClockModel:
    """Build one relaxed lognormal clock model under one explicit rate policy."""
    normalized_rate_policy = rate_policy.strip().casefold()
    if normalized_rate_policy not in RELAXED_CLOCK_RATE_POLICIES:
        raise PhylogeneticsError(
            "relaxed lognormal clock model requires a supported rate policy",
            code="relaxed_lognormal_clock_model_invalid_rate_policy",
            details={
                "rate_policy": rate_policy,
                "allowed_rate_policies": list(RELAXED_CLOCK_RATE_POLICIES),
            },
        )
    if not math.isfinite(mean_clock_rate) or mean_clock_rate <= 0.0:
        raise PhylogeneticsError(
            "relaxed lognormal clock model requires a strictly positive finite mean clock rate",
            code="relaxed_lognormal_clock_model_invalid_mean_clock_rate",
            details={"mean_clock_rate": mean_clock_rate},
        )
    if not math.isfinite(log_standard_deviation) or log_standard_deviation <= 0.0:
        raise PhylogeneticsError(
            "relaxed lognormal clock model requires a strictly positive finite log standard deviation",
            code="relaxed_lognormal_clock_model_invalid_log_standard_deviation",
            details={"log_standard_deviation": log_standard_deviation},
        )
    return RelaxedLognormalClockModel(
        family="relaxed-lognormal",
        rate_policy=normalized_rate_policy,
        mean_clock_rate=float(format(mean_clock_rate, ".15g")),
        log_standard_deviation=float(format(log_standard_deviation, ".15g")),
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


def evaluate_relaxed_lognormal_clock_tree_log_prior(
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    rate_model: RelaxedLognormalClockModel,
) -> RelaxedLognormalClockEvaluationReport:
    """Evaluate one rooted tree pair under one relaxed lognormal clock prior."""
    if rate_model.family != "relaxed-lognormal":
        raise PhylogeneticsError(
            "relaxed lognormal clock model family is unsupported",
            code="relaxed_lognormal_clock_model_family_invalid",
            details={"family": rate_model.family},
        )
    matched_branches = _collect_matched_clock_branches(
        substitution_tree,
        dated_tree,
        workflow_name="relaxed lognormal clock model",
        require_positive_dated_durations=True,
    )
    branch_rate_by_branch_id: dict[str, float] = {}
    branch_rows: list[RelaxedLognormalClockBranchRow] = []
    for matched_branch in matched_branches:
        branch_rate = (
            matched_branch.observed_substitution_branch_length
            / matched_branch.dated_time_duration
        )
        if rate_model.rate_policy == "independent":
            anchor_branch_id = None
            anchor_rate = rate_model.mean_clock_rate
        else:
            anchor_branch_id = matched_branch.parent_branch_id
            anchor_rate = (
                rate_model.mean_clock_rate
                if anchor_branch_id is None
                else branch_rate_by_branch_id[anchor_branch_id]
            )
        log_prior_contribution = _relaxed_clock_lognormal_log_density(
            branch_rate=branch_rate,
            anchor_rate=anchor_rate,
            log_standard_deviation=rate_model.log_standard_deviation,
        )
        branch_rows.append(
            RelaxedLognormalClockBranchRow(
                branch_id=matched_branch.branch_id,
                child_name=matched_branch.child_name,
                descendant_taxa=matched_branch.descendant_taxa,
                anchor_branch_id=anchor_branch_id,
                dated_time_duration=float(
                    format(matched_branch.dated_time_duration, ".15g")
                ),
                observed_substitution_branch_length=float(
                    format(
                        matched_branch.observed_substitution_branch_length,
                        ".15g",
                    )
                ),
                anchor_rate=float(format(anchor_rate, ".15g")),
                branch_rate=float(format(branch_rate, ".15g")),
                expected_substitution_branch_length=float(
                    format(anchor_rate * matched_branch.dated_time_duration, ".15g")
                ),
                branch_rate_deviation=float(
                    format(branch_rate - anchor_rate, ".15g")
                ),
                log_prior_contribution=log_prior_contribution,
            )
        )
        branch_rate_by_branch_id[matched_branch.branch_id] = branch_rate

    total_log_prior = sum(row.log_prior_contribution for row in branch_rows)
    branch_rates = [row.branch_rate for row in branch_rows]
    return RelaxedLognormalClockEvaluationReport(
        family=rate_model.family,
        rate_policy=rate_model.rate_policy,
        mean_clock_rate=rate_model.mean_clock_rate,
        log_standard_deviation=rate_model.log_standard_deviation,
        dated_tree_newick=dumps_newick(dated_tree),
        substitution_tree_newick=dumps_newick(substitution_tree),
        taxa=sorted(substitution_tree.tip_names),
        tip_count=substitution_tree.tip_count,
        internal_node_count=substitution_tree.internal_node_count,
        branch_count=len(branch_rows),
        minimum_branch_rate=min(branch_rates),
        maximum_branch_rate=max(branch_rates),
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


def _collect_matched_clock_branches(
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    *,
    workflow_name: str,
    require_positive_dated_durations: bool,
) -> list[_MatchedClockBranch]:
    _validate_rooted_tree(
        substitution_tree,
        workflow_name=workflow_name,
        tree_role="substitution",
    )
    _validate_rooted_tree(
        dated_tree,
        workflow_name=workflow_name,
        tree_role="dated",
    )
    _require_complete_branch_lengths(
        substitution_tree,
        message=f"{workflow_name} requires complete substitution branch lengths",
    )
    _require_complete_branch_lengths(
        dated_tree,
        message=f"{workflow_name} requires complete dated branch durations",
    )
    _require_nonnegative_branch_lengths(
        substitution_tree,
        message=f"{workflow_name} requires non-negative substitution branch lengths",
    )
    _require_nonnegative_branch_lengths(
        dated_tree,
        message=f"{workflow_name} requires non-negative dated branch durations",
    )
    _require_matching_topology(substitution_tree, dated_tree)

    dated_branch_by_signature = _branch_lookup_by_descendant_taxa(dated_tree)
    matched_branches: list[_MatchedClockBranch] = []
    for parent, child in substitution_tree.iter_edges():
        if child.node_id is None:
            raise PhylogeneticsError(
                f"{workflow_name} requires stable branch ids on substitution trees",
                code="clock_rate_model_requires_stable_branch_ids",
            )
        parent_branch_id = None if parent is substitution_tree.root else parent.node_id
        dated_child = dated_branch_by_signature[tuple(child.descendant_taxa)]
        dated_time_duration = float(dated_child.branch_length or 0.0)
        if require_positive_dated_durations and dated_time_duration <= 0.0:
            raise InvalidBranchLengthError(
                f"{workflow_name} requires strictly positive dated branch durations"
            )
        matched_branches.append(
            _MatchedClockBranch(
                branch_id=child.node_id,
                parent_branch_id=parent_branch_id,
                child_name=child.name,
                descendant_taxa=child.descendant_taxa,
                dated_time_duration=dated_time_duration,
                observed_substitution_branch_length=float(child.branch_length or 0.0),
            )
        )
    return matched_branches


def _copy_with_scaled_branch_lengths(
    dated_tree: PhyloTree,
    *,
    clock_rate: float,
) -> PhyloTree:
    scaled_tree = dated_tree.copy()
    for _parent, child in scaled_tree.iter_edges():
        child.branch_length = float(child.branch_length or 0.0) * clock_rate
    return scaled_tree


def _relaxed_clock_lognormal_log_density(
    *,
    branch_rate: float,
    anchor_rate: float,
    log_standard_deviation: float,
) -> float:
    if branch_rate <= 0.0 or anchor_rate <= 0.0:
        return -math.inf
    mean_log_rate = math.log(anchor_rate) - (
        (log_standard_deviation * log_standard_deviation) / 2.0
    )
    log_rate = math.log(branch_rate)
    z_score = (log_rate - mean_log_rate) / log_standard_deviation
    return (
        -math.log(branch_rate)
        - math.log(log_standard_deviation)
        - (math.log(2.0 * math.pi) / 2.0)
        - ((z_score * z_score) / 2.0)
    )
