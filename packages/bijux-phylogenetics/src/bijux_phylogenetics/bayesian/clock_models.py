from __future__ import annotations

from collections.abc import Mapping
import csv
from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
    UnrootedTreeError,
)

CLOCK_RATE_MODEL_FAMILIES = ("strict-clock", "relaxed-lognormal")
RELAXED_CLOCK_RATE_POLICIES = ("independent", "autocorrelated")
LOCAL_CLOCK_RATE_MODEL_FAMILIES = ("local-clock",)
LOCAL_CLOCK_TARGET_KINDS = ("background", "branch", "clade")
_LOCAL_CLOCK_BACKGROUND_REGIME_ID = "background"


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
class LocalClockRegimeDefinition:
    """One local-clock regime selector over rooted tree branches."""

    regime_id: str
    target_kind: str
    target_label: str | None
    descendant_taxa: list[str]


@dataclass(frozen=True, slots=True)
class LocalClockRateModel:
    """One validated local-clock prior parameterization."""

    family: str
    background_clock_rate: float
    regime_clock_rates: dict[str, float]
    log_standard_deviation: float


@dataclass(frozen=True, slots=True)
class LocalClockRateBranchRow:
    """One branch-level local-clock prior contribution row."""

    branch_id: str
    child_name: str | None
    descendant_taxa: list[str]
    regime_id: str
    target_kind: str
    dated_time_duration: float
    observed_substitution_branch_length: float
    class_clock_rate: float
    branch_rate: float
    expected_substitution_branch_length: float
    branch_rate_deviation: float
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class LocalClockRateRegimeRow:
    """One local-clock regime summary row."""

    regime_id: str
    target_kind: str
    target_label: str | None
    descendant_taxa: list[str]
    branch_count: int
    class_clock_rate: float
    log_prior_contribution: float


@dataclass(frozen=True, slots=True)
class LocalClockRateEvaluationReport:
    """One local-clock prior evaluation report."""

    family: str
    background_clock_rate: float
    log_standard_deviation: float
    dated_tree_newick: str
    substitution_tree_newick: str
    taxa: list[str]
    tip_count: int
    internal_node_count: int
    regime_count: int
    branch_count: int
    total_log_prior: float
    branch_rows: list[LocalClockRateBranchRow]
    regime_rows: list[LocalClockRateRegimeRow]


@dataclass(frozen=True, slots=True)
class _MatchedClockBranch:
    """One matched rooted branch pair between dated and substitution trees."""

    branch_id: str
    parent_branch_id: str | None
    child_name: str | None
    descendant_taxa: list[str]
    dated_time_duration: float
    observed_substitution_branch_length: float


@dataclass(frozen=True, slots=True)
class _ResolvedLocalClockSelector:
    """One tree-resolved local-clock selector with concrete branch ownership."""

    regime_id: str
    target_kind: str
    target_label: str | None
    descendant_taxa: tuple[str, ...]
    matching_branch_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _LocalClockBranchAssignment:
    """One resolved local-clock class assignment for one branch."""

    regime_id: str
    target_kind: str


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


def build_local_clock_rate_model(
    *,
    background_clock_rate: float,
    regime_clock_rates: Mapping[str, float],
    log_standard_deviation: float,
) -> LocalClockRateModel:
    """Build one local-clock prior model with one background class and named regimes."""
    validated_background_clock_rate = _validate_positive_clock_parameter(
        parameter_name="background_clock_rate",
        value=background_clock_rate,
        workflow_name="local clock prior",
        code="local_clock_rate_model_invalid_background_clock_rate",
    )
    if not math.isfinite(log_standard_deviation) or log_standard_deviation <= 0.0:
        raise PhylogeneticsError(
            "local clock prior requires a strictly positive finite log standard deviation",
            code="local_clock_rate_model_invalid_log_standard_deviation",
            details={"log_standard_deviation": log_standard_deviation},
        )
    validated_regime_clock_rates: dict[str, float] = {}
    for regime_id, regime_clock_rate in regime_clock_rates.items():
        normalized_regime_id = regime_id.strip()
        if not normalized_regime_id:
            raise PhylogeneticsError(
                "local clock prior requires non-empty regime ids",
                code="local_clock_rate_model_invalid_regime_id",
            )
        if normalized_regime_id.casefold() == _LOCAL_CLOCK_BACKGROUND_REGIME_ID:
            raise PhylogeneticsError(
                "local clock prior reserves regime_id 'background' for the background class",
                code="local_clock_rate_model_reserved_background_regime_id",
            )
        validated_regime_clock_rates[normalized_regime_id] = (
            _validate_positive_clock_parameter(
                parameter_name=f"regime_clock_rates[{normalized_regime_id!r}]",
                value=regime_clock_rate,
                workflow_name="local clock prior",
                code="local_clock_rate_model_invalid_regime_clock_rate",
            )
        )
    return LocalClockRateModel(
        family="local-clock",
        background_clock_rate=validated_background_clock_rate,
        regime_clock_rates=validated_regime_clock_rates,
        log_standard_deviation=float(format(log_standard_deviation, ".15g")),
    )


def load_local_clock_regime_definitions(
    tree: PhyloTree,
    regime_path: Path,
) -> list[LocalClockRegimeDefinition]:
    """Load one governed local-clock regime table against one rooted tree."""
    _validate_rooted_tree(
        tree,
        workflow_name="local clock prior",
        tree_role="selector",
    )
    with regime_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        required_columns = {"regime_id", "target_kind", "target_label", "taxa"}
        missing_columns = sorted(required_columns - set(fieldnames))
        if missing_columns:
            raise PhylogeneticsError(
                "local clock regime table is missing required columns: "
                + ", ".join(missing_columns),
                code="local_clock_regime_table_missing_columns",
            )
        rows = list(reader)

    if not rows:
        raise PhylogeneticsError(
            "local clock regime table must contain at least one regime row",
            code="local_clock_regime_table_empty",
        )

    known_descendant_taxa = {
        frozenset(node.descendant_taxa) for node in tree.iter_nodes(order="preorder")
    }
    regime_definitions: list[LocalClockRegimeDefinition] = []
    seen_regime_ids: set[str] = set()
    for row in rows:
        regime_id = (row.get("regime_id") or "").strip()
        if not regime_id:
            raise PhylogeneticsError(
                "local clock regime rows require a non-empty regime_id",
                code="local_clock_regime_table_invalid_regime_id",
            )
        if regime_id.casefold() == _LOCAL_CLOCK_BACKGROUND_REGIME_ID:
            raise PhylogeneticsError(
                "local clock regime_id 'background' is reserved for unassigned branches",
                code="local_clock_regime_table_reserved_background_regime_id",
            )
        if regime_id in seen_regime_ids:
            raise PhylogeneticsError(
                f"local clock regime_id '{regime_id}' is duplicated",
                code="local_clock_regime_table_duplicate_regime_id",
            )
        seen_regime_ids.add(regime_id)
        target_kind = (row.get("target_kind") or "").strip().casefold()
        if target_kind not in {"branch", "clade"}:
            raise PhylogeneticsError(
                "local clock target_kind must be either 'branch' or 'clade'",
                code="local_clock_regime_table_invalid_target_kind",
            )
        descendant_taxa = _parse_local_clock_selector_taxa(
            row.get("taxa"),
            regime_id=regime_id,
        )
        if frozenset(descendant_taxa) not in known_descendant_taxa:
            raise PhylogeneticsError(
                f"local clock regime '{regime_id}' does not match one stable tree clade",
                code="local_clock_regime_table_unknown_target",
            )
        regime_definitions.append(
            LocalClockRegimeDefinition(
                regime_id=regime_id,
                target_kind=target_kind,
                target_label=(row.get("target_label") or "").strip() or None,
                descendant_taxa=descendant_taxa,
            )
        )
    return regime_definitions


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
                branch_rate_deviation=float(format(branch_rate - anchor_rate, ".15g")),
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


def evaluate_local_clock_tree_log_prior(
    substitution_tree: PhyloTree,
    dated_tree: PhyloTree,
    regime_definitions: list[LocalClockRegimeDefinition],
    rate_model: LocalClockRateModel,
) -> LocalClockRateEvaluationReport:
    """Evaluate one rooted tree pair under one user-defined local-clock prior."""
    if rate_model.family != "local-clock":
        raise PhylogeneticsError(
            "local clock prior family is unsupported",
            code="local_clock_rate_model_family_invalid",
            details={"family": rate_model.family},
        )
    matched_branches = _collect_matched_clock_branches(
        substitution_tree,
        dated_tree,
        workflow_name="local clock prior",
        require_positive_dated_durations=True,
    )
    resolved_selectors = _resolve_local_clock_selectors(
        substitution_tree,
        regime_definitions,
    )
    branch_assignment_by_branch_id, branch_count_by_regime_id = (
        _resolve_local_clock_branch_assignments(
            substitution_tree,
            resolved_selectors,
        )
    )
    expected_regime_ids = {selector.regime_id for selector in resolved_selectors}
    provided_regime_ids = set(rate_model.regime_clock_rates)
    if expected_regime_ids != provided_regime_ids:
        raise PhylogeneticsError(
            "local clock prior requires one clock rate for every named regime",
            code="local_clock_rate_model_missing_regime_rates",
            details={
                "expected_regime_ids": sorted(expected_regime_ids),
                "provided_regime_ids": sorted(provided_regime_ids),
            },
        )

    class_clock_rate_by_regime_id = {
        _LOCAL_CLOCK_BACKGROUND_REGIME_ID: rate_model.background_clock_rate,
        **rate_model.regime_clock_rates,
    }
    branch_rows: list[LocalClockRateBranchRow] = []
    log_prior_by_regime_id: dict[str, float] = dict.fromkeys(
        class_clock_rate_by_regime_id,
        0.0,
    )
    for matched_branch in matched_branches:
        assignment = branch_assignment_by_branch_id[matched_branch.branch_id]
        class_clock_rate = class_clock_rate_by_regime_id[assignment.regime_id]
        branch_rate = (
            matched_branch.observed_substitution_branch_length
            / matched_branch.dated_time_duration
        )
        log_prior_contribution = _relaxed_clock_lognormal_log_density(
            branch_rate=branch_rate,
            anchor_rate=class_clock_rate,
            log_standard_deviation=rate_model.log_standard_deviation,
        )
        log_prior_by_regime_id[assignment.regime_id] += log_prior_contribution
        branch_rows.append(
            LocalClockRateBranchRow(
                branch_id=matched_branch.branch_id,
                child_name=matched_branch.child_name,
                descendant_taxa=matched_branch.descendant_taxa,
                regime_id=assignment.regime_id,
                target_kind=assignment.target_kind,
                dated_time_duration=float(
                    format(matched_branch.dated_time_duration, ".15g")
                ),
                observed_substitution_branch_length=float(
                    format(
                        matched_branch.observed_substitution_branch_length,
                        ".15g",
                    )
                ),
                class_clock_rate=float(format(class_clock_rate, ".15g")),
                branch_rate=float(format(branch_rate, ".15g")),
                expected_substitution_branch_length=float(
                    format(
                        class_clock_rate * matched_branch.dated_time_duration, ".15g"
                    )
                ),
                branch_rate_deviation=float(
                    format(branch_rate - class_clock_rate, ".15g")
                ),
                log_prior_contribution=log_prior_contribution,
            )
        )

    regime_rows = _build_local_clock_regime_rows(
        resolved_selectors=resolved_selectors,
        branch_count_by_regime_id=branch_count_by_regime_id,
        class_clock_rate_by_regime_id=class_clock_rate_by_regime_id,
        log_prior_by_regime_id=log_prior_by_regime_id,
    )
    return LocalClockRateEvaluationReport(
        family=rate_model.family,
        background_clock_rate=rate_model.background_clock_rate,
        log_standard_deviation=rate_model.log_standard_deviation,
        dated_tree_newick=dumps_newick(dated_tree),
        substitution_tree_newick=dumps_newick(substitution_tree),
        taxa=sorted(substitution_tree.tip_names),
        tip_count=substitution_tree.tip_count,
        internal_node_count=substitution_tree.internal_node_count,
        regime_count=len(regime_rows),
        branch_count=len(branch_rows),
        total_log_prior=sum(row.log_prior_contribution for row in branch_rows),
        branch_rows=branch_rows,
        regime_rows=regime_rows,
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
    return {tuple(child.descendant_taxa): child for _parent, child in tree.iter_edges()}


def _validate_positive_clock_parameter(
    *,
    parameter_name: str,
    value: float,
    workflow_name: str,
    code: str,
) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise PhylogeneticsError(
            f"{workflow_name} requires a strictly positive finite {parameter_name}",
            code=code,
            details={parameter_name: value},
        )
    return float(format(value, ".15g"))


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


def _parse_local_clock_selector_taxa(
    raw_taxa: str | None,
    *,
    regime_id: str,
) -> list[str]:
    descendant_taxa = [
        token.strip() for token in (raw_taxa or "").split("|") if token.strip()
    ]
    if not descendant_taxa:
        raise PhylogeneticsError(
            f"local clock regime '{regime_id}' requires at least one taxon",
            code="local_clock_regime_table_empty_taxa",
        )
    if len(set(descendant_taxa)) != len(descendant_taxa):
        raise PhylogeneticsError(
            f"local clock regime '{regime_id}' contains duplicate taxa",
            code="local_clock_regime_table_duplicate_taxa",
        )
    return sorted(descendant_taxa)


def _resolve_local_clock_selectors(
    tree: PhyloTree,
    regime_definitions: list[LocalClockRegimeDefinition],
) -> list[_ResolvedLocalClockSelector]:
    node_by_descendant_taxa = {
        frozenset(node.descendant_taxa): node
        for node in tree.iter_nodes(order="preorder")
    }
    resolved_selectors: list[_ResolvedLocalClockSelector] = []
    for regime_definition in regime_definitions:
        node = node_by_descendant_taxa.get(frozenset(regime_definition.descendant_taxa))
        if node is None or node.node_id is None:
            raise PhylogeneticsError(
                f"local clock regime '{regime_definition.regime_id}' does not match one stable tree clade",
                code="local_clock_regime_table_unknown_target",
            )
        if regime_definition.target_kind == "branch":
            if node is tree.root:
                raise PhylogeneticsError(
                    f"local clock branch regime '{regime_definition.regime_id}' cannot target the root",
                    code="local_clock_regime_table_root_branch_target",
                )
            matching_branch_ids = (node.node_id,)
        else:
            if node is tree.root:
                raise PhylogeneticsError(
                    f"local clock clade regime '{regime_definition.regime_id}' cannot target the whole tree; use the background rate instead",
                    code="local_clock_regime_table_root_clade_target",
                )
            target_taxa = frozenset(regime_definition.descendant_taxa)
            matching_branch_ids = tuple(
                child.node_id
                for _parent, child in tree.iter_edges()
                if child.node_id is not None
                and frozenset(child.descendant_taxa).issubset(target_taxa)
            )
        resolved_selectors.append(
            _ResolvedLocalClockSelector(
                regime_id=regime_definition.regime_id,
                target_kind=regime_definition.target_kind,
                target_label=regime_definition.target_label,
                descendant_taxa=tuple(regime_definition.descendant_taxa),
                matching_branch_ids=matching_branch_ids,
            )
        )
    return resolved_selectors


def _resolve_local_clock_branch_assignments(
    tree: PhyloTree,
    resolved_selectors: list[_ResolvedLocalClockSelector],
) -> tuple[dict[str, _LocalClockBranchAssignment], dict[str, int]]:
    branch_assignment_by_branch_id: dict[str, _LocalClockBranchAssignment] = {}
    branch_count_by_regime_id: dict[str, int] = {_LOCAL_CLOCK_BACKGROUND_REGIME_ID: 0}
    for selector in resolved_selectors:
        branch_count_by_regime_id[selector.regime_id] = 0
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise PhylogeneticsError(
                "local clock prior requires stable branch ids",
                code="local_clock_rate_model_requires_stable_branch_ids",
            )
        candidates = [
            selector
            for selector in resolved_selectors
            if child.node_id in selector.matching_branch_ids
        ]
        if not candidates:
            assignment = _LocalClockBranchAssignment(
                regime_id=_LOCAL_CLOCK_BACKGROUND_REGIME_ID,
                target_kind="background",
            )
        else:
            ranked_candidates = sorted(
                candidates,
                key=lambda selector: (
                    0 if selector.target_kind == "branch" else 1,
                    len(selector.descendant_taxa),
                    selector.regime_id,
                ),
            )
            top_selector = ranked_candidates[0]
            top_priority = (
                0 if top_selector.target_kind == "branch" else 1,
                len(top_selector.descendant_taxa),
            )
            tied_regime_ids = [
                selector.regime_id
                for selector in ranked_candidates
                if (
                    0 if selector.target_kind == "branch" else 1,
                    len(selector.descendant_taxa),
                )
                == top_priority
            ]
            if len(tied_regime_ids) > 1:
                raise PhylogeneticsError(
                    "local clock branch assignment is ambiguous for clade "
                    + "|".join(child.descendant_taxa)
                    + ": "
                    + ", ".join(tied_regime_ids),
                    code="local_clock_rate_model_ambiguous_branch_assignment",
                )
            assignment = _LocalClockBranchAssignment(
                regime_id=top_selector.regime_id,
                target_kind=top_selector.target_kind,
            )
        branch_assignment_by_branch_id[child.node_id] = assignment
        branch_count_by_regime_id[assignment.regime_id] = (
            branch_count_by_regime_id.get(assignment.regime_id, 0) + 1
        )
    empty_regime_ids = [
        selector.regime_id
        for selector in resolved_selectors
        if branch_count_by_regime_id[selector.regime_id] == 0
    ]
    if empty_regime_ids:
        raise PhylogeneticsError(
            "local clock selector resolution left one or more regimes without branches: "
            + ", ".join(sorted(empty_regime_ids)),
            code="local_clock_rate_model_empty_regime_assignment",
        )
    return branch_assignment_by_branch_id, branch_count_by_regime_id


def _build_local_clock_regime_rows(
    *,
    resolved_selectors: list[_ResolvedLocalClockSelector],
    branch_count_by_regime_id: Mapping[str, int],
    class_clock_rate_by_regime_id: Mapping[str, float],
    log_prior_by_regime_id: Mapping[str, float],
) -> list[LocalClockRateRegimeRow]:
    regime_rows = [
        LocalClockRateRegimeRow(
            regime_id=_LOCAL_CLOCK_BACKGROUND_REGIME_ID,
            target_kind="background",
            target_label=None,
            descendant_taxa=[],
            branch_count=branch_count_by_regime_id[_LOCAL_CLOCK_BACKGROUND_REGIME_ID],
            class_clock_rate=class_clock_rate_by_regime_id[
                _LOCAL_CLOCK_BACKGROUND_REGIME_ID
            ],
            log_prior_contribution=log_prior_by_regime_id[
                _LOCAL_CLOCK_BACKGROUND_REGIME_ID
            ],
        )
    ]
    for selector in resolved_selectors:
        regime_rows.append(
            LocalClockRateRegimeRow(
                regime_id=selector.regime_id,
                target_kind=selector.target_kind,
                target_label=selector.target_label,
                descendant_taxa=list(selector.descendant_taxa),
                branch_count=branch_count_by_regime_id[selector.regime_id],
                class_clock_rate=class_clock_rate_by_regime_id[selector.regime_id],
                log_prior_contribution=log_prior_by_regime_id[selector.regime_id],
            )
        )
    return regime_rows


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
