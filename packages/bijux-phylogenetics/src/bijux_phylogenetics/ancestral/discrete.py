from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.ancestral.common import (
    dump_pruned_tree,
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
    write_ancestral_rows,
)

_STATE_SELECTION_REL_TOL = 1e-9
_STATE_SELECTION_ABS_TOL = 1e-12
_DISCRETE_LOG_PARAMETER_LOWER_BOUND = -10.0
_DISCRETE_LOG_PARAMETER_UPPER_BOUND = 5.0
_DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE = 1e-10


@dataclass(slots=True)
class DiscreteAncestralEstimate:
    """One discrete ancestral-state estimate for a tree node."""

    node: str
    node_name: str | None
    is_tip: bool
    descendant_taxa: list[str]
    state_set: list[str]
    most_likely_state: str
    state_probabilities: dict[str, float]
    ambiguous: bool
    confidence: float
    interpretation: str
    unstable: bool
    downstream_risks: list[str]


@dataclass(slots=True)
class DiscreteTransitionRateRow:
    """One directed transition rate from a fitted discrete likelihood model."""

    source_state: str
    target_state: str
    transition_allowed: bool
    step_distance: int
    rate: float


@dataclass(slots=True)
class DiscreteOptimizerDiagnostics:
    """Optimizer state for one discrete likelihood fit."""

    optimizer_name: str
    parameter_count: int
    initial_candidate_count: int
    best_initial_scale: float
    converged: bool
    iteration_count: int
    function_evaluation_count: int
    simplex_shrink_count: int
    hit_lower_parameter_bound: bool
    hit_upper_parameter_bound: bool


@dataclass(slots=True)
class DiscreteModelBaselineComparison:
    """Likelihood-model comparison against the equal-rates baseline."""

    baseline_model: str
    baseline_log_likelihood: float
    baseline_parameter_count: int
    baseline_aic: float
    delta_log_likelihood: float
    delta_aic: float
    preferred_model_by_aic: str


@dataclass(slots=True)
class DiscreteLikelihoodFitResult:
    """Internal likelihood fit details for one discrete ancestral reconstruction."""

    estimates: list[DiscreteAncestralEstimate]
    ordered_states: list[str]
    state_order: list[str]
    log_likelihood: float
    parameter_count: int
    aic: float
    transition_rate_rows: list[DiscreteTransitionRateRow]
    allowed_transition_pairs: list[tuple[str, str]]
    optimizer_diagnostics: DiscreteOptimizerDiagnostics
    overparameterized: bool
    baseline_comparison: DiscreteModelBaselineComparison | None


@dataclass(slots=True)
class DiscreteAncestralReport:
    """Discrete ancestral-state reconstruction report."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    model: str
    state_ordering: str
    root_prior_mode: str | None
    fixed_root_state: str | None
    ordered_states: list[str]
    taxon_count: int
    observed_states: list[str]
    state_counts: dict[str, int]
    sparse_states: list[str]
    analysis_tree_newick: str
    dropped_missing_taxa: list[str]
    minimal_change_count: int | None
    parsimonious_root_state_count: int | None
    warnings: list[str]
    unstable_nodes: list[str]
    weak_support_nodes: list[str]
    estimates: list[DiscreteAncestralEstimate]
    log_likelihood: float | None
    parameter_count: int | None
    aic: float | None
    transition_rate_rows: list[DiscreteTransitionRateRow]
    allowed_transition_pairs: list[tuple[str, str]]
    optimizer_diagnostics: DiscreteOptimizerDiagnostics | None
    overparameterized: bool
    baseline_comparison: DiscreteModelBaselineComparison | None


@dataclass(slots=True)
class DiscreteAncestralSummary:
    """Reviewer-facing summary for one discrete ancestral reconstruction."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    root_prior_mode: str | None
    fixed_root_state: str | None
    analyzed_taxon_count: int
    excluded_taxon_count: int
    internal_node_count: int
    ambiguous_internal_node_count: int
    unstable_node_count: int
    weak_support_node_count: int
    observed_state_count: int
    sparse_state_count: int
    minimal_change_count: int | None
    parsimonious_root_state_count: int | None
    root_node: str
    root_most_likely_state: str
    root_confidence: float
    log_likelihood: float | None
    parameter_count: int | None
    aic: float | None
    optimizer_converged: bool | None
    optimizer_iteration_count: int | None
    optimizer_function_evaluation_count: int | None
    overparameterized: bool
    baseline_model: str | None
    baseline_delta_aic: float | None
    preferred_model_by_aic: str | None
    warning_count: int


@dataclass(slots=True)
class DiscreteAncestralExclusion:
    """One excluded tip from a discrete ancestral reconstruction."""

    taxon: str
    reason: str


def reconstruct_discrete_ancestral_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
) -> DiscreteAncestralReport:
    """Reconstruct discrete ancestral states under Fitch or Mk likelihood models."""
    resolved_model = _resolve_discrete_model_name(model)
    if resolved_model == "fitch" and state_ordering != "unordered":
        raise ValueError(
            "ordered discrete ancestral reconstruction requires a likelihood model"
        )
    if resolved_model == "fitch" and (
        root_prior_mode != "equal" or fixed_root_state is not None
    ):
        raise ValueError(
            "fitch discrete ancestral reconstruction does not support root-prior assumptions"
        )
    if resolved_model == "fitch" and allowed_transition_pairs is not None:
        raise ValueError(
            "fitch discrete ancestral reconstruction does not support allowed-transition constraints"
        )
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if resolved_model != "fitch":
        fit_result = _reconstruct_likelihood_estimates(
            dataset,
            model=resolved_model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
            root_prior_mode=root_prior_mode,
            fixed_root_state=fixed_root_state,
            allowed_transition_pairs=allowed_transition_pairs,
        )
        unstable_nodes = [
            estimate.node
            for estimate in fit_result.estimates
            if estimate.unstable and not estimate.is_tip
        ]
        weak_support_nodes = [
            estimate.node
            for estimate in fit_result.estimates
            if not estimate.is_tip and estimate.confidence < 0.75
        ]
        warnings = list(dataset.warnings)
        if unstable_nodes:
            warnings.append(
                "one or more discrete ancestral nodes remain unstable across candidate states"
            )
        if weak_support_nodes:
            warnings.append(
                "low-confidence ancestral state assignments should not be overinterpreted as definitive transitions"
            )
        if fit_result.overparameterized:
            warnings.append(
                "the discrete likelihood fit is likely overparameterized relative to the analyzed taxon count"
            )
        if (
            fit_result.baseline_comparison is not None
            and fit_result.baseline_comparison.preferred_model_by_aic
            == "equal-rates"
        ):
            warnings.append(
                "the equal-rates baseline remains preferred by AIC over the requested discrete likelihood model"
            )
        return DiscreteAncestralReport(
            tree_path=tree_path,
            traits_path=traits_path,
            taxon_column=dataset.taxon_column,
            trait=trait,
            model=resolved_model,
            state_ordering=state_ordering,
            root_prior_mode=root_prior_mode,
            fixed_root_state=fixed_root_state,
            ordered_states=fit_result.ordered_states,
            taxon_count=len(dataset.taxa),
            observed_states=dataset.observed_states,
            state_counts=dataset.state_counts,
            sparse_states=dataset.sparse_states,
            analysis_tree_newick=dump_pruned_tree(dataset.tree),
            dropped_missing_taxa=dataset.dropped_missing_taxa,
            minimal_change_count=None,
            parsimonious_root_state_count=None,
            warnings=warnings,
            unstable_nodes=unstable_nodes,
            weak_support_nodes=weak_support_nodes,
            estimates=fit_result.estimates,
            log_likelihood=fit_result.log_likelihood,
            parameter_count=fit_result.parameter_count,
            aic=fit_result.aic,
            transition_rate_rows=fit_result.transition_rate_rows,
            allowed_transition_pairs=fit_result.allowed_transition_pairs,
            optimizer_diagnostics=fit_result.optimizer_diagnostics,
            overparameterized=fit_result.overparameterized,
            baseline_comparison=fit_result.baseline_comparison,
        )
    estimates: list[DiscreteAncestralEstimate] = []

    candidate_sets: dict[str, set[str]] = {}
    minimal_change_count = 0

    def record_candidate_sets(node) -> tuple[set[str], int]:
        if node.is_leaf():
            state = dataset.states_by_taxon[node.name]
            candidate_set = {state}
            candidate_sets[node_signature(node)] = candidate_set
            return candidate_set, 0
        child_results = [record_candidate_sets(child) for child in node.children]
        candidate = set(child_results[0][0])
        minimal_changes = sum(result[1] for result in child_results)
        for child_set, _ in child_results[1:]:
            intersection = candidate & child_set
            if intersection:
                candidate = intersection
            else:
                candidate |= child_set
                minimal_changes += 1
        candidate_sets[node_signature(node)] = candidate
        return candidate, minimal_changes

    _, minimal_change_count = record_candidate_sets(dataset.tree.root)

    for node in dataset.tree.iter_nodes():
        signature = node_signature(node)
        if node.is_leaf():
            resolved_state = dataset.states_by_taxon[node.name]
            probabilities = {resolved_state: 1.0}
            state_set = [resolved_state]
        else:
            state_set = sorted(candidate_sets[signature])
            probabilities = {state: 1.0 / len(state_set) for state in state_set}
            resolved_state = state_set[0]
        estimates.append(
            _build_discrete_estimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                state_set=state_set,
                most_likely_state=resolved_state,
                state_probabilities=probabilities,
            )
        )
    unstable_nodes = [
        estimate.node
        for estimate in estimates
        if estimate.unstable and not estimate.is_tip
    ]
    weak_support_nodes = [
        estimate.node
        for estimate in estimates
        if not estimate.is_tip and estimate.confidence < 0.75
    ]
    warnings = list(dataset.warnings)
    if unstable_nodes:
        warnings.append(
            "one or more discrete ancestral nodes remain unstable across candidate states"
        )
    if weak_support_nodes:
        warnings.append(
            "low-confidence ancestral state assignments should not be overinterpreted as definitive transitions"
        )

    return DiscreteAncestralReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=resolved_model,
        state_ordering=state_ordering,
        root_prior_mode=None,
        fixed_root_state=None,
        ordered_states=list(ordered_states or []),
        taxon_count=len(dataset.taxa),
        observed_states=dataset.observed_states,
        state_counts=dataset.state_counts,
        sparse_states=dataset.sparse_states,
        analysis_tree_newick=dump_pruned_tree(dataset.tree),
        dropped_missing_taxa=dataset.dropped_missing_taxa,
        minimal_change_count=minimal_change_count,
        parsimonious_root_state_count=len(
            candidate_sets[node_signature(dataset.tree.root)]
        ),
        warnings=warnings,
        unstable_nodes=unstable_nodes,
        weak_support_nodes=weak_support_nodes,
        estimates=estimates,
        log_likelihood=None,
        parameter_count=None,
        aic=None,
        transition_rate_rows=[],
        allowed_transition_pairs=[],
        optimizer_diagnostics=None,
        overparameterized=False,
        baseline_comparison=None,
    )


def summarize_discrete_ancestral_report(
    report: DiscreteAncestralReport,
) -> DiscreteAncestralSummary:
    """Summarize the main review facts for one discrete ancestral report."""
    internal_estimates = [
        estimate for estimate in report.estimates if not estimate.is_tip
    ]
    if not internal_estimates:
        raise ValueError(
            "discrete ancestral summary requires at least one internal-node estimate"
        )
    ambiguous_internal_node_count = sum(
        1 for estimate in internal_estimates if estimate.ambiguous
    )
    root_estimate = max(
        internal_estimates,
        key=lambda estimate: (
            len(estimate.descendant_taxa),
            estimate.node,
        ),
    )
    return DiscreteAncestralSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        root_prior_mode=report.root_prior_mode,
        fixed_root_state=report.fixed_root_state,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.dropped_missing_taxa),
        internal_node_count=len(internal_estimates),
        ambiguous_internal_node_count=ambiguous_internal_node_count,
        unstable_node_count=len(report.unstable_nodes),
        weak_support_node_count=len(report.weak_support_nodes),
        observed_state_count=len(report.observed_states),
        sparse_state_count=len(report.sparse_states),
        minimal_change_count=report.minimal_change_count,
        parsimonious_root_state_count=report.parsimonious_root_state_count,
        root_node=root_estimate.node,
        root_most_likely_state=root_estimate.most_likely_state,
        root_confidence=root_estimate.confidence,
        log_likelihood=report.log_likelihood,
        parameter_count=report.parameter_count,
        aic=report.aic,
        optimizer_converged=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.converged
        ),
        optimizer_iteration_count=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.iteration_count
        ),
        optimizer_function_evaluation_count=(
            None
            if report.optimizer_diagnostics is None
            else report.optimizer_diagnostics.function_evaluation_count
        ),
        overparameterized=report.overparameterized,
        baseline_model=(
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.baseline_model
        ),
        baseline_delta_aic=(
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.delta_aic
        ),
        preferred_model_by_aic=(
            None
            if report.baseline_comparison is None
            else report.baseline_comparison.preferred_model_by_aic
        ),
        warning_count=len(report.warnings),
    )


def discrete_ancestral_exclusions(
    report: DiscreteAncestralReport,
) -> list[DiscreteAncestralExclusion]:
    """Return one explicit exclusion row per dropped tip taxon."""
    return [
        DiscreteAncestralExclusion(
            taxon=taxon,
            reason="missing_discrete_trait_state",
        )
        for taxon in report.dropped_missing_taxa
    ]


def write_discrete_ancestral_summary_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one summary ledger for a discrete ancestral reconstruction."""
    summary = summarize_discrete_ancestral_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "root_prior_mode",
            "fixed_root_state",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "unstable_node_count",
            "weak_support_node_count",
            "observed_state_count",
            "sparse_state_count",
            "minimal_change_count",
            "parsimonious_root_state_count",
            "root_node",
            "root_most_likely_state",
            "root_confidence",
            "log_likelihood",
            "parameter_count",
            "aic",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "overparameterized",
            "baseline_model",
            "baseline_delta_aic",
            "preferred_model_by_aic",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "root_prior_mode": summary.root_prior_mode or "",
                "fixed_root_state": summary.fixed_root_state or "",
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "unstable_node_count": str(summary.unstable_node_count),
                "weak_support_node_count": str(summary.weak_support_node_count),
                "observed_state_count": str(summary.observed_state_count),
                "sparse_state_count": str(summary.sparse_state_count),
                "minimal_change_count": _format_optional_int(
                    summary.minimal_change_count
                ),
                "parsimonious_root_state_count": _format_optional_int(
                    summary.parsimonious_root_state_count
                ),
                "root_node": summary.root_node,
                "root_most_likely_state": summary.root_most_likely_state,
                "root_confidence": str(summary.root_confidence),
                "log_likelihood": _format_optional_float(summary.log_likelihood),
                "parameter_count": _format_optional_int(summary.parameter_count),
                "aic": _format_optional_float(summary.aic),
                "optimizer_converged": _format_optional_bool(
                    summary.optimizer_converged
                ),
                "optimizer_iteration_count": _format_optional_int(
                    summary.optimizer_iteration_count
                ),
                "optimizer_function_evaluation_count": _format_optional_int(
                    summary.optimizer_function_evaluation_count
                ),
                "overparameterized": str(summary.overparameterized).lower(),
                "baseline_model": summary.baseline_model or "",
                "baseline_delta_aic": _format_optional_float(
                    summary.baseline_delta_aic
                ),
                "preferred_model_by_aic": summary.preferred_model_by_aic or "",
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_discrete_ancestral_probability_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one internal-node marginal-probability ledger for a discrete reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_state",
            "state_set",
            "state_probabilities",
            "confidence",
            "ambiguous",
            "unstable",
            "interpretation",
        ],
        rows=[
            {
                "node": estimate.node,
                "node_name": estimate.node_name or "",
                "descendant_taxa": ",".join(estimate.descendant_taxa),
                "most_likely_state": estimate.most_likely_state,
                "state_set": ",".join(estimate.state_set),
                "state_probabilities": json.dumps(
                    estimate.state_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(estimate.confidence),
                "ambiguous": str(estimate.ambiguous).lower(),
                "unstable": str(estimate.unstable).lower(),
                "interpretation": estimate.interpretation,
            }
            for estimate in report.estimates
            if not estimate.is_tip
        ],
    )


def write_discrete_ancestral_exclusion_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one explicit excluded-tip ledger for a discrete reconstruction."""
    exclusions = discrete_ancestral_exclusions(report)
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in exclusions
        ],
    )


def write_discrete_ancestral_transition_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one fitted transition-rate ledger for a discrete likelihood reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_state",
            "target_state",
            "transition_allowed",
            "step_distance",
            "rate",
        ],
        rows=[
            {
                "source_state": row.source_state,
                "target_state": row.target_state,
                "transition_allowed": str(row.transition_allowed).lower(),
                "step_distance": str(row.step_distance),
                "rate": str(row.rate),
            }
            for row in report.transition_rate_rows
        ],
    )


def write_discrete_ancestral_fit_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one optimizer and model-fit ledger for a discrete likelihood reconstruction."""
    diagnostics = report.optimizer_diagnostics
    baseline = report.baseline_comparison
    return write_ancestral_rows(
        path,
        columns=[
            "model",
            "taxon_count",
            "state_count",
            "parameter_count",
            "log_likelihood",
            "aic",
            "overparameterized",
            "optimizer_name",
            "optimizer_converged",
            "optimizer_iteration_count",
            "optimizer_function_evaluation_count",
            "simplex_shrink_count",
            "initial_candidate_count",
            "best_initial_scale",
            "hit_lower_parameter_bound",
            "hit_upper_parameter_bound",
            "baseline_model",
            "baseline_parameter_count",
            "baseline_log_likelihood",
            "baseline_aic",
            "delta_log_likelihood",
            "delta_aic",
            "preferred_model_by_aic",
        ],
        rows=[
            {
                "model": report.model,
                "taxon_count": str(report.taxon_count),
                "state_count": str(len(report.observed_states)),
                "parameter_count": _format_optional_int(report.parameter_count),
                "log_likelihood": _format_optional_float(report.log_likelihood),
                "aic": _format_optional_float(report.aic),
                "overparameterized": str(report.overparameterized).lower(),
                "optimizer_name": "" if diagnostics is None else diagnostics.optimizer_name,
                "optimizer_converged": _format_optional_bool(
                    None if diagnostics is None else diagnostics.converged
                ),
                "optimizer_iteration_count": _format_optional_int(
                    None if diagnostics is None else diagnostics.iteration_count
                ),
                "optimizer_function_evaluation_count": _format_optional_int(
                    None
                    if diagnostics is None
                    else diagnostics.function_evaluation_count
                ),
                "simplex_shrink_count": _format_optional_int(
                    None if diagnostics is None else diagnostics.simplex_shrink_count
                ),
                "initial_candidate_count": _format_optional_int(
                    None
                    if diagnostics is None
                    else diagnostics.initial_candidate_count
                ),
                "best_initial_scale": _format_optional_float(
                    None if diagnostics is None else diagnostics.best_initial_scale
                ),
                "hit_lower_parameter_bound": _format_optional_bool(
                    None
                    if diagnostics is None
                    else diagnostics.hit_lower_parameter_bound
                ),
                "hit_upper_parameter_bound": _format_optional_bool(
                    None
                    if diagnostics is None
                    else diagnostics.hit_upper_parameter_bound
                ),
                "baseline_model": "" if baseline is None else baseline.baseline_model,
                "baseline_parameter_count": _format_optional_int(
                    None if baseline is None else baseline.baseline_parameter_count
                ),
                "baseline_log_likelihood": _format_optional_float(
                    None if baseline is None else baseline.baseline_log_likelihood
                ),
                "baseline_aic": _format_optional_float(
                    None if baseline is None else baseline.baseline_aic
                ),
                "delta_log_likelihood": _format_optional_float(
                    None if baseline is None else baseline.delta_log_likelihood
                ),
                "delta_aic": _format_optional_float(
                    None if baseline is None else baseline.delta_aic
                ),
                "preferred_model_by_aic": (
                    "" if baseline is None else baseline.preferred_model_by_aic
                ),
            }
        ],
    )


def _resolve_discrete_model_name(model: str) -> str:
    aliases = {
        "fitch": "fitch",
        "equal-rates": "equal-rates",
        "er": "equal-rates",
        "symmetric": "symmetric",
        "sym": "symmetric",
        "all-rates-different": "all-rates-different",
        "ard": "all-rates-different",
    }
    resolved = aliases.get(model)
    if resolved is None:
        raise ValueError(f"unsupported discrete ancestral model: {model}")
    return resolved


def _reconstruct_likelihood_estimates(
    dataset,
    *,
    model: str,
    state_ordering: str,
    ordered_states: list[str] | None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    include_equal_rates_baseline: bool = True,
) -> DiscreteLikelihoodFitResult:
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    resolved_allowed_transition_pairs = _resolve_allowed_transition_pairs(
        state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    rate_matrix, default_root_prior, optimizer_diagnostics = _fit_discrete_mk_model(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    root_prior = _resolve_root_prior(
        state_order,
        state_counts=dataset.state_counts,
        mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
        default_root_prior=default_root_prior,
    )
    log_likelihood = _tree_log_likelihood(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    reported_log_likelihood = _reported_discrete_log_likelihood(
        log_likelihood,
        root_prior_mode=root_prior_mode,
        state_count=len(state_order),
    )
    posterior_by_node = _estimate_marginal_state_probabilities(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    estimates: list[DiscreteAncestralEstimate] = []
    for node in dataset.tree.iter_nodes():
        signature = node_signature(node)
        raw_probabilities = {
            state: float(probability)
            for state, probability in posterior_by_node[signature].items()
        }
        probabilities = _stable_probability_mapping(raw_probabilities)
        material_states = _material_state_set(
            probabilities,
            preferred_order=state_order,
        )
        estimates.append(
            _build_discrete_estimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=node_descendant_taxa(node),
                state_set=material_states,
                most_likely_state=_select_most_likely_state(
                    raw_probabilities,
                    preferred_order=state_order,
                ),
                state_probabilities=probabilities,
            )
        )
    parameter_count = _parameter_count(
        len(state_order),
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    aic = (2.0 * parameter_count) - (2.0 * reported_log_likelihood)
    baseline_comparison: DiscreteModelBaselineComparison | None = None
    if (
        include_equal_rates_baseline
        and model != "equal-rates"
        and allowed_transition_pairs is None
    ):
        try:
            baseline_fit = _reconstruct_likelihood_estimates(
                dataset,
                model="equal-rates",
                state_ordering=state_ordering,
                ordered_states=state_order,
                root_prior_mode=root_prior_mode,
                fixed_root_state=fixed_root_state,
                allowed_transition_pairs=None,
                include_equal_rates_baseline=False,
            )
        except (RuntimeError, ValueError):
            baseline_fit = None
        if baseline_fit is not None:
            baseline_comparison = DiscreteModelBaselineComparison(
                baseline_model="equal-rates",
                baseline_log_likelihood=baseline_fit.log_likelihood,
                baseline_parameter_count=baseline_fit.parameter_count,
                baseline_aic=baseline_fit.aic,
                delta_log_likelihood=(
                    reported_log_likelihood - baseline_fit.log_likelihood
                ),
                delta_aic=aic - baseline_fit.aic,
                preferred_model_by_aic=(
                    model if aic <= baseline_fit.aic else "equal-rates"
                ),
            )
    overparameterized = _detect_discrete_overparameterization(
        taxon_count=len(dataset.taxa),
        parameter_count=parameter_count,
    )
    return DiscreteLikelihoodFitResult(
        estimates=estimates,
        ordered_states=(state_order if state_ordering == "ordered" else []),
        state_order=state_order,
        log_likelihood=reported_log_likelihood,
        parameter_count=parameter_count,
        aic=aic,
        transition_rate_rows=_build_transition_rate_rows(
            state_order=state_order,
            state_ordering=state_ordering,
            rate_matrix=rate_matrix,
            allowed_transition_pairs=resolved_allowed_transition_pairs,
        ),
        allowed_transition_pairs=[
            (state_order[left_index], state_order[right_index])
            for left_index, right_index in sorted(resolved_allowed_transition_pairs)
        ],
        optimizer_diagnostics=optimizer_diagnostics,
        overparameterized=overparameterized,
        baseline_comparison=baseline_comparison,
    )


def _resolve_state_order(
    observed_states: list[str],
    *,
    state_ordering: str,
    ordered_states: list[str] | None,
) -> list[str]:
    if state_ordering == "unordered":
        return list(observed_states)
    if ordered_states is None:
        return list(observed_states)
    missing_states = sorted(set(observed_states) - set(ordered_states))
    if missing_states:
        raise ValueError(
            "ordered discrete ancestral reconstruction is missing observed states: "
            + ", ".join(missing_states)
        )
    return list(ordered_states)


def _fit_discrete_mk_model(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> tuple[numpy.ndarray, numpy.ndarray, DiscreteOptimizerDiagnostics]:
    parameter_count = _parameter_count(
        len(state_order),
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    if parameter_count < 1:
        raise ValueError(
            "discrete ancestral reconstruction requires at least one allowed transition"
        )
    initial_scales = (0.1, 1.0, 3.0)
    initial_candidates = [
        numpy.full(parameter_count, math.log(scale), dtype=float)
        for scale in initial_scales
    ]
    best_log_parameters: numpy.ndarray | None = None
    best_log_likelihood = float("-inf")
    best_run: _DiscreteOptimizationRun | None = None
    for initial_scale, initial in zip(initial_scales, initial_candidates, strict=True):
        run = _optimize_log_parameters(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            initial_log_parameters=initial,
            initial_scale=initial_scale,
        )
        if run.log_likelihood > best_log_likelihood:
            best_run = run
            best_log_parameters = run.log_parameters
            best_log_likelihood = run.log_likelihood
    if best_log_parameters is None or best_run is None:
        raise RuntimeError(
            "discrete ancestral optimization did not produce rate parameters"
        )
    if parameter_count > 1:
        # Plateau regularization is only intended to tame multi-parameter ridges
        # where several rate vectors induce effectively identical likelihood.
        # Applying it to identified one-parameter fits moves the estimate away
        # from the live `ape` maximum without improving stability.
        best_log_parameters = _regularize_plateau_log_parameters(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            log_parameters=best_log_parameters,
            reference_log_likelihood=best_log_likelihood,
        )
    rate_matrix = _rate_matrix_from_log_parameters(
        best_log_parameters,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    root_prior = _uniform_root_prior(len(state_order))
    return rate_matrix, root_prior, DiscreteOptimizerDiagnostics(
        optimizer_name="nelder-mead",
        parameter_count=parameter_count,
        initial_candidate_count=len(initial_candidates),
        best_initial_scale=best_run.initial_scale,
        converged=best_run.converged,
        iteration_count=best_run.iteration_count,
        function_evaluation_count=best_run.function_evaluation_count,
        simplex_shrink_count=best_run.simplex_shrink_count,
        hit_lower_parameter_bound=bool(
            numpy.any(
                best_log_parameters
                <= (_DISCRETE_LOG_PARAMETER_LOWER_BOUND + 1e-9)
            )
        ),
        hit_upper_parameter_bound=bool(
            numpy.any(
                best_log_parameters
                >= (_DISCRETE_LOG_PARAMETER_UPPER_BOUND - 1e-9)
            )
        ),
    )


@dataclass(slots=True)
class _DiscreteOptimizationRun:
    log_parameters: numpy.ndarray
    log_likelihood: float
    initial_scale: float
    converged: bool
    iteration_count: int
    function_evaluation_count: int
    simplex_shrink_count: int


def _optimize_log_parameters(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    initial_log_parameters: numpy.ndarray,
    initial_scale: float,
) -> _DiscreteOptimizationRun:
    simplex = [
        numpy.clip(
            initial_log_parameters.copy(),
            _DISCRETE_LOG_PARAMETER_LOWER_BOUND,
            _DISCRETE_LOG_PARAMETER_UPPER_BOUND,
        )
    ]
    for index in range(initial_log_parameters.size):
        vertex = simplex[0].copy()
        vertex[index] += 0.75
        simplex.append(
            numpy.clip(
                vertex,
                _DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                _DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
        )
    scores = [
        _evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            log_parameters=vertex,
        )
        for vertex in simplex
    ]
    function_evaluation_count = len(scores)
    alpha = 1.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5
    converged = False
    iteration_count = 0
    simplex_shrink_count = 0
    for iteration_count in range(1, 601):
        ordering = sorted(
            range(len(simplex)),
            key=lambda index: scores[index],
            reverse=True,
        )
        simplex = [simplex[index] for index in ordering]
        scores = [scores[index] for index in ordering]
        if (
            max(numpy.linalg.norm(vertex - simplex[0]) for vertex in simplex[1:]) < 1e-6
            and max(abs(score - scores[0]) for score in scores[1:]) < 1e-9
        ):
            converged = True
            break
        centroid = numpy.mean(simplex[:-1], axis=0)
        reflected = numpy.clip(
            centroid + alpha * (centroid - simplex[-1]),
            _DISCRETE_LOG_PARAMETER_LOWER_BOUND,
            _DISCRETE_LOG_PARAMETER_UPPER_BOUND,
        )
        reflected_score = _evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            log_parameters=reflected,
        )
        function_evaluation_count += 1
        if scores[0] >= reflected_score > scores[-2]:
            simplex[-1] = reflected
            scores[-1] = reflected_score
            continue
        if reflected_score > scores[0]:
            expanded = numpy.clip(
                centroid + gamma * (reflected - centroid),
                _DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                _DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
            expanded_score = _evaluate_log_likelihood(
                tree,
                states_by_taxon,
                state_order=state_order,
                model=model,
                state_ordering=state_ordering,
                allowed_transition_pairs=allowed_transition_pairs,
                log_parameters=expanded,
            )
            function_evaluation_count += 1
            if expanded_score > reflected_score:
                simplex[-1] = expanded
                scores[-1] = expanded_score
            else:
                simplex[-1] = reflected
                scores[-1] = reflected_score
            continue
        if reflected_score > scores[-1]:
            contracted = numpy.clip(
                centroid + rho * (reflected - centroid),
                _DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                _DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
        else:
            contracted = numpy.clip(
                centroid + rho * (simplex[-1] - centroid),
                _DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                _DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
        contracted_score = _evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            log_parameters=contracted,
        )
        function_evaluation_count += 1
        if contracted_score > scores[-1]:
            simplex[-1] = contracted
            scores[-1] = contracted_score
            continue
        best_vertex = simplex[0]
        new_simplex = [best_vertex]
        new_scores = [scores[0]]
        simplex_shrink_count += 1
        for vertex in simplex[1:]:
            shrunk = numpy.clip(
                best_vertex + sigma * (vertex - best_vertex),
                _DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                _DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
            new_simplex.append(shrunk)
            new_scores.append(
                _evaluate_log_likelihood(
                    tree,
                    states_by_taxon,
                    state_order=state_order,
                    model=model,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                    log_parameters=shrunk,
                )
            )
            function_evaluation_count += 1
        simplex = new_simplex
        scores = new_scores
    best_index = max(range(len(scores)), key=lambda index: scores[index])
    return _DiscreteOptimizationRun(
        log_parameters=simplex[best_index],
        log_likelihood=scores[best_index],
        initial_scale=initial_scale,
        converged=converged,
        iteration_count=iteration_count,
        function_evaluation_count=function_evaluation_count,
        simplex_shrink_count=simplex_shrink_count,
    )


def _evaluate_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    log_parameters: numpy.ndarray,
) -> float:
    rate_matrix = _rate_matrix_from_log_parameters(
        log_parameters,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    return _tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=_uniform_root_prior(len(state_order)),
    )


def _regularize_plateau_log_parameters(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    log_parameters: numpy.ndarray,
    reference_log_likelihood: float,
) -> numpy.ndarray:
    regularized = log_parameters.copy()
    for parameter_index in range(regularized.size):
        current_value = float(regularized[parameter_index])
        if current_value <= (_DISCRETE_LOG_PARAMETER_LOWER_BOUND + 1e-9):
            continue
        current_reference_log_likelihood = _evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            log_parameters=regularized,
        )
        if current_reference_log_likelihood < (
            reference_log_likelihood - _DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE
        ):
            current_reference_log_likelihood = reference_log_likelihood
        low = _DISCRETE_LOG_PARAMETER_LOWER_BOUND
        high = current_value
        for _ in range(80):
            midpoint = (low + high) / 2.0
            candidate = regularized.copy()
            candidate[parameter_index] = midpoint
            candidate_log_likelihood = _evaluate_log_likelihood(
                tree,
                states_by_taxon,
                state_order=state_order,
                model=model,
                state_ordering=state_ordering,
                allowed_transition_pairs=allowed_transition_pairs,
                log_parameters=candidate,
            )
            if (
                current_reference_log_likelihood - candidate_log_likelihood
                <= _DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE
            ):
                high = midpoint
            else:
                low = midpoint
        regularized[parameter_index] = high
    return regularized


def _parameter_count(
    state_count: int,
    *,
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> int:
    allowed_pair_count = len(allowed_transition_pairs)
    if model == "equal-rates":
        return 1 if allowed_pair_count > 0 else 0
    if state_ordering == "ordered":
        edge_count = max(
            sum(
                1
                for left_index in range(state_count)
                for right_index in range(left_index + 1, state_count)
                if _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
                and _transition_allowed(
                    right_index,
                    left_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
            ),
            0,
        )
        if model == "symmetric":
            return edge_count
        return allowed_pair_count
    if model == "symmetric":
        return max(
            sum(
                1
                for left_index in range(state_count)
                for right_index in range(left_index + 1, state_count)
                if _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
                and _transition_allowed(
                    right_index,
                    left_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                )
            ),
            0,
        )
    return allowed_pair_count


def _rate_matrix_from_log_parameters(
    log_parameters: numpy.ndarray,
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> numpy.ndarray:
    state_count = len(state_order)
    rate_matrix = numpy.zeros((state_count, state_count), dtype=float)
    parameter_index = 0
    if model == "equal-rates":
        rate = math.exp(float(log_parameters[0]))
        for left_index in range(state_count):
            for right_index in range(state_count):
                if _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    rate_matrix[left_index, right_index] = rate
    elif model == "symmetric":
        for left_index in range(state_count):
            for right_index in range(left_index + 1, state_count):
                if not _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                if not _transition_allowed(
                    right_index,
                    left_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                rate = math.exp(float(log_parameters[parameter_index]))
                parameter_index += 1
                rate_matrix[left_index, right_index] = rate
                rate_matrix[right_index, left_index] = rate
    else:
        for left_index in range(state_count):
            for right_index in range(state_count):
                if not _transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                rate_matrix[left_index, right_index] = math.exp(
                    float(log_parameters[parameter_index])
                )
                parameter_index += 1
    for state_index in range(state_count):
        rate_matrix[state_index, state_index] = -float(
            numpy.sum(rate_matrix[state_index, :])
        )
    return rate_matrix


def _transition_allowed(
    left_index: int,
    right_index: int,
    *,
    state_count: int,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> bool:
    if left_index == right_index:
        return False
    allowed_by_order = state_ordering == "unordered" or (
        abs(left_index - right_index) == 1
        and max(left_index, right_index) < state_count
    )
    if not allowed_by_order:
        return False
    return (left_index, right_index) in allowed_transition_pairs


def _resolve_allowed_transition_pairs(
    state_order: list[str],
    *,
    model: str,
    state_ordering: str,
    allowed_transition_pairs: list[tuple[str, str]] | None,
) -> set[tuple[int, int]]:
    state_to_index = {state: index for index, state in enumerate(state_order)}
    if allowed_transition_pairs is None:
        pairs = {
            (left_index, right_index)
            for left_index in range(len(state_order))
            for right_index in range(len(state_order))
            if left_index != right_index
        }
    else:
        pairs: set[tuple[int, int]] = set()
        for source_state, target_state in allowed_transition_pairs:
            if source_state not in state_to_index:
                raise ValueError(
                    "allowed transition source state is not present in the analyzed state vocabulary: "
                    f"{source_state}"
                )
            if target_state not in state_to_index:
                raise ValueError(
                    "allowed transition target state is not present in the analyzed state vocabulary: "
                    f"{target_state}"
                )
            if source_state == target_state:
                raise ValueError(
                    "allowed transition pairs must connect distinct states"
                )
            pairs.add((state_to_index[source_state], state_to_index[target_state]))
    if model == "symmetric":
        asymmetric_pairs = [
            (left_index, right_index)
            for left_index, right_index in sorted(pairs)
            if (right_index, left_index) not in pairs
        ]
        if asymmetric_pairs:
            raise ValueError(
                "symmetric discrete ancestral reconstruction requires bidirectional allowed transitions"
            )
    filtered_pairs = {
        (left_index, right_index)
        for left_index, right_index in pairs
        if left_index != right_index
        and (state_ordering == "unordered" or abs(left_index - right_index) == 1)
    }
    if not filtered_pairs:
        raise ValueError(
            "discrete ancestral reconstruction requires at least one allowed transition after applying constraints"
        )
    return filtered_pairs


def _build_transition_rate_rows(
    *,
    state_order: list[str],
    state_ordering: str,
    rate_matrix: numpy.ndarray,
    allowed_transition_pairs: set[tuple[int, int]],
) -> list[DiscreteTransitionRateRow]:
    rows: list[DiscreteTransitionRateRow] = []
    state_count = len(state_order)
    for left_index, source_state in enumerate(state_order):
        for right_index, target_state in enumerate(state_order):
            if left_index == right_index:
                continue
            rows.append(
                DiscreteTransitionRateRow(
                    source_state=source_state,
                    target_state=target_state,
                    transition_allowed=_transition_allowed(
                        left_index,
                        right_index,
                        state_count=state_count,
                        state_ordering=state_ordering,
                        allowed_transition_pairs=allowed_transition_pairs,
                    ),
                    step_distance=abs(left_index - right_index),
                    rate=float(rate_matrix[left_index, right_index]),
                )
            )
    return rows


def _tree_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> float:
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = _transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    def visit(node) -> tuple[numpy.ndarray, float]:
        if node.is_leaf():
            likelihood = numpy.zeros(len(state_order), dtype=float)
            likelihood[state_index[states_by_taxon[node.name]]] = 1.0
            return likelihood, 0.0
        partial = numpy.ones(len(state_order), dtype=float)
        log_scale = 0.0
        for child in node.children:
            child_partial, child_scale = visit(child)
            branch_transition = transition(_branch_length(child))
            partial *= branch_transition @ child_partial
            log_scale += child_scale
        scale = float(partial.sum())
        if scale <= 0.0:
            return partial, float("-inf")
        partial /= scale
        return partial, log_scale + math.log(scale)

    root_partial, subtree_log_scale = visit(tree.root)
    root_weight = root_prior * root_partial
    root_scale = float(root_weight.sum())
    if root_scale <= 0.0:
        return float("-inf")
    return subtree_log_scale + math.log(root_scale)


def _estimate_marginal_state_probabilities(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    state_index = {state: index for index, state in enumerate(state_order)}
    posterior_by_node: dict[str, numpy.ndarray] = {}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = _transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    def postorder(node) -> numpy.ndarray:
        signature = node_signature(node)
        if node.is_leaf():
            partial = numpy.zeros(len(state_order), dtype=float)
            partial[state_index[states_by_taxon[node.name]]] = 1.0
            posterior_by_node[signature] = partial
            return partial
        partial = numpy.ones(len(state_order), dtype=float)
        for child in node.children:
            child_partial = postorder(child)
            partial *= transition(_branch_length(child)) @ child_partial
        partial = _normalize_array(partial)
        posterior_by_node[signature] = partial
        return partial

    postorder(tree.root)
    root_signature = node_signature(tree.root)
    posterior_by_node[root_signature] = _normalize_array(
        root_prior * posterior_by_node[root_signature]
    )

    def preorder(node) -> None:
        parent_signature = node_signature(node)
        if node.is_leaf():
            return
        parent_probabilities = posterior_by_node[parent_signature]
        for child in node.children:
            if child.is_leaf():
                continue
            child_signature = node_signature(child)
            branch_transition = transition(_branch_length(child))
            child_probabilities = posterior_by_node[child_signature]
            denominator = child_probabilities @ branch_transition
            updated = (parent_probabilities / denominator) @ branch_transition
            posterior_by_node[child_signature] = _normalize_array(
                updated * child_probabilities
            )
            preorder(child)

    preorder(tree.root)
    return {
        node: {
            state: float(format(probability, ".15g"))
            for state, probability in zip(
                state_order,
                _normalize_array(probabilities),
                strict=True,
            )
        }
        for node, probabilities in posterior_by_node.items()
    }


def _transition_probability_matrix(
    rate_matrix: numpy.ndarray,
    branch_length: float,
) -> numpy.ndarray:
    if branch_length <= 0.0:
        return numpy.eye(rate_matrix.shape[0], dtype=float)
    eigenvalues, eigenvectors = numpy.linalg.eig(rate_matrix)
    inverse_vectors = numpy.linalg.inv(eigenvectors)
    diagonal = numpy.diag(numpy.exp(eigenvalues * branch_length))
    transition = eigenvectors @ diagonal @ inverse_vectors
    transition = numpy.real_if_close(transition, tol=1000).astype(float)
    transition[transition < 0.0] = 0.0
    row_sums = transition.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return transition / row_sums


def _normalize_array(values: numpy.ndarray) -> numpy.ndarray:
    total = float(values.sum())
    if total <= 0.0:
        return numpy.full(values.shape[0], 1.0 / values.shape[0], dtype=float)
    return values / total


def _uniform_root_prior(state_count: int) -> numpy.ndarray:
    return numpy.full(state_count, 1.0 / state_count, dtype=float)


def _empirical_root_prior(
    state_order: list[str], state_counts: dict[str, int]
) -> numpy.ndarray:
    return _normalize_array(
        numpy.array(
            [float(state_counts.get(state, 0)) for state in state_order],
            dtype=float,
        )
    )


def _fixed_root_prior(state_order: list[str], fixed_root_state: str) -> numpy.ndarray:
    if fixed_root_state not in state_order:
        raise ValueError(
            "fixed root state is not available in the analyzed state vocabulary: "
            f"{fixed_root_state}"
        )
    prior = numpy.zeros(len(state_order), dtype=float)
    prior[state_order.index(fixed_root_state)] = 1.0
    return prior


def _resolve_root_prior(
    state_order: list[str],
    *,
    state_counts: dict[str, int],
    mode: str,
    fixed_root_state: str | None,
    default_root_prior: numpy.ndarray | None = None,
) -> numpy.ndarray:
    if mode == "equal":
        if fixed_root_state is not None:
            raise ValueError("fixed_root_state requires root_prior_mode 'fixed'")
        if default_root_prior is not None:
            return default_root_prior
        return _uniform_root_prior(len(state_order))
    if mode == "empirical":
        if fixed_root_state is not None:
            raise ValueError("fixed_root_state requires root_prior_mode 'fixed'")
        return _empirical_root_prior(state_order, state_counts)
    if mode == "fixed":
        if fixed_root_state is None:
            raise ValueError("root_prior_mode 'fixed' requires a fixed_root_state")
        return _fixed_root_prior(state_order, fixed_root_state)
    raise ValueError(f"unsupported discrete ancestral root prior mode: {mode}")


def _reported_discrete_log_likelihood(
    log_likelihood: float,
    *,
    root_prior_mode: str,
    state_count: int,
) -> float:
    if root_prior_mode != "equal" or state_count < 2:
        return log_likelihood
    return log_likelihood + math.log(state_count)


def _branch_length(node) -> float:
    if node.branch_length is None:
        return 1.0
    return max(float(node.branch_length), 0.0)


def _material_state_set(
    state_probabilities: dict[str, float],
    *,
    preferred_order: list[str] | None = None,
) -> list[str]:
    return sorted(
        state
        for state, probability in state_probabilities.items()
        if probability >= 0.1
    ) or [
        _select_most_likely_state(
            state_probabilities,
            preferred_order=preferred_order,
        )
    ]


def _stable_probability_mapping(
    state_probabilities: dict[str, float],
) -> dict[str, float]:
    return {
        state: stable_value(probability)
        for state, probability in state_probabilities.items()
    }


def _select_most_likely_state(
    state_probabilities: dict[str, float],
    *,
    preferred_order: list[str] | None = None,
) -> str:
    if not state_probabilities:
        raise ValueError("state probabilities must not be empty")
    best_probability = max(state_probabilities.values())
    tied_states = [
        state
        for state, probability in state_probabilities.items()
        if math.isclose(
            probability,
            best_probability,
            rel_tol=_STATE_SELECTION_REL_TOL,
            abs_tol=_STATE_SELECTION_ABS_TOL,
        )
    ]
    if preferred_order is not None:
        order_lookup = {state: index for index, state in enumerate(preferred_order)}
        return max(
            tied_states,
            key=lambda state: (order_lookup.get(state, len(order_lookup)), state),
        )
    return max(tied_states)


def _build_discrete_estimate(
    *,
    node: str,
    node_name: str | None,
    is_tip: bool,
    descendant_taxa: list[str],
    state_set: list[str] | None = None,
    most_likely_state: str,
    state_probabilities: dict[str, float],
) -> DiscreteAncestralEstimate:
    resolved_state_set = sorted(state_set or state_probabilities)
    ordered_probabilities = sorted(state_probabilities.values(), reverse=True)
    confidence = ordered_probabilities[0] if ordered_probabilities else 0.0
    runner_up = ordered_probabilities[1] if len(ordered_probabilities) > 1 else 0.0
    unstable = not is_tip and ((confidence - runner_up) < 0.15 or confidence < 0.7)
    if is_tip:
        interpretation = "observed tip state"
    elif unstable:
        interpretation = "unstable node state"
    elif confidence >= 0.9:
        interpretation = "strongly supported node state"
    else:
        interpretation = "moderately supported node state"
    return DiscreteAncestralEstimate(
        node=node,
        node_name=node_name,
        is_tip=is_tip,
        descendant_taxa=descendant_taxa,
        state_set=resolved_state_set,
        most_likely_state=most_likely_state,
        state_probabilities=state_probabilities,
        ambiguous=len(resolved_state_set) > 1,
        confidence=confidence,
        interpretation=interpretation,
        unstable=unstable,
        downstream_risks=_discrete_downstream_risks(unstable),
    )


def _discrete_downstream_risks(unstable: bool) -> list[str]:
    if not unstable:
        return []
    return [
        "transition counts and inferred ancestral geography may change under alternative state models",
        "biological narratives about ancestral states should be treated as provisional for this node",
    ]


def _detect_discrete_overparameterization(
    *,
    taxon_count: int,
    parameter_count: int,
) -> bool:
    return parameter_count >= taxon_count


def _format_optional_int(value: int | None) -> str:
    if value is None:
        return ""
    return str(value)


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)


def _format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return str(value).lower()
