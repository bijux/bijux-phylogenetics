from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import load_discrete_dataset, node_signature
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .models import (
    DiscreteStateEvolutionReport,
    GeographicAnalysisReadinessReport,
    TransitionModelReport,
    TransitionSummaryReport,
)
from .state_coding import (
    _detect_sparse_state_instability,
    _summarize_dominant_state_bias,
    detect_state_imbalance_problems,
    validate_discrete_state_coding,
)
from .transition_engine import (
    _build_transition_count_matrix,
    _estimate_node_states,
    _estimate_transition_rate_uncertainty,
    _estimate_transition_support_rows,
    _fit_transition_matrix,
    _fitch_candidate_sets,
    _normalize_probabilities,
    _pseudo_log_likelihood,
    _resolve_er_states,
    _resolve_state_order,
    _root_prior,
    _stationary_frequencies,
    _transition_events,
)


def assess_geographic_state_analysis_readiness(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> GeographicAnalysisReadinessReport:
    """Decide whether one geographic discrete-state analysis is credible enough to run."""
    tree_validation = validate_tree_path(tree_path)
    coding = validate_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    imbalance = detect_state_imbalance_problems(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    dominant_state_bias = _summarize_dominant_state_bias(imbalance.state_counts)
    blockers: list[str] = []
    warnings: list[str] = []

    if not tree_validation.syntax_valid or not tree_validation.biologically_safe:
        blockers.append(
            "tree validation failed and the geographic analysis is not safe to interpret"
        )
    if not tree_validation.rooted:
        blockers.append("geographic ancestral-state analysis requires a rooted tree")
    if not coding.valid:
        blockers.append(
            "discrete geographic states contain unsupported labels or coding patterns"
        )
    if any(warning.code == "single-state-dataset" for warning in imbalance.warnings):
        blockers.append(
            "geographic analysis requires at least two observed states after matching taxa to the tree"
        )
    rare_state_count = sum(1 for count in imbalance.state_counts.values() if count < 2)
    if imbalance.state_counts and rare_state_count == len(imbalance.state_counts):
        blockers.append(
            "one or more geographic states are too sparse to estimate transitions credibly"
        )
    if dominant_state_bias.biased:
        blockers.append(
            "observed geographic states are dominated by one state and the sampling is too biased for credible transition inference"
        )

    warnings.extend(tree_validation.warnings)
    warnings.extend(warning.message for warning in imbalance.warnings)
    if (
        dominant_state_bias.message is not None
        and dominant_state_bias.message not in warnings
    ):
        warnings.append(dominant_state_bias.message)

    return GeographicAnalysisReadinessReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=coding.taxon_column,
        trait=trait,
        valid=not blockers,
        blockers=blockers,
        warnings=warnings,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        coding_validation=coding,
        imbalance=imbalance,
        dominant_state_bias=dominant_state_bias,
        tree_validation_decision=tree_validation.validity_decision,
    )


def run_discrete_state_transition_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteStateEvolutionReport:
    """Run a deterministic discrete-state evolution workflow on one tree and trait."""
    if model == "meristic":
        _resolve_discrete_model_name(model)
    if model not in {"equal-rates", "symmetric", "all-rates-different"}:
        raise ValueError(f"unsupported discrete-state model: {model}")
    coding = validate_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    if not coding.valid:
        raise AncestralReconstructionError(
            "discrete-state evolution input contains unsupported state labels"
        )
    imbalance = detect_state_imbalance_problems(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if any(warning.code == "single-state-dataset" for warning in imbalance.warnings):
        raise AncestralReconstructionError(
            "discrete-state evolution requires at least two observed states"
        )
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    state_order = _resolve_state_order(
        dataset.observed_states,
        allowed_states=allowed_states,
        ordered_states=ordered_states,
        state_ordering=state_ordering,
    )
    candidate_sets = _fitch_candidate_sets(dataset.tree, dataset.states_by_taxon)
    stationary = _stationary_frequencies(dataset.states_by_taxon, state_order)
    er_resolved = _resolve_er_states(
        dataset.tree, candidate_sets, dataset.states_by_taxon, state_order
    )
    er_events = _transition_events(dataset.tree, er_resolved)
    count_matrix = _build_transition_count_matrix(
        state_order,
        er_events,
        model=model,
        state_ordering=state_ordering,
    )
    matrix = _fit_transition_matrix(
        model, state_order, stationary, er_events, state_ordering=state_ordering
    )
    root_prior = _root_prior(
        model, stationary, candidate_sets[node_signature(dataset.tree.root)]
    )
    estimates = _estimate_node_states(
        dataset.tree,
        candidate_sets,
        dataset.states_by_taxon,
        state_order,
        matrix,
        root_prior,
        state_ordering=state_ordering,
    )
    resolved_states = {
        estimate.node: estimate.most_likely_state for estimate in estimates
    }
    events = _transition_events(dataset.tree, resolved_states)
    transition_counts: dict[str, int] = {}
    for event in events:
        key = f"{event.source_state}->{event.target_state}"
        transition_counts[key] = transition_counts.get(key, 0) + 1
    support_rows = _estimate_transition_support_rows(
        estimates=estimates,
        events=events,
        transition_matrix=matrix,
    )
    branch_count = len(events)
    transition_count = sum(1 for event in events if event.changed)
    strongly_supported_transition_count = sum(
        1 for row in support_rows if row.strongly_supported
    )
    strongly_supported_transition_counts: dict[str, int] = {}
    for row in support_rows:
        if row.strongly_supported:
            strongly_supported_transition_counts[row.inferred_transition] = (
                strongly_supported_transition_counts.get(row.inferred_transition, 0) + 1
            )
    uncertainty = _estimate_transition_rate_uncertainty(
        model=model,
        state_ordering=state_ordering,
        transition_matrix=matrix,
        count_matrix=count_matrix,
    )
    instability = _detect_sparse_state_instability(
        state_counts=dataset.state_counts,
        count_matrix=count_matrix,
    )
    dominant_state_bias = _summarize_dominant_state_bias(dataset.state_counts)
    transition_model = TransitionModelReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        likelihood_method="deterministic-node-probability",
        state_ordering=state_ordering,
        ordered_states=state_order if state_ordering == "ordered" else [],
        state_order=state_order,
        parameter_count=(
            1
            if model == "equal-rates"
            else (
                len(state_order) * max(len(state_order) - 1, 0) // 2
                if model == "symmetric"
                else len(state_order) * max(len(state_order) - 1, 0)
            )
        ),
        pseudo_log_likelihood=0.0,
        aic=0.0,
        stationary_frequencies=stationary,
        transition_matrix=matrix,
        uncertainty=uncertainty,
        root_state_probabilities=_normalize_probabilities(root_prior),
    )
    transition_model.pseudo_log_likelihood = _pseudo_log_likelihood(
        estimates, events, transition_model
    )
    transition_model.aic = float(
        format(
            2.0 * transition_model.parameter_count
            - 2.0 * transition_model.pseudo_log_likelihood,
            ".15g",
        )
    )
    summary = TransitionSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        branch_count=branch_count,
        transition_count=transition_count,
        strongly_supported_transition_count=strongly_supported_transition_count,
        transition_counts=dict(sorted(transition_counts.items())),
        strongly_supported_transition_counts=dict(
            sorted(strongly_supported_transition_counts.items())
        ),
        support_rows=support_rows,
        events=events,
    )
    warnings = list(dataset.warnings)
    warnings.extend(warning.message for warning in imbalance.warnings)
    if instability.unstable:
        warnings.append(
            "sparse-state transition estimates may be unstable for one or more source-target paths"
        )
    if dominant_state_bias.biased and dominant_state_bias.message is not None:
        warnings.append(dominant_state_bias.message)
    return DiscreteStateEvolutionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        model=model,
        likelihood_method="deterministic-node-probability",
        state_ordering=state_ordering,
        ordered_states=state_order if state_ordering == "ordered" else [],
        analysis_tree_newick=dumps_newick(dataset.tree),
        taxon_count=len(dataset.taxa),
        observed_states=state_order,
        state_counts=dataset.state_counts,
        coding_validation=coding,
        imbalance=imbalance,
        instability=instability,
        dominant_state_bias=dominant_state_bias,
        transition_model=transition_model,
        estimates=estimates,
        transition_summary=summary,
        warnings=warnings,
    )


def estimate_ancestral_geographic_states(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteStateEvolutionReport:
    """Estimate ancestral geographic states over a rooted tree."""
    readiness = assess_geographic_state_analysis_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    if not readiness.valid:
        raise AncestralReconstructionError(
            "geographic state analysis is inappropriate: "
            + "; ".join(readiness.blockers)
        )
    return run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
