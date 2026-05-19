from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    node_signature,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.core.traits import load_tsv_summary
from bijux_phylogenetics.diagnostics.validation import validate_tree_path
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from .models import (
    BiogeographicComputedResult as BiogeographicComputedResult,
    BiogeographicInterpretationReport as BiogeographicInterpretationReport,
    DiscreteEvolutionNarrative as DiscreteEvolutionNarrative,
    DiscreteEvolutionReportBuildResult as DiscreteEvolutionReportBuildResult,
    DiscreteModelComparisonReport as DiscreteModelComparisonReport,
    DiscreteModelComparisonRow as DiscreteModelComparisonRow,
    DiscreteStateEvolutionReport as DiscreteStateEvolutionReport,
    DiscreteTransitionReferenceObservation as DiscreteTransitionReferenceObservation,
    DiscreteTransitionReferenceRate as DiscreteTransitionReferenceRate,
    DiscreteTransitionReferenceValidationReport as DiscreteTransitionReferenceValidationReport,
    DominantStateBiasReport as DominantStateBiasReport,
    GeographicAnalysisReadinessReport as GeographicAnalysisReadinessReport,
    ModelSensitiveRegionRow as ModelSensitiveRegionRow,
    NodeStateDifference as NodeStateDifference,
    NodeStateEstimate as NodeStateEstimate,
    SparseStateInstabilityReport as SparseStateInstabilityReport,
    StateCodingAuditReport as StateCodingAuditReport,
    StateCodingAuditRow as StateCodingAuditRow,
    StateCodingIssue as StateCodingIssue,
    StateCodingValidationReport as StateCodingValidationReport,
    StateImbalanceReport as StateImbalanceReport,
    StateImbalanceWarning as StateImbalanceWarning,
    TransitionEvent as TransitionEvent,
    TransitionModelReport as TransitionModelReport,
    TransitionRateRow as TransitionRateRow,
    TransitionRateUncertaintyReport as TransitionRateUncertaintyReport,
    TransitionRateUncertaintyRow as TransitionRateUncertaintyRow,
    TransitionSummaryReport as TransitionSummaryReport,
    TransitionSupportRow as TransitionSupportRow,
)
from .transition_engine import (
    _build_transition_count_matrix,
    _estimate_node_states,
    _estimate_transition_rate_uncertainty,
    _estimate_transition_support_rows,
    _fitch_candidate_sets,
    _fit_transition_matrix,
    _normalize_probabilities,
    _pseudo_log_likelihood,
    _resolve_er_states,
    _resolve_state_order,
    _root_prior,
    _stationary_frequencies,
    _transition_events,
)

_DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST = ("|", "/", ";")
_DEFAULT_STATE_COLORS = (
    "#0f766e",
    "#1d4ed8",
    "#c2410c",
    "#7c3aed",
    "#b91c1c",
    "#047857",
    "#a16207",
    "#0f172a",
)


def _detect_sparse_state_instability(
    *,
    state_counts: dict[str, int],
    count_matrix: dict[str, dict[str, float]],
) -> SparseStateInstabilityReport:
    sparse_states = sorted(state for state, count in state_counts.items() if count < 2)
    zero_support_transitions = sorted(
        f"{source}->{target}"
        for source, targets in count_matrix.items()
        for target, value in targets.items()
        if value == 0.0
    )
    warning_count = int(bool(sparse_states)) + int(bool(zero_support_transitions))
    return SparseStateInstabilityReport(
        sparse_states=sparse_states,
        zero_support_transitions=zero_support_transitions,
        warning_count=warning_count,
        unstable=bool(sparse_states or zero_support_transitions),
    )


def _summarize_dominant_state_bias(
    state_counts: dict[str, int],
) -> DominantStateBiasReport:
    if not state_counts:
        return DominantStateBiasReport(
            dominant_states=[],
            dominant_fraction=0.0,
            biased=False,
            message=None,
        )
    total = sum(state_counts.values())
    max_count = max(state_counts.values())
    dominant_states = sorted(
        state for state, count in state_counts.items() if count == max_count
    )
    dominant_fraction = float(format(max_count / max(total, 1), ".15g"))
    biased = dominant_fraction >= 0.8
    return DominantStateBiasReport(
        dominant_states=dominant_states,
        dominant_fraction=dominant_fraction,
        biased=biased,
        message=(
            "one state dominates most taxa and may compress minority-state transition evidence"
            if biased
            else None
        ),
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


def _quantile(sorted_values: list[float], fraction: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return float(format(sorted_values[0], ".15g"))
    index = max(
        0, min(len(sorted_values) - 1, int(round(fraction * (len(sorted_values) - 1))))
    )
    return float(format(sorted_values[index], ".15g"))


def _build_model_sensitive_regions(
    differences: list[NodeStateDifference],
) -> list[ModelSensitiveRegionRow]:
    rows: list[ModelSensitiveRegionRow] = []
    for difference in differences:
        if not difference.differs:
            continue
        left_probability = difference.left_probabilities.get(difference.left_state, 0.0)
        right_probability = difference.right_probabilities.get(
            difference.right_state, 0.0
        )
        rows.append(
            ModelSensitiveRegionRow(
                node=difference.node,
                descendant_taxa=difference.descendant_taxa,
                left_state=difference.left_state,
                right_state=difference.right_state,
                sensitivity_score=float(
                    format(abs(left_probability - right_probability), ".15g")
                ),
            )
        )
    return sorted(rows, key=lambda row: (-row.sensitivity_score, row.node))


def validate_discrete_state_coding(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> StateCodingValidationReport:
    """Detect impossible or unsupported discrete-state labels."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    if state_ordering not in {"unordered", "ordered"}:
        raise ValueError(f"unsupported state ordering: {state_ordering}")
    ordered = list(ordered_states or [])
    if ordered and len(set(ordered)) != len(ordered):
        raise AncestralReconstructionError(
            "ordered state vocabulary contains duplicate labels"
        )
    allowed = list(allowed_states or ordered or [])
    allowed_set = set(allowed)
    issues: list[StateCodingIssue] = []
    usable_taxa: list[str] = []
    observed_states: set[str] = set()
    for row in table.rows:
        taxon = row[table.taxon_column]
        if taxon not in tree_taxa:
            continue
        raw_state = row[trait].strip()
        if not raw_state:
            continue
        invalid_delimiter = next(
            (
                token
                for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST
                if token in raw_state
            ),
            None,
        )
        if invalid_delimiter is not None:
            issues.append(
                StateCodingIssue(
                    taxon=taxon,
                    raw_state=raw_state,
                    code="unsupported-state-delimiter",
                    message=f"state label contains reserved delimiter '{invalid_delimiter}'",
                )
            )
            continue
        if allowed and raw_state not in allowed_set:
            issues.append(
                StateCodingIssue(
                    taxon=taxon,
                    raw_state=raw_state,
                    code="unordered-state-vocabulary"
                    if state_ordering == "ordered" and ordered
                    else "unsupported-state-label",
                    message=(
                        "state label is not present in the declared ordered vocabulary"
                        if state_ordering == "ordered" and ordered
                        else "state label is not present in the allowed state vocabulary"
                    ),
                )
            )
            continue
        observed_states.add(raw_state)
        usable_taxa.append(taxon)
    if state_ordering == "ordered" and ordered:
        missing_from_order = sorted(observed_states - set(ordered))
        issues.extend(
            StateCodingIssue(
                taxon="",
                raw_state=state,
                code="unordered-state-vocabulary",
                message="observed state is missing from the declared ordered vocabulary",
            )
            for state in missing_from_order
        )
    return StateCodingValidationReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        allowed_states=allowed,
        state_ordering=state_ordering,
        ordered_states=ordered,
        valid=not issues,
        issues=issues,
        observed_states=sorted(observed_states),
        usable_taxa=sorted(usable_taxa),
    )


def audit_discrete_state_coding(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    coding_map: dict[str, str] | None = None,
) -> StateCodingAuditReport:
    """Show how raw metadata states become model states or get excluded."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    mapping = dict(coding_map or {})
    ordered = list(ordered_states or [])
    allowed = list(allowed_states or ordered or [])
    allowed_set = set(allowed)
    rows: list[StateCodingAuditRow] = []
    for row in table.rows:
        taxon = row[table.taxon_column]
        raw_state = row[trait].strip()
        in_tree = taxon in tree_taxa
        if not raw_state:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=None,
                    in_tree=in_tree,
                    included=False,
                    issue_code="missing-state",
                    note="state is blank and cannot be used",
                )
            )
            continue
        normalized_state = mapping.get(raw_state, raw_state)
        invalid_delimiter = next(
            (
                token
                for token in _DEFAULT_ALLOWED_STATE_PATTERN_BLOCKLIST
                if token in normalized_state
            ),
            None,
        )
        if not in_tree:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=False,
                    included=False,
                    issue_code="taxon-not-in-tree",
                    note="taxon does not overlap the tree and is excluded before state modeling",
                )
            )
            continue
        if invalid_delimiter is not None:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=True,
                    included=False,
                    issue_code="unsupported-state-delimiter",
                    note=f"normalized state contains reserved delimiter '{invalid_delimiter}'",
                )
            )
            continue
        if allowed and normalized_state not in allowed_set:
            rows.append(
                StateCodingAuditRow(
                    taxon=taxon,
                    raw_state=raw_state,
                    normalized_state=normalized_state,
                    in_tree=True,
                    included=False,
                    issue_code="unsupported-state-label"
                    if state_ordering != "ordered"
                    else "unordered-state-vocabulary",
                    note=(
                        "normalized state is absent from the declared ordered vocabulary"
                        if state_ordering == "ordered" and ordered
                        else "normalized state is absent from the allowed state vocabulary"
                    ),
                )
            )
            continue
        rows.append(
            StateCodingAuditRow(
                taxon=taxon,
                raw_state=raw_state,
                normalized_state=normalized_state,
                in_tree=True,
                included=True,
                issue_code=None,
                note="state is retained for discrete-state modeling",
            )
        )
    included_row_count = sum(1 for row in rows if row.included)
    return StateCodingAuditReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        state_ordering=state_ordering,
        ordered_states=ordered,
        coding_map=mapping,
        row_count=len(rows),
        included_row_count=included_row_count,
        excluded_row_count=len(rows) - included_row_count,
        rows=rows,
    )


def detect_state_imbalance_problems(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> StateImbalanceReport:
    """Flag rare or degenerate discrete-state inputs over the tree overlap."""
    tree = load_tree(tree_path)
    table = (
        load_tsv_summary(traits_path)
        if taxon_column is None
        else load_taxon_table(traits_path, taxon_column=taxon_column)
    )
    if trait not in table.columns:
        raise AncestralReconstructionError(
            f"trait table does not contain column '{trait}'"
        )
    tree_taxa = set(tree.tip_names)
    state_counts: dict[str, int] = {}
    for row in table.rows:
        taxon = row[table.taxon_column]
        state = row[trait].strip()
        if taxon in tree_taxa and state:
            state_counts[state] = state_counts.get(state, 0) + 1
    observed_states = sorted(state_counts)
    warnings: list[StateImbalanceWarning] = []
    rare_states = sorted(state for state, count in state_counts.items() if count < 2)
    if len(observed_states) < 2:
        warnings.append(
            StateImbalanceWarning(
                code="single-state-dataset",
                message="only one observed state remains after pruning to usable tree taxa",
                affected_states=observed_states,
            )
        )
    if rare_states:
        warnings.append(
            StateImbalanceWarning(
                code="rare-states",
                message="one or more states are represented by fewer than two taxa",
                affected_states=rare_states,
            )
        )
    dominant_fraction = (
        max(state_counts.values()) / max(sum(state_counts.values()), 1)
        if state_counts
        else 0.0
    )
    if dominant_fraction >= 0.8 and observed_states:
        dominant_states = [
            state
            for state, count in state_counts.items()
            if count == max(state_counts.values())
        ]
        warnings.append(
            StateImbalanceWarning(
                code="dominant-state-skew",
                message="one state dominates most observed taxa and may overwhelm transition inference",
                affected_states=sorted(dominant_states),
            )
        )
    return StateImbalanceReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=table.taxon_column,
        trait=trait,
        taxon_count=sum(state_counts.values()),
        observed_states=observed_states,
        state_counts=state_counts,
        warnings=warnings,
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


def compare_discrete_state_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    left_model: str = "equal-rates",
    right_model: str = "all-rates-different",
    allowed_states: list[str] | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> DiscreteModelComparisonReport:
    """Compare discrete-state reconstructions across two supported models."""
    left = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=left_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right = run_discrete_state_transition_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=right_model,
        allowed_states=allowed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    right_by_node = {estimate.node: estimate for estimate in right.estimates}
    differences: list[NodeStateDifference] = []
    for left_estimate in left.estimates:
        right_estimate = right_by_node[left_estimate.node]
        differences.append(
            NodeStateDifference(
                node=left_estimate.node,
                descendant_taxa=left_estimate.descendant_taxa,
                left_state=left_estimate.most_likely_state,
                right_state=right_estimate.most_likely_state,
                differs=left_estimate.most_likely_state
                != right_estimate.most_likely_state,
                left_probabilities=left_estimate.state_probabilities,
                right_probabilities=right_estimate.state_probabilities,
            )
        )
    rows = [
        DiscreteModelComparisonRow(
            model=left_model,
            parameter_count=left.transition_model.parameter_count,
            pseudo_log_likelihood=left.transition_model.pseudo_log_likelihood,
            aic=left.transition_model.aic,
            transition_count=left.transition_summary.transition_count,
        ),
        DiscreteModelComparisonRow(
            model=right_model,
            parameter_count=right.transition_model.parameter_count,
            pseudo_log_likelihood=right.transition_model.pseudo_log_likelihood,
            aic=right.transition_model.aic,
            transition_count=right.transition_summary.transition_count,
        ),
    ]
    better_model = min(rows, key=lambda row: row.aic).model
    sensitive_regions = _build_model_sensitive_regions(differences)
    return DiscreteModelComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=left.taxon_column,
        trait=trait,
        left_model=left_model,
        right_model=right_model,
        better_model=better_model,
        rows=rows,
        node_differences=differences,
        sensitive_region_count=len(sensitive_regions),
        sensitive_regions=sensitive_regions,
    )
