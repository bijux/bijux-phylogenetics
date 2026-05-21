from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralDiscreteDataset,
    dump_pruned_tree,
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
)

from .likelihood.fit_workflow import reconstruct_likelihood_estimates
from .models import DiscreteAncestralEstimate, DiscreteAncestralReport
from .policy import rerooting_method_compatibility, resolve_discrete_model_name

_STATE_SELECTION_REL_TOL = 1e-9
_STATE_SELECTION_ABS_TOL = 1e-12


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
    resolved_model = resolve_discrete_model_name(model)
    _validate_discrete_reconstruction_request(
        model=resolved_model,
        state_ordering=state_ordering,
        root_prior_mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return reconstruct_discrete_ancestral_states_from_dataset(
        dataset,
        model=resolved_model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        root_prior_mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
        allowed_transition_pairs=allowed_transition_pairs,
    )


def reconstruct_discrete_ancestral_states_from_dataset(
    dataset: AncestralDiscreteDataset,
    *,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
) -> DiscreteAncestralReport:
    """Reconstruct discrete ancestral states from one native ancestral dataset."""
    resolved_model = resolve_discrete_model_name(model)
    _validate_discrete_reconstruction_request(
        model=resolved_model,
        state_ordering=state_ordering,
        root_prior_mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    if resolved_model != "fitch":
        fit_result = reconstruct_likelihood_estimates(
            dataset,
            model=resolved_model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
            root_prior_mode=root_prior_mode,
            fixed_root_state=fixed_root_state,
            allowed_transition_pairs=allowed_transition_pairs,
            stable_probability_mapping=_stable_probability_mapping,
            material_state_set=_material_state_set,
            select_most_likely_state=_select_most_likely_state,
            build_discrete_estimate=_build_discrete_estimate,
            detect_discrete_overparameterization=_detect_discrete_overparameterization,
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
            fit_result.optimizer_diagnostics is not None
            and not fit_result.optimizer_diagnostics.converged
        ):
            warnings.append(
                "the discrete likelihood optimizer did not converge and should be interpreted cautiously"
            )
        if fit_result.optimizer_diagnostics is not None and (
            fit_result.optimizer_diagnostics.hit_lower_parameter_bound
            or fit_result.optimizer_diagnostics.hit_upper_parameter_bound
        ):
            warnings.append(
                "one or more discrete rate parameters hit an optimizer bound and should be interpreted as weakly identified"
            )
        if (
            fit_result.baseline_comparison is not None
            and fit_result.baseline_comparison.preferred_model_by_aic == "equal-rates"
        ):
            warnings.append(
                "the equal-rates baseline remains preferred by AIC over the requested discrete likelihood model"
            )
        warnings.extend(fit_result.rerooting_method_compatibility.notes)
        return DiscreteAncestralReport(
            tree_path=dataset.tree_path,
            traits_path=dataset.traits_path,
            taxon_column=dataset.taxon_column,
            trait=dataset.trait,
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
            rerooting_method_compatibility=fit_result.rerooting_method_compatibility,
            log_likelihood=fit_result.log_likelihood,
            parameter_count=fit_result.parameter_count,
            aic=fit_result.aic,
            transition_rate_rows=fit_result.transition_rate_rows,
            allowed_transition_pairs=fit_result.allowed_transition_pairs,
            optimizer_diagnostics=fit_result.optimizer_diagnostics,
            overparameterized=fit_result.overparameterized,
            baseline_comparison=fit_result.baseline_comparison,
        )
    return _reconstruct_fitch_report(
        dataset,
        model=resolved_model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )


def _validate_discrete_reconstruction_request(
    *,
    model: str,
    state_ordering: str,
    root_prior_mode: str,
    fixed_root_state: str | None,
    allowed_transition_pairs: list[tuple[str, str]] | None,
) -> None:
    if model == "fitch" and state_ordering != "unordered":
        raise ValueError(
            "ordered discrete ancestral reconstruction requires a likelihood model"
        )
    if model == "fitch" and (
        root_prior_mode != "equal" or fixed_root_state is not None
    ):
        raise ValueError(
            "fitch discrete ancestral reconstruction does not support root-prior assumptions"
        )
    if model == "fitch" and allowed_transition_pairs is not None:
        raise ValueError(
            "fitch discrete ancestral reconstruction does not support allowed-transition constraints"
        )


def _reconstruct_fitch_report(
    dataset: AncestralDiscreteDataset,
    *,
    model: str,
    state_ordering: str,
    ordered_states: list[str] | None,
) -> DiscreteAncestralReport:
    estimates: list[DiscreteAncestralEstimate] = []
    candidate_sets: dict[str, set[str]] = {}

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
    warnings.extend(
        rerooting_method_compatibility(
            model=model,
            state_ordering=state_ordering,
            root_prior_mode="equal",
        ).notes
    )

    return DiscreteAncestralReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=model,
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
        rerooting_method_compatibility=rerooting_method_compatibility(
            model=model,
            state_ordering=state_ordering,
            root_prior_mode="equal",
        ),
        log_likelihood=None,
        parameter_count=None,
        aic=None,
        transition_rate_rows=[],
        allowed_transition_pairs=[],
        optimizer_diagnostics=None,
        overparameterized=False,
        baseline_comparison=None,
    )


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
