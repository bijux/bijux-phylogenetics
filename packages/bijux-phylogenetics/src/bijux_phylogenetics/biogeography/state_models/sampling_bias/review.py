from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    node_descendant_taxa,
    node_signature,
    stable_value,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    NodeStateEstimate,
    TransitionEvent,
    TransitionSupportRow,
    _estimate_node_states,
    _estimate_transition_support_rows,
    _fit_transition_matrix,
    _fitch_candidate_sets,
    _resolve_state_order,
    _root_prior,
    _stationary_frequencies,
    _transition_events,
    audit_discrete_state_coding,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .comparison_assembly import (
    build_exclusion_rows,
    build_node_rows,
    build_summary,
    build_transition_rows,
)
from .contracts import (
    GeographicSamplingBiasNodeRow,
    GeographicSamplingBiasReport,
    GeographicSamplingBiasSummary,
    GeographicSamplingBiasTransitionRow,
    GeographicSamplingCountRow,
)
from .weighting_policy import (
    build_count_rows,
    build_warnings,
    included_region_counts,
    resolve_region_weights,
)

_MODEL_ALIAS_TO_INTERNAL = {
    "er": "equal-rates",
    "sym": "symmetric",
    "ard": "all-rates-different",
    "equal-rates": "equal-rates",
    "symmetric": "symmetric",
    "all-rates-different": "all-rates-different",
}

_INTERNAL_MODEL_TO_ALIAS = {
    "equal-rates": "er",
    "symmetric": "sym",
    "all-rates-different": "ard",
}

@dataclass(frozen=True, slots=True)
class _SamplingBiasModelSurface:
    taxon_column: str
    estimates: list[NodeStateEstimate]
    events: list[TransitionEvent]
    support_rows: list[TransitionSupportRow]


def summarize_geographic_sampling_bias(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "er",
    allowed_regions: list[str] | None = None,
    weights_path: Path | None = None,
    region_column: str = "region",
    weight_column: str = "weight",
) -> GeographicSamplingBiasReport:
    """Review how explicit region weights change biogeographic state inference."""
    internal_model = _resolve_internal_model(model)
    audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        allowed_states=allowed_regions,
    )
    exclusion_rows = build_exclusion_rows(audit)
    included_counts = included_region_counts(audit)
    weights, weighting_mode = resolve_region_weights(
        included_counts,
        weights_path=weights_path,
        region_column=region_column,
        weight_column=weight_column,
    )
    baseline = _run_sampling_bias_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_regions=allowed_regions,
        region_weights=None,
    )
    weighted = _run_sampling_bias_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=internal_model,
        allowed_regions=allowed_regions,
        region_weights=weights,
    )
    count_rows = build_count_rows(included_counts, weights)
    node_rows = build_node_rows(baseline, weighted)
    transition_rows = build_transition_rows(baseline, weighted)
    warnings = build_warnings(count_rows, weighting_mode, node_rows)
    summary = build_summary(
        trait=trait,
        taxon_column=baseline.taxon_column,
        model=_INTERNAL_MODEL_TO_ALIAS[internal_model],
        internal_model=internal_model,
        analyzed_taxon_count=sum(included_counts.values()),
        excluded_taxon_count=len(exclusion_rows),
        weighting_mode=weighting_mode,
        count_rows=count_rows,
        node_rows=node_rows,
        transition_rows=transition_rows,
        warning_count=len(warnings),
    )
    return GeographicSamplingBiasReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=baseline.taxon_column,
        model=summary.model,
        internal_model=internal_model,
        weighting_mode=weighting_mode,
        summary=summary,
        count_rows=count_rows,
        node_rows=node_rows,
        transition_rows=transition_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )
def _normalize_probabilities(probabilities: dict[str, float]) -> dict[str, float]:
    total = sum(probabilities.values())
    if total <= 0.0:
        uniform = 1.0 / max(len(probabilities), 1)
        return {region: stable_value(uniform) for region in sorted(probabilities)}
    return {
        region: stable_value(probability / total)
        for region, probability in sorted(probabilities.items())
    }


def _most_likely_region(probabilities: dict[str, float]) -> str:
    return max(sorted(probabilities), key=lambda region: probabilities[region])


def _resolve_internal_model(model: str) -> str:
    if model not in _MODEL_ALIAS_TO_INTERNAL:
        raise ValueError(f"unsupported biogeography model alias: {model}")
    return _MODEL_ALIAS_TO_INTERNAL[model]


def _run_sampling_bias_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
    model: str,
    allowed_regions: list[str] | None,
    region_weights: dict[str, float] | None,
) -> _SamplingBiasModelSurface:
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    if len(dataset.observed_states) < 2:
        raise AncestralReconstructionError(
            "geographic sampling-bias review requires at least two observed regions"
        )
    state_order = _resolve_state_order(
        dataset.observed_states,
        allowed_states=allowed_regions,
        ordered_states=None,
        state_ordering="unordered",
    )
    candidate_sets = _fitch_candidate_sets(dataset.tree, dataset.states_by_taxon)
    stationary = (
        _weighted_stationary_frequencies(
            dataset.states_by_taxon, state_order, region_weights
        )
        if region_weights is not None
        else _stationary_frequencies(dataset.states_by_taxon, state_order)
    )
    priority_weights = (
        {state: region_weights.get(state, 1.0) for state in state_order}
        if region_weights is not None
        else dict.fromkeys(state_order, 1.0)
    )
    er_resolved = _resolve_biased_er_states(
        dataset.tree,
        candidate_sets,
        dataset.states_by_taxon,
        state_order,
        priority_weights,
    )
    er_events = _transition_events(dataset.tree, er_resolved)
    matrix = _fit_transition_matrix(
        model,
        state_order,
        stationary,
        er_events,
        state_ordering="unordered",
    )
    root_prior = _root_prior(
        model,
        stationary,
        candidate_sets[node_signature(dataset.tree.root)],
    )
    estimates = _estimate_node_states(
        dataset.tree,
        candidate_sets,
        dataset.states_by_taxon,
        state_order,
        matrix,
        root_prior,
        state_ordering="unordered",
    )
    resolved_states = {
        estimate.node: estimate.most_likely_state for estimate in estimates
    }
    events = _transition_events(dataset.tree, resolved_states)
    support_rows = _estimate_transition_support_rows(
        estimates=estimates,
        events=events,
        transition_matrix=matrix,
    )
    return _SamplingBiasModelSurface(
        taxon_column=dataset.taxon_column,
        estimates=estimates,
        events=events,
        support_rows=support_rows,
    )


def _weighted_stationary_frequencies(
    states_by_taxon: dict[str, str],
    state_order: list[str],
    region_weights: dict[str, float],
) -> dict[str, float]:
    weighted_counts = {
        state: stable_value(
            sum(
                region_weights.get(observed_state, 1.0)
                for observed_state in states_by_taxon.values()
                if observed_state == state
            )
        )
        for state in state_order
    }
    return _normalize_probabilities(weighted_counts)


def _resolve_biased_er_states(
    tree,
    candidate_sets: dict[str, list[str]],
    states_by_taxon: dict[str, str],
    state_order: list[str],
    priority_weights: dict[str, float],
) -> dict[str, str]:
    support_by_node = {
        node_signature(node): _state_support(node, states_by_taxon, state_order)
        for node in tree.iter_nodes()
    }
    resolved: dict[str, str] = {}
    root_signature = node_signature(tree.root)
    root_candidates = candidate_sets[root_signature]
    resolved[root_signature] = _best_supported_state(
        root_candidates,
        support_by_node[root_signature],
        priority_weights,
    )

    def visit(node, parent_state: str) -> None:
        for child in node.children:
            signature = node_signature(child)
            candidates = candidate_sets[signature]
            if parent_state in candidates:
                chosen = parent_state
            else:
                chosen = _best_supported_state(
                    candidates,
                    support_by_node[signature],
                    priority_weights,
                )
            resolved[signature] = chosen
            if not child.is_leaf():
                visit(child, chosen)

    visit(tree.root, resolved[root_signature])
    return resolved


def _state_support(
    node,
    states_by_taxon: dict[str, str],
    state_order: list[str],
) -> dict[str, int]:
    counts = dict.fromkeys(state_order, 0)
    for taxon in node_descendant_taxa(node):
        state = states_by_taxon.get(taxon)
        if state is not None:
            counts[state] += 1
    return counts


def _best_supported_state(
    candidate_states: list[str],
    support_counts: dict[str, int],
    priority_weights: dict[str, float],
) -> str:
    return max(
        sorted(candidate_states),
        key=lambda state: (
            support_counts.get(state, 0),
            priority_weights.get(state, 0.0),
            state,
        ),
    )
