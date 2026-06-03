from __future__ import annotations

from bijux_phylogenetics.ancestral.common import stable_value

from ..likelihood import GeographicExcludedTaxonRow
from .contracts import (
    GeographicSamplingBiasNodeRow,
    GeographicSamplingBiasSummary,
    GeographicSamplingBiasTransitionRow,
    GeographicSamplingCountRow,
)
from .weighting_policy import DOMINANT_REGION_THRESHOLD


def build_node_rows(
    baseline,
    weighted,
) -> list[GeographicSamplingBiasNodeRow]:
    weighted_by_node = {estimate.node: estimate for estimate in weighted.estimates}
    internal_estimates = [
        estimate for estimate in baseline.estimates if not estimate.is_tip
    ]
    root_node = max(
        internal_estimates,
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    ).node
    rows: list[GeographicSamplingBiasNodeRow] = []
    for estimate in internal_estimates:
        weighted_estimate = weighted_by_node[estimate.node]
        weighted_region = weighted_estimate.most_likely_state
        weighted_confidence = stable_value(
            weighted_estimate.state_probabilities[weighted_region]
        )
        unweighted_confidence = stable_value(
            estimate.state_probabilities[estimate.most_likely_state]
        )
        rows.append(
            GeographicSamplingBiasNodeRow(
                node=estimate.node,
                node_name=estimate.node_name,
                descendant_taxa=list(estimate.descendant_taxa),
                is_root=estimate.node == root_node,
                unweighted_region=estimate.most_likely_state,
                weighted_region=weighted_region,
                unweighted_confidence=unweighted_confidence,
                weighted_confidence=weighted_confidence,
                confidence_delta=stable_value(
                    weighted_confidence - unweighted_confidence
                ),
                changed=weighted_region != estimate.most_likely_state,
                unweighted_region_probabilities=dict(estimate.state_probabilities),
                weighted_region_probabilities=dict(
                    weighted_estimate.state_probabilities
                ),
            )
        )
    return rows


def build_transition_rows(
    baseline,
    weighted,
) -> list[GeographicSamplingBiasTransitionRow]:
    baseline_support = {
        (row.parent_node, row.child_node): row for row in baseline.support_rows
    }
    weighted_support = {
        (row.parent_node, row.child_node): row for row in weighted.support_rows
    }
    weighted_event_by_branch = {
        (event.parent_node, event.child_node): event for event in weighted.events
    }
    rows: list[GeographicSamplingBiasTransitionRow] = []
    for event in baseline.events:
        weighted_event = weighted_event_by_branch[(event.parent_node, event.child_node)]
        rows.append(
            GeographicSamplingBiasTransitionRow(
                parent_node=event.parent_node,
                child_node=event.child_node,
                child_descendant_taxa=event.child_node.split("|"),
                unweighted_source_region=event.source_state,
                unweighted_target_region=event.target_state,
                weighted_source_region=weighted_event.source_state,
                weighted_target_region=weighted_event.target_state,
                unweighted_transition=f"{event.source_state}->{event.target_state}",
                weighted_transition=(
                    f"{weighted_event.source_state}->{weighted_event.target_state}"
                ),
                unweighted_changed=event.changed,
                weighted_changed=weighted_event.changed,
                changed_by_weighting=(
                    event.source_state != weighted_event.source_state
                    or event.target_state != weighted_event.target_state
                ),
                unweighted_support=baseline_support[
                    (event.parent_node, event.child_node)
                ].support,
                weighted_support=weighted_support[
                    (event.parent_node, event.child_node)
                ].support,
            )
        )
    return rows


def build_summary(
    *,
    trait: str,
    taxon_column: str,
    model: str,
    internal_model: str,
    analyzed_taxon_count: int,
    excluded_taxon_count: int,
    weighting_mode: str,
    count_rows: list[GeographicSamplingCountRow],
    node_rows: list[GeographicSamplingBiasNodeRow],
    transition_rows: list[GeographicSamplingBiasTransitionRow],
    warning_count: int,
) -> GeographicSamplingBiasSummary:
    root_row = next(row for row in node_rows if row.is_root)
    dominant_row = max(count_rows, key=lambda row: (row.sample_fraction, row.region))
    weighted_dominant_row = max(
        count_rows,
        key=lambda row: (row.weighted_sample_fraction, row.region),
    )
    changed_node_count = sum(1 for row in node_rows if row.changed)
    changed_transition_count = sum(
        1
        for row in transition_rows
        if row.unweighted_transition != row.weighted_transition
    )
    return GeographicSamplingBiasSummary(
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        internal_model=internal_model,
        weighting_mode=weighting_mode,
        analyzed_taxon_count=analyzed_taxon_count,
        excluded_taxon_count=excluded_taxon_count,
        observed_region_count=len(count_rows),
        region_dominated=dominant_row.sample_fraction >= DOMINANT_REGION_THRESHOLD,
        dominant_region=dominant_row.region,
        dominant_region_fraction=dominant_row.sample_fraction,
        weighted_region_dominated=(
            weighted_dominant_row.weighted_sample_fraction >= DOMINANT_REGION_THRESHOLD
        ),
        weighted_dominant_region=weighted_dominant_row.region,
        weighted_dominant_region_fraction=weighted_dominant_row.weighted_sample_fraction,
        root_region_unweighted=root_row.unweighted_region,
        root_region_weighted=root_row.weighted_region,
        root_region_changed=root_row.changed,
        compared_internal_node_count=len(node_rows),
        changed_internal_node_count=changed_node_count,
        compared_transition_count=len(transition_rows),
        changed_transition_count=changed_transition_count,
        warning_count=warning_count,
    )


def build_exclusion_rows(audit) -> list[GeographicExcludedTaxonRow]:
    return [
        GeographicExcludedTaxonRow(
            taxon=row.taxon,
            raw_region=row.raw_state,
            normalized_region=row.normalized_state,
            reason=row.issue_code or "excluded",
            note=row.note,
        )
        for row in audit.rows
        if not row.included
    ]
